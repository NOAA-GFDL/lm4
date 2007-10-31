#define __DEBUG1__(x) write(*,'(a12,g)')#x,x
#define __DEBUG2__(x1,x2) write(*,'(99(a12,g))')#x1,x1,#x2,x2 
#define __DEBUG3__(x1,x2,x3) write(*,'(99(a12,g))')#x1,x1,#x2,x2,#x3,x3 
#define __DEBUG4__(x1,x2,x3,x4) write(*,'(99(a12,g))')#x1,x1,#x2,x2,#x3,x3,#x4,x4 
! ============================================================================
! vegetation (fast) biophysics
! ============================================================================
module vegetation_mod

use fms_mod, only: write_version_number, error_mesg, NOTE, file_exist, close_file, &
     open_namelist_file, check_nml_error, stdlog
use mpp_mod, only: mpp_sum, mpp_pe, mpp_root_pe
use time_manager_mod, only: time_type, time_type_to_real, get_date, operator(-)
use constants_mod,    only: tfreeze, rdgas, rvgas, hlv, hlf, cp_air, PI
use sphum_mod, only: qscomp
use nf_utils_mod, only: nfu_def_var, nfu_get_var_int, nfu_put_var_int

use vegn_tile_mod, only: vegn_tile_type, update_derived_vegn_data, &
     vegn_seed_demand, vegn_seed_supply, vegn_add_bliving
use soil_tile_mod, only: soil_tile_type, soil_ave_temp, soil_ave_theta
     
use land_constants_mod, only : NBANDS, BAND_VIS, d608, seconds_per_year
use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=), &
     get_elmt_indices
use land_tile_diag_mod, only : &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, diag_buff_type
use land_data_mod,      only : land_state_type, lnd
use land_tile_io_mod, only : create_tile_out_file, write_tile_data_r0d_fptr, &
     write_tile_data_i0d_fptr, read_tile_data_r0d_fptr, read_tile_data_i0d_fptr,&
     print_netcdf_error
use vegn_data_mod, only : SP_C4GRASS, read_vegn_data_namelist, &
     tau_drip_l, tau_drip_s, T_transp_min, cold_month_threshold, soil_carbon_depth_scale
use vegn_cohort_mod, only : vegn_cohort_type, vegn_phys_prog_type, &
     height_from_biomass, &
     vegn_data_heat_capacity, vegn_data_intrcptn_cap, &
     get_vegn_wet_frac, vegn_data_cover
use canopy_air_mod, only : cana_turbulence
     
use cohort_list_mod, only : insert_cohort, first_cohort, current_cohort, n_cohorts
use cohort_io_mod, only :  read_create_cohorts, create_cohort_dimension, &
     read_cohort_data_r0d_fptr,  read_cohort_data_i0d_fptr,&
     write_cohort_data_r0d_fptr, write_cohort_data_i0d_fptr
use land_debug_mod, only : is_watch_point, set_current_point
use vegn_radiation_mod, only : vegn_radiation_init, vegn_radiation
use vegn_photosynthesis_mod, only : vegn_photosynthesis_init, vegn_photosynthesis
use static_vegn_mod, only : static_vegn_init, static_vegn_end, static_vegn_override
use vegn_dynamics_mod, only : vegn_dynamics_init, vegn_carbon_int, vegn_growth, &
     vegn_daily_npp, vegn_phenology, vegn_biogeography
use vegn_disturbance_mod, only : vegn_nat_mortality, vegn_disturbance, update_fuel
implicit none
private

! ==== public interfaces =====================================================
public :: read_vegn_namelist
public :: vegn_init
public :: vegn_end

public :: vegn_get_cover
public :: vegn_radiation
public :: vegn_diffusion

public :: vegn_step_1
public :: vegn_step_2
public :: vegn_step_3

public :: update_vegn_slow
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: vegetation.F90,v 15.1.2.2 2007/09/16 22:14:23 slm Exp $', &
   tagname = '$Name: omsk_2007_10 $', &
   module_name = 'vegn'

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2               = .false.
real    :: init_Wl           = 0
real    :: init_Ws           = 0
real    :: init_Tv           = 288.
real    :: cpw = 1952.  ! specific heat of water vapor at constant pressure
real    :: clw = 4218.  ! specific heat of water (liquid)
real    :: csw = 2106.  ! specific heat of water (ice)
real    :: init_cohort_bl    = 0.05 ! initial biomass of leaves, kg C/m2
real    :: init_cohort_blv   = 0.0  ! initial biomass of labile store, kg C/m2
real    :: init_cohort_br    = 0.05 ! initial biomass of fine roots, kg C/m2
real    :: init_cohort_bsw   = 0.05 ! initial biomass of sapwood, kg C/m2
real    :: init_cohort_bwood = 0.05 ! initial biomass of heartwood, kg C/m2
real    :: init_cohort_cmc   = 0.0  ! initial intercepted water
character(32) :: rad_to_use = 'big-leaf' ! or 'two-stream'
character(32) :: snow_rad_to_use = 'ignore' ! or 'paint-leaves'
character(32) :: photosynthesis_to_use = 'simple' ! or 'leuning'
logical :: write_soil_carbon_restart = .FALSE. ! indicates whether to write
                        ! information for soil carbon acceleration
logical :: do_cohort_dynamics   = .TRUE. ! if true, do vegetation growth
logical :: do_patch_disturbance = .TRUE. ! 
logical :: do_phenology         = .TRUE. 
logical :: do_biogeography      = .TRUE.
logical :: do_seed_transport    = .TRUE.
namelist /vegn_nml/ &
    lm2, init_Wl, init_Ws, init_Tv, cpw, clw, csw, &
    init_cohort_bl, init_cohort_blv, init_cohort_br, init_cohort_bsw, &
    init_cohort_bwood, init_cohort_cmc, &
    rad_to_use, snow_rad_to_use, photosynthesis_to_use, &
    write_soil_carbon_restart, &
    do_cohort_dynamics, do_patch_disturbance, do_phenology, &
    do_biogeography, do_seed_transport
    
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
type(time_type) :: time ! *** NOT YET USED
real            :: delta_time      ! fast time step

! diagnostic field ids
integer :: id_vegn_type, id_temp, id_wl, id_ws, id_height, id_lai, id_sai, id_leaf_size, &
   id_root_density, id_root_zeta, id_rs_min, id_leaf_refl, id_leaf_tran,&
   id_leaf_emis, id_snow_crit, id_stomatal, id_an_op, id_an_cl, &
   id_bl, id_blv, id_br, id_bsw, id_bwood, id_species, id_status, &
   id_con_v_h, id_con_v_v
! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_vegn_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call read_vegn_data_namelist()

  call write_version_number(version, tagname)
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=vegn_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'vegn_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  if (mpp_pe() == mpp_root_pe()) then
     write (stdlog(), nml=vegn_nml)
  endif

  ! ---- initialize vegetation radiation options
  call vegn_radiation_init(rad_to_use, snow_rad_to_use)

  ! ---- initialize vegetation photosynthesis options
  call vegn_photosynthesis_init(photosynthesis_to_use)

end subroutine read_vegn_namelist


! ============================================================================
! initialize vegetation
subroutine vegn_init ( id_lon, id_lat, id_band )
  integer, intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in) :: id_lat  ! ID of land latitude (Y) axis
  integer, intent(in) :: id_band ! ID of spectral band axis

  ! ---- local vars
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce ! current and tail tile list elements
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  type(vegn_cohort_type), pointer :: cohort! pointer to initial cohort for cold-start
  integer :: n_accum
  integer :: nmn_acm

  module_is_initialized = .TRUE.

  ! ---- make module copy of time and calculate time step ------------------
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)


  ! ---- initialize vegn state ---------------------------------------------
  n_accum = 0
  nmn_acm = 0
  if (file_exist('INPUT/vegn1.res.nc')) then
     __NF_ASRT__(nf_open('INPUT/vegn1.res.nc',NF_NOWRITE,unit))
     ! read the cohort index and generate appropriate number of cohorts
     ! for each vegetation tile
     call read_create_cohorts(unit)
     
     ! read cohort data
     call read_cohort_data_r0d_fptr(unit, 'tv', cohort_tv_ptr )
     call read_cohort_data_r0d_fptr(unit, 'wl', cohort_wl_ptr )
     call read_cohort_data_r0d_fptr(unit, 'ws', cohort_ws_ptr )
     __NF_ASRT__(nf_close(unit))     

     __NF_ASRT__(nf_open('INPUT/vegn2.res.nc',NF_NOWRITE,unit))
     ! read global variables
     __NF_ASRT__(nfu_get_var_int(unit,'n_accum',n_accum))
     __NF_ASRT__(nfu_get_var_int(unit,'nmn_acm',nmn_acm))

     call read_cohort_data_i0d_fptr(unit, 'species', cohort_species_ptr )
     call read_cohort_data_r0d_fptr(unit, 'hite', cohort_height_ptr )
     call read_cohort_data_r0d_fptr(unit, 'bl', cohort_bl_ptr )
     call read_cohort_data_r0d_fptr(unit, 'blv', cohort_blv_ptr )
     call read_cohort_data_r0d_fptr(unit, 'br', cohort_br_ptr )
     call read_cohort_data_r0d_fptr(unit, 'bsw', cohort_bsw_ptr )
     call read_cohort_data_r0d_fptr(unit, 'bwood', cohort_bwood_ptr )
     call read_cohort_data_r0d_fptr(unit, 'bliving', cohort_bliving_ptr )
     call read_cohort_data_i0d_fptr(unit, 'status', cohort_status_ptr )
     call read_cohort_data_r0d_fptr(unit, 'npp_prev_day', cohort_npp_previous_day_ptr )

     call read_tile_data_r0d_fptr(unit,'age',vegn_age_ptr)
     call read_tile_data_r0d_fptr(unit,'fsc',vegn_fast_soil_C_ptr)
     call read_tile_data_r0d_fptr(unit,'ssc',vegn_slow_soil_C_ptr)
     call read_tile_data_r0d_fptr(unit,'fsc_pool',vegn_fsc_pool_ptr)
     call read_tile_data_r0d_fptr(unit,'fsc_rate',vegn_fsc_rate_ptr)
     call read_tile_data_r0d_fptr(unit,'ssc_pool',vegn_fsc_pool_ptr)
     call read_tile_data_r0d_fptr(unit,'ssc_rate',vegn_fsc_rate_ptr)
     ! monthly-mean values
     call read_tile_data_r0d_fptr(unit,'tc_av', vegn_tc_av_ptr)
     call read_tile_data_r0d_fptr(unit,'theta_av', vegn_theta_av_ptr)
     call read_tile_data_r0d_fptr(unit,'tsoil_av', vegn_tsoil_av_ptr)
     call read_tile_data_r0d_fptr(unit,'precip_av', vegn_precip_av_ptr)
     call read_tile_data_r0d_fptr(unit,'lambda', vegn_lambda_ptr)
     call read_tile_data_r0d_fptr(unit,'fuel', vegn_fuel_ptr)
     ! annual-mean values
     call read_tile_data_r0d_fptr(unit,'t_ann', vegn_t_ann_ptr)
     call read_tile_data_r0d_fptr(unit,'t_cold', vegn_t_ann_ptr)
     call read_tile_data_r0d_fptr(unit,'p_ann', vegn_p_ann_ptr)
     call read_tile_data_r0d_fptr(unit,'ncm', vegn_ncm_ptr)
     ! accumulated values for annual averaging
     call read_tile_data_r0d_fptr(unit,'t_ann_acm', vegn_t_ann_acm_ptr)
     call read_tile_data_r0d_fptr(unit,'t_cold_acm', vegn_t_cold_acm_ptr)
     call read_tile_data_r0d_fptr(unit,'p_ann_acm', vegn_p_ann_acm_ptr)
     call read_tile_data_r0d_fptr(unit,'ncm_acm', vegn_ncm_acm_ptr)
     __NF_ASRT__(nf_close(unit))     
  endif
  ! Go through all tiles and initilize the cohorts that have not been initilialized yet --
  ! this allows to read partial restarts. Also initialize accumulation counters to zero
  ! or the values from the restarts
  te  = tail_elmt(lnd%tile_map)
  ce = first_elmt(lnd%tile_map)
  do while(ce /= te)
     tile=>current_tile(ce)  ! get pointer to current tile
     ce=next_elmt(ce)       ! advance position to the next tile
     if (.not.associated(tile%vegn)) cycle

     tile%vegn%n_accum = n_accum
     tile%vegn%nmn_acm = nmn_acm

     if (n_cohorts(tile%vegn%cohorts)>0) cycle ! skip initialized tiles
     
     ! create and initialize cohorts for this vegetation tile
     ! for now, just create a new cohort with default values of biomasses
     allocate(cohort)
     cohort%prog%Wl = init_Wl
     cohort%prog%Ws = init_Ws
     cohort%prog%Tv = init_Tv
     
     !        cohort%species = SP_C4GRASS ! to be determined based on input later?
     cohort%species = tile%vegn%tag
     cohort%bl      = init_cohort_bl
     cohort%blv     = init_cohort_blv
     cohort%br      = init_cohort_br
     cohort%bsw     = init_cohort_bsw
     cohort%bwood   = init_cohort_bwood
     cohort%bliving = cohort%bl+cohort%br+cohort%blv+cohort%bsw
     cohort%npp_previous_day = 0.0
     cohort%height = height_from_biomass(cohort%bliving+cohort%bwood)
     cohort%status  = 0
     call insert_cohort(cohort,tile%vegn%cohorts)
     cohort=>NULL()
  enddo
  
  if (file_exist('INPUT/soil_carbon.res.nc')) then
     __NF_ASRT__(nf_open('INPUT/soil_carbon.res.nc',NF_NOWRITE,unit))
     call error_mesg('veg_data_init','reading soil_carbon restart',NOTE)
     call read_tile_data_r0d_fptr(unit,'asoil_in',vegn_asoil_in_ptr)
     call read_tile_data_r0d_fptr(unit,'fsc_in',vegn_fsc_in_ptr)
     call read_tile_data_r0d_fptr(unit,'ssc_in',vegn_ssc_in_ptr)
     __NF_ASRT__(nf_close(unit))     
  endif
  
  ! initialize carbon integrator
  call vegn_dynamics_init ( id_lon, id_lat, lnd%time, delta_time )

  ! initialize static vegetation
  call static_vegn_init ()
  call static_vegn_override ( lnd%time )

  ! initialize vegetation diagnostic fields
  call vegn_diag_init ( id_lon, id_lat, id_band, lnd%time )

  ! ---- diagnostic section
  ce = first_elmt(lnd%tile_map, is=lnd%is, js=lnd%js)
  te  = tail_elmt(lnd%tile_map)
  do while(ce /= te)
     tile => current_tile(ce)
     ce=next_elmt(ce)     
     if (.not.associated(tile%vegn)) cycle ! skip non-vegetetion tiles
     ! send the data
     call send_tile_data(id_vegn_type,  real(tile%vegn%tag), tile%diag)
  enddo

end subroutine vegn_init

! ============================================================================
subroutine vegn_diag_init ( id_lon, id_lat, id_band, time )
  integer        , intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer        , intent(in) :: id_lat  ! ID of land latitude (Y) axis
  integer        , intent(in) :: id_band ! ID of spectral band axis
  type(time_type), intent(in) :: time    ! initial time for diagnostic fields
  
  id_vegn_type = register_tiled_static_field ( module_name, 'vegn_type',  &
       (/id_lon,id_lat/), 'vegetation type', missing_value=-1.0 )

  id_temp = register_tiled_diag_field ( module_name, 'temp',  &
       (/id_lon,id_lat/), time, 'canopy temperature', 'degK', missing_value=-1.0 )
  id_wl = register_tiled_diag_field ( module_name, 'wl',  &
       (/id_lon,id_lat/), time, 'canopy liquid water content', 'kg/m2', missing_value=-1.0 )
  id_ws = register_tiled_diag_field ( module_name, 'ws',  &
       (/id_lon,id_lat/), time, 'canopy solid water content', 'kg/m2', missing_value=-1.0 )

  id_height = register_tiled_diag_field ( module_name, 'height',  &
       (/id_lon,id_lat/), time, 'vegetation height', 'm', missing_value=-1.0 )
  id_lai    = register_tiled_diag_field ( module_name, 'lai',  &
       (/id_lon,id_lat/), time, 'leaf area index', 'm2/m2', missing_value=-1.0 )
  id_sai    = register_tiled_diag_field ( module_name, 'sai',  &
       (/id_lon,id_lat/), time, 'stem area index', 'm2/m2', missing_value=-1.0 )
  id_leaf_size = register_tiled_diag_field ( module_name, 'leaf_size',  &
       (/id_lon,id_lat/), time, missing_value=-1.0 )
  id_root_density = register_tiled_diag_field ( module_name, 'root_density',  &
       (/id_lon,id_lat/), time, 'total biomass below ground', 'kg/m2', missing_value=-1.0 )
  id_root_zeta = register_tiled_diag_field ( module_name, 'root_zeta',  &
       (/id_lon,id_lat/), time, 'e-folding depth of root biomass', 'm',missing_value=-1.0 )
  id_rs_min = register_tiled_diag_field ( module_name, 'rs_min',  &
       (/id_lon,id_lat/), time, missing_value=-1.0 )
  id_leaf_refl = register_tiled_diag_field ( module_name, 'leaf_refl',  &
       (/id_lon,id_lat,id_band/), time, 'reflectivity of leaf', missing_value=-1.0 )
  id_leaf_tran = register_tiled_diag_field ( module_name, 'leaf_tran',  &
       (/id_lon,id_lat,id_band/), time, 'transmittance of leaf', missing_value=-1.0 )
  id_leaf_emis = register_tiled_diag_field ( module_name, 'leaf_emis',  &
       (/id_lon,id_lat/), time, 'leaf_emissivity', missing_value=-1.0 )
  id_snow_crit = register_tiled_diag_field ( module_name, 'snow_crit',  &
       (/id_lon,id_lat/), time, missing_value=-1.0 )
  id_stomatal = register_tiled_diag_field ( module_name, 'stomatal_cond',  &
       (/id_lon,id_lat/), time, 'vegetation stomatal conductance', missing_value=-1.0 )
  id_an_op = register_tiled_diag_field ( module_name, 'An_op',  &
       (/id_lon,id_lat/), time, 'net photosynthesis with open stomata', &
       '(mol CO2)(m2 of leaf)^-1 year^-1', missing_value=-1e20 )
  id_an_cl = register_tiled_diag_field ( module_name, 'An_cl',  &
       (/id_lon,id_lat/), time, 'net photosynthesis with closed stomata', &
       '(mol CO2)(m2 of leaf)^-1 year^-1', missing_value=-1e20 )
  
  id_bl = register_tiled_diag_field ( module_name, 'bl',  &
       (/id_lon,id_lat/), time, 'biomass of leaves', 'kg C/m2', missing_value=-1.0 )
  id_blv = register_tiled_diag_field ( module_name, 'blv',  &
       (/id_lon,id_lat/), time, 'biomass in labile store', 'kg C/m2', missing_value=-1.0 )
  id_br = register_tiled_diag_field ( module_name, 'br',  &
       (/id_lon,id_lat/), time, 'biomass of fine roots', 'kg C/m2', missing_value=-1.0 )
  id_bsw = register_tiled_diag_field ( module_name, 'bsw',  &
       (/id_lon,id_lat/), time, 'biomass of sapwood', 'kg C/m2', missing_value=-1.0 )
  id_bwood = register_tiled_diag_field ( module_name, 'bwood',  &
       (/id_lon,id_lat/), time, 'biomass of heartwood', 'kg C/m2', missing_value=-1.0 )

  id_species = register_tiled_diag_field ( module_name, 'species',  &
       (/id_lon,id_lat/), time, 'species type', missing_value=-1.0 )
  id_status = register_tiled_diag_field ( module_name, 'status',  &
       (/id_lon,id_lat/), time, 'status of leaves', missing_value=-1.0 )

  id_con_v_h = register_tiled_diag_field ( module_name, 'con_v_h', (/id_lon,id_lat/), &
       time, 'conductance for sensible heat between canopy and canopy air', &
       'm/s', missing_value=-1.0 )
  id_con_v_v = register_tiled_diag_field ( module_name, 'con_v_v', (/id_lon,id_lat/), &
       time, 'conductance for water vapor between canopy and canopy air', &
       'm/s', missing_value=-1.0 )

end subroutine


! ============================================================================
! write restart file and release memory
subroutine vegn_end (tile_dim_length)
  integer, intent(in) :: tile_dim_length ! length of tile dimension in the 
                                         ! output file
  ! ---- local vars ----------------------------------------------------------
  integer :: unit ! restart file unit 
  integer :: ierr
  logical :: restart_created ! flag indicating that the restart file was created
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile

  module_is_initialized =.FALSE.

  ! ---- write restart file --------------------------------------------------

  ! create output file, including internal structure necessary for tile output
  ! NOTE: restart is not created if the number of the tiles in the domain is  
  ! equal to zero; in this case nothing needs to be written
  call error_mesg('vegn_end','writing NetCDF restart',NOTE)
  call create_tile_out_file(unit,'RESTART/vegn1.res.nc', &
          lnd%glon*180.0/PI, lnd%glat*180/PI, vegn_tile_exists, tile_dim_length, &
          created=restart_created)
  ! create compressed dimension for vegetation cohorts -- must be called even
  ! if restart has not been created, because it calls mpp_max and that should 
  ! be called on all PEs to work
  call create_cohort_dimension(unit)

  if (restart_created) then
     call write_cohort_data_r0d_fptr(unit,'tv',cohort_tv_ptr,'vegetation temperature','degrees_K')
     call write_cohort_data_r0d_fptr(unit,'wl',cohort_wl_ptr,'vegetation liquid water content','kg/m2')
     call write_cohort_data_r0d_fptr(unit,'ws',cohort_ws_ptr,'vegetation solid water content','kg/m2')
     ! close output file
     __NF_ASRT__(nf_close(unit))
  endif

  call create_tile_out_file(unit,'RESTART/vegn2.res.nc', &
          lnd%glon*180.0/PI, lnd%glat*180/PI, vegn_tile_exists, tile_dim_length, &
          created=restart_created)
  ! create compressed dimension for vegetation cohorts -- see note above
  call create_cohort_dimension(unit)

  if (restart_created) then
     ! store global variables
     ierr = nf_redef(unit)
     __NF_ASRT__(nfu_def_var(unit,'n_accum',NF_INT,long_name='number of accumulated steps'))
     __NF_ASRT__(nfu_def_var(unit,'nmn_acm',NF_INT,long_name='number of accumulated months'))
     ierr = nf_enddef(unit)
     ! find first tile
     ce = first_elmt(lnd%tile_map) ; te = tail_elmt(lnd%tile_map)
     do while ( ce /= te )
        tile => current_tile(ce) ; ce=next_elmt(ce)
        if(associated(tile%vegn)) exit
     enddo
     __NF_ASRT__(nfu_put_var_int(unit,'n_accum',tile%vegn%n_accum))
     __NF_ASRT__(nfu_put_var_int(unit,'nmn_acm',tile%vegn%nmn_acm))
     
     call write_cohort_data_i0d_fptr(unit,'species', cohort_species_ptr, 'vegetation species')
     call write_cohort_data_r0d_fptr(unit,'hite', cohort_height_ptr, 'vegetation height','m')
     call write_cohort_data_r0d_fptr(unit,'bl', cohort_bl_ptr, 'biomass of leaves per individual','kg C/m2')
     call write_cohort_data_r0d_fptr(unit,'blv', cohort_blv_ptr, 'biomass of virtual leaves (labile store) per individual','kg C/m2')
     call write_cohort_data_r0d_fptr(unit,'br', cohort_br_ptr, 'biomass of fine roots per individual','kg C/m2')
     call write_cohort_data_r0d_fptr(unit,'bsw', cohort_bsw_ptr, 'biomass of sapwood per individual','kg C/m2')
     call write_cohort_data_r0d_fptr(unit,'bwood', cohort_bwood_ptr, 'biomass of heartwood per individual','kg C/m2')
     call write_cohort_data_r0d_fptr(unit,'bliving', cohort_bliving_ptr, 'total living biomass per individual','kg C/m2')
!     call write_cohort_data_r0d_fptr(unit,'tleaf', cohort_tleaf_ptr, 'leaf temperature','degK')
     call write_cohort_data_i0d_fptr(unit,'status', cohort_status_ptr, 'leaf status')
!     call write_cohort_data_r0d_fptr(unit,'intercept_l', cohort_cmc_ptr, 'intercepted water per cohort','kg/m2')
     call write_cohort_data_r0d_fptr(unit,'npp_prev_day', cohort_npp_previous_day_ptr, 'previous day NPP','kg C/(m2 year)')

     call write_tile_data_r0d_fptr(unit,'age',vegn_age_ptr,'vegetation age', 'yr')
     call write_tile_data_r0d_fptr(unit,'fsc',vegn_fast_soil_C_ptr,'fast soil carbon', 'kg C/m2')
     call write_tile_data_r0d_fptr(unit,'ssc',vegn_slow_soil_C_ptr,'slow soil carbon', 'kg C/m2')
     call write_tile_data_r0d_fptr(unit,'fsc_pool',vegn_fsc_pool_ptr,'intermediate pool for fast soil carbon input', 'kg C/m2')
     call write_tile_data_r0d_fptr(unit,'fsc_rate',vegn_fsc_rate_ptr,'conversion rate of fsc_pool to fast soil carbon', 'kg C/(m2 yr)')
     call write_tile_data_r0d_fptr(unit,'ssc_pool',vegn_fsc_pool_ptr,'intermediate pool for slow soil carbon input', 'kg C/m2')
     call write_tile_data_r0d_fptr(unit,'ssc_rate',vegn_fsc_rate_ptr,'conversion rate of ssc_pool to slow soil carbon', 'kg C/(m2 yr)')

     ! monthly-mean values
     call write_tile_data_r0d_fptr(unit,'tc_av', vegn_tc_av_ptr,'average canopy air temperature','degK')
     call write_tile_data_r0d_fptr(unit,'theta_av', vegn_theta_av_ptr,'average soil moisture')
     call write_tile_data_r0d_fptr(unit,'tsoil_av', vegn_tsoil_av_ptr,'average bulk soil temperature for soil carbon','degK')
     call write_tile_data_r0d_fptr(unit,'precip_av', vegn_precip_av_ptr,'average total precipitation','kg/(m2 s)')
     call write_tile_data_r0d_fptr(unit,'lambda', vegn_lambda_ptr,'dryness parameter')
     call write_tile_data_r0d_fptr(unit,'fuel', vegn_fuel_ptr,'fuel density','kg C/m2')
     ! annual-mean values
     call write_tile_data_r0d_fptr(unit,'t_ann', vegn_t_ann_ptr,'average annual canopy air temperature','degK')
     call write_tile_data_r0d_fptr(unit,'t_cold', vegn_t_ann_ptr,'average canopy air temperature of coldest month','degK')
     call write_tile_data_r0d_fptr(unit,'p_ann', vegn_p_ann_ptr,'average annual precipitation','kg/(m2 s)')
     call write_tile_data_r0d_fptr(unit,'ncm', vegn_ncm_ptr,'number of cold months')
     ! accumulated values for annual averaging
     call write_tile_data_r0d_fptr(unit,'t_ann_acm', vegn_t_ann_acm_ptr,'accumulated annual canopy air temperature','degK')
     call write_tile_data_r0d_fptr(unit,'t_cold_acm', vegn_t_cold_acm_ptr,'accumulated temperature of coldest month','degK')
     call write_tile_data_r0d_fptr(unit,'p_ann_acm', vegn_p_ann_acm_ptr,'accumulated precipitation','kg/(m2 s)')
     call write_tile_data_r0d_fptr(unit,'ncm_acm', vegn_ncm_acm_ptr,'accumulated number of cold months')
     

     __NF_ASRT__(nf_close(unit))
  endif

  if (write_soil_carbon_restart) then
     call create_tile_out_file(unit,'RESTART/soil_carbon.res.nc', &
          lnd%glon*180.0/PI, lnd%glat*180/PI, vegn_tile_exists, tile_dim_length, created=restart_created)
     if (restart_created) then
        call write_tile_data_r0d_fptr(unit,'asoil_in',vegn_asoil_in_ptr,'aerobic activity modifier', 'unitless')
        call write_tile_data_r0d_fptr(unit,'fsc_in',vegn_fsc_in_ptr,'fast soil carbon input', 'kg C/m2')
        call write_tile_data_r0d_fptr(unit,'ssc_in',vegn_ssc_in_ptr,'slow soil carbon input', 'kg C/m2')
        __NF_ASRT__(nf_close(unit))
     endif
  endif

  ! finalize static vegetation, if necessary
  call static_vegn_end ()

end subroutine vegn_end


! ============================================================================
subroutine vegn_get_cover(vegn, snow_depth, vegn_cover)
  type(vegn_tile_type), intent(in)  :: vegn
  real,                 intent(in)  :: snow_depth
  real,                 intent(out) :: vegn_cover

  type(vegn_cohort_type), pointer :: cohort
  real :: vegn_cover_snow_factor

  ! get the pointer to the first (and, currently, the only) cohort
  cohort => current_cohort(first_cohort(vegn%cohorts))

  call vegn_data_cover(cohort, snow_depth, vegn_cover, vegn_cover_snow_factor)
  
end subroutine vegn_get_cover


! ============================================================================
subroutine vegn_diffusion ( vegn, vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf)
  type(vegn_tile_type), intent(in) :: vegn
  real,                intent(out) :: &
       vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf
  
  type(vegn_cohort_type), pointer :: cohort

  ! get the pointer to the first (and, currently, the only) cohort
  cohort => current_cohort(first_cohort(vegn%cohorts))

  vegn_cover  = cohort%cover
  vegn_lai    = cohort%lai
  vegn_sai    = cohort%sai
  vegn_height = cohort%height
  vegn_d_leaf = cohort%leaf_size

end subroutine vegn_diffusion


! ============================================================================
subroutine vegn_step_1 ( vegn, diag, &
        p_surf, ustar, drag_q, &
        SWdn, RSv, precip_l, precip_s, &
        land_d, land_z0s, land_z0m, grnd_z0s, &
        soil_beta, soil_water_supply, soil_uptake_T, &
        cana_T, cana_q, cana_co2, &
        ! output
        con_g_h, con_g_v, & ! aerodynamic conductance between canopy air and canopy, for heat and vapor flux
        vegn_T,vegn_Wl,  vegn_Ws,           & ! temperature, water and snow mass of the canopy
        vegn_ifrac,                         & ! intercepted fraction of liquid and frozen precipitation
        vegn_lai,                           & ! leaf area index
        drip_l, drip_s,                     & ! water and snow drip rate from precipitation, kg/(m2 s)
        vegn_hcap,                          & ! vegetation heat capacity
        Hv0,   DHvDTv,   DHvDTc,            & ! sens heat flux
        Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
        Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
        Efi0,  DEfiDTv,  DEfiDqc,  DEfiDwl,  DEfiDwf  ) ! sublimation of intercepted snow
  type(vegn_tile_type), intent(inout) :: vegn ! vegetation data
  type(diag_buff_type), intent(inout) :: diag ! diagnostic buffer
  real, intent(in) :: &
       p_surf,    & ! surface pressure, N/m2
       ustar,     & ! friction velocity, m/s
       drag_q,    & ! bulk drag coefficient for specific humidity
       SWdn(NBANDS), & ! downward SW radiation at the top of the canopy, W/m2
       RSv (NBANDS), & ! net SW radiation balance of the canopy, W/m2
       precip_l, precip_s, & ! liquid and solid precipitation rates, kg/(m2 s)
       land_d, land_z0s, land_z0m, & ! land displacement height and roughness, m
       grnd_z0s, & ! roughness of ground surface (including snow effect)
       soil_beta, & ! relative water availability
       soil_water_supply, & ! supply of water to roots per unit active root biomass, kg/m2
       soil_uptake_T, & ! temperature of the water uptake from soil
       cana_T,    & ! temperature of canopy air
       cana_q,    & ! specific humidity of canopy air
       cana_co2     ! co2 concentration in canopy air
  ! output -- coefficients of linearized expressions for fluxes
  real, intent(out) ::   &
       vegn_T,vegn_Wl,  vegn_Ws,& ! temperature, water and snow mass of the canopy
       vegn_ifrac, & ! intercepted fraction of liquid and frozen precipitation
       vegn_lai, & ! vegetation leaf area index
       drip_l, drip_s, & ! water and snow drip rate from precipitation, kg/(m2 s)
       vegn_hcap, & ! total vegetation heat capacity, including intercepted water and snow
       con_g_h, con_g_v, & ! aerodynamic conductance between ground and canopy air
       Hv0,   DHvDTv,   DHvDTc, & ! sens heat flux
       Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
       Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
       Efi0,  DEfiDTv,  DEfiDqc,  DEfiDwl,  DEfiDwf    ! sublimation of intercepted snow
  
  ! ---- local vars 
  real :: &
       ft,DftDwl,DftDwf, & ! fraction of canopy not covered by intercepted water/snow, and its' 
                    ! derivatives w.r.t. intercepted water masses 
       fw,DfwDwl,DfwDwf, & ! fraction of canopy covered by intercepted water, and its' 
                    ! derivatives w.r.t. intercepted water masses 
       fs,DfsDwl,DfsDwf, & ! fraction of canopy covered by intercepted snow, and its' 
                    ! derivatives w.r.t. intercepted water masses
       stomatal_cond, & ! integral stomatal conductance of canopy
       con_v_h, con_v_v, & ! aerodyn. conductance between canopy and CAS, for heat and vapor
       total_cond, &! overall conductance from inside stomata to canopy air 
       qvsat,     & ! sat. specific humidity at the leaf T
       DqvsatDTv, & ! derivative of qvsat w.r.t. leaf T
       rho,       & ! density of canopy air
       photosynt, & ! photosynthesis
       photoresp    ! photo-respiration
  type(vegn_cohort_type), pointer :: cohort
  
  ! get the pointer to the first (and, currently, the only) cohort
  cohort => current_cohort(first_cohort(vegn%cohorts))

  ! calculate the fractions of intercepted precipitation
  vegn_ifrac = cohort%cover

  ! get the lai
  vegn_lai = cohort%lai

  ! calculate the aerodynamic conductance coeficients
  call cana_turbulence(ustar, &
     cohort%cover, cohort%height, cohort%lai, cohort%sai, cohort%leaf_size, &
     land_d, land_z0m, land_z0s, grnd_z0s, &
     con_v_h, con_v_v, con_g_h, con_g_v)

  ! calculate the vegetation photosynthesis and associated stomatal conductance
  call vegn_photosynthesis ( vegn, &
     SWdn(BAND_VIS), RSv(BAND_VIS), cana_q, cana_co2, p_surf, drag_q, &
     soil_beta, soil_water_supply/delta_time, &
     stomatal_cond, photosynt, photoresp )

  call get_vegn_wet_frac ( cohort, fw, DfwDwl, DfwDwf, fs, DfsDwl, DfsDwf )
  ! transpiring fraction and its derivatives
  ft     = 1 - fw - fs
  DftDwl = - DfwDwl - DfsDwl
  DftDwf = - DfwDwf - DfsDwf
  call qscomp(cohort%prog%Tv, p_surf, qvsat, DqvsatDTv)

  rho = p_surf/(rdgas*cana_T *(1+d608*cana_q))
  
  ! get the vegetation temperature
  vegn_T  =  cohort%prog%Tv
  ! get the amount of intercepted water and snow
  vegn_Wl =  cohort%prog%Wl
  vegn_Ws =  cohort%prog%Ws
  ! calculate the drip rates
  drip_l  = max(vegn_Wl,0.0)/tau_drip_l
  drip_s  = max(vegn_Ws,0.0)/tau_drip_s
  ! correct the drip rates so that the amount of water and snow accumulated over time step 
  ! is no larger then the canopy water-holding capacity
  drip_l = max((vegn_Wl+precip_l*delta_time*vegn_ifrac-cohort%W_max)/delta_time,drip_l)
  drip_s = max((vegn_Ws+precip_s*delta_time*vegn_ifrac-cohort%W_max)/delta_time,drip_s)

  ! calculate the total heat capacity
  call vegn_data_heat_capacity (cohort, vegn_hcap)
  if(cohort%lai == 0 ) vegn_hcap = 1.0 ! set some heat capacity for degenerate case
  vegn_hcap = vegn_hcap + clw*cohort%prog%Wl + csw*cohort%prog%Ws
  ! calculate the coefficient of sensible heat flux linearization
  Hv0     =  2*rho*cp_air*con_v_h*(cohort%prog%Tv - cana_T)
  DHvDTv  =  2*rho*cp_air*con_v_h
  DHvDTc  = -2*rho*cp_air*con_v_h
  ! calculate the coefficients of the transpiration linearization
  if(con_v_v==0.and.stomatal_cond==0) then
     total_cond = 0.0
  else
     total_cond = stomatal_cond*con_v_v/(stomatal_cond+con_v_v)
  endif

  if(qvsat>cana_q)then
     ! flux is directed from the surface: transpiration is possible, and the
     ! evaporation of intercepted water depends on the fraction of wet/snow
     ! covered canopy.

     ! prohibit transpiration if leaf temperature below some predefined minimum
     ! typically (268K, but check namelist)
     if(cohort%prog%Tv < T_transp_min) total_cond = 0 
     ! calculate the transpiration linearization coefficients
     Et0     =  rho*total_cond*ft*(qvsat - cana_q)
     DEtDTv  =  rho*total_cond*ft*DqvsatDTv
     DEtDqc  = -rho*total_cond*ft
     DEtDwl  =  rho*total_cond*DftDwl*(qvsat - cana_q)
     DEtDwf  =  rho*total_cond*DftDwf*(qvsat - cana_q)
     ! calculate the coefficients of the intercepted liquid evaporation linearization
     Eli0    =  rho*con_v_v*fw*(qvsat - cana_q)
     DEliDTv =  rho*con_v_v*fw*DqvsatDTv
     DEliDqc = -rho*con_v_v*fw
     DEliDwl =  rho*con_v_v*DfwDwl*(qvsat-cana_q)
     DEliDwf =  rho*con_v_v*DfwDwf*(qvsat-cana_q)
     ! calculate the coefficients of the intercepted snow evaporation linearization
     Efi0    =  rho*con_v_v*fs*(qvsat - cana_q)
     DEfiDTv =  rho*con_v_v*fs*DqvsatDTv
     DEfiDqc = -rho*con_v_v*fs
     DEfiDwl =  rho*con_v_v*DfsDwl*(qvsat-cana_q)
     DEfiDwf =  rho*con_v_v*DfsDwf*(qvsat-cana_q)
  else
     ! Flux is ditected TOWARD the surface: no transpration (assuming plants do not
     ! take water through stomata), and condensation does not depend on the fraction
     ! of wet canopy -- dew formation occurs on the entire surface

     ! prohibit transpiration:
     Et0     = 0
     DEtDTv  = 0; DEtDwl = 0; DEtDwf = 0;
     DEtDqc  = 0
     ! calculate dew or frost formation rates, depending on the temperature
     Eli0    = 0; Efi0    = 0
     DEliDTv = 0; DEfiDTv = 0
     DEliDqc = 0; DEfiDqc = 0
     DEliDwl = 0; DEfiDwl = 0
     DEliDwf = 0; DEfiDwf = 0
     ! calculate the coefficients of the intercepted liquid condensation linearization
     if(vegn_T >= tfreeze) then
        Eli0    =  rho*con_v_v*(qvsat - cana_q)
        DEliDTv =  rho*con_v_v*DqvsatDTv
        DEliDqc = -rho*con_v_v
     else
        ! calculate the coefficients of the intercepted snow condensation linearization
        Efi0    =  rho*con_v_v*(qvsat - cana_q)
        DEfiDTv =  rho*con_v_v*DqvsatDTv
        DEfiDqc = -rho*con_v_v
     endif
  endif
  ! ---- diagnostic section
  call send_tile_data(id_stomatal, stomatal_cond, diag)
  call send_tile_data(id_an_op, cohort%An_op, diag)
  call send_tile_data(id_an_cl, cohort%An_cl, diag)
  call send_tile_data(id_con_v_h, con_v_h, diag)
  call send_tile_data(id_con_v_v, con_v_v, diag)

end subroutine vegn_step_1


! ============================================================================
subroutine vegn_step_2 ( vegn, diag, &
     delta_Tv, delta_wl, delta_wf, &
     vegn_melt, &
     vegn_ovfl_l,  vegn_ovfl_s,  & ! overflow of liquid and solid water from the canopy, kg/(m2 s)
     vegn_ovfl_Hl, vegn_ovfl_Hs, & ! heat flux carried from canopy by overflow, W/(m2 s)
     vegn_LMASS, vegn_FMASS, vegn_HEAT )
! ============================================================================
! Given the surface solution, substitute it back into the vegetation
! equations to determine corresponding new vegetation state.
! Perform balance of intercepted water.
! ----------------------------------------------------------------------------

  ! ---- arguments -----------------------------------------------------------
  type(vegn_tile_type) , intent(inout) :: vegn
  type(diag_buff_type) , intent(inout) :: diag
  real, intent(in) :: &
       delta_Tv, & ! change in vegetation temperature, degK
       delta_wl, & ! change in intercepted liquid water mass, kg/m2
       delta_wf    ! change in intercepted frozen water mass, kg/m2 
  real, intent(out) :: &
       vegn_melt, &
       vegn_ovfl_l,   vegn_ovfl_s,   & ! overflow of liquid and solid water from the canopy
       vegn_ovfl_Hl, vegn_ovfl_Hs, & ! heat flux from canopy due to overflow
       vegn_LMASS, &! mass of liquid water in canopy, kg/m2
       vegn_FMASS, &! mass of frozen water in canopy, kg/m2
       vegn_HEAT    ! total heat content of canopy, J/m2

  ! ---- local variables
  real :: &
     vegn_W_max, &  ! max. possible amount of water (solid+liquid) in the canopy
     mcv, &
     cap0, melt_per_deg, &
     Wl, Ws,  & ! positively defined amounts of water and snow on canopy
     overflow   ! excess of total canopy water over allowed maximum at the end of 
                ! implicit time step
  type(vegn_cohort_type), pointer :: cohort
  
  if (is_watch_point()) then
     write(*,*)'#### vegn_step_2 input ####'
     __DEBUG3__(delta_Tv, delta_wl, delta_wf)
     __DEBUG1__(cohort%prog%Tv)
  endif

  ! get the pointer to the first (and, currently, the only) cohort
  cohort => current_cohort(first_cohort(vegn%cohorts))

  ! update vegetation state
  cohort%prog%Tv = cohort%prog%Tv + delta_Tv
  cohort%prog%Wl = cohort%prog%Wl + delta_wl
  cohort%prog%Ws = cohort%prog%Ws + delta_wf 

  call vegn_data_intrcptn_cap(cohort, vegn_W_max)
  call vegn_data_heat_capacity(cohort, mcv)


  ! ---- update for evaporation and interception -----------------------------
  cap0 = mcv + clw*cohort%prog%Wl + csw*cohort%prog%Ws

  if(is_watch_point()) then
     write (*,*)'#### vegn_step_2 #### 1'
     __DEBUG1__(cap0)
     __DEBUG1__(cohort%prog%Tv)
     __DEBUG2__(cohort%prog%Wl, cohort%prog%Ws)
  endif
  ! melt on the vegetation should probably be prohibited altogether, since
  ! the amount of melt or freeze calculated this way is severly underestimated 
  ! (depending on the overall vegetation heat capacity) which leads to extended 
  ! periods when the canopy temperature is fixed at freezing point.
  if (lm2) then 
     vegn_melt = 0
  else
     ! ---- freeze/melt of intercepted water
     ! heat capacity of leaf + intercepted water/snow _can_ go below zero if the 
     ! total water content goes below zero as a result of implicit time step.
     ! If it does, we just prohibit melt, setting it to zero.
     if(cap0 > 0)then
        melt_per_deg = cap0 / hlf
        if (cohort%prog%Ws>0 .and. cohort%prog%Tv>tfreeze) then
           vegn_melt =  min(cohort%prog%Ws, (cohort%prog%Tv-tfreeze)*melt_per_deg)
        else if (cohort%prog%Wl>0 .and. cohort%prog%Tv<tfreeze) then
           vegn_melt = -min(cohort%prog%Wl, (tfreeze-cohort%prog%Tv)*melt_per_deg)
        else
           vegn_melt = 0
        endif
        cohort%prog%Ws = cohort%prog%Ws - vegn_melt
        cohort%prog%Wl = cohort%prog%Wl + vegn_melt
        if (vegn_melt/=0) &
             cohort%prog%Tv = tfreeze + (cap0*(cohort%prog%Tv-tfreeze) - hlf*vegn_melt) &
             / ( cap0 + (clw-csw)*vegn_melt )
        vegn_melt = vegn_melt / delta_time
     else
        vegn_melt = 0
     endif
  endif

  if(is_watch_point()) then
     write (*,*)'#### vegn_step_2 #### 1'
     __DEBUG1__(cap0)
     __DEBUG1__(cohort%prog%Tv)
     __DEBUG3__(vegn_melt, cohort%prog%Wl, cohort%prog%Ws)
  endif

  ! ---- update for overflow -------------------------------------------------
  Wl = max(cohort%prog%Wl,0.0); Ws = max(cohort%prog%Ws,0.0)
  overflow = max (0.,Wl+Ws-vegn_W_max)/delta_time
  if (overflow>0) then
     vegn_ovfl_l = overflow*Wl/(Wl+Ws) ; vegn_ovfl_s = overflow-vegn_ovfl_l
  else
     vegn_ovfl_l = 0                   ; vegn_ovfl_s = 0
  endif
  vegn_ovfl_Hl = clw*vegn_ovfl_l*(cohort%prog%Tv-tfreeze)
  vegn_ovfl_Hs = csw*vegn_ovfl_s*(cohort%prog%Tv-tfreeze)

  cohort%prog%Wl = cohort%prog%Wl - vegn_ovfl_l*delta_time
  cohort%prog%Ws = cohort%prog%Ws - vegn_ovfl_s*delta_time
  vegn_LMASS = cohort%prog%Wl
  vegn_FMASS = cohort%prog%Ws
  vegn_HEAT  = (mcv + clw*cohort%prog%Wl+ csw*cohort%prog%Ws)*(cohort%prog%Tv-tfreeze)

  if(is_watch_point()) then
     write(*,*)'#### vegn_step_2 output #####'
     __DEBUG3__(vegn_melt, vegn_ovfl_l, vegn_ovfl_s)
     __DEBUG2__(vegn_ovfl_Hl,vegn_ovfl_Hs)
     __DEBUG3__(vegn_LMASS, vegn_FMASS, vegn_HEAT)
  endif

  ! ---- diagnostic section
  call send_tile_data(id_temp,   cohort%prog%Tv, diag)
  call send_tile_data(id_wl,     cohort%prog%Wl, diag)
  call send_tile_data(id_ws,     cohort%prog%Ws, diag)

  call send_tile_data(id_height, cohort%height, diag)
  call send_tile_data(id_lai, cohort%lai, diag)
  call send_tile_data(id_sai, cohort%sai, diag)
  call send_tile_data(id_leaf_size, cohort%leaf_size, diag)
  call send_tile_data(id_root_density, cohort%root_density, diag)
  call send_tile_data(id_root_zeta, cohort%root_zeta, diag)
  call send_tile_data(id_rs_min, cohort%rs_min, diag)
  call send_tile_data(id_leaf_refl, cohort%leaf_refl, diag)
  call send_tile_data(id_leaf_tran, cohort%leaf_tran, diag)
  call send_tile_data(id_leaf_emis, cohort%leaf_emis, diag)
  call send_tile_data(id_snow_crit, cohort%snow_crit, diag)

  ! slow variables -- may be need to be moved to update_vegn_slow later 
  call send_tile_data(id_bl, cohort%bl, diag)
  call send_tile_data(id_blv, cohort%blv, diag)
  call send_tile_data(id_br, cohort%br, diag)
  call send_tile_data(id_bsw, cohort%bsw, diag)
  call send_tile_data(id_bwood, cohort%bwood, diag)
  call send_tile_data(id_species, real(cohort%species), diag)
  call send_tile_data(id_status, real(cohort%status), diag)
  
end subroutine vegn_step_2


! ============================================================================
! do the vegetation calculations that require updated (end-of-timestep) values 
! of progtnostic land variables
subroutine vegn_step_3(vegn, soil, cana_T, precip, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(in)    :: soil
  real, intent(in) :: cana_T ! canopy temperature, deg K
  real, intent(in) :: precip ! total (rain+snow) precipitation, kg/(m2 s)
  type(diag_buff_type), intent(inout) :: diag
  
  ! ---- local vars
  real :: tsoil ! average temperature of soil for soil carbon decomposition, deg K
  real :: theta ! average soil wetness, unitless
  type(vegn_cohort_type), pointer :: cohort

  tsoil = soil_ave_temp (soil,soil_carbon_depth_scale)
  theta = soil_ave_theta(soil,soil_carbon_depth_scale)

  call vegn_carbon_int(vegn, tsoil, theta, diag)

  ! --- accumulate values for climatological averages
  vegn%tc_av     = vegn%tc_av + cana_T
  vegn%tsoil_av  = vegn%tsoil_av + tsoil
  vegn%precip_av = vegn%precip_av + precip
  cohort => current_cohort(first_cohort(vegn%cohorts))
  vegn%theta_av  = vegn%theta_av + soil_ave_theta(soil,cohort%root_zeta)

  vegn%n_accum   = vegn%n_accum+1
  
end subroutine vegn_step_3


! ============================================================================
! update slow components of the vegetation model
subroutine update_vegn_slow( )

  ! ---- local vars ----------------------------------------------------------
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  integer :: i,j,k ! current point indices

  ! get components of calendar dates for this and previous time step
  call get_date(lnd%time,             year0,month0,day0,hour,minute,second)
  call get_date(lnd%time-lnd%dt_slow, year1,month1,day1,hour,minute,second)

  ce = first_elmt(lnd%tile_map, lnd%is, lnd%js) ; te = tail_elmt(lnd%tile_map)
  do while ( ce /= te )
     call get_elmt_indices(ce,i,j,k) ; call set_current_point(i,j,k) ! this is for debug output only
     tile => current_tile(ce) ; ce=next_elmt(ce)
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

     if (day1 /= day0) then
        call vegn_daily_npp(tile%vegn)
     endif

     ! monthly averaging
     if (month1 /= month0) then
        ! compute averages from accumulated monthly values 
        tile%vegn%tc_av     = tile%vegn%tc_av     / tile%vegn%n_accum
        tile%vegn%tsoil_av  = tile%vegn%tsoil_av  / tile%vegn%n_accum
        tile%vegn%theta_av  = tile%vegn%theta_av  / tile%vegn%n_accum
        tile%vegn%precip_av = tile%vegn%precip_av / tile%vegn%n_accum
        ! accumulate annual values
        tile%vegn%p_ann_acm = tile%vegn%p_ann_acm+tile%vegn%precip_av
        tile%vegn%t_ann_acm = tile%vegn%t_ann_acm+tile%vegn%tc_av
        if ( tile%vegn%tc_av < cold_month_threshold ) & 
             tile%vegn%ncm_acm = tile%vegn%ncm_acm+1
        tile%vegn%t_cold_acm = min(tile%vegn%t_cold_acm, tile%vegn%tc_av)

        tile%vegn%nmn_acm = tile%vegn%nmn_acm+1 ! increase the number of accumulated months
     endif

     ! annual averaging
     if (year1 /= year0) then
        if(tile%vegn%nmn_acm /= 0) then
           ! calculate annual averages from accumulated values
           tile%vegn%p_ann  = tile%vegn%p_ann_acm/tile%vegn%nmn_acm
           tile%vegn%t_ann  = tile%vegn%t_ann_acm/tile%vegn%nmn_acm
           tile%vegn%t_cold = tile%vegn%t_cold_acm
           tile%vegn%ncm    = tile%vegn%ncm_acm
           ! reset accumulated values
           tile%vegn%ncm_acm    = 0
           tile%vegn%p_ann_acm  = 0
           tile%vegn%t_ann_acm  = 0
           tile%vegn%t_cold_acm = HUGE(tile%vegn%t_cold_acm)
        endif
!!$        call calc_miami_npp(tile%vegn)
        tile%vegn%nmn_acm = 0
     endif

     if (year1 /= year0 .and. do_biogeography) then
        call vegn_biogeography(tile%vegn)
     endif

     if (month1 /= month0.and.do_patch_disturbance) then
        call update_fuel(tile%vegn)
     endif

     if (day1 /= day0 .and. do_cohort_dynamics) then
        call vegn_growth(tile%vegn)
        call vegn_nat_mortality(tile%vegn,86400.0)
     endif

     if  (month1 /= month0 .and. do_phenology) then
        call vegn_phenology (tile%vegn)
     endif

     if (year1 /= year0 .and. do_patch_disturbance) then
        call vegn_disturbance(tile%vegn, seconds_per_year)
     endif
     ! reset averages and number of steps to 0 before the start of new month
     if (month1 /= month0) then
        tile%vegn%n_accum  = 0
        tile%vegn%tc_av    = 0.
        tile%vegn%tsoil_av = 0.
        tile%vegn%theta_av = 0.
        tile%vegn%precip_av= 0.
     endif

     !reset fuel and drought months before the start of new year
     if (year1 /= year0) then
        tile%vegn%lambda     = 0
        tile%vegn%fuel       = 0
     endif

  enddo

  ! seed transport
  if (year1 /= year0 .and. do_seed_transport) then
     call vegn_seed_transport()
  endif

  ! override biomasses with static vegetation
  if(month1/=month0) &
     call static_vegn_override(lnd%time)

end subroutine update_vegn_slow


! ============================================================================
subroutine vegn_seed_transport()

  ! local vars
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  integer :: i,j ! current point indices
  real :: total_seed_supply
  real :: total_seed_demand
  real :: f_supply ! fraction of the supply that gets spent
  real :: f_demand ! fraction of the demand thet gets satisfied

  ce = first_elmt(lnd%tile_map, lnd%is, lnd%js) ; te = tail_elmt(lnd%tile_map)
  total_seed_supply = 0.0; total_seed_demand = 0.0
  do while ( ce /= te )
     call get_elmt_indices(ce,i,j)
     tile => current_tile(ce) ; ce=next_elmt(ce)
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

     total_seed_supply = total_seed_supply + vegn_seed_supply(tile%vegn)*tile%frac*lnd%garea(i,j)
     total_seed_demand = total_seed_demand + vegn_seed_demand(tile%vegn)*tile%frac*lnd%garea(i,j)
  enddo
  ! sum totals globally
  call mpp_sum(total_seed_demand)
  call mpp_sum(total_seed_supply)
  ! if either demand or supply are zeros we dont need (or can't) transport anythig
  if (total_seed_demand==0.or.total_seed_supply==0)then
     return
  end if

  ! calculate the fraction of the supply that's going to be used
  f_supply = MIN(total_seed_demand/total_seed_supply, 1.0)
  ! calculate the fraction of the demand that's going to be satisfied
  f_demand = MIN(total_seed_supply/total_seed_demand, 1.0)
  ! note that either f_supply or f_demand is 1; the mass conservation law in the
  ! following calculations is satisfied since 
  ! f_demand*total_seed_demand - f_supply*total_seed_supply == 0

  ! redistribute part (or possibly all) of the supply to satisfy part (or possibly all) 
  ! of the demand
  ce = first_elmt(lnd%tile_map) ; te = tail_elmt(lnd%tile_map)
  do while ( ce /= te )
     call get_elmt_indices(ce,i,j)
     tile => current_tile(ce) ; ce=next_elmt(ce)
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body
     
     call vegn_add_bliving(tile%vegn, &
          f_demand*vegn_seed_demand(tile%vegn)-f_supply*vegn_seed_supply(tile%vegn))
  enddo
end subroutine vegn_seed_transport


! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
end function vegn_tile_exists


! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_VEGN_ACCESSOR_0D(xtype,x) subroutine vegn_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x;endif;end subroutine

#define DEFINE_COHORT_ACCESSOR(xtype,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;end subroutine

#define DEFINE_COHORT_COMPONENT_ACCESSOR(xtype,component,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%component%x;end subroutine

DEFINE_VEGN_ACCESSOR_0D(real,age)
DEFINE_VEGN_ACCESSOR_0D(real,fast_soil_C)
DEFINE_VEGN_ACCESSOR_0D(real,slow_soil_C)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_pool)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_rate)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_pool)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_rate)
DEFINE_VEGN_ACCESSOR_0D(real,asoil_in)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_in)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_in)
DEFINE_VEGN_ACCESSOR_0D(real,tc_av)
DEFINE_VEGN_ACCESSOR_0D(real,theta_av)
DEFINE_VEGN_ACCESSOR_0D(real,tsoil_av)
DEFINE_VEGN_ACCESSOR_0D(real,precip_av)
DEFINE_VEGN_ACCESSOR_0D(real,fuel)
DEFINE_VEGN_ACCESSOR_0D(real,lambda)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold)
DEFINE_VEGN_ACCESSOR_0D(real,ncm)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold_acm)
DEFINE_VEGN_ACCESSOR_0D(real,ncm_acm)

DEFINE_COHORT_ACCESSOR(integer,species)
DEFINE_COHORT_ACCESSOR(real,bl)
DEFINE_COHORT_ACCESSOR(real,br)
DEFINE_COHORT_ACCESSOR(real,blv)
DEFINE_COHORT_ACCESSOR(real,bsw)
DEFINE_COHORT_ACCESSOR(real,bwood)
DEFINE_COHORT_ACCESSOR(real,bliving)
DEFINE_COHORT_ACCESSOR(integer,status)
DEFINE_COHORT_ACCESSOR(real,npp_previous_day)

DEFINE_COHORT_COMPONENT_ACCESSOR(real,prog,tv)
DEFINE_COHORT_COMPONENT_ACCESSOR(real,prog,wl)
DEFINE_COHORT_COMPONENT_ACCESSOR(real,prog,ws)

DEFINE_COHORT_ACCESSOR(real,height)

end module vegetation_mod
