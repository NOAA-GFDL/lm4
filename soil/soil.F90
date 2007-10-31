! ============================================================================
! soil model module
! ============================================================================
module soil_mod

use fms_mod,            only: error_mesg, file_exist,  open_namelist_file,    &
                              check_nml_error, stdlog, write_version_number, &
                              close_file, mpp_pe, mpp_root_pe, FATAL, NOTE
use mpp_io_mod,         only: mpp_open, MPP_RDONLY
use time_manager_mod,   only: time_type, increment_time, time_type_to_real
use diag_manager_mod,   only: diag_axis_init
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o, PI
use horiz_interp_mod,   only: horiz_interp

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR
use soil_tile_mod, only : &
     soil_tile_type, soil_pars_type, soil_prog_type, read_soil_data_namelist, &
     soil_data_beta, soil_data_radiation, soil_data_diffusion, soil_data_thermodynamics, &
     soil_data_hydraulics, soil_data_w_sat, max_lev

use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, get_elmt_indices, &
     operator(/=)
use land_utils_mod, only : put_to_tiles_r0d_fptr, put_to_tiles_r1d_fptr
use land_tile_diag_mod, only : diag_buff_type, &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, send_tile_data_r0d_fptr, send_tile_data_r1d_fptr, &
     send_tile_data_i0d_fptr
use land_data_mod,      only : land_state_type, lnd
use land_io_mod, only : read_field
use land_tile_io_mod, only : create_tile_out_file, write_tile_data_r1d_fptr, &
     read_tile_data_r1d_fptr, print_netcdf_error
use nf_utils_mod, only : nfu_def_dim, nfu_put_att_text
use vegn_tile_mod, only : vegn_tile_type
use land_debug_mod, only : is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public :: read_soil_namelist
public :: soil_init
public :: soil_end

public :: soil_get_sfc_temp
public :: soil_radiation
public :: soil_diffusion
public :: soil_step_1
public :: soil_step_2
! =====end of public interfaces ==============================================



! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'soil',&
    version = '',&
    tagname = ''


! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2                  = .false.
logical :: use_bucket           = .false.    ! single-layer soil water
logical :: bifurcate            = .false.    ! consider direct evap from bucket
real    :: init_temp            = 288.        ! cold-start soil T
real    :: init_w               = 150.      ! cold-start w(l)/dz(l)
real    :: init_groundwater     =   0.        ! cold-start gw storage
   real :: cpw = 1952. ! specific heat of water vapor at constant pressure
   real :: clw = 4218. ! specific heat of water (liquid)
   real :: csw = 2106. ! specific heat of water (ice)
character(len=16) :: albedo_to_use = '' ! or 'albedo-map'

namelist /soil_nml/ lm2, use_bucket,             bifurcate,             &
                    init_temp,      &
                    init_w,       &
                    init_groundwater, cpw, clw, csw, &
                    albedo_to_use
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
type(time_type) :: time
real            :: delta_time

logical         :: use_beta, use_soil_rh, use_single_geo

integer         :: num_l              ! # of water layers
real            :: dz    (max_lev)    ! thicknesses of layers
real            :: zfull (max_lev)
real            :: zhalf (max_lev+1)

! ---- diagnostic field IDs
integer :: id_lwc, id_swc, id_temp, id_ie, id_sn, id_bf, id_hie, id_hsn, &
    id_hbf, id_heat_cap, id_type, id_tau_gw, id_w_fc, &
    id_refl_dry_dif, id_refl_dry_dir, id_refl_sat_dif, id_refl_sat_dir

! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_soil_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l

  call read_soil_data_namelist(num_l,dz,use_single_geo)

  call write_version_number(version, tagname)
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=soil_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'soil_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  if (mpp_pe() == mpp_root_pe()) then
     write (stdlog(), nml=soil_nml)
  endif

  ! ---- set up vertical discretization
  zhalf(1) = 0
  do l = 1, num_l;   
     zhalf(l+1) = zhalf(l) + dz(l)
     zfull(l) = 0.5*(zhalf(l+1) + zhalf(l))
  enddo

end subroutine read_soil_namelist


! ============================================================================
! initialize soil model
subroutine soil_init ( id_lon, id_lat, id_band, use_E_max   )
  integer, intent(in)  :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in)  :: id_lat  ! ID of land latitude (Y) axis
  integer, intent(in)  :: id_band ! ID of spectral band axis
  logical, intent(out) :: use_E_max

  ! ---- local vars
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce  ! tail and current tile list elements
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  real, allocatable :: tau_gw(:,:), albedo(:,:,:) ! input data buffers for respective variables

  module_is_initialized = .TRUE.
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)

  ! ---- set flags for use in water-availability constraints -----------------
  use_beta    = .false.
  use_E_max   = .false.
  use_soil_rh = .false.
!  if (use_bucket .and. .not.bifurcate) then    !************* TEMPORARY!!!!
!      use_beta = .true.
!    else if (use_bucket .and. bifurcate) then
!      use_E_max = .true.
!    else
!      use_soil_rh = .true.
!    endif
  use_soil_rh = .true.
  use_beta    = .true.


  ! -------- initialize soil state --------
  if (file_exist('INPUT/soil.res.nc')) then
     call error_mesg('soil_init','reading NetCDF restart',NOTE)
     __NF_ASRT__(nf_open('INPUT/soil.res.nc',NF_NOWRITE,unit))
     call read_tile_data_r1d_fptr(unit, 'temp'         , soil_T_ptr  )
     call read_tile_data_r1d_fptr(unit, 'wl'           , soil_wl_ptr )
     call read_tile_data_r1d_fptr(unit, 'ws'           , soil_ws_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater'  , soil_groundwater_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater_T', soil_groundwater_T_ptr)
     __NF_ASRT__(nf_close(unit))     
  else
     te = tail_elmt (lnd%tile_map)
     ce = first_elmt(lnd%tile_map)
     do while(ce /= te)
        tile=>current_tile(ce)  ! get pointer to current tile
        ce=next_elmt(ce)       ! advance position to the next tile

        if (.not.associated(tile%soil)) cycle
        
        if (init_temp.ge.tfreeze) then      ! USE SOIL TFREEZE HERE
           tile%soil%prog(1:num_l)%wl = init_w*dz(1:num_l)
           tile%soil%prog(1:num_l)%ws = 0
        else
           tile%soil%prog(1:num_l)%wl = 0
           tile%soil%prog(1:num_l)%ws = init_w*dz(1:num_l)
        endif
        tile%soil%prog%T             = init_temp
        tile%soil%prog%groundwater   = init_groundwater
        tile%soil%prog%groundwater_T = init_temp
     enddo
  endif
  
  ! initialize soil model diagnostic fields
  call soil_diag_init ( id_lon, id_lat, id_band )

  ! read groundwater residence time field, if requested
  if (.not.use_single_geo) then
     allocate(tau_gw(lnd%is:lnd%ie,lnd%js:lnd%je))
     call read_field( 'INPUT/groundwater_residence.nc','tau',lnd%lonb, lnd%latb, tau_gw, &
          interp='bilinear' )
     call put_to_tiles_r0d_fptr( tau_gw, lnd%tile_map, soil_tau_groundwater_ptr )
     deallocate(tau_gw)
  endif

  ! set dry soil albedo values, if requested
  if (trim(albedo_to_use)=='albedo-map') then
     allocate(albedo(lnd%is:lnd%ie,lnd%js:lnd%je,NBANDS))
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_VIS',lnd%lonb, lnd%latb, albedo(:,:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_NIR',lnd%lonb, lnd%latb, albedo(:,:,BAND_NIR),'bilinear')
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_dry_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_dry_dif_ptr )
     ! for now, put the same value into the saturated soil albedo, so that
     ! the albedo doesn't depend on soil wetness
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_sat_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, lnd%tile_map, soil_refl_sat_dif_ptr )
     deallocate(albedo)
  endif
  
  ! ---- static diagnostic section
  call send_tile_data_r0d_fptr(id_tau_gw,       lnd%tile_map, soil_tau_groundwater_ptr)
  call send_tile_data_r1d_fptr(id_w_fc,         lnd%tile_map, soil_w_fc_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dir, lnd%tile_map, soil_refl_dry_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dif, lnd%tile_map, soil_refl_dry_dif_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dir, lnd%tile_map, soil_refl_sat_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dif, lnd%tile_map, soil_refl_sat_dif_ptr)
  call send_tile_data_i0d_fptr(id_type,         lnd%tile_map, soil_tag_ptr)
 
end subroutine soil_init


! ============================================================================
subroutine soil_diag_init ( id_lon, id_lat, id_band )
  integer, intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in) :: id_lat  ! ID of land latitude (Y) axis
  integer, intent(in) :: id_band ! ID of spectral band axis  

  ! ---- local vars
  integer :: axes(3)
  integer :: id_zhalf, id_zfull

  ! define vertical axis and its' edges
  id_zhalf = diag_axis_init ( &
       'zhalf_soil', zhalf(1:num_l+1), 'meters', 'z', 'half level',  -1, set_name='soil' )
  id_zfull = diag_axis_init ( &
       'zfull_soil', zfull(1:num_l),   'meters', 'z', 'full level',  -1, set_name='soil', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/ id_lon, id_lat, id_zfull /)

  ! define diagnostic fields
  id_lwc = register_tiled_diag_field ( module_name, 'soil_liq', axes,  &
       Time, 'bulk density of liquid water', 'kg/m3', missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'soil_ice',  axes,  &
       Time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'soil_T',  axes,       &
       Time, 'temperature',            'degK',  missing_value=-100.0 )
  id_ie  = register_tiled_diag_field ( module_name, 'soil_rie',  axes(1:2),  &
       Time, 'inf exc runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_sn  = register_tiled_diag_field ( module_name, 'soil_rsn',  axes(1:2),  &
       Time, 'satn runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_bf  = register_tiled_diag_field ( module_name, 'soil_rbf',  axes(1:2),  &
       Time, 'baseflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_hie  = register_tiled_diag_field ( module_name, 'soil_hie',  axes(1:2), &
       Time, 'heat ie runf',            'W/m2',  missing_value=-100.0 )
  id_hsn  = register_tiled_diag_field ( module_name, 'soil_hsn',  axes(1:2), &
       Time, 'heat sn runf',            'W/m2',  missing_value=-100.0 )
  id_hbf  = register_tiled_diag_field ( module_name, 'soil_hbf',  axes(1:2), &
       Time, 'heat bf runf',            'W/m2',  missing_value=-100.0 )

  id_heat_cap = register_tiled_diag_field ( module_name, 'soil_heat_cap',  &
       axes, Time, 'heat capacity of dry soil','J/(kg K)', missing_value=-100.0 )
  
  id_type = register_tiled_static_field ( module_name, 'soil_type',  &
       axes(1:2), 'soil type', missing_value=-1.0 )
  id_tau_gw = register_tiled_static_field ( module_name, 'tau_gw',  &
       axes(1:2), 'groundwater residence time', 's', missing_value=-100.0 )
  id_w_fc = register_tiled_static_field ( module_name, 'w_fc',  &
       axes, 'soil field capacity', missing_value=-1.0 )
  id_refl_dry_dir = register_tiled_static_field ( module_name, 'refl_dry_dir',  &
       (/id_lon, id_lat, id_band/), 'reflectance of dry soil for direct light', missing_value=-1.0 )
  id_refl_dry_dif = register_tiled_static_field ( module_name, 'refl_dry_dif',  &
       (/id_lon, id_lat, id_band/), 'reflectance of dry soil for diffuse light', missing_value=-1.0 )
  id_refl_sat_dir = register_tiled_static_field ( module_name, 'refl_sat_dir',  &
       (/id_lon, id_lat, id_band/), 'reflectance of saturated soil for direct light', missing_value=-1.0 )
  id_refl_sat_dif = register_tiled_static_field ( module_name, 'refl_sat_dif',  &
       (/id_lon, id_lat, id_band/), 'reflectance of saturated soil for diffuse light', missing_value=-1.0 )

end subroutine soil_diag_init


! ============================================================================
subroutine soil_end (tile_dim_length)
  integer, intent(in) :: tile_dim_length ! length of tile dimension in the 
                                         ! output file

  integer :: unit
  logical :: restart_created ! flag indicating that the restart file was created

  module_is_initialized =.FALSE.

  ! ---- write restart file --------------------------------------------------

  ! create output file, including internal structure necessary for output
  call error_mesg('soil_end','writing NetCDF restart',NOTE)
  call create_tile_out_file(unit,'RESTART/soil.res.nc', &
          lnd%glon*180.0/PI, lnd%glat*180/PI, soil_tile_exists, tile_dim_length,&
          created=restart_created)

  if (restart_created) then
     ! in addition, define vertical coordinate
     __NF_ASRT__(nfu_def_dim(unit,'zfull',zfull(1:num_l),'full level','m'))
     __NF_ASRT__(nfu_put_att_text(unit,'zfull','positive','down'))
        
     ! write out fields
     call write_tile_data_r1d_fptr(unit,'temp'         ,soil_T_ptr   ,'zfull','soil temperature','degrees_K')
     call write_tile_data_r1d_fptr(unit,'wl'           ,soil_wl_ptr  ,'zfull','liquid water content','kg/m2')
     call write_tile_data_r1d_fptr(unit,'ws'           ,soil_ws_ptr  ,'zfull','solid water content','kg/m2')
     call write_tile_data_r1d_fptr(unit,'groundwater'  ,soil_groundwater_ptr  ,'zfull')
     call write_tile_data_r1d_fptr(unit,'groundwater_T',soil_groundwater_T_ptr ,'zfull')
   
     ! close file
     __NF_ASRT__(nf_close(unit))
  endif

end subroutine soil_end


! ============================================================================
subroutine soil_get_sfc_temp ( soil, soil_T )
  type(soil_tile_type), intent(in) :: soil
  real, intent(out) :: soil_T

  soil_T= soil%prog(1)%T
end subroutine soil_get_sfc_temp


! ============================================================================
! compute soil radiative properties
subroutine soil_radiation ( soil,&
     soil_refl_dir, soil_refl_dif, soil_refl_lw, soil_emis )
  type(soil_tile_type), intent(in) :: soil
  real, intent(out) :: soil_refl_dir(NBANDS), soil_refl_dif(NBANDS), soil_refl_lw, soil_emis

  call soil_data_radiation ( soil, soil_refl_dir, soil_refl_dif, soil_emis )
  soil_refl_lw = 1 - soil_emis
  if(any(soil_refl_dif<0).or.any(soil_refl_dif>1).or.&
     any(soil_refl_dir<0).or.any(soil_refl_dir>1)) then
    write(*,*)'soil_refl is out of range'
    write(*,*)'soil_refl_dif=',soil_refl_dif
    write(*,*)'soil_refl_dir=',soil_refl_dir
  endif
end subroutine soil_radiation


! ============================================================================
! compute soil roughness
subroutine soil_diffusion ( soil, soil_z0s, soil_z0m )
  type(soil_tile_type), intent(in) :: soil
  real, intent(out) :: soil_z0s, soil_z0m

  call soil_data_diffusion ( soil, soil_z0s, soil_z0m )
end subroutine soil_diffusion


! ============================================================================
! update soil properties explicitly for time step.
! MAY WISH TO INTRODUCE 'UNIVERSAL' SENSITIVITIES FOR SIMPLICITY.
! T-DEPENDENCE OF HYDRAULIC PROPERTIES COULD BE DONE LESS FREQUENTLY.
! integrate soil-heat conduction equation upward from bottom of soil
! to surface, delivering linearization of surface ground heat flux.
subroutine soil_step_1 ( soil, vegn, &
                         soil_T, soil_uptake_T, soil_beta, soil_water_supply, soil_E_max, &
                         soil_rh, soil_liq, soil_ice, soil_subl, soil_tf, &
                         soil_G0, soil_DGDT )
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  real, intent(out) :: &
       soil_T, &    ! temperature of the upper layer of the soil, degK
       soil_uptake_T, & ! effective temperature for vegetation water uptake, degK
       soil_beta, &
       soil_water_supply, & ! supply of water to vegetation per unit total active root biomass, kg/m2 
       soil_E_max, &
       soil_rh,   &
       soil_liq,  & ! amount of liquid water available for evaporation 
       soil_ice,  & ! amount of ice available for sublimation
       soil_subl, & ! part of sublimation in water vapor flux, dimensionless [0,1]
       soil_tf,   & ! soil freezing temperature, degK
       soil_G0, soil_DGDT ! linearization of ground heat flux
  ! ---- local vars
  real :: bbb, denom, dt_e, vlc_sfc, vsc_sfc
  real, dimension(num_l) :: aaa, ccc, thermal_cond, heat_capacity
  integer :: l

  if(is_watch_point()) then
     write(*,*) 'soil%tag', soil%tag
     write(*,*) 'soil%pars%k_sat_ref', soil%pars%k_sat_ref 
     write(*,*) 'soil%pars%psi_sat_ref', soil%pars%psi_sat_ref
     write(*,*) 'soil%pars%chb', soil%pars%chb
     write(*,*) 'soil%pars%w_sa', soil%pars%w_sat
  endif
! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  soil_T = soil%prog(1)%T
  call soil_data_beta (soil, vegn, soil_beta, soil_water_supply)
  soil_uptake_T     = sum(soil%uptake_frac*soil%prog%T)
  if (.not.use_beta) soil_beta = 1

  vlc_sfc = max(0.0, soil%prog(1)%wl / (dens_h2o * dz(1)))
  vsc_sfc = max(0.0, soil%prog(1)%ws / (dens_h2o * dz(1)))
  call soil_data_thermodynamics ( soil, vlc_sfc, vsc_sfc,  &  
                                  soil_E_max, soil_rh, &
                                  thermal_cond )
  do l = 1, num_l
     heat_capacity(l) = soil%heat_capacity_dry(l) *dz(l) &
          + clw*soil%prog(l)%wl + csw*soil%prog(l)%ws
  enddo
!!! CHECK (bucket depth needed?) ***  soil_E_max = soil_E_max/dz(1)
  if (.not.use_soil_rh) soil_rh=1

  soil_liq  = max(soil%prog(1)%wl, 0.)
  soil_ice  = max(soil%prog(1)%ws, 0.)
  if (soil_liq + soil_ice > 0) then
     soil_subl = soil_ice / (soil_liq + soil_ice)
  else
     soil_subl = 0
  endif

  if(num_l > 1) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz(l+1)/thermal_cond(l+1) &
                     + dz(l)/thermal_cond(l)   )
        aaa(l+1) = - dt_e * delta_time / heat_capacity(l+1)
        ccc(l)   = - dt_e * delta_time / heat_capacity(l)
     enddo

     bbb = 1.0 - aaa(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(soil%prog(num_l)%T - soil%prog(num_l-1)%T)
     soil%e(num_l-1) = -aaa(num_l)/denom
     soil%f(num_l-1) = dt_e/denom
     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*soil%e(l)
        dt_e = - ( ccc(l)*(soil%prog(l+1)%T - soil%prog(l)%T  ) &
                  -aaa(l)*(soil%prog(l)%T   - soil%prog(l-1)%T) )
        soil%e(l-1) = -aaa(l)/denom
        soil%f(l-1) = (dt_e - ccc(l)*soil%f(l))/denom
     end do
     denom = delta_time/(heat_capacity(1) )
     soil_G0   = ccc(1)*(soil%prog(2)%T- soil%prog(1)%T + soil%f(1)) / denom
     soil_DGDT = (1 - ccc(1)*(1-soil%e(1))) / denom   
  else  ! one-level case
     denom = delta_time/heat_capacity(1)
     soil_G0    = 0.
     soil_DGDT  = 1. / denom
  end if
  
  ! set soil freezing temperature
  soil_tf = soil%pars%tfreeze

  if(is_watch_point()) then
     write(*,*) '#### soil_step_1 checkpoint 1 ####'
     write(*,*) 'mask    ', .true.
     write(*,*) 'T       ', soil_T
     write(*,*) 'uptake_T', soil_uptake_T
     write(*,*) 'beta    ', soil_beta
     write(*,*) 'E_max   ', soil_E_max
     write(*,*) 'rh      ', soil_rh
     write(*,*) 'liq     ', soil_liq
     write(*,*) 'ice     ', soil_ice
     write(*,*) 'subl    ', soil_subl
     write(*,*) 'G0      ', soil_G0
     write(*,*) 'DGDT    ', soil_DGDT
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g))') 'level=',l,&
             'uptake_frac=', soil%uptake_frac(l),&
             'T=', soil%prog(l)%T
     enddo
  endif
end subroutine soil_step_1


! ============================================================================
! apply boundary flows to soil water and move soil water vertically.
  subroutine soil_step_2 ( soil, diag, soil_subl, snow_lprec, snow_hlprec,  &
                           vegn_uptk, &
                           subs_DT, subs_M_imp, subs_evap, &
                           soil_levap, soil_fevap, soil_hadvec, soil_melt, &
                           soil_lrunf, soil_hlrunf, soil_Ttop, soil_Ctop, &
                           soil_LMASS, soil_FMASS, soil_HEAT )
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: &
     soil_subl     !
  real, intent(in) :: &
     snow_lprec, &
     snow_hlprec, &
     vegn_uptk, &
     subs_DT,       &!
     subs_M_imp,       &! rate of phase change of non-evaporated soil water
     subs_evap
  real, intent(out) :: &
     soil_levap, soil_fevap, soil_hadvec, soil_melt, &
     soil_lrunf, soil_hlrunf, soil_Ttop, soil_Ctop, &
     soil_LMASS, soil_FMASS, soil_HEAT

  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l) :: del_t, eee, fff, &
             psi, DThDP, hyd_cond, DKDP, K, DKDPm, DKDPp, grad, &
             vlc, vsc, dW_l, u_minus, u_plus, DPsi, soil_w_fc, soil_w_sat
  real, dimension(num_l+1) :: flow, infilt
  real, dimension(num_l  ) :: div
  real      :: &
     lprec_eff, hlprec_eff, tflow, hcap,cap_flow, &
     melt_per_deg, melt, adj,&
     lrunf_sn,lrunf_ie,lrunf_bf, hlrunf_sn,hlrunf_ie,hlrunf_bf, &
     Qout, DQoutDP,&
     tau_gw, c0, c1, c2, x, aaa, bbb, ccc, ddd, xxx, Dpsi_min, Dpsi_max, &
     liq_frac, excess_wat, excess_liq, excess_ice, h1, h2, summax, &
     space_avail, liq_placed, ice_placed, excess_t
  logical :: stiff
  real, dimension(num_l-1) :: del_z
  integer :: l
  real :: jj
  
  jj = 1.
  
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 1 #####'
     write(*,*) 'mask    ', .true.
     write(*,*) 'subs_evap    ', subs_evap
     write(*,*) 'snow_lprec   ', snow_lprec
     write(*,*) 'uptake  ', vegn_uptk
     write(*,*) 'subs_M_imp   ', subs_M_imp
     write(*,*) 'theta_s ', soil%pars%w_sat
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g))') 'level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws+soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws,&
             ' gw=', soil%prog(l)%groundwater
     enddo
  endif

  ! ---- record fluxes ---------
  soil_levap  = subs_evap*(1-soil_subl)
  soil_fevap  = subs_evap*   soil_subl
  soil_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution --------------
  del_t(1) = subs_DT
  soil%prog(1)%T = soil%prog(1)%T + del_t(1)
  if ( num_l > 1) then
    do l = 1, num_l-1
      del_t(l+1) = soil%e(l) * del_t(l) + soil%f(l)
      soil%prog(l+1)%T = soil%prog(l+1)%T + del_t(l+1)
    end do
  end if

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 2 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g))') 'level=',l, 'T=', soil%prog(l)%T, &
             'del_t=', del_t(l), 'e=', soil%e(l), 'f=', soil%f(l)
     enddo
  endif

  soil_hadvec = 0.
  ! ---- extract evap from soil and do implicit melt --------------------
  IF(LM2) THEN
    do l = 1, num_l
      soil%prog(l)%wl = soil%prog(l)%wl &
                      - soil%uptake_frac(l)*soil_levap*delta_time
      soil_hadvec = soil_hadvec - soil%uptake_frac(l)*soil_levap &
                                      *clw*(soil%prog(l)%T-tfreeze)
    enddo
!      soil%prog(1)%wl = soil%prog(1)%wl - soil_levap*delta_time
!      soil%prog(1)%ws = soil%prog(1)%ws - soil_fevap*delta_time
!      soil_hadvec = -cpw*subs_evap*(soil%prog(1)%T-tfreeze)
    hcap = soil%heat_capacity_dry(1)*dz(1) &
                    + clw*soil%prog(1)%wl + csw*soil%prog(1)%ws
    soil%prog(1)%T = soil%prog(1)%T + (   &
               +((clw-cpw)*soil_levap                              &
               + (csw-cpw)*soil_fevap)*(soil%prog(1)%T  -tfreeze) &
                                            )*delta_time/ hcap
    soil%prog(1)%wl = soil%prog(1)%wl + subs_M_imp
    soil%prog(1)%ws = soil%prog(1)%ws - subs_M_imp
    soil%prog(1)%T  = tfreeze + (hcap*(soil%prog(1)%T-tfreeze) ) &
                           / ( hcap + (clw-csw)*subs_M_imp )
  ELSE
    soil%prog(1)%wl = soil%prog(1)%wl - soil_levap*delta_time
    soil%prog(1)%ws = soil%prog(1)%ws - soil_fevap*delta_time
    soil_hadvec = -cpw*subs_evap*(soil%prog(1)%T-tfreeze)
    hcap = soil%heat_capacity_dry(1)*dz(1) &
                       + clw*soil%prog(1)%wl + csw*soil%prog(1)%ws
    soil%prog(1)%T = soil%prog(1)%T + (   &
                  +((clw-cpw)*soil_levap                              &
                  + (csw-cpw)*soil_fevap)*(soil%prog(1)%T  -tfreeze) &
                                               )*delta_time/ hcap
    soil%prog(1)%wl = soil%prog(1)%wl + subs_M_imp
    soil%prog(1)%ws = soil%prog(1)%ws - subs_M_imp
    soil%prog(1)%T  = tfreeze + (hcap*(soil%prog(1)%T-tfreeze) ) &
                              / ( hcap + (clw-csw)*subs_M_imp )
  ENDIF

  do l = 1, num_l
    soil%prog(l)%wl = soil%prog(l)%wl &
                      - soil%uptake_frac(l)*vegn_uptk*delta_time
    soil_hadvec = soil_hadvec - soil%uptake_frac(l)*vegn_uptk &
                                      *clw*(soil%prog(l)%T-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g))') ' level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws+soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws
     enddo
  endif
  hlrunf_ie=0;hlrunf_sn=0;hlrunf_bf=0;lrunf_ie=0;lrunf_sn = 0;lrunf_bf = 0
  ! ---- push down any excess surface water, with heat ---------------------
  call soil_data_w_sat(soil, soil_w_sat)
  liq_frac=0;excess_wat=0;excess_liq=0;excess_ice=0;h1=0;h2=0;liq_frac=0
  l = 1
  summax = max(0.,soil%prog(l)%wl)+max(0.,soil%prog(l)%ws)
  if (summax > 0) then
     liq_frac = max(0.,soil%prog(l)%wl) / summax
  else
     liq_frac = 1
  endif
  excess_wat = max(0., soil%prog(l)%wl + soil%prog(l)%ws &
       - dens_h2o*dz(l)*soil_w_sat(l) )
  excess_liq = excess_wat*liq_frac
  excess_ice = excess_wat-excess_liq
  excess_t   = soil%prog(l)%T
  soil%prog(l)%wl = soil%prog(l)%wl - excess_liq
  soil%prog(l)%ws = soil%prog(l)%ws - excess_ice

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.001 #####'
     write(*,*) ' level=', l,&
          ' summax =', summax,&
          ' liq_frac =', liq_frac,&
          ' soil_w_sat =', soil_w_sat(l),&
          ' excess_liq =', excess_liq,&
          ' excess_ice =', excess_ice, &
          ' dens_h2o=', dens_h2o, &
          ' dz(l)=',dz(l),&
          'friday am'
  endif

  do l = 2, num_l-1
     if (excess_liq+excess_ice>0) then
        space_avail = dens_h2o*dz(l)*soil_w_sat(l) &
             - soil%prog(l)%wl + soil%prog(l)%ws
        liq_placed = max(min(space_avail, excess_liq), 0.)
        ice_placed = max(min(space_avail-liq_placed, excess_ice), 0.)
        h1 = (soil%heat_capacity_dry(l)*dz(l) &
             + csw*soil%prog(l)%ws + clw*soil%prog(l)%wl)
        h2 = liq_placed*clw+ice_placed*csw
        soil%prog(l)%T = (h1 * soil%prog(l)%T &
             + h2 * excess_T )  / (h1+h2)
        soil%prog(l)%wl = soil%prog(l)%wl + liq_placed
        soil%prog(l)%ws = soil%prog(l)%ws + ice_placed
        excess_liq = excess_liq - liq_placed
        excess_ice = excess_ice - ice_placed
     endif
  enddo

  l = num_l
  if (excess_liq+excess_ice > 0) then
     liq_placed = excess_liq
     ice_placed = excess_ice
     h1 = (soil%heat_capacity_dry(l)*dz(l) &
          + csw*soil%prog(l)%ws + clw*soil%prog(l)%wl)
     h2 = liq_placed*clw+ice_placed*csw
     soil%prog(l)%T = (h1 * soil%prog(l)%T &
          + h2 * excess_T )  / (h1+h2)
     soil%prog(l)%wl = soil%prog(l)%wl + liq_placed
     soil%prog(l)%ws = soil%prog(l)%ws + ice_placed
  endif

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.01 #####'
     do l = 1, num_l
        write(*,'(x,a,x,i2.2,100(x,a,g))') ' level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws +soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws
     enddo
  endif

  ! ---- fetch soil hydraulic properties -------------------------------------
  vlc=0;vsc=0
  do l = 1, num_l
    vlc(l) = max(0., soil%prog(l)%wl / (dens_h2o*dz(l)))
    vsc(l) = max(0., soil%prog(l)%ws / (dens_h2o*dz(l)))
  enddo
  call soil_data_hydraulics (soil, vlc, vsc, &
                   psi, DThDP, hyd_cond, DKDP, Dpsi_min, Dpsi_max, tau_gw, &
                   soil_w_fc )

  IF (lm2) THEN ! ********************************

     if(is_watch_point()) then
        write(*,*) ' ##### soil_step_2 checkpoint 3.1 #####'
        do l = 1, num_l
           write(*,'(x,a,x,i2.2,100(x,a,g))')'level=', l, 'vlc', vlc(l), 'K  ', hyd_cond(l)
        enddo
     endif
  ! ---- remainder of mass fluxes and associated sensible heat fluxes --------
    flow=1
    flow(1)  = 0
    do l = 1, num_l
      infilt(l) = soil%uptake_frac(l)*snow_lprec *delta_time
      flow(l+1) = max(0., soil%prog(l)%wl + flow(l) &
            + infilt(l) - soil_w_fc(l)*dz(l)*dens_h2o)
      dW_l(l) = flow(l) - flow(l+1) + infilt(l)
      soil%prog(l)%wl = soil%prog(l)%wl + dW_l(l)
    enddo
    do l = 1, num_l
      flow(l) = flow(l) + infilt(l)
    enddo
    dW_l=0
    lrunf_bf = lrunf_bf + flow(num_l)/delta_time
  ELSE   ! ********************************
    div = 0.; flow=0
    do l = 1, num_l
      if (vsc(l).eq.0. .and. psi(l).gt.0.) then
        div(l) = 0.15*dens_h2o*dz(l)/tau_gw
      endif
    enddo
    lrunf_bf = lrunf_bf + sum(div)

    if(is_watch_point()) then
       do l = 1, num_l
          write(*,'(a,1x,i2.2,100(2x,g))')'div,vsc,psi,dz',l,div(l),vsc(l),psi(l),dz(l)
       enddo
       write(*,*)'lrunf_bf',lrunf_bf
       write(*,*)'tau_gw',tau_gw
       write(*,*)'dens_h2o',dens_h2o
    endif
  ! ---- soil-water flow ----------------------------------------------------
    stiff = all(DThDP.eq.0)
    if (snow_lprec.ne.0. .and. psi(num_l).gt.0.) then
      lrunf_sn = snow_lprec*min((psi(num_l)/zhalf(num_l))**soil%pars%rsa_exp,1.)
      hlrunf_sn = lrunf_sn*snow_hlprec/snow_lprec
    else
      lrunf_sn = 0.
      hlrunf_sn = 0.
    endif
    lprec_eff = snow_lprec - lrunf_sn
    hlprec_eff = snow_hlprec - hlrunf_sn
    flow(1) = delta_time*lprec_eff
    do l = 1, num_l-1
      del_z(l) = zfull(l+1)-zfull(l)
      K(l) = 0.5*(hyd_cond(l)+hyd_cond(l+1))
      DKDPm(l) = 0. !0.5*DKDP(l)
      DKDPp(l) = 0. ! 0.5*DKDP(l+1)
!        K(l) = hyd_cond(l)
!        DKDPm(l) = DKDP(l)
!        DKDPp(l) = 0
      grad(l)  = jj*(psi(l+1)-psi(l))/del_z(l) - 1
    enddo

    if(is_watch_point()) then
       write(*,*) '##### soil_step_2 checkpoint 3.1 #####'
       do l = 1, num_l
          write(*,'(x,a,x,i2.2,x,a,100(x,g))') 'level=', l, 'DThDP,hyd_cond,psi,DKDP', &
               DThDP(l),&
               hyd_cond(l),&
               psi(l),&
               DKDP(l)
       enddo
       do l = 1, num_l-1
          write(*,'(a,i2.2,1x,a,100(2x,g))') 'interface=', l, 'K,DKDPm,DKDPp,grad,del_z', &
               K(l),&
               DKDPm(l),&
               DKDPp(l),&
               grad(l)
       enddo
    endif


    l = num_l
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    aaa =     - ( jj* K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
!      where (stiff)
    bbb = xxx - (- jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
    ddd = - K(l-1) *grad(l-1) - div(l)
!        elsewhere
!          Qout = hyd_cond(l) ! gravity drainage
!          DQoutDP = DKDP(l)  ! gravity drainage
!          Qout = 0.                ! no drainage
!          DQoutDP = 0.             ! no drainage
!          where (psi(l).gt.0.) ! linear baseflow from gw
!              Qout = 0.15*psi(l)/tau_gw
!              DQoutDP = 0.15/tau_gw
!            elsewhere
!              Qout = 0.
!              DQoutDP = 0.
!            endwhere
!          bbb = xxx - (- jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
!                      -DQoutDP )
!          ddd = -Qout - K(l-1) *grad(l-1)
!        endwhere
    eee(l-1) = -aaa/bbb
    fff(l-1) =  ddd/bbb
  
    if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g))') 'l,a,b, ,d', l,aaa, bbb,ddd
    endif

    do l = num_l-1, 2, -1
      xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
      aaa = - ( jj*K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
      bbb = xxx-( -jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
                  -jj*K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
      ccc =   - (  jj*K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
      ddd =       K(l)*grad(l) - K(l-1)*grad(l-1) &
                            - div(l)
      eee(l-1) =                    -aaa/(bbb+ccc*eee(l))
      fff(l-1) =  (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
      if(is_watch_point()) then
         write(*,'(a,i2.2,100(2x,g))') 'l,a,b,c,d', l,aaa, bbb,ccc,ddd
      endif
    enddo
  
    lrunf_ie=0
    l = 1
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    bbb = xxx - ( -jj*K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
    ccc =     - (  jj*K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
    ddd =          flow(1)/delta_time +    K(l)     *grad(l) &
                            - div(l)
    if (stiff) then
      dPsi(l) =  - psi(l)
    else
      dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
      dPsi(l) = min (dPsi(l), Dpsi_max)
      dPsi(l) = max (dPsi(l), Dpsi_min)
    endif
    flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                      - K(l)*grad(l))*delta_time
    lrunf_ie         = lprec_eff - flow(l)/delta_time
  
    if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g))') 'l,  b,c,d', l, bbb,ccc,ddd
       write(*,*) ' ##### soil_step_2 checkpoint 3.2 #####'
       write(*,*) 'ie,sn,bf:', lrunf_ie,lrunf_sn,lrunf_bf
       do l = 1, num_l-1
          write(*,'(a,i2.2,100(2x,g))') 'l,eee(l),fff(l)',l,eee(l),fff(l)
       enddo
       write(*,*) 'DThDP(1)', DThDP(1)
       write(*,*) 'K(1)', K(1)
       write(*,*) 'grad(1)', grad(1)
       write(*,*) 'ddd(1)', ddd
       write(*,*) 'ccc(1)', ccc
       write(*,*) 'bbb(1)', bbb
       write(*,*) 'dPsi(1)', dPsi(1)
       write(*,*) 'Psi(1)', Psi(1)
       write(*,*) 'div(1)', div(1)
    endif

    do l = 2, num_l
      dPsi(l) = eee(l-1)*dPsi(l-1) + fff(l-1)
    enddo
  
    do l = 1, num_l-1
      flow(l+1) = delta_time*( &
           -K(l)*(grad(l)&
           +jj*(DPsi(l+1)-DPsi(l))/ del_z(l)) &
           -grad(l)*(DKDPp(l)*Dpsi(l+1)+ &
                           DKDPm(l)*Dpsi(l) )  )
      dW_l(l) = flow(l) - flow(l+1) - div(l)*delta_time
!      soil%prog(l)%wl = soil%prog(l)%wl + dW_l(l)
    enddo
!  where (stiff)
    flow(num_l+1) = 0.
!    elsewhere
!      flow(num_l+1) = (Qout &
!                         + DQoutDP*DPsi(num_l)) * delta_time
!    endwhere
    dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                            - div(num_l)*delta_time
!    soil%prog(num_l)%wl = soil%prog(num_l)%wl + dW_l(num_l)

    if(is_watch_point()) then
       write(*,*) ' ##### soil_step_2 checkpoint 3.21 #####'
       do l = 1, num_l
          write(*,'(i2.2,100(2x,a,g))') l,&
               ' dW_l=', dW_l(l),&
               ' flow=', flow(l),&
               ' div=', div(l)
       enddo
    endif

! Squeegee back out of soil any excess liquid saturation in those rare situations
! where lrunf_ie is negative:
    if (lrunf_ie < -1.0e-4) then
       do l = num_l, 1, -1
          adj = max(dW_l(l)+soil%prog(l)%ws+soil%prog(l)%wl &
               - soil_w_sat(l)*dz(l)*dens_h2o, 0. )

          if(is_watch_point()) then
             write(*,*) '3.22 l=', l,&
                  ' soil_prog%wl=',soil%prog(l)%wl,  &
                  ' soil_prog%ws=',soil%prog(l)%ws , &
                  ' soil_w_sat=', soil_w_sat(l), &
                  ' dz=', dz(l), &
                  ' adj=', adj
          endif

          adj = min(adj, max(0.,soil%prog(l)%wl))

          if(is_watch_point()) then
             write(*,*) '3.23 l=', l, ' adj=', adj
          endif

          dW_l(l) = dW_l(l) - adj
          flow(l) = flow(l+1) + dW_l(l) + div(l)*delta_time
       enddo
       lrunf_ie = lprec_eff - flow(1)/delta_time
    endif
!
    do l = 1, num_l
       soil%prog(l)%wl = soil%prog(l)%wl + dW_l(l)
    enddo

  ENDIF ! ************************************

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',soil%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g))') l, &
             'Th=', (soil%prog(l)%ws +soil%prog(l)%wl)/(dens_h2o*dz(l)), &
             'wl=', soil%prog(l)%wl, &
             'ws=', soil%prog(l)%ws, &
             'dW_l=', dW_l(l), &
             'dPsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif


  soil_hadvec = soil_hadvec + snow_hlprec
  if  (snow_lprec.ne.0.) then
    tflow = tfreeze + snow_hlprec/(clw*snow_lprec)
  else
    tflow = tfreeze
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.4 ***** '
     write(*,*) ' tfreeze', tfreeze
     write(*,*) ' snow_hlprec', snow_hlprec
  endif


! For initial testing, use top-down-flow weights to advect heat.
  u_minus = 1.
  u_plus  = 0.
  if (flow(1).lt.0.) u_minus(1) = 0.
  hcap = (soil%heat_capacity_dry(num_l)*dz(num_l) &
                              + csw*soil%prog(num_l)%ws)/clw
  aaa = -flow(num_l) * u_minus(num_l)
  bbb =  hcap + soil%prog(num_l)%wl - dW_l(num_l) - aaa
  eee(num_l-1) = -aaa/bbb
  fff(num_l-1) = aaa*(soil%prog(num_l)%T-soil%prog(num_l-1)%T) / bbb

  do l = num_l-1, 2, -1
    hcap = (soil%heat_capacity_dry(l)*dz(l) &
                              + csw*soil%prog(l)%ws)/clw
    aaa = -flow(l)   * u_minus(l)
    ccc =  flow(l+1) * u_plus (l)
    bbb =  hcap + soil%prog(l)%wl - dW_l(l) - aaa - ccc
    eee(l-1) = -aaa / ( bbb +ccc*eee(l) )
    fff(l-1) = (   aaa*(soil%prog(l)%T-soil%prog(l-1)%T)    &
                       + ccc*(soil%prog(l)%T-soil%prog(l+1)%T)    &
                       - ccc*fff(l) ) / ( bbb +ccc*eee(l) )
  enddo
    
  hcap = (soil%heat_capacity_dry(1)*dz(1) + csw*soil%prog(1)%ws)/clw
  aaa = -flow(1) * u_minus(1)
  ccc =  flow(2) * u_plus (1)
  bbb =  hcap + soil%prog(1)%wl - dW_l(1) - aaa - ccc

  del_t(1) =  (  aaa*(soil%prog(1)%T-tflow          ) &
                     + ccc*(soil%prog(1)%T-soil%prog(2)%T) &
                     - ccc*fff(1) ) / (bbb+ccc*eee(1))
  soil%prog(1)%T = soil%prog(1)%T + del_t(1)

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.4.1 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'aaa', aaa
     write(*,*) 'bbb', bbb
     write(*,*) 'ccc', ccc
     write(*,*) 'del_t(1)', del_t(1)
     write(*,*) ' T(1)', soil%prog(1)%T
  endif

  do l = 1, num_l-1
    del_t(l+1) = eee(l)*del_t(l) + fff(l)
    soil%prog(l+1)%T = soil%prog(l+1)%T + del_t(l+1)
  enddo

  tflow = soil%prog(num_l)%T

!  do l = 1, num_l
!    where (mask)
!        hcap = soil%heat_capacity_dry(l)*dz(l) &
!                 + clw*(soil%prog(l)%wl-dW_l(l)) + csw*soil%prog(l)%ws
!        cap_flow = clw*flow(l)
!        soil%prog(l)%T = (hcap*soil%prog(l)%T + cap_flow*tflow) &
!                         /(hcap                 + cap_flow      )
!        tflow  = soil%prog(l)%T
!      endwhere
!    enddo

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.5 ***** '
     write(*,*) 'hcap', hcap
!     write(*,*) 'cap_flow', cap_flow
     do l = 1, num_l
        write(*,*) 'level=', l, ' T', soil%prog(l)%T
     enddo
  endif

  ! ---- groundwater ---------------------------------------------------------
  ! THIS T AVERAGING IS WRONG, BECAUSE IT NEGLECTS THE MEDIUM  ***
  ! ALSO, FREEZE-THAW IS NEEDED!
  ! PROBABLY THIS SECTION WILL BE DELETED ANYWAY, WITH GW TREATED ABOVE.
  IF (lm2) THEN
    do l = 1, 1      !TEMPORARY LAYER THING !!!!***
      if (soil%prog(l)%groundwater + flow(num_l+1) .ne. 0.) then ! TEMP FIX
          soil%prog(l)%groundwater_T =    &
           + (soil%prog(l)%groundwater*soil%prog(l)%groundwater_T &
              + flow(num_l+1)*tflow) &
            /(soil%prog(l)%groundwater + flow(num_l+1))
      endif
      c0 = delta_time/tau_gw
      c1 = exp(-c0)
      c2 = (1-c1)/c0
      x  = (1-c1)*soil%prog(l)%groundwater/delta_time &
                          + (1-c2)*flow(num_l+1)/delta_time
      soil%prog(l)%groundwater = c1 * soil%prog(l)%groundwater &
                                + c2 * flow(num_l+1)
      soil_lrunf  = x
      soil_hlrunf = x*clw*(soil%prog(l)%groundwater_T-tfreeze)
      soil_hadvec = soil_hadvec - soil_hlrunf
    enddo
  ELSE
    if (lprec_eff.ne.0. .and. flow(1).ge.0. ) then
      hlrunf_ie = lrunf_ie*hlprec_eff/lprec_eff
    else if (flow(1).lt.0. ) then
      hlrunf_ie = hlprec_eff - (flow(1)/delta_time)*clw &
                         *(soil%prog(1)%T-tfreeze)
    else
      hlrunf_ie = 0.
    endif
    hlrunf_bf = hlrunf_bf + clw*sum(div*(soil%prog%T-tfreeze))
    soil_lrunf  = lrunf_sn + lrunf_ie + lrunf_bf
    soil_hlrunf = hlrunf_sn + hlrunf_bf + hlrunf_ie 
    soil_hadvec = soil_hadvec - soil_hlrunf
  ENDIF

  do l = 1, num_l
    ! ---- compute explicit melt/freeze --------------------------------------
    hcap = soil%heat_capacity_dry(l)*dz(l) &
             + clw*soil%prog(l)%wl + csw*soil%prog(l)%ws
    melt_per_deg = hcap/hlf
    if       (soil%prog(l)%ws>0 .and. soil%prog(l)%T>soil%pars%tfreeze) then
      melt =  min(soil%prog(l)%ws, (soil%prog(l)%T-soil%pars%tfreeze)*melt_per_deg)
    else if (soil%prog(l)%wl>0 .and. soil%prog(l)%T<soil%pars%tfreeze) then
      melt = -min(soil%prog(l)%wl, (soil%pars%tfreeze-soil%prog(l)%T)*melt_per_deg)
    else
      melt = 0
    endif
    soil%prog(l)%wl = soil%prog(l)%wl + melt
    soil%prog(l)%ws = soil%prog(l)%ws - melt
    soil%prog(l)%T = tfreeze &
       + (hcap*(soil%prog(l)%T-tfreeze) - hlf*melt) &
                            / ( hcap + (clw-csw)*melt )
    soil_melt = soil_melt + melt / delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 5 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g))') ' level=', l,&
             ' T =', soil%prog(l)%T,&
             ' Th=', (soil%prog(l)%ws +soil%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', soil%prog(l)%wl,&
             ' ws=', soil%prog(l)%ws,&
             ' gw=', soil%prog(l)%groundwater
     enddo
  endif


  soil_LMASS = 0
  soil_FMASS = 0
  soil_HEAT = 0
  do l = 1, num_l
    soil_LMASS = soil_LMASS &
         +      soil%prog(l)%wl
    soil_FMASS = soil_FMASS &
         +      soil%prog(l)%ws
    soil_HEAT = soil_HEAT &
         + (soil%heat_capacity_dry(l)*dz(l) &
             + clw*soil%prog(l)%wl + csw*soil%prog(l)%ws)  &
                        * (soil%prog(l)%T-tfreeze)
  enddo
  soil_LMASS = soil_LMASS +     soil%prog(1)%groundwater
  soil_HEAT = soil_HEAT + clw*soil%prog(1)%groundwater    &
                                  *(soil%prog(1)%groundwater_T-tfreeze)
  soil_Ttop = soil%prog(1)%T
  soil_Ctop = soil%heat_capacity_dry(1)*dz(1) &
    + clw*soil%prog(1)%wl + csw*soil%prog(1)%ws


! ----------------------------------------------------------------------------
! given solution for surface energy balance, write diagnostic output.
!  

  ! ---- increment time and do diagnostics -----------------------------------
  time = increment_time(time, int(delta_time), 0)
  
  ! ---- diagnostic section
   call send_tile_data(id_temp, soil%prog(1:num_l)%T, diag)
   call send_tile_data(id_lwc,  soil%prog(1:num_l)%wl/dz(1:num_l), diag)
   call send_tile_data(id_swc,  soil%prog(1:num_l)%ws/dz(1:num_l), diag)
   call send_tile_data(id_ie,   lrunf_ie, diag)
   call send_tile_data(id_sn,   lrunf_sn, diag)
   call send_tile_data(id_bf,   lrunf_bf, diag)
   call send_tile_data(id_hie,  hlrunf_ie, diag)
   call send_tile_data(id_hsn,  hlrunf_sn, diag)
   call send_tile_data(id_hbf,  hlrunf_bf, diag)

   call send_tile_data(id_heat_cap, soil%heat_capacity_dry(1:num_l), diag)

end subroutine soil_step_2


! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function soil_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   soil_tile_exists = associated(tile%soil)
end function soil_tile_exists


! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_SOIL_COMPONENT_ACCESSOR_0D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x;endif;end subroutine
#define DEFINE_SOIL_COMPONENT_ACCESSOR_1D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p(:);p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x;endif;end subroutine

DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau_groundwater)

subroutine soil_w_fc_ptr(t,p);
  type(land_tile_type),pointer::t;
  real,pointer::p(:);
  p=>NULL();
  if(associated(t))then;
     if(associated(t%soil))p=>t%soil%w_fc;
  endif;
end subroutine

subroutine soil_tag_ptr(t,p);
  type(land_tile_type),pointer::t;
  integer,pointer::p;
  p=>NULL();
  if(associated(t))then;
     if(associated(t%soil))p=>t%soil%tag;
  endif;
end subroutine

DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dir)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dif)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dir)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dif)

DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,T)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,wl)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,ws)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,groundwater)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,prog,groundwater_T)

end module soil_mod



