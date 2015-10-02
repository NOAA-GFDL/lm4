! ============================================================================
! canopy air
! ============================================================================
module canopy_air_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : write_version_number, error_mesg, FATAL, NOTE, file_exist, &
     close_file, check_nml_error, mpp_pe, mpp_root_pe, stdlog, string
use fms_io_mod, only : read_compressed, restart_file_type, free_restart_type
use fms_io_mod, only : field_exist, get_field_size, save_restart
use fms_io_mod, only : register_restart_axis, register_restart_field, set_domain, nullify_domain
use time_manager_mod, only : time_type, time_type_to_real
use constants_mod, only : rdgas, rvgas, cp_air, PI, VONKARM
use sphum_mod, only : qscomp
use field_manager_mod, only : parse, MODEL_ATMOS, MODEL_LAND
use tracer_manager_mod, only : get_number_tracers, get_tracer_index, get_tracer_names
use tracer_manager_mod, only : get_tracer_index, query_method, NO_TRACER

use nf_utils_mod, only : nfu_inq_var
use land_constants_mod, only : NBANDS,d608,mol_CO2,mol_air
use land_tracers_mod, only : ntcana, isphum, ico2
use cana_tile_mod, only : cana_tile_type, &
     canopy_air_mass, canopy_air_mass_for_tracers, cpw
use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=)
use land_data_mod,      only : land_state_type, lnd
use land_tile_io_mod, only : create_tile_out_file, &
     read_tile_data_r0d_fptr, write_tile_data_r0d_fptr, &
     read_tile_data_r1d_fptr, write_tile_data_r1d_fptr, &
     get_input_restart_name, print_netcdf_error, &
     gather_tile_data, assemble_tiles
use land_debug_mod, only : is_watch_point, check_temp_range

implicit none
private

! ==== public interfaces =====================================================
public :: read_cana_namelist
public :: cana_init
public :: cana_end
public :: save_cana_restart
public :: save_cana_restart_new
public :: cana_radiation
public :: cana_turbulence
public :: cana_roughness
public :: cana_state
public :: cana_step_1
public :: cana_step_2
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
  version = '$Id: canopy_air.F90,v 21.0.2.7 2015/04/09 19:29:13 pjp Exp $', &
  tagname = '$Name: ulm_pjp $', &
  module_name = 'canopy_air_mod'

! options for turbulence parameter calculations
integer, parameter :: TURB_LM3W = 1, TURB_LM3V = 2
! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
real :: init_T           = 288.
real :: init_T_cold      = 260.
real :: init_q           = 0.
real :: init_co2         = 350.0e-6 ! ppmv = mol co2/mol of dry air
character(len=32) :: turbulence_to_use = 'lm3w' ! or lm3v
logical :: save_qco2     = .TRUE.
logical :: sfc_dir_albedo_bug = .FALSE. ! if true, reverts to buggy behavior
! where direct albedo was mistakenly used for part of sub-canopy diffuse light
logical :: allow_small_z0 = .FALSE. ! to use z0 provided by lake and glac modules
real :: lai_min_turb = 0.0 ! fudge to desensitize Tv to SW/cosz inconsistency
real :: bare_rah_sca     = 0.01     ! bare-ground resistance between ground and canopy air, s/m
namelist /cana_nml/ &
  init_T, init_T_cold, init_q, init_co2, turbulence_to_use, &
  canopy_air_mass, canopy_air_mass_for_tracers, cpw, save_qco2, &
  sfc_dir_albedo_bug, allow_small_z0, lai_min_turb, bare_rah_sca
!---- end of namelist --------------------------------------------------------

logical            :: module_is_initialized =.FALSE.
integer :: turbulence_option ! selected option of turbulence parameters 
     ! calculations

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)


contains


! ============================================================================
subroutine read_cana_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=cana_nml, iostat=io)
     ierr = check_nml_error(io, 'cana_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=cana_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'cana_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=cana_nml)
  endif
end subroutine read_cana_namelist

! ============================================================================
! initialize canopy air
subroutine cana_init ( id_lon, id_lat, new_land_io )
  integer, intent(in)          :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in)          :: id_lat  ! ID of land latitude (Y) axis
  logical, intent(in)          :: new_land_io  ! This is a transition var and will be removed

  ! ---- local vars ----------------------------------------------------------
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce ! last and current tile
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  character(len=256) :: restart_file_name
  character(len=17) :: restart_base_name='INPUT/cana.res.nc'
  integer :: siz(4)
  integer, allocatable :: idx(:)         ! I/O domain vector of compressed indices
  real,    allocatable :: r0d(:)         ! I/O domain vector of real data
  logical :: found,restart_exists

  character(len=32)  :: name  ! name of the tracer
  integer            :: tr, i ! tracer indices
  real               :: init_tr(ntcana) ! initial (cold-start) values of tracers
  real               :: value ! used for parameter parsing
  character(len=32)  :: scheme
  character(len=128) :: parameters

  module_is_initialized = .TRUE.

  ! ---- initialize cana state -----------------------------------------------
  ! get the initial conditions for tracers

  ! For  now, we get the initial value from the surface_value parameter of the
  ! *atmospheric* tracer vertical profile, except for co2 and sphum. init_q in 
  ! the cana_nml sets the initial value of the specific humidity, and init_co2
  ! sets the initial value of the dry volumetric mixing ration for co2.
  ! If surface_value is not defined in tracer table, then initial condition is zero
  init_tr(:) = 0.0
  do i = 1, ntcana
     call get_tracer_names(MODEL_LAND, i, name=name)
     tr = get_tracer_index(MODEL_ATMOS, name)
     if (tr==NO_TRACER) cycle ! no such tracer in the atmos
     ! TODO: possibly we need to add an ability to read init value from some parameter
     ! in the land tracer table?
     scheme = ''; parameters = ''
     if (query_method('profile_type', MODEL_ATMOS, tr, scheme, parameters)) then
        if (parse(parameters,'surface_value',value)>0) init_tr(i) = value
     endif
  enddo
  init_tr(isphum) = init_q
  init_tr(ico2)   = init_co2*mol_CO2/mol_air*(1-init_tr(isphum)) ! convert to kg CO2/kg wet air

  ! first, set the initial values
  te = tail_elmt (lnd%tile_map)
  ce = first_elmt(lnd%tile_map)
  do while(ce /= te)
     tile=>current_tile(ce)  ! get pointer to current tile
     ce=next_elmt(ce)       ! advance position to the next tile
     
     if (.not.associated(tile%cana)) cycle
     
     if (associated(tile%glac)) then
        tile%cana%T = init_T_cold
     else
        tile%cana%T = init_T
     endif
     tile%cana%tr(:) = init_tr(:)
  enddo

  ! then read the restart if it exists
  call get_input_restart_name(restart_base_name,restart_exists,restart_file_name)
  if (restart_exists) then
     call error_mesg('cana_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     if(new_land_io)then
        call error_mesg('cana_init', 'Using new canopy restart read', NOTE)

        call get_field_size(restart_base_name,'tile_index',siz, field_found=found, domain=lnd%domain)
        if ( .not. found ) call error_mesg(trim(module_name), &
                       'tile_index axis not found in '//trim(restart_file_name), FATAL)
        allocate(idx(siz(1)),r0d(siz(1)))
        call read_compressed(restart_base_name,'tile_index',idx, domain=lnd%domain, timelevel=1)

        call read_compressed(restart_base_name,'temp',r0d, domain=lnd%domain, timelevel=1)
        call assemble_tiles(cana_T_ptr,idx,r0d)

        do tr = 1, ntcana
          call get_tracer_names(MODEL_LAND, tr, name=name)
          if(field_exist(restart_base_name,name,domain=lnd%domain)) then
            call read_compressed(restart_base_name,name,r0d, domain=lnd%domain, timelevel=1)
            call assemble_tiles(cana_tr_ptr,idx,r0d,tr)
          endif
        enddo
        deallocate(idx,r0d)
     else
        __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
        call read_tile_data_r0d_fptr(unit, 'temp', cana_T_ptr  )
        do tr = 1, ntcana
           call get_tracer_names(MODEL_LAND, tr, name=name)
           if(nfu_inq_var(unit,trim(name))==NF_NOERR) then
              call error_mesg('cana_init','reading tracer "'//trim(name)//'"',NOTE)
              call read_tile_data_r1d_fptr(unit,name,cana_tr_ptr,tr)
           else
              call error_mesg('cana_init', 'tracer "'//trim(name)// &
                   '" was set to initial value '//string(init_tr(tr)), NOTE)
           endif
        enddo
        __NF_ASRT__(nf_close(unit))     
     endif
  else
     call error_mesg('cana_init',&
          'cold-starting canopy air',&
          NOTE)
  endif

  ! initialize options, to avoid expensive character comparisons during 
  ! run-time
  if (trim(turbulence_to_use)=='lm3v') then
     turbulence_option = TURB_LM3V
  else if (trim(turbulence_to_use)=='lm3w') then
     turbulence_option = TURB_LM3W
  else
     call error_mesg('cana_init', 'canopy air turbulence option turbulence_to_use="'// &
          trim(turbulence_to_use)//'" is invalid, use "lm3w" or "lm3v"', FATAL)
  endif
  
end subroutine cana_init


! ============================================================================
! release memory
subroutine cana_end ()

  module_is_initialized =.FALSE.

end subroutine cana_end


! ============================================================================
subroutine save_cana_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  integer :: unit            ! restart file i/o unit
  character(len=32)  :: name,units
  character(len=128) :: longname
  integer :: tr

  call error_mesg('cana_end','writing NetCDF restart',NOTE)
  call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'cana.res.nc',&
          lnd%coord_glon, lnd%coord_glat, cana_tile_exists, tile_dim_length)

     ! write temperature
     call write_tile_data_r0d_fptr(unit,'temp' ,cana_T_ptr,'canopy air temperature','degrees_K')
     ! write all tracers
     do tr = 1,ntcana
        call get_tracer_names(MODEL_LAND, tr, name, longname, units)
        if (tr==ico2.and..not.save_qco2) cycle
        call write_tile_data_r1d_fptr(unit,name,cana_tr_ptr,tr,'canopy air '//trim(longname),trim(units))
     enddo
     ! close output file
     __NF_ASRT__(nf_close(unit))
end subroutine save_cana_restart

! ============================================================================
subroutine save_cana_restart_new(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  character(len=267) :: fname
  type(restart_file_type) :: cana_restart ! restart file i/o object
  integer, allocatable :: idx(:)
  real, allocatable :: cana_T(:)
  real, allocatable :: cana_tr(:,:)
  integer :: id_restart, isize, tr
  character(len=32)  :: name,units
  character(len=128) :: longname

  call error_mesg('cana_end','writing new format NetCDF restart',NOTE)
! must set domain so that io_domain is available
  call set_domain(lnd%domain)
! Note that fname is updated for tile & rank numbers during file creation
  fname = trim(timestamp)//'cana.res.nc'
  call create_tile_out_file(cana_restart,idx,fname,cana_tile_exists,tile_dim_length)
  isize = size(idx)

  allocate(cana_T(isize))
  allocate(cana_tr(isize,ntcana))

  ! Output data provides signature
  call gather_tile_data(cana_T_ptr,idx,cana_T)
  id_restart = register_restart_field(cana_restart,fname,'temp',cana_T, &
                                      longname='canopy air temperature',units='degrees_K')

  do tr = 1,ntcana
     call get_tracer_names(MODEL_LAND, tr, name, longname, units)
     if (tr==ico2.and..not.save_qco2) cycle
     call gather_tile_data(cana_tr_ptr,idx,tr,cana_tr(:,tr))
     id_restart = register_restart_field(cana_restart,fname,name,cana_tr(:,tr), &
                                         longname='canopy air '//trim(longname), units=trim(units))
  enddo

  ! save performs io domain aggregation through mpp_io as with regular domain data
  call save_restart(cana_restart)

  deallocate(idx,cana_T,cana_tr)

  call free_restart_type(cana_restart)
  call nullify_domain()
end subroutine save_cana_restart_new


! ============================================================================
! set up constants for linearization of radiative transfer, using information
! provided by soil, snow and vegetation modules.
subroutine cana_radiation (lm2, &
     subs_refl_dir, subs_refl_dif, subs_refl_lw, & 
     snow_refl_dir, snow_refl_dif, snow_refl_lw, & 
     snow_area, &
     vegn_refl_dir, vegn_refl_dif, vegn_tran_dir, vegn_tran_dif, &
     vegn_tran_dir_dir, vegn_refl_lw, vegn_tran_lw,  &
     vegn_cover, &
     Sg_dir, Sg_dif, Sv_dir, Sv_dif, &
     land_albedo_dir, land_albedo_dif )

  logical, intent(in) :: lm2
  real, intent(in) :: &
       subs_refl_dir(NBANDS), subs_refl_dif(NBANDS), subs_refl_lw,  & ! sub-snow reflectances for direct, diffuse, and LW radiation respectively
       snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_refl_lw,  & ! snow reflectances for direct, diffuse, and LW radiation respectively
       snow_area, &
       vegn_refl_dir(NBANDS), vegn_tran_dir(NBANDS), & ! vegn reflectance & transmittance for direct light
       vegn_tran_dir_dir(NBANDS), & !
       vegn_refl_dif(NBANDS), vegn_tran_dif(NBANDS), & ! vegn reflectance & transmittance for diffuse light 
       vegn_refl_lw,  vegn_tran_lw,  & ! vegn reflectance & transmittance for thermal radiation 
       vegn_cover

  real, intent(out) :: &
     Sg_dir(NBANDS), Sg_dif(NBANDS), & ! fraction of downward short-wave absorbed by ground and snow
     Sv_dir(NBANDS), Sv_dif(NBANDS), & ! fraction of downward short-wave absorbed by vegetation
     land_albedo_dir(NBANDS), land_albedo_dif(NBANDS)

  ! ---- local vars
  real :: &
       grnd_refl_dir(NBANDS), & ! SW reflectances of ground surface (by spectral band)
       grnd_refl_dif(NBANDS), & ! SW reflectances of ground surface (by spectral band)
       grnd_refl_lw             ! LW reflectance of ground surface
  real :: &
       subs_up_from_dir(NBANDS), subs_up_from_dif(NBANDS), &
       subs_dn_dir_from_dir(NBANDS),  subs_dn_dif_from_dif(NBANDS), subs_dn_dif_from_dir(NBANDS)

  grnd_refl_dir = subs_refl_dir + (snow_refl_dir - subs_refl_dir) * snow_area
  grnd_refl_dif = subs_refl_dif + (snow_refl_dif - subs_refl_dif) * snow_area
  grnd_refl_lw  = subs_refl_lw  + (snow_refl_lw  - subs_refl_lw ) * snow_area

  ! ---- shortwave -----------------------------------------------------------
  ! allocation to canopy and ground, based on solution for single
  ! vegetation layer of limited cover. both ground and vegetation are gray.
  Sv_dir = 0; Sv_dif = 0
  IF (LM2) THEN

     subs_dn_dir_from_dir = vegn_tran_dir_dir
     subs_dn_dif_from_dir = vegn_tran_dir
     subs_dn_dif_from_dif = vegn_tran_dif
     if (sfc_dir_albedo_bug) then
        subs_up_from_dir = grnd_refl_dir &
             * (subs_dn_dir_from_dir + subs_dn_dif_from_dir)
     else
        subs_up_from_dir = grnd_refl_dir*subs_dn_dir_from_dir + &
                           grnd_refl_dif*subs_dn_dif_from_dir
     endif
     subs_up_from_dif = grnd_refl_dif*subs_dn_dif_from_dif
     land_albedo_dir = subs_up_from_dir+vegn_refl_dir
     land_albedo_dif = subs_up_from_dif+vegn_refl_dif

  ELSE

     subs_dn_dir_from_dir = vegn_tran_dir_dir
     subs_dn_dif_from_dir = (vegn_tran_dir + vegn_refl_dif*grnd_refl_dir*vegn_tran_dir_dir)&
                          / (1 - grnd_refl_dif*vegn_refl_dif)
     subs_dn_dif_from_dif = vegn_tran_dif &
                          / (1 - grnd_refl_dif*vegn_refl_dif)
     if (sfc_dir_albedo_bug) then
        subs_up_from_dir = grnd_refl_dir * (subs_dn_dir_from_dir + subs_dn_dif_from_dir)
     else
        subs_up_from_dir = grnd_refl_dir*subs_dn_dir_from_dir + &
                           grnd_refl_dif*subs_dn_dif_from_dir
     endif
     subs_up_from_dif = grnd_refl_dif*subs_dn_dif_from_dif
     land_albedo_dir  = subs_up_from_dir*vegn_tran_dif + vegn_refl_dir
     land_albedo_dif  = subs_up_from_dif*vegn_tran_dif + vegn_refl_dif

  ENDIF

  Sg_dir = subs_dn_dir_from_dir + subs_dn_dif_from_dir - subs_up_from_dir
  Sg_dif = subs_dn_dif_from_dif - subs_up_from_dif
  Sv_dir = 1 - Sg_dir - land_albedo_dir
  Sv_dif = 1 - Sg_dif - land_albedo_dif

end subroutine cana_radiation


! ============================================================================
subroutine cana_turbulence (u_star,&
     vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf, &
     land_d, land_z0m, land_z0s, grnd_z0s, &
     con_v_h, con_v_v, con_g_h, con_g_v )
  real, intent(in) :: &
       u_star, & ! friction velocity, m/s
       land_d, land_z0m, land_z0s, grnd_z0s, & 
       vegn_cover, vegn_height, &
       vegn_lai, vegn_sai, vegn_d_leaf
  real, intent(out) :: &
       con_v_h, con_v_v, & ! one-sided foliage-cas conductance per unit of ground area
       con_g_h, con_g_v    ! ground-CAS conductance per unit ground area

  !---- local constants
  real, parameter :: a_max = 3
  real, parameter :: leaf_co = 0.01 ! meters per second^(1/2)
                                    ! leaf_co = g_b(z)/sqrt(wind(z)/d_leaf)
  ! ---- local vars 
  real :: a        ! parameter of exponential wind profile within canopy:
                   ! u = u(ztop)*exp(-a*(1-z/ztop))
  real :: height   ! effective height of vegetation
  real :: wind     ! normalized wind on top of canopy, m/s
  real :: Kh_top   ! turbulent exchange coefficient on top of the canopy
  real :: vegn_idx ! total vegetation index = LAI+SAI
  real :: rah_sca  ! ground-SCA resistance

  vegn_idx = vegn_lai+vegn_sai  ! total vegetation index
  select case(turbulence_option)
  case(TURB_LM3W)
     if(vegn_cover > 0) then
        wind  = u_star/VONKARM*log((vegn_height-land_d)/land_z0m) ! normalized wind on top of the canopy
        a     = vegn_cover*a_max
        con_v_h = (2*max(vegn_lai,lai_min_turb)*leaf_co*(1-exp(-a/2))/a)*sqrt(wind/vegn_d_leaf)
        con_g_h = u_star*a*VONKARM*(1-land_d/vegn_height) &
             / (exp(a*(1-grnd_z0s/vegn_height)) - exp(a*(1-(land_z0s+land_d)/vegn_height)))
     else
        con_v_h = 0
        con_g_h = 0
     endif
  case(TURB_LM3V)
     height = max(vegn_height,0.1) ! effective height of the vegetation
     a = a_max
     wind=u_star/VONKARM*log((height-land_d)/land_z0m) ! normalized wind on top of the canopy
  
     con_v_h = (2*max(vegn_lai,lai_min_turb)*leaf_co*(1-exp(-a/2))/a)*sqrt(wind/vegn_d_leaf)

     if (land_d > 0.06 .and. vegn_idx > 0.25) then
        Kh_top = VONKARM*u_star*(height-land_d)
        rah_sca = height/a/Kh_top * &
             (exp(a*(1-grnd_z0s/height)) - exp(a*(1-(land_z0m+land_d)/height)))
        rah_sca = min(rah_sca,1250.0)
     else
        rah_sca = bare_rah_sca
     endif
     con_g_h = 1.0/rah_sca
  end select
  con_g_v = con_g_h
  con_v_v = con_v_h
end subroutine

! ============================================================================
! update effective surface roughness lengths for CAS-to-atmosphere fluxes
! and conductances for canopy-to-CAS and ground-to-CAS fluxes
!
! Strategy: Always define a canopy present. Non-vegetated situation is simply
! a limit as vegetation density approaches (but isn't allowed to reach) zero.
! Create expressions for the outputs that reduce to the special
! cases of full canopy cover and no canopy. Full canopy solution is that
! from Bonan (NCAR/TN-417+STR, 1996, p. 63). Thus, setting cover=1 in
! recovers Bonan. Letting cover approach 0 makes con_v_coef go to zero,
! preventing exchange with canopy, and makes con_g_coef go infinite,
! removing sub-canopy resistance and putting all resistance above the
! canopy, where it can be affected by stability adjustments.
!
! ** However, there is still a problem with this formulation when a
! canopy is present, because surface flux (I think) is not told to
! subtract out the resistances associated with con_v_coef and con_g_coef,
! which thus seem to be double-counted. For testing LM2, we should set them
! to zero anyway.
subroutine cana_roughness(lm2, &
     subs_z0m, subs_z0s, &
     snow_z0m, snow_z0s, snow_area, &
     vegn_cover, vegn_height, vegn_lai, vegn_sai, &
     land_d, land_z0m, land_z0s, is_lake_or_glac )
  logical, intent(in) :: lm2, is_lake_or_glac
  real, intent(in) :: &
       subs_z0m, subs_z0s, snow_z0m, snow_z0s, snow_area, vegn_cover, vegn_height, &
       vegn_lai, vegn_sai
  real, intent(out) :: &
       land_d    ,&
       land_z0m  ,&
       land_z0s

  !---- local constants
  real, parameter :: d_h_max = 2./3.
  real, parameter :: z0m_h_max = 1/7.35

  ! ---- local vars 
  real :: d_h      ! ratio of displacement height to vegetation height
  real :: z0m_h    ! ratio of roughness length to vegetation height
  real :: grnd_z0m, grnd_z0s
  real :: z0s_h, z0s_h_max
  real :: vegn_idx ! total vegetation index = LAI+SAI
  real :: height   ! effective vegetation height

  grnd_z0m = exp( (1-snow_area)*log(subs_z0m) + snow_area*log(snow_z0m))
  grnd_z0s = exp( (1-snow_area)*log(subs_z0s) + snow_area*log(snow_z0s))

  select case(turbulence_option)
  case(TURB_LM3W)
     if(vegn_cover > 0) then
        z0s_h_max = z0m_h_max*grnd_z0s/grnd_z0m ! to ensure cover->0 limit works
        d_h = vegn_cover*d_h_max
        if (lm2) then
           if (vegn_lai.gt.1) then  ! TEMP ***
              z0m_h = z0m_h_max
              z0s_h = z0s_h_max
           else
              z0m_h = grnd_z0m/vegn_height
              z0s_h = grnd_z0s/vegn_height
           endif
        else
           z0m_h = exp( vegn_cover*log(z0m_h_max) + (1-vegn_cover)*log(grnd_z0m/vegn_height))
           z0s_h = exp( vegn_cover*log(z0s_h_max) + (1-vegn_cover)*log(grnd_z0s/vegn_height))
        endif
        land_d   = d_h*vegn_height
        land_z0m = z0m_h*vegn_height
        land_z0s = z0s_h*vegn_height
     else
        land_d   = 0
        land_z0m = grnd_z0m
        land_z0s = grnd_z0s
     endif
     
  case(TURB_LM3V)
     height = max(vegn_height,0.1) ! effective height of the vegetation
     vegn_idx = vegn_lai+vegn_sai  ! total vegetation index
     if(vegn_idx>1e-4) then
        land_d = 1.1*height*log(1+(0.07*vegn_idx)**0.25)
        if(vegn_idx>2.85) then
           land_z0m = 0.3*(height-land_d)
        else
           land_z0m = grnd_z0m + 0.3*height*sqrt(0.07*vegn_idx)
        endif
     else 
        ! bare soil or leaf off
        land_z0m = 0.1 *height
        land_d   = 0.66*height
     endif
     land_z0s = land_z0m*exp(-2.0) 

     if (allow_small_z0.and.is_lake_or_glac) then
         land_d = 0
	 land_z0m = grnd_z0m
	 land_z0s = grnd_z0s
       endif
     
  end select

end subroutine cana_roughness

! ============================================================================
subroutine cana_state ( cana, cana_T, cana_q, cana_co2 )
  type(cana_tile_type), intent(in)  :: cana
  real, optional      , intent(out) :: cana_T, cana_q, cana_co2

  if (present(cana_T))   cana_T   = cana%T
  if (present(cana_q))   cana_q   = cana%tr(isphum)
  if (present(cana_co2)) cana_co2 = cana%tr(ico2)
  
end subroutine

! ============================================================================
subroutine cana_step_1 ( cana,&
     p_surf, con_g_h, con_g_v, grnd_T, grnd_rh, grnd_rh_psi, &
     Hge,  DHgDTg, DHgDTc,    &
     Ege,  DEgDTg, DEgDqc, DEgDpsig     )
  type(cana_tile_type), intent(in) :: cana
  real, intent(in) :: &
     p_surf,  & ! surface pressure, Pa
     con_g_h, & ! conductivity between ground and CAS for heat
     con_g_v, & ! conductivity between ground and CAS for vapor
     grnd_T,  & ! ground temperature, degK
     grnd_rh, & ! ground relative humidity
     grnd_rh_psi ! psi derivative of ground relative humidity
  real, intent(out) ::   &
     Hge,  DHgDTg, DHgDTc, & ! linearization of the sensible heat flux from ground
     Ege,  DEgDTg, DEgDqc, DEgDpsig    ! linearization of evaporation from ground

  ! ---- local vars
  real :: rho, grnd_q, qsat, DqsatDTg

  call check_temp_range(grnd_T,'cana_step_1','grnd_T')

  call qscomp(grnd_T,p_surf,qsat,DqsatDTg)
  grnd_q = grnd_rh * qsat

  rho      =  p_surf/(rdgas*cana%T*(1+d608*cana%tr(isphum)))
  Hge      =  rho*cp_air*con_g_h*(grnd_T - cana%T)
  DHgDTg   =  rho*cp_air*con_g_h
  DHgDTc   = -rho*cp_air*con_g_h
  Ege      =  rho*con_g_v*(grnd_q  - cana%tr(isphum))
  DEgDTg   =  rho*con_g_v*DqsatDTg*grnd_rh
  DEgDqc   = -rho*con_g_v
  DEgDpsig =  rho*con_g_v*qsat*grnd_rh_psi
  if(is_watch_point())then
     write(*,*)'#### cana_step_1 input ####'
     __DEBUG1__(p_surf)
     __DEBUG2__(con_g_h,con_g_v)
     __DEBUG2__(grnd_T,grnd_rh)
     write(*,*)'#### cana_step_1 internals ####'
     __DEBUG4__(rho, grnd_q, qsat, DqsatDTg)
     __DEBUG2__(cana%T,cana%tr(isphum))
     write(*,*)'#### cana_step_1 output ####'
     __DEBUG3__(Hge,  DHgDTg, DHgDTc)
     __DEBUG4__(Ege,  DEgDTg, DEgDqc, DEgDpsig)
  endif
end subroutine cana_step_1 


! ============================================================================
subroutine cana_step_2 ( cana, delta_Tc, delta_qc )
  type(cana_tile_type), intent(inout) :: cana
  real, intent(in) ::  &
     delta_Tc, & ! change in canopy air temperature
     delta_qc    ! change in canopy air humidity

  cana%T = cana%T + delta_Tc
  cana%tr(isphum) = cana%tr(isphum) + delta_qc
end subroutine cana_step_2

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function cana_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   cana_tile_exists = associated(tile%cana)
end function cana_tile_exists

! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
#define DEFINE_CANA_ACCESSOR_0D(xtype,x) subroutine cana_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%cana))p=>t%cana%x;endif;end subroutine
#define DEFINE_CANA_ACCESSOR_1D(xtype,x) subroutine cana_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p(:);p=>NULL();if(associated(t))then;if(associated(t%cana))p=>t%cana%x;endif;end subroutine

DEFINE_CANA_ACCESSOR_0D(real,T)
DEFINE_CANA_ACCESSOR_1D(real,tr)

end module canopy_air_mod
