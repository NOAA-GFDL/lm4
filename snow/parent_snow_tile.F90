module parent_snow_tile_mod
#include <fms_platform.h>
#include "shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : file_exist, check_nml_error, close_file, stdlog, FATAL, NOTE
use constants_mod,only: tfreeze, hlf
use land_constants_mod, only : NBANDS
use land_tile_selectors_mod, only : tile_selector_type
use land_data_mod, only : log_version
use land_debug_mod, only : is_watch_point
use snowpack_mod, only : snowpack_t, use_mcm_masking, depth_crit

implicit none
private

! ==== public interfaces =====================================================
public :: read_snow_data_namelist 
public :: read_snow_data_namelist_brief 
public :: snow_data_thermodynamics 
public :: snow_data_hydraulics 
public :: snow_data_area 
public :: snow_radiation 
public :: mc_fict, z0_momentum, k_over_B, num_l, dz, distinct_snow_on_glacier
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'parent_snow_tile_mod'
#include "../shared/version_variable.inc"

integer, parameter, public :: max_lev = 10

! ! from the modis brdf/albedo product user's guide:
real, parameter :: g_iso  = 1.
real, parameter :: g_vol  = 0.189184
real, parameter :: g_geo  = -1.377622
real, parameter :: g0_iso = 1.0
real, parameter :: g1_iso = 0.0
real, parameter :: g2_iso = 0.0
real, parameter :: g0_vol = -0.007574
real, parameter :: g1_vol = -0.070987
real, parameter :: g2_vol =  0.307588
real, parameter :: g0_geo = -1.284909
real, parameter :: g1_geo = -0.166314
real, parameter :: g2_geo =  0.041840

! range of temperatures for ramp between "warm" and "cold" albedo
real, parameter :: t_range = 10.0 ! degK

! ==== types =================================================================




type, abstract, public :: snow_tile_type
  ! variables common to the two snow models:
  integer :: tag ! kind of the tile
  integer :: nlayers !< number of snow layers
  ! variables needed for old snow model only:
  real, allocatable :: wl(:)
  real, allocatable :: ws(:)
  real, allocatable :: T(:)
  real, allocatable :: e(:), f(:)
  type(snowpack_t) :: sp ! structure with data for glass snow model
  contains
  procedure(func_snow_is_selected),   deferred :: snow_is_selected
  procedure(func_snow_roughness),   deferred :: snow_roughness
  procedure(func_stock_pe),   deferred :: stock_pe
  procedure(func_snow_active),   deferred :: snow_active
  procedure(func_snow_tile_heat),   deferred :: snow_tile_heat
  procedure(func_snow_get_sfc_temp),   deferred :: snow_get_sfc_temp
  procedure(func_merge_snow_tiles), deferred :: merge_snow_tiles
  procedure(func_get_snow_tile_tag), deferred :: get_snow_tile_tag
  procedure(func_snow_get_wsi), deferred :: get_wsi
  procedure(func_snow_get_wli), deferred :: get_wli
  procedure(func_snow_get_Ti), deferred :: get_Ti
  procedure(func_snow_set_wsi), deferred :: set_wsi
  procedure(func_snow_set_wli), deferred :: set_wli
  procedure(func_snow_set_Ti), deferred :: set_Ti
  procedure(func_get_snow_total_ice), deferred :: ice
  procedure(func_get_snow_total_liq), deferred :: liq
end type snow_tile_type

abstract interface
  ! module procedures
  subroutine func_merge_snow_tiles(snow2, w2, snow1, w1)
    import :: snow_tile_type
    real, intent(in) :: w1
    real, intent(in) :: w2
    class(snow_tile_type), intent(in) :: snow1
    class(snow_tile_type), intent(inout) :: snow2
  end subroutine func_merge_snow_tiles

  integer function func_get_snow_tile_tag(snow) result(tag)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
  end function func_get_snow_tile_tag

  logical function func_snow_is_selected(snow, sel) result(cm1_snow_is_selected)
    import :: snow_tile_type
    import :: tile_selector_type
    type(tile_selector_type),  intent(in) :: sel
    class(snow_tile_type), intent(in) :: snow
  end function func_snow_is_selected

  subroutine func_snow_roughness(snow, snow_z0s, snow_z0m)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
    real, intent(out):: snow_z0s, snow_z0m
  end subroutine func_snow_roughness

  subroutine func_stock_pe(snow, twd_liq, twd_sol  )
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
    real,                  intent(out)   :: twd_liq, twd_sol
  end subroutine func_stock_pe

  real function func_snow_tile_heat(snow) result(heat)
    import :: snow_tile_type
    class(snow_tile_type), intent(in)  :: snow
  end function func_snow_tile_heat

  logical function func_snow_active(snow) result(snow_active)
    import :: snow_tile_type
    class(snow_tile_type), intent(in)  :: snow
  end function func_snow_active

  subroutine func_snow_get_sfc_temp(snow, snow_T)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
    real,                  intent(out)   :: snow_T
  end subroutine func_snow_get_sfc_temp

  real function func_snow_get_wli(snow, i) result(res)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
    integer, intent(in) :: i
  end function func_snow_get_wli

  real function func_snow_get_wsi(snow, i) result(res)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
    integer, intent(in) :: i
  end function func_snow_get_wsi

  real function func_snow_get_Ti(snow, i) result(res)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
    integer, intent(in) :: i
  end function func_snow_get_Ti

  subroutine func_snow_set_wli(snow, i, v) 
    import :: snow_tile_type
    class(snow_tile_type), intent(inout) :: snow
    integer, intent(in) :: i
    real, intent(in) :: v
  end subroutine func_snow_set_wli

  subroutine func_snow_set_wsi(snow, i, v) 
    import :: snow_tile_type
    class(snow_tile_type), intent(inout) :: snow
    integer, intent(in) :: i
    real, intent(in) :: v
  end subroutine func_snow_set_wsi

  subroutine func_snow_set_Ti(snow, i, v) 
    import :: snow_tile_type
    class(snow_tile_type), intent(inout) :: snow
    integer, intent(in) :: i
    real, intent(in) :: v
  end subroutine func_snow_set_Ti

  real function func_get_snow_total_ice(snow) result(ice)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
  end function func_get_snow_total_ice

  real function func_get_snow_total_liq(snow) result(liq)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
  end function func_get_snow_total_liq

end interface



! ==== module data ===========================================================
logical, public :: use_brdf ! not protected because it is set in snow.F90

!---- namelist ---------------------------------------------------------------
character(len=16), PUBLIC:: snow_option = 'cm'  ! or 'gl' later on
! logical :: use_mcm_masking       = .false.   ! MCM snow mask fn
real    :: w_sat                 = 670.
real    :: psi_sat               = -0.06
real    :: k_sat                 = 0.02
real    :: chb                   = 3.5
real    :: thermal_cond_ref      = 0.3
! real    :: depth_crit            = 0.0167
real    :: z0_momentum           = 0.001
real    :: refl_snow_max_dir(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
real    :: refl_snow_max_dif(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
real    :: refl_snow_min_dir(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
real    :: refl_snow_min_dif(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
real    :: emis_snow_max         = 0.95      ! reset to 1 for MCM
real    :: emis_snow_min         = 0.90      ! reset to 1 for MCM
real    :: k_over_B              = 2         ! reset to 0 for MCM
integer :: num_l                 = 3         ! number of snow levels
real    :: dz(max_lev)           = (/0.1,0.8,0.1,0.,0.,0.,0.,0.,0.,0./)
                                              ! rel. thickness of model layers,
                                              ! from top down
! real, protected, public :: &
!    cpw = 1952.0, &  ! specific heat of water vapor at constant pressure
!    clw = 4218.0, &  ! specific heat of water (liquid)
!    csw = 2106.0     ! specific heat of water (ice)
real    :: mc_fict = 10. * 4218 ! additional (fictitious) soil heat capacity (for numerical stability?).
! from analysis of modis data (ignoring temperature dependence):
  real :: f_iso_cold(NBANDS) = (/ 0.354, 0.530 /)
  real :: f_vol_cold(NBANDS) = (/ 0.200, 0.252 /)
  real :: f_geo_cold(NBANDS) = (/ 0.054, 0.064 /)
  real :: f_iso_warm(NBANDS) = (/ 0.354, 0.530 /)
  real :: f_vol_warm(NBANDS) = (/ 0.200, 0.252 /)
  real :: f_geo_warm(NBANDS) = (/ 0.054, 0.064 /)

logical :: distinct_snow_on_glacier = .FALSE. ! if TRUE, the following parameters define
           ! reflectance of snow on glaciers, otherwise snow reflectance does not depend
           ! on the underlying surface (except overlap).
real :: f_iso_cold_on_glacier(NBANDS) = (/ 0.354, 0.530 /)
real :: f_vol_cold_on_glacier(NBANDS) = (/ 0.200, 0.252 /)
real :: f_geo_cold_on_glacier(NBANDS) = (/ 0.054, 0.064 /)
real :: f_iso_warm_on_glacier(NBANDS) = (/ 0.354, 0.530 /)
real :: f_vol_warm_on_glacier(NBANDS) = (/ 0.200, 0.252 /)
real :: f_geo_warm_on_glacier(NBANDS) = (/ 0.054, 0.064 /)
real :: refl_snow_max_dir_on_glacier(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
real :: refl_snow_max_dif_on_glacier(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
real :: refl_snow_min_dir_on_glacier(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
real :: refl_snow_min_dif_on_glacier(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM

namelist /snow_data_nml/  w_sat,                    &
     psi_sat,                k_sat,                 &
     chb,                                           &
     thermal_cond_ref,                              &
     z0_momentum,                                   &
     f_iso_cold, f_vol_cold, f_geo_cold, &
     f_iso_warm, f_vol_warm, f_geo_warm, &
     refl_snow_max_dir,    refl_snow_min_dir,   &
     refl_snow_max_dif,    refl_snow_min_dif,   &
     emis_snow_max,          emis_snow_min,         &
     k_over_B,             &
     num_l,                   dz, mc_fict, &
! snow radiative parameters on glacier
     distinct_snow_on_glacier, &
     f_iso_cold_on_glacier, f_vol_cold_on_glacier, f_geo_cold_on_glacier, &
     f_iso_warm_on_glacier, f_vol_warm_on_glacier, f_geo_warm_on_glacier, &
     refl_snow_max_dir_on_glacier,    refl_snow_min_dir_on_glacier,   &
     refl_snow_max_dif_on_glacier,    refl_snow_min_dif_on_glacier, snow_option

! ---- end of namelist --------------------------------------------------------

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



! ============================================================================
subroutine read_snow_data_namelist_brief()
! subroutine read_snow_data_namelist()
  ! integer, intent(out) :: snow_num_l
  ! real,    intent(out) :: snow_dz(:)
  ! real,    intent(out) :: snow_mc_fict

  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=snow_data_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=snow_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'snow_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  unit=stdlog()
  write(unit, nml=snow_data_nml)

  ! initialize global module data here

  ! set up output arguments
  ! snow_num_l = num_l
  ! snow_dz    = dz
  ! snow_mc_fict = mc_fict

end subroutine read_snow_data_namelist_brief


! ============================================================================
subroutine read_snow_data_namelist(snow_num_l, snow_dz, snow_mc_fict)
! subroutine read_snow_data_namelist()
  integer, intent(out) :: snow_num_l
  real,    intent(out) :: snow_dz(:)
  real,    intent(out) :: snow_mc_fict

  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=snow_data_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=snow_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'snow_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  unit=stdlog()
  write(unit, nml=snow_data_nml)

  ! initialize global module data here

  ! set up output arguments
  snow_num_l = num_l
  snow_dz    = dz
  snow_mc_fict = mc_fict


    ! EZSNOW - check nml values here
  if(is_watch_point()) then
    write(*,*) "EZNML CHECK - READ_SNOW_DATA_NAMELIST"
      __DEBUG1__(w_sat) 
      __DEBUG1__(psi_sat) 
      __DEBUG1__(k_sat) 
      __DEBUG1__(chb) 
      __DEBUG1__(thermal_cond_ref) 
      __DEBUG1__(z0_momentum) 
      __DEBUG1__(refl_snow_max_dir)
      __DEBUG1__(refl_snow_max_dif)
      __DEBUG1__(refl_snow_min_dir)
      __DEBUG1__(refl_snow_min_dif)
      __DEBUG1__(emis_snow_max) 
      __DEBUG1__(emis_snow_min) 
      __DEBUG1__(k_over_B) 
      __DEBUG1__(num_l) 
      __DEBUG1__(dz)
      __DEBUG1__(mc_fict)
      __DEBUG1__(f_iso_cold)
      __DEBUG1__(f_vol_cold)
      __DEBUG1__(f_geo_cold)
      __DEBUG1__(f_iso_warm)
      __DEBUG1__(f_vol_warm)
      __DEBUG1__(f_geo_warm)
      __DEBUG1__(distinct_snow_on_glacier)
      __DEBUG1__(f_iso_cold_on_glacier)
      __DEBUG1__(f_vol_cold_on_glacier)
      __DEBUG1__(f_geo_cold_on_glacier)
      __DEBUG1__(f_iso_warm_on_glacier)
      __DEBUG1__(f_vol_warm_on_glacier)
      __DEBUG1__(f_geo_warm_on_glacier)
      __DEBUG1__(refl_snow_max_dir_on_glacier)
      __DEBUG1__(refl_snow_max_dif_on_glacier)
      __DEBUG1__(refl_snow_min_dir_on_glacier)
      __DEBUG1__(refl_snow_min_dif_on_glacier)
    endif


end subroutine read_snow_data_namelist

! ============================================================================
! compute snow thermodynmamic properties.
subroutine snow_data_thermodynamics ( snow_rh, thermal_cond)
  real, intent(out) :: snow_rh
  real, intent(out) :: thermal_cond(:)

  ! snow surface assumed to have air at saturation
  snow_rh = 1

  ! these will eventually be functions of water contents and T.
  thermal_cond  = thermal_cond_ref

end subroutine snow_data_thermodynamics


! ============================================================================
! compute snow hydraulic properties (assumed dependent only on wl)
subroutine snow_data_hydraulics (wl, ws, psi, hyd_cond )
  real, intent(in),  dimension(:) :: wl, ws
  real, intent(out), dimension(:) :: psi, hyd_cond

  ! ---- local vars
  integer :: l

  do l = 1, num_l
    psi     (l) = psi_sat *(w_sat/(wl(l)+ws(l)))**chb
    hyd_cond(l) = k_sat*(wl(l)/w_sat)**(3+2*chb)
  enddo

end subroutine snow_data_hydraulics


! ============================================================================
! compute snow area
subroutine snow_data_area ( snow_depth, snow_area )
    real, intent(in)  :: snow_depth
    real, intent(out) :: snow_area

  snow_area = 0.
  if (use_mcm_masking) then
     snow_area = min(1., 0.5*sqrt(max(0.,snow_depth)/depth_crit))
  else
     snow_area = max(0.,snow_depth) / (max(0.,snow_depth) + depth_crit)
  endif

end subroutine snow_data_area

! ============================================================================
! compute snow properties needed to do soil-canopy-atmos energy balance
subroutine snow_radiation ( snow_T, cosz, on_glacier,&
     snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis )
  real, intent(in) :: snow_T  ! snow temperature, deg K
  real, intent(in) :: cosz ! cosine of zenith angle
  logical, intent(in) :: on_glacier ! TRUE if snow is on glacier
  real, intent(out) :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_refl_lw, snow_emis

  if (on_glacier.and.distinct_snow_on_glacier) then
     call snow_rad_calculations ( snow_T, cosz, &
        f_iso_warm_on_glacier, f_vol_warm_on_glacier, f_geo_warm_on_glacier, &
        f_iso_cold_on_glacier, f_vol_cold_on_glacier, f_geo_cold_on_glacier, &
        refl_snow_min_dir_on_glacier, refl_snow_max_dir_on_glacier, &
        refl_snow_min_dif_on_glacier, refl_snow_max_dif_on_glacier, &
        snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis )
  else
     call snow_rad_calculations ( snow_T, cosz, &
        f_iso_warm, f_vol_warm, f_geo_warm, &
        f_iso_cold, f_vol_cold, f_geo_cold, &
        refl_snow_min_dir, refl_snow_max_dir, &
        refl_snow_min_dif, refl_snow_max_dif, &
        snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis )
  endif
end subroutine snow_radiation

! ============================================================================
subroutine snow_rad_calculations ( snow_T, cosz, &
     f_iso_warm, f_vol_warm, f_geo_warm, &
     f_iso_cold, f_vol_cold, f_geo_cold, &
     refl_snow_min_dir, refl_snow_max_dir, &
     refl_snow_min_dif, refl_snow_max_dif, &
     snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis )
  real, intent(in) :: snow_T  ! snow temperature, deg K
  real, intent(in) :: cosz ! cosine of zenith angle
  real, intent(in), dimension(NBANDS) :: &
     f_iso_warm, f_vol_warm, f_geo_warm, &
     f_iso_cold, f_vol_cold, f_geo_cold, &
     refl_snow_min_dir, refl_snow_max_dir, refl_snow_min_dif, refl_snow_max_dif
  real, intent(out) :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_refl_lw, snow_emis

  ! ---- local vars
  real :: blend
  real :: warm_value_dir(NBANDS), cold_value_dir(NBANDS)
  real :: warm_value_dif(NBANDS), cold_value_dif(NBANDS)
  real :: zenith_angle, zsq, zcu

  blend = max(0.,min(1.,1.-(tfreeze-snow_T)/t_range))
  if (use_brdf) then
     zenith_angle = acos(cosz)
     zsq = zenith_angle*zenith_angle
     zcu = zenith_angle*zsq
     warm_value_dir = f_iso_warm*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                    + f_vol_warm*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                    + f_geo_warm*(g0_geo+g1_geo*zsq+g2_geo*zcu)
     cold_value_dir = f_iso_cold*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                    + f_vol_cold*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                    + f_geo_cold*(g0_geo+g1_geo*zsq+g2_geo*zcu)
     cold_value_dif = g_iso*f_iso_cold + g_vol*f_vol_cold + g_geo*f_geo_cold
     warm_value_dif = g_iso*f_iso_warm + g_vol*f_vol_warm + g_geo*f_geo_warm
  else
     warm_value_dir = refl_snow_min_dir
     cold_value_dir = refl_snow_max_dir
     warm_value_dif = refl_snow_min_dif
     cold_value_dif = refl_snow_max_dif
  endif
  snow_refl_dir = cold_value_dir + blend*(warm_value_dir-cold_value_dir)
  snow_refl_dif = cold_value_dif + blend*(warm_value_dif-cold_value_dif)
  snow_emis     = emis_snow_max + blend*(emis_snow_min-emis_snow_max  )
  snow_refl_lw  = 1 - snow_emis
end subroutine snow_rad_calculations


end module parent_snow_tile_mod
