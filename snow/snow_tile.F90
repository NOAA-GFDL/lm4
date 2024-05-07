module snow_tile_mod
#include <fms_platform.h>
#include "../shared/debug.inc"

use mpp_mod, only: input_nml_file

use fms_mod, only : check_nml_error, stdlog, FATAL, NOTE
use constants_mod,only: tfreeze, hlf
use land_constants_mod, only : NBANDS, &
! MODIS BRDF model parameters
    g_iso, g0_iso, g1_iso, g2_iso, &
    g_vol, g0_vol, g1_vol, g2_vol, &
    g_geo, g0_geo, g1_geo, g2_geo
use land_tile_selectors_mod, only : tile_selector_type
use land_data_mod, only : log_version
use land_debug_mod, only : is_watch_point
use tile_diag_buff_mod, only : diag_buff_type

use snowpack_mod, only : snowpack_t, use_mcm_masking, depth_crit

implicit none
private

! ==== public interfaces =====================================================
public :: read_snow_data_namelist
public :: snow_data_thermodynamics
public :: snow_data_area
public :: snow_radiation
public :: mc_fict, z0_momentum, k_over_B, num_l, dz, distinct_snow_on_glacier
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snow_tile_mod'
#include "../shared/version_variable.inc"

integer, parameter, public :: max_lev = 10

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
  procedure(func_snow_is_selected),    deferred :: snow_is_selected
  procedure(func_snow_roughness),      deferred :: snow_roughness
  procedure(func_snow_rad_prop),       deferred :: radiative_properties
  procedure(func_stock_pe),            deferred :: stock_pe
  procedure(func_snow_active),         deferred :: snow_active
  procedure(func_get_real_0D),         deferred :: snow_tile_heat
  procedure(func_get_real_0D),         deferred :: sfc_temp
  procedure(func_merge_snow_tiles),    deferred :: merge_snow_tiles
  procedure(func_get_snow_tile_tag),   deferred :: get_snow_tile_tag
  procedure(func_get_real_1D),         deferred :: get_wsi
  procedure(func_get_real_1D),         deferred :: get_wli
  procedure(func_get_real_1D),         deferred :: get_Ti
  procedure(func_set_real_1D),         deferred :: set_wsi
  procedure(func_set_real_1D),         deferred :: set_wli
  procedure(func_set_real_1D),         deferred :: set_Ti
  procedure(func_get_real_0D),         deferred :: ice
  procedure(func_get_real_0D),         deferred :: liq
  procedure(func_get_depth_area),      deferred :: get_depth_area

  procedure(func_sweep_tiny),          deferred :: sweep_tiny
  procedure(func_partition_sw),        deferred :: partition_sw

  procedure(func_step1),               deferred :: step1
  procedure(func_step2),               deferred :: step2

  procedure :: send_diag !< Send model-specific output to diagnostics.
                         !! Default implementation does nothing.
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

  subroutine func_snow_rad_prop(snow, &
    cosz, subs_refl_dif, p_atm, on_glacier, &
    snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis)
    import :: snow_tile_type
    class(snow_tile_type), intent(inout) :: snow ! slm: it is only inout because gl_snow stores something (intermediate results?)
    real, intent(in) :: cosz
    real, intent(in) :: subs_refl_dif(:)
    real, intent(in) :: p_atm
    logical, intent(in) :: on_glacier
    real, intent(out) :: snow_refl_dir(:), snow_refl_dif(:)
    real, intent(out) :: snow_refl_lw, snow_emis
  end subroutine

  subroutine func_stock_pe(snow, twd_liq, twd_sol  )
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
    real,                  intent(out)   :: twd_liq, twd_sol
  end subroutine func_stock_pe

  real function func_get_real_0D(snow)
    import :: snow_tile_type
    class(snow_tile_type), intent(in)  :: snow
  end function

  logical function func_snow_active(snow) result(snow_active)
    import :: snow_tile_type
    class(snow_tile_type), intent(in)  :: snow
  end function func_snow_active

  real function func_get_real_1D(snow, i)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
    integer,               intent(in) :: i
  end function

  subroutine func_set_real_1D(snow, i, v)
    import :: snow_tile_type
    class(snow_tile_type), intent(inout) :: snow
    integer, intent(in) :: i
    real,    intent(in) :: v
  end subroutine

  subroutine func_get_depth_area(snow, snow_depth, snow_area)
    import :: snow_tile_type
    class(snow_tile_type), intent(in) :: snow
    real, intent(out) :: snow_depth, snow_area
  end subroutine

  ! removes snow if its amount is tiny
  subroutine func_sweep_tiny(snow, lrunf, frunf, hlrunf, hfrunf, lost_wc_em, lost_wc_im)
    import :: snow_tile_type
    class(snow_tile_type), intent(inout) :: snow
    real, intent(out) :: lrunf, frunf   ! liquid and solid runoff generated by sweeper
    real, intent(out) :: hlrunf, hfrunf ! heat carried by liquid ans solid runoff
    real, intent(out) :: lost_wc_em(:), lost_wc_im(:) ! tracer losses
  end subroutine

  ! Given the shortwave radiation absorbed by snow + substrate (fswg) [W/m2]
  ! as well its direct and diffuse components (fswg_dir, fswg_dif)
  ! partition it between surface of snow (where it was absorbed entirely in old cm snow model)
  ! and, if requested, absorption within the snowpack
  ! and absoirption in the underlying substrate (lake/soil/glacier)
  subroutine func_partition_sw( snow, &
    fswg, fswg_dir, fswg_dif,    & ! input
    fswg_substrate, fswg_surface ) ! output
    import :: snow_tile_type
    class(snow_tile_type), intent(inout) :: snow !< state of snowpack
    real, intent(in)  :: fswg ! total sw absorbed by snow + substrate [W/m2]
    real, intent(in)  :: fswg_dir(:), fswg_dif(:) ! total sw absorbed by snow + substrate (dir only, dif only) [W/m2]
    ! logical, intent(IN) :: assign_substrate_sw_to_surface ! if true, override code and assign excess heat to surface instead that passing it to underlying substrate
    real, intent(out) :: fswg_substrate ! sw radiation passed to substrate [W/m2]
    real, intent(out) :: fswg_surface   ! sw radiation to be absorbed at the surface [W/m2]
  end subroutine

  subroutine func_step1( snow, p_surf, grnd_T, snow_G_Z, snow_G_TZ, &
         snow_active, snow_T, snow_rh, snow_liq, snow_ice, &
         snow_subl, snow_area, snow_G0, snow_DGDT, snow_E_max )
    import :: snow_tile_type
    class(snow_tile_type), intent(inout) :: snow
    real,    intent(in) :: p_surf
    real,    intent(in) :: grnd_T
    real,    intent(in) :: snow_G_Z
    real,    intent(in) :: snow_G_TZ
    logical, intent(out):: snow_active
    real,    intent(out):: &
         snow_T, snow_rh, snow_liq, snow_ice, &
         snow_subl, snow_area, snow_G0, snow_DGDT, &
         snow_E_max
  end subroutine

  subroutine func_step2 ( snow, snow_subl,                     &
       vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
       DTg,  Mg_imp,  evapg,  fswg,  flwg,  sensg,  &
       use_tfreeze_in_grnd_latent, &
       ! output
       subs_DT, &
       subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens,  &
       snow_fsw, snow_flw, snow_sens, &
       snow_levap, snow_fevap, snow_melt, &
       snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
       snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, &
       snow_avrg_T , &
       ! additional input/output added by Enrico for standalone model only
       !    snow_rho, snow_age, snow_sph, snow_optd, & ! average snow properties
       ! heat1, verbose, hfevap, dt, wind_atm, t_atm, &
       dt, wind_atm, t_atm, p_surf, &
       wetdep, drydep, grnd_T_preprec, &
       ! for conservation checks only :
       begw_check, begh_check, &
       G0, DGDTg, snow_G_Z, snow_G_TZ, &
       mass_lai_em_1, mass_lai_im_1, &
       lost_wc_em_st, lost_wc_im_st, &
       lost_wc_em, lost_wc_im)
    import :: snow_tile_type
    class(snow_tile_type), intent(inout) :: snow
    real, intent(in) :: &
       snow_subl, vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec
    real, intent(in) :: &
       DTg, Mg_imp, evapg, fswg, flwg, sensg
    logical, intent(in) :: use_tfreeze_in_grnd_latent
    real, intent(out) :: &
           subs_DT, subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens, &
           snow_fsw, snow_flw, snow_sens, &
           snow_levap, snow_fevap, snow_melt, &
           snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
           snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, snow_avrg_T
     real, intent(out) :: grnd_T_preprec
     !  additional in-out variables
     !  real, intent(out) :: heat1
     !  real, intent(out) :: hfevap ! heat released by subl [kg m^-2 s^-1]
     !  logical, intent(in) :: verbose
     !  real, intent(out) :: delta_heat_DTg
      real, intent(in) :: dt ! delta time step
      real, intent(in) :: wind_atm, t_atm, p_surf
      real, intent(in) :: wetdep(:) ! wet deposition rate of tracers from atmosphere [ppm]
      real, intent(in) :: drydep(:) ! dry deposition rate of tracers from atmosphere [mg/m2/s]
      real, intent(out):: lost_wc_em(:), lost_wc_im(:)
      real, intent(in) :: mass_lai_em_1(:), mass_lai_im_1(:) ! mass of LAIs at beginning of step, for mass cons checks
      real, intent(in) :: lost_wc_em_st(:), lost_wc_im_st(:)
      real, intent(in) :: begw_check, begh_check
      real, intent(in) :: G0, DGDTg, snow_G_Z, snow_G_TZ
  end subroutine

end interface



! ==== module data ===========================================================
logical, public :: use_brdf ! not protected because it is set in snow.F90

!---- namelist ---------------------------------------------------------------
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
     refl_snow_max_dif_on_glacier,    refl_snow_min_dif_on_glacier

! ---- end of namelist --------------------------------------------------------

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
subroutine read_snow_data_namelist(snow_num_l, snow_dz, snow_mc_fict)
  integer, intent(out) :: snow_num_l
  real,    intent(out) :: snow_dz(:)
  real,    intent(out) :: snow_mc_fict

  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
  read (input_nml_file, nml=snow_data_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_data_nml')
  unit=stdlog()
  write(unit, nml=snow_data_nml)

  ! initialize global module data here

  ! set up output arguments
  snow_num_l = num_l
  snow_dz    = dz
  snow_mc_fict = mc_fict
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

subroutine send_diag(snow, diag)
  class(snow_tile_type), intent(inout) :: snow
  type(diag_buff_type),  intent(inout) :: diag !< diagnostic buffer
  ! this default implementation does nothing
end subroutine

end module snow_tile_mod
