! ============================================================================
! snow model module
! ============================================================================
module snow_mod

#include "../shared/debug.inc"

use fms_mod, only : error_mesg, FATAL, NOTE, lowercase

use land_constants_mod, only : NBANDS

use land_data_mod, only : log_version
use land_debug_mod, only : is_watch_point, land_error_message

use cm_snow_mod, only: cm_read_snow_namelist, cm_snow_init, cm_snow_end, &
    cm_save_snow_restart, cm_snow_get_depth_area, &
    cm_snow_step_1, cm_snow_step_2

use gl_snow_mod, only: gl_read_snow_namelist, gl_snow_init, gl_snow_end, &
    gl_save_snow_restart, gl_snow_get_depth_area

use snow_tile_mod, only : snow_tile_type
use snow_base_mod, only : read_snow_model_namelist, snow_option, SNOW_CM, SNOW_GL

use cm_snow_tile_mod, only: cm_snow_tile_type

use gl_snow_tile_mod, only: gl_snow_tile_type

use snow_evolution_mod, only: gl_snow_step_2, gl_sweep_tiny_snow, gl_compute_snow_albedo, &
                              albedo_to_use, use_internal_sources, thresh_snow_depth_swheat, &
                              assign_substrate_sw_to_surface

use snow_constants_mod, only: NTRACERS

implicit none
private

! ==== public interfaces =====================================================
public :: read_snow_namelist
public :: snow_init
public :: snow_end
public :: save_snow_restart
public :: snow_get_depth_area ! interface
public :: snow_step_1 ! interface
public :: snow_step_2 ! interface
public :: compute_snow_albedo

! re-export snow model selector
public :: snow_option, SNOW_CM, SNOW_GL

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snow_mod'
#include "../shared/version_variable.inc"

contains


subroutine read_snow_namelist()
  call log_version(version, module_name, &
  __FILE__)

  call read_snow_model_namelist()

  select case(snow_option)
  case(SNOW_CM)
     call cm_read_snow_namelist()
  case(SNOW_GL)
     call gl_read_snow_namelist()
  case default
     call land_error_message('read_snow_namelist: The value of snow_option is invalid. This should never happen. See developer', FATAL)
  end select
end subroutine read_snow_namelist


! initialize snow model
subroutine snow_init()
  select case(snow_option)
  case(SNOW_CM)
     call cm_snow_init()
  case(SNOW_GL)
     call gl_snow_init()
  case default
     call land_error_message('snow_init: The value of snow_option is invalid. This should never happen. See developer', FATAL)
  end select
end subroutine snow_init


! initialize snow model
subroutine snow_end()
  select case(snow_option)
  case(SNOW_CM)
     call cm_snow_end()
  case(SNOW_GL)
     call gl_snow_end()
  case default
     call land_error_message('snow_end: The value of snow_option is invalid. This should never happen. See developer', FATAL)
  end select
end subroutine snow_end


! save snow model restart file
subroutine save_snow_restart(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  select case(snow_option)
  case(SNOW_CM)
     call cm_save_snow_restart(tile_dim_length, timestamp)
  case(SNOW_GL)
     call gl_save_snow_restart(tile_dim_length, timestamp)
  case default
     call land_error_message('snow_end: The value of snow_option is invalid. This should never happen. See developer', FATAL)
  end select
end subroutine save_snow_restart


subroutine snow_get_depth_area(snow, snow_depth, snow_area)
  class(snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_depth, snow_area

  select type (snow)
  type is (cm_snow_tile_type)
      call cm_snow_get_depth_area(snow, snow_depth, snow_area)
   type is (gl_snow_tile_type)
      call gl_snow_get_depth_area(snow, snow_depth, snow_area)
  class default
      call error_mesg( &
        'snow_get_depth_area in snow_tile_mod', &
        'type is incorrect!', FATAL)
  end select
end subroutine snow_get_depth_area


subroutine snow_step_1( snow, snow_G_Z, snow_G_TZ, &
         snow_active, snow_T, snow_rh, snow_liq, snow_ice, &
         snow_subl, snow_area, snow_G0, snow_DGDT )
  class(snow_tile_type), intent(inout) :: snow
  real,                 intent(in) :: snow_G_Z
  real,                 intent(in) :: snow_G_TZ
  logical,              intent(out):: snow_active
  real,                 intent(out):: &
       snow_T, snow_rh, snow_liq, snow_ice, &
       snow_subl, snow_area, snow_G0, snow_DGDT

  select type (snow)
  type is (cm_snow_tile_type)
      call cm_snow_step_1 ( snow, snow_G_Z, snow_G_TZ, &
                         snow_active, snow_T, snow_rh, snow_liq, snow_ice, &
                         snow_subl, snow_area, snow_G0, snow_DGDT )
  class default
      call error_mesg( &
        'snow_step_1 in snow_tile_mod', &
        'type is incorrect!', FATAL)
  end select
end subroutine snow_step_1


subroutine compute_snow_albedo(snow, snow_T, cosz, on_glacier, p_atm, subs_refl_dif, & ! input
                snow_refl_dir, snow_refl_dif)
  class(snow_tile_type), intent(inout) :: snow !< state of snowpack
  real, intent(in) :: p_atm  ! ! atm pressure [Pa] from forcing
  real, intent(in) :: snow_T  ! snow temperature, deg K [get it from s instead?]
  real, intent(in) :: cosz ! cosine of zenith angle
  logical, intent(in) :: on_glacier ! TRUE if snow is on glacier
  real, dimension(NBANDS), intent(IN) :: subs_refl_dif
  real, dimension(NBANDS), intent(OUT) :: snow_refl_dir
  real, dimension(NBANDS), intent(OUT) :: snow_refl_dif
!   real, intent(OUT) :: snow_refl_lw, snow_emis

  select type (snow)
  type is (gl_snow_tile_type)
      call gl_compute_snow_albedo ( snow%sp, snow_T, cosz, on_glacier, p_atm, subs_refl_dif, & ! input
                snow_refl_dir, snow_refl_dif)

  class default
      call error_mesg( &
        'compute snow albedo in snow_mod', &
        'type is incorrect: works only with gl_snow_tile_type!', FATAL)
  end select
end subroutine compute_snow_albedo


subroutine snow_step_2 ( snow, snow_subl,                     &
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
    real, intent(in) :: wetdep(NTRACERS) ! wet deposition rate of tracers from atmosphere [ppm]
    real, intent(in) :: drydep(NTRACERS) ! dry deposition rate of tracers from atmosphere [mg/m2/s]
    real, intent(out), DIMENSION(NTRACERS) :: lost_wc_em, lost_wc_im
    real, intent(in), DIMENSION(NTRACERS) :: mass_lai_em_1, mass_lai_im_1 ! mass of LAIs at beginning of step, for mass cons checks
    real, intent(in), DIMENSION(NTRACERS) :: lost_wc_em_st, lost_wc_im_st
    real, intent(in) :: begw_check, begh_check
    real, intent(in) :: G0, DGDTg, snow_G_Z, snow_G_TZ

  select type (snow)

  type is (cm_snow_tile_type)
      call cm_snow_step_2 ( snow, snow_subl,                     &
                           vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
                           DTg,  Mg_imp,  evapg,  fswg,  flwg,  sensg,  &
                           use_tfreeze_in_grnd_latent, subs_DT, &
                           subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens,  &
                           snow_fsw, snow_flw, snow_sens, &
                           snow_levap, snow_fevap, snow_melt, &
                           snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
                           snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, &
                           snow_avrg_T )
  type is (gl_snow_tile_type)
      call gl_snow_step_2 ( snow%sp, snow_subl,                     &
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
  class default
      call error_mesg( &
        'snow_step_2 in snow_tile_mod', &
        'type is incorrect: must be cm_snow_tile_type, or gl_snow_tile_type!', FATAL)
  end select
end subroutine snow_step_2

end module snow_mod
