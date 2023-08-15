! ============================================================================
! snow model module
! ============================================================================
module snow_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : error_mesg, file_exist, check_nml_error, &
     stdlog, close_file, mpp_pe, mpp_root_pe, FATAL, NOTE, lowercase
use time_manager_mod,   only: time_type_to_real
use constants_mod,      only: tfreeze, hlv, hlf, PI

use land_constants_mod, only : NBANDS

use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_data_mod, only : lnd, log_version
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, get_tile_data
use land_debug_mod, only : is_watch_point

use cm_snow_mod, only: cm_read_snow_namelist, cm_snow_init, cm_snow_end, &
    cm_save_snow_restart, cm_snow_get_depth_area, &
    cm_sweep_tiny_snow, cm_snow_step_1, cm_snow_step_2

use gl_snow_mod, only: gl_read_snow_namelist, gl_snow_init, gl_snow_end, &
    gl_save_snow_restart, gl_snow_get_depth_area

use parent_snow_tile_mod, only : &
     snow_tile_type, read_snow_data_namelist, &
     read_snow_data_namelist_brief, &
     snow_option

use cm_snow_tile_mod, only: cm_snow_tile_type

use gl_snow_tile_mod, only: gl_snow_tile_type

use snow_evolution_mod, only: gl_snow_step_2, gl_sweep_tiny_snow, gl_compute_snow_albedo

use snow_constants_mod, only: NTRACERS 


implicit none
private

! ==== public interfaces =====================================================
public :: read_snow_namelist
public :: snow_init
public :: snow_end
public :: save_snow_restart
public :: snow_get_depth_area ! interface
public :: sweep_tiny_snow ! interface
public :: snow_step_1 ! interface
public :: snow_step_2 ! interface
public :: compute_snow_albedo
public :: partition_sw_heat_in_snow





! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snow_mod'
#include "../shared/version_variable.inc"

contains


subroutine read_snow_namelist()

  call read_snow_data_namelist_brief()

  select case(trim(lowercase(snow_option)))
      case('cm')
         call cm_read_snow_namelist()
      case('gl')
         call gl_read_snow_namelist()
      case default
         call error_mesg( &
            'read_snow_namelist in snow_mod', &
            'snow_option = "'//trim(snow_option)//'" is incorrect, use "cm", or "gl"', FATAL)
  end select
end subroutine read_snow_namelist


! initialize snow model
subroutine snow_init()
  select case(trim(lowercase(snow_option)))
      case('cm')
         call cm_snow_init()
      case('gl')
         call gl_snow_init()
      case default
         call error_mesg( &
            'snow_init in snow_mod', &
            'snow_option = "'//trim(snow_option)//'" is incorrect, use "cm", or "gl"', FATAL)
  end select
end subroutine snow_init


! initialize snow model
subroutine snow_end()
  select case(trim(lowercase(snow_option)))
      case('cm')
         call cm_snow_end()
      case('gl')
         call gl_snow_end()
      case default
         call error_mesg( &
            'snow_end in snow_mod', &
            'snow_option = "'//trim(snow_option)//'" is incorrect, use "cm", or "gl"', FATAL)
  end select
end subroutine snow_end



! initialize snow model
subroutine save_snow_restart(tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name


  select case(trim(lowercase(snow_option)))
      case('cm')
         call cm_save_snow_restart(tile_dim_length, timestamp)
      case('gl')
         call gl_save_snow_restart(tile_dim_length, timestamp)
      case default
         call error_mesg( &
            'save_snow_restart in snow_mod', &
            'snow_option = "'//trim(snow_option)//'" is incorrect, use "cm", or "gl"', FATAL)
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


subroutine sweep_tiny_snow(snow, lrunf, frunf, hlrunf, hfrunf, lost_wc_em, lost_wc_im)
  class(snow_tile_type), intent(inout) :: snow
  real, intent(out) :: lrunf, frunf, hlrunf, hfrunf
  real, intent(out), dimension(NTRACERS) :: lost_wc_em, lost_wc_im

  select type (snow)
  type is (cm_snow_tile_type)
      call cm_sweep_tiny_snow(snow, lrunf, frunf, hlrunf, hfrunf)
      lost_wc_em = 0.0
      lost_wc_im = 0.0
  type is (gl_snow_tile_type)
      call gl_sweep_tiny_snow(snow%sp, lrunf, frunf, hlrunf, hfrunf,lost_wc_em, lost_wc_im)
  class default 
      call error_mesg( &
        'sweep_tiny_snow in snow_tile_mod', &
        'snow tile type is incorrect!', FATAL)
  end select
end subroutine sweep_tiny_snow


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


subroutine partition_sw_heat_in_snow( &
   snow, fswg, fswg_dir, fswg_dif, & ! input
   snow_option_passed, albedo_to_use, use_internal_sources, thresh_snow_depth_swheat, &  ! input
   fswg_substrate, fswg_surface) ! output
   !
   ! Given the shortwave radiation absorbed by snow + substrate (fswg) [W/m2]
   ! as well its direct and diffuse components (fswg_dir, fswg_dif)
   ! partition it between surface of snow (where it was absorbed entirely in old cm snow model)
   ! and, if requested, absorption within the snowpack
   ! andabsoirption in the underlying substrate (lake/soil/glacier)
   !
   class(snow_tile_type), intent(inout) :: snow !< state of snowpack
   real, intent(IN), dimension(NBANDS) :: fswg, fswg_dir, fswg_dif ! total sw absorbed by snow + substrate (total, dir only, dif only) [W/m2]
   real, intent(IN) :: snow_option_passed, albedo_to_use, use_internal_sources, thresh_snow_depth_swheat ! model options
   real, intent(OUT) :: fswg_substrate ! sw radiation passed to substrate [W/m2]
   real, intent(OUT) :: fswg_surface ! sw radiation to be absorbed at surface [W/m2]

   if (ALLOCATED(snow%sp%swheat)) DEALLOCATE(snow%sp%swheat) 
   
   if (trim(lowercase(snow_option)) == 'gl') then
      if (trim(lowercase(albedo_to_use))=='snicar') then
         if ((use_internal_sources) .and. ((snow%sp%depth() > thresh_snow_depth_swheat) & 
                                    .and. (snow%sp%nlayers > 0))) then
            ALLOCATE(snow%sp%swheat(snow%sp%nlayers))
            sum_sw_frac_dir = 0.0 ! init total fractions of sw down absorbed by snowpack
            sum_sw_frac_dif = 0.0 ! init total fractions of sw down absorbed by snowpack
            do il=1,snow%sp%nlayers
               snow%sp%swheat(il) = &     
                     fswg_dir(1) * snow%sp%sw_frac_dir(il, 1) + & 
                     fswg_dif(1) * snow%sp%sw_frac_dif(il, 1) + &
                     fswg_dir(2) * snow%sp%sw_frac_dir(il, 2) + & 
                     fswg_dif(2) * snow%sp%sw_frac_dif(il, 2)    
               sum_sw_frac_dir(1)  = sum_sw_frac_dir(1) + snow%sp%sw_frac_dir(il, 1) 
               sum_sw_frac_dif(1)  = sum_sw_frac_dif(1) + snow%sp%sw_frac_dif(il, 1) 
               sum_sw_frac_dir(2)  = sum_sw_frac_dir(2) + snow%sp%sw_frac_dir(il, 2) 
               sum_sw_frac_dif(2)  = sum_sw_frac_dif(2) + snow%sp%sw_frac_dif(il, 2) 
            enddo
            if ((sum_sw_frac_dir(1)>1.0+1E-7).or. (sum_sw_frac_dif(1)>1.0+1E-7) .or. & 
               (sum_sw_frac_dir(2) >1.0+1E-7).or. (sum_sw_frac_dif(2) >1.0+1E-7)  ) then
               write(*,*) "sum of sw_frac_dir(1):", sum_sw_frac_dir(1)
               write(*,*) "sum of sw_frac_dir(2):", sum_sw_frac_dir(2)
               write(*,*) "sum of sw_frac_dif(1):", sum_sw_frac_dif(1)
               write(*,*) "sum of sw_frac_dif(2):", sum_sw_frac_dif(2)
               call land_error_message("Error in sw sources from SNICAR: a total is larger than 1!", severity=FATAL)
            endif
            if ((sum_sw_frac_dir(1)<0.0-1E-7).or. (sum_sw_frac_dif(1)<0.0-1E-7) .or. & 
               (sum_sw_frac_dir(2) <0.0-1E-7).or. (sum_sw_frac_dif(2) <0.0-1E-7)  ) then
               write(*,*) "sum of sw_frac_dir(1):", sum_sw_frac_dir(1)
               write(*,*) "sum of sw_frac_dir(2):", sum_sw_frac_dir(2)
               write(*,*) "sum of sw_frac_dif(1):", sum_sw_frac_dif(1)
               write(*,*) "sum of sw_frac_dif(2):", sum_sw_frac_dif(2)
               call land_error_message("Error in sw sources from SNICAR: a total is below 0!", severity=FATAL)
            endif
            fswg_surface = 0.0
            fswg_substrate = & 
                  fswg_dir(1) * (1.0 - sum_sw_frac_dir(1)) + &
                  fswg_dif(1) * (1.0 - sum_sw_frac_dif(1)) + &
                  fswg_dir(2) * (1.0 - sum_sw_frac_dir(2)) + &
                  fswg_dif(2) * (1.0 - sum_sw_frac_dif(2))
         else ! albedo = snicar, but do not use internal sw sources
            if (snow%sp%nlayers>0) then
               ALLOCATE(snow%sp%swheat(snow%sp%nlayers))
               snow%sp%swheat = 0.0 ! don't change fswg in this case 
            else
               ALLOCATE(snow%sp%swheat(1))
               snow%sp%swheat = 0.0 ! don't change fswg in this case 
            endif
            fswg_surface = fswg
            fswg_substrate = 0.0
         endif
      else ! snow option = GL but albedo model not SNICAR
         if ((use_internal_sources) .and. ((snow%sp%depth() > thresh_snow_depth_swheat) & 
                                    .and. (snow%sp%nlayers > 0))) then
            call snow%sp%sw_sources(fswg_dir, fswg_dif, fswg_substrate)
            fswg_surface = 0.0
         else
            ALLOCATE(snow%sp%swheat(snow%sp%nlayers))
            snow%sp%swheat = 0.0 
            fswg_surface = fswg
            fswg_substrate = 0.0
         endif
      endif ! end albedo choice for GL snow model option
   else ! case of CM snow model: all sw absorption occurs at the surface (part of surface energy balance)
      fswg_surface=fswg
      fswg_substrate = 0.0
   endif
end subroutine partition_sw_heat_in_snow


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
    real, intent(in) :: wetdep(NTRACERS) ! wet deposition of tracers from atmosphere [ppm]
    real, intent(in) :: drydep(NTRACERS) ! dry 
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



