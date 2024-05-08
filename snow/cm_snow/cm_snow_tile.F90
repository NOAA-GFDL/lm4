module cm_snow_tile_mod
#include <fms_platform.h>
#include "../../shared/debug.inc"


use fms_mod, only : input_nml_file, check_nml_error, stdlog, mpp_pe, mpp_root_pe, FATAL
use time_manager_mod, only: time_type_to_real
use constants_mod,only: tfreeze, hlf
use land_tile_selectors_mod, only : tile_selector_type
use land_data_mod, only : lnd, log_version
use land_debug_mod, only : is_watch_point, land_error_message

use snow_tile_mod, only : snow_tile_type, mc_fict, z0_momentum, k_over_B, num_l, dz, &
      snow_data_area, snow_data_thermodynamics, snow_radiation
use snowpack_mod, only : cpw, clw, csw

implicit none
private

! ==== public interfaces =====================================================
public :: cm_snow_tile_type
public :: cm_snow_tile_ctor

public :: read_snow_cm_namelist
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'cm_snow_tile_mod'
#include "../../shared/version_variable.inc"
real, parameter :: heat_capacity_retro = 1.6e6


! ==== types =================================================================

type, extends(snow_tile_type) :: cm_snow_tile_type
   ! data structure already defined in parent snow type
contains
    procedure :: merge_snow_tiles => cm_merge_snow_tiles_wrapper
    procedure :: get_snow_tile_tag => cm_get_snow_tile_tag
    procedure :: snow_is_selected => cm_snow_is_selected
    procedure :: n_layers => cm_snow_nlayers
    procedure :: snow_roughness => cm_snow_roughness
    procedure :: radiative_properties => cm_snow_rad_prop
    procedure :: stock_pe => cm_snow_tile_stock_pe
    procedure :: snow_active => cm_snow_active
    procedure :: snow_tile_heat => cm_snow_tile_heat
    procedure :: sfc_temp => cm_snow_get_sfc_temp
    procedure :: get_Ti => cm_snow_get_Ti
    procedure :: get_wli => cm_snow_get_wli
    procedure :: get_wsi => cm_snow_get_wsi
    procedure :: set_Ti => cm_snow_set_Ti
    procedure :: set_wli => cm_snow_set_wli
    procedure :: set_wsi => cm_snow_set_wsi
    procedure :: ice => cm_snow_get_total_ice
    procedure :: liq => cm_snow_get_total_liq
    procedure :: get_depth_area => cm_snow_get_depth_area
    procedure :: lai_im => cm_snow_lai  ! both lai_im and lai_em return zero
    procedure :: lai_em => cm_snow_lai

    procedure :: sweep => cm_sweep_tiny_snow
    procedure :: partition_sw => cm_partition_sw

    procedure :: step1 => cm_snow_step_1
    procedure :: step2 => cm_snow_step_2
end type cm_snow_tile_type

! ---- module data
real :: delta_time ! fast (physical) time step [s]


!---- namelist ---------------------------------------------------------------
logical, public, protected :: retro_heat_capacity  = .false.
logical, public, protected :: lm2  = .false.
logical, public, protected :: steal = .false.
character(16), public, protected :: albedo_to_use = ''  ! or 'brdf-params'
real, public, protected :: max_snow       = 1000.
real, public, protected :: wet_max        = 0.0  ! TEMP, move to snow_data
real, public, protected :: snow_density   = 300. ! TEMP, move to snow_data and generalize
real, public, protected :: init_temp = 260.   ! cold-start snow T
real, public, protected :: init_pack_ws   =   0.
real, public, protected :: init_pack_wl   =   0.
real, public, protected :: min_snow_mass = 0.
logical, public, protected :: prevent_tiny_snow = .FALSE. ! if true, tiny snow is removed at the
   ! beginning of fast time step to avoid numerical issues. There is no harm
   ! in doing that, but it changes answers, so for compatibility with older code
   ! turn it off.

namelist /cm_snow_nml/ retro_heat_capacity, lm2, steal, albedo_to_use, &
                    max_snow, wet_max, snow_density, &
                    init_temp, init_pack_ws, init_pack_wl, &
                    min_snow_mass, prevent_tiny_snow
!---- end of namelist --------------------------------------------------------

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine read_snow_cm_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l            ! layer iterator

  call log_version(version, module_name, &
  __FILE__)
  read (input_nml_file, nml=cm_snow_nml, iostat=io)
  ierr = check_nml_error(io, 'cm_snow_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=cm_snow_nml)
  endif

  delta_time = time_type_to_real(lnd%dt_fast) ! [s]
end subroutine

! ============================================================================
function cm_snow_tile_ctor() result(ptr)
  type(cm_snow_tile_type), pointer :: ptr ! return value
  ! class(snow_tile_type), pointer :: ptr ! return value

  allocate(ptr)
  allocate(ptr%ws(num_l))
  allocate(ptr%wl(num_l))
  allocate(ptr%T(num_l))
  allocate(ptr%e(num_l))
  allocate(ptr%f(num_l))
  ptr%nlayers = num_l

end function cm_snow_tile_ctor

! ============================================================================
function cm_snow_tile_copy_ctor(snow) result(ptr)
  type(cm_snow_tile_type), pointer :: ptr ! return value
  type(cm_snow_tile_type), intent(in) :: snow ! tile to copy

  allocate(ptr)
  ! copy all non-pointer members
  ptr = snow
  ! no need to allocate storage for allocatable components of the type, because
  ! F2003 takes care of that, and also takes care of copying data
end function cm_snow_tile_copy_ctor

! ============================================================================
subroutine cm_delete_snow_tile(snow)
  ! type(cm_snow_tile_type), pointer :: snow
  class(cm_snow_tile_type), pointer :: snow

  ! no need to deallocate components of tile, because F2003 takes care of
  ! allocatable components deallocation when tile is deallocated
  deallocate(snow)
end subroutine cm_delete_snow_tile

! =============================================================================
function cm_snow_tiles_can_be_merged(snow1,snow2) result(response)
  logical :: response
  ! type(cm_snow_tile_type), intent(in) :: snow1,snow2
  class(cm_snow_tile_type), intent(in) :: snow1,snow2

  response = .TRUE.
end function cm_snow_tiles_can_be_merged


! =============================================================================
subroutine cm_merge_snow_tiles_wrapper(snow2, w2, snow1, w1)
  class(snow_tile_type), intent(in)    :: snow1
  class(cm_snow_tile_type), intent(inout) :: snow2
  real                , intent(in)    :: w1, w2 ! relative weights
  ! write(*,*) "CM_SNOW MERGING TILES WRAPPER"
  select type (snow1)
  type is (cm_snow_tile_type)
      call cm_merge_snow_tiles(snow2, w2, snow1, w1)
  class default
   call land_error_message( &
        'cm_merge_snow_tiles_wrapper in cm_snow_tile_mod: type is incorrect!', FATAL)
  end select
end subroutine cm_merge_snow_tiles_wrapper

! =============================================================================
subroutine cm_merge_snow_tiles(snow2, w2, snow1, w1)
  ! type(cm_snow_tile_type), intent(in)    :: snow1
  ! type(cm_snow_tile_type), intent(inout) :: snow2
  class(cm_snow_tile_type), intent(in)    :: snow1
  class(cm_snow_tile_type), intent(inout) :: snow2
  real                , intent(in)    :: w1, w2 ! relative weights

  ! ---- local vars
  real    :: x1, x2 ! normalized weights
  real    :: HEAT1, HEAT2
  integer :: i

  ! write(*,*) "CM_SNOW MERGING TILES..."
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1-x1

  do i = 1, num_l
    HEAT1 = (mc_fict*dz(i)+clw*snow1%wl(i)+csw*snow1%ws(i))*(snow1%T(i)-tfreeze)
    HEAT2 = (mc_fict*dz(i)+clw*snow2%wl(i)+csw*snow2%ws(i))*(snow2%T(i)-tfreeze)
    snow2%wl(i) = snow1%wl(i)*x1 + snow2%wl(i)*x2
    snow2%ws(i) = snow1%ws(i)*x1 + snow2%ws(i)*x2
    if (snow2%wl(i)/=0.or.snow2%ws(i)/=0) then
       snow2%T(i)  = (HEAT1*x1+HEAT2*x2)/&
            (mc_fict*dz(i)+clw*snow2%wl(i)+csw*snow2%ws(i))+tfreeze
    else
       snow2%T(i)  = snow1%T(i)*x1 + snow2%T(i)*x2
    endif
  enddo
end subroutine cm_merge_snow_tiles

! =============================================================================
! returns true if tile fits the specified selector
function cm_snow_is_selected(snow, sel)
  logical cm_snow_is_selected
  type(tile_selector_type),  intent(in) :: sel
  ! type(cm_snow_tile_type),      intent(in) :: snow
  class(cm_snow_tile_type),      intent(in) :: snow

  cm_snow_is_selected = .TRUE.
end function cm_snow_is_selected

! ============================================================================
! retruns tag of the tile
function cm_get_snow_tile_tag(snow) result(tag)
  integer :: tag
  ! type(cm_snow_tile_type), intent(in) :: snow
  class(cm_snow_tile_type), intent(in) :: snow

  tag = snow%tag
end function cm_get_snow_tile_tag

integer function cm_snow_nlayers(snow)
  class(cm_snow_tile_type), intent(in) :: snow
  cm_snow_nlayers = snow%nlayers
end function

! ============================================================================
subroutine cm_snow_roughness(snow, snow_z0s, snow_z0m)
  ! type(cm_snow_tile_type), intent(in) :: snow ! not used
  class(cm_snow_tile_type), intent(in) :: snow ! not used
  real, intent(out):: snow_z0s, snow_z0m

  snow_z0m =  z0_momentum
  snow_z0s =  z0_momentum * exp(-k_over_B)
end subroutine cm_snow_roughness

! ============================================================================
subroutine cm_snow_tile_stock_pe (snow, twd_liq, twd_sol  )
  ! type(cm_snow_tile_type),  intent(in)    :: snow
  class(cm_snow_tile_type),  intent(in)    :: snow
  real,                  intent(out)   :: twd_liq, twd_sol
  integer n

  twd_liq = 0.
  twd_sol = 0.
  do n=1, size(snow%wl)
    twd_liq = twd_liq + snow%wl(n)
    twd_sol = twd_sol + snow%ws(n)
    enddo

end subroutine cm_snow_tile_stock_pe

! returns snow radiative properties: short-wave refletances (by spectral band),
! long-wave reflecatanc, emissivity
subroutine cm_snow_rad_prop (snow, cosz, subs_refl_dif, p_atm, on_glacier, &
                             snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis)
  class(cm_snow_tile_type), intent(inout) :: snow
  real, intent(in) :: cosz
  real, intent(in) :: subs_refl_dif(:) ! slm: not used?
  real, intent(in) :: p_atm            ! not used
  logical, intent(in) :: on_glacier
  real, intent(out) :: snow_refl_dir(:), snow_refl_dif(:)
  real, intent(out) :: snow_refl_lw, snow_emis

  real :: snow_top_temp

  if (snow%snow_active()) then
      snow_top_temp = snow%sfc_temp()
  else
      snow_top_temp = TFREEZE ! NOT used in this case
  endif
  call snow_radiation ( snow_top_temp, cosz, on_glacier, &
      snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis)
end subroutine

! ============================================================================
! returns snow heat content, J/m2
function cm_snow_tile_heat (snow) result(heat) ; real heat
  class(cm_snow_tile_type), intent(in)  :: snow

  integer :: i
  heat = 0
  do i = 1,num_l
     heat = heat - snow%ws(i)*hlf &
        + (mc_fict*dz(i) + clw*snow%wl(i) + csw*snow%ws(i))  &
                                      * (snow%T(i)-tfreeze)
  enddo
end function cm_snow_tile_heat

! ============================================================================
! returns true if snow plays a role
function cm_snow_active(snow) ; logical cm_snow_active
  ! type(cm_snow_tile_type), intent(in)  :: snow
  class(cm_snow_tile_type), intent(in)  :: snow
  cm_snow_active = ( sum(snow%ws(1:num_l)) > 0 )
end function cm_snow_active

! ============================================================================
real function cm_snow_get_sfc_temp(snow)
  ! type(cm_snow_tile_type), intent(in) :: snow
  class(cm_snow_tile_type), intent(in) :: snow

  cm_snow_get_sfc_temp = snow%T(1)
end function

real function cm_snow_get_Ti(snow, i) result(res)
  class(cm_snow_tile_type), intent(in) :: snow
  integer, intent(in) :: i
  res = snow%T(i)
end function

real function cm_snow_get_wsi(snow, i) result(res)
  class(cm_snow_tile_type), intent(in) :: snow
  integer, intent(in) :: i
  res = snow%ws(i)
end function

real function cm_snow_get_wli(snow, i) result(res)
  class(cm_snow_tile_type), intent(in) :: snow
  integer, intent(in) :: i
  res = snow%wl(i)
end function

subroutine cm_snow_set_Ti(snow, i, v)
  class(cm_snow_tile_type), intent(inout) :: snow
  integer, intent(in) :: i
  real, intent(in) :: v
  snow%T(i) = v
end subroutine

subroutine cm_snow_set_wsi(snow, i, v)
  class(cm_snow_tile_type), intent(inout) :: snow
  integer, intent(in) :: i
  real, intent(in) :: v
  snow%ws(i) = v
end subroutine

subroutine cm_snow_set_wli(snow, i, v)
  class(cm_snow_tile_type), intent(inout) :: snow
  integer, intent(in) :: i
  real, intent(in) :: v
  snow%wl(i) = v
end subroutine


real function cm_snow_get_total_ice(snow) result(ice)
  class(cm_snow_tile_type), intent(in) :: snow
  ice = sum(snow%ws(:))
end function


real function cm_snow_get_total_liq(snow) result(liq)
  class(cm_snow_tile_type), intent(in) :: snow
  liq = sum(snow%wl(:))
end function


subroutine cm_snow_get_depth_area(snow, snow_depth, snow_area)
  class(cm_snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_depth, snow_area

  integer :: l

  snow_depth= 0.0
  do l = 1, num_l
     snow_depth = snow_depth + snow%ws(l)
  enddo
  snow_depth = snow_depth / snow_density
  call snow_data_area (snow_depth, snow_area )
end subroutine

subroutine cm_snow_lai(snow, tracers)
  class(cm_snow_tile_type), intent(in) :: snow
  real,                     intent(out) :: tracers(:)
  tracers(:) = 0.0
end subroutine

! ============================================================================
! if snow amount is below specified limit, sweeps it into runoff
subroutine cm_sweep_tiny_snow(snow, lrunf, frunf, hlrunf, hfrunf, lost_wc_em, lost_wc_im)
  class(cm_snow_tile_type), intent(inout) :: snow
  real, intent(out) :: lrunf, frunf   ! liquid and solid runoff generated by sweeper
  real, intent(out) :: hlrunf, hfrunf ! heat carried by liquid ans solid runoff
  real, intent(out) :: lost_wc_em(:), lost_wc_im(:) ! tracer losses

  real :: snow_mass
  integer :: l

  lost_wc_em(:)=0; lost_wc_im(:) = 0

  lrunf=0 ; frunf=0 ; hlrunf=0 ; hfrunf=0
  if (.not.prevent_tiny_snow) return ! do nothing, return zeros

  snow_mass  = sum(snow%ws)
  ! check if the snow is small enough to warrant sweeping
  if ( snow_mass<0 .or. snow_mass >= min_snow_mass ) return

  lrunf  = sum(snow%wl) ; frunf  = snow_mass
  hlrunf = 0.0 ; hfrunf = 0.0
  do l = 1, num_l
     hlrunf = hlrunf + clw*snow%wl(l)*(snow%T(l)-tfreeze)
     hfrunf = hfrunf + csw*snow%ws(l)*(snow%T(l)-tfreeze)
  enddo
  snow%ws = 0 ; snow%wl = 0
end subroutine cm_sweep_tiny_snow

subroutine cm_partition_sw( &
   snow, fswg, fswg_dir, fswg_dif, & ! input
   fswg_substrate, fswg_surface) ! output
   !
   ! Given the shortwave radiation absorbed by snow + substrate (fswg) [W/m2]
   ! as well its direct and diffuse components (fswg_dir, fswg_dif)
   ! partition it between surface of snow (where it was absorbed entirely in old cm snow model)
   ! and, if requested, absorption within the snowpack
   ! andabsoirption in the underlying substrate (lake/soil/glacier)
   !
   class(cm_snow_tile_type), intent(inout) :: snow !< state of snowpack
   real, intent(in)  :: fswg ! total sw absorbed by snow + substrate [W/m2]
   real, intent(in)  :: fswg_dir(:), fswg_dif(:) ! total sw absorbed by snow + substrate (dir only, dif only) [W/m2]
   real, intent(out) :: fswg_substrate ! sw radiation passed to substrate [W/m2]
   real, intent(out) :: fswg_surface   ! sw radiation to be absorbed at the surface [W/m2]

   ! case of CM snow model: all sw absorption occurs at the surface (part of surface energy balance)
   fswg_surface=fswg
   fswg_substrate = 0.0

   if (is_watch_point()) then
      write(*,*) "##### cm_partition_sw checkpoint 1: #####"
      __DEBUG1__(fswg)
      __DEBUG1__(fswg_surface)
      __DEBUG1__(fswg_substrate)
   endif
end subroutine cm_partition_sw

! ============================================================================
! update snow properties explicitly for time step.
! integrate snow-heat conduction equation upward from bottom of snow
! to surface, delivering linearization of surface ground heat flux.
subroutine cm_snow_step_1 ( snow, p_surf, grnd_T, snow_G_Z, snow_G_TZ, &
                         snow_active, snow_T, snow_rh, snow_liq, snow_ice, &
                         snow_subl, snow_area, snow_G0, snow_DGDT, snow_E_max )
  class(cm_snow_tile_type), intent(inout) :: snow
  real,    intent(in) :: p_surf
  real,    intent(in) :: grnd_T
  real,    intent(in) :: snow_G_Z
  real,    intent(in) :: snow_G_TZ
  logical, intent(out):: snow_active
  real,    intent(out):: &
       snow_T, snow_rh, snow_liq, snow_ice, &
       snow_subl, snow_area, snow_G0, snow_DGDT, &
       snow_E_max

  ! ---- local vars
  real :: snow_depth, bbb, denom, dt_e
  real, dimension(num_l):: aaa, ccc, thermal_cond, dz_phys, heat_capacity
  integer :: l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  snow_T = tfreeze
  snow_T = snow%T(1)

  call snow_data_thermodynamics ( snow_rh, thermal_cond )
  snow_depth= 0.0
  do l = 1, num_l
     snow_depth = snow_depth + snow%ws(l)
  enddo
  snow_depth = snow_depth / snow_density
  call snow_data_area (snow_depth, snow_area )
  ! ---- only liquid in the top snow layer is available to freeze implicitly
  snow_liq =     snow%wl(1)
  ! ---- snow in any layer can be melted implicitly
  snow_ice = sum(snow%ws(:))

! ---- fractionate evaporation/sublimation according to sfc phase ratios
!  where (max(snow%ws(1),0.)+max(snow%wl(1),0.)>0)
!      snow_subl = max(snow%ws(1),0.) &
!       /(max(snow%ws(1),0.)+max(snow%wl(1),0.))
!    elsewhere
!      snow_subl = 0
!    endwhere
!  snow_active = snow_subl>0.
  if (snow_depth>0) then
     snow_subl = 1.
  else
     snow_subl = 0
  endif
  snow_active = snow_subl>0.

  do l = 1, num_l
     dz_phys(l) = dz(l)*snow_depth
  enddo

  if (retro_heat_capacity) then
     do l = 1, num_l
        heat_capacity(l) = heat_capacity_retro*dz_phys(l)
     enddo
  else
     do l = 1, num_l
        heat_capacity(l) = mc_fict*dz(l) + &
             clw*snow%wl(l) + csw*snow%ws(l)
     enddo
  endif

!  if(num_l > 1) then
  if (snow_depth > 0) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz_phys(l+1)/thermal_cond(l+1) &
                     + dz_phys(l)/thermal_cond(l)   )
        aaa(l+1) = - dt_e * delta_time / heat_capacity(l+1)
        ccc(l)   = - dt_e * delta_time / heat_capacity(l)
     enddo

     bbb = 1.0 - aaa(num_l) + delta_time*snow_G_TZ/heat_capacity(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(snow%T(num_l) - snow%T(num_l-1)) &
          - delta_time*snow_G_Z/heat_capacity(num_l)
     snow%e(num_l-1) = -aaa(num_l)/denom
     snow%f(num_l-1) = dt_e/denom

     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*snow%e(l)
        dt_e = - ( ccc(l)*(snow%T(l+1) - snow%T(l)  ) &
                  -aaa(l)*(snow%T(l)   - snow%T(l-1)) )
        snow%e(l-1) = -aaa(l)/denom
        snow%f(l-1) = (dt_e - ccc(l)*snow%f(l))/denom
     enddo

     denom = delta_time/heat_capacity(1)
     snow_G0    = ccc(1)*(snow%T(2)- snow%T(1) &
          + snow%f(1)) / denom
     snow_DGDT  = (1 - ccc(1)*(1-snow%e(1))) / denom
  endif

!    else  ! one-level case
!      denom = delta_time/heat_capacity(1)
!      snow_G0    = 0.
!      snow_DGDT  = 1. / denom
!    end if

  if (snow_depth <= 0) then
     snow_G0   = snow_G_Z
     snow_DGDT = snow_G_TZ
  endif

  snow_E_max = HUGE(1.0)

  if(is_watch_point()) then
     write(*,*) 'snow_depth', snow_depth
     write(*,*) '############ snow_step_1 output'
     write(*,*) 'mask      ', .true.
     write(*,*) 'snow_T    ', snow_T
     write(*,*) 'snow_rh   ', snow_rh
     write(*,*) 'snow_liq  ', snow_liq
     write(*,*) 'snow_ice  ', snow_ice
     write(*,*) 'snow_subl ', snow_subl
     write(*,*) 'snow_area ', snow_area
     write(*,*) 'snow_G_Z  ', snow_G_Z
     write(*,*) 'snow_G_TZ ', snow_G_TZ
     write(*,*) 'snow_G0   ', snow_G0
     write(*,*) 'snow_DGDT ', snow_DGDT
     write(*,*) '############ end of snow_step_1 output'
  endif

end subroutine cm_snow_step_1

subroutine cm_snow_step_2 ( snow, snow_subl,                     &
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
    class(cm_snow_tile_type), intent(inout) :: snow
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

  ! ---- local vars
  real, dimension(num_l) :: del_t, M_layer
  real :: depth, &
         cap0, dW_l, dW_s, dcap, dheat,&
         melt, melt_per_deg, drain,       &
         snow_mass, sum_liq, &
         sum_heat, sum_sno, &
         snow_transfer, frac,&
         liq_rate, hliq_rate,&
         sno_rate, hsno_rate, fict_heat, &
         evapg_lm2, vegn_fprec_lm2, &
         snow_LMASS, snow_FMASS, snow_HEAT
  integer :: l, l_old
  real :: new_ws(num_l)
  real :: new_wl(num_l)
  real :: new_T(num_l)
  ! --------------------------------------------------------------------------


  depth= 0.
  do l = 1, num_l
    depth = depth + snow%ws(l)
  enddo
  depth = depth / snow_density

  if(is_watch_point()) then
     write(*,*) '############ snow_step_2 input'
     write(*,*) 'mask       ', .TRUE.
     write(*,*) 'snow_subl  ', snow_subl
     write(*,*) 'vegn_lprec ', vegn_lprec
     write(*,*) 'vegn_fprec ', vegn_fprec
     write(*,*) 'vegn_hlprec', vegn_hlprec
     write(*,*) 'vegn_hfprec', vegn_hfprec
     write(*,*) 'DTg        ', DTg
     write(*,*) 'Mg_imp     ', Mg_imp
     write(*,*) 'evapg      ', evapg
     write(*,*) 'fswg       ', fswg
     write(*,*) 'flwg       ', flwg
     write(*,*) 'sensg      ', sensg
     write(*,*) '############ end of snow_step_2 input'

     write(*,*) 'depth   ', depth
     do l = 1, num_l
        write(*,'(i2,3(x,a,g23.16))') l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             ' T =', snow%T(l)
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l;
        snow_LMASS = snow_LMASS + snow%wl(l)
        snow_FMASS = snow_FMASS + snow%ws(l)
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
                                                * (snow%T(l)-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 1.01 ***** '
     call cm_print_snow_integrals(snow)
  endif

  ! ---- record fluxes -------------------------------------------------------
  call cm_get_snow_integrals(snow, snow_LMASS, snow_FMASS, snow_HEAT)
  if (lm2.and.steal) then
    if (snow_FMASS-Mg_imp > 0.) then
      if (evapg <= (snow_FMASS-Mg_imp)/delta_time) then
          evapg_lm2 = evapg
          vegn_fprec_lm2 = vegn_fprec
      else if (evapg <= (snow_FMASS-Mg_imp)/delta_time+vegn_fprec) then
          evapg_lm2 = evapg
          vegn_fprec_lm2 = vegn_fprec - evapg + (snow_FMASS-Mg_imp)/delta_time
      else
          evapg_lm2 = (snow_FMASS-Mg_imp)/delta_time+vegn_fprec
          vegn_fprec_lm2 = 0.
      endif
    else
      evapg_lm2 = 0.
      vegn_fprec_lm2 = vegn_fprec
    endif
  else
     evapg_lm2 = evapg
     vegn_fprec_lm2 = vegn_fprec
  endif
  vegn_fprec_lm2 = vegn_fprec
  if (depth>0) then
        snow_fsw   = fswg
        snow_flw   = flwg
        snow_sens  = sensg
        snow_levap = evapg_lm2*(1-snow_subl)
        snow_fevap = evapg_lm2*   snow_subl
  else
        snow_fsw    = 0
        snow_flw    = 0
        snow_sens   = 0
        snow_levap  = 0
        snow_fevap  = 0
  endif
  subs_fsw = fswg - snow_fsw
  subs_flw = flwg - snow_flw
  subs_evap = evapg - snow_levap - snow_fevap
  subs_sens = sensg - snow_sens

  ! ---- load surface temp change and perform back substitution --------------
  if (depth>0) then
      del_t(1) = DTg
      snow%T(1)  = snow%T(1) + del_t(1)
  endif
  if ( num_l > 1) then
     do l = 1, num_l-1
        if (depth>0) then
            del_t(l+1) = snow%e(l) * del_t(l) + snow%f(l)
            snow%T(l+1) = snow%T(l+1) + del_t(l+1)
        endif
     enddo
  endif
  if (depth>0) then
    subs_DT = del_t(num_l)
  else
    subs_DT = DTg
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2 ***** '
     do l = 1, num_l
        write(*,'(i2,a,g23.16)') l,' T =', snow%T(l)
     enddo
     call cm_print_snow_integrals(snow)
  endif

  ! ---- evaporation and sublimation -----------------------------------------
  if (depth>0) then
      snow%wl(1) = snow%wl(1) - snow_levap*delta_time
      snow%ws(1) = snow%ws(1) - snow_fevap*delta_time
      cap0 = mc_fict*dz(1) + clw*snow%wl(1) + csw*snow%ws(1)
      ! T adjustment for nonlinear terms (del_T)*(del_W)
      dheat = delta_time*(clw*snow_levap+csw*snow_fevap)*del_T(1)
      ! take out extra heat not claimed in advance for evaporation
      if (use_tfreeze_in_grnd_latent) dheat = dheat &
            - delta_time*((cpw-clw)*snow_levap+(cpw-csw)*snow_fevap) &
                               *(snow%T(1)-del_T(1)-tfreeze)
      snow%T(1)  = snow%T(1)  + dheat/cap0
    endif

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2.5 ***** '
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))')l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             ' T =', snow%T(l)
     enddo
     call cm_print_snow_integrals(snow)
  endif

  ! ---- distribute implicit phase change downward through snow layers -------
  if (depth>0) then
      snow_melt = Mg_imp/delta_time
  else
      snow_melt = 0
  endif
  M_layer = 0.
  subs_M_imp = Mg_imp
  do l = 1, num_l
    if (depth>0 .and. subs_M_imp.gt.0) then
        M_layer(l) =  min( subs_M_imp, max(0.,snow%ws(l)) )
        subs_M_imp = subs_M_imp - M_layer(l)
    endif
  enddo
  if (depth>0) then
      M_layer(1) = M_layer(1) + subs_M_imp
      subs_M_imp = 0.
  endif
  do l = 1, num_l
    if (depth>0) then
          cap0 = mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l)
          snow%wl(l) = snow%wl(l) + M_layer(l)
          snow%ws(l) = snow%ws(l) - M_layer(l)
          snow%T(l)  = tfreeze + (cap0*(snow%T(l)-tfreeze) ) &
                                                          / ( cap0 + (clw-csw)*M_layer(l) )
    endif
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 3 ***** '
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))') l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             '  T=', snow%T(l)
     enddo
     call cm_print_snow_integrals(snow)
  endif

! ----------------------------------------------------------------------------
!   call snow_data_hydraulics (pars, snow%wl, psi, hyd_cond )

! ---- remainder of mass fluxes and associated sensible heat fluxes ----------
  liq_rate = vegn_lprec
  sno_rate = vegn_fprec_lm2
  hliq_rate = vegn_hlprec
  if (vegn_fprec.ne.0.) then
          hsno_rate = vegn_hfprec*(vegn_fprec_lm2/vegn_fprec)
  else
          hsno_rate = 0.
  endif

  do l = 1, num_l
    if(depth>0 .or. vegn_fprec_lm2>0) then
    ! ---- mix inflow with existing snow and water ---------------------------
          cap0 = mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l)
          dW_l = liq_rate*delta_time
          dW_s = sno_rate*delta_time
          dcap = clw*dW_l + csw*dW_s
          snow%ws(l) = snow%ws(l) + dW_s
          snow%wl(l) = snow%wl(l) + dW_l
          snow%T(l)  = tfreeze + (cap0*(snow%T(l)-tfreeze) &
                               + (hsno_rate+hliq_rate)*delta_time) /(cap0 + dcap)
    endif

    if(is_watch_point()) then
       write(*,*) ' ***** snow_step_2 checkpoint 4a ***** '
       write(*,'(i2,3(a,g23.16))') l,&
            ' wl=', snow%wl(l),&
            ' ws=', snow%ws(l),&
            '  T=', snow%T(l)
    endif

    if (depth>0 .or. vegn_fprec_lm2>0) then
    ! ---- compute explicit melt/freeze --------------------------------------
          melt_per_deg = (cap0+dcap)/hlf
          if (snow%ws(l)>0 .and. snow%T(l)>tfreeze) then
                  melt =  min(snow%ws(l), (snow%T(l)-tfreeze)*melt_per_deg)
      elseif (snow%wl(l)>0 .and. snow%T(l)<tfreeze) then
                  melt = -min(snow%wl(l), (tfreeze-snow%T(l))*melt_per_deg)
      else
                  melt = 0
          endif
          snow_melt = snow_melt + melt/delta_time
          snow%wl(l) = snow%wl(l) + melt
          snow%ws(l) = snow%ws(l) - melt
!        where (cap0+dcap.ne.0.) &
!        snow%T(l)  = snow%T(l)  - melt/melt_per_deg
          snow%T(l) = tfreeze &
                 + ((cap0+dcap)*(snow%T(l)-tfreeze) - hlf*melt) &
                                                          / ( cap0+dcap + (clw-csw)*melt )
    endif

   if(is_watch_point()) then
      write(*,*) ' ***** snow_step_2 checkpoint 4b ***** '
      write(*,'(i2,3(a,g23.16))')l,&
            ' wl=', snow%wl(l),&
            ' ws=', snow%ws(l),&
            '  T=', snow%T(l)
   endif

   if (depth>0 .or. vegn_fprec_lm2>0) then
    ! ---- compute drainage from this layer to next --------------------------
        drain = max (0., snow%wl(l) - wet_max*snow%ws(l))
        snow%wl(l) = snow%wl(l) - drain
        liq_rate = drain / delta_time
        hliq_rate = clw*liq_rate*(snow%T(l)-tfreeze)
        sno_rate = 0
        hsno_rate = 0
    endif
  enddo

  snow_lprec  = liq_rate
  snow_hlprec = hliq_rate

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 4c ***** '
     call cm_print_snow_integrals(snow)
  endif

! ---- conceptually remove fictitious mass/heat for the moment ---------------
  fict_heat = 0.
  do l = 1, num_l
    fict_heat = fict_heat + dz(l)*snow%T(l)     ! (*mc_fict)
  enddo

  snow_mass  = sum(snow%ws)
  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 4d ***** '
     write(*,*) 'max_snow    ', max_snow
     write(*,*) 'snow_mass   ', snow_mass
  endif

! ---- remove any isolated snow molecules (!) or sweep any excess snow from top of pack ----

  snow_lrunf  = 0.
  snow_frunf  = 0.
  snow_hlrunf = 0.
  snow_hfrunf = 0.
  if (0. < snow_mass .and. snow_mass < min_snow_mass ) then
     do l = 1, num_l
       snow_hlrunf = snow_hlrunf  &
         + clw*snow%wl(l)*(snow%T(l)-tfreeze)
       snow_hfrunf = snow_hfrunf  &
         + csw*snow%ws(l)*(snow%T(l)-tfreeze)
       enddo
     snow_lrunf  = sum(snow%wl)
     snow_frunf  = snow_mass
     snow_mass   = 0.
     snow%ws = 0.
     snow%wl = 0.
  else if (max_snow < snow_mass) then
     snow_frunf  = snow_mass - max_snow
     snow_mass  = max_snow
     sum_sno  = 0
     snow_transfer = 0
     do l = 1, num_l
        if (sum_sno + snow%ws(l) > snow_frunf) then
           snow_transfer = snow_frunf - sum_sno
        else
           snow_transfer = snow%ws(l)
        endif
        if (snow%ws(l) > 0) then
           frac = snow_transfer / snow%ws(l)
        else
           frac = 1.
        endif
        sum_sno  = sum_sno  + snow_transfer
        snow_lrunf  = snow_lrunf  +     frac*snow%wl(l)
        snow_hlrunf = snow_hlrunf + clw*frac*snow%wl(l)*(snow%T(l)-tfreeze)
        snow_hfrunf = snow_hfrunf + csw*frac*snow%ws(l)*(snow%T(l)-tfreeze)
        snow%ws(l) = (1-frac)*snow%ws(l)
        snow%wl(l) = (1-frac)*snow%wl(l)
     enddo
  endif
  snow_lrunf  = snow_lrunf  / delta_time
  snow_frunf  = snow_frunf  / delta_time
  snow_hlrunf = snow_hlrunf / delta_time
  snow_hfrunf = snow_hfrunf / delta_time

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 5 ***** '
     write(*,*) 'fict_heat         ', fict_heat
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))')l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             ' T =', snow%T(l)
     enddo
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 5.01 ***** '
     call cm_print_snow_integrals(snow)
  endif

  depth= 0.
  new_ws=0
  new_wl=0
  new_T=0
  do l = 1, num_l
     depth = depth + snow%ws(l)
  enddo
  depth = depth / snow_density

!************************** fudge to avoid T=NaN from too-small mass **
!   if(depth*snow_density < min_snow_mass .and. depth>0.) then
!       depth = 0
!       snow%ws = 0
!       snow%wl = 0
!     endif

! ---- re-layer the snowpack ------------------------------------------------
  do l = 1, num_l
     if (depth > 0) then
        new_ws(l) = snow_mass*dz(l)
        sum_sno = 0
        sum_liq = 0
        sum_heat = 0
     endif
     do l_old = 1, num_l
        if (depth > 0) then
           if (sum_sno + snow%ws(l_old) > new_ws(l)) then
              snow_transfer = new_ws(l) - sum_sno
           else
              snow_transfer = snow%ws(l_old)
           endif
           if (snow%ws(l_old) .ne. 0.) then
              frac = snow_transfer / snow%ws(l_old)
           else
              frac = 1
           endif
           sum_sno  = sum_sno  + snow_transfer
           sum_liq  = sum_liq  + frac*     snow%wl(l_old)
           sum_heat = sum_heat + frac*&
                (clw*snow%wl(l_old) + csw*snow%ws(l_old))&
                *snow%T(l_old)
           snow%ws(l_old) = (1.-frac)*snow%ws(l_old)
           snow%wl(l_old) = (1.-frac)*snow%wl(l_old)
           if(is_watch_point()) then
              write(*,'(i2,2x,a,i2,99(2x,a,g23.16))')l,&
                  'l_old=',l_old, 'snow_transfer=',snow_transfer,'frac=',frac,&
                   'sum_sno=',sum_sno,'sum_liq=',sum_liq,'sum_heat=',sum_heat
           endif
        endif

     enddo
     if (depth > 0) then
        new_wl(l) = sum_liq
        new_T(l)  = sum_heat / (clw*new_wl(l) + csw*new_ws(l))
     endif
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 5.1 ***** '
     write(*,*) 'depth             ', depth
     write(*,*) 'fict_heat         ', fict_heat
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))')l,&
             ' new_wl=', new_wl(l),&
             ' new_ws=', new_ws(l),&
             ' new_T =', new_T(l)
     enddo
  endif

! add back fictional mass/heat
  do l = 1, num_l
    if (depth > 0) &
    new_T(l) = ( &
    (clw*new_wl(l) + csw*new_ws(l))*new_T(l)  &
      + mc_fict*dz(l)*fict_heat ) &
      / (clw*new_wl(l) + csw*new_ws(l) + dz(l)*mc_fict)
  enddo

  do l = 1, num_l
    if (depth > 0) then
      snow%ws(l) = new_ws(l)
      snow%wl(l) = new_wl(l)
      snow%T(l)  = new_T(l)
    endif
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 6 ***** '
     write(*,*) 'evap         ', subs_evap
     write(*,*) 'snow_lprec', snow_lprec
     write(*,*) 'depth        ', depth
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))')l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             ' T =', snow%T(l)
     enddo
  endif

  call cm_get_snow_integrals(snow, snow_LMASS, snow_FMASS, snow_HEAT)
  snow_Tbot = snow%T(num_l)
  snow_Cbot = mc_fict*dz(num_l) &
        + clw*snow%wl(num_l) + csw*snow%ws(num_l)
  snow_C = sum(mc_fict*dz(1:num_l) &
        + clw*snow%wl(1:num_l) + csw*snow%ws(1:num_l))
  snow_avrg_T = snow_HEAT/snow_C+tfreeze

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 7 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif
end subroutine cm_snow_step_2

! ============================================================================
subroutine cm_get_snow_integrals(snow, snow_LMASS, snow_FMASS, snow_HEAT)
  type(cm_snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_LMASS, snow_FMASS, snow_HEAT

  integer :: l

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l;
    snow_LMASS = snow_LMASS + snow%wl(l)
    snow_FMASS = snow_FMASS + snow%ws(l)
    snow_HEAT = snow_HEAT + &
      (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
                                            * (snow%T(l)-tfreeze)
  enddo
end subroutine cm_get_snow_integrals

! ============================================================================
subroutine cm_print_snow_integrals(snow)
  type(cm_snow_tile_type), intent(in) :: snow

  real    :: snow_LMASS, snow_FMASS, snow_HEAT
  call cm_get_snow_integrals(snow, snow_LMASS, snow_FMASS, snow_HEAT)
  __DEBUG3__(snow_LMASS, snow_FMASS, snow_HEAT)
end subroutine cm_print_snow_integrals

end module cm_snow_tile_mod
