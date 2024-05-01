module cm_snow_tile_mod
#include <fms_platform.h>
#include "../../shared/debug.inc"


use fms_mod, only : input_nml_file, check_nml_error, stdlog, mpp_pe, mpp_root_pe, FATAL
use constants_mod,only: tfreeze, hlf
use land_tile_selectors_mod, only : tile_selector_type
use land_data_mod, only : log_version
use snow_tile_mod, only : snow_tile_type, mc_fict, z0_momentum, k_over_B, num_l, dz
use land_debug_mod, only : is_watch_point, land_error_message
use snowpack_mod, only : cpw, clw, csw

implicit none
private

! ==== public interfaces =====================================================
public :: cm_snow_tile_type
public :: cm_snow_tile_ctor
public :: cm_snow_tile_copy_ctor
public :: cm_delete_snow_tile
public :: cm_snow_tiles_can_be_merged
public :: cm_merge_snow_tiles
public :: cm_snow_is_selected
public :: cm_get_snow_tile_tag
public :: cm_snow_tile_stock_pe
public :: cm_snow_tile_heat
public :: cm_snow_active
public :: cm_snow_roughness
public :: cm_snow_get_sfc_temp

public :: read_snow_cm_namelist
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'cm_snow_tile_mod'
#include "../../shared/version_variable.inc"


! ==== types =================================================================

type, extends(snow_tile_type) :: cm_snow_tile_type
   ! data structure already defined in parent snow type
   contains
    procedure :: merge_snow_tiles => cm_merge_snow_tiles_wrapper
    procedure :: get_snow_tile_tag => cm_get_snow_tile_tag
    procedure :: snow_is_selected => cm_snow_is_selected
    procedure :: snow_roughness => cm_snow_roughness
    procedure :: stock_pe => cm_snow_tile_stock_pe
    procedure :: snow_active => cm_snow_active
    procedure :: snow_tile_heat => cm_snow_tile_heat
    procedure :: snow_get_sfc_temp => cm_snow_get_sfc_temp
    procedure :: get_Ti => cm_snow_get_Ti
    procedure :: get_wli => cm_snow_get_wli
    procedure :: get_wsi => cm_snow_get_wsi
    procedure :: set_Ti => cm_snow_set_Ti
    procedure :: set_wli => cm_snow_set_wli
    procedure :: set_wsi => cm_snow_set_wsi
    procedure :: ice => cm_snow_get_total_ice
    procedure :: liq => cm_snow_get_total_liq

    procedure :: sweep_tiny => cm_sweep_tiny_snow
    procedure :: partition_sw => cm_partition_sw

end type cm_snow_tile_type

! ==== module data ===========================================================

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

! ============================================================================
! returns snow heat content, J/m2
function cm_snow_tile_heat (snow) result(heat) ; real heat
  ! type(cm_snow_tile_type), intent(in)  :: snow
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
subroutine cm_snow_get_sfc_temp(snow, snow_T)
  ! type(cm_snow_tile_type), intent(in) :: snow
  class(cm_snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_T

  snow_T = snow%T(1)
end subroutine

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

end module cm_snow_tile_mod
