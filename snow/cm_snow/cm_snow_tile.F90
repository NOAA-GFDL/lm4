module cm_snow_tile_mod
#include <fms_platform.h>

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : file_exist, check_nml_error, close_file, stdlog, error_mesg, FATAL, NOTE
use constants_mod,only: tfreeze, hlf
use land_constants_mod, only : NBANDS
use land_tile_selectors_mod, only : tile_selector_type
use land_data_mod, only : log_version
use parent_snow_tile_mod, only : snow_tile_type, mc_fict, z0_momentum, k_over_B, num_l, dz
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
end type cm_snow_tile_type

! ==== module data ===========================================================

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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


end module cm_snow_tile_mod
