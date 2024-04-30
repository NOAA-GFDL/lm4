module snow_tile_mod
#include <fms_platform.h>

use fms_mod, only : error_mesg, FATAL, NOTE, lowercase

use land_data_mod, only : log_version
use land_debug_mod, only : land_error_message

use cm_snow_tile_mod, only : cm_snow_tile_ctor
use gl_snow_tile_mod, only : gl_snow_tile_ctor
use parent_snow_tile_mod, only: snow_tile_type, snow_option

implicit none
private

! ==== public interfaces =====================================================
public :: new_snow_tile
public :: delete_snow_tile
public :: snow_tiles_can_be_merged
! ==== end of public interfaces ==============================================

interface new_snow_tile
   module procedure snow_tile_ctor
   module procedure snow_tile_copy_ctor
end interface

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snow_tile_mod'
#include "../shared/version_variable.inc"
! ==== types =================================================================

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
function snow_tile_ctor(tag) result(ptr)
  class(snow_tile_type), pointer :: ptr ! return value
  integer, optional, intent(in) :: tag ! kind of tile
  select case(trim(lowercase(snow_option)))
  case('cm')
    ptr => cm_snow_tile_ctor()
  case('gl')
    ptr => gl_snow_tile_ctor()
  case default
    call error_mesg( &
        'snow_tile_ctor in snow_tile_mod', &
        'snow_option = "'//trim(snow_option)//'" is incorrect, use "cm", or "gl"', FATAL)
   end select
end function snow_tile_ctor


function snow_tile_copy_ctor(snow) result(ptr)
  class(snow_tile_type), pointer :: ptr ! return value
  class(snow_tile_type), intent(in) :: snow ! tile to copy

  real liq1, liq2, ice1, ice2, heat1, heat2, dheat, dwat

  allocate(ptr, source=snow)
  ptr%sp = snow%sp
  ptr%sp%snow = snow%sp%snow

  call snow%stock_pe(liq1, ice1)
  call ptr%stock_pe(liq2, ice2)
  heat1 =  snow%snow_tile_heat()
  heat2 = ptr%snow_tile_heat()

  if (abs(liq1-liq2)>1E-2) call land_error_message("snow_tile_copy_ctor in snow_tile_mod: liquid water non conserved!", FATAL)
  if (abs(ice1-ice2)>1E-2) call land_error_message("snow_tile_copy_ctor in snow_tile_mod: frozen water non conserved!", FATAL)
  if (abs(heat1-heat2)>1E-2) call land_error_message("snow_tile_copy_ctor in snow_tile_mod: heat non conserved!", FATAL)
end function snow_tile_copy_ctor


subroutine delete_snow_tile(snow)
  class(snow_tile_type), pointer :: snow
  deallocate(snow)
end subroutine delete_snow_tile


function snow_tiles_can_be_merged(snow1,snow2) result(response)
  logical :: response
  class(snow_tile_type), intent(in) :: snow1,snow2
! select case(trim(snow_option))
! case('cm')
!    response = cm1_snow_tiles_can_be_merged(snow1, snow2)
! case('gl')
!    response = cm2_snow_tiles_can_be_merged(snow1, snow2)
! case default
!    call error_mesg( &
!         'snow_tiles_can_be_merged in snow_tile_mod', &
!         'snow_option = "'//trim(snow_option)//'" is incorrect, use "cm1", "cm2", or "gl"', FATAL)
! end select
! // TODO to make it type-specific if necessary
  response = .TRUE.
end function snow_tiles_can_be_merged


end module snow_tile_mod
