module snow_tile_mod
#include <fms_platform.h>

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : error_mesg, file_exist, check_nml_error, close_file, stdlog, FATAL, NOTE, lowercase
use constants_mod,only: tfreeze, hlf
use land_constants_mod, only : NBANDS
use land_tile_selectors_mod, only : tile_selector_type
use land_data_mod, only : log_version
use land_debug_mod, only : is_watch_point, land_error_message

use cm_snow_tile_mod, only : cm_snow_tile_type, cm_snow_tile_ctor, & 
  cm_snow_tile_copy_ctor, cm_delete_snow_tile, &
  cm_snow_tiles_can_be_merged, cm_merge_snow_tiles, cm_snow_is_selected, &
  cm_get_snow_tile_tag, cm_snow_tile_stock_pe, cm_snow_tile_heat, &
  cm_snow_active, &
  cm_snow_roughness, cm_snow_get_sfc_temp

use gl_snow_tile_mod, only : gl_snow_tile_type, gl_snow_tile_ctor, & 
  gl_snow_tile_copy_ctor, gl_delete_snow_tile, &
  gl_snow_tiles_can_be_merged, gl_merge_snow_tiles, gl_snow_is_selected, &
  gl_get_snow_tile_tag, gl_snow_tile_stock_pe, gl_snow_tile_heat, &
  gl_snow_active, &
  gl_snow_roughness, gl_snow_get_sfc_temp

use parent_snow_tile_mod, only: snow_tile_type, snow_option, num_l


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
