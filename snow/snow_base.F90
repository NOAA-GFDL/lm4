module snow_base_mod
#include <fms_platform.h>

use mpp_mod, only: input_nml_file
use fms_mod, only: check_nml_error, lowercase, &
            stdlog, mpp_pe, mpp_root_pe, error_mesg, FATAL, NOTE

use land_data_mod, only : log_version
use land_debug_mod, only : land_error_message

use snow_tile_mod, only: snow_tile_type
use cm_snow_tile_mod, only : cm_snow_tile_type, new_cm_snow_tile
use gl_snow_tile_mod, only : gl_snow_tile_type, new_gl_snow_tile

implicit none
private

! ==== public interfaces =====================================================
public :: read_snow_model_namelist
public :: new_snow_tile
public :: delete_snow_tile
public :: snow_tiles_can_be_merged
! ==== end of public interfaces ==============================================

interface new_snow_tile
   module procedure snow_tile_ctor
   module procedure snow_tile_copy
end interface

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snow_base_mod'
#include "../shared/version_variable.inc"

!---- namelist ---------------------------------------------------------------
character(32) :: model_to_use = 'cm' ! name of the snow model to use
namelist /snow_nml/ model_to_use

integer, public, protected :: snow_option = -1
integer, public, parameter :: &
    SNOW_CM = 1, &  ! Milly model
    SNOW_GL = 2     ! GLASS snow model, Zorzetto et al. 2024

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine read_snow_model_namelist()
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, __FILE__)
  read (input_nml_file, nml=snow_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=snow_nml)
  endif

  ! parse snow model options
  select case (lowercase(model_to_use))
  case('cm')
    snow_option = SNOW_CM
  case('gl')
    snow_option = SNOW_GL
  case default
    call error_mesg('read_snow_model_namelist', &
        '"'//trim(model_to_use)//'" is an invalid option for model_to_use', FATAL)
  end select
end subroutine

! ============================================================================
function snow_tile_ctor(tag) result(ptr)
  class(snow_tile_type), pointer :: ptr ! return value
  integer, optional, intent(in) :: tag ! kind of tile
  select case(snow_option)
  case(SNOW_CM)
     ptr => new_cm_snow_tile()
  case(SNOW_GL)
     ptr => new_gl_snow_tile()
  case default
     call land_error_message('snow_tile_ctor: The value of snow_option is invalid. This should never happen. See developer', FATAL)
  end select
end function snow_tile_ctor


function snow_tile_copy(snow) result(ptr)
  class(snow_tile_type), pointer    :: ptr  ! return value
  class(snow_tile_type), intent(in) :: snow ! tile to copy

  real liq1, liq2, ice1, ice2, heat1, heat2, dheat, dwat

  select type (snow)
  type is (cm_snow_tile_type)
     ptr => new_cm_snow_tile(snow)
  type is (gl_snow_tile_type)
     ptr => new_gl_snow_tile(snow)
  class default
     call land_error_message('snow_tile_copy: The type of snow tile is invalid. This should never happen. See developer', FATAL)
  end select

  ! slm: why are these checks here? Why tolerances are so high? Copies must
  ! be identical: ani discreapanc y is an error
  liq1 = snow%liq(); ice1 = snow%ice(); heat1 = snow%snow_tile_heat()
  liq2 = ptr%liq();  ice2 = ptr%ice();  heat2 = ptr%snow_tile_heat()

  if (abs(liq1-liq2)>1E-2) call land_error_message("snow_tile_copy_ctor in snow_tile_mod: liquid water non conserved!", FATAL)
  if (abs(ice1-ice2)>1E-2) call land_error_message("snow_tile_copy_ctor in snow_tile_mod: frozen water non conserved!", FATAL)
  if (abs(heat1-heat2)>1E-2) call land_error_message("snow_tile_copy_ctor in snow_tile_mod: heat non conserved!", FATAL)
end function


subroutine delete_snow_tile(snow)
  class(snow_tile_type), pointer :: snow
  deallocate(snow)
end subroutine delete_snow_tile


function snow_tiles_can_be_merged(snow1,snow2) result(response)
  logical :: response
  class(snow_tile_type), intent(in) :: snow1,snow2
!  select case(snow_option)
!  case(SNOW_CM)
!     response = cm1_snow_tiles_can_be_merged(snow1, snow2)
!  case(SNOW_GL)
!     response = cm2_snow_tiles_can_be_merged(snow1, snow2)
!  case default
!     call land_error_message('snow_tiles_can_be_merged: The value of snow_option is invalid. This should never happen. See developer', FATAL)
!  end select
!  // TODO to make it type-specific if necessary
  response = .TRUE.
end function snow_tiles_can_be_merged


end module snow_base_mod
