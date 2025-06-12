module soil_BGC_base_mod

use fms_mod, only: check_nml_error, input_nml_file, &
            stdlog, mpp_pe, mpp_root_pe, error_mesg, FATAL
use land_data_mod, only: log_version
use land_debug_mod, only: land_error_message

use soil_tile_mod, only: soil_tile_type
use soil_BGC_type_mod, only: soil_BGC_t
use soil_BGC_SIMPLE_type_mod, only: soil_BGC_SIMPLE_t, new_soilc_SIMPLE, read_soil_BGC_SIMPLE_namelist
use soil_BGC_CORPSE_type_mod, only: soil_BGC_CORPSE_t, new_soilc_CORPSE, read_soil_BGC_CORPSE_namelist
use soil_BGC_GIMICS_type_mod, only: soil_BGC_GIMICS_t, new_soilc_GIMICS, read_soil_BGC_GIMICS_namelist

implicit none; private

public :: read_soil_BGC_namelist
public :: new_soilc, delete_soilc

public :: soil_BGC_option

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'soil_BGC_base_mod'
#include "../shared/version_variable.inc"

! ==== module interfaces ======================================================
interface new_soilc
   module procedure soilc_ctor
   module procedure soilc_copy
end interface

!---- namelist ---------------------------------------------------------------
character(32) :: model_to_use = 'SIMPLE' ! or 'CORPSE'
namelist /soil_BGC_nml/ model_to_use

! soil carbon options
integer, protected :: soil_BGC_option
integer, public, parameter :: &
    SOIL_BGC_SIMPLE        = 1, & ! SIMPLE decomposition
    SOIL_BGC_CORPSE        = 2, & ! CORPSE model
    SOIL_BGC_GIMICS        = 3    ! GIMICS model

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_soil_BGC_namelist()
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, __FILE__)
  read (input_nml_file, nml=soil_BGC_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_BGC_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=soil_BGC_nml)
  endif

  ! parse soil carbon option
  select case (model_to_use)
  case('SIMPLE')
    soil_BGC_option = SOIL_BGC_SIMPLE
  case('CORPSE')
    soil_BGC_option = SOIL_BGC_CORPSE
  case('GIMICS')
    soil_BGC_option = SOIL_BGC_GIMICS
  case default
    call error_mesg('read_soil_BGC_namelist', &
        '"'//trim(model_to_use)//'" is an invalid option for model_to_use', FATAL)
  end select

  select case (soil_BGC_option)
  case (SOIL_BGC_SIMPLE)
    call read_soil_BGC_SIMPLE_namelist()
  case (SOIL_BGC_CORPSE)
    call read_soil_BGC_CORPSE_namelist()
  case (SOIL_BGC_GIMICS)
    call read_soil_BGC_GIMICS_namelist()
  end select
end subroutine read_soil_BGC_namelist

!> @brief Create new empty soil carbon container
!! @return Pointer to new allocated and initialized soil carbon container
function soilc_ctor(soil) result(ptr)
  class(soil_BGC_t), pointer :: ptr
  type(soil_tile_type), intent(in) :: soil

  select case (soil_BGC_option)
  case (SOIL_BGC_SIMPLE)
    ptr => new_soilc_SIMPLE(soil)
  case (SOIL_BGC_CORPSE)
    ptr => new_soilc_CORPSE(soil)
  case (SOIL_BGC_GIMICS)
    ptr => new_soilc_GIMICS(soil)
  case default
    call land_error_message('soilc_ctor: The value of soil_BGC_option is invalid. This should never happen. See developer', FATAL)
  end select
end function soilc_ctor

!> @brief Create a copy of soil carbon container
!! @return Pointer to a copy of given soil carbon container
function soilc_copy(soilc) result(ptr)
  class(soil_BGC_t), pointer :: ptr
  class(soil_BGC_t), intent(in) :: soilc

  allocate(ptr, source=soilc)
  ! copy all non-pointer members
  select type(soilc)
  type is (soil_BGC_SIMPLE_t)
     ptr => new_soilc_SIMPLE(soilc)
  type is (soil_BGC_CORPSE_t)
     ptr => new_soilc_CORPSE(soilc)
  type is (soil_BGC_GIMICS_t)
     ptr => new_soilc_GIMICS(soilc)
  class default
    call land_error_message('soilc_copy: The type of soilc is invalid. This should never happen. See developer', FATAL)
  end select
end function soilc_copy

!> @brief Deallocate soil carbon contaner
subroutine delete_soilc(ptr)
  class(soil_BGC_t), pointer :: ptr

  ! no need to deallocate components of soil_tile, because F2003 takes care of
  ! allocatable components deallocation when soil_tile is deallocated
  deallocate(ptr)
end subroutine delete_soilc

end module