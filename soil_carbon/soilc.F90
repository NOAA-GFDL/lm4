module soilc_mod

use fms_mod, only: check_nml_error, input_nml_file, &
            stdlog, mpp_pe, mpp_root_pe, error_mesg, FATAL
use land_data_mod, only: log_version
use land_debug_mod, only: land_error_message

use soil_tile_mod, only: soil_tile_type
use soilc_type_mod, only: soilc_t
use soilc_CENT_type_mod, only: soilc_CENT_t, new_soilc_CENT, read_soilc_CENT_namelist
use soil_carbon_mod, only:soilc_CORPSE_t, new_soilc_CORPSE, read_soilc_CORPSE_namelist

implicit none; private

public :: read_soil_carbon_namelist
public :: new_soilc, delete_soilc

public :: save_soilc_equilibration_data

public :: soil_carbon_option

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'soilc_mod'
#include "../shared/version_variable.inc"

! ==== module interfaces ======================================================
interface new_soilc
   module procedure soilc_ctor
   module procedure soilc_copy
end interface

!---- namelist ---------------------------------------------------------------
character(32) :: soil_carbon_model_to_use = 'CENTURY-like' ! or 'CENTURY-like-by-layer', or 'CORPSE', or 'CORPSE-N'
logical, protected :: save_soilc_equilibration_data = .FALSE. ! indicates whether to write
                        ! information for soil carbon acceleration
namelist /soil_carbon_nml/ soil_carbon_model_to_use, save_soilc_equilibration_data

! soil carbon options
integer, protected :: soil_carbon_option
integer, public, parameter :: &
    SOILC_CENTURY          = 1, & ! CENTURY-like decomposition
    SOILC_CORPSE           = 3    ! CORPSE model

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_soil_carbon_namelist()
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
  read (input_nml_file, nml=soil_carbon_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_carbon_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=soil_carbon_nml)
  endif

  ! parse soil carbon option
  select case (soil_carbon_model_to_use)
  case('CENTURY-like')
    soil_carbon_option = SOILC_CENTURY
  case('CORPSE')
    soil_carbon_option = SOILC_CORPSE
  case default
    call error_mesg('read_soil_carbon_namelist', &
        '"'//trim(soil_carbon_model_to_use)//'" is an invalid option for soil_carbon_model_to_use', FATAL)
  end select

  select case (soil_carbon_option)
  case (SOILC_CENTURY)
    call read_soilc_CENT_namelist()
  case (SOILC_CORPSE)
    call read_soilc_CORPSE_namelist()
  end select
end subroutine read_soil_carbon_namelist

!> @brief Create new empty soil carbon container
!! @return Pointer to new allocated and initialized soil carbon container
function soilc_ctor(soil) result(ptr)
  class(soilc_t), pointer :: ptr
  type(soil_tile_type), intent(in) :: soil

  select case (soil_carbon_option)
  case (SOILC_CENTURY)
    ptr => new_soilc_CENT(soil)
  case (SOILC_CORPSE)
    ptr => new_soilc_CORPSE(soil)
  case default
    call land_error_message('soilc_ctor: The value of soil_carbon_option is invalid. This should never happen. See developer', FATAL)
  end select
end function soilc_ctor

!> @brief Create a copy of soil carbon container
!! @return Pointer to a copy of given soil carbon container
function soilc_copy(soilc) result(ptr)
  class(soilc_t), pointer :: ptr
  class(soilc_t), intent(in) :: soilc

  allocate(ptr, source=soilc)
  ! copy all non-pointer members
  select type(soilc)
  type is (soilc_CENT_t)
     ptr => new_soilc_CENT(soilc)
  type is (soilc_CORPSE_t)
     ptr => new_soilc_CORPSE(soilc)
  class default
    call land_error_message('soilc_copy: The type of soilc is invalid. This should never happen. See developer', FATAL)
  end select
end function soilc_copy

!> @brief Deallocate soil carbon contaner
subroutine delete_soilc(ptr)
  class(soilc_t), pointer :: ptr

  ! no need to deallocate components of soil_tile, because F2003 takes care of
  ! allocatable components deallocation when soil_tile is deallocated
  deallocate(ptr)
end subroutine delete_soilc

end module