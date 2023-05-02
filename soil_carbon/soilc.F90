module soilc_mod

use fms_mod, only: check_nml_error, input_nml_file, &
            stdlog, mpp_pe, mpp_root_pe, error_mesg, FATAL
use land_data_mod, only: log_version

use soil_tile_mod, only: soil_tile_type
use soil_carbon_mod, only: soil_carbon_option, &
    SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, SOILC_CORPSE_N, &
    read_soilc_CORPSE_namelist, &
    soilc_t, soilc_CENT_t, soilc_CORPSE_t, &
    soilc_CENT_copy, soilc_CORPSE_copy, soilc_ALL_ctor

implicit none; private

public :: read_soil_carbon_namelist
public :: new_soilc, delete_soilc

public :: save_soilc_equilibration_data

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
  case('CENTURY-like-by-layer')
    soil_carbon_option = SOILC_CENTURY_BY_LAYER
  case('CORPSE')
    soil_carbon_option = SOILC_CORPSE
  case('CORPSE-N')
    soil_carbon_option = SOILC_CORPSE_N
  case default
    call error_mesg('read_soil_carbon_namelist', &
        '"'//trim(soil_carbon_model_to_use)//'" is an invalid option for soil_carbon_model_to_use', FATAL)
  end select

  select case (soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
!    call read_soilc_CENT_namelist()
  case (SOILC_CORPSE, SOILC_CORPSE_N)
    call read_soilc_CORPSE_namelist()
  end select
end subroutine read_soil_carbon_namelist

! ============================================================================
function soilc_ctor(soil) result(ptr)
  class(soilc_t), pointer :: ptr
  type(soil_tile_type), intent(in) :: soil

  select case (soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
    ptr => soilc_ALL_ctor(soil)
  case (SOILC_CORPSE, SOILC_CORPSE_N)
    ptr => soilc_ALL_ctor(soil)
  end select
end function soilc_ctor

! ============================================================================
function soilc_copy(soilc) result(ptr)
  class(soilc_t), pointer :: ptr
  class(soilc_t), intent(in) :: soilc

  allocate(ptr, source=soilc)
  ! copy all non-pointer members
  select type(soilc)
  type is (soilc_CORPSE_t)
      ptr => soilc_CORPSE_copy(soilc)
  type is (soilc_CENT_t)
      ptr => soilc_CENT_copy(soilc)
  end select
end function soilc_copy

! ============================================================================
subroutine delete_soilc(ptr)
  type(soilc_t), pointer :: ptr

  ! no need to deallocate components of soil_tile, because F2003 takes care of
  ! allocatable components deallocation when soil_tile is deallocated
  deallocate(ptr)
end subroutine delete_soilc

end module