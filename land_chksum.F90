module land_chksum_mod

#include <fms_platform.h>

use iso_fortran_env
use mpp_mod, only : mpp_chksum, mpp_npes, mpp_get_current_pelist
use netcdf, only: NF90_FILL_DOUBLE, NF90_FILL_INT

implicit none
private

public :: get_land_chksum

interface get_land_chksum
  module procedure get_land_chksum_is
  module procedure get_land_chksum_i0d
  module procedure get_land_chksum_i1d
  module procedure get_land_chksum_r0d
  module procedure get_land_chksum_r1d
  module procedure get_land_chksum_r2d
end interface get_land_chksum

contains

subroutine get_land_chksum_is(variable_data, chksum)
  character(len=32), intent(out) :: chksum !< Checksum value converted to a string
  integer, intent(in) :: variable_data !< Input variable data

  integer, allocatable :: all_pe(:) !< Array with the pelist
  integer(LONG_KIND) :: chksum_val !< Checksum value

  ! Create a vector of pes
  allocate(all_pe(mpp_npes()))
  call mpp_get_current_pelist(all_pe)

  chksum_val = variable_data !mpp_chksum(variable_data, all_pe)

  ! Convert to string
  chksum = ""
  write(chksum, "(Z16)") chksum_val

  deallocate(all_pe)
end subroutine get_land_chksum_is

subroutine get_land_chksum_i0d(variable_data, chksum)
  character(len=32), intent(out) :: chksum
  integer, pointer, intent(in) :: variable_data(:)

  integer, allocatable :: all_pe(:)
  integer(LONG_KIND) :: chksum_val

  ! Create a vector of pes
  allocate(all_pe(mpp_npes()))
  call mpp_get_current_pelist(all_pe)

  chksum_val = mpp_chksum(variable_data, all_pe, mask_val=NF90_FILL_INT)

  ! Convert to string
  chksum = ""
  write(chksum, "(Z16)") chksum_val
end subroutine get_land_chksum_i0d

subroutine get_land_chksum_i1d(variable_data, chksum)
  character(len=32), intent(out) :: chksum
  integer, pointer, intent(in) :: variable_data(:,:)

  integer, allocatable :: all_pe(:)
  integer(LONG_KIND) :: chksum_val

  ! Create a vector of pes
  allocate(all_pe(mpp_npes()))
  call mpp_get_current_pelist(all_pe)

  chksum_val = mpp_chksum(variable_data, all_pe, mask_val=NF90_FILL_INT)

  ! Convert to string
  chksum = ""
  write(chksum, "(Z16)") chksum_val
end subroutine get_land_chksum_i1d

subroutine get_land_chksum_r0d(variable_data, chksum)
  character(len=32), intent(out) :: chksum
  real, pointer, intent(in) :: variable_data(:)

  integer, allocatable :: all_pe(:)
  integer(LONG_KIND) :: chksum_val

  ! Create a vector of pes
  allocate(all_pe(mpp_npes()))
  call mpp_get_current_pelist(all_pe)

  chksum_val = mpp_chksum(variable_data, all_pe, mask_val=NF90_FILL_DOUBLE)

  ! Convert to string
  chksum = ""
  write(chksum, "(Z16)") chksum_val
end subroutine get_land_chksum_r0d

subroutine get_land_chksum_r1d(variable_data, chksum)
  character(len=32), intent(out) :: chksum
  real, pointer, intent(in) :: variable_data(:,:)

  integer, allocatable :: all_pe(:)
  integer(LONG_KIND) :: chksum_val

  ! Create a vector of pes
  allocate(all_pe(mpp_npes()))
  call mpp_get_current_pelist(all_pe)

  chksum_val = mpp_chksum(variable_data, all_pe, mask_val=NF90_FILL_DOUBLE)

  ! Convert to string
  chksum = ""
  write(chksum, "(Z16)") chksum_val
end subroutine get_land_chksum_r1d

subroutine get_land_chksum_r2d(variable_data, chksum)
  character(len=32), intent(out) :: chksum
  real, pointer, intent(in) :: variable_data(:,:,:)

  integer, allocatable :: all_pe(:)
  integer(LONG_KIND) :: chksum_val

  ! Create a vector of pes
  allocate(all_pe(mpp_npes()))
  call mpp_get_current_pelist(all_pe)

  chksum_val = mpp_chksum(variable_data, all_pe, mask_val=NF90_FILL_DOUBLE)

  ! Convert to string
  chksum = ""
  write(chksum, "(Z16)") chksum_val
end subroutine get_land_chksum_r2d

end module
