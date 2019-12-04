module land_chksum_mod

#include <fms_platform.h>

use iso_fortran_env
use mpp_mod, only : mpp_chksum, mpp_npes
use netcdf, only: NF90_FILL_DOUBLE, NF90_FILL_INT

contains

subroutine get_land_chksum_is(variable_data, chksum)
	character(len=32), intent(out) :: chksum
	integer, intent(in) :: variable_data

   integer, allocatable :: all_pe(:)
   integer(LONG_KIND) :: chksum_val
	integer :: i
	integer :: npes

	! Get the total number of pes
	npes = mpp_npes()

	! Create a vector of pes
  	allocate(all_pe(npes))
  	DO i = 1,npes
    	all_pe(i) = i-1
  	enddo

   chksum_val = variable_data !mpp_chksum(variable_data, all_pe)

	! Convert to string 
  	chksum = ""
   write(chksum, "(Z16)") chksum_val

end subroutine get_land_chksum_is

subroutine get_land_chksum_i0d(variable_data, chksum)
	character(len=32), intent(out) :: chksum
	integer, pointer, intent(in) :: variable_data(:)

   integer, allocatable :: all_pe(:)
   integer(LONG_KIND) :: chksum_val
	integer :: i
	integer :: npes

	! Get the total number of pes
	npes = mpp_npes()

	! Create a vector of pes
  	allocate(all_pe(npes))
  	DO i = 1,npes
    	all_pe(i) = i-1
  	enddo

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
	integer :: i
	integer :: npes

	! Get the total number of pes
	npes = mpp_npes()

	! Create a vector of pes
  	allocate(all_pe(npes))
  	DO i = 1,npes
    	all_pe(i) = i-1
  	enddo

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
	integer :: i
	integer :: npes

	! Get the total number of pes
	npes = mpp_npes()

	! Create a vector of pes
  	allocate(all_pe(npes))
  	DO i = 1,npes
    	all_pe(i) = i-1
  	enddo

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
	integer :: i
	integer :: npes

	! Get the total number of pes
	npes = mpp_npes()

	! Create a vector of pes
  	allocate(all_pe(npes))
  	DO i = 1,npes
    	all_pe(i) = i-1
  	enddo

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
	integer :: i
	integer :: npes

	! Get the total number of pes
	npes = mpp_npes()

	! Create a vector of pes
  	allocate(all_pe(npes))
  	DO i = 1,npes
    	all_pe(i) = i-1
  	enddo

   chksum_val = mpp_chksum(variable_data, all_pe, mask_val=NF90_FILL_DOUBLE)

	! Convert to string 
  	chksum = ""
   write(chksum, "(Z16)") chksum_val

end subroutine get_land_chksum_r2d

end module
