!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Land Model 4 (LM4).
!*
!* LM4 is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* LM4 is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with LM4.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
module land_utils_mod

use land_tile_mod, only : land_tile_type, land_tile_enum_type, land_tile_list_type, &
     first_elmt, loop_over_tiles, fptr_r0, fptr_r0i

implicit none
private
! ==== public interfaces =====================================================
public :: put_to_tiles_r0d_fptr
public :: put_to_tiles_r1d_fptr
! ==== end of public interfaces ==============================================

contains

! ============================================================================
subroutine put_to_tiles_r0d_fptr(x2d, tile_map, fptr)
  real, intent(in)                         :: x2d     (:)
  type(land_tile_list_type), intent(inout) :: tile_map(:)
  procedure(fptr_r0)                       :: fptr ! subroutine returning the pointer to the data

  integer :: l
  type(land_tile_enum_type)     :: ce      ! tile list elements enumerator
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  real                , pointer :: ptr     ! pointer to the data element within a tile

  ce = first_elmt( tile_map )
  do while(loop_over_tiles(ce,tileptr,l))
     call fptr(tileptr,ptr)
     if (associated(ptr)) ptr=x2d(l)
  enddo
end subroutine


! ============================================================================
subroutine put_to_tiles_r1d_fptr(x2d, tile_map, fptr)
  real, intent(in)                         :: x2d     (:,:)
  type(land_tile_list_type), intent(inout) :: tile_map(:)
  procedure(fptr_r0i)                      :: fptr ! subroutine returning the pointer to the data

  integer :: l, k
  type(land_tile_enum_type)     :: ce      ! tile list elements enumerator
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  real                , pointer :: ptr     ! pointer to the data element within a tile

  ce = first_elmt( tile_map )
  do while(loop_over_tiles(ce,tileptr,l))
     do k = 1, size(x2d,2)
        call fptr(tileptr,k,ptr)
        if (associated(ptr)) ptr=x2d(l,k)
     enddo
  enddo
end subroutine

end module land_utils_mod
