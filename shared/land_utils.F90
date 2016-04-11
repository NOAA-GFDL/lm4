module land_utils_mod

use land_tile_mod, only : land_tile_type, land_tile_enum_type, land_tile_list_type, &
     first_elmt, tail_elmt, next_elmt, operator(/=), get_elmt_indices, current_tile, &
     fptr_r0, fptr_r0i

implicit none
private
! ==== public interfaces =====================================================
public :: put_to_tiles_r0d_fptr
public :: put_to_tiles_r1d_fptr
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
#include "../shared/version_variable.inc"

contains

! ============================================================================
subroutine put_to_tiles_r0d_fptr(x2d, tile_map, fptr)
  real, intent(in)                         :: x2d     (:,:)
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  procedure(fptr_r0)                       :: fptr ! subroutine returning the pointer to the data

  integer :: i,j
  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  real                , pointer :: ptr     ! pointer to the data element within a tile

  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     call get_elmt_indices(ce,i,j)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if (associated(ptr)) ptr=x2d(i,j)
     ce=next_elmt(ce)
  enddo
end subroutine


! ============================================================================
subroutine put_to_tiles_r1d_fptr(x2d, tile_map, fptr)
  real, intent(in)                         :: x2d     (:,:,:)
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  procedure(fptr_r0i)                      :: fptr ! subroutine returning the pointer to the data

  integer :: i,j,k
  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  real                , pointer :: ptr     ! pointer to the data element within a tile

  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     call get_elmt_indices(ce,i,j)
     tileptr => current_tile(ce)
     do k = 1, size(x2d,3)
        call fptr(tileptr,k,ptr)
        if (associated(ptr)) ptr=x2d(i,j,k)
     enddo
     ce=next_elmt(ce)
  enddo
end subroutine

end module land_utils_mod
