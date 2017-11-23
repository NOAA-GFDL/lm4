module land_utils_mod

use land_debug_mod, only : check_conservation, water_cons_tol, carbon_cons_tol, &
     nitrogen_cons_tol, heat_cons_tol, do_check_conservation
use soil_carbon_mod, only : soil_carbon_option, SOILC_CORPSE_N
use land_tile_mod, only : land_tile_type, land_tile_enum_type, land_tile_list_type, &
     first_elmt, loop_over_tiles, fptr_r0, fptr_r0i, &
     get_tile_water, land_tile_carbon, land_tile_nitrogen, land_tile_heat

implicit none
private
! ==== public interfaces =====================================================
public :: put_to_tiles_r0d_fptr
public :: put_to_tiles_r1d_fptr
public :: check_conservation_1, check_conservation_2
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

! ============================================================================
! + conservation check, part 1: calculate totals. It cannot be in land_debug_mod 
! because it uses land tile type.
subroutine check_conservation_1(tile,lmass,fmass,cmass,nmass,heat)
  type(land_tile_type), intent(in) :: tile
  real, intent(out), optional :: lmass,fmass,cmass,nmass,heat ! stocks to check against

  real :: lmass1,fmass1

  if (.not.do_check_conservation) return

  call get_tile_water(tile,lmass1,fmass1)
  if (present(lmass)) lmass = lmass1
  if (present(fmass)) fmass = fmass1
  if (present(cmass)) cmass = land_tile_carbon(tile)
  if (present(nmass).and.soil_carbon_option==SOILC_CORPSE_N) then
     nmass  = land_tile_nitrogen(tile)
  endif
  if (present(heat)) heat  = land_tile_heat(tile)
end subroutine check_conservation_1

! ============================================================================
! + conservation check, part 2: calculate totals in final state, and compare
! with previous totals. It cannot be in land_debug_mod because it uses land tile type.
subroutine check_conservation_2(tile,tag,lmass,fmass,cmass,nmass,heat)
  type(land_tile_type), intent(in) :: tile
  character(*), intent(in) :: tag
  real, intent(in), optional :: lmass,fmass,cmass,nmass,heat ! stocks to check against

  real :: lmass1,fmass1,cmass1,nmass1,heat1
  if (.not.do_check_conservation) return

  if (present(lmass).or.present(fmass)) then
     call get_tile_water(tile,lmass1,fmass1)
     if(present(lmass)) call check_conservation (tag,'liquid water', lmass, lmass1, water_cons_tol)
     if(present(lmass)) call check_conservation (tag,'frozen water', fmass, fmass1, water_cons_tol)
  endif
  if (present(cmass)) then
     cmass1 = land_tile_carbon(tile)
     call check_conservation (tag,'carbon', cmass, cmass1, carbon_cons_tol)
  endif
  if (present(nmass).and.soil_carbon_option==SOILC_CORPSE_N) then
     nmass1  = land_tile_nitrogen(tile)
     call check_conservation (tag,'nitrogen', nmass, nmass1, nitrogen_cons_tol)
  endif
  if (present(heat)) then
     heat1  = land_tile_heat(tile)
     call check_conservation (tag,'heat content', heat, heat1, heat_cons_tol)
  endif
end subroutine check_conservation_2

end module land_utils_mod
