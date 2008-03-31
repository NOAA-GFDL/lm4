module cana_tile_mod

use land_tile_selectors_mod, only : &
     tile_selector_type

implicit none
private

! ==== public interfaces =====================================================
public :: cana_prog_type
public :: cana_tile_type

public :: new_cana_tile, delete_cana_tile
public :: cana_tiles_can_be_merged, merge_cana_tiles
public :: get_cana_tile_tag
public :: cana_is_selected

! ==== end of public interfaces ==============================================
interface new_cana_tile
   module procedure cana_tile_ctor
   module procedure cana_tile_copy_ctor
end interface

type :: cana_prog_type
  real T
  real q
  real :: co2 ! co2 concentration in canopy air, mol/mol
end type cana_prog_type

type :: cana_tile_type
   type(cana_prog_type) :: prog
end type cana_tile_type


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! =============================================================================
function cana_tile_ctor() result(ptr)
  type(cana_tile_type), pointer :: ptr ! return value

  allocate(ptr)
end function cana_tile_ctor

! =============================================================================
function cana_tile_copy_ctor(cana) result(ptr)
  type(cana_tile_type), pointer :: ptr ! return value
  type(cana_tile_type), intent(in) :: cana ! return value

  allocate(ptr)
  ptr = cana
end function cana_tile_copy_ctor

! =============================================================================
subroutine delete_cana_tile(cana)
  type(cana_tile_type), pointer :: cana

  deallocate(cana)
end subroutine delete_cana_tile

! =============================================================================
function cana_tiles_can_be_merged(cana1,cana2) result(response)
  logical :: response
  type(cana_tile_type), intent(in) :: cana1,cana2

  response = .TRUE.
end function

! =============================================================================
subroutine merge_cana_tiles(cana1,w1,cana2,w2)
  type(cana_tile_type), intent(in)    :: cana1
  type(cana_tile_type), intent(inout) :: cana2
  real                , intent(in)    :: w1, w2
  
  ! ---- local vars
  real :: x1,x2 ! normalized weights
  
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1-x1
  
  cana2%prog%T = cana1%prog%T*x1+cana2%prog%T*x2
  cana2%prog%q = cana1%prog%q*x1+cana2%prog%q*x2
end subroutine

! =============================================================================
! returns tag of the tile
function get_cana_tile_tag(cana) result(tag)
  integer :: tag
  type(cana_tile_type), intent(in) :: cana
  
  tag = 1
end function

! =============================================================================
! returns true if tile fits the specified selector
function cana_is_selected(cana, sel)
  logical cana_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(cana_tile_type),      intent(in) :: cana

  cana_is_selected = .TRUE.
end function


end module cana_tile_mod
