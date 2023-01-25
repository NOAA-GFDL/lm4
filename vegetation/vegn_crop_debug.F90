module crop_debug_mod

#include "../shared/debug.inc"

use land_debug_mod, only : is_watch_cell, set_current_point, log_date
use vegn_data_mod, only : FORM_GRASS, LU_CROP, spdata
use vegn_tile_mod, only: vegn_tile_type
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_data_mod, only : lnd, log_version

implicit none
private

public :: debug_crop
public :: debug_crop_1

contains

! ============================================================================
subroutine debug_crop(vegn, tag)
  type(vegn_tile_type), intent(in) :: vegn
  character(*),         intent(in) :: tag

  integer :: k
  logical :: do_debug

  if (.not.is_watch_cell())         return
  if (vegn%landuse.ne.LU_CROP) return

  do_debug = .FALSE.
  do k = 1,vegn%n_cohorts
    associate(cc=>vegn%cohorts(k))
    do_debug = do_debug.or.((spdata(cc%species)%lifeform.ne.FORM_GRASS).and.(cc%nindivs>0))
    end associate ! cc
  enddo

  if (do_debug) then
     call log_date('#### debug_crop: '//trim(tag)//' ',lnd%time)
     do k = 1, vegn%n_cohorts
        associate(cc=>vegn%cohorts(k))
        write(*,'(i2.2," : layer ",i2.2)',advance='NO') k, cc%layer
        call dpri('frac',cc%layerfrac)
        call dpri('height',cc%height)
        call dpri('zbot',cc%zbot)
        call dpri('LAI',cc%lai)
        ! call dpri('bl',cc%bl)
        ! call dpri('leafarea',cc%leafarea)
        call dpri('crownarea',cc%crownarea)
        call dpri('nindivs',cc%nindivs)
        ! call dpri('gapfrac',spdata(sp)%internal_gap_frac)
        ! call dpri('layerarea',layer_area(cc%layer))
        call dpri('species',spdata(cc%species)%name)
        write(*,*)
        end associate ! cc
     enddo
  endif

end subroutine debug_crop

! ============================================================================
subroutine debug_crop_1(tag)
  character(*), intent(in) :: tag

  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  integer :: k,l ! current point indices

  ce = first_elmt(land_tile_map, lnd%ls)
  do while (loop_over_tiles(ce,tile,l,k))
     call set_current_point(l,k)
     if (.not.associated(tile%vegn)) cycle
     call debug_crop(tile%vegn,tag)
  enddo

end subroutine debug_crop_1

end module crop_debug_mod