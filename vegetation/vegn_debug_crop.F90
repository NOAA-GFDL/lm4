module debug_crop_mod

#include "../shared/debug.inc"

use land_debug_mod, only : is_watch_cell, set_current_point, log_date
use vegn_data_mod, only : FORM_GRASS, LU_CROP, spdata
use vegn_tile_mod, only: vegn_tile_type, vegn_tile_LAI
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
  real :: LAI

  if (.not.is_watch_cell())    return
  if (vegn%landuse.ne.LU_CROP) return

  call log_date('#### debug_crop: '//trim(tag)//' ',lnd%time)
! do k = 1, vegn%n_cohorts
!    associate(cc=>vegn%cohorts(k))
!    write(*,'(i2.2," : layer ",i2.2)',advance='NO') k, cc%layer
!    call dpri('frac',cc%layerfrac)
!    call dpri('height',cc%height)
!    call dpri('zbot',cc%zbot)
!    call dpri('LAI',cc%lai)
!    call dpri('crownarea',cc%crownarea)
!    call dpri('nindivs',cc%nindivs)
!    call dpri('species',spdata(cc%species)%name)
!    end associate ! cc
!    write(*,*)
! enddo
  call dpri('Maize_pday',   vegn%Crop%crop_cal_Maize(8))
  call dpri('Maize_hday',   vegn%Crop%crop_cal_Maize(11))
  call dpri('Soy_pday',     vegn%Crop%crop_cal_Soy(8))
  call dpri('Soy_hday',     vegn%Crop%crop_cal_Soy(11))
  call dpri('SW_pday',      vegn%Crop%crop_cal_SW(8))
  call dpri('SW_hday',      vegn%Crop%crop_cal_SW(11))
  call dpri('WW_pday',      vegn%Crop%crop_cal_WW(8))
  call dpri('WW_hday',      vegn%Crop%crop_cal_WW(11))
  call dpri('Rice_pday',    vegn%Crop%crop_cal_Rice_1(8))
  call dpri('Rice_hday',    vegn%Crop%crop_cal_Rice_1(11))
  call dpri('current_crop', vegn%Crop%current_crop)
  call dpri('plant_opt',    vegn%Crop%plant_opt)
  call dpri('harvest_opt',  vegn%Crop%harvest_opt)
  call dpri('status',       vegn%Crop%status)
  call dpri('LAI',          vegn_tile_LAI(vegn))
  write(*,*)

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

end module debug_crop_mod
