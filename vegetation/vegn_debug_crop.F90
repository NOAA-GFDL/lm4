module vegn_debug_crop_mod

#include "../shared/debug.inc"

use land_debug_mod, only : is_watch_cell, set_current_point, log_date
use vegn_data_mod, only : FORM_GRASS, LU_CROP, spdata
use vegn_tile_mod, only : vegn_tile_type, vegn_tile_LAI
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, first_elmt, loop_over_tiles
use land_data_mod, only : lnd, log_version
use time_manager_mod, only: get_date

implicit none
private

public :: debug_crop
public :: debug_crop_1
public :: date_string ! watchpoint_code
public :: is_crop_watchpoint ! watchpoint_code

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
  call dpri('potential_crop(1)', vegn%Crop%potential_crop(1))
  call dpri('potential_crop(2)', vegn%Crop%potential_crop(2))
  call dpri('potential_crop(3)', vegn%Crop%potential_crop(3))
  call dpri('potential_crop(4)', vegn%Crop%potential_crop(4))
  call dpri('potential_crop(5)', vegn%Crop%potential_crop(5))
  call dpri('potential_crop(6)', vegn%Crop%potential_crop(6))
  call dpri('potential_crop(7)', vegn%Crop%potential_crop(7))
  call dpri('potential_crop(8)', vegn%Crop%potential_crop(8))
  call dpri('potential_crop(9)', vegn%Crop%potential_crop(9))
  call dpri('potential_crop(10)',vegn%Crop%potential_crop(10))
  call dpri('chosen_calendars(1,1)',vegn%Crop%chosen_calendars(1,1))
  call dpri('chosen_calendars(2,1)',vegn%Crop%chosen_calendars(2,1))
  call dpri('chosen_crop(1)',vegn%Crop%chosen_crop(1))
  call dpri('chosen_crop(2)',vegn%Crop%chosen_crop(2))
  call dpri('crop_calendars(1,1,1,1)', vegn%Crop%crop_calendars(1,1,1,1))
  call dpri('crop_calendars(2,1,1,1)', vegn%Crop%crop_calendars(2,1,1,1))
  call dpri('crop_calendars(1,1,2,1)', vegn%Crop%crop_calendars(1,1,2,1))
  call dpri('crop_calendars(2,1,2,1)', vegn%Crop%crop_calendars(2,1,2,1))
  call dpri('crop_calendars(1,1,1,2)', vegn%Crop%crop_calendars(1,1,1,2))
  call dpri('crop_calendars(2,1,1,2)', vegn%Crop%crop_calendars(2,1,1,2))
  call dpri('crop_calendars(1,1,2,2)', vegn%Crop%crop_calendars(1,1,2,2))
  call dpri('crop_calendars(2,1,2,2)', vegn%Crop%crop_calendars(2,1,2,2))
  call dpri('crop_calendars(1,1,1,3)', vegn%Crop%crop_calendars(1,1,1,3))
  call dpri('crop_calendars(2,1,1,3)', vegn%Crop%crop_calendars(2,1,1,3))
  call dpri('crop_calendars(1,1,2,3)', vegn%Crop%crop_calendars(1,1,2,3))
  call dpri('crop_calendars(2,1,2,3)', vegn%Crop%crop_calendars(2,1,2,3))
  call dpri('crop_calendars(1,1,1,4)', vegn%Crop%crop_calendars(1,1,1,4))
  call dpri('crop_calendars(2,1,1,4)', vegn%Crop%crop_calendars(2,1,1,4))
  call dpri('crop_calendars(1,1,2,4)', vegn%Crop%crop_calendars(1,1,2,4))
  call dpri('crop_calendars(2,1,2,4)', vegn%Crop%crop_calendars(2,1,2,4))
  call dpri('crop_calendars(1,1,1,5)', vegn%Crop%crop_calendars(1,1,1,5))
  call dpri('crop_calendars(2,1,1,5)', vegn%Crop%crop_calendars(2,1,1,5))
  call dpri('crop_calendars(1,1,2,5)', vegn%Crop%crop_calendars(1,1,2,5))
  call dpri('crop_calendars(2,1,2,5)', vegn%Crop%crop_calendars(2,1,2,5))
  call dpri('crop_calendars(1,1,1,6)', vegn%Crop%crop_calendars(1,1,1,6))
  call dpri('crop_calendars(2,1,1,6)', vegn%Crop%crop_calendars(2,1,1,6))
  call dpri('crop_calendars(1,1,2,6)', vegn%Crop%crop_calendars(1,1,2,6))
  call dpri('crop_calendars(2,1,2,6)', vegn%Crop%crop_calendars(2,1,2,6))
  call dpri('crop_calendars(1,1,1,7)', vegn%Crop%crop_calendars(1,1,1,7))
  call dpri('crop_calendars(2,1,1,7)', vegn%Crop%crop_calendars(2,1,1,7))
  call dpri('crop_calendars(1,1,2,7)', vegn%Crop%crop_calendars(1,1,2,7))
  call dpri('crop_calendars(2,1,2,7)', vegn%Crop%crop_calendars(2,1,2,7))
  call dpri('crop_calendars(1,1,1,8)', vegn%Crop%crop_calendars(1,1,1,8))
  call dpri('crop_calendars(2,1,1,8)', vegn%Crop%crop_calendars(2,1,1,8))
  call dpri('crop_calendars(1,1,2,8)', vegn%Crop%crop_calendars(1,1,2,8))
  call dpri('crop_calendars(2,1,2,8)', vegn%Crop%crop_calendars(2,1,2,8))
  call dpri('crop_calendars(1,1,1,9)', vegn%Crop%crop_calendars(1,1,1,9))
  call dpri('crop_calendars(2,1,1,9)', vegn%Crop%crop_calendars(2,1,1,9))
  call dpri('crop_calendars(1,1,2,9)', vegn%Crop%crop_calendars(1,1,2,9))
  call dpri('crop_calendars(2,1,2,9)', vegn%Crop%crop_calendars(2,1,2,9))
  call dpri('crop_calendars(1,1,1,10)',vegn%Crop%crop_calendars(1,1,1,10))
  call dpri('crop_calendars(2,1,1,10)',vegn%Crop%crop_calendars(2,1,1,10))
  call dpri('crop_calendars(1,1,2,10)',vegn%Crop%crop_calendars(1,1,2,10))
  call dpri('crop_calendars(2,1,2,10)',vegn%Crop%crop_calendars(2,1,2,10))
  call dpri('LAI',vegn_tile_LAI(vegn))
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
! ============================================================================
  function date_string() result(string) ! watchpoint_code
  character(len=19) :: string
  integer :: year,month,day,hour,minute,second

  string = '    -  -     :  :  '
  call get_date(lnd%time, year,month,day,hour,minute,second)
  write(string( 1: 4),'(i4.4)') year
  write(string( 6: 7),'(i2.2)') month
  write(string( 9:10),'(i2.2)') day
  write(string(12:13),'(i2.2)') hour
  write(string(15:16),'(i2.2)') minute
  write(string(18:19),'(i2.2)') second
  end function date_string
! ============================================================================
  logical function is_crop_watchpoint(L,vegn,crop_watchpoint) ! watchpoint_code
  integer, intent(in) :: L
  type(vegn_tile_type), intent(in) :: vegn
  integer, intent(in) :: crop_watchpoint(3)

  if(vegn%landuse /= LU_CROP) then
    is_crop_watchpoint = .FALSE.
    return
  endif
  if(crop_watchpoint(1) == lnd%i_index(L) .and. crop_watchpoint(2) == lnd%j_index(L) .and. crop_watchpoint(3) == lnd%sg_face) then
    is_crop_watchpoint = .TRUE.
  else
    is_crop_watchpoint = .FALSE.
  endif
  end function is_crop_watchpoint
! ============================================================================

end module vegn_debug_crop_mod
