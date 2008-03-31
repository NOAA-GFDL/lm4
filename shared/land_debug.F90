module land_debug_mod

use fms_mod, only: &
     error_mesg, file_exist, open_namelist_file, check_nml_error, stdlog, &
     write_version_number, close_file, mpp_pe, mpp_npes, mpp_root_pe, FATAL, NOTE
use grid_mod, only: &
     get_grid_ntiles
implicit none
private

! ==== public interfaces =====================================================
public :: land_debug_init
public :: land_debug_end

public :: set_current_point
public :: is_watch_point

! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'land_debug',&
    version     = '$Id: land_debug.F90,v 15.0.2.3 2007/12/05 19:41:35 slm Exp $',&
    tagname     = '$Name: omsk_2008_03 $'

! ==== module variables ======================================================
integer :: current_debug_level = 0
integer :: mosaic_tile = 0

!---- namelist ---------------------------------------------------------------
integer :: watch_point(4)=(/0,0,0,1/) ! coordinates of the point of interest, i,j,tile,mosaic_tile

namelist/land_debug_nml/ watch_point


contains

! ============================================================================
subroutine land_debug_init()
  ! ---- local vars
  integer :: unit, ierr, io, ntiles

  call write_version_number(version, tagname)
  
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=land_debug_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'soil_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  if (mpp_pe() == mpp_root_pe()) then
     write (stdlog(), nml=land_debug_nml)
  endif
  ! set number of our mosaic tile 
  call get_grid_ntiles('LND',ntiles)
  mosaic_tile = ntiles*mpp_pe()/mpp_npes() + 1  ! assumption

end subroutine land_debug_init

! ============================================================================
subroutine land_debug_end()
end subroutine land_debug_end

! ============================================================================
subroutine set_current_point(i,j,k)
  integer, intent(in) :: i,j,k

  current_debug_level = 0
  if ( watch_point(1)==i.and. &
       watch_point(2)==j.and. &
       watch_point(3)==k.and. &
       watch_point(4)==mosaic_tile) then
     current_debug_level = 1
  endif
end subroutine set_current_point

! ============================================================================
function is_watch_point()
  logical :: is_watch_point
  is_watch_point = (current_debug_level > 0)
end function is_watch_point

end module land_debug_mod
