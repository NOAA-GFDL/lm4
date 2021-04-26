module land_debug_mod

use ieee_arithmetic

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif
use mpp_mod, only: mpp_max
use constants_mod, only: PI
use fms_mod, only: error_mesg, file_exist, check_nml_error, stdlog, lowercase, &
     close_file, mpp_pe, mpp_root_pe, string, FATAL, WARNING, NOTE
use time_manager_mod, only : &
     time_type, get_date, set_date, operator(<=), operator(>=)
use land_data_mod, only: lnd, log_version

! NOTE TO SELF: the "!$" sentinels are not comments: they are compiled if OpenMP
! support is turned on
!$ use omp_lib, only: OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM

implicit none
private

! ==== public interfaces =====================================================
public :: land_debug_init
public :: land_debug_end

public :: set_current_point, set_current_point_sg
public :: get_current_point
public :: current_i, current_j, current_k, current_face
public :: is_watch_point
public :: is_watch_cell
public :: is_watch_time
public :: get_watch_point

public :: check_temp_range
public :: check_var_range
public :: check_conservation

public :: land_error_message
public :: log_date
public :: string_from_time
public :: dpri

interface dpri
   module procedure debug_printout_r0d
   module procedure debug_printout_i0d
   module procedure debug_printout_l0d
   module procedure debug_printout_t0d
   module procedure debug_printout_r1d
   module procedure debug_printout_i1d
   module procedure debug_printout_t1d
   module procedure debug_printout_r2d
end interface dpri

interface check_var_range
   module procedure check_var_range_0d
   module procedure check_var_range_1d
end interface check_var_range

interface check_temp_range
   module procedure check_temp_range_0d
   module procedure check_temp_range_1d
end interface check_temp_range

! checksum output trigger: if TRUE, checksums of land fields are printed each
! time step. Useful for reproducibility debugging.
public :: do_checksums
! conservation tolerances for use across the code. This module does not use
! them, just serves as a convenient place to share them across all land code
public :: water_cons_tol
public :: carbon_cons_tol, nitrogen_cons_tol
public :: heat_cons_tol
public :: do_check_conservation

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_debug_mod'
#include "../shared/version_variable.inc"

! ==== module variables ======================================================
integer              :: mosaic_tile_sg = 0, mosaic_tile_ug = 0
integer, allocatable :: curr_i(:), curr_j(:), curr_k(:), curr_l(:)
logical, allocatable :: watched_tile(:), watched_cell(:)
type(time_type)      :: start_watch_time, stop_watch_time
character(128)       :: fixed_label_format
integer              :: watch_point_lindex = 0  ! watch point index in unstructured grid.

!---- namelist ---------------------------------------------------------------
integer :: watch_point(4)=[0,0,0,1] ! coordinates of the point of interest,
           ! i, j, subgrid tile, cubic sphere face
integer :: start_watching(6) = [    1, 1, 1, 0, 0, 0 ] ! YYYY, MM, DD, HH, MM, SS
integer :: stop_watching(6)  = [ 9999, 1, 1, 0, 0, 0 ] ! YYYY, MM, DD, HH, MM, SS
logical :: watch_conservation = .FALSE. ! if true, conservation check reports are
           ! printed for watch_point, in addition to regular debug output
logical :: print_hex_debug = .FALSE. ! if TRUE, hex representation of debug
           ! values is also printed
integer :: label_len = 12  ! minimum length of text labels for debug output
logical :: trim_labels = .FALSE. ! if TRUE, the length of text labels in debug
           ! printout is never allowed to exceed label_len, resulting in
           ! trimming of the labels. Set it to TRUE to match earlier debug
           ! printout
character(64) :: value_format = 'short'
logical, protected :: do_checksums = .FALSE. ! if TRUE, checksums of land fields are printed each
           ! time step. Useful for reproducibility debugging.
real    :: temp_lo = 120.0 ! lower limit of "reasonable" temperature range, deg K
real    :: temp_hi = 373.0 ! upper limit of "reasonable" temperature range, deg K

namelist/land_debug_nml/ watch_point, &
   start_watching, stop_watching, watch_conservation, &
   value_format, print_hex_debug, label_len, trim_labels, &
   do_checksums, &
   temp_lo, temp_hi

logical, protected :: do_check_conservation = .FALSE.
real, protected    :: water_cons_tol  = 1e-11 ! tolerance of water conservation checks
real, protected    :: carbon_cons_tol = 1e-13 ! tolerance of carbon conservation checks
real, protected    :: heat_cons_tol   = 1e-8  ! tolerance of heat conservation checks
real, protected    :: nitrogen_cons_tol = 1e-13 ! tolerance of nitrogen conservation checks, kgN/m2
namelist/land_conservation_nml/ do_check_conservation, &
      water_cons_tol, carbon_cons_tol, heat_cons_tol, nitrogen_cons_tol

contains

! ============================================================================
subroutine land_debug_init()
  ! ---- local vars
  integer :: unit, ierr, io, ntiles, l
  integer :: max_threads

  call log_version(version, module_name, &
  __FILE__)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=land_debug_nml, iostat=io)
  ierr = check_nml_error(io, 'land_debug_nml')
  read (input_nml_file, nml=land_conservation_nml, iostat=io)
  ierr = check_nml_error(io, 'land_conservation_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=land_debug_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_debug_nml')
     enddo
10   continue
     call close_file (unit)
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=land_conservation_nml, iostat=io, end=11)
        ierr = check_nml_error (io, 'land_conservation_nml')
     enddo
11   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=land_debug_nml)
     write(unit, nml=land_conservation_nml)
  endif

  ! set number of threads and allocate by-thread arrays
  max_threads = 1
!$  max_threads = OMP_GET_MAX_THREADS()
  allocate(curr_i(max_threads),curr_j(max_threads),curr_k(max_threads),curr_l(max_threads))
  allocate(watched_tile(max_threads), watched_cell(max_threads))
  watched_tile(:) = .FALSE. ; watched_cell(:) = .FALSE.

  ! construct the label format string for output
  fixed_label_format = '(a'//trim(string(label_len))//')'

  ! construct value format
  if (trim(lowercase(value_format))=='short') then
     value_format = '(g13.6)'
  else if (trim(lowercase(value_format))=='long') then
     value_format = '(g23.16)'
  else if (trim(lowercase(value_format))=='full') then
     value_format = '(g23.16)'
  else
     ! do nothing: use value_format provided in the namelis for printing
  endif

  start_watch_time = set_date(start_watching(1), start_watching(2), start_watching(3), &
                              start_watching(4), start_watching(5), start_watching(6)  )
  stop_watch_time  = set_date( stop_watching(1),  stop_watching(2),  stop_watching(3), &
                               stop_watching(4),  stop_watching(5),  stop_watching(6)  )

  ! Set up the unstructured grid index of the watch point.
  watch_point_lindex = 0
  do l = lnd%ls, lnd%le
     if ( watch_point(1) == lnd%i_index(l) .AND. &
          watch_point(2) == lnd%j_index(l)       ) then
        watch_point_lindex = l
     endif
  enddo

end subroutine land_debug_init

! ============================================================================
subroutine land_debug_end()
  deallocate(curr_i,curr_j,curr_k,curr_l)
  deallocate(watched_tile, watched_cell)
end subroutine

! ============================================================================
! The primary purpose of set_current_point and set_current_point_sg is to set the the
! triggers for debug output on a specific subgrid tile or grid cell, controlled by
! parameters in land_debug_nml. The coordinates of current point can also be retrieved,
! which is used in conservation-checking and land error message subroutines
!
! Since some parts of the land model work on unstructured grid, and some other parts on
! structured grid, there are two subroutines, with set_current_point used in unstructured
! and set_current_point_sg -- in structured grid parts.
!
! Typical usage pattern on unstructured grid:
!
! do while(loop_over_tiles(ce,tile,l,k))
!     call set_current_point(l,k)
!     ...
!     if (is_watch_point()) then
!         ... print out debug information for this subgrid tile
!     endif
!     if (is_watch_cell()) then
!         ... print out debug information for the cell, regardless of tile index
!     endif
! enddo
!
! Typical usage pattern on structured grid (e.g. in river code):
!
! do j = jsc, jec
! do i = isc, iec
! do k = 1, n_subgrid_tiles
!     call set_current_point_sg(i,j,k)
!     if (is_watch_point()) then
!         ... print out debug information for this subgrid tile
!     endif
!     if (is_watch_cell()) then
!         ... print out debug information for the cell, regardless of tile index
!     endif
! enddo
! enddo
! enddo
!
! Note that the difference in more subtle than just the input parameters: on the same
! processor, the cubic sphere face number may be different for structured and unstructured
! grid.

subroutine set_current_point(l,k) ! for unstructured grid
  integer, intent(in) :: l ! grid cell index on unstructured grid
  integer, intent(in) :: k ! subgrid tile index

  integer :: thread
  ! set the coordinates of the current point
  thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  curr_i(thread) = lnd%i_index(l) ; curr_j(thread) = lnd%j_index(l)
  curr_k(thread) = k; curr_l(thread) = l

  ! trigger watch tile/cell indicators
  watched_tile(thread) = .FALSE. ; watched_cell(thread) = .FALSE.
  if ( watch_point_lindex /= l ) return
  if ( watch_point(4) /= lnd%ug_face ) return
  watched_cell(thread) = .TRUE.
  if ( watch_point(3) /= curr_k(thread) ) return
  watched_tile(thread) = .TRUE.
end subroutine set_current_point


subroutine set_current_point_sg(i,j,k)  ! for unstructured grid
  integer, intent(in) :: i,j  ! grid cell horizontal indices
  integer, intent(in) :: k    ! subgrid tile index

  integer :: thread
  ! set the coordinates of the current point
  thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  curr_i(thread) = i ; curr_j(thread) = j ; curr_k(thread) = k; curr_l(thread) = 0

  ! trigger watch tile/cell indicators
  watched_tile(thread) = .FALSE. ; watched_cell(thread) = .FALSE.
  if ( watch_point(1)/=i ) return
  if ( watch_point(2)/=j ) return
  if ( watch_point(4)/=lnd%sg_face) return
  watched_cell(thread) = .TRUE.
  if ( watch_point(3)/=k ) return
  watched_tile(thread) = .TRUE.
end subroutine set_current_point_sg

! ============================================================================
subroutine get_current_point(i,j,k,face)
  integer, intent(out), optional :: i,j  ! grid cell indices
  integer, intent(out), optional :: k    ! subgrid tile index
  integer, intent(out), optional :: face ! cubic sphere face; could be unstructured on structured
    ! grid face, depending on what the current grid is (more specifically, was current
    ! point set by set_current_point or set_current_point_sg)

  integer :: thread
  thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1

  if (present(i)) i = curr_i(thread)
  if (present(j)) j = curr_j(thread)
  if (present(k)) k = curr_k(thread)
  if (present(face)) then
     ! curr_l used as an indicator of the structured/unstructured grid : set_current_point_sg
     ! sets it to zero, while the set_current_point (which takes unstructured grid index
     ! and is supposed to be used on unstructured grid) sets it to a non-zero value.
     ! Assuming lnd%ls:lnd%le range does not include zero, which is currently always true.
     if (curr_l(thread) > 0) then
        face = lnd%ug_face
     else
        face = lnd%sg_face
     endif
  endif
end subroutine get_current_point

! ============================================================================
integer function current_i()
  integer :: thread
  thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  current_i = curr_i(thread)
end function

integer function current_j()
  integer :: thread
  thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  current_j = curr_j(thread)
end function

integer function current_k()
  integer :: thread
  thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  current_k = curr_k(thread)
end function

integer function current_face()
  integer :: thread
  thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1

  ! curr_l used as an indicator of the structured/unstructured grid : set_current_point_sg
  ! sets it to zero, while the set_current_point (which takes unstructured grid index
  ! and is supposed to be used on unstructured grid) sets it to non-zero value.
  ! This assumes that lnd%ls:lnd%le range never includes zero.
  if (curr_l(thread) > 0) then
     current_face = lnd%ug_face
  else
     current_face = lnd%sg_face
  endif
end function

! ============================================================================
logical function is_watch_time()
   is_watch_time = lnd%time >= start_watch_time &
             .and. lnd%time <= stop_watch_time
end function is_watch_time

! ============================================================================
logical function is_watch_point()
  integer :: thread
  thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1
  is_watch_point=.FALSE.
  if (.not.watched_tile(thread)) return
  if (.not.is_watch_time()) return
  is_watch_point=.TRUE.
end function is_watch_point

! ============================================================================
! returns true, if the watch point is within the grid cell, regardless of
! the tile number
logical function is_watch_cell()
  integer :: thread
  thread = 1
!$  thread = OMP_GET_THREAD_NUM()+1

  is_watch_cell = .FALSE.
  if (.not.watched_cell(thread)) return
  if (.not.is_watch_time()) return
  is_watch_cell = .TRUE.
end function is_watch_cell


! ============================================================================
subroutine get_watch_point(i,j,k,face,l)
  integer, intent(out), optional :: i,j,k,face,l
  if (present(i)) i = watch_point(1)
  if (present(j)) j = watch_point(2)
  if (present(k)) k = watch_point(3)
  if (present(face)) face = watch_point(4)
  if (present(l)) l = watch_point_lindex
end subroutine get_watch_point

! ============================================================================
! checks if the temperature within reasonable range, and prints a message
! if it is not
subroutine check_temp_range_0d(temp, tag, varname)
  real, intent(in) :: temp ! temperature to check
  character(*), intent(in) :: tag ! tag to print
  character(*), intent(in) :: varname ! name of the variable for printout

  call check_var_range(temp,temp_lo,temp_hi,tag,varname,WARNING)
end subroutine check_temp_range_0d

subroutine check_temp_range_1d(temp, tag, varname)
  real, intent(in) :: temp(:) ! temperature to check
  character(*), intent(in) :: tag ! tag to print
  character(*), intent(in) :: varname ! name of the variable for printout

  call check_var_range(temp,temp_lo,temp_hi,tag,varname,WARNING)
end subroutine check_temp_range_1d

! ============================================================================
! checks if the value is within specified range, and prints a message
! if it is not
subroutine check_var_range_0d(value, lo, hi, tag, varname, severity)
  real        , intent(in) :: value    ! value to check
  real        , intent(in) :: lo,hi    ! lower and upper bounds of acceptable range
  character(*), intent(in) :: tag      ! tag to print
  character(*), intent(in) :: varname  ! name of the variable for printout
  integer     , intent(in) :: severity ! severity of the non-conservation error:
         ! Can be WARNING, FATAL, or negative. Negative means check is not done.

  ! ---- local vars
  integer :: y,mo,d,h,m,s ! components of date
  integer :: thread, face
  real    :: lon, lat ! current coordinates, degree
  character(512) :: message
  character :: ew,ns  ! hemisphere indicators

  if (severity<0) return

  if(ieee_is_finite(value)) then
      if(lo<=value.and.value<=hi) return
  endif

  thread = 1
!$   thread = OMP_GET_THREAD_NUM()+1
  call get_date(lnd%time,y,mo,d,h,m,s)
  call get_current_coordinates(thread, lon, lat, face)
  ew = 'E'; if (lon<0) ew = 'W'
  ns = 'N'; if (lat<0) ns = 'S'
  write(message,'(a,g23.16,a,2(f7.2,a),a,3(i4,",")i4,a,i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))')&
       trim(varname)//' out of range: value=', value, ' at', abs(lon),ew,abs(lat),ns, &
       '  (i,j,tile,face) = (',curr_i(thread),curr_j(thread),curr_k(thread),face, &
       ') time=',y,mo,d,h,m,s
  call error_mesg(trim(tag),trim(message),severity)
end subroutine check_var_range_0d


! ============================================================================
subroutine check_var_range_1d(value, lo, hi, tag, varname, severity)
  real        , intent(in) :: value(:) ! value to check
  real        , intent(in) :: lo,hi    ! lower and upper bounds of acceptable range
  character(*), intent(in) :: tag      ! tag to print
  character(*), intent(in) :: varname  ! name of the variable for printout
  integer     , intent(in) :: severity ! severity of the non-conservation error:
         ! Can be WARNING, FATAL, or negative. Negative means check is not done.

  ! ---- local vars
  integer :: i
  integer :: y,mo,d,h,m,s ! components of date
  integer :: thread, face
  real    :: lon, lat ! current coordinates, degree
  character :: ew,ns  ! hemisphere indicators
  character(512) :: message

  if (severity<0) return

  do i = 1,size(value)
     if(ieee_is_finite(value(i))) then
        if(lo<=value(i).and.value(i)<=hi) cycle
     endif
     thread = 1
!$   thread = OMP_GET_THREAD_NUM()+1
     call get_date(lnd%time,y,mo,d,h,m,s)
     call get_current_coordinates(thread, lon, lat, face)
     ew = 'E'; if (lon<0) ew = 'W'
     ns = 'N'; if (lat<0) ns = 'S'
     write(message,'(a,g23.16,a,2(f7.2,a),a,3(i4,",")i4,a,i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))')&
          trim(varname)//'('//trim(string(i))//')'//' out of range: value=', value(i),&
          ' at', abs(lon),ew,abs(lat),ns, &
          '  (i,j,tile,face) = (',curr_i(thread),curr_j(thread),curr_k(thread),face, &
          ') time=',y,mo,d,h,m,s
     call error_mesg(trim(tag),trim(message),severity)
  enddo
end subroutine check_var_range_1d


! ============================================================================
subroutine print_label(description)
  character(*), intent(in) :: description

  if (trim_labels.or.len_trim(description)<label_len) then
     write(*,fixed_label_format,advance='NO')trim(description)
  else
     write(*,'(x,a,g23.16)',advance='NO')trim(description)
  endif
end subroutine print_label

! ============================================================================
! debug printout procedures
subroutine debug_printout_r0d(description,value)
  character(*), intent(in) :: description
  real        , intent(in) :: value

  call print_label(description)
  write(*,value_format,advance='NO')value
  if(print_hex_debug) write(*,'(z17)',advance='NO')value
end subroutine


subroutine debug_printout_i0d(description,value)
  character(*), intent(in) :: description
  integer     , intent(in) :: value

  call print_label(description)
  write(*,value_format,advance='NO')value
end subroutine


subroutine debug_printout_l0d(description,value)
  character(*), intent(in) :: description
  logical     , intent(in) :: value

  call print_label(description)
  write(*,value_format,advance='NO')value
end subroutine

subroutine debug_printout_t0d(description,value)
  character(*), intent(in) :: description
  character(*), intent(in) :: value

  call print_label(description)
  write(*,value_format,advance='NO')'"'//trim(value)//'"'
end subroutine

subroutine debug_printout_r1d(description,values)
  character(*), intent(in) :: description
  real        , intent(in) :: values(:)

  integer :: i

  call print_label(description)
  do i = 1,size(values)
     write(*,value_format,advance='NO')values(i)
     if(print_hex_debug) write(*,'(z17)',advance='NO')values(i)
  enddo
end subroutine

subroutine debug_printout_i1d(description,values)
  character(*), intent(in) :: description
  integer     , intent(in) :: values(:)

  integer :: i

  call print_label(description)
  do i = 1,size(values)
     write(*,value_format,advance='NO')values(i)
  enddo
end subroutine

subroutine debug_printout_t1d(description,values)
  character(*), intent(in) :: description
  character(*), intent(in) :: values(:)

  integer :: i

  call print_label(description)
  do i = 1,size(values)
     write(*,value_format,advance='NO')'"'//trim(values(i))//'"'
  enddo
end subroutine

subroutine debug_printout_r2d(description,values)
  character(*), intent(in) :: description
  real        , intent(in) :: values(:,:)

  integer :: i,j

  ! TODO: print 2D value as a matrix
  call print_label(description)
  do i = 1,size(values,1)
  do j = 1,size(values,2)
     write(*,value_format,advance='NO')values(i,j)
  enddo
  enddo
end subroutine


! ============================================================================
! checks the conservation of a substance and issues a message with specified
! severity if the difference is not within tolerance.
subroutine check_conservation(tag, substance, d1, d2, tolerance, severity)
  character(*), intent(in) :: tag ! message tag (subroutine name or some such)
  character(*), intent(in) :: substance ! name of the substance for printout
  real, intent(in) :: d1,d2 ! values to check
  real, intent(in) :: tolerance ! tolerance of the test
  integer, intent(in), optional :: severity ! severity of the non-conservation error:
         ! Can be WARNING, FATAL, or negative. Negative means check is not done.

  ! ---- local vars
  integer :: y,mo,d,h,m,s ! components of date
  integer :: thread, face
  real    :: lon, lat ! current coordinates, degree
  character(512) :: message
  integer :: severity_

  if(.not.do_check_conservation) return

  severity_=FATAL
  if (present(severity))severity_=severity

  if (severity_<0) return

  if (abs(d2-d1)<tolerance) then
     if (is_watch_point().and.watch_conservation) then
     write(*,'(3(x,a,g23.16))')&
          trim(tag)//': conservation of '//trim(substance)//'; before=', d1, 'after=', d2, 'diff=',d2-d1
     endif
  else
     thread = 1
!$   thread = OMP_GET_THREAD_NUM()+1
     call get_date(lnd%time,y,mo,d,h,m,s)
     call get_current_coordinates(thread, lon, lat, face)
     write(message,'(3(x,a,g23.16),2(x,a,f9.4),4(x,a,i4),x,a,i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))')&
          'conservation of '//trim(substance)//' is violated; before=', d1, 'after=', d2, 'diff=',d2-d1,&
          'at lon=',lon, 'lat=',lat, &
          'i=',curr_i(thread),'j=',curr_j(thread),'tile=',curr_k(thread),'face=',lnd%ug_face, &
          'time=',y,mo,d,h,m,s
     call error_mesg(tag,message,severity_)
  endif
end subroutine check_conservation

! returns coordinates and face of current point
subroutine get_current_coordinates(thread, lon, lat, face)
   integer, intent(in)  :: thread   ! our thread
   real, intent(out)    :: lon, lat ! coordinates of current point, degree
   integer, intent(out) :: face     ! number of cubic sphere face

   if(curr_l(thread) == 0) then
      lon = lnd%sg_lon(curr_i(thread),curr_j(thread))*180.0/PI
      lat = lnd%sg_lat(curr_i(thread),curr_j(thread))*180.0/PI
      face = lnd%sg_face
   else
      lon = lnd%ug_lon(curr_l(thread))*180.0/PI
      lat = lnd%ug_lat(curr_l(thread))*180.0/PI
      face = lnd%ug_face
   endif
   do while (lon > 180.0)
      lon = lon - 360.0
   enddo
   do while (lon < -180.0)
      lon = lon + 360.0
   enddo
end subroutine get_current_coordinates

! ============================================================================
! print a message with current coordinates and time
subroutine land_error_message(text,severity)
  character(*), intent(in) :: text
  integer, intent(in), optional :: severity

  integer :: y,mo,d,h,m,s ! components of date
  real    :: lon, lat ! current coordinates, degree
  integer :: thread, face
  character(512) :: message
  integer :: severity_

  severity_=WARNING
  if (present(severity))severity_=severity

  thread = 1
!$   thread = OMP_GET_THREAD_NUM()+1
  call get_date(lnd%time,y,mo,d,h,m,s)
  call get_current_coordinates(thread, lon, lat, face)
  write(message,'(2(x,a,f9.4),4(x,a,i4),x,a,i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))') &
       'at lon=',lon, 'lat=',lat, &
       'i=',curr_i(thread),'j=',curr_j(thread),'tile=',curr_k(thread),'face=',face, &
       'time=',y,mo,d,h,m,s
  call error_mesg(text,message,severity_)

end subroutine land_error_message

! ============================================================================
function string_from_time(time) result(str)
  character(19) :: str  ! YYYY-MM-DD HH:MM:YY
  type(time_type), intent(in) :: time

  integer :: y,mo,d,h,m,s ! components of date for debug printout

  call get_date(lnd%time,y,mo,d,h,m,s)
  write(str,'(i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))')y,mo,d,h,m,s
end function string_from_time

! ============================================================================
! print time in the debug output
subroutine log_date(tag,time)
  character(*),    intent(in) :: tag
  type(time_type), intent(in) :: time
  integer :: y,mo,d,h,m,s ! components of date for debug printout

  call get_date(lnd%time,y,mo,d,h,m,s)
  write(*,'(a,i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))') tag,y,mo,d,h,m,s
end subroutine log_date

end module land_debug_mod
