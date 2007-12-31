module static_vegn_mod

use time_manager_mod,   only : time_type, set_date
use get_cal_time_mod,   only : get_cal_time
use fms_mod,            only : write_version_number, error_mesg, FATAL, NOTE, &
     mpp_pe, file_exist, open_namelist_file, close_file, check_nml_error, stdlog, &
     mpp_root_pe
use time_interp_mod,    only : time_interp

use land_io_mod,        only : print_netcdf_error
use vegn_cohort_mod,    only : vegn_cohort_type
use cohort_io_mod,      only : read_cohort_data_i0d_fptr, read_cohort_data_r0d_fptr 


implicit none
private

! ==== public interface =====================================================
public :: static_vegn_init
public :: static_vegn_end

public :: static_vegn_override
! ==== end of public interface ==============================================

! ==== module constants =====================================================
character(len=*), parameter :: &
     module_name = 'static_vegn_mod', &
     version     = '$Id: vegn_static_override.F90,v 15.0.2.1 2007/09/16 22:14:24 slm Exp $', &
     tagname     = '$Name: omsk_2007_12 $'

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.
integer :: ncid ! netcdf id of the input file
type(time_type),allocatable :: time_line(:) ! time line of input data
type(time_type)             :: ts,te        ! beginning and end of time interval

! ---- namelist variables ---------------------------------------------------
logical :: use_static_veg = .FALSE.
character(len=512) :: input_file = & ! name of input file for static vegetation
     "INPUT/static_veg_data.nc"
character(len=10)  :: timeline   = 'normal' ! type of timeline ('normal' or 'loop')
integer, dimension(6) :: &
     start_loop = (/1,1,1,0,0,0/), & ! beginning of the time loop
     end_loop   = (/1,1,1,0,0,0/)    ! end of the time loop

namelist/static_veg_nml/use_static_veg,input_file,timeline,start_loop,end_loop

! ==== NetCDF declarations ==================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

contains

! ===========================================================================
subroutine static_vegn_init()

  ! ---- local vars
  integer :: unit, ierr, io, unlimdim, timelen, timeid
  integer :: i
  character(len=NF_MAX_NAME) :: timename ! name of the time variable
  real, allocatable          :: t(:)     ! temporary real timeline
  character(len=256)         :: units    ! units of time in the file
  character(len=256)         :: calendar ! calendar of the data

  if(module_is_initialized) return

  call write_version_number(version, tagname)

  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=static_veg_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'static_veg_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  
  if (mpp_pe() == mpp_root_pe()) then
     write (stdlog(), nml=static_veg_nml)
  endif

  if(use_static_veg) then

     ! SET UP LOOP BOUNDARIES
     ts = set_date(start_loop(1),start_loop(2),start_loop(3), start_loop(4),start_loop(5),start_loop(6))
     te = set_date(end_loop(1)  ,end_loop(2)  ,end_loop(3)  , end_loop(4)  ,end_loop(5)  ,end_loop(6)  )
     
     ! OPEN INPUT FILE
     __NF_ASRT__(nf_open(input_file,NF_NOWRITE,ncid))
     
     ! READ TIME AXIS DATA
     __NF_ASRT__(nf_inq_unlimdim( ncid, unlimdim ))
     __NF_ASRT__(nf_inq_dimname ( ncid, unlimdim, timename ))
     __NF_ASRT__(nf_inq_varid   ( ncid, timename, timeid ))
     __NF_ASRT__(nf_inq_dimlen( ncid, unlimdim, timelen ))
     allocate (time_line(timelen), t(timelen))
     __NF_ASRT__(nf_get_var_double (ncid, timeid, t ))
     
     ! GET UNITS OF THE TIME
     units = ' '
     __NF_ASRT__(nf_get_att_text(ncid, timeid,'units',units))
     
     ! GET CALENDAR OF THE DATA
     calendar = ' '
     ierr = nf_get_att_text(ncid, timeid, 'calendar',calendar)
     if(ierr/=NF_NOERR) &
          ierr = nf_get_att_text(ncid, timeid,'calendar_type',calendar)
     if(ierr/=NF_NOERR) &
          calendar='JULIAN' ! use model calendar? how to get the name of the model calendar?
       
     ! CONVERT TIME TO THE FMS TIME_TYPE AND STORE IT IN THE TIMELINE FOR THE
     ! DATA SET
     do i = 1, size(t)
        ! set the respective value in the timeline
        time_line(i) = get_cal_time(t(i),units,calendar)
     enddo
     
     deallocate(t)
  endif

  module_is_initialized = .true.

end subroutine static_vegn_init

! ===========================================================================
subroutine static_vegn_end()
  if(.not.use_static_veg) return;

  __NF_ASRT__(nf_close(ncid))

  deallocate(time_line)
  module_is_initialized = .false.

end subroutine static_vegn_end

! ===========================================================================
subroutine static_vegn_override (time)
  type(time_type), intent(in)    :: time

  ! ---- local vars 
  integer :: index1, index2 ! result of time interpolation (only index1 is used)
  real    :: weight         ! another result of time interp, not used

  if(.not.use_static_veg)return;

  !   time_interp to find out the index of the current time interval
  if (timeline == 'loop') then
     call time_interp(time, ts, te, time_line, weight, index1, index2, &
                      correct_leap_year_inconsistency=.true.)
  else if (timeline == 'normal') then
     call time_interp(time, time_line, weight, index1, index2)
  else
     call error_mesg(module_name,'timeline option "'//trim(timeline)// &
          '" is incorrect, use "normal" or "loop"', FATAL)
  endif

  ! read the data into cohort variables
  call read_cohort_data_i0d_fptr(ncid, 'species' , cohort_species_ptr , index1)
  call read_cohort_data_r0d_fptr(ncid, 'bl'      , cohort_bl_ptr      , index1)
  call read_cohort_data_r0d_fptr(ncid, 'blv'     , cohort_blv_ptr     , index1)
  call read_cohort_data_r0d_fptr(ncid, 'br'      , cohort_br_ptr      , index1)
  call read_cohort_data_r0d_fptr(ncid, 'bsw'     , cohort_bsw_ptr     , index1)
  call read_cohort_data_r0d_fptr(ncid, 'bwood'   , cohort_bwood_ptr   , index1)
  call read_cohort_data_r0d_fptr(ncid, 'bliving' , cohort_bliving_ptr , index1)
  call read_cohort_data_i0d_fptr(ncid, 'status'  , cohort_status_ptr  , index1)

  ! derived variables will be updated in update_land_bc_fast
end subroutine static_vegn_override

! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_COHORT_ACCESSOR(xtype,x)\
subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;\
end subroutine

DEFINE_COHORT_ACCESSOR(integer,species)
DEFINE_COHORT_ACCESSOR(real,bl)
DEFINE_COHORT_ACCESSOR(real,br)
DEFINE_COHORT_ACCESSOR(real,blv)
DEFINE_COHORT_ACCESSOR(real,bsw)
DEFINE_COHORT_ACCESSOR(real,bwood)
DEFINE_COHORT_ACCESSOR(real,bliving)
DEFINE_COHORT_ACCESSOR(integer,status)

end module static_vegn_mod
