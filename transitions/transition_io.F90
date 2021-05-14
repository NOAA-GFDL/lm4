module transition_io_mod

use constants_mod, only : PI
use mpp_mod, only : mpp_error, FATAL
use fms_mod, only : string, error_mesg, FATAL, WARNING, NOTE, &
     mpp_pe, lowercase, file_exist, close_file, &
     check_nml_error, stdlog, mpp_root_pe, fms_error_handler

use time_manager_mod, only : time_type, set_date, get_date, set_time, &
     operator(+), operator(-), operator(>), operator(<), operator(<=), operator(/), &
     operator(//), operator(==), days_in_year, print_date, increment_date, get_time, &
     valid_calendar_types, get_calendar_type
use get_cal_time_mod, only : get_cal_time
use horiz_interp_mod, only : horiz_interp_type, horiz_interp_init, &
     horiz_interp_new, horiz_interp_del
use nfu_mod, only : nfu_validtype, nfu_inq_var, nfu_get_dim_bounds, nfu_get_rec, &
     nfu_get_dim, nfu_get_var, nfu_get_valid_range, nfu_is_valid

use land_tile_io_mod, only : print_netcdf_error
use land_data_mod, only : lnd, log_version, horiz_interp_ug

implicit none
private

! ==== public interface =====================================================
public :: transition_io_init
public :: infile_T
public :: varset_T
! ==== end of public interface ==============================================

! usage:

! type(infile_t), pointer => f
! type(var_set_t) :: varset()
! f => open(path)
! do k1 = ...
! do k2 = ...
!     ....
!     input_tran(k1,k2)%add(f,fieldname)
! enddo
! perhaps make input_tran an associative array, to get rid of reallocation for irrigation?

! ==== module constants =====================================================
character(len=*), parameter :: module_name = 'transitions_io_mod'
#include "../shared/version_variable.inc"

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

! ==== data types ===========================================================

!> container for information about input file and grid information
!!
!! We assume that all variables from this file are on the same grid (horizontal and time),
!! that the valid values mask is the same does not change in time, and that the same
!! normalization factor (if any) must be applied to all of them
type :: infile_T
  character(1024) :: path      = '' !< file path
  character(1024) :: static    = '' !< static file path
  character(16)   :: data_type = '' !< type of input data (LUH1 or LUH2). Due to differences
      !! in the normalization in the input data sets (per unit area of land or per unit area
      !! of grid cell) and differences in definition of valid values mask, interpolation is set
      !! up slightly differently depending on the type of the data set.
  integer         :: ncid = -1 !< netcdf file ID; switch to io2 file structure in future

  type(time_type), allocatable :: time_in(:)   !< input data time axis

  logical :: grid_initialized = .FALSE. !< set to TRUE when horizontal interpolator is set up
  integer                      :: nlon_in=-1, nlat_in=-1 !< sizes of input data horizontal grid
  real, allocatable            :: norm_in(:,:) !< normalizing factor to convert input data to
      !! units of [fractions of vegetated area per year]
  type(horiz_interp_type)      :: interp       !< horizontal interpolator
contains
  procedure :: init         => infile_init
  procedure :: destroy      => infile_destroy
  procedure :: inq_var      => infile_inq_var
end type infile_T

!> structure that represents a set of variables
!!
!! Since in general there is no on-to-one correspondence between LUH database land use
!! types and land model land use types, some of the variables need to be aggregated
!! on input. For example, in LUH2 data base there are several types of agriculture, while
!! the model only has one.
!!
!! This structure holds  a set of variables that come from the same file,
!! and are added together to get one model field
type :: varset_T
  type(infile_T), pointer :: file => NULL() !< pointer to input file object
  character(NF_MAX_NAME) :: name  = '' !< internal name of the field
  integer       :: nvars = 0  !< number of variable ids
  character(NF_MAX_NAME), allocatable :: varname(:) !< names of the input fields
contains
  procedure :: addvar   => varset_add_var
  procedure :: descr    => varset_descr
  procedure :: destroy  => varset_destroy
  procedure :: get_data => varset_get_data
end type varset_T

contains

subroutine transition_io_init()
  call log_version(version, module_name, &
  __FILE__)
end subroutine transition_io_init

! ==== infile_T member functions =============================================

! ============================================================================
!> open input file and set up basic grid information
!! \returns file ID: index of file in the table of input files
!! \throws FATAL "file ... could not be opened"
!!
!! opens file; reads time axis from the file
subroutine infile_init(this, path, static, data_type)
  class(infile_T), intent(inout) :: this
  character(*),    intent(in)    :: path   !< file path
  character(*),    intent(in)    :: static !< static file path
  character(*),    intent(in)    :: data_type !< data type, LUH1 or LUH2

  integer :: ierr

  this%path      = path
  this%static    = static
  this%data_type = data_type
  ierr = nf_open(path,NF_NOWRITE,this%ncid)
  if(ierr/=NF_NOERR) call mpp_error(FATAL, &
      'file "'//trim(path)//'" could not be opened because '//nf_strerror(ierr), FATAL)
  ! get time axis
  call get_time_axis(this%ncid,this%time_in)
end subroutine infile_init

! ============================================================================
!> free memory and close input files
subroutine infile_destroy(this)
  class(infile_T), intent(inout) :: this

  integer :: ierr
  ! close input file
  if (this%ncid > 0) ierr = nf_close(this%ncid)
  ! deallocate timeline
  if (allocated(this%time_in)) deallocate(this%time_in)
  ! deallocate interpolator, if it exists
  call horiz_interp_del(this%interp)
  ! deallocate norm, if it exist
  if (allocated(this%norm_in)) deallocate(this%norm_in)
  this%nlon_in = -1; this%nlat_in = -1
  this%grid_initialized = .FALSE.
end subroutine infile_destroy

! ============================================================================
subroutine infile_inq_var(this, varname, found)
  class(infile_T), intent(inout) :: this
  character(*),    intent(in)    :: varname
  logical,         intent(out)   :: found

  integer :: ierr
  integer :: dimids(NF_MAX_VAR_DIMS), dimlens(NF_MAX_VAR_DIMS)

  ierr = nfu_inq_var(this%ncid, trim(varname), dimids=dimids, dimlens=dimlens)

  select case(ierr)
  case (NF_NOERR)
     found = .TRUE.
  case (NF_ENOTVAR)
     found = .FALSE.
!       call error_mesg('land_transitions_init',&
!            'field "'//trim(varname)//'" not found in file "'//trim(filename)//'"',&
!            NOTE)
     ! do nothing in this case, it is OK for only subset of variables
     ! to be present in the file
     return
  case default
     call mpp_error(FATAL,&
          'error initializing field "'//varname//&
          '" from file "'//trim(this%path)//'" : '//nf_strerror(ierr))
  end select
end subroutine infile_inq_var

! ==== end of infile_T member functions ======================================

! ============================================================================
!> read time axis from a file
subroutine get_time_axis(ncid, time_in)
  integer, intent(in) :: ncid
  type(time_type), allocatable :: time_in(:)

  integer :: timedim ! id of the record (time) dimension
  integer :: timevar ! id of the time variable
  character(len=NF_MAX_NAME) :: timename  ! name of the time variable
  character(len=256)         :: timeunits ! units ot time in the file
  character(len=32) :: calendar ! model calendar
  real, allocatable :: time(:)  ! real values of time coordinate
  integer :: i, nrec

  ! get the time axis
  __NF_ASRT__(nf_inq_unlimdim(ncid, timedim))
  __NF_ASRT__(nf_inq_dimlen(ncid, timedim, nrec))
  allocate(time(nrec), time_in(nrec))
  __NF_ASRT__(nfu_get_dim(ncid, timedim, time))
  ! get units of time
  __NF_ASRT__(nf_inq_dimname(ncid, timedim, timename))
  __NF_ASRT__(nf_inq_varid(ncid, timename, timevar))
  timeunits = ' '
  __NF_ASRT__(nf_get_att_text(ncid,timevar,'units',timeunits))
  ! get model calendar
  calendar=valid_calendar_types(get_calendar_type())

  ! loop through the time axis and get time_type values in time_in
  if (index(lowercase(timeunits),'calendar_year')>0) then
     do i = 1,size(time)
        time_in(i) = set_date(nint(time(i)),1,1,0,0,0) ! uses model calendar
     end do
  else
     do i = 1,size(time)
        time_in(i) = get_cal_time(time(i),timeunits,calendar)
     end do
  endif
  deallocate(time)
end subroutine get_time_axis

! ==================================================================
!> set up horizontal grid information for a file
subroutine setup_hgrid(this,varname)
  class(infile_T), intent(inout) :: this
  character(*),    intent(in)    :: varname !< name of a variable whose horizontal grid is
      !! used as for horizontal interpolator setup. It is assumed that all variables
      !! read from a file are on the same grid.

  real, allocatable :: lon_in(:,:),lat_in(:,:) ! horizontal grid of input data
  real, allocatable :: buffer_in(:,:) ! buffers for input data reading
  real, allocatable :: mask_in  (:,:) ! valid data mask on the input data grid
  integer :: dimids(NF_MAX_VAR_DIMS), dimlens(NF_MAX_VAR_DIMS)
  type(nfu_validtype) :: v ! valid values range
  integer :: ncid1 ! ID of static file
  integer :: ierr

  if (this%grid_initialized) return ! do nothing if grid is already set up
  ! TODO: possibly check that variable size is the same

  __NF_ASRT__(nfu_inq_var(this%ncid, trim(varname), dimids=dimids, dimlens=dimlens))
  this%nlon_in = dimlens(1); this%nlat_in=dimlens(2)

  ! get the boundaries of the horizontal axes
  allocate(lon_in(this%nlon_in+1,1), lat_in(1,this%nlat_in+1) )
  __NF_ASRT__(nfu_get_dim_bounds(this%ncid, dimids(1), lon_in(:,1)))
  __NF_ASRT__(nfu_get_dim_bounds(this%ncid, dimids(2), lat_in(1,:)))

  ! allocate temporary variables
  allocate(buffer_in    (this%nlon_in,this%nlat_in), &
           mask_in      (this%nlon_in,this%nlat_in), &
           this%norm_in (this%nlon_in,this%nlat_in)  )

  ! get the first record from variable and obtain the mask of valid data
  ! assume that valid mask does not change with time
  __NF_ASRT__(nfu_get_rec(this%ncid,varname,1,buffer_in))
  ! get the valid range for the variable
  __NF_ASRT__(nfu_get_valid_range(this%ncid,varname,v))
  ! get the mask
  where (nfu_is_valid(buffer_in,v))
     mask_in = 1
  elsewhere
     mask_in = 0
  end where

  ! calculate the normalizing factor to convert input data to units of
  ! [fraction of vegetated area per year]
  select case (trim(lowercase(this%data_type)))
  case ('luh1')
     ! LUH1 (CMIP5) data were converted on pre-processing
     this%norm_in = 1.0
  case ('luh2')
     ! read static file and calculate normalizing factor
     ! LUH2 data are in [fraction of cell area per year]
     if (trim(this%static)=='') call mpp_error(FATAL, &
          'using LUH2 data set, but static data file is not specified')
     ierr=nf_open(this%static,NF_NOWRITE,ncid1)
     if(ierr/=NF_NOERR) call error_mesg('land_transitions_init', &
          'using LUH2 data set, but static data file "'// &
          trim(this%static)//'" could not be opened because '//nf_strerror(ierr), FATAL)
     __NF_ASRT__(nfu_get_var(ncid1,'landfrac',buffer_in))
     where (buffer_in > 0.0)
        this%norm_in = 1.0/buffer_in
     elsewhere
        this%norm_in = 0.0
        mask_in = 0
     end where
     ierr = nf_close(ncid1)
  case default
     call error_mesg('land_transitions_init','unknown data_type "'&
                    //trim(this%data_type)//'", use "luh1" or "luh2"', FATAL)
  end select

  ! initialize horizontal interpolator
  call horiz_interp_new(this%interp, lon_in*PI/180,lat_in*PI/180, &
       lnd%sg_lonb, lnd%sg_latb, &
       interp_method='conservative',&
       mask_in=mask_in, is_latlon_in=.TRUE. )

  ! get rid of temporary allocated data
  deallocate(buffer_in, mask_in,lon_in,lat_in)

  this%grid_initialized = .TRUE.
end subroutine setup_hgrid

! ==== varset members ========================================================

! ============================================================================
!> add variable to a variable set
subroutine varset_add_var(this,infile,varname)
   class(varset_T), intent(inout) :: this
   class(infile_T), target        :: infile   !< input file
   character(*),    intent(in)    :: varname  !< name of the variable in input file

   character(NF_MAX_NAME), allocatable :: varname_(:)
   logical :: found

   if (.not.associated(this%file)) then
      this%file => infile
   else if (.not.associated(this%file,infile)) then
      call mpp_error(FATAL, 'variable set already associated with different file')
   endif

   ! allocate space for variable names on the first call
   if (.not.allocated(this%varname)) then
      allocate(this%varname(10))
      this%varname(:) = ''
   endif

   ! expand space for variable names if necessary
   if (this%nvars >= size(this%varname)) then
      ! make space for new variables
      allocate(varname_(size(this%varname)+10))
      varname_(:) = ''
      varname_(1:this%nvars) = this%varname(1:this%nvars)
      call move_alloc(varname_,this%varname)
   endif

   call this%file%inq_var(varname,found)
   if (found) then
      call error_mesg('land_transitions_init',&
           'adding field "'//trim(varname)//'" from file "'//trim(this%file%path)//'"'//&
           ' to transition "'//trim(this%name)//'"',&
           NOTE)
      this%nvars = this%nvars+1
      this%varname(this%nvars) = trim(varname)
      ! set up grid in the input file
      call setup_hgrid (this%file,varname)
   else
      call error_mesg('land_transitions_init',&
           'did not find field "'//trim(varname)//'" in file "'//trim(this%file%path)//'"'//&
           ' for transition "'//trim(this%name)//'"',&
           NOTE)
   endif
end subroutine varset_add_var

! ============================================================================
!> read, aggregate, and interpolate set of transitions
subroutine varset_get_data(this,rec,frac)
   class(varset_T), intent(in) :: this
   integer, intent(in) :: rec
   real, intent(out) :: frac(:)

   real, allocatable :: buff0(:,:),buff1(:,:)
   integer :: i

   frac = 0.0
   if (this%nvars == 0) return

   if (.not.associated(this%file)) call mpp_error(FATAL, &
       'variable set "'//trim(this%name)//'" has no associated file')

   allocate(buff0(this%file%nlon_in,this%file%nlat_in), &
            buff1(this%file%nlon_in,this%file%nlat_in)  )
   buff1 = 0.0
   do i = 1,this%nvars
      __NF_ASRT__(nfu_get_rec(this%file%ncid,this%varname(i),rec,buff0))
      buff1 = buff1 + buff0
   enddo
   call horiz_interp_ug(this%file%interp,buff1*this%file%norm_in,frac)
   deallocate(buff0,buff1)
end subroutine varset_get_data

! ============================================================================
!> string representation of variable set
!! \returns a string representing all the variables included in this set
function varset_descr(this) result(str)
  character(:), allocatable :: str
  class(varset_T), intent(in) :: this

  character(NF_MAX_NAME) :: varname
  integer :: i

  str = trim(this%name)//' = '
  if (this%nvars == 0) then
     str = str//'0'
  else
     do i = 1, this%nvars
        if (i==1) then
           str = str//trim(this%varname(i))
        else
           str = str//' + '//trim(this%varname(i))
        endif
     enddo
  endif
end function varset_descr

! ============================================================================
!> destructor
subroutine varset_destroy(this)
  class(varset_T), intent(inout) :: this
  if (allocated(this%varname)) deallocate(this%varname)
  this%file => NULL()
end subroutine varset_destroy
! ==== end of varset members ====================================================

end module transition_io_mod
