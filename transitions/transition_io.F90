module transition_io_mod

use netcdf, only: nf90_max_name

use constants_mod, only : PI
use mpp_mod, only : mpp_error, FATAL
use mpp_domains_mod, only : mpp_pass_sg_to_ug
use fms_mod, only : string, lowercase, error_mesg, FATAL, WARNING, NOTE
use fms_io_mod, only : get_file_name

use time_manager_mod, only : time_type, set_date, get_date, set_time, &
     operator(+), operator(-), operator(>), operator(<), operator(<=), operator(/), &
     operator(//), operator(==), days_in_year, print_date, increment_date, get_time, &
     valid_calendar_types, get_calendar_type
use time_interp_mod, only : time_interp
use get_cal_time_mod, only : get_cal_time
use horiz_interp_mod, only : horiz_interp_type, horiz_interp_init, &
     horiz_interp_new, horiz_interp_del
use fms2_io_mod, only: FmsNetcdfFile_t, Valid_t, read_data, open_file, close_file, &
    get_valid, is_valid, variable_exists, get_variable_size, &
    get_unlimited_dimension_name, get_dimension_size, get_variable_attribute, &
    get_variable_dimension_names, get_variable_num_dimensions
use axis_utils2_mod, only: axis_edges

use land_data_mod, only : lnd, log_version, horiz_interp_ug
use land_debug_mod, only : set_current_point, is_watch_cell, &
     get_current_point, check_var_range, log_date

implicit none
private

! ==== public interface =====================================================
public :: transition_io_init
public :: infile_T
public :: varset_T
public :: new_infile_LUH1, new_infile_LUH2, new_infile_CS
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

!> container for information about input file and grid
!!
!! We assume that all variables from this file are on the same grid (horizontal and time),
!! that the valid values mask is the same does not change in time, and that the same
!! normalization factor (if any) must be applied to all of them
!!
!! This is not an abstract class, but many of its methods should never be called
type :: infile_t
  character(1024)       :: path = '' !< file path
  type(FmsNetcdfFile_t) :: ncobj     !< netcdf fms_io2 file object

  type(time_type), allocatable :: time_in(:)   !< input data time axis
  integer         :: nlon_in=-1, nlat_in=-1 !< sizes of input data horizontal grid

contains
  procedure :: setup_hgrid => infile_t_setup_hgrid ! set up a horizontal grid and interpolator for conversion to UG
  procedure :: var_exists  => infile_t_var_exists  ! returns TRUE if variable is found in the file
  procedure :: get_record  => infile_t_get_record  ! reads a single record of given variable, on variable's native grid
  procedure :: to_ug       => infile_t_to_ug       ! converts from variable native grid to model grid, on unstructured domain
  final     :: infile_t_destroy ! destructor
end type infile_t

!> container for information about input file on regular lat-lon grid
type, extends(infile_t) :: infile_latlon_t
  character(1024) :: static    = '' !< static file path
  type(FmsNetcdfFile_t) :: statobj  !< netcdf fms_io2 file object for static file

  character(16)   :: data_type = '' !< type of input data (LUH1 or LUH2). Due to differences
      !! in the normalization in the input data sets (per unit area of land or per unit area
      !! of grid cell) and differences in definition of valid values mask, interpolation is set
      !! up slightly differently depending on the type of the data set.

  logical :: grid_initialized = .FALSE. !< set to TRUE when horizontal interpolator is set up
  real, allocatable            :: norm_in(:,:) !< normalizing factor to convert input data to
      !! units of [fractions of vegetated area per year]
  type(horiz_interp_type)      :: interp       !< horizontal interpolator
contains
  final     :: infile_latlon_destroy
  procedure :: setup_hgrid => infile_latlon_setup_hgrid ! set up a horizontal grid and interpolator for conversion to UG
  procedure :: to_ug       => infile_latlon_to_ug
end type infile_latlon_t

!> container for information about input file on model native (cubic sphere)
type, extends(infile_t) :: infile_cs_t
contains
  procedure :: setup_hgrid => infile_cs_setup_hgrid
  procedure :: to_ug       => infile_cs_to_ug
  procedure :: get_record  => infile_cs_get_record
end type infile_cs_t

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
  class(infile_T), pointer :: file => NULL() !< pointer to input file object
  character(NF_MAX_NAME) :: name  = '' !< internal name of the field
  integer       :: nvars = 0  !< number of variable ids
  character(NF_MAX_NAME), allocatable :: varname(:) !< names of the input fields
contains
  procedure :: addvar   => varset_add_var
  procedure :: descr    => varset_descr
  procedure :: destroy  => varset_destroy
  procedure :: get_data => varset_get_data
  procedure :: integrate   => varset_integrate   ! integrate over time between t1 and t2
  procedure :: interpolate => varset_interpolate ! interpolate in time
end type varset_T

! ---- module variables
logical :: module_is_initialized = .FALSE.

contains

subroutine transition_io_init()
  if(module_is_initialized) return
  call log_version(version, module_name, &
  __FILE__)
  module_is_initialized = .TRUE.
end subroutine transition_io_init

! ==== infile_T member functions =============================================

! ============================================================================
!> open input file and set up basic grid information
!! \returns file ID: index of file in the table of input files
!! \throws FATAL "file ... could not be opened"
!!
!! opens file; reads time axis from the file
subroutine infile_t_open(this, path,static,data_type)
  class(infile_latlon_t), intent(inout) :: this
  character(*),    intent(in)    :: path   !< file path
  character(*),    intent(in)    :: static !< static file path
  character(*),    intent(in)    :: data_type !< data type, LUH1 or LUH2

  logical :: exists,static_exists

  this%path = path ! store path for future reference
  exists = open_file(this%ncobj, this%path, mode="read")
  if(.not.exists) call mpp_error(FATAL, &
      'file "'//trim(this%path)//'" could not be opened because it does not exist', FATAL)
  ! get the time axis from file
  static_exists = open_file(this%statobj, static, "read")
  write (*,*) 'This datatype:',trim(lowercase(this%data_type))
  if(trim(lowercase(this%data_type)) == 'luh2' .and. .not. static_exists) call &
      error_mesg('land_transition_io_infile_init', &
      trim(static)//'" could not be opened.', FATAL)
  call get_time_axis(this%ncobj,this%time_in)
end subroutine infile_t_open

! ============================================================================
!> destructor: free memory and close input files
subroutine infile_t_destroy(this)
  type(infile_t), intent(inout) :: this

  ! close input file
  call close_file(this%ncobj)
  ! deallocate timeline
  if (allocated(this%time_in)) deallocate(this%time_in)

  this%nlon_in = -1; this%nlat_in = -1
end subroutine infile_t_destroy

! ============================================================================
subroutine infile_t_setup_hgrid(this, varname)
  class(infile_t), intent(inout) :: this
  character(*),    intent(in)    :: varname

  call mpp_error(FATAL, 'infile_t_setup_hgrid should never be called')
end subroutine infile_t_setup_hgrid

! ============================================================================
logical function infile_t_var_exists(this, varname)
  class(infile_T), intent(inout) :: this
  character(*),    intent(in)    :: varname

  infile_t_var_exists = variable_exists(this%ncobj,trim(varname))
end function infile_t_var_exists

! ============================================================================
!> given variable name and record number, read the variable data for this record
subroutine infile_t_get_record(this, varname, rec, buff)
  class(infile_t), intent(in)  :: this !< file object
  character(*),    intent(in)  :: varname !< name of the variable
  integer,         intent(in)  :: rec  !< record number
  real,            intent(out) :: buff(:,:) !< buffer for output data

  !TODO: check buffer size
  call read_data(this%ncobj, varname, buff, unlim_dim_level=rec )
end subroutine infile_t_get_record

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!> convert input data from native grid to model grid, on unstructured domain
subroutine infile_t_to_ug(this,data2,data1)
  class(infile_t), intent(in)  :: this !< file object
  real,            intent(in)  :: data2(:,:) !< input 2D data
  real,            intent(out) :: data1(:) !< output data, on model's unstructured grid

  call mpp_error(FATAL, 'infile_t_to_ug should never be called')
end subroutine infile_t_to_ug
! ==== end of infile_T member functions ======================================

! ============================================================================
!> read time axis from a file
subroutine get_time_axis(ncobj, time_in)
  type(FmsNetcdfFile_t), intent(in) :: ncobj
  type(time_type), allocatable :: time_in(:)

  character(len=nf90_max_name) :: timename  ! name of the time variable
  character(len=256)         :: timeunits ! units ot time in the file
  character(len=32) :: calendar ! model calendar
  real, allocatable :: time(:)  ! real values of time coordinate
  integer :: i, nrec

  ! get the time axis
  call get_unlimited_dimension_name(ncobj, timename)
  call get_dimension_size(ncobj, timename, nrec)
  allocate(time(nrec), time_in(nrec))
  ! get units of time
  call read_data(ncobj, timename, time)
  timeunits = ' '
  call get_variable_attribute(ncobj, timename, "units", timeunits)
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

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function new_infile_LUH1(path) result(ptr)
  class(infile_latlon_t), pointer :: ptr
  character(*), intent(in) :: path

  allocate(ptr)
  !call infile_t_open(ptr,path,)
!   ptr%static    = ''
!   ptr%data_type = 'luh1'
end function new_infile_LUH1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function new_infile_LUH2(path,static) result(ptr)
  class(infile_latlon_t), pointer :: ptr
  character(*), intent(in) :: path
  character(*), intent(in) :: static

  allocate(ptr)
  call infile_t_open(ptr,path,static,'luh2')
  ptr%static    = static
  ptr%data_type = 'luh2'
end function new_infile_LUH2

! ==================================================================
!> set up horizontal grid information for a file
subroutine infile_latlon_setup_hgrid(this,varname)
  class(infile_latlon_t), intent(inout) :: this
  character(*),    intent(in)    :: varname !< name of a variable whose horizontal grid is
      !! used as for horizontal interpolator setup. It is assumed that all variables
      !! read from a file are on the same grid.

  real, allocatable :: lon_in(:,:),lat_in(:,:) ! horizontal grid of input data
  real, allocatable :: buffer_in(:,:) ! buffers for input data reading
  real, allocatable :: mask_in  (:,:) ! valid data mask on the input data grid
  character(len=256), allocatable :: dimnames(:)
  integer, allocatable :: dimlens(:)
  integer :: ndims

  type(Valid_t) :: v

  if (this%grid_initialized) return ! do nothing if grid is already set up
  ! TODO: possibly check that variable size is the same

  ndims = get_variable_num_dimensions(this%ncobj, varname)
  allocate(dimlens(ndims))
  call get_variable_size(this%ncobj, varname, dimlens)
  this%nlon_in = dimlens(1); this%nlat_in=dimlens(2)
  deallocate(dimlens)

  ! allocate temporary variables
  allocate(buffer_in    (this%nlon_in,this%nlat_in), &
           mask_in      (this%nlon_in,this%nlat_in), &
           this%norm_in (this%nlon_in,this%nlat_in)  )
  allocate(lon_in(this%nlon_in+1,1), lat_in(1,this%nlat_in+1) )

  ! get the boundaries of the horizontal axes
  allocate(dimnames(ndims))
  call get_variable_dimension_names(this%ncobj, varname, dimnames)
  call axis_edges(this%ncobj, dimnames(1), lon_in(:,1))
  call axis_edges(this%ncobj, dimnames(2), lat_in(1,:))
  deallocate(dimnames)

  ! get the first record from variable and obtain the mask of valid data
  ! assume that valid mask does not change with time
  call read_data(this%ncobj, varname, buffer_in, unlim_dim_level=1)
  ! get the valid range for the variable
  v = get_valid(this%ncobj, varname)
  ! get the mask
  where (is_valid(buffer_in,v))
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
     call read_data(this%statobj, 'landfrac', buffer_in)
     write(*,*) 'landfrac statobj is: ',this%statobj%path
     where (buffer_in > 0.0)
        this%norm_in = 1.0/buffer_in
     elsewhere
        this%norm_in = 0.0
        mask_in = 0
     end where
     call close_file(this%statobj)
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
end subroutine infile_latlon_setup_hgrid

subroutine infile_latlon_to_ug(this, data2, data1)
  class(infile_latlon_t), intent(in)  :: this !< file object
  real,            intent(in)  :: data2(:,:) !< input 2D data
  real,            intent(out) :: data1(:) !< output data, on model's unstructured grid
  call horiz_interp_ug(this%interp,data2*this%norm_in,data1)
end subroutine infile_latlon_to_ug


subroutine infile_latlon_destroy(this)
  type(infile_latlon_t), intent(inout) :: this

  ! deallocate interpolator, if it exists
  call horiz_interp_del(this%interp)
  ! deallocate norm, if it exist
  if (allocated(this%norm_in)) deallocate(this%norm_in)
end subroutine infile_latlon_destroy



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function new_infile_CS(path) result(ptr)
  class(infile_cs_t), pointer :: ptr
  character(*), intent(in) :: path

  logical :: read_dist, io_domain_exist, found_file

  allocate(ptr)

  found_file = get_file_name(path, ptr%path, read_dist, io_domain_exist, domain=lnd%sg_domain)
  if (.not.found_file) call mpp_error(FATAL, &
     'file "'//trim(path)//'" not found')

!   call infile_t_open(ptr,ptr%path)

!   ! data are suposed to be on SG grid compute domain
!   ptr%nlon_in   = lnd%ie-lnd%is+1
!   ptr%nlat_in   = lnd%je-lnd%js+1
end function new_infile_CS

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine infile_cs_setup_hgrid(this,varname)
  class(infile_cs_t), intent(inout) :: this
  character(*),       intent(in)    :: varname

  ! do nothing here
end subroutine infile_cs_setup_hgrid

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine infile_cs_to_ug(this,data2,data1)
  class(infile_cs_t), intent(in)  :: this !< file object
  real,            intent(in)  :: data2(:,:) !< input 2D data, on SG domain
  real,            intent(out) :: data1(:) !< output data, on model's unstructured grid

  call mpp_pass_SG_to_UG(lnd%ug_domain, data2, data1)
end subroutine infile_cs_to_ug

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine infile_cs_get_record(this,varname,rec,buff)
  class(infile_cs_t), intent(in)  :: this !< file object
  character(*),    intent(in)  :: varname !< name of the variable
  integer,         intent(in)  :: rec  !< record number
  real,            intent(out) :: buff(:,:) !< buffer for output data, supposed to be on SG domain

!   call read_data(this%path, varname, buff, lnd%sg_domain, rec)
! slm: I am not entirely sure this fms2_io call is equivalent to the older version above:
! see signature of the function below. Do we need to create FmsNetcdfDomainFile_t for CS files?
  call read_data(this%ncobj, varname, buff, unlim_dim_level=rec)

!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
! subroutine domain_read_2d(fileobj, variable_name, vdata, unlim_dim_level, &
!                           corner, edge_lengths)
!
!   type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
!   character(len=*), intent(in) :: variable_name !< Variable name.
!   class(*), dimension(:,:), intent(inout) :: vdata !< Data that will
!                                                    !! be written out
!                                                    !! to the netcdf file.
!   integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
!                                                    !! dimension.
!   integer, dimension(2), intent(in), optional :: corner !< Array of starting
!                                                         !! indices describing
!                                                         !! where the data
!                                                         !! will be written to.
!   integer, dimension(2), intent(in), optional :: edge_lengths !< The number of
!                                                               !! elements that
!                                                               !! will be written
!                                                               !! in each dimension.

end subroutine infile_cs_get_record




! ==== varset members ========================================================

! ============================================================================
!> add variable to a variable set
subroutine varset_add_var(this,infile,varname)
   class(varset_T), intent(inout) :: this
   class(infile_T), target        :: infile   !< input file
   character(*),    intent(in)    :: varname  !< name of the variable in input file

   character(NF_MAX_NAME), allocatable :: varname_(:)

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

   if (this%file%var_exists(varname)) then
      call error_mesg('land_transitions_init',&
           'adding field "'//trim(varname)//'" from file "'//trim(this%file%path)//'"'//&
           ' to transition "'//trim(this%name)//'"',&
           NOTE)
      this%nvars = this%nvars+1
      this%varname(this%nvars) = trim(varname)
      ! set up grid in the input file
      call this%file%setup_hgrid (varname)
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
      call this%file%get_record(this%varname(i),rec,buff0)
      buff1 = buff1 + buff0
   enddo
   call this%file%to_ug(buff1,frac)

   deallocate(buff0,buff1)
end subroutine varset_get_data

! ==============================================================================
! given boundaries of time interval [t1,t2], calculates total transition (time
! integral of transition rates) over the specified interval
subroutine varset_integrate(tran, t1, t2, frac)
  class(varset_T), intent(in)  :: tran
  type(time_type), intent(in)  :: t1,t2 ! time boundaries
  real           , intent(out) :: frac(:)

  ! ---- local vars
  integer :: n ! size of time axis
  type(time_type) :: ts,te
  integer         :: i1,i2
  real :: w  ! time interpolation weight
  real :: dt ! current time interval, in years
  real :: sum(size(frac(:)))
  integer :: l

  ! adjust the integration limits, in case they are out of range
  associate(time_in => tran%file%time_in)
  n = size(time_in)
  ts = t1;
  if (ts<time_in(1)) ts = time_in(1)
  if (ts>time_in(n)) ts = time_in(n)
  te = t2
  if (te<time_in(1)) te = time_in(1)
  if (te>time_in(n)) te = time_in(n)

  call time_interp(ts, time_in, w, i1,i2)
  call tran%get_data(i1,frac)

  dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
  sum = -frac*w*dt
  do while(time_in(i2)<=te)
     call tran%get_data(i1,frac)
     dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
     sum = sum+frac*dt
     i2 = i2+1
     i1 = i2-1
     if(i2>size(time_in)) exit ! from loop
  enddo

  call time_interp(te,time_in,w,i1,i2)
  call tran%get_data(i1,frac)
  dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
  frac = sum+frac*w*dt
  end associate
  ! check the transition rate validity
  do l = 1,size(frac(:))
     call set_current_point(l+lnd%ls-1,1)
     call check_var_range(frac(l),0.0,HUGE(1.0),'integral_transition',tran%name, FATAL)
  enddo
end subroutine varset_integrate

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine varset_interpolate(this, t2, irr_area, interp)
  class(varset_T), intent(in)  :: this ! id of the field
  type(time_type), intent(in)  :: t2 ! time boundaries
  real           , intent(out) :: irr_area(:)
  character(*)   , intent(in), optional :: interp ! 'exact', 'before', or 'after'

  real :: irr_area1(lnd%ls:lnd%le)
  real :: irr_area2(lnd%ls:lnd%le)
  integer :: i1,i2, n
  real :: w  ! time interpolation weight
  type(time_type) :: time_adjust
  character(16) :: interp_

  interp_ = 'exact'
  if (present(interp)) interp_ = interp

  if (.not.associated(this%file)) call mpp_error(FATAL, &
       'variable set "'//trim(this%name)//'" has no associated file')
  n = size(this%file%time_in)
  time_adjust = t2
  if (time_adjust<this%file%time_in(1)) time_adjust = this%file%time_in(1)
  if (time_adjust>this%file%time_in(n)) time_adjust = this%file%time_in(n)

  call time_interp(time_adjust, this%file%time_in, w, i1,i2)

  select case (trim(lowercase(interp_)))
  case('exact')
     call this%get_data(i1,irr_area1)
     call this%get_data(i2,irr_area2)

     irr_area = irr_area1*(1-w)+irr_area2*w
  case('before')
     call this%get_data(i1,irr_area)
  case('after')
     call this%get_data(i2,irr_area)
  case default
     call mpp_error(FATAL, &
         'interpolation type "'//trim(interp)//'" is incorrect, mus be "before", "after", or "exact"')
  end select
end subroutine varset_interpolate

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
