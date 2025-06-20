module transition_io_mod

use netcdf, only: nf90_max_name
use constants_mod, only : PI
use fms_mod, only : string, lowercase, error_mesg, FATAL, WARNING, NOTE

use time_manager_mod, only : time_type, set_date, valid_calendar_types, get_calendar_type, &
     operator(+), operator(-), operator(>), operator(<), operator(<=), operator(/), &
     operator(//), operator(==)
use time_interp_mod, only : time_interp ! used for irrigation
use get_cal_time_mod, only : get_cal_time
use horiz_interp_mod, only : horiz_interp_type, horiz_interp_new, horiz_interp_del
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

! ==== data types ===========================================================

!> container for information about input file and grid
!!
!! We assume that all variables from this file are on the same grid (horizontal and time),
!! that the valid values mask is the same does not change in time, and that the same
!! normalization factor (if any) must be applied to all of them
type :: infile_T
  character(1024) :: path      = '' !< file path
  character(1024) :: static    = '' !< static path
  character(16)   :: data_type = '' !< type of input data (LUH1 or LUH2). Due to differences
      !! in the normalization in the input data sets (per unit area of land or per unit area
      !! of grid cell) and differences in definition of valid values mask, interpolation is set
      !! up slightly differently depending on the type of the data set.
  type(FmsNetcdfFile_t)        :: statobj !< static fms2_io file object
  type(FmsNetcdfFile_t)        :: ncobj !< netcdf fms_io2 file object

  type(time_type), allocatable :: time_in(:)   !< input data time axis

  logical :: grid_initialized = .FALSE. !< set to TRUE when horizontal interpolator is set up
  integer                      :: nlon_in=-1, nlat_in=-1 !< sizes of input data horizontal grid
  real, allocatable            :: norm_in(:,:) !< normalizing factor to convert input data to
      !! units of [fractions of vegetated area per year]
  type(horiz_interp_type)      :: interp       !< horizontal interpolator
contains
  procedure :: init         => infile_init
  procedure :: destroy      => infile_destroy
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
  character(len=nf90_max_name) :: name  = '' !< internal name of the field
  integer       :: nvars = 0  !< number of variable ids
  character(len=nf90_max_name), dimension(:), allocatable :: varname !< names of the input fields
contains
  procedure :: addvar   => varset_add_var
  procedure :: descr    => varset_descr
  procedure :: destroy  => varset_destroy
  procedure :: get_data => varset_get_data
  procedure :: interpolate => varset_interpolate
end type varset_T

! ---- module variables
logical :: module_is_initialized = .FALSE.
integer :: ndims
integer, dimension(:), allocatable :: dimlens
character(len=nf90_max_name), dimension(:), allocatable :: dimnames

contains

subroutine transition_io_init()
  if(module_is_initialized) return

  call log_version(version, module_name, __FILE__)
  module_is_initialized = .TRUE.
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

  logical :: path_exists, static_exists

  this%path      = path
  this%static    = static
  this%data_type = data_type
  path_exists = open_file(this%ncobj, this%path, "read")
  if(.not. path_exists) call error_mesg('land_transition_io_infile_init', &
      trim(path)//'" could not be opened.', FATAL)
  static_exists = open_file(this%statobj, this%static, "read")
  if(trim(lowercase(this%data_type)) == 'luh2' .and. .not. static_exists) call &
      error_mesg('land_transition_io_infile_init', &
      trim(static)//'" could not be opened.', FATAL)
  ! get time axis
  call get_time_axis(this%ncobj,this%time_in)
end subroutine infile_init

! ============================================================================
!> free memory and close input files
subroutine infile_destroy(this)
  class(infile_T), intent(inout) :: this

  ! close input file
  call close_file(this%ncobj)
  ! deallocate timeline
  if (allocated(this%time_in)) deallocate(this%time_in)
  ! deallocate interpolator, if it exists
  call horiz_interp_del(this%interp)
  ! deallocate norm, if it exist
  if (allocated(this%norm_in)) deallocate(this%norm_in)
  this%nlon_in = -1; this%nlat_in = -1
  this%grid_initialized = .FALSE.
end subroutine infile_destroy

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
end subroutine setup_hgrid

! ==== varset members ========================================================

! ============================================================================
!> add variable to a variable set
subroutine varset_add_var(this,infile,varname)
   class(varset_T), intent(inout) :: this
   class(infile_T), target        :: infile   !< input file
   character(*),    intent(in)    :: varname  !< name of the variable in input file

   character(len=nf90_max_name), allocatable :: varname_(:)

   if (.not.associated(this%file)) then
      this%file => infile
   else if (.not.associated(this%file,infile)) then
      call error_mesg('transition_io_varset_add_var', 'variable set already associated with different file', FATAL)
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

   if (variable_exists(this%file%ncobj, varname)) then
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
!> read, aggregate, and interpolate (in horizontal dimensions) the set of transitions
subroutine varset_get_data(this,rec,frac)
   class(varset_T), intent(in) :: this
   integer, intent(in)  :: rec     !< 1-based index of time record to read from input fields
   real,    intent(out) :: frac(:) !< aggregated value of this variable set

   real, allocatable :: buff0(:,:),buff1(:,:)
   integer :: i

   frac = 0.0
   if (this%nvars == 0) return

   if (.not.associated(this%file)) call error_mesg('transition_io_varset_get_data', &
       'variable set "'//trim(this%name)//'" has no associated file', FATAL)

   allocate(buff0(this%file%nlon_in,this%file%nlat_in), &
            buff1(this%file%nlon_in,this%file%nlat_in)  )
   buff1 = 0.0
   do i = 1,this%nvars
      call read_data(this%file%ncobj, this%varname(i), buff0, unlim_dim_level=rec)
      buff1 = buff1 + buff0
   enddo
   call horiz_interp_ug(this%file%interp,buff1*this%file%norm_in,frac)
   deallocate(buff0,buff1)
end subroutine varset_get_data

! ============================================================================
!> interpolate variables that belong to the set in time
subroutine varset_interpolate(this, time, frac, interp)
  class(varset_T), intent(in)  :: this
  type(time_type), intent(in)  :: time    !< time to interpolate to
  real           , intent(out) :: frac(:) !< value of the aggregated and interpolated input fields, on unstructured grid
  character(*)   , intent(in), optional :: interp !< time interpolation method
    !! 'exact' interpolates linearly in time,
    !! 'before' takes data from the beginning of the time interval,
    !! 'after' takes data from the beginning of the time interval,
    !! Default is 'exact'

  real :: frac1(lnd%ls:lnd%le)
  real :: frac2(lnd%ls:lnd%le)
  integer :: i1,i2, n
  real :: w  ! time interpolation weight
  type(time_type) :: time_adjust
  character(16) :: interp_

  interp_ = 'exact'
  if (present(interp)) interp_ = interp

  if (.not.associated(this%file)) call error_mesg('transition_io',&
       'variable set "'//trim(this%name)//'" has no associated file', FATAL)

  n = size(this%file%time_in)
  time_adjust = time
  if (time_adjust<this%file%time_in(1)) time_adjust = this%file%time_in(1)
  if (time_adjust>this%file%time_in(n)) time_adjust = this%file%time_in(n)

  call time_interp(time_adjust, this%file%time_in, w, i1,i2)

  select case (trim(lowercase(interp_)))
  case('exact')
     call this%get_data(i1,frac1)
     call this%get_data(i2,frac2)

     frac = frac1*(1-w)+frac2*w
  case('before')
     call this%get_data(i1,frac)
  case('after')
     call this%get_data(i2,frac)
  case default
     call error_mesg('transition_io',&
         'time interpolation method "'//trim(interp)//'" is incorrect, must be "before", "after", or "exact"',FATAL)
  end select
end subroutine varset_interpolate

! ============================================================================
!> string representation of variable set
!! \returns a string representing all the variables included in this set
function varset_descr(this) result(str)
  character(:), allocatable :: str
  class(varset_T), intent(in) :: this

  character(len=nf90_max_name) :: varname
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
