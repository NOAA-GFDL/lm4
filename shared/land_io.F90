module land_io_mod

use netcdf, only: nf90_max_name
use mpp_domains_mod, only : mpp_pass_sg_to_ug

use constants_mod,     only : PI
use fms_mod, only: error_mesg, FATAL, stdlog, mpp_pe, &
     mpp_root_pe, string, check_nml_error, input_nml_file
use mpp_io_mod, only: axistype, mpp_get_axis_data
use axis_utils_mod, only : get_axis_bounds
use horiz_interp_mod,  only : horiz_interp_type, &
     horiz_interp_new, horiz_interp_del, horiz_interp
use land_numerics_mod, only : nearest, bisect
use land_data_mod, only : log_version, lnd, horiz_interp_ug
use time_interp_external_mod, only: time_interp_external_init, &
     time_interp_external, init_external_field
use time_manager_mod, only: time_type
use mpp_domains_mod, only : domain2d

use fms2_io_mod, only: close_file, FmsNetcdfFile_t, get_valid, get_variable_attribute, &
                       get_variable_num_dimensions, get_variable_dimension_names, get_variable_size, &
                       is_valid, open_file, read_data, Valid_t, variable_att_exists, variable_exists
use legacy_mod, only: axis_edges

implicit none
private

! ==== public interface ======================================================
public :: init_cover_field
public :: read_field
public :: read_land_io_namelist
public :: external_ts_type
public :: init_external_ts, del_external_ts
public :: read_external_ts
public :: input_buf_size

! ==== end of public interface ===============================================

interface read_field
   module procedure read_field_N_2D, read_field_N_3D
   module procedure read_field_N_2D_int, read_field_N_3D_int
end interface

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_io_mod'
#include "../shared/version_variable.inc"

real, parameter :: DEFAULT_FILL_INT  = -HUGE(1)
real, parameter :: DEFAULT_FILL_REAL = -HUGE(1.0)

! ==== module types ==========================================================
type external_ts_type
   character(256) :: filename
   character(64)  :: fieldname
   integer :: id ! ID of external field
   type(horiz_interp_type) :: interp ! interpolator
   real :: fill ! fill value for missing data
end type external_ts_type

! ==== module data ===========================================================
logical :: module_is_initialized = .false.
character(len=64)  :: interp_method = "conservative"
integer, protected :: input_buf_size = 65536 ! input buffer size for tile and cohort reading
namelist /land_io_nml/ interp_method, input_buf_size
contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine read_land_io_namelist()
  integer :: io, ierr, unit


  module_is_initialized = .TRUE.

  call log_version (version, module_name, &
  __FILE__)

  read (input_nml_file, nml=land_io_nml, iostat=io)
  ierr = check_nml_error(io, 'land_io_nml')

  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_io_nml)
  endif

  if(trim(interp_method) .NE. "conservative" .AND. trim(interp_method) .NE. "conserve_great_circle") then
     call error_mesg ( module_name,'interp_method should be "conservative" or "conserve_great_circle"', FATAL)
  endif

  if (input_buf_size <= 0) then
     call error_mesg ( module_name,'input_buf_size must be larger than zero', FATAL)
  endif

end subroutine read_land_io_namelist


! ============================================================================
! This procedure creates and initializes a field of fractional coverage.
subroutine init_cover_field( &
     cover_to_use, filename, cover_field_name, frac_field_name, &
     lonb, latb, uniform_cover, input_cover_types, frac)
  character(len=*), intent(in) :: cover_to_use
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: cover_field_name, frac_field_name
  real            , intent(in) :: lonb(:,:), latb(:,:) ! boundaries of the grid cells
  integer         , intent(in) :: uniform_cover
  integer         , intent(in) :: input_cover_types(:)
  real            , intent(out):: frac(:,:) ! output-global map of soil fractional coverage

  ! ---- local vars ---------------------------------------------------------
  integer :: l, k     ! iterators
  integer :: cover_id
  real    :: maxfrac, total

  if( .not. module_is_initialized ) &
       call error_mesg(module_name,'land_io_init is not called', FATAL)

  frac = 0

  if (cover_to_use == 'multi-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
  else if (cover_to_use=='single-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
     do l = 1,size(frac,1)
        total = sum(frac(l,:))
        if (total <= 0) cycle ! operate on valid input data points only
        maxfrac=0 ; cover_id=1
        do k = 1,size(frac,2)
           if(frac(l,k).gt.maxfrac) then
              maxfrac=frac(l,k)
              cover_id=k
           endif
        enddo
        ! set all fractions except dominant fraction to zero
        frac(l,:) = 0.0
        frac(l,cover_id) = total
     enddo
  else if (cover_to_use == 'uniform') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
     do l = 1,size(frac,1)
        total = sum(frac(l,:))
        if (total <= 0) cycle ! operate on valid input data points only
        ! set all fractions except dominant fraction to zero
        frac(l,:) = 0.0
        frac(l,uniform_cover) = total
     enddo
  else
     call error_mesg ( module_name,'illegal value of cover_to_use '//cover_to_use, FATAL )
  endif

end subroutine init_cover_field


! ============================================================================
subroutine read_cover_field(file, cover_field_name, frac_field_name,&
     lonb, latb, input_cover_types, frac)
  character(len=*)  , intent(in)  :: file            ! file to read from
  character(len=*)  , intent(in)  :: cover_field_name, frac_field_name
  real              , intent(in)  :: lonb(:,:),latb(:,:) ! boundaries of the model grid
  real              , intent(out) :: frac(:,:)     ! resulting fractions
  integer, optional , intent(in)  :: input_cover_types(:)

  ! --- local vars
  type(FmsNetcdfFile_t) :: fileobj
  logical :: exists

  exists = open_file(fileobj, trim(file), "read")
  if (.not. exists) then
    call error_mesg(module_name, 'input file "'//trim(file)//'" does not exist', &
                    fatal)
  endif

! If field named 'cover' does not exist in file then read field named 'frac'
! 'cover' does not exist in either ground_type.nc or cover_type.nc
! The extent of the third dimension is 10 in ground_type.nc and 11 in cover_type.nc
 
  if (variable_exists(fileobj, cover_field_name)) then
    call do_read_cover_field(fileobj, cover_field_name, lonb, latb, input_cover_types, frac)
  elseif (variable_exists(fileobj, frac_field_name)) then
    call do_read_fraction_field(fileobj, frac_field_name, lonb, latb, input_cover_types, frac)
  else
    call error_mesg(module_name, &
                    'neither "'//trim(cover_field_name)//'" nor "'//&
                    frac_field_name//'" is present in input file "'//trim(file)//'"' , &
                    fatal)
  endif
  call close_file(fileobj)

end subroutine read_cover_field

! ============================================================================
subroutine do_read_cover_field(fileobj, name, lonb, latb, input_cover_types, frac)

  type(FmsNetcdfFile_t), intent(in)  :: fileobj
  character(len=*), intent(in) :: name

  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:)

  ! ---- local vars
  integer, dimension(:), allocatable :: dimlens
  character(len=nf90_max_name), dimension(:), allocatable :: dimnames
  integer :: ndims
  integer :: nlon, nlat, k
  integer, allocatable :: in_cover(:,:)
  real, allocatable    :: in_lonb(:), in_latb(:), x(:,:), r_in_cover(:,:)
  type(horiz_interp_type) :: interp
  type(Valid_t) :: v
  integer :: in_j_start, in_j_end, in_j_count ! limits of the latitude belt we read
  integer :: start(4), count(4)
  real :: min_in_latb, max_in_latb, y

  ! check the field dimenions and get size of the longitude and latitude axes
  ndims = get_variable_num_dimensions(fileobj, name)
  if (ndims .ne. 2) then
    call error_mesg('do_read_cover_field', &
                    'cover field "'//trim(name)//'" in file "'//trim(fileobj%path)// &
                    '" must be two-dimensional (lon,lat)', fatal)
  endif
  allocate(dimlens(ndims))
  allocate(dimnames(ndims))
  call get_variable_size(fileobj, name, dimlens)
  call get_variable_dimension_names(fileobj, name, dimnames)
  nlon = dimlens(1)
  nlat = dimlens(2)
  allocate(in_lonb(nlon+1), in_latb(nlat+1))
  call axis_edges(fileobj, dimnames(1), in_lonb)
  in_lonb = in_lonb*PI/180
  call axis_edges(fileobj, dimnames(2), in_latb)
  in_latb = in_latb*PI/180
  deallocate(dimlens)
  deallocate(dimnames)

  ! to minimize the i/o and work done by horiz_interp, find the boundaries
  ! of latitude belt in input data that covers the entire latb array
  min_in_latb = minval(in_latb); max_in_latb = maxval(in_latb)
  y = minval(latb)
  if (y<min_in_latb) then
     in_j_start = 1
  else if (y>max_in_latb) then
     in_j_start = nlat
  else
     in_j_start=bisect(in_latb, y)
  endif

  y = maxval(latb)
  if (y<min_in_latb) then
     in_j_end = 1
  else if (y>max_in_latb) then
     in_j_end = nlat
  else
     in_j_end = bisect(in_latb, y)
  endif
  in_j_count = in_j_end - in_j_start + 1

  ! check for unreasonable values
  if (in_j_start<1) &
     call error_mesg('do_read_cover_field','reading field "'//trim(name)//'" from file "'&
                     //trim(fileobj%path)//'" input latitude start index ('&
                     //trim(string(in_j_start))//') is out of bounds', FATAL)
  if (in_j_count<1) &
     call error_mesg('do_read_cover_field','reading field "'//trim(name)//'" from file "'&
                     //trim(fileobj%path)//'" computed input latitude count for domain'&
                     //' is not positive, perhaps input data do not cover entire globe', FATAL)
  if (in_j_start+in_j_count-1>nlat) &
     call error_mesg('do_read_cover_field','reading field "'//trim(name)//'" from file "'&
                     //trim(fileobj%path)//'input latitude count ('&
                     //trim(string(in_j_count))//') is too large (start index='&
                     //trim(string(in_j_start))//')', FATAL)

  ! allocate input data buffers
  allocate ( x(nlon,in_j_count), in_cover(nlon,in_j_count), r_in_cover(nlon,in_j_count) )

  ! read input data
  start = (/1,in_j_start,1,1/)
  count(1:2) = shape(in_cover)
  count(3:4) = 1
  call read_data(fileobj, name, r_in_cover, corner=start, edge_lengths=count)
  in_cover = r_in_cover
  v = get_valid(fileobj, name)

  call horiz_interp_new(interp, in_lonb,in_latb(in_j_start:in_j_start+in_j_count), &
       lonb,latb, interp_method=trim(interp_method))
  frac=0
  do k = 1,size(input_cover_types(:))
     x=0
     where(is_valid(r_in_cover,v).and.in_cover==input_cover_types(k)) x = 1
     call horiz_interp_ug(interp,x,frac(:,k))
  enddo

  call horiz_interp_del(interp)

  ! clean up memory
  deallocate(in_lonb, in_latb, in_cover, x)

end subroutine do_read_cover_field


! ============================================================================
 subroutine do_read_fraction_field(fileobj, name, lonb, latb, input_cover_types, frac)

  type(FmsNetcdfFile_t), intent(in)  :: fileobj
  character(len=*), intent(in) :: name
  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:)

  ! ---- local vars
  integer, dimension(:), allocatable :: dimlens
  character(len=nf90_max_name), dimension(:), allocatable :: dimnames
  integer :: ndims
  integer :: nlon, nlat, ntypes, k, cover
  real, allocatable :: in_frac(:,:,:)
  real, allocatable :: in_lonb(:), in_latb(:)
  real, allocatable :: in_mask(:,:)
  type(horiz_interp_type) :: interp
  type(Valid_t) :: v
  integer :: in_j_end, in_j_start, in_j_count ! limits of the latitude belt we read
  integer :: start(4), count(4)
  real :: min_in_latb, max_in_latb, y

  ! check the field dimenions and get size of the longitude and latitude axes
  ndims = get_variable_num_dimensions(fileobj, name)
  if (ndims .ne. 3) then
    call error_mesg('do_read_cover_field', &
                    'cover field "'//trim(name)//'" in file "'//trim(fileobj%path)// &
                    '" must be two-dimensional (lon,lat,_)', fatal)
  endif
  allocate(dimlens(ndims))
  allocate(dimnames(ndims))
  call get_variable_size(fileobj, name, dimlens)
  call get_variable_dimension_names(fileobj, name, dimnames)
  nlon = dimlens(1)
  nlat = dimlens(2)
  ntypes = dimlens(3)
  allocate(in_lonb(nlon+1), in_latb(nlat+1))
  call axis_edges(fileobj, dimnames(1), in_lonb)
  in_lonb = in_lonb*PI/180
  call axis_edges(fileobj, dimnames(2), in_latb)
  in_latb = in_latb*PI/180
  deallocate(dimlens)
  deallocate(dimnames)

  ! find the boundaries of latitude belt in input data that covers the
  ! entire latb array
  min_in_latb = minval(in_latb); max_in_latb = maxval(in_latb)
  y = minval(latb)
  if (y<min_in_latb) then
     in_j_start = 1
  else if (y>max_in_latb) then
     in_j_start = nlat
  else
     in_j_start=bisect(in_latb, y)
  endif

  y = maxval(latb)
  if (y<min_in_latb) then
     in_j_end = 1
  else if (y>max_in_latb) then
     in_j_end = size(in_latb)-1
  else
     in_j_end = bisect(in_latb, y)
  endif
  in_j_count = in_j_end - in_j_start + 1

  ! check for unreasonable values
if (in_j_start<1) &
     call error_mesg('do_read_fraction_field','reading field "'//trim(name)//'" from file "'&
                     //trim(fileobj%path)//'input latitude start index ('&
                     //trim(string(in_j_start))//') is out of bounds', FATAL)
  if (in_j_count<1) &
     call error_mesg('do_read_fraction_field','reading field "'//trim(name)//'" from file "'&
                     //trim(fileobj%path)//'computed input latitude count for domain'&
                     //' is not positive, perhaps input data do not cover entire globe', FATAL)
  if (in_j_start+in_j_count-1>nlat) &
     call error_mesg('do_read_fraction_field','reading field "'//trim(name)//'" from file "'&
                     //trim(fileobj%path)//'input latitude count ('&
                     //trim(string(in_j_count))//') is too large (start index='&
                     //trim(string(in_j_start))//')', FATAL)

  allocate( in_mask(nlon,in_j_count), in_frac(nlon,in_j_count,ntypes) )

  ! read input data
! iret = nf_get_vara_double(ncid, varid, (/1,in_j_start,1/), (/nlon,in_j_count,ntypes/), in_frac)
  start = (/1,in_j_start,1,1/)
  count(1:3) = shape(in_frac)
  count(4) = 1
  call read_data(fileobj, name, in_frac, corner=start, edge_lengths=count)
  v = get_valid(fileobj, name)


  ! Initialize horizontal interpolator; we assume that the valid data mask is
  ! the same for all levels in input frac array. This is probably a good assumption
  ! in all cases.
  where(is_valid(in_frac(:,:,1),v))
     in_mask = 1.0
  elsewhere
     in_mask = 0.0
  end where
  call horiz_interp_new(interp, &
       in_lonb,in_latb(in_j_start:in_j_start+in_j_count), lonb,latb,&
       interp_method=trim(interp_method), mask_in=in_mask)

  frac = 0
  do k = 1,size(input_cover_types)
     cover = input_cover_types(k)
     if (cover<1.or.cover>ntypes) then
        cycle ! skip all invalid indices in the array of input cover types
     endif

     call horiz_interp_ug(interp,in_frac(:,:,cover),frac(:,k))
  enddo

  ! clean up memory
  call horiz_interp_del(interp)
  deallocate(in_lonb, in_latb, in_frac, in_mask)

end subroutine do_read_fraction_field

! ============================================================================
subroutine read_field_N_2D_int(fileobj, varname, data_ug, interp, fill)
  type(FmsNetcdfFile_t), intent(in) :: fileobj
  character(*), intent(in)  :: varname
  integer,      intent(out) :: data_ug(:)
  character(*), intent(in), optional :: interp  ! kind of interpolation
  integer,      intent(in), optional :: fill    ! fill value for missing pints
  ! ---- local vars
  real :: data3(size(data_ug,1),1)
  real :: fill_

  fill_ = DEFAULT_FILL_INT
  if (present(fill)) fill_ = fill

  call read_field_N_3D(fileobj, varname, data3, interp, fill_)
  data_ug = nint(data3(:,1))
end subroutine read_field_N_2D_int

! ============================================================================
subroutine read_field_N_3D_int(fileobj, varname, data_ug, interp, fill)
  type(FmsNetcdfFile_t), intent(in) :: fileobj
  character(*), intent(in)  :: varname
  integer,      intent(out) :: data_ug(:,:)
  character(*), intent(in), optional :: interp
  integer,      intent(in), optional :: fill
  ! ---- local vars
  real :: data3(size(data_ug,1),size(data_ug,2))
  real :: fill_

  fill_ = DEFAULT_FILL_INT
  if (present(fill)) fill_ = fill

  call read_field_N_3D(fileobj, varname, data3, interp, fill_)
  data_ug = nint(data3(:,:))
end subroutine read_field_N_3D_int

! ============================================================================
subroutine read_field_N_2D(fileobj, varname, data_ug, interp, fill)
  type(FmsNetcdfFile_t), intent(in) :: fileobj
  character(*), intent(in)  :: varname
  real,         intent(out) :: data_ug(:)
  character(*), intent(in),  optional :: interp
  real,         intent(in), optional :: fill
  ! ---- local vars
  real    :: data3(size(data_ug,1),1)

  call read_field_N_3D(fileobj, varname, data3, interp, fill)
  data_ug = data3(:,1)
end subroutine read_field_N_2D

! ============================================================================
subroutine read_field_N_3D(fileobj, varname, data_ug, interp, fill)
  type(FmsNetcdfFile_t), intent(in) :: fileobj
  character(*), intent(in)  :: varname
  real,         intent(out) :: data_ug(lnd%ls:,:)
  character(*), intent(in), optional :: interp
  real,         intent(in), optional :: fill
 ! ---- local vars
  character(len=20) :: interp_
  real    :: fill_
  integer, dimension(:), allocatable :: dimlens
  integer :: varndims ! number of variable dimension
  integer :: nlon, nlat, nlev ! size of input grid
  real,    allocatable :: in_lonb(:), in_latb(:), in_lon(:), in_lat(:)
  character(len=nf90_max_name), dimension(:), allocatable :: dimnames
  type(Valid_t) :: v
  real,    allocatable :: in_data(:,:,:) ! input buffer
  logical, allocatable :: lmask(:,:,:) ! mask of valid input values
  real,    allocatable :: rmask(:,:,:) ! real mask for interpolator
  real,    allocatable :: data_sg(:,:,:) ! data on structured grid
  real,    allocatable :: omask(:,:) ! mask of valid output data
  real,    allocatable :: data2(:,:)
  integer :: k,imap,jmap,l
  type(horiz_interp_type) :: hinterp
  integer :: ndim,nvar,natt,nrec
  integer :: jstart, jend

  interp_ = 'bilinear'
  if(present(interp)) interp_ = interp

  fill_ = DEFAULT_FILL_REAL
  if (present(fill)) fill_=fill

  if (.not. variable_exists(fileobj, varname)) then
    call error_mesg('read_field', &
                    'variable "'//trim(varname)//'" not found in file "'//trim(fileobj%path)//'"', &
                    FATAL)
  endif

  ! get the dimensions of our variable
  varndims = get_variable_num_dimensions(fileobj, varname)
  if (varndims .lt. 2 .or. varndims .gt. 3) then
     call error_mesg('read_field','variable "'//trim(varname)//'" in file "'//trim(fileobj%path)// &
          '" is '//string(varndims)//'D, but only reading 2D or 3D variables is supported', FATAL)
  endif
  allocate(dimlens(varndims))
  call get_variable_size(fileobj, varname, dimlens)
  nlon = dimlens(1) ; nlat = dimlens(2)
  nlev = 1;
  if (varndims==3) nlev=dimlens(3)
  if(nlev/=size(data_ug,2)) then
     call error_mesg('read_field','3rd dimension length of the variable "'&
          //trim(varname)//'" ('//trim(string(nlev))//') in file "'//trim(fileobj%path)//&
          '" is different from the expected size of data ('// trim(string(size(data_ug,2)))//')', &
          FATAL)
  endif
  deallocate(dimlens)
  v = get_valid(fileobj, varname)

  ! read boundaries of the grid cells in longitudinal direction
  allocate (in_lon  (nlon),   in_lat  (nlat),  &
            in_lonb (nlon+1), in_latb (nlat+1) )
  allocate(dimnames(varndims))
  call get_variable_dimension_names(fileobj, varname, dimnames)
  call read_data(fileobj, dimnames(1), in_lon)
  in_lon = in_lon*PI/180
  call axis_edges(fileobj, dimnames(1), in_lonb)
  in_lonb = in_lonb*PI/180
  call read_data(fileobj, dimnames(2), in_lat)
  in_lat = in_lat*PI/180
  call axis_edges(fileobj, dimnames(2), in_latb)
  in_latb = in_latb*PI/180
  deallocate(dimnames)

  select case(trim(interp_))
  case('nearest')
     allocate (in_data(nlon, nlat, nlev), lmask(nlon, nlat, nlev))
     ! read input data. In case of nearest interpolation we need global fields.
     call read_data(fileobj, varname, in_data)
     lmask = is_valid(in_data,v)
     do k = 1,size(data_ug,2)
        do l = lnd%ls,lnd%le
           call nearest (lmask(:,:,k), in_lon, in_lat, lnd%ug_lon(l), lnd%ug_lat(l), imap, jmap)
           if(imap<0 .or. jmap<0) then
              call error_mesg('read_field', 'imap or jamp is negative' ,FATAL)
           endif
           data_ug(l,k) = in_data(imap,jmap,k)
        enddo
     enddo
     deallocate(in_data,lmask)
  case('bilinear')
     ! we do bilinear interpolation directly on unstructured grid
     call jlimits(minval(lnd%ug_lat),maxval(lnd%ug_lat),in_lat,jstart,jend)
     call read_input()
     call horiz_interp_new(hinterp, in_lonb, in_latb(jstart:jend+1), &
            reshape(lnd%ug_lon,[lnd%le-lnd%ls+1,1]), & ! reshape converts 1D array (N) to array of shape (N,1)
            reshape(lnd%ug_lat,[lnd%le-lnd%ls+1,1]), &
            interp_method='bilinear')
     allocate(omask(size(data_ug,1),1), data2(size(data_ug,1),1))
     do k = 1,size(data_ug,2)
        call horiz_interp(hinterp, in_data(:,:,k), data2(:,:), mask_in=rmask(:,:,k), mask_out=omask(:,:))
        data_ug(:,k) = data2(:,1)
        where (omask(:,1)==0.0) data_ug(:,k) = fill_
     enddo
     call horiz_interp_del(hinterp)
     deallocate(in_data, rmask, omask, data2)
  case('conservative')
     ! conservative interpolation is done on structured grid, and then interpolated values
     ! are passed to unstructured grid
     call jlimits(minval(lnd%sg_latb),maxval(lnd%sg_latb),in_lat,jstart,jend)
     call read_input() ! allocates and fills in_data, rmask
     ! we create horiz interpolator inside the loop, because data masks may be different for
     ! different levels
     allocate (data_sg(lnd%is:lnd%ie,lnd%js:lnd%je,size(data_ug,2))  ) ! data on structured grid
     do k = 1,size(data_ug,2)
        call horiz_interp_new(hinterp, in_lonb, in_latb(jstart:jend+1), lnd%sg_lonb, lnd%sg_latb, &
             mask_in=rmask(:,:,k), interp_method='conservative')
        data_sg(:,:,k) = fill_
        call horiz_interp(hinterp,in_data(:,:,k),data_sg(:,:,k))
        call horiz_interp_del(hinterp)
     enddo
     call mpp_pass_sg_to_ug(lnd%ug_domain, data_sg, data_ug)
     deallocate(in_data, rmask, data_sg)
  case default
     call error_mesg('read_field','Unknown interpolation method "'//trim(interp_)//'". use "nearest", "bilinear", or "conservative"', FATAL)
  end select

  deallocate(in_lonb, in_latb, in_lon, in_lat)

  contains ! internal subroutines

  subroutine read_input
    ! Note that it changes data in host subroutine, so it must be internal
    ! ---- local vars
    integer :: start(4), count(4)

    allocate(in_data(nlon,jstart:jend,nlev), rmask(nlon,jstart:jend,nlev))
    start(1)  = 1;       count(1)  = nlon
    start(2)  = jstart;  count(2)  = jend-jstart+1
    start(3)  = 1;       count(3)  = nlev
    start(4:) = 1;       count(4:) = 1
    ! read input data
    call read_data(fileobj, varname, in_data, corner=start, edge_lengths=count)
    where (is_valid(in_data,v))
       rmask = 1.0
    elsewhere
       rmask = 0.0
    end where
  end subroutine read_input

  subroutine jlimits(minlat, maxlat, in_lat, jstart, jend)
    ! to minimize memory footprint, find the latitudinal boundaries in input
    ! data grid that cover our domain.

    ! This subroutine doesn't have to be internal, but it is not used (and perhaps not
    ! useful) anywhere else
    real,    intent(in)  :: minlat, maxlat
    real,    intent(in)  :: in_lat(:)
    integer, intent(out) :: jstart,jend

    integer :: nlat, j
    nlat = size(in_lat)

    jstart = 1; jend = nlat
    do j = 1, nlat
       if(minlat < in_lat(j)) then
          jstart = j-1
          exit
       endif
    enddo
    jstart = max(jstart-1,1)

    do j = 1, nlat
       if(maxlat < in_lat(j)) then
          jend = j
          exit
       endif
    enddo
    jend = min(jend+1,nlat)
  end subroutine jlimits
end subroutine read_field_N_3D

! ==============================================================================
! simplified interface for the time_inerp_external: takes care of creating
! the horizntal interpolator
! ==============================================================================
subroutine init_external_ts(ts, filename, fieldname, interp, fill)
  type(external_ts_type), intent(inout) :: ts
  character(*), intent(in) :: filename, fieldname
  character(*), intent(in) :: interp ! interpolation method
  real,         intent(in), optional :: fill ! fill value for missing data

! NOTE: filling missing data is not really implemented yet. It is not clear how to get
! the input data to determine the input valid data mask. Besides, missing input data mask
! may be different at different times -- not sure if time_interp_external can handle that
! at all.
! TODO: really implement missing data masking and filling

  integer :: axis_sizes(4)
  type(axistype) :: axis_centers(4), axis_bounds(4)
  real, allocatable :: lon_in(:), lat_in(:)

  ! initialize external field
  ts%filename = filename
  ts%fieldname = fieldname
  ts%id = init_external_field(filename,fieldname, domain=lnd%sg_domain, &
       axis_centers=axis_centers, axis_sizes=axis_sizes, &
       use_comp_domain=.TRUE., override=.TRUE.)
  !  get lon and lat of the input (source) grid, assuming that axis%data contains
  !  lat and lon of the input grid (in degrees)
  call get_axis_bounds(axis_centers(1),axis_bounds(1),axis_centers)
  call get_axis_bounds(axis_centers(2),axis_bounds(2),axis_centers)
  allocate(lon_in(axis_sizes(1)+1))
  allocate(lat_in(axis_sizes(2)+1))
  call mpp_get_axis_data(axis_bounds(1),lon_in)
  call mpp_get_axis_data(axis_bounds(2),lat_in)

  select case (trim(interp))
  case ('bilinear')
     call horiz_interp_new(ts%interp, lon_in*PI/180, lat_in*PI/180, lnd%sg_lon, lnd%sg_lat, &
          interp_method='bilinear')
  case ('conservative')
     call horiz_interp_new(ts%interp, lon_in*PI/180, lat_in*PI/180, lnd%sg_lonb, lnd%sg_latb, &
          interp_method='conservative')
  case default
     call error_mesg('init_external_ts','Unknown interpolation method "'//trim(interp)//'". use "bilinear" or "conservative"', FATAL)
  end select
  deallocate(lon_in,lat_in)
  ts%fill = DEFAULT_FILL_REAL
  if (present(fill)) ts%fill = fill
end subroutine init_external_ts


! ==============================================================================
subroutine del_external_ts(ts)
  type(external_ts_type), intent(inout) :: ts
  call horiz_interp_del(ts%interp)
end subroutine del_external_ts


! ==============================================================================
subroutine read_external_ts(ts,time,data_ug)
  type(external_ts_type), intent(in)  :: ts
  type(time_type),        intent(in)  :: time
  real,                   intent(out) :: data_ug(:)

  real :: data_sg(lnd%is:lnd%ie,lnd%js:lnd%je)

  call time_interp_external(ts%id, time, data_sg, horz_interp=ts%interp)
  call mpp_pass_sg_to_ug(lnd%ug_domain, data_sg, data_ug)
end subroutine read_external_ts

end module
