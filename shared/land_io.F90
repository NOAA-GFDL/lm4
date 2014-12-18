module land_io_mod

use mpp_io_mod, only : fieldtype, mpp_get_info, mpp_get_fields
use mpp_io_mod, only : mpp_get_axes, mpp_get_axis_data, mpp_read, validtype
use mpp_io_mod, only : mpp_get_atts, MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE
use mpp_io_mod, only : axistype, mpp_open, mpp_close, mpp_is_valid, mpp_get_file_name
use mpp_io_mod, only : mpp_get_field_index
use mpp_io_mod, only : mpp_get_file_name

use axis_utils_mod, only : get_axis_bounds

use constants_mod,     only : PI
use fms_mod,           only : file_exist, error_mesg, FATAL, stdlog, mpp_pe, &
     mpp_root_pe, write_version_number, string, check_nml_error, close_file

use mpp_mod, only: mpp_sync
#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif


use horiz_interp_mod,  only : horiz_interp_type, &
     horiz_interp_new, horiz_interp_del, &
     horiz_interp

use land_numerics_mod, only : nearest, bisect
use nf_utils_mod,      only : nfu_validtype, nfu_get_dim, nfu_get_dim_bounds, &
     nfu_get_valid_range, nfu_is_valid, nfu_inq_var, nfu_get_var

implicit none
private

! ==== public interface ======================================================
public :: init_cover_field
public :: read_field
public :: read_land_io_namelist

public :: print_netcdf_error

public :: input_buf_size
! ==== end of public interface ===============================================

interface read_field
   module procedure read_field_N_2D, read_field_N_3D
   module procedure read_field_I_2D, read_field_I_3D
   module procedure read_field_N_2D_int, read_field_N_3D_int
   module procedure read_field_I_2D_int, read_field_I_3D_int
end interface

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'land_io_mod', &
     version     = '$Id: land_io.F90,v 21.0 2014/12/15 21:50:52 fms Exp $', &
     tagname     = '$Name: ulm $'

logical :: module_is_initialized = .false.
character(len=64)  :: interp_method = "conservative"
integer :: input_buf_size = 65536 ! input buffer size for tile and cohort reading
namelist /land_io_nml/ interp_method, input_buf_size

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine read_land_io_namelist()
  integer :: io, ierr, unit


  module_is_initialized = .TRUE.

  ! [1] print out version number
  call write_version_number (version, tagname)

#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=land_io_nml, iostat=io)
     ierr = check_nml_error(io, 'land_io_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=land_io_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_io_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif   
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_io_nml)
     call close_file (unit)
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
  real            , intent(out):: frac(:,:,:) ! output-global map of soil fractional coverage

  ! ---- local vars ---------------------------------------------------------
  integer :: i,j,k     ! iterators
  integer :: cover_id
  real    :: maxfrac, total

  if( .not. module_is_initialized ) &
       call error_mesg(module_name,'land_io_init is not called', FATAL)
 
  frac = 0
  
  if (cover_to_use == 'multi-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
  else if (cover_to_use=='single-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
     do j = 1,size(frac,2)
     do i = 1,size(frac,1)
        total = sum(frac(i,j,:))
        if (total <= 0) cycle ! operate on valid input data points only
        maxfrac=0 ; cover_id=1
        do k = 1,size(frac,3)
           if(frac(i,j,k).gt.maxfrac) then
              maxfrac=frac(i,j,k)
              cover_id=k
           endif
        enddo
        ! set all fractions except dominant fraction to zero
        frac(i,j,:) = 0.0
        frac(i,j,cover_id) = total
     enddo
     enddo
  else if (cover_to_use == 'uniform') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
     do j = 1,size(frac,2)
     do i = 1,size(frac,1)
        total = sum(frac(i,j,:))
        if (total <= 0) cycle ! operate on valid input data points only
        ! set all fractions except dominant fraction to zero
        frac(i,j,:) = 0.0
        frac(i,j,uniform_cover) = total
     enddo
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
  real              , intent(out) :: frac(:,:,:)     ! resulting fractions
  integer, optional , intent(in)  :: input_cover_types(:)

  ! --- local vars
! integer :: ncid, varid
  integer :: input_unit , ndim , nvar , natt , nrec , iret
  type(fieldtype), allocatable, dimension(:) :: fields
  type(fieldtype) :: field

  if (.not.file_exist(file)) call error_mesg(module_name,'input file "'//trim(file)//'" does not exist',FATAL)

! If field named 'cover' does not exist in file then read field named 'frac'
! 'cover' does not exist in either ground_type.nc or cover_type.nc
! The extent of the third dimension is 10 in ground_type.nc and 11 in cover_type.nc
  call mpp_open(input_unit, trim(file), action=MPP_RDONLY, form=MPP_NETCDF, &
       threading=MPP_MULTI, fileset=MPP_SINGLE, iostat=iret)
  call mpp_get_info(input_unit,ndim,nvar,natt,nrec)
  allocate(fields(nvar))
  call mpp_get_fields(input_unit,fields)

! if(nf_inq_varid(ncid,cover_field_name,varid)==NF_NOERR) then
!    call do_read_cover_field(ncid,varid,lonb,latb,input_cover_types,frac)
! else if ( nf_inq_varid(ncid,frac_field_name,varid)==NF_NOERR) then
!     call do_read_fraction_field(ncid,varid,lonb,latb,input_cover_types,frac)

  if(get_field(fields,cover_field_name,field)==0) then
     call do_read_cover_field(input_unit,field,lonb,latb,input_cover_types,frac)
  else if ( get_field(fields,frac_field_name,field)==0) then
     call do_read_fraction_field(input_unit,field,lonb,latb,input_cover_types,frac)
  else
     call error_mesg(module_name,&
          'neither "'//trim(cover_field_name)//'" nor "'//&
          frac_field_name//'" is present in input file "'//trim(file)//'"' ,&
          FATAL)
  endif
  call mpp_close(input_unit)

end subroutine read_cover_field

! ============================================================================
function get_field(fields,field_name,field)
  type(fieldtype), intent(in) :: fields(:)
  character(len=*), intent(in) :: field_name
  type(fieldtype), intent(out) :: field
  integer :: get_field, n
  character(len=256) :: name_out

  n =  mpp_get_field_index(fields,trim(field_name))
  if ( n > 0 ) then
    get_field = 0
    field = fields(n)
  else
    get_field = 1
  endif

end function get_field
! ============================================================================
subroutine do_read_cover_field(input_unit, field, lonb, latb, input_cover_types, frac)
  integer, intent(in)  :: input_unit
  type(fieldtype), intent(in) :: field
  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:,:)

  ! ---- local vars
  integer :: nlon, nlat ! size of input map
  integer :: k
  integer, allocatable :: in_cover(:,:)
  real, allocatable    :: in_lonb(:), in_latb(:), x(:,:), r_in_cover(:,:)
  type(horiz_interp_type) :: interp
  integer :: vardims(1024)
  type(validtype) :: v
  integer :: in_j_start, in_j_count ! limits of the latitude belt we read
  integer :: num_dims, dimlens(1024)
  type(fieldtype), allocatable :: fields(:)
  type(axistype),  allocatable :: axes(:)
  type(axistype) :: axes_bnd
  integer :: ndim,nvar,natt,nrec
  integer :: start(4), count(4)

  ! find out dimensions, etc
  call mpp_get_info(input_unit,ndim,nvar,natt,nrec)
  allocate(axes(ndim), fields(nvar))
  call mpp_get_axes(input_unit, axes)
  call mpp_get_fields(input_unit, fields)
  ! get size of the longitude and latitude axes
  call mpp_get_atts(field, ndim=num_dims, siz=dimlens) ! field is of type fieldtype and replaces varid.
  nlon = dimlens(1); nlat = dimlens(2)
  allocate ( in_lonb(nlon+1), in_latb(nlat+1) )
  call get_axis_bounds(axes(1), axes_bnd, axes)
  call mpp_get_axis_data(axes_bnd, in_lonb(:))
  call get_axis_bounds(axes(2), axes_bnd, axes)
  call mpp_get_axis_data(axes_bnd, in_latb(:))
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0

  ! to minimize the i/o and work done by horiz_interp, find the boundaries
  ! of latitude belt in input data that covers the entire latb array
  in_j_start=bisect(in_latb, minval(latb))
  in_j_count=bisect(in_latb, maxval(latb))-in_j_start+1

  ! check for unreasonable values
  if (in_j_start<1) &
     call error_mesg('do_read_cover_field','input latitude start index ('&
                     //string(in_j_start)//') is out of bounds', FATAL)
  if (in_j_start+in_j_count-1>nlat) &
     call error_mesg('do_read_cover_field','input latitude count ('&
                     //string(in_j_count)//') is too large (start index='&
                     //string(in_j_start)//')', FATAL)

  ! allocate input data buffers
  allocate ( x(nlon,in_j_count), in_cover(nlon,in_j_count), r_in_cover(nlon,in_j_count) )

  ! read input data
! iret = nf_get_vara_int(ncid,varid, (/1,in_j_start/), (/nlon,in_j_count/), in_cover)
  start = (/1,in_j_start,1,1/)
  count(1:2) = shape(in_cover)
  count(3:4) = 1
  call mpp_read( input_unit, field, r_in_cover, start, count )
  in_cover = r_in_cover
  call mpp_get_atts(field,valid=v)

  call horiz_interp_new(interp, in_lonb,in_latb(in_j_start:in_j_start+in_j_count), &
       lonb,latb, interp_method=trim(interp_method))
  frac=0
  do k = 1,size(input_cover_types(:))
     x=0
     where(mpp_is_valid(r_in_cover,v).and.in_cover==input_cover_types(k)) x = 1
     call horiz_interp(interp,x,frac(:,:,k))
  enddo

  call horiz_interp_del(interp)

  ! clean up memory
  deallocate(in_lonb, in_latb, in_cover, x)

end subroutine do_read_cover_field


! ============================================================================
 subroutine do_read_fraction_field(input_unit,field,lonb,latb,input_cover_types,frac)
  integer, intent(in)  :: input_unit
  type(fieldtype), intent(in) :: field
  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:,:)

  ! ---- local vars
  integer :: nlon, nlat, ntypes, k, cover
  real, allocatable :: in_frac(:,:,:)
  real, allocatable :: in_lonb(:), in_latb(:)
  real, allocatable :: in_mask(:,:)
  type(horiz_interp_type) :: interp
  type(validtype) :: v
  integer :: vardims(1024)
  integer :: in_j_end, in_j_start, in_j_count ! limits of the latitude belt we read
  integer :: num_dims, dimlens(1024)
  type(fieldtype), allocatable :: fields(:)
  type(axistype),  allocatable :: axes(:)
  type(axistype) :: axes_bnd
  integer :: ndim,nvar,natt,nrec
  integer :: start(4), count(4)
  character(len=256) :: bnd_name, err_msg

  ! find out dimensions, etc
  call mpp_get_info(input_unit,ndim,nvar,natt,nrec)
  allocate(axes(ndim), fields(nvar))
  call mpp_get_axes(input_unit, axes)
  call mpp_get_fields(input_unit, fields)
  ! get size of the longitude and latitude axes
  call mpp_get_atts(field, ndim=num_dims, siz=dimlens)
  nlon = dimlens(1); nlat = dimlens(2) ; ntypes = dimlens(3)
  allocate ( in_lonb(nlon+1), in_latb(nlat+1) )

  call get_axis_bounds(axes(1), axes_bnd, axes) ! This routine requires some explanation. See below.
! get_axis_bounds works like this:
! It looks at the axis appearing as the first argument to see if the name of a bounds axis exists.
! If it does then: It looks for the bounds axis amoung the array of axes supplied as the third argument and returns it as the second argument.
! If it does not then: It creates a bounds axis using data from the axis appearing as the first argument and returns this as the second argument.

  call mpp_get_axis_data(axes_bnd, in_lonb)
  err_msg = ''
  call get_axis_bounds(axes(2), axes_bnd, axes, bnd_name=bnd_name, err_msg=err_msg)
  if(len_trim(err_msg) > 0 ) then
    err_msg = trim(err_msg)//' File: '//trim(mpp_get_file_name(input_unit))
    call error_mesg('do_read_fraction_field',err_msg, FATAL)
  endif
  call mpp_get_axis_data(axes_bnd, in_latb)
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0
  ! find the boundaries of latitude belt in input data that covers the 
  ! entire latb array
  in_j_start=bisect(in_latb, minval(latb))
  in_j_count=bisect(in_latb, maxval(latb))-in_j_start+1
  in_j_end = in_j_start + in_j_count - 1

  ! check for unreasonable values
  if (in_j_start<1) &
     call error_mesg('do_read_fraction_field','input latitude start index ('&
                     //string(in_j_start)//') is out of bounds', FATAL)
  if (in_j_start+in_j_count-1>nlat) &
     call error_mesg('do_read_fraction_field','input latitude count ('&
                     //string(in_j_count)//') is too large (start index='&
                     //string(in_j_start)//')', FATAL)

  allocate( in_mask(nlon,in_j_count), in_frac(nlon,in_j_count,ntypes) )

  ! read input data
! iret = nf_get_vara_double(ncid, varid, (/1,in_j_start,1/), (/nlon,in_j_count,ntypes/), in_frac)
  start = (/1,in_j_start,1,1/)
  count(1:3) = shape(in_frac)
  count(4) = 1
  call mpp_read( input_unit, field, in_frac, start, count ) ! interface called here is mpp_read_region_r3D
  call mpp_get_atts(field,valid=v)

  frac = 0
  do k = 1,size(input_cover_types)
     cover = input_cover_types(k)
     if (cover<1.or.cover>ntypes) then
        cycle ! skip all invalid indices in the array of input cover types
     endif

     where(mpp_is_valid(in_frac(:,:,cover),v))
        in_mask = 1.0
     elsewhere
        in_mask = 0
     end where
     call horiz_interp_new(interp, &
          in_lonb,in_latb(in_j_start:in_j_start+in_j_count), lonb,latb,&
          interp_method=trim(interp_method), mask_in=in_mask)
     call horiz_interp(interp,in_frac(:,:,cover),frac(:,:,k))
     call horiz_interp_del(interp)
  enddo

  ! clean up memory
  deallocate(in_lonb, in_latb, in_frac, in_mask)

end subroutine do_read_fraction_field

! ============================================================================
subroutine read_field_N_2D_int(filename, varname, lon, lat, data, interp, mask)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  real, intent(in)  :: lon(:,:),lat(:,:)
  integer, intent(out) :: data(:,:)
  character(len=*), intent(in), optional :: interp
  logical, intent(out), optional :: mask(:,:)

  ! ---- local vars ----------------------------------------------------------
  real    :: data3(size(data,1),size(data,2),1)
  logical :: mask3(size(data,1),size(data,2),1)

  call read_field_N_3D(filename, varname, lon, lat, data3, interp, mask3)
  data = nint(data3(:,:,1))
  if (present(mask)) &
     mask = mask3(:,:,1)

end subroutine read_field_N_2D_int

! ============================================================================
subroutine read_field_N_3D_int(filename, varname, lon, lat, data, interp, mask)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  real, intent(in)  :: lon(:,:),lat(:,:)
  integer, intent(out) :: data(:,:,:)
  character(len=*), intent(in), optional :: interp
  logical, intent(out), optional :: mask(:,:,:)

  ! ---- local vars ----------------------------------------------------------
  real    :: data3(size(data,1),size(data,2),size(data,3))

  call read_field_N_3D(filename, varname, lon, lat, data3, interp, mask)
  data = nint(data3(:,:,:))

end subroutine read_field_N_3D_int

! ============================================================================
subroutine read_field_I_2D_int(ncid, varname, lon, lat, data, interp, mask)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real, intent(in) :: lon(:,:),lat(:,:)
  integer, intent(out) :: data(:,:)
  character(len=*), intent(in), optional  :: interp
  logical, intent(out), optional :: mask(:,:)
  ! ---- local vars
  real    :: data3(size(data,1),size(data,2),1)
  logical :: mask3(size(data,1),size(data,2),1)

  call read_field_I_3D(ncid, varname, lon, lat, data3, interp, mask3)
  data = nint(data3(:,:,1))
  if (present(mask)) mask = mask3(:,:,1)

end subroutine read_field_I_2D_int

! ============================================================================
subroutine read_field_I_3D_int(ncid, varname, lon, lat, data, interp, mask)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real, intent(in) :: lon(:,:),lat(:,:)
  integer, intent(out) :: data(:,:,:)
  character(len=*), intent(in), optional  :: interp
  logical, intent(out), optional :: mask(:,:,:)
  ! ---- local vars
  real    :: data3(size(data,1),size(data,2),size(data,3))

  call read_field_I_3D(ncid, varname, lon, lat, data3, interp, mask)
  data = nint(data3(:,:,:))

end subroutine read_field_I_3D_int

! ============================================================================
subroutine read_field_N_2D(filename, varname, lon, lat, data, interp, mask)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  real, intent(in)  :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:)
  character(len=*), intent(in), optional :: interp
  logical, intent(out), optional :: mask(:,:)

  ! ---- local vars ----------------------------------------------------------
  real    :: data3(size(data,1),size(data,2),1)
  logical :: mask3(size(data,1),size(data,2),1)

  call read_field_N_3D(filename, varname, lon, lat, data3, interp, mask3)
  data = data3(:,:,1)
  if (present(mask)) mask = mask3(:,:,1)

end subroutine read_field_N_2D

! ============================================================================
subroutine read_field_N_3D(filename, varname, lon, lat, data, interp, mask)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  real, intent(in)  :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:,:)
  character(len=*), intent(in), optional :: interp
  logical, intent(out), optional :: mask(:,:,:)

  ! ---- local vars ----------------------------------------------------------
  integer :: ierr, input_unit

  ! Files read: biodata.nc, geohydrology.nc, soil_brdf.nc
  input_unit = -9999
  call mpp_open(input_unit, trim(filename), action=MPP_RDONLY, form=MPP_NETCDF, &
           threading=MPP_MULTI, fileset=MPP_SINGLE, iostat=ierr)
  call read_field_I_3D(input_unit, varname, lon, lat, data, interp, mask)
  call mpp_sync()
  call mpp_close(input_unit)

end subroutine read_field_N_3D

! ============================================================================
subroutine read_field_I_2D(ncid, varname, lon, lat, data, interp, mask)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real, intent(in) :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:)
  character(len=*), intent(in), optional  :: interp
  logical, intent(out), optional :: mask(:,:)

  ! ---- local vars
  real    :: data3(size(data,1),size(data,2),1)
  logical :: mask3(size(data,1),size(data,2),1)

  call read_field_I_3D(ncid, varname, lon, lat, data3, interp, mask3)
  data = data3(:,:,1)
  if (present(mask)) mask = mask3(:,:,1)

end subroutine read_field_I_2D

! ============================================================================
subroutine read_field_I_3D(input_unit, varname, lon, lat, data, interp, mask)
  integer, intent(in) :: input_unit
  character(len=*), intent(in) :: varname
  real, intent(in) :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:,:)
  character(len=*), intent(in), optional  :: interp
  logical, intent(out), optional :: mask(:,:,:)

  ! ---- local vars ----------------------------------------------------------
  integer :: nlon, nlat, nlev ! size of input grid
  integer :: varndims ! number of variable dimension
  integer :: vardims(1024) ! IDs of variable dimension
  integer :: dimlens(1024) ! sizes of respective dimensions
  real,    allocatable :: in_lonb(:), in_latb(:), in_lon(:), in_lat(:)
  real,    allocatable :: x(:,:,:) ! input buffer
  logical, allocatable :: imask(:,:,:) ! mask of valid input values
  real,    allocatable :: rmask(:,:,:) ! real mask for interpolator
  real    :: omask(size(data,1),size(data,2),size(data,3)) ! mask of valid output data
  character(len=20) :: interpolation 
  integer :: i,j,k,imap,jmap !
  type(validtype) :: v
  type(horiz_interp_type) :: hinterp
  integer :: ndim,nvar,natt,nrec
  type(axistype), allocatable :: all_axes(:), varaxes(:)
  type(axistype):: axis_bnd
  type(fieldtype), allocatable :: fields(:)
  type(fieldtype) :: fld
  character(len=256) :: bnd_name, axis_name, field_name, file_name
  character(len=512) :: err_msg
  logical :: found
  integer :: bnd_ndims, bnd_index
  real, allocatable  :: bnd_data(:,:)
  integer :: bnd_dimlens(4)

  interpolation = "bilinear"
  if(present(interp)) interpolation = interp
  
  call mpp_get_info(input_unit,ndim,nvar,natt,nrec)
  allocate(fields(nvar), all_axes(ndim))
  call mpp_get_axes(input_unit, all_axes)
  call mpp_get_fields(input_unit, fields)
  k = mpp_get_field_index(fields,trim(varname))
  if(k > 0) then
    fld = fields(k)
  else
    file_name = mpp_get_file_name(input_unit)
    call error_mesg('read_field','variable "'//trim(varname)//'" not found in file "'//trim(file_name)//'"',FATAL)
  endif
  ! get the dimensions of our variable

  call mpp_get_atts(fld, ndim=varndims, siz=dimlens, valid=v)
  if(varndims<2.or.varndims>3) then
     call error_mesg('read_field','variable "'//trim(varname)//'" is '//string(varndims)//&
          'D, but only reading 2D or 3D variables is supported', FATAL)
  endif
  allocate(varaxes(varndims))
  call mpp_get_atts(fld, axes=varaxes)
  nlon = dimlens(1) ; nlat = dimlens(2)
  nlev = 1; 
  if (varndims==3) nlev=dimlens(3)
  if(nlev/=size(data,3)) then
     call error_mesg('read_field','3rd dimension length of the variable "'&
          //trim(varname)//'" ('//trim(string(nlev))//') is different from the expected size of data ('// &
          trim(string(size(data,3)))//')', FATAL)
  endif

  allocate (                 &
       in_lon  (nlon),   in_lat  (nlat),   &
       in_lonb (nlon+1), in_latb (nlat+1), &
       x       (nlon, nlat, nlev) ,&
       imask    (nlon, nlat, nlev) , rmask(nlon, nlat, nlev) )

  ! read boundaries of the grid cells in longitudinal direction
  call mpp_get_axis_data(varaxes(1), in_lon)
  call mpp_get_axis_data(varaxes(2), in_lat)
  in_lon = in_lon*PI/180.0; in_lat = in_lat*PI/180.0
  call get_axis_bounds(varaxes(1), axis_bnd, all_axes)
  call mpp_get_axis_data(axis_bnd, in_lonb)
  err_msg = ''
  call get_axis_bounds(varaxes(2), axis_bnd, all_axes, bnd_name=bnd_name, err_msg=err_msg)
  if(len_trim(err_msg) > 0 ) then
    ! The boundary data for varaxes(2) was not found amoung the axes passed in as "all_axes".
    ! Look for a field (as opposed to as axis) that has the name of the boundary data.
    bnd_index = mpp_get_field_index(fields,bnd_name)
    if(bnd_index > 0) then
      found = .TRUE.
    else
      found = .FALSE.
    endif

    ! There are multiple error checks in the code below.
    ! Each requires a long error message. The bulk of the error
    ! message is assigned here in case it is needed later.
    file_name = mpp_get_file_name(input_unit)
    call mpp_get_atts(varaxes(2), name=axis_name)
    err_msg = 'axis bound "'//trim(bnd_name)//'" of axis "'//trim(axis_name)
    err_msg = trim(err_msg)//'" of variable "'//trim(varname)//'" of file "'//trim(file_name)//'"'

    if(.not.found) then
      ! Neither an axis or a field was found that has the name of the boundary data.
      err_msg = trim(err_msg)//'  The axis bound name is not found in this file.'
      call error_mesg('read_field',trim(err_msg),FATAL)
    endif

    if(found) then
      ! A field was found that has the name of the boundary data.
      ! This field could have one or two dimensions.
      ! Lets handle both cases.
      call mpp_get_atts(fields(bnd_index), ndim=bnd_ndims, siz=bnd_dimlens)
      if(bnd_ndims < 1 .or. bnd_ndims > 2) then
        err_msg = trim(err_msg)//'  The axis bound data has more that 2 dimensions. Only 1 or 2 is allowed.'
        call error_mesg('read_field',trim(err_msg),FATAL)
      endif
      if(bnd_dimlens(2) < nlat .or. bnd_dimlens(2) > nlat+1) then
         err_msg = trim(err_msg)//'  The size of the axis bound data does not conform to the size of the axis data.'
         call error_mesg('read_field',trim(err_msg),FATAL)
      endif
      if(bnd_ndims == 1) then
        allocate(bnd_data(nlat+1,1))
      else if(bnd_ndims == 2) then
        if(bnd_dimlens(1) /= 2) then
          err_msg = trim(err_msg)//'  The first dimension of the bound data must be 2 or nlat'
          call error_mesg('read_field',trim(err_msg),FATAL)
        endif
        allocate(bnd_data(2,nlat))
      endif
      call mpp_read(input_unit, fields(bnd_index), bnd_data)
      if(bnd_ndims == 1) then
        in_latb = bnd_data(:,1)
      else
        in_latb(1:nlat) = bnd_data(1,:)
        in_latb(nlat+1) = bnd_data(2,nlat)
      endif
      deallocate(bnd_data)
    endif
  else
    call mpp_get_axis_data(axis_bnd, in_latb)
  endif
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0
  ! read input data
  call mpp_read(input_unit, fld, x)
  
  imask = mpp_is_valid(x,v)
  rmask = 1.0
  where(.not.imask) rmask = 0.0

  select case(trim(interpolation))
  case ("nearest")
     do k = 1,size(data,3)
     do j = 1,size(data,2)
     do i = 1,size(data,1)
        call nearest (imask(:,:,k), in_lon, in_lat, lon(i,j), lat(i,j), imap, jmap)
        data(i,j,k) = x(imap,jmap,k)
     enddo
     enddo
     enddo
  case default
     call horiz_interp_new(hinterp, in_lonb, in_latb, lon, lat, interp_method=interpolation)
     do k = 1,size(data,3)
        call horiz_interp(hinterp,x(:,:,k),data(:,:,k),mask_in=rmask(:,:,k), mask_out=omask(:,:,k))
     enddo
     if (present(mask)) mask(:,:,:) = (omask(:,:,:)/=0.0)
     call horiz_interp_del(hinterp)
  end select

  deallocate(in_lonb, in_latb, in_lon, in_lat, x, imask, rmask)
  deallocate(all_axes, varaxes, fields)

end subroutine read_field_I_3D

! ============================================================================
subroutine print_netcdf_error(ierr, file, line)
  ! prints out NetCDF library error message, including file name and line number
  integer,          intent(in) :: ierr ! error code
  character(len=*), intent(in) :: file ! name of the file
  integer,          intent(in) :: line ! number of line in the file

  ! ---- local vars
  character(len=1024) :: mesg

  if (ierr.ne.NF_NOERR) then
     write(mesg, "('File ',a,' Line ',i4.4,' :: ',a)") &
          trim(file),line,trim(NF_STRERROR(ierr))
     call error_mesg('NetCDF', mesg, FATAL)
  endif
end subroutine print_netcdf_error

end module
