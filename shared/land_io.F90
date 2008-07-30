module land_io_mod

use constants_mod,     only : PI
use fms_mod,           only : file_exist, error_mesg, FATAL, stdlog, mpp_pe, &
     mpp_root_pe, write_version_number
use horiz_interp_mod,  only : horiz_interp_type, &
     horiz_interp_new, horiz_interp_del, &
     horiz_interp

use nf_utils_mod,      only : nfu_validtype, nfu_get_dim, nfu_get_dim_bounds, &
     nfu_get_valid_range, nfu_is_valid

implicit none
private

! ==== public interface ======================================================
public :: nearest
public :: init_cover_field
public :: read_field

public :: print_netcdf_error
! ==== end of public interface ===============================================

interface read_field
   module procedure read_field_0
   module procedure read_field_1
end interface

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'land_io_mod', &
     version     = '$Id: land_io.F90,v 16.0 2008/07/30 22:13:09 fms Exp $', &
     tagname     = '$Name: perth $'


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine nearest(mask, lon, lat, plon, plat, iout, jout, aspect)
! finds nearest point that is not masked out in input data
  logical, intent(in) :: mask(:,:)  ! mask of valid input points (.true. if valid point)
  real,    intent(in) :: lon(:)     ! longitudes of input grd central points, radian
  real,    intent(in) :: lat(:)     ! latitudes of input grd central points, radian
  real,    intent(in) :: plon, plat ! coordinates of destination point, radian
  integer, intent(out):: iout, jout ! indices of nearest valid (unmasked) point
  real,    intent(in), optional :: &
       aspect ! aspect ratio of metrics; >1 means longitudinal direction is
              ! preferred

  ! ---- local vars ----------------------------------------------------------
  integer :: i,j
  real    :: r,r1
  real    :: zsize

  zsize = 1.0
  if(present(aspect)) zsize = aspect 

  r = 2*PI*max(zsize,1.0)  ! some value larger than any possible distance

  do j = 1, size(mask,2)
  do i = 1, size(mask,1)
     if (.not.mask(i,j)) cycle
     r1 = distance(plon,plat,lon(i),lat(j))
     if ( r1 < r ) then
        iout = i
        jout = j
        r = r1
     endif
  enddo
  enddo

  ! ---- NOTES ---------------------------------------------------------------
  ! implemented in very naive and inefficient approach
  ! ---- BUGS ----------------------------------------------------------------
  ! aspect ratio is not used in this formulation: one can fix it by multiplying
  ! longitude difference by aspect (and confining the result to [-PI, +PI]). 
  ! however, in this case the actual aspect ratio (as measured by raitio of 
  ! distances in two directions) would depend on latitude.

contains 
  real function distance(lon1, lat1, lon2, lat2)
    ! calculates distance between points on unit square
    real, intent(in) :: lon1,lat1,lon2,lat2
    
    real :: x1,y1,z1, x2,y2,z2
    real :: dlon
    dlon = (lon2-lon1)
    
    z1 = sin(lat1)*zsize ;  z2 = sin(lat2)*zsize
    y1 = 0.0             ;  y2 = cos(lat2)*sin(dlon)
    x1 = cos(lat1)       ;  x2 = cos(lat2)*cos(dlon)
    
    ! distance = acos(x1*x2 + z1*z2)
    distance = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
  end function distance

end subroutine nearest


! ============================================================================
! This procedure creates and initializes a field of fractional coverage.
! It should be called for global grid; otherwise the part that fills
! missing data points may fail to find any good data, or do it in nproc-
! dependent way.
subroutine init_cover_field( &
     cover_to_use, filename, cover_field_name, frac_field_name, &
     glonb, glatb, uniform_cover, input_cover_types, frac)
  character(len=*), intent(in) :: cover_to_use
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: cover_field_name, frac_field_name
  real            , intent(in) :: glonb(:,:), glatb(:,:) ! boundaries of the grid cells
  integer         , intent(in) :: uniform_cover
  integer         , intent(in) :: input_cover_types(:)
  real            , intent(out):: frac(:,:,:) ! output-global map of soil fractional coverage

  ! ---- local vars ---------------------------------------------------------
  integer :: i,j,k     ! iterators
  integer :: cover_id
  real    :: maxfrac, total

  frac = 0
  
  if (cover_to_use == 'multi-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,glonb,glatb,input_cover_types,frac)
  else if (cover_to_use=='single-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,glonb,glatb,input_cover_types,frac)
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
     call read_cover_field(filename,cover_field_name,frac_field_name,glonb,glatb,input_cover_types,frac)
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
  integer :: ncid, varid

  if (.not.file_exist(file)) &
       call error_mesg(module_name,'input file "'//trim(file)//'" does not exist',FATAL)

  __NF_ASRT__( nf_open(file, NF_NOWRITE, ncid) )
  if(nf_inq_varid(ncid,cover_field_name,varid)==NF_NOERR) then
     call do_read_cover_field(ncid,varid,lonb,latb,input_cover_types,frac)
  else if ( nf_inq_varid(ncid,frac_field_name,varid)==NF_NOERR) then
     call do_read_fraction_field(ncid,varid,lonb,latb,input_cover_types,frac)
  else
     call error_mesg(module_name,&
          'neither "'//trim(cover_field_name)//'" nor "'//&
          frac_field_name//'" is present in input file "'//trim(file)//'"' ,&
          FATAL)
  endif
  __NF_ASRT__( nf_close(ncid) )

end subroutine read_cover_field

! ============================================================================
subroutine do_read_cover_field(ncid,varid,lonb,latb,input_cover_types,frac)
  integer, intent(in)  :: ncid, varid
  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:,:)

  ! ---- local vars
  integer :: nlon, nlat ! size of input map
  integer :: k
  integer, allocatable :: in_cover(:,:)
  real, allocatable    :: in_lonb(:), in_latb(:), x(:,:)
  type(horiz_interp_type) :: interp
  integer :: vardims(NF_MAX_VAR_DIMS)
  type(nfu_validtype) :: v

  ! find out dimensions, etc
  __NF_ASRT__( nf_inq_vardimid(ncid,varid,vardims) )
  ! get size of the longitude and latitude axes
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(1), nlon) )
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(2), nlat) )
  allocate (                         &
       in_lonb (nlon+1)             ,&
       in_latb (nlat+1)             ,&
       x       (0:nlon-1, 0:nlat-1) ,&
       in_cover(0:nlon-1, 0:nlat-1)  )
  __NF_ASRT__( nfu_get_dim_bounds(ncid, vardims(1), in_lonb) )
  __NF_ASRT__( nfu_get_dim_bounds(ncid, vardims(2), in_latb) )
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0

  ! read input data
  __NF_ASRT__( nf_get_var_int(ncid,varid,in_cover) )
  __NF_ASRT__( nfu_get_valid_range(ncid,varid,v) )

  call horiz_interp_new(interp, in_lonb,in_latb, lonb,latb, &
       interp_method='conservative')
  frac=0
  do k = 1,size(input_cover_types(:))
     x=0
     where(nfu_is_valid(in_cover,v).and.in_cover==input_cover_types(k)) x = 1

     call horiz_interp(interp,x,frac(:,:,k))
  enddo

  call horiz_interp_del(interp)

  ! clean up memory
  deallocate(in_lonb, in_latb, in_cover, x)

end subroutine do_read_cover_field


! ============================================================================
subroutine do_read_fraction_field(ncid,varid,lonb,latb,input_cover_types,frac)
  integer, intent(in)  :: ncid, varid
  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:,:)

  ! ---- local vars
  integer :: nlon, nlat, ntypes, k, cover
  real, allocatable :: in_frac(:,:,:)
  real, allocatable :: in_lonb(:), in_latb(:)
  real, allocatable :: in_mask(:,:)
  type(horiz_interp_type) :: interp
  type(nfu_validtype) :: v
  integer :: vardims(NF_MAX_VAR_DIMS)

  ! find out dimensions, etc
  __NF_ASRT__( nf_inq_vardimid(ncid,varid,vardims) )
  ! get size of the longitude and latitude axes
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(1), nlon) )
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(2), nlat) )
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(3), ntypes))
  allocate (                         &
       in_lonb (nlon+1)             ,&
       in_latb (nlat+1)             ,&
       in_mask (0:nlon-1, 0:nlat-1) ,&
       in_frac (0:nlon-1, 0:nlat-1, ntypes) )
  __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(1), in_lonb))
  __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(2), in_latb))
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0

  ! read input data
  __NF_ASRT__( nf_get_var_double(ncid,varid,in_frac) )
  __NF_ASRT__( nfu_get_valid_range(ncid,varid,v) )

  call horiz_interp_new(interp, in_lonb,in_latb, lonb,latb,&
       interp_method='conservative')
  frac = 0
  do k = 1,size(input_cover_types)
     cover = input_cover_types(k)
     in_mask = 0.0
     where(nfu_is_valid(in_frac(:,:,cover),v)) in_mask = 1.0
     call horiz_interp(interp,in_frac(:,:,cover),frac(:,:,k),mask_in=in_mask)
  enddo
  call horiz_interp_del(interp)

  ! clean up memory
  deallocate(in_lonb, in_latb, in_frac, in_mask)

end subroutine do_read_fraction_field


! ============================================================================
subroutine read_field_1(filename, varname, lon, lat, data, interp)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  real, intent(in)  :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:)
  character(len=*), intent(in), optional :: interp

  ! ---- local vars ----------------------------------------------------------
  integer :: ncid
  integer :: iret

  iret = nf_open(filename,NF_NOWRITE,ncid)
  if(iret/=NF_NOERR) then
     call error_mesg('read_field','Can''t open netcdf file "'//trim(filename)//'"',FATAL)
  endif
  call read_field_0(ncid, varname, lon, lat, data, interp)
  __NF_ASRT__( nf_close(ncid) )

end subroutine read_field_1

! ============================================================================
subroutine read_field_0(ncid, varname, lon, lat, data, interp)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real, intent(in) :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:)
  character(len=*), intent(in), optional  :: interp

  ! ---- local vars ----------------------------------------------------------
  integer :: varid
  integer :: nlon, nlat ! size of input grid
  integer :: vardims(NF_MAX_VAR_DIMS)
  real,    allocatable :: in_lonb(:), in_latb(:), in_lon(:), in_lat(:), x(:,:), rmask(:,:)
  logical, allocatable :: mask(:,:)
  character(len=20) :: interpolation 
  integer :: i,j,imap,jmap !
  real    :: missing

  interpolation = "bilinear"
  if(present(interp)) interpolation = interp
  
  ! get the dimensions of our variable
  __NF_ASRT__( nf_inq_varid(ncid,varname,varid) )
  __NF_ASRT__( nf_inq_vardimid(ncid,varid,vardims) )
  ! get size of the longitude and latitude axes
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(1), nlon) )
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(2), nlat) )

  allocate (                 &
       in_lon  (nlon),   in_lat  (nlat),   &
       in_lonb (nlon+1), in_latb (nlat+1), &
       x       (nlon, nlat) ,&
       mask    (nlon, nlat) , rmask(nlon, nlat) )

  ! read boundaries of the grid cells in longitudinal direction
  __NF_ASRT__(nfu_get_dim(ncid, vardims(1), in_lon))
  __NF_ASRT__(nfu_get_dim(ncid, vardims(2), in_lat))
  in_lon = in_lon*PI/180.0; in_lat = in_lat*PI/180.0
  __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(1), in_lonb))
  __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(2), in_latb))
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0
  ! read input data
  __NF_ASRT__( nf_inq_varid(ncid,varname,varid) )
  __NF_ASRT__( nf_get_var_double(ncid,varid,x) ) ! assuming real is real*8
  mask = .true.
  rmask = 1.0
  if(nf_get_att_double(ncid,varid,"missing_value",missing)==NF_NOERR) then
     mask = x/=missing
     where(.not.mask) rmask = 0.0
  endif

  select case(trim(interpolation))
  case ("bilinear")
     call horiz_interp(x, in_lonb, in_latb, lon,lat, data, mask_in=rmask, interp_method='bilinear')
  case ("nearest")
     do j = 1,size(data,2)
     do i = 1,size(data,1)
        call nearest (mask, in_lon, in_lat, lon(i,j), lat(i,j), imap, jmap)
        data(i,j) = x(imap,jmap)
     enddo
     enddo
  case default
     call error_mesg(module_name, interpolation//" is not a valid interpolation method",FATAL)
  end select

  deallocate(in_lonb, in_latb, in_lon, in_lat, x, mask, rmask)

end subroutine read_field_0

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
