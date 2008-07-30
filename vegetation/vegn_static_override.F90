module static_vegn_mod

use constants_mod,      only : pi
use time_manager_mod,   only : time_type, set_date
use get_cal_time_mod,   only : get_cal_time
use fms_mod,            only : write_version_number, error_mesg, FATAL, NOTE, &
     mpp_pe, file_exist, open_namelist_file, close_file, check_nml_error, stdlog, &
     mpp_root_pe
use time_interp_mod,    only : time_interp

use nf_utils_mod,       only : nfu_inq_dim, nfu_get_dim, nfu_inq_compressed_var, &
     nfu_get_compressed_rec
use land_data_mod,      only : lnd
use land_io_mod,        only : print_netcdf_error, nearest
use land_tile_mod,      only : land_tile_type, land_tile_enum_type, first_elmt, &
     tail_elmt, next_elmt, current_tile, operator(/=)
use vegn_cohort_mod,    only : vegn_cohort_type


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
     version     = '$Id: vegn_static_override.F90,v 16.0 2008/07/30 22:30:22 fms Exp $', &
     tagname     = '$Name: perth $'

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.
integer :: ncid ! netcdf id of the input file
type(time_type),allocatable :: time_line(:) ! time line of input data
type(time_type)             :: ts,te        ! beginning and end of time interval
integer, allocatable :: map_i(:,:), map_j(:,:)! remapping arrays: for each of the
     ! land grid cells in current domain they hold indices of corresponding points 
     ! in the input grid.
! ---- namelist variables ---------------------------------------------------
logical :: use_static_veg = .FALSE.
character(len=512) :: input_file = & ! name of input file for static vegetation
     "INPUT/static_veg_data.nc"
character(len=10)  :: timeline   = 'normal' ! type of timeline ('normal' or 'loop')
integer, dimension(6) :: &
     start_loop = (/1,1,1,0,0,0/), & ! beginning of the time loop
     end_loop   = (/1,1,1,0,0,0/)    ! end of the time loop
logical :: fill_land_mask = .FALSE. ! if true, all the vegetation points on the
     ! map are filled with the information from static vegetation data, using
     ! nearest point remap; otherwise only the points that overlap with valid
     ! static vegetation data are overriden.
namelist/static_veg_nml/use_static_veg,input_file,timeline,start_loop,end_loop,&
     fill_land_mask

! ==== NetCDF declarations ==================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

contains

! ===========================================================================
subroutine static_vegn_init()

  ! ---- local vars
  integer :: unit, ierr, io, unlimdim, timelen, timeid
  integer :: i,j
  character(len=NF_MAX_NAME) :: dimname  ! name of the dimension variable : time, lon, and lat
  integer                    :: ndims    ! rank of input vars
  integer                    :: dimids (NF_MAX_VAR_DIMS) ! netcdf IDs of input var dimensions
  integer                    :: dimlens(NF_MAX_VAR_DIMS) ! sizes of respective dimensions
  real, allocatable          :: t(:)     ! temporary real timeline
  character(len=256)         :: units    ! units of time in the file
  character(len=256)         :: calendar ! calendar of the data
  real, allocatable          :: in_lon(:)! longitude coordinates in input file
  real, allocatable          :: in_lat(:)! latitude coordinates in input file
  logical, allocatable       :: mask(:,:)! mask of valid points in input data 
  integer, allocatable       :: data(:,:,:,:) ! temprary array used to calculate the mask of
                                         ! valid input data

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
     __NF_ASRT__(nf_inq_dimname ( ncid, unlimdim, dimname ))
     __NF_ASRT__(nf_inq_varid   ( ncid, dimname, timeid ))
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

     ! READ HORIZONTAL COORDINATES
     __NF_ASRT__(nfu_inq_compressed_var(ncid,'species',ndims=ndims,dimids=dimids,dimlens=dimlens))
     allocate(in_lon(dimlens(1)),in_lat(dimlens(2)))
     __NF_ASRT__(nfu_get_dim(ncid,dimids(1),in_lon)) ! get longitude
     __NF_ASRT__(nfu_get_dim(ncid,dimids(2),in_lat)) ! get latitude
     in_lon = in_lon*PI/180.0 ; in_lat = in_lat*PI/180.0

     ! COMPUTE INDEX REMAPPING ARRAY
     allocate(map_i(lnd%is:lnd%ie,lnd%js:lnd%je))
     allocate(map_j(lnd%is:lnd%ie,lnd%js:lnd%je))
     allocate(mask(size(in_lon),size(in_lat)))

     if(fill_land_mask) then
        ! READ THE FIRST RECORD AND CALCULTE THE MASK OF THE VALID INPUT DATA
        allocate(data(dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
        !               lon        lat        tile       cohort
        data(:,:,:,:) = -1
        __NF_ASRT__(nfu_get_compressed_rec(ncid,'species',1,data))
        mask(:,:) = (data(:,:,1,1)>=0)
        deallocate(data)
     else
        mask(:,:) = .TRUE.
     endif
   
     do j = lnd%js,lnd%je
     do i = lnd%is,lnd%ie
        call nearest(mask,in_lon,in_lat,lnd%glon(i,j),lnd%glat(i,j),map_i(i,j),map_j(i,j))
     enddo
     enddo

     deallocate (in_lon,in_lat,mask)
     deallocate(t)
  endif

  module_is_initialized = .true.

end subroutine static_vegn_init

! ===========================================================================
subroutine static_vegn_end()
  if(.not.use_static_veg) return;

  __NF_ASRT__(nf_close(ncid))

  deallocate(time_line)
  deallocate(map_i,map_j)
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
  call read_remap_cohort_data_i0d_fptr(ncid, 'species' , cohort_species_ptr , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'bl'      , cohort_bl_ptr      , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'blv'     , cohort_blv_ptr     , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'br'      , cohort_br_ptr      , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'bsw'     , cohort_bsw_ptr     , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'bwood'   , cohort_bwood_ptr   , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'bliving' , cohort_bliving_ptr , map_i, map_j, index1)
  call read_remap_cohort_data_i0d_fptr(ncid, 'status'  , cohort_status_ptr  , map_i, map_j, index1)

  ! derived variables will be updated in update_land_bc_fast
end subroutine static_vegn_override


! ============================================================================
subroutine read_remap_cohort_data_i0d_fptr(ncid,name,fptr,map_i,map_j,rec)
  integer           , intent(in) :: ncid ! netcdf id
  character(len=*)  , intent(in) :: name ! name of the variable to read
  integer           , intent(in) :: map_i(lnd%is:,lnd%js:) ! re-mapping index
  integer           , intent(in) :: map_j(lnd%is:,lnd%js:) ! re-mapping index
  integer, optional , intent(in) :: rec  ! record number (in case there are 
                                         ! several in the file) 
  ! subroutine returning the pointer to the data to be written
  interface
     subroutine fptr(cohort, ptr)
       use vegn_cohort_mod, only : vegn_cohort_type
       type(vegn_cohort_type), pointer :: cohort ! input
       integer               , pointer :: ptr    ! returned pointer to the data
     end subroutine fptr
  end interface

  ! ---- local constants
  character(*), parameter :: module_name = 'read_remap_cohort_data_i0d_fptr'

  ! ---- local vars
  integer :: i,j,k,n,ii,jj,ndims, iret
  integer :: rec_     ! record number
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type)   , pointer :: tile
  type(vegn_cohort_type) , pointer :: cohort
  integer, pointer :: ptr ! pointer to the individual cohort data
  integer, allocatable :: data(:,:,:,:) ! buffer for input data
  logical, allocatable :: mask(:,:,:,:) ! validity mask for input data
  logical :: has_records
  integer :: dimlens(NF_MAX_VAR_DIMS)
  integer :: recsize

  ! assign the internal record number
  if(present(rec)) then
     rec_ = rec
  else
     rec_ = 1
  endif

  ! get the size of dimensions
  iret=nfu_inq_compressed_var(ncid,name,ndims=ndims,dimlens=dimlens,recsize=recsize,has_records=has_records)
  __NF_ASRT__(iret)
  ! check the number of dimensions of the variable
  if (has_records) then
     if(ndims/=5) call error_mesg(module_name,&
          'time-dependent variable "'//trim(name)//'" has incorrect number of dimensions',&
          FATAL)
  else
     if(ndims/=4) call error_mesg(module_name,&
          'time-independent variable "'//trim(name)//'" has incorrect number of dimensions',&
          FATAL)
  endif

  ! allocate input buffers
  allocate(data(dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
  allocate(mask(dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
  !             lon        lat        tile       cohort

  mask = .FALSE.
  __NF_ASRT__(nfu_get_compressed_rec(ncid,name,rec_,data,mask))

  ! distribute data over cohorts. NOTE that this is slightly different from the restart
  ! reading procedure. On reading the restart, all the tiles are counted in sequence,
  ! while here only tne vegetation tiles.
  do j = lnd%js, lnd%je
  do i = lnd%is, lnd%ie
     ii = map_i(i,j); jj = map_j(i,j)
     if ((ii.le.0).or.(jj.le.0)) cycle ! skip un-mapped points
     if (.not.any(mask(ii,jj,:,:))) cycle ! skip points where there is no data 

     ce = first_elmt (lnd%tile_map(i,j))
     te = tail_elmt  (lnd%tile_map(i,j))
     k = 1
     do while(ce/=te.and.k<=dimlens(3))
        tile=>current_tile(ce); ce=next_elmt(ce);
        if (.not.associated(tile%vegn)) cycle
        do n = 1,min(size(tile%vegn%cohorts(:)),dimlens(4))
           cohort=>tile%vegn%cohorts(n)
           call fptr(cohort,ptr)
           if(associated(ptr).and.mask(ii,jj,k,n)) ptr = data(ii,jj,k,n)
        enddo
        k = k+1
     enddo
  enddo
  enddo
  
  ! free allocated memory
  deallocate(data,mask)

end subroutine read_remap_cohort_data_i0d_fptr

! ============================================================================
subroutine read_remap_cohort_data_r0d_fptr(ncid,name,fptr,map_i,map_j,rec)
  integer           , intent(in) :: ncid ! netcdf id
  character(len=*)  , intent(in) :: name ! name of the variable to read
  integer           , intent(in) :: map_i(lnd%is:,lnd%js:) ! re-mapping index
  integer           , intent(in) :: map_j(lnd%is:,lnd%js:) ! re-mapping index
  integer, optional , intent(in) :: rec  ! record number (in case there are 
                                         ! several in the file) 
  ! subroutine returning the pointer to the data to be written
  interface
     subroutine fptr(cohort, ptr)
       use vegn_cohort_mod, only : vegn_cohort_type
       type(vegn_cohort_type), pointer :: cohort ! input
       real                  , pointer :: ptr    ! returned pointer to the data
     end subroutine fptr
  end interface

  ! ---- local constants
  character(*), parameter :: module_name = 'read_remap_cohort_data_i0d_fptr'

  ! ---- local vars
  integer :: i,j,k,n,ii,jj,ndims, iret
  integer :: rec_     ! record number
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type)   , pointer :: tile
  type(vegn_cohort_type) , pointer :: cohort
  real, pointer :: ptr ! pointer to the individual cohort data
  real,    allocatable :: data(:,:,:,:) ! buffer for input data
  logical, allocatable :: mask(:,:,:,:) ! validity mask for input data
  logical :: has_records
  integer :: dimlens(NF_MAX_VAR_DIMS)
  integer :: recsize

  ! assign the internal record number
  if(present(rec)) then
     rec_ = rec
  else
     rec_ = 1
  endif

  ! get the size of dimensions
  iret=nfu_inq_compressed_var(ncid,name,ndims=ndims,dimlens=dimlens,recsize=recsize,has_records=has_records)
  __NF_ASRT__(iret)
  ! check the number of dimensions of the variable
  if (has_records) then
     if(ndims/=5) call error_mesg(module_name,&
          'time-dependent variable "'//trim(name)//'" has incorrect number of dimensions',&
          FATAL)
  else
     if(ndims/=4) call error_mesg(module_name,&
          'time-independent variable "'//trim(name)//'" has incorrect number of dimensions',&
          FATAL)
  endif

  ! allocate input buffers
  allocate(data(dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
  allocate(mask(dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
  !             lon        lat        tile       cohort

  mask = .FALSE.
  __NF_ASRT__(nfu_get_compressed_rec(ncid,name,rec_,data,mask))

  ! distribute data over cohorts. NOTE that this is slightly different from the restart
  ! reading procedure. On reading the restart, all the tiles are enumerated in sequence,
  ! while here only the vegetation tiles -- non-vegetatet tiles are simply skipped.
  do j = lnd%js, lnd%je
  do i = lnd%is, lnd%ie
     ii = map_i(i,j); jj = map_j(i,j)
     if ((ii.le.0).or.(jj.le.0)) cycle ! skip un-mapped points
     if (.not.any(mask(ii,jj,:,:))) cycle ! skip points where there is no data 

     ce = first_elmt (lnd%tile_map(i,j))
     te = tail_elmt  (lnd%tile_map(i,j))
     k = 1
     do while(ce/=te.and.k<=dimlens(3))
        tile=>current_tile(ce); ce=next_elmt(ce);
        if (.not.associated(tile%vegn)) cycle
        do n = 1,min(size(tile%vegn%cohorts(:)),dimlens(4))
           cohort=>tile%vegn%cohorts(n)
           call fptr(cohort,ptr)
           if(associated(ptr).and.mask(ii,jj,k,n)) ptr = data(ii,jj,k,n)
        enddo
        k = k+1
     enddo
  enddo
  enddo
  
  ! free allocated memory
  deallocate(data,mask)

end subroutine read_remap_cohort_data_r0d_fptr

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
