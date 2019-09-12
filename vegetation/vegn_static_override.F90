module static_vegn_mod

use constants_mod,      only : pi
use mpp_mod,            only : mpp_max, mpp_sum
use time_manager_mod,   only : time_type, set_date, time_type_to_real, &
     get_calendar_type, valid_calendar_types, operator(-), get_date
use get_cal_time_mod,   only : get_cal_time

use fms_mod,            only : error_mesg, FATAL, NOTE, &
     mpp_pe, input_nml_file, check_nml_error, stdlog, lowercase, &
     mpp_root_pe, fms_error_handler
use time_interp_mod,    only : time_interp
use diag_manager_mod,   only : get_base_date

use land_data_mod,      only : log_version, lnd
use land_numerics_mod,  only : nearest
use land_tile_io_mod,   only : create_tile_out_file, gather_tile_index
use land_tile_mod,      only : land_tile_map, land_tile_type, land_tile_enum_type, first_elmt, &
     tail_elmt, next_elmt, current_tile, operator(/=), nitems
use vegn_cohort_mod,    only : vegn_cohort_type
use cohort_io_mod,      only : create_cohort_dimension_new, gather_cohort_data, &
     gather_cohort_index

use fms2_io_mod, only: FmsNetcdfUnstructuredDomainFile_t, register_axis, &
                       register_field, register_variable_attribute, unlimited, &
                       register_restart_field, write_restart, get_dimension_size, &
                       get_variable_size, read_data, FmsNetcdfFile_t, open_file, &
                       close_file, variable_exists, variable_att_exists, &
                       get_variable_attribute, get_variable_num_dimensions, get_unlimited_dimension_name

implicit none
private

! ==== public interface =====================================================
public :: read_static_vegn_namelist
public :: static_vegn_init
public :: static_vegn_end

public :: read_static_vegn
public :: write_static_vegn
! ==== end of public interface ==============================================

! ==== module constants =====================================================
character(len=*), parameter :: module_name = 'static_vegn_mod'
#include "../shared/version_variable.inc"

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.
integer :: ncid  ! netcdf id of the input file
integer, allocatable, dimension(:) :: species, status
real,    allocatable, dimension(:) :: bl, blv, br, bsw, bwood, bliving
type(time_type),allocatable :: time_line(:) ! time line of input data
type(time_type)             :: ts,te        ! beginning and end of time interval
integer, allocatable :: map_i(:), map_j(:)! remapping arrays: for each of the
     ! land grid cells in current domain they hold indices of corresponding points
     ! in the input grid.
type(time_type) :: base_time ! model base time for static vegetation output
integer :: ispecies, ibl, iblv, ibr, ibsw, ibwood, ibliving, istatus

type(FmsNetcdfUnstructuredDomainFile_t) :: static_veg_file ! handle of output file, for new IO

integer :: ncid2 ! netcdf id of the output file, for old IO
integer :: tile_dim_length ! length of tile dimension in output files. global max of number of tiles per gridcell
integer, allocatable :: cidx(:) ! cohort compression index, local for current PE

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
     ! static vegetation data are overridden.
logical :: write_static_veg = .FALSE. ! if true, the state of vegetation is saved
     ! periodically for future use as static vegetation input
character(16) :: static_veg_freq = 'daily' ! or 'monthly', or 'annual'
     ! specifies the frequency for writing the static vegetation data file

namelist/static_veg_nml/use_static_veg,input_file,timeline,start_loop,end_loop,&
     fill_land_mask, write_static_veg, static_veg_freq

logical :: input_is_multiface ! TRUE if the input files are face-specific
type(FmsNetcdfFile_t) :: fileobj
type(FmsNetcdfUnstructuredDomainFile_t) :: fileobj_domainug

contains

! ===========================================================================
subroutine read_static_vegn_namelist(static_veg_used)
  logical, intent(out) :: static_veg_used

  ! ---- local vars
  integer :: unit, ierr, io

  call log_version(version, module_name, &
  __FILE__)

  read (input_nml_file, nml=static_veg_nml, iostat=io)
  ierr = check_nml_error(io, 'static_veg_nml')

  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=static_veg_nml)
  endif

  if (    (trim(static_veg_freq)=='daily') &
      .or.(trim(static_veg_freq)=='monthly') &
      .or.(trim(static_veg_freq)=='annual') ) then
     ! static_veg_freq is OK -- do nothing
  else
     call error_mesg('static_vegn_init','option static_veg_freq="'&
          //trim(static_veg_freq)&
          //'" is invalid, use "daily", "monthly", or "annual"', FATAL)
  endif

  static_veg_used = use_static_veg
end subroutine read_static_vegn_namelist


! ===========================================================================
subroutine static_vegn_init( )

  ! ---- local vars
  integer :: i,j,k
  integer, dimension(5) :: dimlens ! sizes of respective dimensions (lon, lat, tile, cohort, time)
  character(len=256)         :: units    ! units of time in the file
  character(len=256)         :: calendar ! calendar of the data
  real, allocatable          :: in_lon(:)! longitude coordinates in input file
  real, allocatable          :: in_lat(:)! latitude coordinates in input file
  logical, allocatable       :: mask(:,:)! mask of valid points in input data
  integer :: m, n, l
  integer, dimension(:), allocatable :: siz

!  logical :: input_is_multiface ! TRUE if the input files are face-specific
  integer, allocatable :: cidx(:), idata(:)
  logical:: exists
  real, dimension(:), allocatable :: t
  integer :: ndims
  character(len=256) :: dimension_name !NAME OF UNLIMITED DIMENSION (i.e "time")

  if(module_is_initialized) return

  if(use_static_veg) then
    ! SET UP LOOP BOUNDARIES
    ts = set_date(start_loop(1),start_loop(2),start_loop(3), start_loop(4),start_loop(5),start_loop(6))
    te = set_date(end_loop(1)  ,end_loop(2)  ,end_loop(3)  , end_loop(4)  ,end_loop(5)  ,end_loop(6)  )

    ! OPEN INPUT FILE
    exists = open_file(fileobj, input_file, "read")
    if (exists) then
       call error_mesg('static_vegn_init','Reading global static vegetation file "'&
            //trim(input_file)//'"', NOTE)
       input_is_multiface = .FALSE.
    else
       if(lnd%nfaces==1) then
          ! for 1-face grid we cannot use multi-face input, even if it exists
          call error_mesg('static_vegn_init','input file "'//trim(input_file)&
                  //'" does not exist', FATAL)
       else
          ! if there is more then one face, try opening face-specific input with consideration of io_layout
          exists = open_file(fileobj_domainug, input_file, "read", lnd%ug_domain)
          if (.not. exists) then
             call error_mesg('static_vegn_init','"'//trim(input_file)// &
             '" and corresponding distributed file are not found', FATAL)

          endif
          call error_mesg('static_vegn_init','Reading face-specific vegetation file "'&
               //trim(input_file)//'"', NOTE)
          input_is_multiface = .TRUE.
       endif
    endif

    units = ' '
    calendar = 'JULIAN'

    ! READ TIME AXIS DATA
    ! GET UNITS OF THE TIME AND CALENDAR OF THE DATA
    ! CONVERT TIME TO THE FMS TIME_TYPE AND STORE IT IN THE TIMELINE FOR THE DATA SET
    ! READ HORIZONTAL COORDINATES
    if (input_is_multiface) then
      ! Get the name of the unlimited dimension (i.e "time")
      call get_unlimited_dimension_name(fileobj_domainug, dimension_name)
      call get_variable_attribute(fileobj_domainug, dimension_name, "units", units)
      call get_variable_attribute(fileobj_domainug, dimension_name, "calendar", calendar)
      call get_dimension_size(fileobj_domainug, dimension_name, dimlens(5))
      call get_dimension_size(fileobj_domainug, "lon", dimlens(1))
      call get_dimension_size(fileobj_domainug, "lat", dimlens(2))
    else
      ! Get the name of the unlimited dimension (i.e "time")
      call get_unlimited_dimension_name(fileobj_domainug, dimension_name)
      call get_variable_attribute(fileobj_domainug, dimension_name, "units", units)
      call get_variable_attribute(fileobj_domainug, dimension_name, "calendar", calendar)
      call get_dimension_size(fileobj_domainug, dimension_name, dimlens(5))
      call get_dimension_size(fileobj, "lon", dimlens(1))
      call get_dimension_size(fileobj, "lat", dimlens(2))
    endif
    allocate(t(dimlens(5)))
    allocate(time_line(dimlens(5)))
    allocate(in_lon(dimlens(1)))
    allocate(in_lat(dimlens(2)))
    if (input_is_multiface) then
      call read_data(fileobj_domainug, dimension_name, t)
      call read_data(fileobj_domainug, "lon", in_lon)
      call read_data(fileobj_domainug, "lat", in_lat)
    else
      call read_data(fileobj, dimension_name, t)
      call read_data(fileobj, "lon", in_lon)
      call read_data(fileobj, "lat", in_lat)
    endif
    do i = 1, dimlens(5)
       ! set the respective value in the timeline
       time_line(i) = get_cal_time(t(i), units, calendar)
    enddo
    in_lon = in_lon*PI/180.0
    in_lat = in_lat*PI/180.0

    ! COMPUTE INDEX REMAPPING ARRAY
    allocate(map_i(lnd%ls:lnd%le))
    allocate(map_j(lnd%ls:lnd%le))
    map_i = -1
    map_j = -1
    if( .not. input_is_multiface ) then
       allocate(mask(size(in_lon),size(in_lat)))
       mask = .false.

       if(fill_land_mask) then
          ! READ THE FIRST RECORD AND CALCULATE THE MASK OF THE VALID INPUT DATA
          call get_dimension_size(fileobj, "tile", dimlens(3))
          call get_dimension_size(fileobj, "cohort", dimlens(4))

          ! Note: The input file used for initial testing had
          ! lon = 144, lat = 90, tile = 2, cohort = 1
          ndims = get_variable_num_dimensions(fileobj, "cohort_index")
          allocate(siz(ndims))
          call get_variable_size(fileobj, "cohort_index", siz)
          allocate(cidx(siz(1)), idata(siz(1)))
          deallocate(siz)
          call read_data(fileobj, 'cohort_index', cidx)
          call read_data(fileobj, 'species', idata, unlim_dim_level=1)
          do n = 1,size(cidx)
             m = cidx(n)
             i = modulo(m,dimlens(1))+1
             m = m/dimlens(1)
             j = modulo(m,dimlens(2))+1
             m = m/dimlens(2)
             ! k = modulo(m,dimlens(3))+1 ! This is how to get tile number, if it were needed.
             m = m/dimlens(3)
                ! L = m+1  ! This is how to get cohort number, if it were needed. No need to do
                           ! modulo with dimlens(4) because at this point m is always < dimlens(4)
             if(idata(n)>=0 .or. mask(i,j)) then
                mask(i,j) = .TRUE. ! If species exists in any cohort of this grid cell then mask is .TRUE.
             endif
          enddo
          deallocate(idata)
       else
          mask(:,:) = .TRUE.
       endif
    endif
    deallocate (in_lon,in_lat)
    if(allocated(mask)) deallocate(mask)
  endif

  if(write_static_veg) &
      call init_writing_static_veg()

  module_is_initialized = .true.
end subroutine static_vegn_init


! create output file for static vegetation
subroutine init_writing_static_veg()
  integer, allocatable :: tidx(:)
  integer :: k, l, csize, iret
  integer :: year, month, day, hour, minute, sec ! components of base date
  character(len=256) :: units ! units of time in the output file

  ! count all land tiles and determine the length of tile dimension
  ! sufficient for the current domain
  tile_dim_length = 0
  do l = lnd%ls, lnd%le
     k = nitems(land_tile_map(l))
     tile_dim_length = max(tile_dim_length,k)
  enddo

  ! [1.1] calculate the tile dimension length by taking the max across all domains
  call mpp_max(tile_dim_length)
  call gather_tile_index(vegn_tile_exists,tidx)
  call gather_cohort_index(tile_dim_length, cidx)

  csize = size(cidx)
  allocate(species(csize),status(csize),bl(csize), blv(csize), br(csize), bsw(csize), &
           bwood(csize), bliving(csize))

  call gather_tile_index(vegn_tile_exists, tidx)

  call get_base_date(year,month,day,hour,minute,sec)
  base_time = set_date(year, month, day, hour, minute, sec)
  write(units, 11) year, month, day, hour, minute, sec

    call create_tile_out_file(static_veg_file, 'static_veg_out.nc', tidx, tile_dim_length)
  call create_cohort_dimension_new(static_veg_file, cidx, 'static_veg_out.nc', tile_dim_length)

  call register_axis(static_veg_file, "time", unlimited)
  call register_field(static_veg_file, "time", "double", (/"time"/))
  call register_variable_attribute(static_veg_file, "time", "units", units)
  call register_variable_attribute(static_veg_file, "time", "calendar", &
                                   valid_calendar_types(get_calendar_type()))

  call register_restart_field(static_veg_file, "species", species, (/"cohort_index"/))
  call register_variable_attribute(static_veg_file, "species", "long_name", &
                                   "vegetation species")

  call register_restart_field(static_veg_file, "bl", bl, (/"cohort_index"/))
  call register_variable_attribute(static_veg_file, "bl", "long_name", &
                                   "biomass of leaves per individual")
  call register_variable_attribute(static_veg_file, "bl", "units", &
                                   "kg C/m2")

  call register_restart_field(static_veg_file, "blv", blv, (/"cohort_index"/))
  call register_variable_attribute(static_veg_file, "blv", "long_name", &
                                   "biomass of virtual leaves (labile store) per individual")
  call register_variable_attribute(static_veg_file, "blv", "units", &
                                   "kg C/m2")

  call register_restart_field(static_veg_file, "br", br, (/"cohort_index"/))
  call register_variable_attribute(static_veg_file, "br", "long_name", &
                                   "biomass of fine roots per individual")
  call register_variable_attribute(static_veg_file, "br", "units", &
                                   "kg C/m2")

  call register_restart_field(static_veg_file, "bsw", bsw, (/"cohort_index"/))
  call register_variable_attribute(static_veg_file, "bsw", "long_name", &
                                   "biomass of sapwood per individual")
  call register_variable_attribute(static_veg_file, "bsw", "units", &
                                   "kg C/m2")

  call register_restart_field(static_veg_file, "bwood", bwood, (/"cohort_index"/))
  call register_variable_attribute(static_veg_file, "bwood", "long_name", &
                                   "biomass of heartwood per individual")
  call register_variable_attribute(static_veg_file, "bwood", "units", &
                                   "kg C/m2")

  call register_restart_field(static_veg_file, "bliving", bliving, (/"cohort_index"/))
  call register_variable_attribute(static_veg_file, "bliving", "long_name", &
                                   "total living biomass per individual")
  call register_variable_attribute(static_veg_file, "bliving", "units", &
                                   "")

  call register_restart_field(static_veg_file, "status", status, (/"cohort_index"/))
  call register_variable_attribute(static_veg_file, "status", "long_name", &
                                   "leaf status")
  call register_variable_attribute(static_veg_file, "status", "units", &
                                   "")

11 format('days since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
end subroutine init_writing_static_veg

! ===========================================================================
subroutine static_vegn_end()

  if (use_static_veg) then
    if (input_is_multiface) then
      call close_file(fileobj_domainug)
    else
      call close_file(fileobj)
    endif
    call close_file(static_veg_file)

  endif
  module_is_initialized = .false.
end subroutine static_vegn_end

! ===========================================================================
subroutine read_static_vegn (time, err_msg)
  type(time_type), intent(in)    :: time
  character(len=*), intent(out), optional :: err_msg

  ! ---- local vars
  integer :: index1, index2 ! result of time interpolation (only index1 is used)
  real    :: weight         ! another result of time interp, not used
  character(len=256) :: msg
  integer,dimension(:), allocatable :: siz
  integer, allocatable :: cidx(:), idata(:)
  real,    allocatable :: rdata(:)
  integer :: ndims

  if(.not.use_static_veg)return;

  msg = ''
  !   time_interp to find out the index of the current time interval
  if (timeline == 'loop') then
     call time_interp(time, ts, te, time_line, weight, index1, index2, &
                      correct_leap_year_inconsistency=.true., err_msg=msg)
  else if (timeline == 'normal') then
     call time_interp(time, time_line, weight, index1, index2, err_msg=msg)
  else
     call error_mesg(module_name,'timeline option "'//trim(timeline)// &
          '" is incorrect, use "normal" or "loop"', FATAL)
  endif
  if(msg /= '') then
    if(fms_error_handler('read_static_vegn','Message from time_interp: '//trim(msg),err_msg)) return
  endif

  ! read the data into cohort variables
  if(input_is_multiface) then
   ndims = get_variable_num_dimensions(fileobj, "cohort_index")
   allocate(siz(ndims))
   call get_variable_size(fileobj_domainug, "cohort_index", siz)
   allocate(cidx(siz(1)), idata(siz(1)), rdata(siz(1)))
   deallocate(siz)
   call read_data(fileobj_domainug, "cohort_index", cidx, unlim_dim_level=index1)
   call read_data(fileobj_domainug, "species", idata, unlim_dim_level=index1)
   call read_remap_cohort_data_i0d_new(fileobj_domainug, "species", cohort_species_ptr, map_i, map_j, cidx, idata)
   call read_data(fileobj_domainug, "bl", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj_domainug, "bl", cohort_bl_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj_domainug, "blv", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj_domainug, "blv", cohort_blv_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj_domainug, "br", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj_domainug, "br", cohort_br_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj_domainug, "bsw", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj_domainug, "bsw", cohort_bsw_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj_domainug, "bwood", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj_domainug, "bwood", cohort_bwood_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj_domainug, "bliving", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj_domainug, "bliving", cohort_bliving_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj_domainug, "status", idata, unlim_dim_level=index1)
   call read_remap_cohort_data_i0d_new(fileobj_domainug, "status", cohort_status_ptr, map_i, map_j, cidx, idata)
   deallocate(cidx, idata, rdata)
  else
   ndims = get_variable_num_dimensions(fileobj, "cohort_index")
   allocate(siz(ndims))
   call get_variable_size(fileobj, "cohort_index", siz)
   allocate(cidx(siz(1)), idata(siz(1)), rdata(siz(1)))
   deallocate(siz)
   call read_data(fileobj, "cohort_index", cidx)
   call read_data(fileobj, "species", idata, unlim_dim_level=index1)
   call read_remap_cohort_data_i0d_new(fileobj, "species", cohort_species_ptr, map_i, map_j, cidx, idata)
   call read_data(fileobj, "bl", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj, "bl", cohort_bl_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj, "blv", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj, "blv", cohort_blv_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj, "br", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj, "br", cohort_br_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj, "bsw", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj, "bsw", cohort_bsw_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj, "bwood", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj, "bwood", cohort_bwood_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj, "bliving", rdata, unlim_dim_level=index1)
   call read_remap_cohort_data_r0d_new(fileobj, "bliving", cohort_bliving_ptr, map_i, map_j, cidx, rdata)
   call read_data(fileobj, "status", idata, unlim_dim_level=index1)
   call read_remap_cohort_data_i0d_new(fileobj, "status", cohort_status_ptr, map_i, map_j, cidx, idata)
   deallocate(cidx, idata, rdata)
  endif

  ! derived variables will be updated in update_land_bc_fast
end subroutine read_static_vegn


! ===========================================================================
subroutine write_static_vegn()

  real :: t ! time in output units
  integer :: rec ! number of record to write
  ! components of the date
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1

  if(.not.write_static_veg) return;

  ! get components of calendar dates for this and previous time step
  call get_date(lnd%time,             year0,month0,day0,hour,minute,second)
  call get_date(lnd%time-lnd%dt_fast, year1,month1,day1,hour,minute,second)

  if (.not.((trim(static_veg_freq)=='daily'  .and.  day1/=day0)   &
        .or.(trim(static_veg_freq)=='monthly'.and.month1/=month0) &
        .or.(trim(static_veg_freq)=='annual' .and. year1/=year0))) return

  t = (time_type_to_real(lnd%time)-time_type_to_real(base_time))/86400
  call gather_cohort_data(cohort_species_ptr,cidx,tile_dim_length,species)
  call gather_cohort_data(cohort_bl_ptr,cidx,tile_dim_length,bl)
  call gather_cohort_data(cohort_blv_ptr,cidx,tile_dim_length,blv)
  call gather_cohort_data(cohort_br_ptr,cidx,tile_dim_length,br)
  call gather_cohort_data(cohort_bsw_ptr,cidx,tile_dim_length,bsw)
  call gather_cohort_data(cohort_bwood_ptr,cidx,tile_dim_length,bwood)
  call gather_cohort_data(cohort_bliving_ptr,cidx,tile_dim_length,bliving)
  call gather_cohort_data(cohort_status_ptr,cidx,tile_dim_length,status)

  call write_restart(static_veg_file)
! call write_restart(static_veg_file, unlim_dim_level=t)

end subroutine write_static_vegn

! ============================================================================
#define F90_TYPE       integer
#define READ_REMAP_SUB read_remap_cohort_data_i0d_new
#include "read_remap_cohort_data_new.inc"

#define F90_TYPE       real
#define READ_REMAP_SUB read_remap_cohort_data_r0d_new
#include "read_remap_cohort_data_new.inc"
! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
end function vegn_tile_exists

! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_COHORT_ACCESSOR(xtype,x) subroutine cohort_ ## x ## _ptr(c,p);\
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
