module predefined_tiles_input_mod

use hdf5
use, intrinsic :: iso_c_binding
use constants_mod, only : pi
use fms_mod, only : string, error_mesg, FATAL
use land_data_mod, only : land_state_type, atmos_land_boundary_type, log_version
use land_debug_mod, only : land_error_message
use land_tile_mod, only : insert, new_land_tile_glac, new_land_tile_lake, new_land_tile_soil
use land_tile_mod, only : land_tile_list_type, land_tile_type
use tiling_input_types_mod, only : tile_parameters_type, metadata_predefined_type, &
       lake_predefined_type, glacier_predefined_type, soil_predefined_type
use time_manager_mod, only : get_date

implicit none
private

! ==== public interfaces =====================================================
public :: open_database_predefined_tiles
public :: close_database_predefined_tiles
public :: land_cover_cold_start_0d_predefined_tiles
public :: land_cover_warm_start_0d_predefined_tiles
public :: downscale_atmos
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'predefined_tiles_input_mod'
#include "../shared/version_variable.inc"


! ==== interfaces ======================================================
interface get_parameter_data
   module procedure get_parameter_data_1d_integer
   module procedure get_parameter_data_1d_real
   module procedure get_parameter_data_2d_integer
   module procedure get_parameter_data_2d_real
end interface

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine check_h5err(status,message,severity)
  integer,      intent(in)           :: status   ! HDF5 status
  character(*), intent(in), optional :: message  ! optional
  integer,      intent(in), optional :: severity ! error severity

  integer :: severity_
  integer :: ierr

  if (status .eq. 0) return

  severity_ = FATAL
  if (present(severity)) severity_ = severity

  call h5eprint_f(ierr)
  if (present(message)) then
     call error_mesg('tiling_input', message, severity_)
  else
     call error_mesg('tiling_input', 'HDF5 ERROR', severity_)
  endif
end subroutine

subroutine open_database_predefined_tiles(h5id,lnd)
  type(land_state_type),intent(in) :: lnd
  integer(hid_t), intent(out) :: h5id
  integer :: status
  character(100) :: filename

  !Initialize the fortran library
  call h5open_f(status)
  !call check_h5err(status)

  !Define the file name
  filename = 'INPUT/ptiles.face'// trim(string(lnd%ug_face)) // '.h5'

  !Open access to the model input database
  CALL h5fopen_f(filename,H5F_ACC_RDONLY_F,h5id, status)
  call check_h5err(status, message='Error opening file "'//trim(filename)//'"')

  !CALL h5fopen_f('INPUT/land_model_input_database.h5',H5F_ACC_RDONLY_F,h5id, status)
  !CALL h5fopen_f('INPUT/land_model_input_database.nc',H5F_ACC_RDONLY_F,h5id, status)

end subroutine open_database_predefined_tiles

subroutine close_database_predefined_tiles(h5id)

  integer(hid_t), intent(in) :: h5id
  integer :: status
  !Close access to the model input database
  call h5fclose_f(h5id,status)
    !call check_h5err(status)

  !Close the hdf5 library
  call h5close_f(status)
    !call check_h5err(status)

end subroutine close_database_predefined_tiles

subroutine load_group_into_memory(tile,is,js,h5id,buf_ptr,buf_len,image_ptr)

  integer,intent(in) :: tile,is,js
  integer(hid_t), intent(in) :: h5id
  type(c_ptr),intent(inout) :: buf_ptr
  integer(size_t),intent(inout) :: buf_len
  character(kind=c_char),intent(inout),allocatable,dimension(:),target :: image_ptr
  integer :: status
  integer(hid_t) :: grpid,cell_grpid,dstid
  character(100) :: cellid_string
  integer(hid_t) :: fapl
  integer(size_t), parameter :: memory_increment = 1000000

 !Open access to the group in the database that contains all the group information
 call h5gopen_f(h5id,"grid_data",grpid,status)
 call check_h5err(status)
 !Write the cell id to string
 !write(cellid_string,'(I10)') cellid
 !cellid_string = trim('g' // trim(adjustl(cellid_string)))
 cellid_string = 'tile:' // trim(string(tile)) // ',is:' // trim(string(js)) // ',js:' // trim(string(is))
 !The goal here is to load the desired group of the cellid into memory. This buffer
 !will then be sent to the land model core. However, there is no direct way to do this
 !with a group instead we have to:
 !1.Create a new file in memory
 !2.Copy the desired group to this new file
 !3.Use the HDF5 api to then load this new file as a buffer
 !Ensure that we are always working in memory
 call h5pcreate_f(H5P_FILE_ACCESS_F,fapl,status)
 call check_h5err(status)
 !Setting the third parameter to false ensures that we never write this file to disk
 call h5pset_fapl_core_f(fapl,memory_increment,.False.,status)
 call check_h5err(status)
 !Although we create this file it is always in memory. It never gets written to disk
 call h5fcreate_f("buffer.hdf5",H5F_ACC_TRUNC_F,dstid,status,access_prp=fapl)
 call check_h5err(status)
 !Close access to the property list
 call h5pclose_f(fapl,status)
 !call check_h5err(status)
 !Copy the group from the original database to the new file
 call h5ocopy_f(grpid,cellid_string,dstid,'data',status)
 call check_h5err(status)
 !Flush the file
 call h5fflush_f(dstid,H5F_SCOPE_GLOBAL_F,status)
 call check_h5err(status)
 !Determine the size of the desired group (which is now a file in memory...)
 buf_len = 0
 buf_ptr = C_NULL_PTR
 call h5fget_file_image_f(dstid,buf_ptr,int(0,size_t),status,buf_len)
 call check_h5err(status)
 !Load the entire new file (i.e., desired group) into memory
 allocate(image_ptr(1:buf_len))
 buf_ptr = c_loc(image_ptr(1)(1:1))
 call h5fget_file_image_f(dstid,buf_ptr,buf_len,status)
 call check_h5err(status)
 !Close the copied file (release memory)
 call h5fclose_f(dstid,status)
 !call check_h5err(status)
 !Close access to the grid data group
 call h5gclose_f(grpid,status)
 !call check_h5err(status)

end subroutine load_group_into_memory

subroutine open_image_file(buf_ptr,buf_len,image_ptr,dstid)

 integer(hid_t),intent(out) :: dstid
 type(c_ptr),intent(inout) :: buf_ptr
 integer(size_t),intent(inout) :: buf_len
 character(kind=c_char),intent(inout),allocatable,dimension(:),target :: image_ptr
 integer(hid_t) :: fapl
 integer :: status

 !This will all be don eon the land model side. The only purpose is to open access to
 !the buffer as if it is a file
 !Open a new fapl
 call h5pcreate_f(H5P_FILE_ACCESS_F,fapl,status)
 !call check_h5err(status)
 call h5pset_fapl_core_f(fapl,buf_len,.False.,status)
 !call check_h5err(status)
 !Assign the buffer to memory
 call h5pset_file_image_f(fapl,buf_ptr,buf_len,status)
 !call check_h5err(status)
 !Discard the buffer
 deallocate(image_ptr)
 buf_ptr = C_NULL_PTR
 !Open the file from memory
 dstid = 0
 call h5fopen_f("test_image",H5F_ACC_RDONLY_F,dstid,status,fapl)
 !call check_h5err(status)

end subroutine open_image_file

subroutine land_cover_cold_start_0d_predefined_tiles(tiles,lnd,l,h5id)

  type(land_tile_list_type) , intent(inout) :: tiles   ! list of tiles to insert new tiles into
  type(land_state_type)     , intent(in)    :: lnd     ! land geometry and indexing information
  integer                   , intent(in)    :: l       ! unstructured grid index of current gridcell
  integer(hid_t)            , intent(in)    :: h5id    ! input data set handle

  type(land_tile_type), pointer :: tile
  integer :: itile,tid,i_index,j_index,status
  integer(hid_t) :: dstid,cid
  type(tile_parameters_type) :: tile_parameters
  type(c_ptr) :: buf_ptr
  integer(size_t) :: buf_len
  character(kind=c_char),allocatable,dimension(:),target :: image_ptr
!   real :: lat,lon

  i_index = lnd%i_index(l)
  j_index = lnd%j_index(l)

  !Print out the current lat and lon
!   lon = 180.0*lnd%ug_lon(l)/pi
!   lat = 180.0*lnd%ug_lat(l)/pi
!   print*,"Initializing: ",lat,lon

  !Retrieve buffer and buffer length of desired group (I/O core)
  call load_group_into_memory(lnd%ug_face,i_index,j_index,h5id,buf_ptr,buf_len,image_ptr)

  !Use buffer and buffer length to open new image file (Land model core)
  call open_image_file(buf_ptr,buf_len,image_ptr,dstid)

  !Retrieve the group
  call h5gopen_f(dstid,'data',cid,status)
  !call check_h5err(status)

  !Metadata
  call retrieve_metadata(tile_parameters,cid)
  call retrieve_soil_parameters(tile_parameters,cid)
  call retrieve_lake_parameters(tile_parameters,cid)
  call retrieve_glacier_parameters(tile_parameters,cid)

  !Create the tiles
  do itile = 1,tile_parameters%metadata%ntile
     if (tile_parameters%metadata%frac(itile) .eq. 0.0)cycle
     tid = tile_parameters%metadata%tid(itile)
     select case (tile_parameters%metadata%ttype(itile))
     case(1)
        tile => new_land_tile_glac(                           &
                  tile_parameters%metadata%frac(itile),       & ! land area fraction
                  tid,                                        & ! kind (tag) of glacier tile
                  tile_parameters%glacier,                    & ! collection of glacier parameter arrays
                  tid,                                        & ! itile -- index in glacier parameters
                  tile_parameters%metadata%tile(itile)+1,     & ! pid
                  tile_parameters%metadata%dws_prec(itile,:), &
                  tile_parameters%metadata%dws_srad(itile,:), &
                  tile_parameters%metadata%dws_tavg(itile,:)  &
              )
     case(2)
        tile => new_land_tile_lake(                           &
                  tile_parameters%metadata%frac(itile),       & ! land area fraction
                  tid,                                        & ! kind (tag) of lake tile
                  tile_parameters%lake,                       & ! collection of lake parameter arrays
                  tid,                                        & ! itile  -- index in lake parameters
                  tile_parameters%metadata%tile(itile)+1,     & ! pid
                  tile_parameters%metadata%dws_prec(itile,:), &
                  tile_parameters%metadata%dws_srad(itile,:), &
                  tile_parameters%metadata%dws_tavg(itile,:)  &
              )
     case(3)
        tile => new_land_tile_soil(                           &
                  tile_parameters%metadata%frac(itile),       & ! land area fraction
                  1,                                          & ! kind (tag) of soil tile
                  tile_parameters%soil%vegn(tid),             & ! kind (tag) of vegn tile
                  tile_parameters%soil,                       & ! collection of soil parameter arrays
                  tid,                                        & ! itile  -- index in soil parameters
                  tile_parameters%metadata%tile(itile)+1,     & ! pid
                  tile_parameters%metadata%dws_prec(itile,:), &
                  tile_parameters%metadata%dws_srad(itile,:), &
                  tile_parameters%metadata%dws_tavg(itile,:), &
                  tile_parameters%soil%hidx_j(tid),           &
                  tile_parameters%soil%hidx_k(tid)            &
              )
     case default
         call land_error_message('land_cover_cold_start_0d_predefined_tiles: unknow tile type',FATAL)
     end select
     call insert(tile,tiles)
  enddo

  !Close access to the grid cell group
  call h5gclose_f(cid,status)
  !call check_h5err(status)
  call h5fclose_f(dstid,status)
  !call check_h5err(status)

end subroutine land_cover_cold_start_0d_predefined_tiles

subroutine land_cover_warm_start_0d_predefined_tiles(lnd, h5id, tiles, l, &
  lidx, frac, pid, vegn)

  type(land_state_type)     , intent(in)    :: lnd     ! land geometry and indexing information
  integer(hid_t)            , intent(in)    :: h5id    ! input data set handle
  type(land_tile_list_type) , intent(inout) :: tiles   ! list of tiles to insert new tiles into
  integer                   , intent(in)    :: l       ! unstructured grid index of current gridcell
  integer                   , intent(in)    :: lidx(:) ! l-indices indices of input data
  real                      , intent(in)    :: frac(:) ! fractions of tiles
  integer                   , intent(in)    :: pid (:) ! "parent IDs" of tiles
  integer                   , intent(in)    :: vegn(:) ! vegetation flags of tiles

  type(land_tile_type), pointer :: tile
  integer :: itile,tid,i_index,j_index,status,k
  integer(hid_t) :: dstid,cid
  type(tile_parameters_type) :: tile_parameters
  type(c_ptr) :: buf_ptr
  integer(size_t) :: buf_len
  character(kind=c_char),allocatable,target :: image_ptr(:)
!   real :: lat,lon

  i_index = lnd%i_index(l)
  j_index = lnd%j_index(l)

  !Print out the current lat and lon
!   lon = 180.0*lnd%ug_lon(l)/pi
!   lat = 180.0*lnd%ug_lat(l)/pi
!   print*,"Initializing: ",lat,lon

  !Retrieve buffer and buffer length of desired group (I/O core)
  call load_group_into_memory(lnd%ug_face,i_index,j_index,h5id,buf_ptr,buf_len,image_ptr)

  !Use buffer and buffer length to open new image file (Land model core)
  call open_image_file(buf_ptr,buf_len,image_ptr,dstid)

  !Retrieve the group
  call h5gopen_f(dstid,'data',cid,status)
  !call check_h5err(status)

  call retrieve_metadata(tile_parameters,cid)
  call retrieve_soil_parameters(tile_parameters,cid)
  call retrieve_lake_parameters(tile_parameters,cid)
  call retrieve_glacier_parameters(tile_parameters,cid)

  ! Create the tiles
  do k = 1,size(frac)
     if (lidx(k)/=l) cycle ! skip all tiles that do not belong to the current point
     itile = pid(k)
     tid = tile_parameters%metadata%tid(itile)
     select case (tile_parameters%metadata%ttype(itile))
     case(1)
        tile => new_land_tile_glac(                           &
                  frac(k),                                    & ! land area fraction
                  tid,                                        & ! kind (tag) of glacier tile
                  tile_parameters%glacier,                    & ! collection of glacier parameter arrays
                  tid,                                        & ! itile -- index in glacier parameters
                  pid(k),                                     & ! pid
                  tile_parameters%metadata%dws_prec(itile,:), &
                  tile_parameters%metadata%dws_srad(itile,:), &
                  tile_parameters%metadata%dws_tavg(itile,:)  &
              )
     case(2)
        tile => new_land_tile_lake(                           &
                  frac(k),                                    & ! land area fraction
                  tid,                                        & ! kind (tag) of lake tile
                  tile_parameters%lake,                       & ! collection of lake parameter arrays
                  tid,                                        & ! itile  -- index in lake parameters
                  pid(k),                                     & ! pid
                  tile_parameters%metadata%dws_prec(itile,:), &
                  tile_parameters%metadata%dws_srad(itile,:), &
                  tile_parameters%metadata%dws_tavg(itile,:)  &
              )
     case(3)
        tile => new_land_tile_soil(                           &
                  frac(k),                                    & ! land area fraction
                  1,                                          & ! kind (tag) of soil tile
                  vegn(k),                                    & ! kind (tag) of vegn tile
                  tile_parameters%soil,                       & ! collection of soil parameter arrays
                  tid,                                        & ! itile  -- index in soil parameters
                  pid(k),                                     & ! pid
                  tile_parameters%metadata%dws_prec(itile,:), &
                  tile_parameters%metadata%dws_srad(itile,:), &
                  tile_parameters%metadata%dws_tavg(itile,:), &
                  tile_parameters%soil%hidx_j(tid),           &
                  tile_parameters%soil%hidx_k(tid)            &
              )
     case default
         call land_error_message('land_cover_cold_start_0d_predefined_tiles: unknow tile type',FATAL)
     end select
     call insert(tile,tiles)
  enddo

  !Close access to the grid cell's group
  call h5gclose_f(cid,status)
  !call check_h5err(status)
  call h5fclose_f(dstid,status)
  !call check_h5err(status)
end subroutine land_cover_warm_start_0d_predefined_tiles

subroutine retrieve_metadata(tile_parameters,cid)

  type(tile_parameters_type),intent(inout),target :: tile_parameters
  integer(hid_t),intent(inout) :: cid
  type(metadata_predefined_type),pointer :: metadata
  integer(hid_t) :: dimid,grpid,dsid,varid
  integer :: ntile,nband,status,im
  integer(hsize_t) :: dims(1),maxdims(1)
  !allocate(tile_parameters%metadata)
  metadata => tile_parameters%metadata

  !Retrieve the group id
  call h5gopen_f(cid,"metadata",grpid,status)
  !call check_h5err(status)

  !Retrieve the number of tiles
  call h5dopen_f(grpid,"tile",varid,status)
  !call check_h5err(status)
  call h5dget_space_f(varid,dsid,status)
  !call check_h5err(status)
  call h5sget_simple_extent_dims_f(dsid,dims,maxdims,status)
  !call check_h5err(status)
  call h5dclose_f(varid,status)
  !call check_h5err(status)
  metadata%ntile = dims(1)

  !Retrieve the info
  call get_parameter_data(grpid,'frac',metadata%ntile,metadata%frac)
  call get_parameter_data(grpid,'tile',metadata%ntile,metadata%tile)
  call get_parameter_data(grpid,'type',metadata%ntile,metadata%ttype)
  call get_parameter_data(grpid,'tid',metadata%ntile,metadata%tid)

  !Retrieve meteorology info
  call get_parameter_data(grpid,'prec',metadata%ntile,12,metadata%dws_prec)
  call get_parameter_data(grpid,'srad',metadata%ntile,12,metadata%dws_srad)
  call get_parameter_data(grpid,'tavg',metadata%ntile,12,metadata%dws_tavg)

  !Curate the weights
  do im = 1,12
   if (isnan(sum(metadata%frac*metadata%dws_prec(:,im))))metadata%dws_prec(:,im)=1
   if (isnan(sum(metadata%frac*metadata%dws_tavg(:,im))))metadata%dws_tavg(:,im)=1
  enddo
  !print*,'prec',sum(metadata%frac*metadata%dws_prec(:,1))
  !print*,'tavg',sum(metadata%frac*metadata%dws_tavg(:,1))

  !Clean up (This should go in the database creation)
  where ((metadata%frac .lt. 1.e-8) .and. (metadata%ttype .eq. 2))metadata%frac = 0.0
  metadata%frac = metadata%frac/sum(metadata%frac)

  !Close access to the group
  call h5gclose_f(grpid,status)
  !call check_h5err(status)

end subroutine retrieve_metadata

subroutine retrieve_glacier_parameters(tile_parameters,cid)

  type(tile_parameters_type),intent(inout),target :: tile_parameters
  integer(hid_t),intent(inout) :: cid
  type(glacier_predefined_type),pointer :: glacier
  integer(hid_t) :: dimid,grpid,dsid,varid
  integer :: nglacier,nband,status
  integer(hsize_t) :: dims(2),maxdims(2)
  !allocate(tile_parameters%glacier)
  glacier => tile_parameters%glacier

  !Retrieve the group id
  call h5gopen_f(cid,"glacier",grpid,status)
  !!call check_h5err(status)

  !If it does not exist then exit
  if (status .eq. -1)return

  !Retrieve the number of lake tiles and bands
  call h5dopen_f(grpid,"refl_min_dir",varid,status)
  !call check_h5err(status)
  call h5dget_space_f(varid,dsid,status)
  !call check_h5err(status)
  call h5sget_simple_extent_dims_f(dsid,dims,maxdims,status)

  !call check_h5err(status)
  nglacier = dims(1)
  nband = dims(2)
  glacier%nglacier = nglacier
  call h5sclose_f(dsid,status)
  !call check_h5err(status)
  call h5dclose_f(varid,status)
  !call check_h5err(status)

  !Retrieve the parameters
  call get_parameter_data(grpid,"w_sat",nglacier,glacier%w_sat)
  call get_parameter_data(grpid,"awc_lm2",nglacier,glacier%awc_lm2)
  call get_parameter_data(grpid,"k_sat_ref",nglacier,glacier%k_sat_ref)
  call get_parameter_data(grpid,"psi_sat_ref",nglacier,glacier%psi_sat_ref)
  call get_parameter_data(grpid,"chb",nglacier,glacier%chb)
  call get_parameter_data(grpid,"alpha",nglacier,glacier%alpha)
  call get_parameter_data(grpid,"heat_capacity_ref",nglacier,glacier%heat_capacity_ref)
  call get_parameter_data(grpid,"thermal_cond_ref",nglacier,glacier%thermal_cond_ref)
  call get_parameter_data(grpid,"emis_dry",nglacier,glacier%emis_dry)
  call get_parameter_data(grpid,"emis_sat",nglacier,glacier%emis_sat)
  call get_parameter_data(grpid,"z0_momentum",nglacier,glacier%z0_momentum)
  call get_parameter_data(grpid,"tfreeze",nglacier,glacier%tfreeze)
  call get_parameter_data(grpid,"refl_max_dir",nglacier,nband,glacier%refl_max_dir)
  call get_parameter_data(grpid,"refl_max_dif",nglacier,nband,glacier%refl_max_dif)
  call get_parameter_data(grpid,"refl_min_dir",nglacier,nband,glacier%refl_min_dir)
  call get_parameter_data(grpid,"refl_min_dif",nglacier,nband,glacier%refl_min_dif)

  !Close access to the group
  call h5gclose_f(grpid,status)
  !call check_h5err(status)

end subroutine retrieve_glacier_parameters

subroutine retrieve_lake_parameters(tile_parameters,cid)

  type(tile_parameters_type),intent(inout),target :: tile_parameters
  integer(hid_t), intent(inout) :: cid
  type(lake_predefined_type),pointer :: lake
  integer(hid_t) :: dimid,grpid,dsid,varid
  integer :: nlake,nband,status
  integer(hsize_t) :: dims(2),maxdims(2)
  !allocate(tile_parameters%lake)
  lake => tile_parameters%lake

  !Retrieve the group id
  call h5gopen_f(cid,"lake",grpid,status)
  !!call check_h5err(status)

  !If it does not exist then exit
  if (status .eq. -1)return

  !Retrieve the number of lake tiles and bands
  call h5dopen_f(grpid,"refl_dry_dir",varid,status)
  !call check_h5err(status)
  call h5dget_space_f(varid,dsid,status)
  !call check_h5err(status)
  call h5sget_simple_extent_dims_f(dsid,dims,maxdims,status)
  !call check_h5err(status)
  nlake = dims(1)
  nband = dims(2)
  lake%nlake = nlake
  call h5sclose_f(dsid,status)
  !call check_h5err(status)
  call h5dclose_f(varid,status)
  !call check_h5err(status)

  !Retrieve the parameters
  call get_parameter_data(grpid,"connected_to_next",nlake,lake%connected_to_next)
  call get_parameter_data(grpid,"whole_lake_area",nlake,lake%whole_area)
  call get_parameter_data(grpid,"lake_depth_sill",nlake,lake%depth_sill)
  call get_parameter_data(grpid,"lake_width_sill",nlake,lake%width_sill)
  call get_parameter_data(grpid,"lake_backwater",nlake,lake%backwater)
  call get_parameter_data(grpid,"lake_backwater_1",nlake,lake%backwater_1)
  call get_parameter_data(grpid,"awc_lm2",nlake,lake%awc_lm2)
  call get_parameter_data(grpid,"w_sat",nlake,lake%w_sat)
  call get_parameter_data(grpid,"chb",nlake,lake%chb)
  call get_parameter_data(grpid,"refl_dry_dif",nlake,nband,lake%refl_dry_dif)
  call get_parameter_data(grpid,"refl_dry_dir",nlake,nband,lake%refl_dry_dir)
  call get_parameter_data(grpid,"refl_sat_dif",nlake,nband,lake%refl_sat_dif)
  call get_parameter_data(grpid,"refl_sat_dir",nlake,nband,lake%refl_sat_dir)
  call get_parameter_data(grpid,"z0_momentum",nlake,lake%z0_momentum)
  call get_parameter_data(grpid,"psi_sat_ref",nlake,lake%psi_sat_ref)
  call get_parameter_data(grpid,"emis_dry",nlake,lake%emis_dry)
  call get_parameter_data(grpid,"z0_momentum_ice",nlake,lake%z0_momentum_ice)
  call get_parameter_data(grpid,"heat_capacity_ref",nlake,lake%heat_capacity_ref)
  call get_parameter_data(grpid,"alpha",nlake,lake%alpha)
  call get_parameter_data(grpid,"thermal_cond_ref",nlake,lake%thermal_cond_ref)
  call get_parameter_data(grpid,"emis_sat",nlake,lake%emis_sat)
  call get_parameter_data(grpid,"k_sat_ref",nlake,lake%k_sat_ref)

  !Close access to the group
  call h5gclose_f(grpid,status)
  !call check_h5err(status)

end subroutine retrieve_lake_parameters

subroutine retrieve_soil_parameters(tile_parameters,cid)

  type(tile_parameters_type), intent(inout),target :: tile_parameters
  integer(hid_t), intent(inout) :: cid
  type(soil_predefined_type),pointer :: soil
  integer(hid_t) :: varid,grpid,dsid
  integer :: nsoil,nband,status
  integer(hsize_t) :: dims(2),maxdims(2)
  !allocate(tile_parameters%soil)
  soil => tile_parameters%soil

  !Retrieve the group id
  call h5gopen_f(cid,"soil",grpid,status)
  !!call check_h5err(status)

  !If it does not exist then exit
  if (status .eq. -1)return

  !Retrieve the number of soil tiles and bands
  call h5dopen_f(grpid,"dat_refl_dry_dir",varid,status)
  !call check_h5err(status)
  call h5dget_space_f(varid,dsid,status)
  !call check_h5err(status)
  call h5sget_simple_extent_dims_f(dsid,dims,maxdims,status)
  !call check_h5err(status)
  nsoil = dims(1)
  nband = dims(2)
  soil%nsoil = nsoil
  call h5sclose_f(dsid,status)
  !call check_h5err(status)
  call h5dclose_f(varid,status)
  !call check_h5err(status)

  !Retrieve the parameters
  call get_parameter_data(grpid,"dat_w_sat",nsoil,soil%dat_w_sat)
  call get_parameter_data(grpid,"dat_awc_lm2",nsoil,soil%dat_awc_lm2)
  call get_parameter_data(grpid,"dat_k_sat_ref",nsoil,soil%dat_k_sat_ref)
  call get_parameter_data(grpid,"dat_psi_sat_ref",nsoil,soil%dat_psi_sat_ref)
  call get_parameter_data(grpid,"dat_chb",nsoil,soil%dat_chb)
  call get_parameter_data(grpid,"dat_heat_capacity_dry",nsoil,soil%dat_heat_capacity_dry)
  call get_parameter_data(grpid,"dat_thermal_cond_dry",nsoil,soil%dat_thermal_cond_dry)
  call get_parameter_data(grpid,"dat_thermal_cond_sat",nsoil,soil%dat_thermal_cond_sat)
  call get_parameter_data(grpid,"dat_thermal_cond_exp",nsoil,soil%dat_thermal_cond_exp)
  call get_parameter_data(grpid,"dat_thermal_cond_scale",nsoil,soil%dat_thermal_cond_scale)
  call get_parameter_data(grpid,"dat_thermal_cond_weight",nsoil,soil%dat_thermal_cond_weight)
  call get_parameter_data(grpid,"dat_refl_dry_dir",nsoil,nband,soil%dat_refl_dry_dir)
  call get_parameter_data(grpid,"dat_refl_dry_dif",nsoil,nband,soil%dat_refl_dry_dif)
  call get_parameter_data(grpid,"dat_refl_sat_dir",nsoil,nband,soil%dat_refl_sat_dir)
  call get_parameter_data(grpid,"dat_refl_sat_dif",nsoil,nband,soil%dat_refl_sat_dif)
  call get_parameter_data(grpid,"dat_emis_dry",nsoil,soil%dat_emis_dry)
  call get_parameter_data(grpid,"dat_emis_sat",nsoil,soil%dat_emis_sat)
  call get_parameter_data(grpid,"dat_z0_momentum",nsoil,soil%dat_z0_momentum)
  call get_parameter_data(grpid,"dat_tf_depr",nsoil,soil%dat_tf_depr)
  call get_parameter_data(grpid,"rsa_exp_global",nsoil,soil%rsa_exp_global)
  call get_parameter_data(grpid,"gw_res_time",nsoil,soil%gw_res_time)
  call get_parameter_data(grpid,"gw_hillslope_length",nsoil,soil%gw_hillslope_length)
  call get_parameter_data(grpid,"gw_scale_length",nsoil,soil%gw_scale_length)
  call get_parameter_data(grpid,"gw_hillslope_zeta_bar",nsoil,soil%gw_hillslope_zeta_bar)
  call get_parameter_data(grpid,"gw_hillslope_relief",nsoil,soil%gw_hillslope_relief)
  call get_parameter_data(grpid,"gw_scale_relief",nsoil,soil%gw_scale_relief)
  call get_parameter_data(grpid,"gw_soil_e_depth",nsoil,soil%gw_soil_e_depth)
  call get_parameter_data(grpid,"gw_scale_soil_depth",nsoil,soil%gw_scale_soil_depth)
  !call get_parameter_data(grpid,"gw_hillslope_a",nsoil,soil%gw_hillslope_a)
  !call get_parameter_data(grpid,"gw_hillslope_n",nsoil,soil%gw_hillslope_n)
  call get_parameter_data(grpid,"gw_perm",nsoil,soil%gw_perm)
  call get_parameter_data(grpid,"gw_scale_perm",nsoil,soil%gw_scale_perm)
  call get_parameter_data(grpid,"gw_res_time",nsoil,soil%gw_res_time)
  call get_parameter_data(grpid,"microtopo",nsoil,soil%microtopo)
  call get_parameter_data(grpid,"tile_hlsp_length",nsoil,soil%tile_hlsp_length)
  call get_parameter_data(grpid,"tile_hlsp_slope",nsoil,soil%tile_hlsp_slope)
  call get_parameter_data(grpid,"tile_hlsp_elev",nsoil,soil%tile_hlsp_elev)
  call get_parameter_data(grpid,"tile_hlsp_hpos",nsoil,soil%tile_hlsp_hpos)
  call get_parameter_data(grpid,"tile_hlsp_width",nsoil,soil%tile_hlsp_width)
  call get_parameter_data(grpid,"tile_hlsp_frac",nsoil,soil%tile_hlsp_frac)
  call get_parameter_data(grpid,"hidx_k",nsoil,soil%hidx_k)
  call get_parameter_data(grpid,"hidx_j",nsoil,soil%hidx_j)
  call get_parameter_data(grpid,"vegn",nsoil,soil%vegn)
  call get_parameter_data(grpid,"landuse",nsoil,soil%landuse)
  call get_parameter_data(grpid,"irrigation",nsoil,soil%irrigation)
  call get_parameter_data(grpid,"bl",nsoil,soil%bl)
  call get_parameter_data(grpid,"bsw",nsoil,soil%bsw)
  call get_parameter_data(grpid,"bwood",nsoil,soil%bwood)
  call get_parameter_data(grpid,"br",nsoil,soil%br)
  call get_parameter_data(grpid,"wtd",nsoil,soil%iwtd)
  call get_parameter_data(grpid,"soil_depth",nsoil,soil%soil_depth)
  call get_parameter_data(grpid,"depth_to_bedrock",nsoil,soil%depth_to_bedrock)
  call get_parameter_data(grpid,"hand_ecdf",nsoil,11,soil%hand_ecdf)
  call get_parameter_data(grpid,"hand_bedges",nsoil,11,soil%hand_bedges)
  call get_parameter_data(grpid,"ksat_0cm",nsoil,soil%ksat0cm)
  call get_parameter_data(grpid,"ksat_200cm",nsoil,soil%ksat200cm)

  !Do some basic QC (should be done in the preprocessing...)
  where(isnan(soil%gw_soil_e_depth))soil%gw_soil_e_depth = 3.0
  !where(soil%gw_soil_e_depth .lt. 1.0)soil%gw_soil_e_depth = 1.0
  where(isnan(soil%gw_perm))soil%gw_perm = 2e-13 !HACK
  where(isnan(soil%soil_depth))soil%soil_depth = 2.0
  where(isnan(soil%depth_to_bedrock))soil%depth_to_bedrock = 10.0
  soil%gw_soil_e_depth = soil%soil_depth !Remove dependence on input soil_e_depth
  !Turn everything below DTB to hard rock...
  !soil%gw_perm = 10**(-15.2)!1e-20
  !Divide soil_e_depth by 3
  !soil%gw_soil_e_depth = soil%gw_soil_e_depth/100.0
  where(soil%dat_k_sat_ref .gt. 0.035)soil%dat_k_sat_ref = 0.035 !HACK
  where(soil%dat_psi_sat_ref .gt. -0.01)soil%dat_psi_sat_ref = -0.01 !HACK
  where(soil%dat_chb .lt. 2.0)soil%dat_chb = 2.0 !HACK
  !where(soil%dat_k_sat_ref .gt. 0.01)soil%dat_k_sat_ref = 0.01 !HACK
  !where(soil%dat_k_sat_ref .gt. 0.025)soil%dat_k_sat_ref = 0.0 !HACK
  !Correct psi_sat_ref
  !soil%dat_psi_sat_ref = 10*soil%dat_psi_sat_ref !Need to move back to original
  !where(soil%dat_psi_sat_ref .gt. -0.01)soil%dat_psi_sat_ref = -0.01 !HACK
  where(isnan(soil%dat_k_sat_ref))soil%dat_k_sat_ref = 0.003 !HACK
  where(isnan(soil%bl))soil%bl = 0.0
  where(isnan(soil%br))soil%br = 0.0
  where(isnan(soil%bsw))soil%bsw = 0.0
  where(soil%bwood .gt. 10**10)soil%bwood = 0.0
  where(soil%bl .gt. 10**10)soil%bl = 0.0
  where(soil%br .gt. 10**10)soil%br = 0.0
  where(soil%bsw .gt. 10**10)soil%bsw = 0.0
  where(soil%bwood .gt. 10**10)soil%bwood = 0.0
  !microtopography
  soil%microtopo = 9e9 !meter (dh/2) !0.5 ! meter

  !Close access to the group
  call h5gclose_f(grpid,status)
  !call check_h5err(status)

end subroutine retrieve_soil_parameters

subroutine get_parameter_data_1d_integer(grpid,var,nx,tmp)

 integer(hid_t), intent(in) :: grpid
 character(len=*),intent(in) :: var
 integer,intent(in) :: nx
 !integer,dimension(:),pointer :: tmp
 integer,dimension(:),allocatable :: tmp
 integer(hid_t) :: varid
 integer :: status
 integer(hsize_t) :: dims(1)
 if (allocated(tmp))deallocate(tmp)
 allocate(tmp(nx))

 dims(1) = nx
 call h5dopen_f(grpid,var,varid,status)
 !call check_h5err(status)
 call h5dread_f(varid,H5T_STD_I32LE,tmp,dims,status)
 !call check_h5err(status)
 call h5dclose_f(varid,status)
 !call check_h5err(status)

end subroutine

subroutine get_parameter_data_1d_real(grpid,var,nx,tmp)

 integer(hid_t), intent(in) :: grpid
 character(len=*),intent(in) :: var
 integer,intent(in) :: nx
 !real,dimension(:),pointer :: tmp
 real,dimension(:),allocatable :: tmp
 integer(hid_t) :: varid
 integer :: status
 integer(hsize_t) :: dims(1)
 !if (associated(tmp))deallocate(tmp)
 if (allocated(tmp))deallocate(tmp)
 allocate(tmp(nx))

 dims(1) = nx
 call h5dopen_f(grpid,var,varid,status)
 !call check_h5err(status)
 call h5dread_f(varid,H5T_IEEE_F64LE,tmp,dims,status)
 !call check_h5err(status)
 call h5dclose_f(varid,status)
 !call check_h5err(status)

end subroutine

subroutine get_parameter_data_2d_integer(grpid,var,nx,ny,tmp)

 integer(hid_t), intent(in) :: grpid
 character(len=*),intent(in) :: var
 integer,intent(in) :: nx,ny
 integer,dimension(:,:), allocatable :: tmp
 integer(hid_t) :: varid
 integer :: status
 integer(hsize_t) :: dims(2)
 if (allocated(tmp))deallocate(tmp)
 allocate(tmp(nx,ny))

 dims(1) = nx
 dims(2) = ny
 call h5dopen_f(grpid,var,varid,status)
 !call check_h5err(status)
 call h5dread_f(varid,H5T_STD_I32LE,tmp,dims,status)
 !call check_h5err(status)
 call h5dclose_f(varid,status)
 !call check_h5err(status)

end subroutine

subroutine get_parameter_data_2d_real(grpid,var,nx,ny,tmp)

 integer(hid_t),   intent(in) :: grpid
 character(len=*), intent(in) :: var
 integer,          intent(in) :: nx,ny
 real,dimension(:,:),allocatable :: tmp
 integer(hid_t) :: varid
 integer :: itile,status
 integer(hsize_t) :: dims(2)
 if (allocated(tmp))deallocate(tmp)
 allocate(tmp(nx,ny))

 dims(1) = nx
 dims(2) = ny
 call h5dopen_f(grpid,var,varid,status)
 !call check_h5err(status)
 call h5dread_f(varid,H5T_IEEE_F64LE,tmp,dims,status)
 !call check_h5err(status)
 call h5dclose_f(varid,status)
 !call check_h5err(status)

end subroutine

subroutine downscale_atmos(tile,cplr2land,l,k,lnd)

  implicit none
  type(land_state_type),intent(in) :: lnd
  type(land_tile_type), intent(inout) :: tile ! pointer to current tile
  type(atmos_land_boundary_type), intent(in)    :: cplr2land
  integer, intent(in) :: l,k
  integer :: year,month,day,hour,minute,second
  real :: tas_ds,lprec_ds,fprec_ds,e_res_ds,hlf,tfreeze
  e_res_ds = 0.0
  hlf = 3.34e5 !< Latent heat of fusion [J/kg]
  tfreeze = 273.16 !< Freezing temperature of fresh water [K]

  call get_date(lnd%time,year,month,day,hour,minute,second)

  !Downscale temperature (precipitation temperature)
  tas_ds = tile%dws_tavg(month)*cplr2land%tprec(l,k)
  !Get canopy temperature
  !tas_ds = tile%cana%T

  !Downscale precipitation
  lprec_ds = tile%dws_prec(month)*cplr2land%lprec(l,k)
  fprec_ds = tile%dws_prec(month)*cplr2land%fprec(l,k)

  !Melt any frozen precipitation
  !e_res_ds = e_res_ds - HLF*fprec_ds
  lprec_ds = lprec_ds + fprec_ds
  fprec_ds = 0.0

  !Calculate the heat necessary to bring up the temperature to tfreeze
  !Calculate the max mass of water that could be frozen before surpassing tfreeze
  !When applicable convert liquid water to frozen
  if (tas_ds .le. tfreeze)then
  ! e_res_ds = e_res_ds + HLF*lprec_ds
   fprec_ds = lprec_ds
   lprec_ds = 0.0
  endif

  !Send values back to the land model
  cplr2land%tprec(l,k) = tas_ds
  cplr2land%lprec(l,k) = lprec_ds
  cplr2land%fprec(l,k) = fprec_ds
  tile%e_res_ds = 0.0!e_res_ds

end subroutine

end module predefined_tiles_input_mod
