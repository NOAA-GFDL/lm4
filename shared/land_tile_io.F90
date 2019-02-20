module land_tile_io_mod

use netcdf, only: NF90_MAX_NAME, NF90_FILL_DOUBLE, NF90_FILL_INT
use fms_mod, only : error_mesg, FATAL, mpp_pe
use fms_io_mod, only : get_instance_filename

use fms2_io_mod, only: FmsNetcdfUnstructuredDomainFile_t, &
                       register_axis, register_field, &
                       register_variable_attribute, write_restart, &
                       close_file, variable_exists, get_variable_size, &
                       read_data, write_data, open_file, &
                       get_variable_num_dimensions, compressed_start_and_count

use time_manager_mod, only : time_type
use data_override_mod, only : data_override_ug
use mpp_domains_mod,   only : mpp_pass_SG_to_UG
use land_io_mod, only : read_field, input_buf_size
use land_tile_mod, only : land_tile_type, land_tile_list_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles, &
     tile_exists_func, fptr_i0, fptr_i0i, fptr_r0, fptr_r0i, fptr_r0ij, fptr_r0ijk, &
     land_tile_map

use land_data_mod, only  : lnd
use land_utils_mod, only : put_to_tiles_r0d_fptr


implicit none
private

! ==== public interfaces =====================================================
! restart i/o subroutines: those use CF "compression by gathering" technique
! to pack tile data.
public :: land_restart_type
public :: init_land_restart, open_land_restart, save_land_restart, free_land_restart
public :: add_restart_axis
public :: add_tile_data, add_int_tile_data, add_scalar_data, add_text_data
public :: get_tile_data, get_int_tile_data, get_scalar_data, get_text_data
public :: field_exists
public :: gather_tile_index

public :: read_field

public :: create_tile_out_file

! auxiliary subroutines
public :: get_tile_by_idx

! ==== end of public interfaces ==============================================
interface create_tile_out_file
   module procedure create_tile_out_file_idx_new
end interface

interface add_tile_data
   module procedure add_tile_data_r0d_fptr_r0
   module procedure add_tile_data_r0d_fptr_r0i
   module procedure add_tile_data_r0d_fptr_r0ij
   module procedure add_tile_data_r1d_fptr_r0i
   module procedure add_tile_data_r1d_fptr_r0ij
   module procedure add_tile_data_r1d_fptr_r0ijk
   module procedure add_tile_data_r2d_fptr_r0ij
   module procedure add_tile_data_r2d_fptr_r0ijk
end interface

interface add_int_tile_data
   module procedure add_tile_data_i0d_fptr_i0
   module procedure add_tile_data_i1d_fptr_i0i
end interface

interface get_tile_data
   module procedure get_tile_data_r0d_fptr_r0
   module procedure get_tile_data_r0d_fptr_r0i
   module procedure get_tile_data_r0d_fptr_r0ij
   module procedure get_tile_data_r1d_fptr_r0i
   module procedure get_tile_data_r1d_fptr_r0ij
   module procedure get_tile_data_r1d_fptr_r0ijk
   module procedure get_tile_data_r2d_fptr_r0ij
   module procedure get_tile_data_r2d_fptr_r0ijk
end interface

interface get_int_tile_data
   module procedure get_tile_data_i0d_fptr_i0
   module procedure get_tile_data_i1d_fptr_i0i
end interface

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_tile_io_mod'

! name of the "compressed" dimension (and dimension variable) in the output
! netcdf files -- that is, the dimensions written out using compression by
! gathering, as described in CF conventions. See subroutines write_tile_data,
! read_tile_data, read_unpack_tile_data, write_cohort_data
character(len=*), parameter :: tile_index_name = 'tile_index'

! ==== module types ==========================================================
type axis
   character(128) :: name = ''  ! name of the axis
   integer        :: len = -1   ! length of the axis
end type axis
! land restart type encapsulates the data needed for the land restarts
type land_restart_type
   type(FmsNetcdfUnstructuredDomainFile_t) :: rhandle ! fms_io restart file data type
   logical :: should_free_rhandle = .FALSE.

   character(267) :: basename ='' ! name of the restart file
   character(267) :: filename ='' ! name of the restart file after adding PE number and such
   integer, allocatable :: tidx(:) ! tile index
   integer, allocatable :: cidx(:) ! vegetation cohort index
   integer :: tile_dim_length = -1! length of tile dimension

   ! axis information
   integer    :: nax=0
   type(axis) :: ax(5)
end type land_restart_type

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ==============================================================================
! given a file name and tile indicator funcrion, creates land restart structure
subroutine init_land_restart(restart,filename,tile_exists,tile_dim_length)
  type(land_restart_type), intent(out) :: restart
  character(*),            intent(in)  :: filename ! name of the file to create
  procedure(tile_exists_func) :: tile_exists ! existence detector function:
      ! returns true if specific tile exists (hence should be written to restart)
  integer,                 intent(in)  :: tile_dim_length ! length of the tile dimension
      ! in the output file

  restart%basename=filename
  ! TODO: determine tile_dim_length inside this subroutine. It is equal to the
  ! max number of tiles per grid cell
  restart%tile_dim_length = tile_dim_length
  ! allocate and fill tile compression index
  call gather_tile_index(tile_exists,restart%tidx)

  call create_tile_out_file_idx_new(restart%rhandle,restart%basename,restart%tidx, &
                                    restart%tile_dim_length)
  restart%should_free_rhandle = .TRUE.
end subroutine init_land_restart

! ==============================================================================
subroutine open_land_restart(restart,filename,restart_exists)
  type(land_restart_type), intent(out) :: restart
  character(*),            intent(in)  :: filename
  logical,                 intent(out) :: restart_exists

  ! ---- local vars
  integer,dimension(:),allocatable :: flen ! length of the index
  integer :: ndims

  restart%basename = filename
  restart_exists = open_file(restart%rhandle, restart%basename, "read", &
                             lnd%ug_domain, is_restart=.true.)
  if (.not. restart_exists) return

  !Get the size of the tile dimension from the file.
  if (.not. field_exists(restart, "tile")) then
      call error_mesg("open_land_restart", "dimension 'tile' not found in file '" &
                      //trim(filename)//"'.", FATAL)
  endif
  ndims = get_variable_num_dimensions(restart%rhandle, "tile")
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, "tile", flen)
  restart%tile_dim_length = flen(1)
  deallocate(flen)

  !Get the size of the tile index dimension from the file.
  if (.not. field_exists(restart, "tile_index")) then
      call error_mesg("open_land_restart", "'tile_index' not found in file '" &
                      //trim(filename)//"'.", FATAL)
  endif
  ndims = get_variable_num_dimensions(restart%rhandle, "tile_index")
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, "tile_index", flen)
  allocate(restart%tidx(flen(1)))
  deallocate(flen)

  !Read in the tile_index field from the file.
  call read_data(restart%rhandle, "tile_index", restart%tidx)

  !Get the size of the cohort_index dimension from the file.
  if (field_exists(restart, "cohort_index")) then
      ndims = get_variable_num_dimensions(restart%rhandle, "cohort_index")
      allocate(flen(ndims))
      call get_variable_size(restart%rhandle, "cohort_index", flen)

      !Read in the cohort_index field from the file.
      allocate(restart%cidx(flen(1)))
      deallocate(flen)
      call read_data(restart%rhandle, "cohort_index", restart%cidx)
  endif
  ! TODO: possibly make tile index and cohort index names parameters in this module
  !       just constants, no sense to make them namelists vars
end subroutine open_land_restart

! ==============================================================================
subroutine save_land_restart(restart)
  type(land_restart_type), intent(inout) :: restart

  if (restart%should_free_rhandle) then
       call write_restart(restart%rhandle)
  endif
end subroutine save_land_restart

! ==============================================================================
! restet land restart data type to its initial state
subroutine free_land_restart(restart)
  type(land_restart_type), intent(inout) :: restart

  if (restart%should_free_rhandle) call close_file(restart%rhandle)
  restart%should_free_rhandle = .FALSE.
  restart%basename = ''
  restart%filename = ''
  if (allocated(restart%tidx)) deallocate(restart%tidx)
  if (allocated(restart%cidx)) deallocate(restart%cidx)
  restart%tile_dim_length = -1
end subroutine free_land_restart

! ==============================================================================
subroutine add_restart_axis(restart,name,data,is_unstructured,cartesian,units,longname,sense)
  type(land_restart_type),    intent(inout) :: restart
  character(len=*),           intent(in)    :: name
  real,                       intent(in)    :: data(:)
  logical,                    intent(in)    :: is_unstructured
  character(len=1), optional, intent(in)    :: cartesian
  character(len=*), optional, intent(in)    :: units, longname
  integer,          optional, intent(in)    :: sense

  integer :: n
  real, pointer :: data_(:)

  allocate(data_(size(data)))
  data_(:) = data(:)
  if (is_unstructured) then
    call register_axis(restart%rhandle, name)
  else
    call register_axis(restart%rhandle, name, size(data))
  endif
  call register_field(restart%rhandle, name, "double", (/name/))
  if (present(cartesian)) then
      call register_variable_attribute(restart%rhandle, name, "cartesian_axis", cartesian)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, name, "units", units)
  endif
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, name, "long_name", longname)
  endif
  if (present(sense)) then
      if (sense .eq. -1) then
          call register_variable_attribute(restart%rhandle, name, "positive", "down")
      else
          call register_variable_attribute(restart%rhandle, name, "positive", "up")
      endif
  endif
  call write_data(restart%rhandle,name,data_)
  deallocate(data_)

  ! record dimension information for future use
  n = restart%nax+1; restart%nax = n
  restart%ax(n)%name = name
  restart%ax(n)%len  = size(data)
end subroutine add_restart_axis

! ==============================================================================
logical function field_exists(restart,name)
  type(land_restart_type), intent(in) :: restart
  character(len=*),        intent(in) :: name

  field_exists = variable_exists(restart%rhandle, name)
end function field_exists

! ==============================================================================
subroutine add_scalar_data(restart,varname,datum,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  integer,          intent(in) :: datum
  character(len=*), intent(in), optional :: units, longname

  call register_field(restart%rhandle, varname, "int")
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_INT)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, datum)
end subroutine add_scalar_data

subroutine add_text_data(restart,varname,dim1,dim2,datum,longname)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: dim1, dim2 ! name of the text dimensions
  character,        intent(in) :: datum(:,:)
  character(len=*), intent(in), optional :: longname

  integer :: id_restart, ierr
  character(NF90_MAX_NAME)::dimnames(2)

  call error_mesg('add_text_data','does not work with new io yet', FATAL)
end subroutine add_text_data

subroutine add_tile_data_i0d_fptr_i0(restart,varname,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_i0)           :: fptr ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  integer, pointer :: data(:)

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)
  allocate(data(size(restart%tidx)))
  call gather_tile_data_i0d(fptr,restart%tidx,data)
  call register_field(restart%rhandle, varname, "int", (/"tile_index"/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_INT)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_i0d_fptr_i0

subroutine add_tile_data_r0d_fptr_r0(restart,varname,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0)           :: fptr ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  real, pointer :: data(:)

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)
  allocate(data(size(restart%tidx)))
  call gather_tile_data_r0d(fptr,restart%tidx,data)
  call register_field(restart%rhandle, varname, "double", (/"tile_index"/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_r0d_fptr_r0

subroutine add_tile_data_r0d_fptr_r0i(restart,varname,fptr,index,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0i)          :: fptr ! subroutine returning pointer to the data
  integer ,         intent(in) :: index ! index of the fptr array element to write
  character(len=*), intent(in), optional :: units, longname

  real, pointer :: data(:)

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)
  allocate(data(size(restart%tidx)))
  call gather_tile_data_r0i(fptr,index,restart%tidx,data)
  call register_field(restart%rhandle, varname, "double", (/"tile_index"/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_r0d_fptr_r0i

subroutine add_tile_data_r0d_fptr_r0ij(restart,varname,fptr,idx1,idx2,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0ij)         :: fptr ! subroutine returning pointer to the data
  integer ,         intent(in) :: idx1,idx2 ! indices of the fptr array element to write
  character(len=*), intent(in), optional :: units, longname

  real, pointer :: data(:)

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0ij', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)
  allocate(data(size(restart%tidx)))
  call gather_tile_data_r0ij(fptr,idx1,idx2,restart%tidx,data)
  call register_field(restart%rhandle, varname, "double", (/"tile_index"/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_r0d_fptr_r0ij

! given restart and name of the dimension, returns the size of this dimension
integer function dimlen(restart,dimname)
  type(land_restart_type), intent(in) :: restart ! restart data structure
  character(*),            intent(in) :: dimname ! name of the dimansion

  integer :: i
  dimlen = -1
  do i = 1,restart%nax
     if (trim(restart%ax(i)%name) == trim(dimname)) then
        dimlen=restart%ax(i)%len
        return
     endif
  enddo
  if (dimlen<1) call error_mesg('dimlen', 'axis "'//trim(dimname)//'" not found', FATAL)
end

subroutine add_tile_data_i1d_fptr_i0i(restart,varname,zdim,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim      ! name of the z-dimension
  procedure(fptr_i0i)          :: fptr    ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  integer, pointer :: data(:,:) ! needs to be pointer; we are passing ownership to restart object
  integer :: i,nlev

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_i1d_fptr_i0i', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  nlev = dimlen(restart,zdim)
  allocate(data(size(restart%tidx),nlev))
  call gather_tile_data_i1d(fptr,restart%tidx,data)
  call register_field(restart%rhandle, varname, "int", (/"tile_index",zdim/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_INT)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_i1d_fptr_i0i

subroutine add_tile_data_r1d_fptr_r0i(restart,varname,zdim,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim      ! name of the z-dimension
  procedure(fptr_r0i)          :: fptr    ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  real, pointer :: data(:,:) ! needs to be pointer; we are passing ownership to restart object
  integer :: i,nlev

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0i', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  nlev = dimlen(restart,zdim)
  allocate(data(size(restart%tidx),nlev))
  call gather_tile_data_r1d(fptr,restart%tidx,data)
  call register_field(restart%rhandle, varname, "double", (/"tile_index",zdim/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_r1d_fptr_r0i

subroutine add_tile_data_r1d_fptr_r0ij(restart,varname,zdim,fptr,index,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim      ! name of the z-dimension
  procedure(fptr_r0ij)         :: fptr    ! subroutine returning pointer to the data
  integer         , intent(in) :: index   ! index of the array element to write
  character(len=*), intent(in), optional :: units, longname

  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real, pointer :: data(:,:) ! needs to be pointer; we are passing ownership to restart object
  real, pointer :: ptr ! pointer to the tile data
  integer :: i,n,nlev

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0i', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  nlev = dimlen(restart,zdim)
  allocate(data(size(restart%tidx),nlev))
  data = NF90_FILL_DOUBLE

  ! gather data into an array along the tile dimension. It is assumed that
  ! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(restart%tidx)
     call get_tile_by_idx(restart%tidx(i), tileptr)
     do n = 1,nlev
        call fptr(tileptr,n,index, ptr)
        if(associated(ptr)) then
           data(i,n) = ptr
        endif
     enddo
  enddo
  call register_field(restart%rhandle, varname, "double", (/"tile_index",zdim/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_r1d_fptr_r0ij

subroutine add_tile_data_r1d_fptr_r0ijk(restart,varname,zdim,fptr,idx1,idx2,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim    ! name of the z-dimension
  procedure(fptr_r0ijk)        :: fptr    ! subroutine returning pointer to the data
  integer         , intent(in) :: idx1,idx2  ! indices of the array element to write
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real, pointer :: data(:,:) ! needs to be pointer; we are passing ownership to restart object
  real, pointer :: ptr ! pointer to the tile data
  integer :: i,n,nlev

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0i', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  nlev = -1
  do i = 1,restart%nax
     if (restart%ax(i)%name == zdim) nlev=restart%ax(i)%len
  enddo
  if (nlev<1) call error_mesg('add_tile_data_r0d_fptr_r0i', 'axis "'//trim(zdim)//'" not found', FATAL)

  allocate(data(size(restart%tidx),nlev))
  data = NF90_FILL_DOUBLE

  ! gather data into an array along the tile dimension. It is assumed that
  ! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(restart%tidx)
     call get_tile_by_idx(restart%tidx(i), tileptr)
     do n = 1,nlev
        call fptr(tileptr,n,idx1,idx2, ptr)
        if(associated(ptr)) then
           data(i,n) = ptr
        endif
     enddo
  enddo
  call register_field(restart%rhandle, varname, "double", (/"tile_index",zdim/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_r1d_fptr_r0ijk

subroutine add_tile_data_r2d_fptr_r0ij(restart,varname,dim1,dim2,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: dim1,dim2 ! names of extra dimensions
  procedure(fptr_r0ij)         :: fptr    ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  real, pointer :: data(:,:,:) ! needs to be pointer; we are passing ownership to restart object
  integer :: i,dim1len,dim2len

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r2d_fptr_r0ij', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  dim1len = dimlen(restart,dim1)
  dim2len = dimlen(restart,dim2)
  allocate(data(size(restart%tidx),dim1len,dim2len))
  call gather_tile_data_r2d(fptr,restart%tidx,data)
  call register_field(restart%rhandle, varname, "double", (/"tile_index",dim1,dim2/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_r2d_fptr_r0ij

subroutine add_tile_data_r2d_fptr_r0ijk(restart,varname,dim1,dim2,fptr,index,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: dim1,dim2 ! names of extra dimensions
  procedure(fptr_r0ijk)        :: fptr    ! subroutine returning pointer to the data
  integer         , intent(in) :: index   ! index of the array element to write
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  real, pointer :: data(:,:,:) ! needs to be pointer; we are passing ownership to restart object
  integer :: i,dim1len,dim2len

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r2d_fptr_r0ijk', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  dim1len = dimlen(restart,dim1)
  dim2len = dimlen(restart,dim2)
  allocate(data(size(restart%tidx),dim1len,dim2len))
  call gather_tile_data_r2d_idx(fptr,index,restart%tidx,data)
  call register_field(restart%rhandle, varname, "double", (/"tile_index",dim1,dim2/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
  if (present(longname)) then
      call register_variable_attribute(restart%rhandle, varname, "long_name", &
                                       longname)
  endif
  if (present(units)) then
      call register_variable_attribute(restart%rhandle, varname, "units", &
                                       units)
  endif
  call write_data(restart%rhandle, varname, data)
  deallocate(data)
end subroutine add_tile_data_r2d_fptr_r0ijk

! =============================================================================
subroutine get_scalar_data(restart,varname,datum)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  integer,          intent(out) :: datum

  call read_data(restart%rhandle,varname,datum)
end subroutine get_scalar_data

subroutine get_text_data(restart,varname,text)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character, allocatable, intent(out) :: text(:,:)

  call error_mesg('get_text_data','does not work with new io yet', FATAL)
end subroutine get_text_data

subroutine get_tile_data_i0d_fptr_i0(restart,varname,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_i0)           :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  integer, allocatable :: r(:) ! input data buffer

  allocate(r(size(restart%tidx)))
  call read_data(restart%rhandle, varname, r)
  call distrib_tile_data_i0d(fptr,restart%tidx,r)
  deallocate(r)
end subroutine get_tile_data_i0d_fptr_i0

subroutine get_tile_data_r0d_fptr_r0(restart,varname,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0)           :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  real, allocatable :: r(:) ! input data buffer

  allocate(r(size(restart%tidx)))
  call read_data(restart%rhandle, varname, r)
  call distrib_tile_data_r0d(fptr,restart%tidx,r)
  deallocate(r)
end subroutine get_tile_data_r0d_fptr_r0

subroutine get_tile_data_r0d_fptr_r0i(restart,varname,fptr,index)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0i)          :: fptr    ! subroutine returning pointer to the data
  integer,          intent(in) :: index   ! index where to read the data

  ! ---- local vars
  real, allocatable :: r(:) ! input data buffer

  allocate(r(size(restart%tidx)))
  call read_data(restart%rhandle, varname, r)
  call distrib_tile_data_r0d_idx(fptr,index,restart%tidx,r)
  deallocate(r)
end subroutine get_tile_data_r0d_fptr_r0i

subroutine get_tile_data_r0d_fptr_r0ij(restart,varname,fptr,i1, i2)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0ij)         :: fptr    ! subroutine returning pointer to the data
  integer,          intent(in) :: i1,i2   ! index where to read the data

  ! ---- local vars
  real, allocatable :: r(:) ! input data buffer

  allocate(r(size(restart%tidx)))
  call read_data(restart%rhandle, varname, r)
  call distrib_tile_data_r0d_ij(fptr,i1,i2,restart%tidx,r)
  deallocate(r)
end subroutine get_tile_data_r0d_fptr_r0ij

subroutine get_tile_data_r1d_fptr_r0i(restart,varname,zdim,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  procedure(fptr_r0i)          :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  integer,dimension(:),allocatable :: flen ! size of the input field
  real, allocatable :: r(:,:) ! input data buffer
  integer :: ndims

  if (.not. field_exists(restart, zdim)) then
      call error_mesg("get_tile_data_r0d_fptr_r0i", &
            "axis '"//trim(zdim)//"' was not found in file '"//trim(restart%basename)//"'.", &
            FATAL)
  endif

  !Get the size of z-dimension from the file.
  ndims = get_variable_num_dimensions(restart%rhandle, zdim)
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, zdim, flen)

  !Read in the field from the file.
  allocate(r(size(restart%tidx),flen(1)))
  call read_data(restart%rhandle, varname, r)
  call distrib_tile_data_r1d(fptr,restart%tidx,r)
  deallocate(r)
  deallocate(flen)
end subroutine get_tile_data_r1d_fptr_r0i

subroutine get_tile_data_i1d_fptr_i0i(restart,varname,zdim,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  procedure(fptr_i0i)          :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  integer,dimension(:),allocatable :: flen ! size of the input field
  integer, allocatable :: r(:,:) ! input data buffer
  integer :: ndims

  if (.not. field_exists(restart, zdim)) then
      call error_mesg("get_tile_data_i1d_fptr_i0i", &
            "axis '"//trim(zdim)//"' was not found in file '"//trim(restart%basename)//"'.", &
            FATAL)
  endif

  !Get the size of z-dimension from the file.
  ndims = get_variable_num_dimensions(restart%rhandle, zdim)
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, zdim, flen)

  !Read in the field from the file.
  allocate(r(size(restart%tidx),flen(1)))
  call read_data(restart%rhandle, varname, r)
  call distrib_tile_data_i1d(fptr,restart%tidx,r)
  deallocate(r)
  deallocate(flen)
end subroutine get_tile_data_i1d_fptr_i0i

subroutine get_tile_data_r1d_fptr_r0ij(restart,varname,zdim,fptr,index)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  procedure(fptr_r0ij)         :: fptr    ! subroutine returning pointer to the data
  integer ,         intent(in) :: index

  ! ---- local vars
  integer,dimension(:),allocatable :: flen ! size of the input field
  real, allocatable :: r(:,:) ! input data buffer
  integer :: ndims

  if (.not. field_exists(restart, zdim)) then
      call error_mesg("get_tile_data_r1d_fptr_r0ij", &
            "axis '"//trim(zdim)//"' was not found in file '"//trim(restart%basename)//"'.", &
            FATAL)
  endif

  !Get the size of z-dimension from the file.
  ndims = get_variable_num_dimensions(restart%rhandle, zdim)
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, zdim, flen)

  !Read in the field from the file.
  allocate(r(size(restart%tidx),flen(1)))
  call read_data(restart%rhandle, varname, r)
  call distrib_tile_data_r1d_idx(fptr,index,restart%tidx,r)
  deallocate(r)
  deallocate(flen)
end subroutine get_tile_data_r1d_fptr_r0ij

subroutine get_tile_data_r1d_fptr_r0ijk(restart,varname,zdim,fptr,idx1,idx2)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  procedure(fptr_r0ijk)        :: fptr    ! subroutine returning pointer to the data
  integer ,         intent(in) :: idx1,idx2

  ! ---- local vars
  integer,dimension(:),allocatable :: flen ! size of the input field
  real, allocatable :: r(:,:) ! input data buffer
  integer :: ndims

  if (.not. field_exists(restart, zdim)) then
      call error_mesg("get_tile_data_r1d_fptr_r0ijk", &
            "axis '"//trim(zdim)//"' was not found in file '"//trim(restart%basename)//"'.", &
            FATAL)
  endif

  !Get the size of z-dimension from the file.
  ndims = get_variable_num_dimensions(restart%rhandle, zdim)
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, zdim, flen)

  !Read in the field from the file.
  allocate(r(size(restart%tidx),flen(1)))
  call read_data(restart%rhandle, varname, r)
! call distrib_tile_data_r1d_idx(fptr,idx1,idx2,restart%tidx,r)
  deallocate(r)
  deallocate(flen)
end subroutine get_tile_data_r1d_fptr_r0ijk

subroutine get_tile_data_r2d_fptr_r0ij(restart,varname,dim1,dim2,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: dim1,dim2 ! names of the dimensions
  procedure(fptr_r0ij)         :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  integer,dimension(:),allocatable :: flen ! size of the input field
  integer :: n,m
  real, allocatable :: r(:,:,:) ! input data buffer
  integer :: ndims

  if (.not. field_exists(restart, dim1)) then
      call error_mesg("get_tile_data_r2d_fptr_r0ij", &
           "axis '"//trim(dim1)//"' was not found in file '"//trim(restart%basename)//"'.", &
           FATAL)
  endif

  !Get the size of the first dimension of the field.
  ndims = get_variable_num_dimensions(restart%rhandle, dim1)
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, dim1, flen)
  n = flen(1)
  deallocate(flen)

  if (.not. field_exists(restart, dim2)) then
      call error_mesg("get_tile_data_r2d_fptr_r0ij", &
           "axis '"//trim(dim2)//"' was not found in file '"//trim(restart%basename)//"'.", &
           FATAL)
  endif

  !Get the size of the second dimension of the field.
  ndims = get_variable_num_dimensions(restart%rhandle, dim2)
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, dim2, flen)
  m = flen(1)
  deallocate(flen)

  !Read in the field data from the file.
  allocate(r(size(restart%tidx),n,m))
  call read_data(restart%rhandle, varname, r)
  call distrib_tile_data_r2d(fptr,restart%tidx,r)
  deallocate(r)
end subroutine get_tile_data_r2d_fptr_r0ij

subroutine get_tile_data_r2d_fptr_r0ijk(restart,varname,dim1,dim2,fptr,index)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: dim1,dim2 ! names of the extra dimensions
  procedure(fptr_r0ijk)        :: fptr    ! subroutine returning pointer to the data
  integer,          intent(in) :: index   ! index where to read the data

  ! ---- local vars
  integer,dimension(:),allocatable :: flen ! size of the input field
  integer ::  m,n
  real, allocatable :: r(:,:,:) ! input data buffer
  integer :: ndims

  if (.not. field_exists(restart, dim1)) then
      call error_mesg("get_tile_data_r2d_fptr_r0ijk", &
           "axis '"//trim(dim1)//"' was not found in file '"//trim(restart%basename)//"'.", &
           FATAL)
  endif

  !Get the size of the first dimension of the field.
  ndims = get_variable_num_dimensions(restart%rhandle, dim1)
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, dim1, flen)
  n = flen(1)
  deallocate(flen)

  if (.not. field_exists(restart, dim2)) then
      call error_mesg("get_tile_data_r2d_fptr_r0ijk", &
           "axis '"//trim(dim2)//"' was not found in file '"//trim(restart%basename)//"'.", &
           FATAL)
  endif

  !Get the size of the second dimension of the field.
  ndims = get_variable_num_dimensions(restart%rhandle, dim2)
  allocate(flen(ndims))
  call get_variable_size(restart%rhandle, dim2, flen)
  m = flen(1)
  deallocate(flen)

  !Read in the field from the file.
  allocate(r(size(restart%tidx),n,m))
  call read_data(restart%rhandle, varname, r)
  call distrib_tile_data_r2d_idx(fptr,index,restart%tidx,r)
  deallocate(r)
end subroutine get_tile_data_r2d_fptr_r0ijk

subroutine create_tile_out_file_idx_new(rhandle,name,tidx,tile_dim_length,zaxis_data,soilCCohort_data)
  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: rhandle     ! restart file handle
  character(len=*),      intent(in)  :: name                ! name of the file to create
  integer              , intent(in)  :: tidx(:)             ! integer compressed index of tiles (local)
  integer              , intent(in)  :: tile_dim_length     ! length of tile axis
  real,        optional, intent(in)  :: zaxis_data(:)       ! data for the Z-axis
  real,        optional, intent(in)  :: soilCCohort_data(:)

  logical :: s
  integer :: ntidx
  integer, dimension(:), allocatable :: npes_tidx !Tile index length of each pe in file's pelist.
  integer, dimension(:), allocatable :: npes_tidx_start !Offset of tile index of each pe in file's pelist.
  integer, dimension(tile_dim_length) :: buffer
  integer :: i

  s = open_file(rhandle, name, "overwrite", lnd%ug_domain, is_restart=.true.)
  call register_axis(rhandle, "lon", size(lnd%coord_glon))
  call register_field(rhandle, "lon", "double", (/"lon"/))
  call register_variable_attribute(rhandle, "lon", "units", "degrees_east")
  call register_variable_attribute(rhandle, "lon", "long_name", "longitude")
  call register_variable_attribute(rhandle, "lon", "cartesian_axis", "X")
  call write_data(rhandle, "lon", lnd%coord_glon)

  call register_axis(rhandle, "lat", size(lnd%coord_glat))
  call register_field(rhandle, "lat", "double", (/"lat"/))
  call register_variable_attribute(rhandle, "lat", "units", "degrees_north")
  call register_variable_attribute(rhandle, "lat", "long_name", "latitude")
  call register_variable_attribute(rhandle, "lat", "cartesian_axis", "Y")
  call write_data(rhandle, "lat", lnd%coord_glat)

  ! the size of tile dimension really does not matter for the output, but it does
  ! matter for uncompressing utility, since it uses it as a size of the array to
  ! unpack to create tile index dimension and variable.
  call register_axis(rhandle, "tile", tile_dim_length)
  call register_field(rhandle, "tile", "int", (/"tile"/))
  call register_variable_attribute(rhandle, "tile", "long_name", "tile number within grid cell")
  do i = 1, tile_dim_length
      buffer(i) = i
  enddo
  call write_data(rhandle, "tile", buffer)

  ntidx = size(tidx)
  call compressed_start_and_count(rhandle, ntidx, npes_tidx_start, npes_tidx)
  call register_axis(rhandle, tile_index_name, npes_corner=npes_tidx_start, npes_nelems=npes_tidx)
  deallocate(npes_tidx)
  deallocate(npes_tidx_start)
  call register_field(rhandle, tile_index_name, "int", (/tile_index_name/))
  call register_variable_attribute(rhandle, tile_index_name, "long_name", "compressed land point index")
  call register_variable_attribute(rhandle, tile_index_name, "compress", "tile lat lon")
  call register_variable_attribute(rhandle, tile_index_name, "units", "")
  call register_variable_attribute(rhandle, tile_index_name, "valid_min", 0)
  call write_data(rhandle, tile_index_name, tidx)

  if (present(zaxis_data)) then
      call register_axis(rhandle, "zfull", size(zaxis_data))
      call register_field(rhandle, "zfull", "double", (/"zfull"/))
      call register_variable_attribute(rhandle, "zfull", "long_name", "full level")
      call register_variable_attribute(rhandle, "zfull", "units", "m")
      call register_variable_attribute(rhandle, "zfull", "positive", "down")
      call write_data(rhandle,"zfull",zaxis_data)
  endif

  if (present(soilCCohort_data)) then
      call register_axis(rhandle, "soilCCohort", size(soilCCohort_data))
      call register_field(rhandle, "soilCCohort", "double", (/"soilCCohort"/))
      call register_variable_attribute(rhandle, "soilCCohort", "long_name", "Soil carbon cohort")
      call write_data(rhandle,"soilCCohort",soilCCohort_data)
  endif
end subroutine create_tile_out_file_idx_new

! ============================================================================
! given a tile existence detection function, allocates and fills the tile index vector
subroutine gather_tile_index(tile_exists,idx)
  procedure(tile_exists_func) :: tile_exists   ! existence detector function:
         ! returns true if specific tile exists (hence should be written to restart)
  integer, allocatable,  intent(out) :: idx(:) ! rank local tile index vector

  ! ---- local vars
  type(land_tile_enum_type) :: ce ! tile list enumerator
  type(land_tile_type), pointer :: tile
  integer :: i,j,k,n

  ! count total number of tiles in this domain
  ce = first_elmt(land_tile_map)
  n  = 0
  do while (loop_over_tiles(ce,tile))
     if (tile_exists(tile)) n = n+1
  end do

  ! calculate compressed tile index to be written to the restart file;
  allocate(idx(max(n,1))); idx(:) = -1 ! set init value to a known invalid index
  ce = first_elmt(land_tile_map, lnd%ls)
  n = 1
  do while (loop_over_tiles(ce,tile,i=i,j=j,k=k))
     if(tile_exists(tile)) then
        idx(n) = (k-1)*lnd%nlon*lnd%nlat + (j-1)*lnd%nlon + (i-1)
        n = n+1
     endif
  end do
end subroutine gather_tile_index

! ============================================================================
! given compressed index, returns a pointer to the tile corresponding to this 
! index, or NULL if the index is outside current domain, or if such tile does 
! not exist.
subroutine get_tile_by_idx(idx,ptr)
   integer, intent(in)           :: idx ! compressed-by-gathering tile index
   type(land_tile_type), pointer :: ptr

   ! ---- local vars
   integer :: g,k,npts,l
   type(land_tile_enum_type) :: ce

   ptr=>null()

   if(idx<0) return ! negative indices do not correspond to any tile

   ! given tile idx, calculate global lon, lat, and tile indices
   k = idx
   npts = lnd%nlon*lnd%nlat
   g = modulo(k,npts)+1 ; k = k/npts
   ! do nothing if the index is outside of our domain
   if (g<lnd%gs.or.g>lnd%ge) return ! skip points outside of domain
   ! loop through the list of tiles at the given point to find k+1st tile
   l = lnd%l_index(g)
   if(l < lnd%ls .OR. l > lnd%le) then
      print*, " l= ", l, g, lnd%ls, lnd%le, mpp_pe()
      call error_mesg("land_tile_io", "l < lnd%ls .OR. l > lnd%le", FATAL)
   endif
   ce = first_elmt(land_tile_map(l))
   do while(loop_over_tiles(ce, ptr))
      k = k-1
      if (k<0) exit ! from loop
   enddo
   ! NOTE that at the end of the loop (that is, if there are less tiles in the list
   ! then requested by the idx), loop_over_tiles(ce,ptr) returns NULL

end subroutine get_tile_by_idx

! ============================================================================
subroutine gather_tile_data_i0d(fptr,idx,data)
  procedure(fptr_i0)  :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  integer, intent(out) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  integer, pointer :: ptr ! pointer to the tile data
  integer :: i

  data = NF90_FILL_INT

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) data(i)=ptr
  enddo
end subroutine gather_tile_data_i0d

subroutine gather_tile_data_r0d(fptr,idx,data)
  procedure(fptr_r0) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

  data = NF90_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) data(i)=ptr
  enddo
end subroutine gather_tile_data_r0d

subroutine gather_tile_data_r0i(fptr,n,idx,data)
  procedure(fptr_r0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: n   ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

  data = NF90_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     call fptr(tileptr, n, ptr)
     if(associated(ptr)) data(i)=ptr
  enddo
end subroutine gather_tile_data_r0i

subroutine gather_tile_data_r0ij(fptr,n,m,idx,data)
  procedure(fptr_r0ij):: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: n, m   ! additional index arguments for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

  data = NF90_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     call fptr(tileptr, n, m, ptr)
     if(associated(ptr)) data(i)=ptr
  enddo
end subroutine gather_tile_data_r0ij

subroutine gather_tile_data_r1d(fptr,idx,data)
  procedure(fptr_r0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,j

  data = NF90_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, ptr)
        if(associated(ptr)) data(i,j)=ptr
     enddo
  enddo
end subroutine gather_tile_data_r1d

subroutine gather_tile_data_i1d(fptr,idx,data)
  procedure(fptr_i0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  integer, intent(out) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  integer, pointer :: ptr ! pointer to the tile data
  integer :: i,j

  data = NF90_FILL_INT

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, ptr)
        if(associated(ptr)) data(i,j)=ptr
     enddo
  enddo
end subroutine gather_tile_data_i1d

subroutine gather_tile_data_r2d(fptr,idx,data)
  procedure(fptr_r0ij):: fptr ! subroutine returning the pointer to the data to be written
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:,:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,k,m

  data = NF90_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     do k=1,size(data,2)
     do m=1,size(data,3)
        call fptr(tileptr, k, m, ptr)
        if(associated(ptr)) data(i,k,m)=ptr
     enddo
     enddo
  enddo
end subroutine gather_tile_data_r2d

subroutine gather_tile_data_r2d_idx(fptr,n,idx,data)
  procedure(fptr_r0ijk) :: fptr ! subroutine returning the pointer to the data
  integer, intent(in) :: n ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:,:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,k,m

  data = NF90_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     do k=1,size(data,2)
     do m=1,size(data,3)
        call fptr(tileptr, k, m, n, ptr)
        if(associated(ptr)) data(i,k,m)=ptr
     enddo
     enddo
  enddo
end subroutine gather_tile_data_r2d_idx

subroutine distrib_tile_data_i0d(fptr,idx,data)
  procedure(fptr_i0) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  integer, intent(in) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  integer, pointer :: ptr ! pointer to the tile data
  integer :: i

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) ptr=data(i)
  enddo
end subroutine distrib_tile_data_i0d

subroutine distrib_tile_data_r0d(fptr,idx,data)
  procedure(fptr_r0) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) ptr=data(i)
  enddo
end subroutine distrib_tile_data_r0d

subroutine distrib_tile_data_r0d_idx(fptr,n,idx,data)
  procedure(fptr_r0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: n ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     call fptr(tileptr, n, ptr)
     if(associated(ptr)) ptr=data(i)
  enddo
end subroutine distrib_tile_data_r0d_idx

subroutine distrib_tile_data_r0d_ij(fptr,n,m,idx,data)
  procedure(fptr_r0ij):: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: n,m ! additional index arguments for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     call fptr(tileptr, n, m, ptr)
     if(associated(ptr)) ptr=data(i)
  enddo
end subroutine distrib_tile_data_r0d_ij

subroutine distrib_tile_data_i1d(fptr,idx,data)
  procedure(fptr_i0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  integer, intent(in) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  integer, pointer :: ptr ! pointer to the tile data
  integer :: i,j

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, ptr)
        if(associated(ptr)) ptr=data(i,j)
     enddo
  enddo
end subroutine distrib_tile_data_i1d

subroutine distrib_tile_data_r1d(fptr,idx,data)
  procedure(fptr_r0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,j

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, ptr)
        if(associated(ptr)) ptr=data(i,j)
     enddo
  enddo
end subroutine distrib_tile_data_r1d

subroutine distrib_tile_data_r1d_idx(fptr,n,idx,data)
  procedure(fptr_r0ij) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: n ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,j

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, n, ptr)
        if(associated(ptr)) ptr=data(i,j)
     enddo
  enddo
end subroutine distrib_tile_data_r1d_idx

subroutine distrib_tile_data_r2d(fptr,idx,data)
  procedure(fptr_r0ij):: fptr ! subroutine returning the pointer to the data to be written
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:,:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real, pointer :: ptr ! pointer to the tile data
  integer :: i,k,m

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     do k=1,size(data,2)
     do m=1,size(data,3)
        call fptr(tileptr, k, m, ptr)
        if(associated(ptr)) ptr=data(i,k,m)
     enddo
     enddo
  enddo
end subroutine distrib_tile_data_r2d

subroutine distrib_tile_data_r2d_idx(fptr,n,idx,data)
  procedure(fptr_r0ijk) :: fptr ! subroutine returning the pointer to the data
  integer, intent(in) :: n ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:,:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real, pointer :: ptr ! pointer to the tile data
  integer :: i,k,m

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i), tileptr)
     do k=1,size(data,2)
     do m=1,size(data,3)
        call fptr(tileptr, k, m, n, ptr)
        if(associated(ptr)) ptr=data(i,k,m)
     enddo
     enddo
  enddo
end subroutine distrib_tile_data_r2d_idx

! ============================================================================
subroutine override_tile_data_r0d_fptr(fieldname,fptr,time,override)
  character(len=*), intent(in)   :: fieldname ! field to override
  procedure(fptr_r0) :: fptr ! subroutine returning pointer to the data
  type(time_type),  intent(in)   :: time      ! model time
  logical, optional, intent(out) :: override  ! true if the field has been
                                              ! overridden successfully

  ! ---- local vars
  real    :: data1D(lnd%ls:lnd%le)  ! storage for the input data
  logical :: override_

  call data_override_ug('LND',fieldname,data1D, time, override_ )
  if(present(override)) override=override_
  if(.not.override_) return ! do nothing if the field was not overridden

  ! distribute the data over the tiles
  call put_to_tiles_r0d_fptr(data1d,land_tile_map,fptr)

end subroutine override_tile_data_r0d_fptr

end module
