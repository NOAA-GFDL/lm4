module cohort_io_mod

use fms_mod,          only : error_mesg, FATAL, WARNING, get_mosaic_tile_file
use fms_io_mod,       only : restart_file_type, get_instance_filename
use mpp_mod,          only : mpp_pe, mpp_max, mpp_send, mpp_recv, mpp_sync, &
                             COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4, &
                             mpp_sync_self, stdout
use nf_utils_mod,     only : nfu_inq_dim, nfu_get_var, nfu_put_var, &
     nfu_get_rec, nfu_put_rec, nfu_def_dim, nfu_def_var, nfu_put_att, &
     nfu_inq_var
use land_io_mod,      only : print_netcdf_error, input_buf_size
use land_tile_mod,    only : land_tile_map, land_tile_type, land_tile_list_type, &
     land_tile_enum_type, first_elmt, tail_elmt, next_elmt, &
     current_tile, operator(/=), nitems, loop_over_tiles

use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, get_tile_data, &
     get_tile_by_idx

use vegn_cohort_mod, only: vegn_cohort_type
use land_data_mod, only : lnd

use fms_io_mod, only: fms_io_unstructured_register_restart_axis
use fms_io_mod, only: fms_io_unstructured_register_restart_field
use fms_io_mod, only: HIDX
use fms_io_mod, only: fms_io_unstructured_read



use fms2_io_mod, only: FmsNetcdfUnstructuredDomainFile_t, &
                       register_axis, register_field, register_variable_attribute, &
                       read_data



implicit none
private

! ==== public interfaces =====================================================
public :: read_create_cohorts
public :: create_cohort_dimension
public :: add_cohort_data, add_int_cohort_data
public :: get_cohort_data, get_int_cohort_data
! remove when cleaning up:
public :: write_cohort_data_r0d, write_cohort_data_i0d
public :: gather_cohort_index, gather_cohort_data
public :: create_cohort_dimension_new
! ==== end of public interfaces ==============================================

interface gather_cohort_data
   module procedure gather_cohort_data_r0d
   module procedure gather_cohort_data_i0d
end interface gather_cohort_data

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'cohort_io_mod'
! name of the "compressed" dimension (and dimension variable) in the output
! netcdf files -- that is, the dimensions written out using compression by
! gathering, as described in CF conventions.
character(len=*),   parameter :: cohort_index_name   = 'cohort_index'

abstract interface
  ! given land cohort, returns pointer to some scalar real data
  ! within this cohort, or an unassociated pointer if there is no data
  subroutine cptr_r0(tile, ptr)
     import vegn_cohort_type
     type(vegn_cohort_type), pointer :: tile ! input
     real                , pointer :: ptr  ! returned pointer to the data
  end subroutine cptr_r0
  ! given land cohort, returns pointer to some scalar real data
  ! within this cohort, or an unassociated pointer if there is no data
  subroutine cptr_i0(tile, ptr)
     import vegn_cohort_type
     type(vegn_cohort_type), pointer :: tile ! input
     integer               , pointer :: ptr  ! returned pointer to the data
  end subroutine cptr_i0
end interface
! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
! given compressed index, sizes of the global grid, 2D array of tile lists
! and the lower boundaries of this array, returns a pointer to the cohort
! corresponding to the compressed index, or NULL is the index is outside
! current domain, or such tile does not exist, or such cohort does not exist.
subroutine get_cohort_by_idx(idx,ntiles,ptr)
   integer, intent(in) :: idx ! index
   integer, intent(in) :: ntiles ! size of tile dimension
   type(vegn_cohort_type), pointer :: ptr

   ! ---- local vars
   integer :: tile_idx, k
   type(land_tile_type), pointer :: tile

   ptr=>NULL()
   if ( idx < 0 ) return
   tile_idx = modulo(idx,lnd%nlon*lnd%nlat*ntiles)
   call get_tile_by_idx(tile_idx,tile)
   if(associated(tile)) then
      if (associated(tile%vegn)) then
         k = idx/(lnd%nlon*lnd%nlat*ntiles) ! calculate cohort index within a tile
         ptr=>tile%vegn%cohorts(k+1)
      endif
   endif
end subroutine get_cohort_by_idx

! ============================================================================
subroutine read_create_cohorts(restart)
  type(land_restart_type), intent(inout) :: restart

  if (.not.allocated(restart%cidx)) call error_mesg('read_create_cohorts', &
      'cohort index not found in file "'//restart%filename//'"',FATAL)
  call read_create_cohorts_new(restart%cidx,restart%tile_dim_length)
end subroutine

! ============================================================================
subroutine read_create_cohorts_new(idx,ntiles)
  integer, intent(in) :: idx(:)
  integer, intent(in) :: ntiles

  integer :: ncohorts ! total number of cohorts in restart file
  integer :: nlon, nlat ! size of respective dimensions

  integer :: i,j,t,k,m, n, npts, g, l
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  character(len=64) :: info ! for error message

  ! get the size of dimensions
  nlon = lnd%nlon
  nlat = lnd%nlat
  ncohorts = size(idx)
  npts = nlon*nlat

  do n = 1,ncohorts
     if(idx(n)<0) cycle ! skip illegal indices
     k = idx(n)
     g = modulo(k,npts)+1
     if(g<lnd%gs.or.g>lnd%ge) cycle ! skip points outside of domain
     l = lnd%l_index(g)
     k = k/npts
     t = modulo(k,ntiles)+1 ; k = k/ntiles
     k = k+1

     ce = first_elmt(land_tile_map(l))
     do m = 1,t-1
        ce=next_elmt(ce)
     enddo
     tile=>current_tile(ce)

     if (.not. associated(tile)) then
         call error_mesg("read_create_cohorts_new", &
                         "current tile returned null pointer", &
                         FATAL)
     endif

     if(.not.associated(tile%vegn)) then
        info = ''
        write(info,'("(",3i3,")")')i,j,t
        call error_mesg('read_create_cohort',&
             'vegn tile'//trim(info)//' does not exist, but is necessary to create a cohort', &
             WARNING)
     else
        tile%vegn%n_cohorts = tile%vegn%n_cohorts + 1
     endif
  enddo

  ! go through all tiles in the domain and allocate requested numner of cohorts
  ce = first_elmt(land_tile_map); te = tail_elmt(land_tile_map)
  do while (ce/=te)
     tile=>current_tile(ce); ce = next_elmt(ce)
     if(.not.associated(tile%vegn))cycle
     allocate(tile%vegn%cohorts(tile%vegn%n_cohorts))
  enddo
end subroutine read_create_cohorts_new

! ============================================================================
! creates cohort dimension, if necessary, in the output restart file. NOTE
subroutine create_cohort_dimension(restart)
  type(land_restart_type), intent(inout) :: restart

  call create_cohort_dimension_new(restart%rhandle,restart%cidx,restart%basename,restart%tile_dim_length)
end subroutine create_cohort_dimension

! ============================================================================
! creates cohort dimension, if necessary, in the output restart file. NOTE
! that this subroutine should be called even if restart has not been created
! (because, for example, there happen to be no vegetation in a certain domain),
! for the reason that it calls mpp_max, and that should be called for each
! processor to work.
subroutine create_cohort_dimension_new(rhandle,cidx,name,tile_dim_length)
  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: rhandle ! fms_io restart file data type
  integer, allocatable,    intent(out)   :: cidx(:) ! rank local tile index vector
  character(len=*),        intent(in)    :: name    ! name of the restart file
  integer,                 intent(in)    :: tile_dim_length ! length of tile axis

  integer :: max_cohorts

  call gather_cohort_index(tile_dim_length,cidx)
  max_cohorts = global_max_cohorts()

  call create_cohort_out_file_idx(rhandle,name,cidx,max(max_cohorts,1))
end subroutine create_cohort_dimension_new

subroutine create_cohort_out_file_idx(rhandle,name,cidx,cohorts_dim_length)
  type(FmsNetcdfUnstructuredDomainFile_t),intent(inout) :: rhandle ! fms_io restart file data type
  character(len=*),      intent(in)  :: name                ! name of the file to create
  integer              , intent(in)  :: cidx(:)             ! integer compressed index of tiles (local)
  integer              , intent(in)  :: cohorts_dim_length  ! length of cohorts axis

  ! ---- local vars
  character(256) :: file_name ! full name of the file, including the processor number

  ! form the full name of the file
  call get_instance_filename(trim(name), file_name)
! call get_mosaic_tile_file(trim(file_name),file_name,lnd%ug_domain)

  ! the size of tile dimension really does not matter for the output, but it does
  ! matter for uncompressing utility, since it uses it as a size of the array to
  ! unpack to create tile index dimension and variable.
  call register_axis(rhandle, cohort_index_name, .true.)
  call register_field(rhandle, cohort_index_name, "int", (/cohort_index_name/))
  call register_variable_attribute(rhandle, cohort_index_name, "compressed", &
                                   "cohort tile lat lon")
  call register_variable_attribute(rhandle, cohort_index_name, "units", &
                                   "none")
  call register_variable_attribute(rhandle, cohort_index_name, "long_name", &
                                   "compressed vegetation cohort index")
end subroutine create_cohort_out_file_idx

subroutine distrib_cohort_data_i0d(fptr,idx,ntiles,data)
  integer, intent(in) :: idx(:) ! local vector of cohort indices
  integer, intent(in) :: ntiles ! size of the tile dimension
  integer, intent(in) :: data(:) ! local cohort data
  procedure(cptr_i0) :: fptr ! subroutine returning pointer to the data

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cohort
  integer, pointer :: ptr ! pointer to the individual cohort data
  integer :: mask(size(data)) ! mask of valid data
  integer :: i

  ! gather data into an array along the cohort dimension
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), ntiles, cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) ptr = data(i)
     endif
  enddo
end subroutine distrib_cohort_data_i0d

subroutine distrib_cohort_data_r0d(fptr,idx,ntiles,data)
  integer, intent(in) :: idx(:) ! local vector of cohort indices
  integer, intent(in) :: ntiles ! size of the tile dimension
  real, intent(in) :: data(:) ! local cohort data
  procedure(cptr_r0) :: fptr ! subroutine returning pointer to the data

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cohort
  real, pointer :: ptr ! pointer to the individual cohort data
  integer :: mask(size(data))
  integer :: i

  ! gather data into an array along the cohort dimension
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), ntiles, cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) ptr = data(i)
     endif
  enddo
end subroutine distrib_cohort_data_r0d

! count max number of cohorts per tile
integer function global_max_cohorts()
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile

  ce = first_elmt(land_tile_map)
  global_max_cohorts = 0
  do while (loop_over_tiles(ce,tile))
     if(associated(tile%vegn)) &
        global_max_cohorts = max(global_max_cohorts,tile%vegn%n_cohorts)
  enddo
  call mpp_max(global_max_cohorts)
end function global_max_cohorts

subroutine gather_cohort_index(ntiles, cidx)
  integer,              intent(in)  :: ntiles
  integer, allocatable, intent(out) :: cidx(:)   ! integer compressed index of tiles

  integer :: i,j,k,c,n
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile

  ! count total number of cohorts in our compute domain
  ce = first_elmt(land_tile_map)
  n = 0
  do while (loop_over_tiles(ce,tile))
     if(associated(tile%vegn)) n = n+tile%vegn%n_cohorts
  enddo

  ! calculate compressed cohort index to be written to the restart file
  allocate(cidx(max(n,1))) ; cidx(:) = -1
  ce = first_elmt(land_tile_map, lnd%ls)
  n = 1
  do while (loop_over_tiles(ce,tile,i=i,j=j,k=k))
     if(associated(tile%vegn)) then
        do c = 1,tile%vegn%n_cohorts
           cidx (n) = &
                (c-1)*lnd%nlon*lnd%nlat*ntiles + &
                (k-1)*lnd%nlon*lnd%nlat + &
                (j-1)*lnd%nlon + &
                (i-1)
           n = n+1
        enddo
     endif
  end do
end subroutine gather_cohort_index

subroutine gather_cohort_data_i0d(fptr,idx,ntiles,data)
  procedure(cptr_i0) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:) ! local vector of cohort indices
  integer, intent(in) :: ntiles ! size of the tile dimension
  integer, intent(out) :: data(:) ! local cohort data

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cohort
  integer, pointer :: ptr ! pointer to the individual cohort data
  integer :: i

  ! gather data into an array along the cohort dimension
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), ntiles, cohort)
     data(i) = NF_FILL_INT
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) data(i) = ptr
     endif
  enddo
end subroutine gather_cohort_data_i0d

subroutine gather_cohort_data_r0d(fptr,idx,ntiles,data)
  procedure(cptr_r0)   :: fptr ! subroutine returning pointer to the data
  integer, intent(in)  :: idx(:) ! local vector of cohort indices
  integer, intent(in)  :: ntiles ! size of the tile dimension
  real,    intent(out) :: data(:) ! local cohort data

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cohort
  real, pointer :: ptr ! pointer to the individual cohort data
  integer :: i

  ! gather data into an array along the cohort dimension
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), ntiles, cohort)
     data(i) = NF_FILL_DOUBLE
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) data(i) = ptr
     endif
  enddo
end subroutine gather_cohort_data_r0d

! ============================================================================
subroutine add_cohort_data(restart,varname,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(cptr_r0)           :: fptr ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  real, pointer :: r(:)

  allocate(r(size(restart%cidx)))
  call gather_cohort_data_r0d(fptr,restart%cidx,restart%tile_dim_length,r)
  call register_axis(restart%rhandle, varname, .true.)
  call register_field(restart%rhandle, varname, "double", (/varname/))
  if (present(units)) then
    call register_variable_attribute(restart%rhandle, varname, "units", units)
  endif
  if (present(longname)) then
    call register_variable_attribute(restart%rhandle, varname, "long_name", longname)
  endif
end subroutine add_cohort_data

! ============================================================================
subroutine add_int_cohort_data(restart,varname,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(cptr_i0)           :: fptr ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  integer, pointer :: r(:)
  integer :: id_restart

  allocate(r(size(restart%cidx)))
  call gather_cohort_data_i0d(fptr,restart%cidx,restart%tile_dim_length,r)
  call register_axis(restart%rhandle, varname, .true.)
  call register_field(restart%rhandle, varname, "int", (/varname/))
  if (present(units)) then
    call register_variable_attribute(restart%rhandle, varname, "units", units)
  endif
  if (present(longname)) then
    call register_variable_attribute(restart%rhandle, varname, "long_name", longname)
  endif
end subroutine add_int_cohort_data

! ============================================================================
subroutine get_cohort_data(restart,varname,fptr)
  type(land_restart_type), intent(in) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(cptr_r0)           :: fptr ! subroutine returning pointer to the data

  real, allocatable :: r(:)

  if (.not.allocated(restart%cidx)) call error_mesg('read_create_cohorts', &
      'cohort index not found in file "'//restart%filename//'"',FATAL)
  allocate(r(size(restart%cidx)))
  call read_data(restart%rhandle, varname, r, unlim_dim_level=1)
  call distrib_cohort_data_r0d(fptr,restart%cidx,restart%tile_dim_length,r)
  deallocate(r)
end subroutine get_cohort_data

! ============================================================================
subroutine get_int_cohort_data(restart,varname,fptr)
  type(land_restart_type), intent(in) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(cptr_i0)           :: fptr ! subroutine returning pointer to the data

  integer, allocatable :: r(:)

  if (.not.allocated(restart%cidx)) call error_mesg('read_create_cohorts', &
      'cohort index not found in file "'//restart%filename//'"',FATAL)
  allocate(r(size(restart%cidx)))
  call read_data(restart%rhandle, varname, r, unlim_dim_level=1)
  call distrib_cohort_data_i0d(fptr,restart%cidx,restart%tile_dim_length,r)
  deallocate(r)
end subroutine get_int_cohort_data

#define F90_TYPE real
#define NF_TYPE NF_DOUBLE
#define NF_FILL_VALUE NF_FILL_DOUBLE
#define READ_0D_FPTR read_cohort_data_r0d_fptr
#define WRITE_0D_FPTR write_cohort_data_r0d_fptr
#define WRITE_0D write_cohort_data_r0d
#include "vegn_cohort_io.inc"

#define F90_TYPE integer
#define NF_TYPE NF_INT
#define NF_FILL_VALUE NF_FILL_INT
#define READ_0D_FPTR read_cohort_data_i0d_fptr
#define WRITE_0D_FPTR write_cohort_data_i0d_fptr
#define WRITE_0D write_cohort_data_i0d
#include "vegn_cohort_io.inc"

end module cohort_io_mod
