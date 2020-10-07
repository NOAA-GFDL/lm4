module cohort_io_mod

use netcdf, only: NF90_FILL_DOUBLE, NF90_FILL_INT

use fms_mod,          only : error_mesg, FATAL, WARNING
use fms_io_mod,       only : get_instance_filename
use mpp_mod,          only : mpp_max
use land_io_mod,      only : input_buf_size
use land_tile_mod,    only : land_tile_map, land_tile_type, land_tile_list_type, &
     land_tile_enum_type, first_elmt, tail_elmt, next_elmt, &
     current_tile, operator(/=), nitems, loop_over_tiles

use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, get_tile_data, &
     get_tile_by_idx

use vegn_cohort_mod, only: vegn_cohort_type
use land_data_mod, only : lnd


use fms2_io_mod, only: compressed_start_and_count, FmsNetcdfUnstructuredDomainFile_t, &
                       register_axis, register_field, register_variable_attribute, &
                       read_data, write_data
use mpp_mod, only : mpp_chksum                       
use land_chksum_mod 

implicit none
private

! ==== public interfaces =====================================================
public :: read_create_cohorts
public :: create_cohort_dimension
public :: add_cohort_data, add_int_cohort_data
public :: get_cohort_data, get_int_cohort_data
! remove when cleaning up:
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
  integer :: ncidx
  integer, dimension(:), allocatable :: npes_cidx !Cohort index length of each pe in file's pelist.
  integer, dimension(:), allocatable :: npes_cidx_start !Offset of cohort index of each pe in file's pelist.
  integer, dimension(cohorts_dim_length) :: buffer
  integer :: i

  ! form the full name of the file
  call get_instance_filename(trim(name), file_name)

  ! the size of tile dimension really does not matter for the output, but it does
  ! matter for uncompressing utility, since it uses it as a size of the array to
  ! unpack to create tile index dimension and variable.
  call register_axis(rhandle, "cohort", cohorts_dim_length)
  call register_field(rhandle, "cohort", "int", (/"cohort"/))
  call register_variable_attribute(rhandle, "cohort", "long_name", "cohort number within tile", &
                                   str_len=trim("cohort number within tile"))
  do i = 1, cohorts_dim_length
    buffer(i) = i
  enddo
  call write_data(rhandle, "cohort", buffer)

  ncidx =  size(cidx)
  call compressed_start_and_count(rhandle, ncidx, npes_cidx_start, npes_cidx)
  call register_axis(rhandle, cohort_index_name, npes_corner=npes_cidx_start, npes_nelems=npes_cidx)
  deallocate(npes_cidx_start)
  deallocate(npes_cidx)
  call register_field(rhandle, cohort_index_name, "int", (/cohort_index_name/))
  call register_variable_attribute(rhandle, cohort_index_name, "compress", "cohort tile lat lon", &
                                   str_len=trim("cohort tile lat lon"))
  call register_variable_attribute(rhandle, cohort_index_name, "units", "none", str_len=trim("none"))
  call register_variable_attribute(rhandle, cohort_index_name, "long_name", "compressed vegetation cohort index", &
                                   str_len=trim("compressed vegetation cohort index"))
  call register_variable_attribute(rhandle, cohort_index_name, "valid_min", 0)
  call write_data(rhandle, cohort_index_name, cidx)
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
     data(i) = NF90_FILL_INT
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
     data(i) = NF90_FILL_DOUBLE
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
  character(len=32) :: chksum
  
  real, pointer :: r(:)

  allocate(r(size(restart%cidx)))
  call gather_cohort_data_r0d(fptr,restart%cidx,restart%tile_dim_length,r)
  call register_field(restart%rhandle, varname, "double", (/cohort_index_name, "Time"/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
  if (present(units)) then
    call register_variable_attribute(restart%rhandle, varname, "units", trim(units), str_len=trim(units))
  endif
  if (present(longname)) then
    call register_variable_attribute(restart%rhandle, varname, "long_name", trim(longname), str_len=trim(longname))
  endif

  call get_land_chksum_r0d(r,chksum)
  call register_variable_attribute(restart%rhandle, varname, "checksum", trim(chksum), str_len=trim(chksum))
  call write_data(restart%rhandle, varname, r)
  deallocate(r)
end subroutine add_cohort_data

! ============================================================================
subroutine add_int_cohort_data(restart,varname,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(cptr_i0)           :: fptr ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  character(len=32) :: chksum
  integer, pointer :: r(:)
  integer :: id_restart

  allocate(r(size(restart%cidx)))
  call gather_cohort_data_i0d(fptr,restart%cidx,restart%tile_dim_length,r)
  call register_field(restart%rhandle, varname, "int", (/cohort_index_name, "Time"/))
  call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_INT)
  if (present(units)) then
    call register_variable_attribute(restart%rhandle, varname, "units", trim(units), str_len=trim(units))
  endif
  if (present(longname)) then
    call register_variable_attribute(restart%rhandle, varname, "long_name", trim(longname), str_len=trim(longname))
  endif
  
  call get_land_chksum_i0d(r,chksum)
  call register_variable_attribute(restart%rhandle, varname, "checksum", trim(chksum), str_len=trim(chksum))
  call write_data(restart%rhandle, varname, r)
  deallocate(r)
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
  call read_data(restart%rhandle, varname, r)
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
  call read_data(restart%rhandle, varname, r)
  call distrib_cohort_data_i0d(fptr,restart%cidx,restart%tile_dim_length,r)
  deallocate(r)
end subroutine get_int_cohort_data

end module cohort_io_mod
