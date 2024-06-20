module snowlayers_io_mod

use netcdf, only: NF90_FILL_DOUBLE, NF90_FILL_INT
use fms_mod,          only : error_mesg, FATAL, WARNING
use fms2_io_mod, only: FmsNetcdfUnstructuredDomainFile_t, compressed_start_and_count, &
     register_axis, register_field, register_variable_attribute, read_data, write_data, get_instance_filename
use mpp_mod,          only : mpp_max
use land_data_mod, only    : lnd
use land_io_mod,      only : register_variable_string_attribute
use land_tile_io_mod, only : land_restart_type, get_tile_by_idx
use land_tile_mod,    only : land_tile_map, land_tile_type, &
     land_tile_enum_type, first_elmt, tail_elmt, next_elmt, &
     current_tile, operator(/=), loop_over_tiles
use gl_snow_tile_mod, only : gl_snow_tile_type
use snowpack_mod, only : snow_layer_type


implicit none
private

! ==== public interfaces =====================================================
public :: read_create_snowlayers
public :: create_snowlayer_dimension
public :: add_snowlayer_data, add_int_snowlayer_data
public :: get_snowlayer_data, get_int_snowlayer_data
! remove when cleaning up:
! public :: gather_snowlayer_index, gather_snowlayer_data
! ==== end of public interfaces ==============================================

interface create_snowlayer_dimension
   module procedure create_snowlayer_dimension1
   module procedure create_snowlayer_dimension2
end interface create_snowlayer_dimension

interface gather_snowlayer_data
   module procedure gather_snowlayer_data_r0d
   module procedure gather_snowlayer_data_i0d
end interface gather_snowlayer_data

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snowlayers_io_mod'
! name of the "compressed" dimension (and dimension variable) in the output
! netcdf files -- that is, the dimensions written out using compression by
! gathering, as described in CF conventions.
character(len=*),   parameter :: snowlayer_index_name   = 'snowlayer_index'

abstract interface
  ! given land snowlayer, returns pointer to some scalar real data
  ! within this snowlayer, or an unassociated pointer if there is no data
  subroutine cptr_r0(tile, ptr)
     import snow_layer_type
     type(snow_layer_type), pointer :: tile ! input
     real                , pointer :: ptr  ! returned pointer to the data
  end subroutine cptr_r0
  ! given land snowlayer, returns pointer to some scalar real data
  ! within this snowlayer, or an unassociated pointer if there is no data
  subroutine cptr_i0(tile, ptr)
     import snow_layer_type
     type(snow_layer_type), pointer :: tile ! input
     integer               , pointer :: ptr  ! returned pointer to the data
  end subroutine cptr_i0
end interface

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
! given compressed index, sizes of the global grid, 2D array of tile lists
! and the lower boundaries of this array, returns a pointer to the snowlayer
! corresponding to the compressed index, or NULL is the index is outside
! current domain, or such tile does not exist, or such snowlayer does not exist.
subroutine get_snowlayer_by_idx(idx,ntiles,ptr)
   integer, intent(in) :: idx ! index
   integer, intent(in) :: ntiles ! size of tile dimension
   type(snow_layer_type), pointer :: ptr

   ! ---- local vars
   integer :: tile_idx, k
   type(land_tile_type), pointer :: tile

   ptr=>NULL()
   if ( idx < 0 ) return
   tile_idx = modulo(idx,lnd%nlon*lnd%nlat*ntiles)
   call get_tile_by_idx(tile_idx,tile)
   if(associated(tile)) then
      if (associated(tile%snow)) then
         k = idx/(lnd%nlon*lnd%nlat*ntiles) ! calculate snowlayer index within a tile
         select type (s=>tile%snow); class is (gl_snow_tile_type)
            ptr=>s%sp%snow(k+1)
         end select
      endif
   endif
end subroutine get_snowlayer_by_idx


! ============================================================================
subroutine read_create_snowlayers(restart)
   type(land_restart_type), intent(inout) :: restart

   integer :: nsnowlayers ! total number of snowlayers in restart file
   integer :: ntiles   ! total number of tiles in restart file
   integer :: nlon, nlat ! size of respective dimensions

   integer :: i,j,t,k,m, n, npts, g, l
   type(land_tile_enum_type) :: ce, te
   type(land_tile_type), pointer :: tile
   character(len=64) :: info ! for error message

   if (.not.allocated(restart%cidx)) call error_mesg('read_create_snowlayers', &
       'snowlayer index not found in file "'//restart%filename//'"',FATAL)

   ! get the size of dimensions
   nlon = lnd%nlon
   nlat = lnd%nlat
   ntiles   = restart%tile_dim_length
   nsnowlayers = size(restart%cidx)
   npts = nlon*nlat

   do n = 1,nsnowlayers
      if(restart%cidx(n)<0) cycle ! skip illegal indices
      k = restart%cidx(n)
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
          call error_mesg("read_create_snowlayers", &
                          "current tile returned null pointer", &
                          FATAL)
      endif

      if(.not.associated(tile%snow)) then
         info = ''
         write(info,'("(",3i3,")")')i,j,t
         call error_mesg('read_create_snowlayer',&
              'snow tile'//trim(info)//' does not exist, but is necessary to create a snowlayer', &
              WARNING)
      else
         select type (s=>tile%snow); class is (gl_snow_tile_type)
            s%sp%nlayers = s%sp%nlayers + 1
         end select
      endif
   enddo

   ! go through all tiles in the domain and allocate requested numner of snowlayers
   ce = first_elmt(land_tile_map); te = tail_elmt(land_tile_map)
   do while (ce/=te)
      tile=>current_tile(ce); ce = next_elmt(ce)
      if(.not.associated(tile%snow))cycle
      select type (s=>tile%snow); class is (gl_snow_tile_type)
         ! slm: we probably should not rely on snow%n_layers() function,
         ! because in general it can be invalid until the tile is fully initialized
         ! perhaps there is a better way, avoiding numerous "select type" constructs?
         allocate(s%sp%snow(s%sp%nlayers))
      end select
   enddo
 end subroutine read_create_snowlayers


! ============================================================================
! creates snowlayer dimension, if necessary, in the output restart file.
subroutine create_snowlayer_dimension1(restart)
   type(land_restart_type), intent(inout) :: restart
   call create_snowlayer_dimension2(restart%rhandle,restart%cidx,restart%basename,restart%tile_dim_length)
 end subroutine create_snowlayer_dimension1

 ! ============================================================================
 ! creates snowlayer dimension, if necessary, in the output restart file. NOTE
 ! that this subroutine should be called even if restart has not been created
 ! (because, for example, there happen to be no snow in a certain domain),
 ! for the reason that it calls mpp_max, and that should be called for each
 ! processor to work.
 subroutine create_snowlayer_dimension2(rhandle,cidx,name,tile_dim_length)
   type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: rhandle ! fms_io restart file data type
   integer, allocatable,    intent(out)   :: cidx(:) ! rank local tile index vector
   character(len=*),        intent(in)    :: name    ! name of the restart file
   integer,                 intent(in)    :: tile_dim_length ! length of tile axis

   integer :: max_snowlayers

   call gather_snowlayer_index(tile_dim_length,cidx)
   max_snowlayers = global_max_snowlayers()

   call create_snowlayer_out_file_idx(rhandle,name,cidx,max(max_snowlayers,1))
 end subroutine create_snowlayer_dimension2


subroutine create_snowlayer_out_file_idx(rhandle,name,cidx,snowlayers_dim_length)
   type(FmsNetcdfUnstructuredDomainFile_t),intent(inout) :: rhandle ! fms_io restart file data type
   character(len=*),      intent(in)  :: name                ! name of the file to create
   integer              , intent(in)  :: cidx(:)             ! integer compressed index of tiles (local)
   integer              , intent(in)  :: snowlayers_dim_length  ! length of snowlayers axis

   ! ---- local vars
   character(256) :: file_name ! full name of the file, including the processor number
   integer :: ncidx
   integer, dimension(:), allocatable :: npes_cidx !snowlayer index length of each pe in file's pelist.
   integer, dimension(:), allocatable :: npes_cidx_start !Offset of snowlayer index of each pe in file's pelist.
   integer, dimension(snowlayers_dim_length) :: buffer
   integer :: i

   ! form the full name of the file
   call get_instance_filename(trim(name), file_name)

   ! the size of tile dimension really does not matter for the output, but it does
   ! matter for uncompressing utility, since it uses it as a size of the array to
   ! unpack to create tile index dimension and variable.

   call register_axis(rhandle, "snowlayer", snowlayers_dim_length)
   call register_field(rhandle, "snowlayer", "int", (/"snowlayer"/))
   call register_variable_string_attribute(rhandle, "snowlayer", "long_name", "snowlayer number within tile")
   do i = 1, snowlayers_dim_length
     buffer(i) = i
   enddo
   call write_data(rhandle, "snowlayer", buffer)

   ncidx =  size(cidx)
   call compressed_start_and_count(rhandle, ncidx, npes_cidx_start, npes_cidx)
   call register_axis(rhandle, trim(snowlayer_index_name), npes_corner=npes_cidx_start, npes_nelems=npes_cidx)
   deallocate(npes_cidx_start)
   deallocate(npes_cidx)
   call register_field(rhandle, trim(snowlayer_index_name), "int", (/trim(snowlayer_index_name)/))
   call register_variable_string_attribute(rhandle, trim(snowlayer_index_name), "compress", "snowlayer tile lat lon")
   call register_variable_string_attribute(rhandle, trim(snowlayer_index_name), "units", "none")
   call register_variable_string_attribute(rhandle, trim(snowlayer_index_name), "long_name", "compressed vegetation snowlayer index")
   call register_variable_attribute(rhandle, trim(snowlayer_index_name), "valid_min", 0)
   call write_data(rhandle, trim(snowlayer_index_name), cidx)

 end subroutine create_snowlayer_out_file_idx

subroutine distrib_snowlayer_data_i0d(fptr,idx,ntiles,data)
  integer, intent(in) :: idx(:) ! local vector of snowlayer indices
  integer, intent(in) :: ntiles ! size of the tile dimension
  integer, intent(in) :: data(:) ! local snowlayer data
  procedure(cptr_i0) :: fptr ! subroutine returning pointer to the data

  ! ---- local vars
  type(snow_layer_type), pointer :: snowlayer
  integer, pointer :: ptr ! pointer to the individual snowlayer data
  integer :: mask(size(data)) ! mask of valid data
  integer :: i

  ! gather data into an array along the snowlayer dimension
  do i = 1, size(idx)
     call get_snowlayer_by_idx ( idx(i), ntiles, snowlayer)
     if (associated(snowlayer)) then
        call fptr(snowlayer, ptr)
        if(associated(ptr)) ptr = data(i)
     endif
  enddo
end subroutine distrib_snowlayer_data_i0d

subroutine distrib_snowlayer_data_r0d(fptr,idx,ntiles,data)
  integer, intent(in) :: idx(:) ! local vector of snowlayer indices
  integer, intent(in) :: ntiles ! size of the tile dimension
  real, intent(in) :: data(:) ! local snowlayer data
  procedure(cptr_r0) :: fptr ! subroutine returning pointer to the data

  ! ---- local vars
  type(snow_layer_type), pointer :: snowlayer
  real, pointer :: ptr ! pointer to the individual snowlayer data
  integer :: mask(size(data))
  integer :: i

  ! gather data into an array along the snowlayer dimension
  do i = 1, size(idx)
     call get_snowlayer_by_idx ( idx(i), ntiles, snowlayer)
     if (associated(snowlayer)) then
        call fptr(snowlayer, ptr)
        if(associated(ptr)) ptr = data(i)
     endif
  enddo
end subroutine distrib_snowlayer_data_r0d

! count max number of snowlayers per tile
integer function global_max_snowlayers()
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile

  ce = first_elmt(land_tile_map)
  global_max_snowlayers = 0
  do while (loop_over_tiles(ce,tile))
     if(associated(tile%snow)) &
        global_max_snowlayers = max(global_max_snowlayers,tile%snow%n_layers())
  enddo
  call mpp_max(global_max_snowlayers)
end function global_max_snowlayers

subroutine gather_snowlayer_index(ntiles, cidx)
  integer,              intent(in)  :: ntiles
  integer, allocatable, intent(out) :: cidx(:)   ! integer compressed index of tiles

  integer :: i,j,k,c,n
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile

  ! count total number of snowlayers in our compute domain
  ce = first_elmt(land_tile_map)
  n = 0
  do while (loop_over_tiles(ce,tile))
     if(associated(tile%snow)) n = n + tile%snow%n_layers()
  enddo

  ! calculate compressed snowlayer index to be written to the restart file
  allocate(cidx(max(n,1))) ; cidx(:) = -1
  ce = first_elmt(land_tile_map, lnd%ls)
  n = 1
  do while (loop_over_tiles(ce,tile,i=i,j=j,k=k))
     if(associated(tile%snow)) then
        do c = 1,tile%snow%n_layers()
           cidx (n) = &
                (c-1)*lnd%nlon*lnd%nlat*ntiles + &
                (k-1)*lnd%nlon*lnd%nlat + &
                (j-1)*lnd%nlon + &
                (i-1)
           n = n+1
        enddo
     endif
  end do
end subroutine gather_snowlayer_index

subroutine gather_snowlayer_data_i0d(fptr,idx,ntiles,data)
  procedure(cptr_i0) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:) ! local vector of snowlayer indices
  integer, intent(in) :: ntiles ! size of the tile dimension
  integer, intent(out) :: data(:) ! local snowlayer data

  ! ---- local vars
  type(snow_layer_type), pointer :: snowlayer
  integer, pointer :: ptr ! pointer to the individual snowlayer data
  integer :: i

  ! gather data into an array along the snowlayer dimension
  do i = 1, size(idx)
     call get_snowlayer_by_idx ( idx(i), ntiles, snowlayer)
     data(i) = NF90_FILL_INT
     if (associated(snowlayer)) then
        call fptr(snowlayer, ptr)
        if(associated(ptr)) data(i) = ptr
     endif
  enddo
end subroutine gather_snowlayer_data_i0d

subroutine gather_snowlayer_data_r0d(fptr,idx,ntiles,data)
  procedure(cptr_r0)   :: fptr ! subroutine returning pointer to the data
  integer, intent(in)  :: idx(:) ! local vector of snowlayer indices
  integer, intent(in)  :: ntiles ! size of the tile dimension
  real,    intent(out) :: data(:) ! local snowlayer data

  ! ---- local vars
  type(snow_layer_type), pointer :: snowlayer
  real, pointer :: ptr ! pointer to the individual snowlayer data
  integer :: i

  ! gather data into an array along the snowlayer dimension
  do i = 1, size(idx)
     call get_snowlayer_by_idx ( idx(i), ntiles, snowlayer)
     data(i) = NF90_FILL_DOUBLE
     if (associated(snowlayer)) then
        call fptr(snowlayer, ptr)
        if(associated(ptr)) data(i) = ptr
     endif
  enddo
end subroutine gather_snowlayer_data_r0d


! ============================================================================
subroutine add_snowlayer_data(restart,varname,fptr,longname,units)
   type(land_restart_type), intent(inout) :: restart
   character(len=*), intent(in) :: varname ! name of the variable to write
   procedure(cptr_r0)           :: fptr ! subroutine returning pointer to the data
   character(len=*), intent(in), optional :: units, longname

   real, pointer :: r(:)

   allocate(r(size(restart%cidx)))
   call gather_snowlayer_data_r0d(fptr,restart%cidx,restart%tile_dim_length,r)
   call register_field(restart%rhandle, varname, "double", (/snowlayer_index_name/))
   call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_DOUBLE)
   if (present(units)) then
     call register_variable_string_attribute(restart%rhandle, varname, "units", units)
   endif

   if (present(longname)) then
     call register_variable_string_attribute(restart%rhandle, varname, "long_name", longname)
   endif
   call write_data(restart%rhandle, varname, r)
   deallocate(r)

 end subroutine add_snowlayer_data


! ============================================================================
subroutine add_int_snowlayer_data(restart,varname,fptr,longname,units)
   type(land_restart_type), intent(inout) :: restart
   character(len=*), intent(in) :: varname ! name of the variable to write
   procedure(cptr_i0)           :: fptr ! subroutine returning pointer to the data
   character(len=*), intent(in), optional :: units, longname

   integer, pointer :: r(:)
   integer :: id_restart

   allocate(r(size(restart%cidx)))
   call gather_snowlayer_data_i0d(fptr,restart%cidx,restart%tile_dim_length,r)
   call register_field(restart%rhandle, varname, "int", (/snowlayer_index_name/))
   call register_variable_attribute(restart%rhandle, varname, "_FillValue", NF90_FILL_INT)
   if (present(units)) then
     call register_variable_string_attribute(restart%rhandle, varname, "units", units)
   endif
   if (present(longname)) then
     call register_variable_string_attribute(restart%rhandle, varname, "long_name", longname)
   endif
   call write_data(restart%rhandle, varname, r)
   deallocate(r)

 end subroutine add_int_snowlayer_data


! ============================================================================
subroutine get_snowlayer_data(restart,varname,fptr)
   type(land_restart_type), intent(in) :: restart
   character(len=*), intent(in) :: varname ! name of the variable to write
   procedure(cptr_r0)           :: fptr ! subroutine returning pointer to the data

   real, allocatable :: r(:)
   if (.not.allocated(restart%cidx)) call error_mesg('read_create_snowlayers', &
       'snowlayer index not found in file "'//restart%filename//'"',FATAL)
   allocate(r(size(restart%cidx)))
   call read_data(restart%rhandle, varname, r)
   call distrib_snowlayer_data_r0d(fptr,restart%cidx,restart%tile_dim_length,r)
   deallocate(r)

 end subroutine get_snowlayer_data


! ============================================================================
subroutine get_int_snowlayer_data(restart,varname,fptr)
   type(land_restart_type), intent(in) :: restart
   character(len=*), intent(in) :: varname ! name of the variable to write
   procedure(cptr_i0)           :: fptr ! subroutine returning pointer to the data

   integer, allocatable :: r(:)

   if (.not.allocated(restart%cidx)) call error_mesg('read_create_snowlayers', &
       'snowlayer index not found in file "'//restart%filename//'"',FATAL)
   allocate(r(size(restart%cidx)))
   call read_data(restart%rhandle, varname, r)
   call distrib_snowlayer_data_i0d(fptr,restart%cidx,restart%tile_dim_length,r)
   deallocate(r)
 end subroutine get_int_snowlayer_data


end module snowlayers_io_mod
