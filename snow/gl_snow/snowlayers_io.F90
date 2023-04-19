module snowlayers_io_mod

use netcdf, only: NF90_FILL_DOUBLE, NF90_FILL_INT

use fms_mod,          only : error_mesg, FATAL, WARNING
use fms_io_mod,       only : get_instance_filename
use fms2_io_mod, only: FmsNetcdfUnstructuredDomainFile_t, compressed_start_and_count, &
     register_axis, register_field, register_variable_attribute, read_data, write_data
use mpp_mod,          only : mpp_max
use land_data_mod, only    : lnd
use land_io_mod,      only : register_variable_string_attribute
use land_tile_io_mod, only : land_restart_type, get_tile_by_idx
use land_tile_mod,    only : land_tile_map, land_tile_type, &
     land_tile_enum_type, first_elmt, tail_elmt, next_elmt, &
     current_tile, operator(/=), loop_over_tiles
use snowpack_mod, only : snow_layer_type


implicit none
private

! ==== public interfaces =====================================================
public :: read_create_snowlayers
public :: create_snowlayer_dimension
public :: add_snowlayer_data, add_int_snowlayer_data
public :: get_snowlayer_data, get_int_snowlayer_data
! remove when cleaning up:
public :: gather_snowlayer_index, gather_snowlayer_data
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
         ptr=>tile%snow%sp%snow(k+1)
      endif
   endif
end subroutine get_snowlayer_by_idx

! ! ============================================================================
! subroutine read_create_snowlayers(restart)
!   type(land_restart_type), intent(inout) :: restart

!   if (new_land_io) then
!      if (.not.allocated(restart%cidx)) call error_mesg('read_create_snowlayers', &
!         'snow layers index not found in file "'//restart%filename//'"',FATAL)
!      call read_create_snowlayers_new(restart%cidx,restart%tile_dim_length)
!   else
!      call read_create_snowlayers_orig(restart%ncid,restart%filename)
!   endif
! end subroutine

! ! ============================================================================
! subroutine read_create_snowlayers_orig(ncid, filename)
!   integer, intent(in) :: ncid
!   character(*), intent(in) :: filename

!   integer :: nsnowlayers ! total number of cohorts in restart file
!   integer :: nlon, nlat, ntiles ! size of respective dimensions

!   integer, allocatable :: idx(:)
!   integer :: i,j,t,k,m, n, nn, idxid, ierr
!   integer :: bufsize, npts,g,l
!   type(land_tile_enum_type) :: ce, te
!   type(land_tile_type), pointer :: tile
!   character(len=64) :: info ! for error message

!   ! get the size of dimensions
!   nlon = lnd%nlon ; nlat = lnd%nlat
!   ierr = nfu_inq_dim(ncid,'tile',len=ntiles)
!   if (ierr/=NF_NOERR) call error_mesg('read_create_snowlayers_orig', &
!               'dimension "tile" not found in file "'//trim(filename)//'"', FATAL)

!   ! read the cohort index
!   ierr = nfu_inq_dim(ncid,snowlayer_index_name,len=nsnowlayers)
!   if (ierr/=NF_NOERR) call error_mesg('read_create_snowlayers_orig', &
!               'dimension "'//trim(snowlayer_index_name)//'" not found in file "'//trim(filename)//'"', FATAL)
!   ierr = nfu_inq_var(ncid,snowlayer_index_name,id=idxid)
!   if (ierr/=NF_NOERR) call error_mesg('read_create_snowlayers_orig', &
!               'variable "'//trim(snowlayer_index_name)//'" not found in file "'//trim(filename)//'"', FATAL)
!   bufsize = min(input_buf_size,nsnowlayers)
!   allocate(idx(bufsize))

!   npts = nlon*nlat
!   do nn = 1, nsnowlayers, bufsize
!      __NF_ASRT__(nf_get_vara_int(ncid,idxid,nn,min(bufsize,nsnowlayers-nn+1),idx))

!      do n = 1,min(bufsize,nsnowlayers-nn+1)
!         k = idx(n)
!         g = modulo(k,npts)+1
!         if(g<lnd%gs.or.g>lnd%ge) cycle ! skip points outside of domain
!         l = lnd%l_index(g)
!         k = k/npts
!         t = modulo(k,ntiles)+1 ; k = k/ntiles
!         k = k+1
!         ce = first_elmt(land_tile_map(l))
!         do m = 1,t-1
!            ce=next_elmt(ce)
!         enddo
!         tile=>current_tile(ce)
!       !   if(.not.associated(tile%vegn)) then
!         if(.not.associated(tile%snow)) then
!            i = lnd%i_index(i)
!            j = lnd%j_index(j)
!            info = ''
!            write(info,'("(",3i3,")")')i,j,t
!            call error_mesg('read_create_snowlayer',&
!                 'snow tile'//trim(info)//' does not exist, but is necessary to create a snowlayer', &
!                 WARNING)
!         else
!            ! tile%snow%nlayers = tile%snow%nlayers + 1
!            tile%snow%sp%nlayers = tile%snow%sp%nlayers + 1
!         endif
!      enddo
!   enddo

!   ! go through all tiles in the domain and allocate requested numner of cohorts
!   ce = first_elmt(land_tile_map); te = tail_elmt(land_tile_map)
!   do while (ce/=te)
!       tile=>current_tile(ce); ce = next_elmt(ce)
!       if(.not.associated(tile%snow))cycle
!       ! write(*,*) "trying to allocate snow tile ... "
!       ! write(*,*) "nlayers = ", tile%snow%sp%nlayers
!       ! write(*,*) "size(tile%snow&sp%snow) = ", size(tile%snow%sp%snow)
!       !   write(*,*) "tile%snow%snow = ", tile%snow%snow
!       ! if(tile%snow%sp%nlayers>0) then ! EZSNOW 
!          allocate(tile%snow%sp%snow(tile%snow%sp%nlayers))
!       ! else ! allocate at least one slot if no snow ! EZSNOW
!          ! allocate(tile%snow%sp%snow(1)) ! EZSNOW
!       ! endif
!   enddo

!   ! clean up memory
!   deallocate(idx)
! end subroutine read_create_snowlayers_orig

! ! ============================================================================
! subroutine read_create_snowlayers_new(idx,ntiles)
!   integer, intent(in) :: idx(:)
!   integer, intent(in) :: ntiles

! !   integer :: ncohorts ! total number of cohorts in restart file
!   integer :: nsnowlayers ! total number of cohorts in restart file
!   integer :: nlon, nlat ! size of respective dimensions

!   integer :: i,j,t,k,m, n, npts, g, l
!   type(land_tile_enum_type) :: ce, te
!   type(land_tile_type), pointer :: tile
!   character(len=64) :: info ! for error message

!   ! get the size of dimensions
!   nlon = lnd%nlon
!   nlat = lnd%nlat
! !   ncohorts = size(idx)
!   nsnowlayers = size(idx)
!   npts = nlon*nlat

!   do n = 1,nsnowlayers
!      if(idx(n)<0) cycle ! skip illegal indices
!      k = idx(n)
!      g = modulo(k,npts)+1
!      if(g<lnd%gs.or.g>lnd%ge) cycle ! skip points outside of domain
!      l = lnd%l_index(g)
!      k = k/npts
!      t = modulo(k,ntiles)+1 ; k = k/ntiles
!      k = k+1

!      ce = first_elmt(land_tile_map(l))
!      do m = 1,t-1
!         ce=next_elmt(ce)
!      enddo
!      tile=>current_tile(ce)

!      if (.not. associated(tile)) then
!          call error_mesg("read_create_snowlayers_new", &
!                          "current tile returned null pointer", &
!                          FATAL)
!      endif

!      if(.not.associated(tile%snow)) then
!         info = ''
!         write(info,'("(",3i3,")")')i,j,t
!         call error_mesg('read_create_snowlayer',&
!              'snow tile'//trim(info)//' does not exist, but is necessary to create a snow layer', &
!              WARNING)
!      else
!         tile%snow%sp%nlayers = tile%snow%sp%nlayers + 1
!       !   tile%vegn%n_cohorts = tile%vegn%n_cohorts + 1
!      endif
!   enddo

!   ! go through all tiles in the domain and allocate requested numner of cohorts
!   ce = first_elmt(land_tile_map); te = tail_elmt(land_tile_map)
!   do while (ce/=te)
!      tile=>current_tile(ce); ce = next_elmt(ce)
!      if(.not.associated(tile%snow))cycle
!    !   allocate(tile%vegn%cohorts(tile%vegn%n_cohorts))
!      allocate(tile%snow%sp%snow(tile%snow%sp%nlayers))
!   enddo
! end subroutine read_create_snowlayers_new


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
         !  tile%vegn%n_cohorts = tile%vegn%n_cohorts + 1
         tile%snow%sp%nlayers = tile%snow%sp%nlayers + 1
      endif
   enddo
 
   ! go through all tiles in the domain and allocate requested numner of snowlayers
   ce = first_elmt(land_tile_map); te = tail_elmt(land_tile_map)
   do while (ce/=te)
      tile=>current_tile(ce); ce = next_elmt(ce)
      if(.not.associated(tile%snow))cycle
!    !   allocate(tile%vegn%cohorts(tile%vegn%n_cohorts))
      allocate(tile%snow%sp%snow(tile%snow%sp%nlayers))
   enddo
 end subroutine read_create_snowlayers

! ! ============================================================================
! ! creates cohort dimension, if necessary, in the output restart file. NOTE
! subroutine create_snowlayer_dimension(restart)
!   type(land_restart_type), intent(inout) :: restart

!   if (new_land_io) then
!      call create_snowlayer_dimension_new(restart%rhandle,restart%cidx,restart%basename,restart%tile_dim_length)
!   else
!      call create_snowlayer_dimension_orig(restart%ncid,restart%cidx,restart%tile_dim_length)
!   endif
! end subroutine create_snowlayer_dimension

! ! ============================================================================
! ! creates cohort dimension, if necessary, in the output restart file. NOTE
! ! that this subroutine should be called even if restart has not been created
! ! (because, for example, there happen to be no vegetation in a certain domain),
! ! for the reason that it calls mpp_max, and that should be called for each
! ! processor to work.
! subroutine create_snowlayer_dimension_orig(ncid,cidx,tile_dim_length)
!   integer, intent(in) :: ncid
!   integer, allocatable, intent(out) :: cidx(:)
!   integer, intent(in) :: tile_dim_length

!   ! ---- local vars
!   integer :: i,k,max_snowlayers,p
!   integer :: iret
!   integer, allocatable :: nsnowlayers(:) ! array of idx sizes from all PEs in io_domain
!   integer, allocatable :: idx2(:) ! array of cohort indices from all PEs in io_domain

!   call gather_snowlayer_index(tile_dim_length,cidx)
!   max_snowlayers = global_max_snowlayers()
! !   max_cohorts = global_max_cohorts()

!   if (mpp_pe()/=lnd%io_pelist(1)) then
!      ! if this processor is not doing io (that is, it's not root io_domain
!      ! processor), simply send the data to the root io_domain PE
!      call mpp_send(size(cidx), plen=1,          to_pe=lnd%io_pelist(1), tag=COMM_TAG_1)
!      call mpp_send(cidx(1),    plen=size(cidx), to_pe=lnd%io_pelist(1), tag=COMM_TAG_2)
!   else
!      ! gather the array of cohort index sizes
!    !   allocate(ncohorts(size(lnd%io_pelist)))
!      allocate(nsnowlayers(size(lnd%io_pelist)))
!      nsnowlayers(1) = size(cidx)
!      do p = 2,size(lnd%io_pelist)
!         call mpp_recv(nsnowlayers(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_1)
!      enddo
!      ! gather cohort index from the processors in our io_domain
!      allocate(idx2(sum(nsnowlayers(:))))
!      idx2(1:nsnowlayers(1))=cidx(:)
!      k=nsnowlayers(1)+1
!      do p = 2,size(lnd%io_pelist)
!         call mpp_recv(idx2(k), from_pe=lnd%io_pelist(p), glen=nsnowlayers(p), tag=COMM_TAG_2)
!         k = k+nsnowlayers(p)
!      enddo
!      ! create cohort dimension in the output file
!      iret = nf_redef(ncid)
!      __NF_ASRT__(nfu_def_dim(ncid,'snowlayer',(/(i,i=1,max_snowlayers)/),'snowlayer number within tile'))
!      ! create cohort index
!      __NF_ASRT__(nfu_def_dim(ncid,snowlayer_index_name,idx2,'compressed snow layer index'))
!      __NF_ASRT__(nfu_put_att(ncid,snowlayer_index_name,'compress','snowlayer tile lat lon'))
!      __NF_ASRT__(nfu_put_att(ncid,snowlayer_index_name,'valid_min',0))

!      ! deallocate the data we no longer need
!      deallocate(nsnowlayers,idx2)
!      ! leave the define mode to commit the new definitions to the disk
!      iret = nf_enddef(ncid)
!   endif
!   call mpp_sync_self()
! end subroutine create_snowlayer_dimension_orig

! subroutine create_snowlayer_dimension_new(rhandle,cidx,name,tile_dim_length)
!   type(restart_file_type), intent(inout) :: rhandle ! restart file handle
!   integer, allocatable,    intent(out)   :: cidx(:) ! rank local tile index vector
!   character(len=*),        intent(in)    :: name    ! name of the restart file
!   integer,                 intent(in)    :: tile_dim_length ! length of tile axis

!   integer :: max_snowlayers

! !   call gather_cohort_index(tile_dim_length,cidx)
! !   max_cohorts = global_max_cohorts()
! !   call create_cohort_out_file_idx(rhandle,name,cidx,max(max_cohorts,1))

!   call gather_snowlayer_index(tile_dim_length,cidx)
!   max_snowlayers = global_max_snowlayers()
!   call create_snowlayer_out_file_idx(rhandle,name,cidx,max(max_snowlayers,1))
!   write(*,*) "creating index out output snowlayers:: cidx = ", cidx

! end subroutine create_snowlayer_dimension_new

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

! subroutine create_snowlayer_out_file_idx(rhandle,name,cidx,snowlayers_dim_length)
!   type(restart_file_type), intent(inout) :: rhandle     ! restart file handle
!   character(len=*),      intent(in)  :: name                ! name of the file to create
!   integer              , intent(in)  :: cidx(:)             ! integer compressed index of tiles (local)
!   integer              , intent(in)  :: snowlayers_dim_length  ! length of cohorts axis

!   ! ---- local vars
!   character(256) :: file_name ! full name of the file, including the processor number

!   ! form the full name of the file
!   call get_instance_filename(trim(name), file_name)
!   call get_mosaic_tile_file(trim(file_name),file_name,lnd%ug_domain)

!   ! the size of tile dimension really does not matter for the output, but it does
!   ! matter for uncompressing utility, since it uses it as a size of the array to
!   ! unpack to create tile index dimension and variable.
!   call fms_io_unstructured_register_restart_axis(rhandle, &
!                                                  name, &
!                                                 !  trim(cohort_index_name), &
!                                                  trim(snowlayer_index_name), &
!                                                  cidx, &
!                                                  "snowlayer tile lat lon", &
!                                                  "H", &
!                                                  snowlayers_dim_length, &
!                                                  lnd%ug_domain, &
!                                                  dimlen_name="snowlayer", &
!                                                  dimlen_lname="snowlayer number within tile", &
!                                                  units="none", &
!                                                  longname="compressed vegetation snowlayer index", &
!                                                  imin=0)

! end subroutine create_snowlayer_out_file_idx

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
        global_max_snowlayers = max(global_max_snowlayers,tile%snow%sp%nlayers)
  enddo
  call mpp_max(global_max_snowlayers)
end function global_max_snowlayers

subroutine gather_snowlayer_index(ntiles, cidx)
  integer,              intent(in)  :: ntiles
  integer, allocatable, intent(out) :: cidx(:)   ! integer compressed index of tiles

  integer :: i,j,k,c,n
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile

  ! count total number of cohorts in our compute domain
  ce = first_elmt(land_tile_map)
  n = 0
  do while (loop_over_tiles(ce,tile))
   !   if(associated(tile%vegn)) n = n+tile%vegn%n_cohorts
     if(associated(tile%snow)) n = n + tile%snow%sp%nlayers
  enddo

  ! calculate compressed cohort index to be written to the restart file
  allocate(cidx(max(n,1))) ; cidx(:) = -1
  ce = first_elmt(land_tile_map, lnd%ls)
  n = 1
  do while (loop_over_tiles(ce,tile,i=i,j=j,k=k))
     if(associated(tile%snow)) then
        do c = 1,tile%snow%sp%nlayers
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
   !   call get_cohort_by_idx ( idx(i), ntiles, snowlayer)
     call get_snowlayer_by_idx ( idx(i), ntiles, snowlayer)
     data(i) = NF90_FILL_DOUBLE
     if (associated(snowlayer)) then
        call fptr(snowlayer, ptr)
        if(associated(ptr)) data(i) = ptr
     endif
  enddo
end subroutine gather_snowlayer_data_r0d

! EZSNOW UPDATED BELOW

! ! ============================================================================
! subroutine add_snowlayer_data(restart,varname,fptr,longname,units)
!   type(land_restart_type), intent(inout) :: restart
!   character(len=*), intent(in) :: varname ! name of the variable to write
!   procedure(cptr_r0)           :: fptr ! subroutine returning pointer to the data
!   character(len=*), intent(in), optional :: units, longname

!   real, pointer :: r(:)
!   integer :: id_restart

!   allocate(r(size(restart%cidx)))
!   call gather_snowlayer_data_r0d(fptr,restart%cidx,restart%tile_dim_length,r)
!   if (new_land_io) then
!      id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, &
!           restart%basename, varname, r, (/HIDX/), lnd%ug_domain, &
!           longname=longname, units=units, restart_owns_data=.true.)
!   else
!      call write_snowlayer_data_r0d(restart%ncid,varname,r,longname,units)
!      deallocate(r)
!   endif
! end subroutine add_snowlayer_data

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

! ! ============================================================================
! subroutine add_int_snowlayer_data(restart,varname,fptr,longname,units)
!   type(land_restart_type), intent(inout) :: restart
!   character(len=*), intent(in) :: varname ! name of the variable to write
!   procedure(cptr_i0)           :: fptr ! subroutine returning pointer to the data
!   character(len=*), intent(in), optional :: units, longname

!   integer, pointer :: r(:)
!   integer :: id_restart

!   allocate(r(size(restart%cidx)))
!   call gather_snowlayer_data_i0d(fptr,restart%cidx,restart%tile_dim_length,r)
!   if (new_land_io) then
!      id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, &
!          restart%basename, varname, r, (/HIDX/), lnd%ug_domain, &
!          longname=longname, units=units, restart_owns_data=.true.)
!   else
!      call write_snowlayer_data_i0d(restart%ncid,varname,r,longname,units)
!      deallocate(r)
!   endif
! end subroutine add_int_snowlayer_data

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

! ! ========================================================================================
! subroutine get_snowlayer_data(restart,varname,fptr)
!   type(land_restart_type), intent(in) :: restart
!   character(len=*), intent(in) :: varname ! name of the variable to write
!   procedure(cptr_r0)           :: fptr ! subroutine returning pointer to the data

!   real, allocatable :: r(:)
!   if (new_land_io) then
!      if (.not.allocated(restart%cidx)) call error_mesg('read_create_snowlayers', &
!         'snowlayer index not found in file "'//restart%filename//'"',FATAL)
!      allocate(r(size(restart%cidx)))
!      call fms_io_unstructured_read(restart%basename, varname, r, lnd%ug_domain, timelevel=1)
!      call distrib_snowlayer_data_r0d(fptr,restart%cidx,restart%tile_dim_length,r)
!      deallocate(r)
!   else
!      call read_snowlayer_data_r0d_fptr(restart%ncid,varname,fptr)
!   endif
! end subroutine get_snowlayer_data


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

! ! ============================================================================
! subroutine get_int_snowlayer_data(restart,varname,fptr)
!   type(land_restart_type), intent(in) :: restart
!   character(len=*), intent(in) :: varname ! name of the variable to write
!   procedure(cptr_i0)           :: fptr ! subroutine returning pointer to the data

!   integer, allocatable :: r(:)
!   if (new_land_io) then
!      if (.not.allocated(restart%cidx)) call error_mesg('read_create_snowlayers', &
!         'snowlayer index not found in file "'//restart%filename//'"',FATAL)
!      allocate(r(size(restart%cidx)))
!      call fms_io_unstructured_read(restart%basename, varname, r, lnd%ug_domain, timelevel=1)
!      call distrib_snowlayer_data_i0d(fptr,restart%cidx,restart%tile_dim_length,r)
!      deallocate(r)
!   else
!      call read_snowlayer_data_i0d_fptr(restart%ncid,varname,fptr)
!   endif
! end subroutine get_int_snowlayer_data

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
