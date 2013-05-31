module cohort_io_mod

use fms_mod,          only : error_mesg, FATAL, WARNING
use fms_io_mod,       only : field_size, read_data
use mpp_mod,          only : mpp_pe, mpp_max, mpp_send, mpp_recv, mpp_sync, &
                             COMM_TAG_1, COMM_TAG_2, COMM_TAG_3
!use mpp_io_mod, only : &
!mpp_open, &
!                       MPP_OVERWR, &
!                       MPP_NETCDF, &
!                       mpp_get_axes, &
!                       mpp_close, &
!                       mpp_get_axis_index, &
!                       mpp_get_axis_length, &
!                       mpp_get_axis_data, &
!                       mpp_write_axis_data, &
!                       

use mpp_io_mod,       only : mpp_read, mpp_write, mpp_write_meta, &
                             mpp_get_axis_by_name, mpp_get_axis_data, &
                             fieldtype, axistype

use nf_utils_mod,     only : nfu_inq_dim, nfu_get_var, nfu_put_var, &
                             nfu_get_rec, nfu_put_rec, nfu_def_dim, &
                             nfu_def_var, nfu_put_att, nfu_inq_var
use land_io_mod,      only : print_netcdf_error
use land_tile_mod,    only : land_tile_type, land_tile_list_type, &
                             land_tile_enum_type, first_elmt, tail_elmt, &
                             next_elmt, get_elmt_indices, current_tile, &
                             operator(/=)
use land_tile_io_mod, only : create_tile_out_file, get_tile_by_idx, sync_nc_files

use vegn_cohort_mod,  only : vegn_cohort_type
use land_data_mod,    only : lnd, land_state_type

implicit none
private

! ==== public interfaces =====================================================
! input
public :: read_create_cohorts
public :: read_create_cohorts_new
public :: read_cohort_data_r0d_fptr, read_cohort_data_r0d_fptr_new
public :: read_cohort_data_i0d_fptr, read_cohort_data_i0d_fptr_new
! output
public :: create_cohort_dimension
public :: create_cohort_out_file
public :: write_cohort_data_r0d_fptr
public :: write_cohort_data_i0d_fptr
public :: write_cohort_data_r0d_fptr_new
public :: write_cohort_data_i0d_fptr_new

! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'cohort_io_mod', &
     version     = '$Id: vegn_cohort_io.F90,v 19.0.4.2.2.4 2013/03/19 14:21:41 William.Cooke Exp $', &
     tagname     = '$Name: siena_201305 $'
! name of the "compressed" dimension (and dimension variable) in the output 
! netcdf files -- that is, the dimensions written out using compression by 
! gathering, as described in CF conventions.
character(len=*),   parameter :: cohort_index_name   = 'cohort_index'
integer, parameter :: INPUT_BUF_SIZE = 2048 ! max size of the input buffer for
                                     ! cohort input

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
! given compressed index, sizes of the global grid, 2D array of tile lists
! and the lower boundaries of this array, returns a pointer to the cohort
! corresponding to the compressed index, or NULL is the index is outside 
! current domain, or such tile does not exist, or such cohort does not exist.
subroutine get_cohort_by_idx(idx,nlon,nlat,ntiles,tiles,is,js,ptr)
   integer, intent(in) :: idx ! index
   integer, intent(in) :: nlon, nlat, ntiles
   integer, intent(in) :: is, js
   type(land_tile_list_type), intent(in) :: tiles(is:,js:)
   type(vegn_cohort_type), pointer :: ptr
   
   ! ---- local vars
   integer :: tile_idx, k
   type(land_tile_type), pointer :: tile
   
   ptr=>NULL()
   
   tile_idx = modulo(idx,nlon*nlat*ntiles)
   call get_tile_by_idx(tile_idx,nlon,nlat,tiles,is,js,tile)
   if(associated(tile)) then
      if (associated(tile%vegn)) then
         k = idx/(nlon*nlat*ntiles) ! calculate cohort index within a tile
         ptr=>tile%vegn%cohorts(k+1)
      endif
   endif

end subroutine

! ============================================================================
subroutine read_create_cohorts(ncid)
  integer, intent(in) :: ncid

  integer :: ncohorts ! total number of cohorts in restart file
  integer :: nlon, nlat, ntiles ! size of respective dimensions
 
  integer, allocatable :: idx(:)
  integer :: i,j,t,k,m, n, nn, idxid
  integer :: bufsize
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  character(len=64) :: info ! for error message

  ! get the size of dimensions
  nlon = lnd%nlon ; nlat = lnd%nlat
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! read the cohort index
  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))
  __NF_ASRT__(nfu_inq_var(ncid,cohort_index_name,id=idxid))
  bufsize = min(INPUT_BUF_SIZE,ncohorts)
  allocate(idx(bufsize))
  
  do nn = 1, ncohorts, bufsize
     __NF_ASRT__(nf_get_vara_int(ncid,idxid,nn,min(bufsize,ncohorts-nn+1),idx))
     
     do n = 1,min(bufsize,ncohorts-nn+1)
        if(idx(n)<0) cycle ! skip illegal indices
        k = idx(n)
        i = modulo(k,nlon)+1   ; k = k/nlon
        j = modulo(k,nlat)+1   ; k = k/nlat
        t = modulo(k,ntiles)+1 ; k = k/ntiles
        k = k+1
        
        if (i<lnd%is.or.i>lnd%ie) cycle ! skip points outside of domain
        if (j<lnd%js.or.j>lnd%je) cycle ! skip points outside of domain
        
        ce = first_elmt(lnd%tile_map(i,j))
        do m = 1,t-1
           ce=next_elmt(ce)
        enddo
        tile=>current_tile(ce)
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
  enddo

  ! go through all tiles in the domain and allocate requested numner of cohorts
  ce = first_elmt(lnd%tile_map); te = tail_elmt(lnd%tile_map)
  do while (ce/=te)
     tile=>current_tile(ce); ce = next_elmt(ce)
     if(.not.associated(tile%vegn))cycle
     allocate(tile%vegn%cohorts(tile%vegn%n_cohorts))
  enddo

  ! clean up memory
  deallocate(idx)
end subroutine read_create_cohorts


! ============================================================================
subroutine read_create_cohorts_new(filename, ntiles)
  character(*), intent(in) :: filename ! name of the netcdf file
  integer, intent(out) :: ntiles
!  integer, intent(in) :: idx(:)

  integer :: ncohorts ! total number of cohorts in restart file
  integer :: nlon, nlat ! size of respective dimensions
 
  integer, allocatable :: idx(:)
  integer :: i,j,t,k,m, n, nn, idxid
!  integer :: bufsize
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  character(len=64) :: info ! for error message
  integer :: siz(4)!, ntiles
  logical :: found

  ! get the size of dimensions
  nlon = lnd%nlon ; nlat = lnd%nlat
!  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! read the cohort index
!  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))
!  __NF_ASRT__(nfu_inq_var(ncid,cohort_index_name,id=idxid))
!  bufsize = min(INPUT_BUF_SIZE,ncohorts)
!  allocate(idx(bufsize))
  
!  do nn = 1, ncohorts, bufsize
!     __NF_ASRT__(nf_get_vara_int(ncid,idxid,nn,min(bufsize,ncohorts-nn+1),idx))

  found = .false.
  call field_size(filename,'tile',siz, field_found=found, domain=lnd%domain)
  if ( .not. found ) &
    call error_mesg(trim(module_name),'tile axis not found in '//trim(filename), FATAL)
  ntiles=siz(1) ! 'tile' is an axis variable but does not have any data.

! Find and read the cohort axis data
  call field_size(filename,'cohort_index',siz, field_found=found, domain=lnd%domain)
  if ( .not. found ) &
    call error_mesg(trim(module_name),'cohort_index axis not found in '//trim(filename), FATAL)
  allocate(idx(siz(1)))
  call read_data(filename,'cohort_index',idx, domain=lnd%domain, is_compressed=.true.)
     
     do n = 1,size(idx)!min(bufsize,ncohorts-nn+1)
        if(idx(n)<0) cycle ! skip illegal indices
        k = idx(n)
        i = modulo(k,nlon)+1   ; k = k/nlon
        j = modulo(k,nlat)+1   ; k = k/nlat
        t = modulo(k,ntiles)+1 ; k = k/ntiles
        k = k+1
        
        if (i<lnd%is.or.i>lnd%ie) cycle ! skip points outside of domain
        if (j<lnd%js.or.j>lnd%je) cycle ! skip points outside of domain
        
        ce = first_elmt(lnd%tile_map(i,j))
        do m = 1,t-1
           ce=next_elmt(ce)
        enddo
        tile=>current_tile(ce)
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
!  enddo

  ! go through all tiles in the domain and allocate requested numner of cohorts
  ce = first_elmt(lnd%tile_map); te = tail_elmt(lnd%tile_map)
  do while (ce/=te)
     tile=>current_tile(ce); ce = next_elmt(ce)
     if(.not.associated(tile%vegn))cycle
     allocate(tile%vegn%cohorts(tile%vegn%n_cohorts))
  enddo

  ! clean up memory
!  deallocate(idx)
end subroutine read_create_cohorts_new

! ============================================================================
! creates cohort dimension, if necessary, in the output restart file. NOTE 
! that this subroutine should be called even if restart has not been created
! (because, for example, there happen to be no vegetation in a certain domain),
! for the reason that it calls mpp_max, and that should be called for each
! processor to work.
subroutine create_cohort_dimension(ncid)
  integer, intent(in) :: ncid


  ! ---- local vars
  type(land_tile_enum_type) :: ce, te ! tile list enumerators
  type(land_tile_type), pointer :: tile
 
  integer, allocatable :: idx(:)   ! integer compressed index of tiles
  integer :: i,j,k,c,n,ntiles,max_cohorts,p
  integer :: iret
  integer, allocatable :: ncohorts(:) ! array of idx sizes from all PEs in io_domain
  integer, allocatable :: idx2(:) ! array of cohort indices from all PEs in io_domain

  ! count total number of cohorts in compute domain and max number of
  ! of cohorts per tile
  ce = first_elmt(lnd%tile_map)
  te = tail_elmt (lnd%tile_map)
  n  = 0
  max_cohorts = 0
  do while (ce/=te)
     tile=>current_tile(ce)
     if(associated(tile%vegn))then
        n = n+tile%vegn%n_cohorts 
        max_cohorts = max(max_cohorts,tile%vegn%n_cohorts)
     endif
     ce=next_elmt(ce)
  enddo

  call mpp_max(max_cohorts)

  ! get the size of the tile dimension from the file
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))
  
  ! calculate compressed cohort index to be written to the restart file
  allocate(idx(max(n,1))) ; idx(:) = -1
  ce = first_elmt(lnd%tile_map, lnd%is, lnd%js)
  n = 1
  do while (ce/=te)
     tile=>current_tile(ce)
     if(associated(tile%vegn)) then
        call get_elmt_indices(ce,i,j,k)
        do c = 1,tile%vegn%n_cohorts
           idx (n) = &
                (c-1)*lnd%nlon*lnd%nlat*ntiles + &
                (k-1)*lnd%nlon*lnd%nlat + &
                (j-1)*lnd%nlon + &
                (i-1)        
           n = n+1
        enddo
     endif
     ce=next_elmt(ce)
  end do

  if (mpp_pe()/=lnd%io_pelist(1)) then
     ! if this processor is not doing io (that is, it's not root io_domain
     ! processor), simply send the data to the root io_domain PE
     call mpp_send(size(idx), plen=1,         to_pe=lnd%io_pelist(1), tag=COMM_TAG_1)
     call mpp_send(idx(1),    plen=size(idx), to_pe=lnd%io_pelist(1), tag=COMM_TAG_2)
  else
     ! gather the array of cohort index sizes
     allocate(ncohorts(size(lnd%io_pelist)))
     ncohorts(1) = size(idx)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(ncohorts(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_1)
     enddo
     ! gather cohort index from the processors in our io_domain
     allocate(idx2(sum(ncohorts(:))))
     idx2(1:ncohorts(1))=idx(:)
     k=ncohorts(1)+1
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(idx2(k), from_pe=lnd%io_pelist(p), glen=ncohorts(p), tag=COMM_TAG_2)
        k = k+ncohorts(p)
     enddo
     ! create cohort dimension in the output file
     iret = nf_redef(ncid)
     __NF_ASRT__(nfu_def_dim(ncid,'cohort',max_cohorts))
     ! create cohort index
     __NF_ASRT__(nfu_def_dim(ncid,cohort_index_name,idx2,'compressed vegetation cohort index'))
     __NF_ASRT__(nfu_put_att(ncid,cohort_index_name,'compress','cohort tile lat lon'))
     __NF_ASRT__(nfu_put_att(ncid,cohort_index_name,'valid_min',0))
     ! deallocate the data we no longer need
     deallocate(ncohorts,idx2)
     ! leave the define mode to commit the new definitions to the disk
     iret = nf_enddef(ncid)
  endif
  call sync_nc_files(ncid)
end subroutine create_cohort_dimension


! ============================================================================
! creates cohort dimension, if necessary, in the output restart file. NOTE 
! that this subroutine should be called even if restart has not been created
! (because, for example, there happen to be no vegetation in a certain domain),
! for the reason that it calls mpp_max, and that should be called for each
! processor to work.
subroutine create_cohort_out_file (ncid, name, land, tile_exists, &
     tile_dim_length, iotiles, num_cohorts, zaxis_data, created)
  integer,               intent(out) :: ncid      ! resulting NetCDF id
  character(len=*),      intent(in)  :: name      ! name of the file to create
  type(land_state_type), intent(in)  :: land
  integer,               intent(in)  :: tile_dim_length ! length of tile axis
  integer,               intent(out) :: iotiles   ! number of tiles on I/O domain
  integer,               intent(out) :: num_cohorts  ! number of cohorts on I/O domain
  real,        optional, intent(in)  :: zaxis_data(:) ! data for the Z-axis
  logical,     optional, intent(out) :: created   ! indicates wether the file was 
      ! created; it is set to false if no restart needs to be written, in case 
      ! the total number of qualifying tiles in this domain is equal to zero
  ! the following interface describes the "detector function", which is passed 
  ! through the argument list and must return true for any tile to be written 
  ! to the specific restart, false otherwise
  interface
     logical function tile_exists(tile)
        use land_tile_mod, only : land_tile_type
        type(land_tile_type), pointer :: tile
     end function tile_exists
  end interface

  ! ---- local vars
  type(land_tile_enum_type) :: ce, te ! tile list enumerators
  type(land_tile_type), pointer :: tile
 
  integer, allocatable :: cohort_idx(:)   ! integer compressed index of tiles
  integer :: i,j,k,c,n,max_cohorts,p
  integer :: iret
  integer, allocatable :: ncohorts(:) ! array of cohort_idx sizes from all PEs in io_domain
  integer, allocatable :: cohort_idx2(:) ! array of cohort indices from all PEs in io_domain
  type(axistype) :: tile_axis(1)
  integer, allocatable :: idx(:)   ! integer compressed index of tiles

  ! First gather the cohort index data

  ! count total number of cohorts in compute domain and max number of
  ! of cohorts per tile
  ce = first_elmt(land%tile_map)
  te = tail_elmt (land%tile_map)
  n  = 0
  max_cohorts = 0
  do while (ce/=te)
     tile=>current_tile(ce)
     if(associated(tile%vegn))then
        n = n+tile%vegn%n_cohorts 
        max_cohorts = max(max_cohorts,tile%vegn%n_cohorts)
     endif
     ce=next_elmt(ce)
  enddo

  call mpp_max(max_cohorts)

  ! get the size of the tile dimension from the file
!  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))
  
  ! calculate compressed cohort index to be written to the restart file
  allocate(cohort_idx(max(n,1))) ; cohort_idx(:) = -1
  ce = first_elmt(land%tile_map, land%is, land%js)
  n = 1
  do while (ce/=te)
     tile=>current_tile(ce)
     if(associated(tile%vegn)) then
        call get_elmt_indices(ce,i,j,k)
        do c = 1,tile%vegn%n_cohorts
           cohort_idx (n) = &
                (c-1)*land%nlon*land%nlat*tile_dim_length + &
                (k-1)*land%nlon*land%nlat + &
                (j-1)*land%nlon + &
                (i-1)        
           n = n+1
        enddo
     endif
     ce=next_elmt(ce)
  end do

  if (mpp_pe()/=land%io_pelist(1)) then
     ! if this processor is not doing io (that is, it's not root io_domain
     ! processor), simply send the data to the root io_domain PE
     call mpp_send(size(cohort_idx), plen=1,         to_pe=land%io_pelist(1), tag=COMM_TAG_1)
!     call mpp_send(cohort_idx(1),    plen=size(cohort_idx), to_pe=land%io_pelist(1), tag=COMM_TAG_2)
     call mpp_recv(num_cohorts,    glen=1, from_pe=land%io_pelist(1), tag=COMM_TAG_2)
  else
     ! gather the array of cohort index sizes
     allocate(ncohorts(size(land%io_pelist)))
     ncohorts(1) = size(cohort_idx)
     do p = 2,size(land%io_pelist)
        call mpp_recv(ncohorts(p), from_pe=land%io_pelist(p), glen=1, tag=COMM_TAG_1)
     enddo
    ! gather cohort index from the processors in our io_domain
!     allocate(cohort_idx2(sum(ncohorts(:))))
     num_cohorts = sum(ncohorts(:))
!     cohort_idx2(1:ncohorts(1))=cohort_idx(:)
!     k=ncohorts(1)+1
     do p = 2,size(land%io_pelist)
        call mpp_send(num_cohorts, plen=1, to_pe=land%io_pelist(p), tag=COMM_TAG_2)
!        call mpp_recv(cohort_idx2(k), from_pe=land%io_pelist(p), glen=ncohorts(p), tag=COMM_TAG_2)
!        k = k+ncohorts(p)
     enddo
     ! deallocate the data we no longer need
     deallocate(ncohorts)!,cohort_idx2)
  endif

  call create_tile_out_file (ncid, name, land, tile_exists, &
     tile_dim_length, iotiles, zaxis_data, max_cohorts, cohort_idx, created)

end subroutine create_cohort_out_file


! ============================================================================
subroutine read_cohort_data_r0d_fptr_new(filename,name,fptr,rec)
!ncid,ntiles,idx,field,fptr,rec)
  character(*), intent(in) :: filename
  character(*), intent(in) :: name
!ncid,ntiles,idx,field,fptr,rec)
!  integer              , intent(in) :: ncid ! netcdf id
!  integer              , intent(in) :: ntiles ! length of "tile" axis
!  integer, dimension(:), intent(in) :: idx ! tile_index axis data
!  type(fieldtype)      , intent(in) :: field ! field containing information about the variable to read
  integer, optional , intent(in) :: rec  ! record number (in case there are 
                                         ! several in the file) 
  ! subroutine returning the pointer to the data to be written
  interface ; subroutine fptr(cohort, ptr)
     use vegn_cohort_mod, only : vegn_cohort_type
     type(vegn_cohort_type), pointer :: cohort ! input
     real, pointer :: ptr ! returned pointer to the data
   end subroutine fptr
  end interface

  ! ---- local constants
  character(*), parameter :: module_name = 'read_cohort_data_r0d_fptr_new'

  ! ---- local vars
  integer :: i, n
  integer :: rec_     ! record number
  integer, allocatable :: idx(:) ! index data to be read
  real,    allocatable :: x1d(:) ! data to be read
  real,    pointer     :: ptr ! pointer to the individual cohort data
  type(vegn_cohort_type), pointer :: cohort
  integer :: siz(4), ntiles
  logical :: found

  ! assign the internal record number
  if(present(rec)) then
     rec_ = rec
  else
     rec_ = 1
  endif


  ! read the input 
  found = .false.
  call field_size(filename,'tile',siz, field_found=found, domain=lnd%domain)
  if ( .not. found ) &
    call error_mesg(trim(module_name),'tile axis not found in '//trim(filename), FATAL)
  ntiles=siz(1) ! 'tile' is an axis variable but does not have any data.

! Find and read the cohort axis data
  call field_size(filename,'cohort_index',siz, field_found=found, domain=lnd%domain)
  if ( .not. found ) &
    call error_mesg(trim(module_name),'cohort_index axis not found in '//trim(filename), FATAL)
  allocate(idx(siz(1)))
  call read_data(filename,'cohort_index',idx, domain=lnd%domain, is_compressed=.true.)

  allocate(x1d(siz(1)))
  call read_data(filename,name,x1d, domain=lnd%domain, is_compressed=.true.)
  ! distribute data over cohorts
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), lnd%nlon, lnd%nlat, ntiles,&
          lnd%tile_map, lnd%is, lnd%js,cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) ptr = x1d(i)
     endif
  enddo
  
end subroutine read_cohort_data_r0d_fptr_new


! ============================================================================
subroutine read_cohort_data_i0d_fptr_new(filename,name,fptr,rec)
!ncid,ntiles,idx,field,fptr,rec)
  character(*), intent(in) :: filename
  character(*), intent(in) :: name
!  integer              , intent(in) :: ncid ! netcdf id
!  integer              , intent(in) :: ntiles ! length of "tile" axis
!  integer, dimension(:), intent(in) :: idx ! tile_index axis data
!  type(fieldtype)      , intent(in) :: field ! field containing information about the variable to read
  integer, optional , intent(in) :: rec  ! record number (in case there are 
                                         ! several in the file) 
  ! subroutine returning the pointer to the data to be written
  interface ; subroutine fptr(cohort, ptr)
     use vegn_cohort_mod, only : vegn_cohort_type
     type(vegn_cohort_type), pointer :: cohort ! input
     integer, pointer :: ptr ! returned pointer to the data
   end subroutine fptr
  end interface

  ! ---- local constants
  character(*), parameter :: module_name = 'read_cohort_data_r0d_fptr_new'

  ! ---- local vars
  integer :: i, n
  integer :: rec_     ! record number
  integer, allocatable :: idx(:) ! index data to be read
  real, allocatable :: x1d(:) ! data to be read
  integer, pointer  :: ptr ! pointer to the individual cohort data
  type(vegn_cohort_type), pointer :: cohort

  integer :: siz(4), ntiles
  logical :: found

  ! assign the internal record number
  if(present(rec)) then
     rec_ = rec
  else
     rec_ = 1
  endif

  ! read the input 
  found = .false.
  call field_size(filename,'tile',siz, field_found=found, domain=lnd%domain)
  if ( .not. found ) &
    call error_mesg(trim(module_name),'tile axis not found in '//trim(filename), FATAL)
  ntiles=siz(1) ! 'tile' is an axis variable but does not have any data.

! Find and read the cohort axis data
  call field_size(filename,'cohort_index',siz, field_found=found, domain=lnd%domain)
  if ( .not. found ) &
    call error_mesg(trim(module_name),'cohort_index axis not found in '//trim(filename), FATAL)
  allocate(idx(siz(1)))
  call read_data(filename,'cohort_index',idx, domain=lnd%domain)

  allocate(x1d(siz(1)))
  call read_data(filename,name,x1d, domain=lnd%domain, is_compressed=.true.)
  ! distribute data over cohorts
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), lnd%nlon, lnd%nlat, ntiles,&
          lnd%tile_map, lnd%is, lnd%js,cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) ptr = x1d(i)
     endif
  enddo
  
end subroutine read_cohort_data_i0d_fptr_new

! ============================================================================
subroutine write_cohort_data_r0d_fptr_new(ncid, ntiles, ncohorts, field, fptr)
  integer         , intent(in) :: ncid ! netcdf id
  integer         , intent(in) :: ntiles ! size of the tile dimension in the output file
  integer         , intent(in) :: ncohorts ! size of the cohort index dimension in the output file
  type(fieldtype) , intent(in) :: field
  ! subroutine returning the pointer to the data to be written
  interface ; subroutine fptr(cohort, ptr)
       use vegn_cohort_mod, only : vegn_cohort_type
       type(vegn_cohort_type), pointer :: cohort ! input
       real, pointer :: ptr ! returned pointer to the data
     end subroutine fptr
  end interface

  ! ---- local vars
  integer :: i, varid, record_, p
  integer,  allocatable :: idx(:) ! index dimension
  real,    allocatable :: idx_r(:)  ! index dimension converted to real
  real, allocatable :: data(:) ! data to be written
  real, allocatable :: buffer(:) ! input buffer for data from other PEs
  integer,  allocatable :: mask(:) ! mask of the valid data
  real, pointer :: ptr ! pointer to the individual cohort data
  type(vegn_cohort_type), pointer :: cohort 
  integer :: dimids(2), ndims
  type(axistype) :: cohort_axis(1)

  ! get the length of cohort compressed index
!  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))

  ! get the length of tile dimension
!  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! allocate data
  allocate(data(ncohorts),idx(ncohorts),idx_r(ncohorts),mask(ncohorts))
  data = NF_FILL_DOUBLE
  mask = 0

  ! read cohort index
!  i = nf_enddef(ncid) ! ignore errors (the file may be in data mode already)
!  __NF_ASRT__(nfu_get_var(ncid,cohort_index_name,idx))

  ! read cohort index
  if (mpp_pe()/=lnd%io_pelist(1)) then
     ! if current PE doesn't do io, we need the number of tiles in this io-domain to be received
     call mpp_recv(idx(1), from_pe=lnd%io_pelist(1), glen=ncohorts, tag=COMM_TAG_1)
  else
     ! Send the tile index to all processors in our io_domain
     cohort_axis = mpp_get_axis_by_name(ncid,'cohort_index')
     ! Axis data is stored as double so need to read and convert to integer.
     call mpp_get_axis_data(cohort_axis(1),idx_r)
     idx(1:ncohorts) = int(idx_r(1:ncohorts))
     do p = 2,size(lnd%io_pelist)
        call mpp_send(idx(1), plen=ncohorts, to_pe=lnd%io_pelist(p), tag=COMM_TAG_1)
     enddo
  endif

  ! gather data into an array along the cohort dimension
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), lnd%nlon, lnd%nlat, ntiles,&
                             lnd%tile_map, lnd%is, lnd%js,cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) then 
           data(i) = ptr
           mask(i) = 1
        endif
     endif
  enddo
  
  ! if this processor isn't the root IO processor, simply send data to the root 
  ! IO processor and return from the subroutine
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(data(1), plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_2)
     call mpp_send(mask(1), plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_3)
  else
     ! gather data from the processors in io_domain
     allocate(buffer(size(data)))
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buffer(1), glen=size(data), from_pe=lnd%io_pelist(p), tag=COMM_TAG_2)
        call mpp_recv(mask(1),   glen=size(data), from_pe=lnd%io_pelist(p), tag=COMM_TAG_3)
        where(mask > 0) data = buffer
     enddo
     deallocate(buffer,mask)

  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to 
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
  ! write data
  if (mpp_pe() .eq. lnd%io_pelist(1)) then
    call mpp_write(ncid,field,data)
  endif

  ! free allocated memory
  deallocate(data,idx)
  
end subroutine write_cohort_data_r0d_fptr_new

! ============================================================================
subroutine write_cohort_data_i0d_fptr_new(ncid, ntiles, ncohorts, field, fptr)
  integer         , intent(in) :: ncid ! netcdf id
  integer         , intent(in) :: ntiles ! size of the tile dimension in the output file
  integer         , intent(in) :: ncohorts ! size of the cohort index dimension in the output file
  type(fieldtype) , intent(in) :: field
  ! subroutine returning the pointer to the data to be written
  interface ; subroutine fptr(cohort, ptr)
       use vegn_cohort_mod, only : vegn_cohort_type
       type(vegn_cohort_type), pointer :: cohort ! input
       integer, pointer :: ptr ! returned pointer to the data
     end subroutine fptr
  end interface

  ! ---- local vars
  integer :: i, varid, record_, p
  integer,  allocatable :: idx(:) ! index dimension
  real,    allocatable :: idx_r(:)  ! index dimension converted to real
  integer, allocatable :: data(:) ! data to be written
  integer, allocatable :: buffer(:) ! input buffer for data from other PEs
  integer,  allocatable :: mask(:) ! mask of the valid data
  integer, pointer :: ptr ! pointer to the individual cohort data
  type(vegn_cohort_type), pointer :: cohort 
  integer :: dimids(2), ndims
  type(axistype) :: cohort_axis(1)

  ! allocate data
  allocate(data(ncohorts),idx(ncohorts),idx_r(ncohorts),mask(ncohorts))
  data = NF_FILL_INT
  mask = 0

  ! read cohort index
  if (mpp_pe()/=lnd%io_pelist(1)) then
     ! if current PE doesn't do io, we need the number of tiles in this io-domain to be received
     call mpp_recv(idx(1), from_pe=lnd%io_pelist(1), glen=ncohorts, tag=COMM_TAG_1)
  else
     ! Send the tile index to all processors in our io_domain
     cohort_axis = mpp_get_axis_by_name(ncid,'cohort_index')
     ! Axis data is stored as double so need to read and convert to integer.
     call mpp_get_axis_data(cohort_axis(1),idx_r)
     idx(1:ncohorts) = int(idx_r(1:ncohorts))
     do p = 2,size(lnd%io_pelist)
        call mpp_send(idx(1), plen=ncohorts, to_pe=lnd%io_pelist(p), tag=COMM_TAG_1)
     enddo
  endif

  ! gather data into an array along the cohort dimension
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), lnd%nlon, lnd%nlat, ntiles,&
                             lnd%tile_map, lnd%is, lnd%js,cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) then 
           data(i) = ptr
           mask(i) = 1
        endif
     endif
  enddo
  
  ! if this processor isn't the root IO processor, simply send data to the root 
  ! IO processor and return from the subroutine
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(data(1), plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_2)
     call mpp_send(mask(1), plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_3)
  else
     ! gather data from the processors in io_domain
     allocate(buffer(size(data)))
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buffer(1), glen=size(data), from_pe=lnd%io_pelist(p), tag=COMM_TAG_2)
        call mpp_recv(mask(1),   glen=size(data), from_pe=lnd%io_pelist(p), tag=COMM_TAG_3)
        where(mask > 0) data = buffer
     enddo
     deallocate(buffer,mask)

  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to 
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
  ! write data
  if (mpp_pe() .eq. lnd%io_pelist(1)) then
    call mpp_write(ncid,field,real(data))
  endif

  ! free allocated memory
  deallocate(data,idx)
  
end subroutine write_cohort_data_i0d_fptr_new


#define F90_TYPE real
#define NF_TYPE NF_DOUBLE
#define NF_FILL_VALUE NF_FILL_DOUBLE
#define READ_0D_FPTR read_cohort_data_r0d_fptr
#define WRITE_0D_FPTR write_cohort_data_r0d_fptr
#include "vegn_cohort_io.inc"

#define F90_TYPE integer
#define NF_TYPE NF_INT
#define NF_FILL_VALUE NF_FILL_INT
#define READ_0D_FPTR read_cohort_data_i0d_fptr
#define WRITE_0D_FPTR write_cohort_data_i0d_fptr
#include "vegn_cohort_io.inc"

end module cohort_io_mod
