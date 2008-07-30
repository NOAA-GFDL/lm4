module cohort_io_mod

use fms_mod,          only : error_mesg, FATAL, WARNING
use mpp_mod,          only : mpp_max

use nf_utils_mod,     only : nfu_inq_dim, nfu_get_var, &
     nfu_get_rec, nfu_def_dim, nfu_def_var, nfu_put_att_text
use land_io_mod,      only : print_netcdf_error
use land_tile_mod,    only : land_tile_type, land_tile_list_type, &
     land_tile_enum_type, first_elmt, tail_elmt, next_elmt, get_elmt_indices, &
     current_tile, operator(/=)
use land_tile_io_mod, only : write_tile_data, get_tile_by_idx 

use vegn_cohort_mod, only: vegn_cohort_type
use land_data_mod, only : lnd

implicit none
private

! ==== public interfaces =====================================================
! input
public :: read_create_cohorts
public :: read_cohort_data_r0d_fptr
public :: read_cohort_data_i0d_fptr
! output
public :: create_cohort_dimension
public :: write_cohort_data_r0d_fptr
public :: write_cohort_data_i0d_fptr
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'cohort_io_mod', &
     version     = '$Id: vegn_cohort_io.F90,v 16.0 2008/07/30 22:30:10 fms Exp $', &
     tagname     = '$Name: perth $'
! name of the "compressed" dimension (and dimension variable) in the output 
! netcdf files -- that is, the dimensions written out using compression by 
! gathering, as described in CF conventions.
character(len=*),   parameter :: cohort_index_name   = 'cohort_index'

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
  integer :: i,j,t,k,m, n
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  character(len=64) :: info ! for error message

  ! get the size of dimensions
  nlon = size(lnd%garea,1) ; nlat = size(lnd%garea,2)
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! read the cohort index
  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))
  allocate(idx(ncohorts))
  __NF_ASRT__(nfu_get_var(ncid,cohort_index_name,idx))
  
  do n = 1,size(idx)
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
subroutine read_cohort_data_r0d_fptr(ncid,name,fptr,rec)
  integer           , intent(in) :: ncid ! netcdf id
  character(len=*)  , intent(in) :: name ! name of the variable to read
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
  character(*), parameter :: module_name = 'read_cohort_data_r0d_fptr'

  ! ---- local vars
  integer :: i
  integer :: rec_     ! record number
  integer :: ntiles   ! size of the tile dimension in restart file
  integer :: ncohorts ! total number of cohorts in restart file
  integer, allocatable :: idx(:)    ! index dimension
  real   , allocatable :: data(:)   ! data to be read
  real   , pointer :: ptr ! pointer to the individual cohort data
  type(vegn_cohort_type), pointer :: cohort

  ! assign the internal record number
  if(present(rec)) then
     rec_ = rec
  else
     rec_ = 1
  endif

  ! get the size of the tile dimension
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! get the length of cohort compressed index
  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))

  ! allocate data
  allocate(data(ncohorts),idx(ncohorts))

  ! read the cohort index
  __NF_ASRT__(nfu_get_var(ncid,cohort_index_name,idx))
  ! read the data
  __NF_ASRT__(nfu_get_rec(ncid,name,rec_,data))

  ! distribute data over cohorts
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), size(lnd%garea,1), size(lnd%garea,2), ntiles,&
                             lnd%tile_map, lnd%is, lnd%js,cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) ptr = data(i)
     endif
  enddo
  
  ! free allocated memory
  deallocate(data,idx)
  
end subroutine read_cohort_data_r0d_fptr
  

! ============================================================================
subroutine read_cohort_data_i0d_fptr(ncid,name,fptr,rec)
  integer           , intent(in) :: ncid ! netcdf id
  character(len=*)  , intent(in) :: name ! name of the variable to read
  integer, optional , intent(in) :: rec  ! record number (in case there are 
                                         ! several in the file) 
  ! subroutine returning the pointer to the data to be written
  interface
     subroutine fptr(cohort, ptr)
       use vegn_cohort_mod, only : vegn_cohort_type
       type(vegn_cohort_type), pointer :: cohort ! input
       integer          , pointer :: ptr    ! returned pointer to the data
     end subroutine fptr
  end interface

  ! ---- local constants
  character(*), parameter :: module_name = 'read_cohort_data_i0d_fptr'

  ! ---- local vars
  integer :: i
  integer :: rec_     ! record number
  integer :: ntiles   ! size of the tile dimension in restart file
  integer :: ncohorts ! total number of cohorts in restart file
  integer, allocatable :: idx(:)    ! index dimension
  integer, allocatable :: data(:)   ! data to be read
  integer, pointer :: ptr ! pointer to the individual cohort data
  type(vegn_cohort_type), pointer :: cohort

  ! assign the internal record number
  if(present(rec)) then
     rec_ = rec
  else
     rec_ = 1
  endif

  ! get the size of dimensions
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! get the length of cohort compressed index
  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))

  ! allocate data
  allocate(data(ncohorts),idx(ncohorts))

  ! read the cohort index
  __NF_ASRT__(nfu_get_var(ncid,cohort_index_name,idx))
  ! read the data
  __NF_ASRT__(nfu_get_rec(ncid,name,rec_,data))

  ! distribute data over cohorts
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), size(lnd%garea,1), size(lnd%garea,2), ntiles,&
                             lnd%tile_map, lnd%is, lnd%js,cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) ptr = data(i)
     endif
  enddo
  
  ! free allocated memory
  deallocate(data,idx)
  
end subroutine read_cohort_data_i0d_fptr


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
  integer :: i,j,k,c,n,ntiles,max_cohorts
  integer :: iret

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

  if (n>0) then
     ! get the size of the tile dimension from the file
     __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

     ! calculate compressed cohort index to be written to the restart file
     allocate(idx(n))
     ce = first_elmt(lnd%tile_map, lnd%is, lnd%js)
     n = 1
     do while (ce/=te)
        tile=>current_tile(ce)
        if(associated(tile%vegn)) then
           call get_elmt_indices(ce,i,j,k)
           do c = 1,tile%vegn%n_cohorts
              idx (n) = &
                   (c-1)*size(lnd%garea,1)*size(lnd%garea,2)*ntiles + &
                   (k-1)*size(lnd%garea,1)*size(lnd%garea,2) + &
                   (j-1)*size(lnd%garea,1) + &
                   (i-1)        
              n = n+1
           enddo
        endif
        ce=next_elmt(ce)
     end do
     ! create tile output file, defining horizontal coordinate and compressed
     ! dimension
     iret = nf_redef(ncid)
     __NF_ASRT__(nfu_def_dim(ncid,'cohort',max_cohorts))
     __NF_ASRT__(nfu_def_dim(ncid,cohort_index_name,idx,'compressed vegetation cohort index'))
     __NF_ASRT__(nfu_put_att_text(ncid,cohort_index_name,'compress','cohort tile lat lon'))
  endif
end subroutine create_cohort_dimension


! ============================================================================
subroutine write_cohort_data_r0d_fptr(ncid,name,fptr,long_name,units)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  character(len=*), intent(in), optional :: units, long_name
  ! subroutine returning the pointer to the data to be written
  interface
     subroutine fptr(cohort, ptr)
       use vegn_cohort_mod, only : vegn_cohort_type
       type(vegn_cohort_type), pointer :: cohort ! input
       real             , pointer :: ptr    ! returned pointer to the data
     end subroutine fptr
  end interface

  ! ---- local vars
  integer :: i, varid
  integer :: ntiles   ! size of the tile dimension in the output file
  integer :: ncohorts ! size of the cohort index dimension in the output file
  integer, allocatable :: idx(:)    ! index dimension
  real   , allocatable :: data(:)   ! data to be written
  real, pointer :: ptr ! pointer to the individual cohort data
  type(vegn_cohort_type), pointer :: cohort 

  ! get the length of cohort compressed index
  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))

  ! get the length of tile dimension
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! allocate data
  allocate(data(ncohorts),idx(ncohorts))

  ! read cohort index
  i = nf_enddef(ncid) ! ignore errors (the file may be in data mode already)
  __NF_ASRT__(nfu_get_var(ncid,cohort_index_name,idx))

  ! gather data into an array along the cohort dimension
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), size(lnd%garea,1), size(lnd%garea,2), ntiles,&
                             lnd%tile_map, lnd%is, lnd%js,cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) data(i) = ptr
     endif
  enddo
  
  ! create variable, if it does not exist
  if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
     __NF_ASRT__(nfu_def_var(ncid,name,NF_DOUBLE,(/cohort_index_name/),long_name,units,varid))
  endif
  ! write data
  i = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
  __NF_ASRT__(nf_put_var_double(ncid,varid,data))
  
  ! free allocated memory
  deallocate(data,idx)
  
end subroutine write_cohort_data_r0d_fptr


! ============================================================================
subroutine write_cohort_data_i0d_fptr(ncid,name,fptr,long_name,units)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  character(len=*), intent(in), optional :: units, long_name
  ! subroutine returning the pointer to the data to be written
  interface
     subroutine fptr(cohort, ptr)
       use vegn_cohort_mod, only : vegn_cohort_type
       type(vegn_cohort_type), pointer :: cohort ! input
       integer          , pointer :: ptr    ! returned pointer to the data
     end subroutine fptr
  end interface

  ! ---- local vars
  integer :: i, varid
  integer :: ntiles   ! size of the tile dimension in the output file
  integer :: ncohorts ! size of the cohort dimension in the output file
  integer, allocatable :: idx(:)    ! index dimension
  integer, allocatable :: data(:)   ! data to be written
  integer, pointer :: ptr ! pointer to the individual cohort data
  type(vegn_cohort_type), pointer :: cohort 

  ! get the length of cohort compressed index
  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))

  ! get the length of tile dimension
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! allocate data
  allocate(data(ncohorts),idx(ncohorts))

  ! read cohort index
  i = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
  __NF_ASRT__(nfu_get_var(ncid,cohort_index_name,idx))

  ! gather data into an array along the cohort dimension
  do i = 1, size(idx)
     call get_cohort_by_idx ( idx(i), size(lnd%garea,1), size(lnd%garea,2), ntiles,&
                             lnd%tile_map, lnd%is, lnd%js,cohort)
     if (associated(cohort)) then
        call fptr(cohort, ptr)
        if(associated(ptr)) data(i) = ptr
     endif
  enddo
  
  ! create variable, if it does not exist
  if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
     __NF_ASRT__(nfu_def_var(ncid,name,NF_INT,(/cohort_index_name/),long_name,units,varid))
  endif
  ! write data
  i = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
  __NF_ASRT__(nf_put_var_int(ncid,varid,data))
  
  ! free allocated memory
  deallocate(data,idx)
  
end subroutine write_cohort_data_i0d_fptr

end module cohort_io_mod
