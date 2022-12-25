! -*-f90-*-
! $Id$

! some sanity checks
#ifndef F90_TYPE
#error F90_TYPE is not defined: must be one of FORTRAN 90 types
#endif

#ifndef READ_REMAP_SUB
#error name of READ_REMAP_SUB is not defined
#endif

! ============================================================================
subroutine READ_REMAP_SUB(fileobj, name, fptr, map_i, map_j, cidx, compressed_data)
  class(FmsNetcdfFile_t), intent(in) :: fileobj
  character(len=*), intent(in) :: name
  integer        , intent(in) :: map_i(lnd%ls:) ! re-mapping index
  integer        , intent(in) :: map_j(lnd%ls:) ! re-mapping index
  integer        , intent(in) :: cidx(:)
  F90_TYPE       , intent(in) :: compressed_data(:)

  ! subroutine returning the pointer to the data to be written
  interface
     subroutine fptr(cohort, ptr)
       use vegn_cohort_mod, only : vegn_cohort_type
       type(vegn_cohort_type), pointer :: cohort ! input
       F90_TYPE              , pointer :: ptr    ! returned pointer to the data
     end subroutine fptr
  end interface

  ! ---- local constants
  character(*), parameter :: module_name = "read_remap_cohort_data"

  ! ---- local vars
  integer :: i,j,k,n,ii,jj,ndims, t, c, npos, l
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type)   , pointer :: tile
  type(vegn_cohort_type) , pointer :: cohort
  F90_TYPE, pointer :: ptr ! pointer to the individual cohort data
  F90_TYPE, allocatable :: expanded_data(:,:,:,:) ! buffer for input data
  logical,  allocatable :: mask(:,:,:,:) ! validity mask for input data
  logical :: is_compressed
  integer, dimension(4) :: dimlens
  character(len=256) :: dim_name
  character(len=256) :: compress_att, string

  ! get the size of dimensions
  is_compressed = .false.
  if (variable_exists(fileobj, "cohort_index")) then
    if (variable_att_exists(fileobj, "cohort_index", "compress")) then
      call get_variable_attribute(fileobj, "cohort_index", "compress", compress_att)
      if (trim(compress_att) .ne. "unspecified") then
        is_compressed = .true.
      endif
    endif
  endif

  if(.not.is_compressed) then
    call error_mesg(module_name, &
    'compress attribute not found for cohort_index. Therefore, do not know how to decompress '//trim(name)//' (pjp)',FATAL)
  endif

  ! Get size of each dimension specified by compress_att
  string = compress_att
  do n=1,4
     npos = scan(string, ' ')
     dim_name = string(1:npos-1)
     call get_dimension_size(fileobj, dim_name, dimlens(n))
     string = string(npos+1:len_trim(string))
     npos = verify(string, ' ')
     if(npos == 0) exit
  enddo

  allocate(expanded_data(dimlens(4),dimlens(3),dimlens(2),dimlens(1)))
  allocate(         mask(dimlens(4),dimlens(3),dimlens(2),dimlens(1)))
  expanded_data = 0.0
  mask = .FALSE.
  do n=1,size(cidx)
     k = cidx(n)
     i = modulo(k,dimlens(4))+1 ; k = k/dimlens(4)
     j = modulo(k,dimlens(3))+1 ; k = k/dimlens(3)
     t = modulo(k,dimlens(2))+1 ; k = k/dimlens(2)
     c = k+1
     expanded_data(i,j,t,c) = compressed_data(n)
     mask(i,j,t,c) = .TRUE.
  enddo

  ! distribute data over cohorts. NOTE that this is slightly different from the restart
  ! reading procedure. On reading the restart, all the tiles are counted in sequence,
  ! while here only the vegetation tiles.
  do l = lnd%ls, lnd%le
     ii = map_i(l); jj = map_j(l)
     if ((ii.le.0).or.(jj.le.0)) cycle ! skip un-mapped points
     if (.not.any(mask(ii,jj,:,:))) cycle ! skip points where there is no data

     ce = first_elmt (land_tile_map(l))
     te = tail_elmt  (land_tile_map(l))
     k = 1
tile_loop:  do while(ce/=te.and.k<=dimlens(2))
        tile=>current_tile(ce); ce=next_elmt(ce);
        if (.not.associated(tile%vegn)) cycle
        ! find index of the next valid tile in the input data
        do while(.not.any(mask(ii,jj,k,:)))
           k=k+1 ! go to the next tile if there's no data (i.e. all mask
                 ! values are false for this tile)
           if(k>dimlens(2)) exit tile_loop
        enddo

        do n = 1,min(size(tile%vegn%cohorts(:)),dimlens(1))
           cohort=>tile%vegn%cohorts(n)
           call fptr(cohort,ptr)
           if(associated(ptr).and.mask(ii,jj,k,n)) ptr = expanded_data(ii,jj,k,n)
        enddo
        k = k+1 ! go to the next tile in input data
     enddo tile_loop
  enddo

  ! free allocated memory
  deallocate(expanded_data,mask)

end subroutine
