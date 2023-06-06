module land_tile_diag_mod

use mpp_mod,            only : mpp_sum
use mpp_efp_mod,        only : mpp_reproducing_sum
use time_manager_mod,   only : time_type
use diag_manager_mod,   only : send_data
use fms_mod,            only : error_mesg, string, FATAL, NOTE

use land_tile_selectors_mod, only : tile_selector_type, n_selectors, selectors
use land_tile_mod,      only : land_tile_type, &
     land_tile_enum_type, first_elmt, loop_over_tiles, &
     land_tile_map, tile_is_selected, fptr_i0, fptr_r0, fptr_r0i
use land_data_mod,      only : lnd, log_version, land_data_type
use land_debug_mod,     only : check_var_range, set_current_point
use tile_diag_buff_mod, only : diag_buff_type
use tile_diag_base_mod, only : & ! use everything except send_tile_data interface, which we extend in this module
        BASE_TILED_FIELD_ID, BASE_COHORT_FIELD_ID, &
        OP_AVERAGE, OP_SUM, OP_MAX, OP_MIN, OP_VAR, OP_STD, &
        cmor_name, cmor_mrsos_depth, &
        n_fields, fields, &
        tiled_diag_field_type, &
        tile_diag_base_init, tile_diag_base_end, set_default_diag_filter, &
        register_tiled_area_fields, register_tiled_diag_field, &
        register_tiled_static_field, register_cohort_diag_field, &
        add_tiled_diag_field_alias, add_tiled_static_field_alias, &
        send_tile_data_0d, send_tile_data_1d, send_cohort_data, &
        get_area_id, get_field_id

implicit none; private


! ==== public interface ======================================================
public :: tile_diag_init
public :: tile_diag_end

public :: set_default_diag_filter ! set the default filter for the consequent
          ! register_tiled_diag_field calls. The default filter will handle the
          ! fields that appear in diag_table without subsampling suffix.

public :: diag_buff_type

public :: register_tiled_area_fields
public :: register_tiled_diag_field
public :: register_tiled_static_field
public :: add_tiled_diag_field_alias
public :: add_tiled_static_field_alias

public :: send_tile_data
public :: send_tile_data_r0d_fptr, send_tile_data_r1d_fptr
public :: send_tile_data_i0d_fptr

public :: register_cohort_diag_field
public :: send_cohort_data

public :: dump_tile_diag_fields

public :: get_area_id, get_field_id

public :: send_global_land_diag ! unused?

public :: OP_AVERAGE, OP_SUM, OP_MAX, OP_MIN, OP_VAR, OP_STD
public :: cmor_name, cmor_mrsos_depth

interface send_tile_data
   module procedure send_tile_data_0d
   module procedure send_tile_data_1d
   module procedure send_tile_data_0d_array
end interface
! ==== end of public interface ===============================================


! ==== module constants ======================================================
character(len=*), parameter :: mod_name = 'land_tile_diag_mod'
#include "../shared/version_variable.inc"


! ==== module data ===========================================================
logical :: module_is_initialized = .false.


contains

! ============================================================================
subroutine tile_diag_init()
  if (module_is_initialized) return

  module_is_initialized = .true.
  call log_version(version, mod_name, &
  __FILE__)

  ! initialize diag selectors
  call tile_diag_base_init()
end subroutine

! ============================================================================
subroutine tile_diag_end()
  call tile_diag_base_end()

  module_is_initialized = .false.
end subroutine tile_diag_end

! ============================================================================
subroutine send_tile_data_0d_array(id, x, send_immediately)
  integer, intent(in) :: id
  real   , intent(in) :: x(:,:)
  logical, intent(in), optional :: send_immediately ! if true, send data to diag_manager
    ! right away, instead of waiting for the next tiled diag fields dump

  integer :: l,k
  type(land_tile_enum_type)     :: ce
  type(land_tile_type), pointer :: tileptr

  ce = first_elmt( land_tile_map )
  do while(loop_over_tiles(ce,tileptr,l,k))
     call send_tile_data(id,x(l,k),tileptr%diag)
  enddo
  if(present(send_immediately)) then
     ! TODO: perhaps need to add time to the arguments, instead of using lnd%time?
     ! not clear if this will have any effect
     if(send_immediately) call dump_tile_diag_field(id, lnd%time)
  endif
end subroutine send_tile_data_0d_array

! ============================================================================
subroutine send_tile_data_r0d_fptr(id, fptr)
  integer, intent(in) :: id
  procedure(fptr_r0)  :: fptr

  type(land_tile_enum_type)     :: ce      ! tile list enumerator
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  real                , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( land_tile_map )
  do while(loop_over_tiles(ce,tileptr))
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,ptr,tileptr%diag)
  enddo
end subroutine send_tile_data_r0d_fptr


! ============================================================================
subroutine send_tile_data_r1d_fptr(id, fptr)
  integer, intent(in) :: id
  procedure(fptr_r0i) :: fptr

  type(land_tile_enum_type)     :: ce      ! tile list enumerator
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  real                , pointer :: ptr     ! pointer to the data element within a tile
  real,             allocatable :: buffer(:) ! buffer for accumulating data
  integer :: i, i2,i3, k
  logical :: have_data

  if(id <= 0) return
  if (id < BASE_TILED_FIELD_ID.or. id >= BASE_COHORT_FIELD_ID ) call error_mesg (mod_name, &
         'tile diag field ID ('//string(id)//') is out of range. Perhaps the field was registred with some other call then register_tile_diag_field?', &
         FATAL)
  i = id - BASE_TILED_FIELD_ID ! index in the array of fields

  allocate(buffer(fields(i)%size))
  ce = first_elmt( land_tile_map, ls=lnd%ls )
  do while(loop_over_tiles(ce, tileptr,i2,i3))
#ifdef DEBUG_LAND_TILE_DIAG
     call set_current_point(i2,i3)
#endif
     have_data = .FALSE.
     do k = 1,fields(i)%size
        call fptr(tileptr,k,ptr)
        if(associated(ptr)) then
            buffer(k) = ptr
            have_data = .TRUE.
        else
            buffer(k) = 0.0
        endif
     enddo
     if (have_data) call send_tile_data(id,buffer,tileptr%diag)
  enddo
  deallocate(buffer)
end subroutine send_tile_data_r1d_fptr


! ============================================================================
subroutine send_tile_data_i0d_fptr(id, fptr)
  integer, intent(in) :: id
  procedure(fptr_i0)  :: fptr

  type(land_tile_enum_type)     :: ce      ! tile list enumerator
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  integer             , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( land_tile_map )
  do while(loop_over_tiles(ce,tileptr))
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,real(ptr),tileptr%diag)
  enddo
end subroutine send_tile_data_i0d_fptr


! ============================================================================
subroutine dump_tile_diag_fields(time)
  type(time_type)          , intent(in) :: time       ! current time

  ! ---- local vars
  integer :: ifld ! field number
  integer :: isel ! selector number
  type(land_tile_enum_type)     :: ce
  type(land_tile_type), pointer :: tile
  integer :: total_n_sends(n_fields)

  total_n_sends(:) = fields(1:n_fields)%n_sends
  call mpp_sum(total_n_sends, n_fields, pelist=lnd%pelist)

!$OMP parallel do default(none) shared(land_tile_map,n_fields,total_n_sends,n_selectors,fields,selectors,time)
  do ifld = 1, n_fields
     if (total_n_sends(ifld) == 0) cycle ! no data to send
     ! write(*,*)trim(fields(ifld)%module),'/',trim(fields(ifld)%name)
     do isel = 1, n_selectors
        if (fields(ifld)%ids(isel) <= 0) cycle
        call dump_diag_field_with_sel (fields(ifld)%ids(isel), &
             fields(ifld), selectors(isel), time )
     enddo
  enddo
  ! zero out the number of data points sent to the field
  fields(1:n_fields)%n_sends=0

  ! all the data are sent to the output, so set the data presence tag to FALSE
  ! in all diag buffers in preparation for the next time step
  ce = first_elmt(land_tile_map)
  do while(loop_over_tiles(ce,tile))
    tile%diag%mask(:) = .FALSE.
  enddo
end subroutine dump_tile_diag_fields

! ============================================================================
! dumps a single field
! TODO: perhaps need dump aliases as well
! TODO: perhaps total_n_sends check can be removed to avoid communication
subroutine dump_tile_diag_field(id, time)
  integer, intent(in) :: id ! diag id of the field
  type(time_type), intent(in) :: time       ! current time

  ! ---- local vars
  integer :: ifld ! field number
  integer :: isel ! selector number
  type(land_tile_enum_type)     :: ce
  type(land_tile_type), pointer :: tile
  integer :: total_n_sends

  if (id<=0) return ! do nothing if field not registered

  ifld = id-BASE_TILED_FIELD_ID
  if (ifld<1.or.ifld>n_fields) &
     call error_mesg(mod_name, 'incorrect field id '//string(id)//' in dump_tile_diag_field ', FATAL)

  total_n_sends = fields(ifld)%n_sends
  call mpp_sum(total_n_sends, pelist=lnd%pelist)

  if (total_n_sends == 0) return ! no data to send
!$OMP parallel do default(none) shared(land_tile_map,n_selectors,fields,ifld,selectors,time) private(isel)
  do isel = 1, n_selectors
     if (fields(ifld)%ids(isel) <= 0) cycle
     call dump_diag_field_with_sel (fields(ifld)%ids(isel), &
          fields(ifld), selectors(isel), time )
  enddo
  ! zero out the number of data points sent to the field
  fields(ifld)%n_sends=0

  ! all the data are sent to the output, so set the data presence tag to FALSE
  ! in all diag buffers in preparation for the next time step
  ce = first_elmt(land_tile_map)
  do while(loop_over_tiles(ce,tile))
    tile%diag%mask(fields(ifld)%offset:fields(ifld)%offset+fields(ifld)%size-1) = .FALSE.
  enddo

end subroutine dump_tile_diag_field

! ============================================================================
subroutine dump_diag_field_with_sel(id, field, sel, time)
  integer                    , intent(in) :: id
  type(tiled_diag_field_type), intent(in) :: field
  type(tile_selector_type)   , intent(in) :: sel
  type(time_type)            , intent(in) :: time ! current time

  ! ---- local vars
  integer :: l ! iterators
  integer :: ks,ke ! array boundaries
  integer :: ls, le
  logical :: used ! value returned from send_data (ignored)
  real, allocatable :: buffer(:,:), weight(:,:), var(:,:)
  logical, allocatable :: mask(:,:)
  type(land_tile_enum_type)     :: ce
  type(land_tile_type), pointer :: tile

  ! calculate array boundaries
  ls = lbound(land_tile_map,1); le = ubound(land_tile_map,1)
  ks = field%offset   ; ke = field%offset + field%size - 1

  ! allocate and initialize temporary buffers
  allocate(buffer(ls:le,ks:ke), weight(ls:le,ks:ke), mask(ls:le,ks:ke))
  weight(:,:) = 0.0
  buffer(:,:) = 0.0

  ! accumulate data
  ce = first_elmt(land_tile_map, ls=ls)
  do while(loop_over_tiles(ce, tile, l))
    if ( size(tile%diag%data) < ke )       cycle ! do nothing if there is no data in the buffer
    if ( .not.tile_is_selected(tile,sel) ) cycle ! do nothing if tile is not selected
    select case (field%opcode)
    case (OP_AVERAGE,OP_VAR,OP_STD)
       where(tile%diag%mask(ks:ke))
          buffer(l,:) = buffer(l,:) + tile%diag%data(ks:ke)*tile%frac
       end where
       weight(l,:) = weight(l,:) + tile%frac
    case (OP_SUM)
       where(tile%diag%mask(ks:ke))
          buffer(l,:) = buffer(l,:) + tile%diag%data(ks:ke)
       end where
       weight(l,:) = 1
    end select
  enddo

  ! normalize accumulated data
  mask = (weight>0)
  where (mask) buffer=buffer/weight

  if (field%opcode == OP_VAR.or.field%opcode == OP_STD) then
     ! second loop to process the variance and standard deviation diagnostics.
     ! it may be possible to calc. var and std in one pass with weighted incremental
     ! algorithm from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
     ! code her is more straightforward. buffer(:,:,:) already contains the mean,
     ! and weight(:,:,:) -- sum of  tile fractions
     allocate(var(ls:le,ks:ke))
     var(:,:) = 0.0
     ! the loop is somewhat different from the first, for no particular reason:
     ! perhaps this way is better for performance?
     do l = ls, le
        ce = first_elmt(land_tile_map(l))
        do while(loop_over_tiles(ce,tile))
           if ( size(tile%diag%data) < ke )       cycle ! do nothing if there is no data in the buffer
           if ( .not.tile_is_selected(tile,sel) ) cycle ! do nothing if tile is not selected
           where(tile%diag%mask(ks:ke))
              var(l,:) = var(l,:) + tile%frac*(tile%diag%data(ks:ke)-buffer(l,:))**2
           end where
        enddo
     enddo
     ! renormalize the variance or standard deviation. note that weight is
     ! calculated in the first loop
     select case (field%opcode)
     case (OP_VAR)
         where (mask) buffer = var/weight
     case (OP_STD)
         where (mask) buffer = sqrt(var/weight)
     end select
     deallocate(var)
  endif

  if (field%fill_missing) then
      where (.not. mask) buffer = 0.0
      mask = .true.
  endif
  ! send diag field
  used = send_data(id,buffer,time,mask=mask)

  ! clean up temporary data
  deallocate(buffer,weight,mask)

end subroutine dump_diag_field_with_sel

! ============================================================================
!> \brief Send out the land model field on unstructured grid for global integral
logical function send_global_land_diag( id, diag, Time, tile, mask, Land )
  integer,                 intent(in) :: id
  real,    dimension(:,:), intent(in) :: diag, tile
  type(time_type),         intent(in) :: Time
  logical, dimension(:,:), intent(in) :: mask
  type(land_data_type),    intent(in) :: Land

  real,    dimension(size(diag,1),1)    :: diag_ug, tile_ug, area_ug
  logical, dimension(size(mask,1))    :: mask_ug
  integer :: k
  real    :: area_sum, diag_sum

  ! sum over tiles on unstructured grid
  diag_ug = 0.0
  tile_ug = 0.0
  do k = 1, size(diag,2)
    where (mask(:,k))
      diag_ug(:,1) = diag_ug(:,1) + diag(:,k)*tile(:,k)
      tile_ug(:,1) = tile_ug(:,1) + tile(:,k)
    endwhere
  enddo
  ! average on unstructured grid
  where (tile_ug > 0.0)
    diag_ug = diag_ug/tile_ug
  endwhere
  mask_ug(:) = ANY(mask,dim=2)

  where(mask_ug)
     diag_ug(:,1) = diag_ug(:,1) * lnd%ug_area
     area_ug(:,1) = lnd%ug_area
  elsewhere
     diag_ug(:,1) = 0.0
     area_ug(:,1) = 0.0
  endwhere

  diag_sum = mpp_reproducing_sum(diag_ug) !, overflow_check=.true.)
  area_sum = mpp_reproducing_sum(area_ug) !, overflow_check=.true.)

  send_global_land_diag = send_data( id, diag_sum/area_sum, Time)

end function send_global_land_diag


end module land_tile_diag_mod
