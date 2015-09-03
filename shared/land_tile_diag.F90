module land_tile_diag_mod

use mpp_mod,            only : mpp_sum
use time_manager_mod,   only : time_type
use diag_axis_mod,      only : get_axis_length
use diag_manager_mod,   only : register_diag_field, register_static_field, &
     send_data
use diag_util_mod,      only : log_diag_field_info
use fms_mod,            only : write_version_number, error_mesg, string, FATAL

use land_tile_selectors_mod, only : tile_selectors_init, tile_selectors_end, &
     tile_selector_type, register_tile_selector, selector_suffix, &
     get_n_selectors, get_selector
use land_tile_mod,      only : land_tile_type, diag_buff_type, &
     land_tile_list_type, first_elmt, tail_elmt, next_elmt, get_elmt_indices, &
     land_tile_enum_type, operator(/=), current_tile, &
     tile_is_selected, fptr_i0, fptr_r0, fptr_r0i
use vegn_cohort_mod,    only : vegn_cohort_type
use land_data_mod,      only : lnd
use tile_diag_buff_mod, only : diag_buff_type, realloc_diag_buff

implicit none
private


! ==== public interface ======================================================
public :: tile_diag_init
public :: tile_diag_end

public :: diag_buff_type

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

interface send_tile_data
   module procedure send_tile_data_0d
   module procedure send_tile_data_1d
end interface
interface send_cohort_data
  module procedure send_cohort_data_with_weight
  module procedure send_cohort_data_without_weight
end interface
! ==== end of public interface ===============================================


! ==== module constants ======================================================
character(len=*), parameter :: &
     mod_name = 'land_tile_diag_mod', &
     version  = '$Id$', &
     tagname  = '$Name$'

integer, parameter :: INIT_FIELDS_SIZE      = 1       ! initial size of the fields array
integer, parameter :: BASE_TILED_FIELD_ID   = 65536   ! base value for tiled field 
integer, parameter :: BASE_COHORT_FIELD_ID  = 65536*2 ! base value for cohort field ids
! ids, to distinguish them from regular diagnostic fields. All IDs of tiled
! (that is, registered by register_*tiled_field functions) are larger than 
! BASE_TILED_FIELD_ID. All cohort field IDs are larger than BASE_COHORT_FIELD_ID

integer, parameter :: MIN_DIAG_BUFFER_SIZE = 1     ! min size of the per-tile diagnostic buffer
! operations used for tile data aggregation
integer, parameter, public :: &
    OP_MEAN = 1, & ! weighted average of tile values
    OP_SUM  = 2, & ! sum of all tile values
    OP_MAX  = 3, & ! maximum of all  values
    OP_MIN  = 4, & ! minimum of all  values
    OP_VAR  = 5, & ! variance of tile values
    OP_STD  = 6, & ! standard deviation of tile values
    OP_DOMINANT = 7 ! dominant value (only for cohorts now)
character(32), parameter :: opstrings(6) = (/ & ! symbolica names of the aggregation operations
   'mean                            ' , &
   'sum                             ' , &
   'maximum                         ' , &
   'minimum                         ' , &
   'variance                        ' , &
   'stdev                           '   /)

! TODO: generalize treatment of cohort filters. Possible filters include: selected 
! species, selected species in the canopy, trees above certain age, etc...
integer, parameter :: N_CHRT_FILTERS = 3 ! number of pssible distinct cohort filters,
! currently : all vegetation, upper canopy layer, and understory
character(2),  parameter :: chrt_filter_suffix(N_CHRT_FILTERS) = (/'  ','_1','_U'/)
character(32), parameter :: chrt_filter_name(N_CHRT_FILTERS)   = (/&
    '                                ', &
    ' in top canopy layer            ', &
    ' in understory                  '  /)


! ==== derived types =========================================================
type :: tiled_diag_field_type
   integer, pointer :: ids(:) => NULL()
   integer :: offset ! offset of the field data in the buffer
   integer :: size   ! size of the field data in the per-tile buffers
   integer :: opcode ! aggregation operation
   logical :: static ! if true, the diag field is static
   integer :: n_sends! number of data points sent to the field since last dump
   integer :: alias = 0 ! ID of the first alias in the chain
   character(32) :: module,name ! for debugging purposes only
end type tiled_diag_field_type

type :: cohort_diag_field_type
   integer :: ids(N_CHRT_FILTERS) = -1 ! IDs of tiled diag fields for each of cohort filters
end type cohort_diag_field_type


! ==== module data ===========================================================
logical :: module_is_initialized = .false.

! list of registered fields
type(tiled_diag_field_type), pointer :: fields(:) => NULL()
integer :: n_fields       = 0 ! current number of diag fields
integer :: current_offset = 1 ! current total size of the diag fields per tile

! list of registered cohort fields
type(cohort_diag_field_type), pointer :: cfields(:) => NULL()
integer :: n_cfields     = 0 ! current number of diag fields


contains



! ============================================================================
subroutine tile_diag_init()

  if (module_is_initialized) return

  module_is_initialized = .true.
  call write_version_number(version, tagname)

  ! initialize diag selectors
  call tile_selectors_init()
  call register_tile_selector('')

  ! initialize global data
  allocate(fields(INIT_FIELDS_SIZE))
  n_fields       = 0
  current_offset = 1

  ! initialize array of cohort fields
  allocate(cfields(INIT_FIELDS_SIZE))
  n_cfields       = 0

end subroutine tile_diag_init



! ============================================================================
subroutine tile_diag_end()

  integer :: i

  ! deallocate global data
  do i = 1, n_fields
     deallocate(fields(i)%ids)
  end do
  deallocate(fields)
  deallocate(cfields)

  ! destroy selectors
  call tile_selectors_end()

  module_is_initialized = .false.

end subroutine tile_diag_end


! ============================================================================
function string2opcode(op) result(code)
  integer :: code ! return value
  character(*), intent(in) :: op ! aggregation operation
  
  integer :: i
  
  code =-1
  do i = 1,size(opstrings)
     if (trim(op)==trim(opstrings(i))) code = i
  enddo
end function string2opcode


! ============================================================================
function register_tiled_diag_field(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op) result (id)

  integer :: id

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(*),     intent(in), optional :: op ! aggregation operation
  
  id = reg_field(.false., module_name, field_name, init_time, axes, long_name, &
         units, missing_value, range, op=op)

end function register_tiled_diag_field

! ============================================================================
function register_tiled_static_field(module_name, field_name, axes, &
     long_name, units, missing_value, range, require, op) result (id)

  integer :: id

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  logical,          intent(in), optional :: require
  character(*),     intent(in), optional :: op ! aggregation operation
  
  ! --- local vars
  type(time_type) :: init_time

  id = reg_field(.true., module_name, field_name, init_time, axes, long_name, &
         units, missing_value, range, require, op)

end function register_tiled_static_field


! ============================================================================
subroutine add_tiled_static_field_alias(id0, module_name, field_name, axes, &
     long_name, units, missing_value, range, op)
  integer,          intent(inout) :: id0 ! id of the original diag field on input;
   ! if negative then it may be replaced with the alias id on output
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(*),     intent(in), optional :: op ! aggregation operation

  ! --- local vars
  type(time_type) :: init_time

  call reg_field_alias(id0, .true., module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op)
end subroutine add_tiled_static_field_alias


! ============================================================================
subroutine add_tiled_diag_field_alias(id0, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op)
  integer,          intent(inout) :: id0 ! id of the original diag field on input;
   ! if negative then it may be replaced with the alias id on output
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(*),     intent(in), optional :: op ! aggregation operation

  call reg_field_alias(id0, .false., module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op)
end subroutine add_tiled_diag_field_alias

! ============================================================================
subroutine reg_field_alias(id0, static, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op)


  integer,          intent(inout) :: id0 ! id of the original diag field on input;
  logical,          intent(in) :: static
   ! if negative then it may be replaced with the alias id on output
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(*),     intent(in), optional :: op ! aggregation operation
  
  ! local vars
  integer :: id1
  integer :: ifld0, ifld1
  
  if (id0>0) then
    ifld0 = id0-BASE_TILED_FIELD_ID
    if (ifld0<1.or.ifld0>n_fields) &
       call error_mesg(mod_name, 'incorrect index ifld0 '//string(ifld0)//  &
                    ' in definition of tiled diag field alias "'//          &
                    trim(module_name)//'/'//trim(field_name)//'"', FATAL)
    id1 = reg_field(static, module_name, field_name, init_time, axes, long_name, &
          units, missing_value, range, op=op, offset=fields(ifld0)%offset)
    if (id1>0) then
      ifld1 = id1-BASE_TILED_FIELD_ID
      ! check that sizes of the fields are identical
      if (fields(ifld0)%size/=fields(ifld1)%size) &
         call error_mesg(mod_name, 'sizes of diag field "'//              &
           trim(fields(ifld0)%module)//'/'//trim(fields(ifld0)%name)//    &
           '" and its alias "'//trim(module_name)//'/'//trim(field_name)//&
           '" are not the same', FATAL)
      ! check that "static" status of the fields is the same
      if(fields(ifld0)%static.and..not.fields(ifld1)%static) &
         call error_mesg(mod_name,                                        &
           'attempt to register non-static alias"'//                      &
           trim(module_name)//'/'//trim(field_name)//                     &
           '" of static field "'//                                        &
           trim(fields(ifld0)%module)//'/'//trim(fields(ifld0)%name)//'"',&
           FATAL)
      if(.not.fields(ifld0)%static.and.fields(ifld1)%static) &
         call error_mesg(mod_name,                                        &
           'attempt to register static alias"'//                          &
           trim(module_name)//'/'//trim(field_name)//                     &
           '" of non-static field "'//                                    &
           trim(fields(ifld0)%module)//'/'//trim(fields(ifld0)%name)//'"',&
           FATAL)

      ! copy alias field from the original into the alias, to preserve the chain
      fields(ifld1)%alias = fields(ifld0)%alias
      ! update alias field in the head of alias chain
      fields(ifld0)%alias = ifld1
    endif
  else
    ! the "main" field has not been registered, so simply register the alias
    ! as a diag field
    id0 = reg_field(static, module_name, field_name, init_time, axes, long_name, &
          units, missing_value, range, op=op)
  endif
end subroutine reg_field_alias

! ============================================================================
! provides unified interface for registering a diagnostic field with full set
! of selectors
function reg_field(static, module_name, field_name, init_time, axes, &
     long_name, units, missing_value, range, require, op, offset) result(id)
 
  integer :: id

  logical,          intent(in) :: static
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  logical,          intent(in), optional :: require
  character(*),     intent(in), optional :: op ! aggregation operation
  integer,          intent(in), optional :: offset

  ! ---- local vars
  integer, pointer :: diag_ids(:) ! ids returned by FMS diag manager for each selector
  integer :: i
  integer :: isel    ! selector iterator
  integer :: opcode  ! code of the tile aggregation operation
  integer :: n_selectors ! number of registered diagnostic tile selectors
  type(tiled_diag_field_type), pointer :: new_fields(:)
  type(tile_selector_type) :: sel
  ! ---- global vars: n_fields, fields, current_offset -- all used and updated

  ! log diagnostic field information
  call log_diag_field_info ( module_name, trim(field_name), axes, long_name, units,&
                             missing_value, range, dynamic=.not.static )
  ! go through all possible selectors and try to register a diagnostic field 
  ! with the name derived from field name and selector; if any of the 
  ! registrations succeeds, return a tiled field id, otherwise return 0.
  ! Note that by design one of the selectors have empty name and selects all
  ! the tiles.
  id = 0
  n_selectors = get_n_selectors()
  allocate(diag_ids(n_selectors))
  
  do isel = 1, n_selectors
     ! try to register field+selector pair with FMS diagnostics manager
     sel = get_selector(isel)
     diag_ids(isel) = reg_field_set(static, sel, module_name, field_name, axes, &
          init_time, long_name, units, missing_value, range, require)

  enddo
  
  if(any(diag_ids>0)) then
     ! if any of the field+selector pairs was found for this field, an entry
     ! must be added to the table of tile diagnostic fields

     ! if there is not enough slots in the field table to add another one,
     ! allocate more space
     if(n_fields>=size(fields)) then
        allocate(new_fields(max(2*n_fields,1)))
        new_fields(1:n_fields) = fields(1:n_fields)
        deallocate(fields)
        fields => new_fields
     endif
     ! add the current field to the field table
     n_fields = n_fields+1
     id       = n_fields
     ! set the array of FMS diagnostic field IDs for each selector
     fields(id)%ids => diag_ids
     ! set the field offset in the diagnostic buffers
     if (present(offset)) then
        fields(id)%offset = offset
     else
        fields(id)%offset = current_offset
     endif  
     ! calculate field size per tile and increment current offset to
     ! reserve space in per-tile buffers. We assume that the first two axes 
     ! are horizontal coordinates, so their size is not taken into account
     fields(id)%size = 1
     do i = 3, size(axes(:))
        fields(id)%size = fields(id)%size * get_axis_length(axes(i))
     enddo
     ! if offset is present in the list of the arguments, it means that we don't
     ! want to increase the current_offset -- this is an alias field
     if (.not.present(offset)) &
        current_offset = current_offset + fields(id)%size
     ! store the code of the requested tile aggregation operation
     if(present(op)) then
        fields(id)%opcode = string2opcode(op)
        if (.not.(fields(id)%opcode > 0)) &
           call error_mesg(mod_name,&
              'tile aggregarion operation "'//trim(op)//'" for field "'&
              //trim(module_name)//'/'//trim(field_name)//'"is incorrect, use "mean", "sum", "variance", or "stdev"',&
              FATAL)
     else
        fields(id)%opcode = OP_MEAN
     endif
     ! store the static field flag
     fields(id)%static = static
     ! zero out the number of data points ent to the field
     fields(id)%n_sends = 0
     ! store the name of the field -- for now, only to be able to see what it is 
     ! in the debugger
     fields(id)%module=module_name 
     fields(id)%name=field_name
     ! increment the field id by some (large) number to distinguish it from the 
     ! IDs of regular FMS diagnostic fields
     id = id + BASE_TILED_FIELD_ID
  else
     deallocate(diag_ids)
  endif

end function reg_field


! ============================================================================
! provides unified interface for registering a diagnostic field with a given
! selector, whether static or time-dependent
function reg_field_set(static, sel, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, require) result (id)

  integer :: id 

  logical,          intent(in) :: static
  type(tile_selector_type), intent(in) :: sel
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  logical,          intent(in), optional :: require

  character(len=128) :: fname
  character(len=128) :: lname

  ! form field name as concatenation of name of the field and selector suffix
  fname = trim(field_name)//trim(selector_suffix(sel))
  ! form long name as concatenation of specified long name (if present) and
  ! selector long name
  lname = ''
  if(present(long_name)) lname=long_name
  if(trim(sel%long_name)/='') &
     lname = trim(lname)//' ('//trim(sel%long_name)//')'

  ! try registering diagnostic field with FMS diagnostic manager.
  if (static) then
     id = register_static_field ( module_name, fname,   &
          axes, lname, units, missing_value, range, require, do_not_log=.TRUE. )
  else
     id = register_diag_field ( module_name,  fname,   &
          axes, init_time, lname, units, missing_value, range, &
          mask_variant=.true., do_not_log=.TRUE. )
  endif

end function reg_field_set


! ============================================================================
subroutine send_tile_data_0d(id, x, buffer)
  integer, intent(in) :: id
  real   , intent(in) :: x
  type(diag_buff_type), intent(inout) :: buffer
  
  integer :: idx, i

  if (id <= 0) return
  if (id < BASE_TILED_FIELD_ID.or. id >= BASE_COHORT_FIELD_ID ) call error_mesg (mod_name, &
         'tile diag field ID is out of range. Perhaps the field was not registred with some other call then register_tile_diag_field?', &
         FATAL)

  ! reallocate diagnostic buffer according to the current number and size of 
  ! tiled diag fields
  call realloc_diag_buff(buffer,current_offset)

  ! calculate offset for the current diagnostic field in the buffer
  i = id - BASE_TILED_FIELD_ID 
  idx = fields(i)%offset
  
  ! store the diagnostic data
  buffer%data(idx) = x
  buffer%mask(idx) = .TRUE.

  ! increment sent data counter
  fields(i)%n_sends = fields(i)%n_sends + 1
  ! increment sent data counter in all aliases
  do while(fields(i)%alias>0)
    i=fields(i)%alias
    fields(i)%n_sends = fields(i)%n_sends + 1
  enddo
end subroutine send_tile_data_0d

! ============================================================================
subroutine send_tile_data_1d(id, x, buffer)
  integer, intent(in) :: id
  real   , intent(in) :: x(:)
  type(diag_buff_type), intent(inout) :: buffer

  integer :: is, ie, i
  if (id <= 0) return
  if (id < BASE_TILED_FIELD_ID.or. id >= BASE_COHORT_FIELD_ID ) call error_mesg (mod_name, &
         'tile diag field ID is out of range. Perhaps the field was not registred with some other call then register_tile_diag_field?', &
         FATAL)

  ! reallocate diagnostic buffer according to the current number and size of 
  ! tiled diag fields
  call realloc_diag_buff(buffer, current_offset)
  
  ! calculate offset for the current diagnostic field in the buffer
  i = id - BASE_TILED_FIELD_ID ! index in the array of fields
  is = fields(i)%offset ; ie = is+fields(i)%size-1

  ! store the data
  buffer%data(is:ie) = x(:)
  buffer%mask(is:ie) = .TRUE.

  ! increment sent data counter
  fields(i)%n_sends = fields(i)%n_sends + 1
  ! increment sent data counter in all aliases
  do while(fields(i)%alias>0)
    i=fields(i)%alias
    fields(i)%n_sends = fields(i)%n_sends + 1
  enddo
end subroutine send_tile_data_1d

! NOTE: 2-d fields can be handled similarly to 1-d with reshape

! ============================================================================
subroutine send_tile_data_r0d_fptr(id, tile_map, fptr)
  integer, intent(in) :: id
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  procedure(fptr_r0)  :: fptr

  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  real                , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,ptr,tileptr%diag)
     ce=next_elmt(ce)
  enddo
end subroutine send_tile_data_r0d_fptr


! ============================================================================
subroutine send_tile_data_r1d_fptr(id, tile_map, fptr)
  integer, intent(in) :: id
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  procedure(fptr_r0i) :: fptr

  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  real                , pointer :: ptr     ! pointer to the data element within a tile
  real,             allocatable :: buffer(:) ! buffer for accumulating data
  integer :: i, k

  if(id <= 0) return
  if (id < BASE_TILED_FIELD_ID.or. id >= BASE_COHORT_FIELD_ID ) call error_mesg (mod_name, &
         'tile diag field ID is out of range. Perhaps the field was not registred with some other call then register_tile_diag_field?', &
         FATAL)
  i = id - BASE_TILED_FIELD_ID ! index in the array of fields

  allocate(buffer(fields(i)%size))
  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     tileptr => current_tile(ce)
     do k = 1,fields(i)%size
        call fptr(tileptr,k,ptr)
        if(associated(ptr)) buffer(k) = ptr
     enddo
     call send_tile_data(id,buffer,tileptr%diag)
     ce=next_elmt(ce)
  enddo
  deallocate(buffer)
end subroutine send_tile_data_r1d_fptr


! ============================================================================
subroutine send_tile_data_i0d_fptr(id, tile_map, fptr)
  integer, intent(in) :: id
  type(land_tile_list_type), intent(inout) :: tile_map(:,:)
  procedure(fptr_i0)  :: fptr

  type(land_tile_enum_type)     :: te,ce   ! tail and current tile list elements
  type(land_tile_type), pointer :: tileptr ! pointer to tile   
  integer             , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( tile_map )
  te = tail_elmt ( tile_map )
  do while(ce /= te)
     tileptr => current_tile(ce)
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,real(ptr),tileptr%diag)
     ce=next_elmt(ce)
  enddo
end subroutine send_tile_data_i0d_fptr


! ============================================================================
subroutine dump_tile_diag_fields(tiles, time)
  type(land_tile_list_type), intent(in) :: tiles(:,:) ! 
  type(time_type)          , intent(in) :: time       ! current time

  ! ---- local vars
  integer :: ifld ! field number
  integer :: isel ! selector number
  type(land_tile_enum_type)     :: ce, te
  type(land_tile_type), pointer :: tile
  integer :: total_n_sends(n_fields)
  ! ---- local static variables -- saved between calls
  logical :: first_dump = .TRUE.

  total_n_sends(:) = fields(1:n_fields)%n_sends
  call mpp_sum(total_n_sends, n_fields, pelist=lnd%pelist)

!$OMP parallel do schedule(dynamic) default(shared) private(ifld,isel)
  do ifld = 1, n_fields
     if (total_n_sends(ifld) == 0) cycle ! no data to send 
     do isel = 1, get_n_selectors()
        if (fields(ifld)%ids(isel) <= 0) cycle
        call dump_diag_field_with_sel ( fields(ifld)%ids(isel), tiles, &
             fields(ifld), get_selector(isel), time )
     enddo
  enddo
  ! zero out the number of data points sent to the field 
  fields(1:n_fields)%n_sends=0

  ! all the data are sent to the output, so set the data presence tag to FALSE 
  ! in all diag buffers in preparation for the next time step
  ce = first_elmt(tiles)
  te = tail_elmt (tiles)
  do while(ce /= te)
    tile => current_tile(ce)       ! get the pointer to the current tile
    tile%diag%mask(:) = .FALSE.
    ce = next_elmt(ce)            ! move to the next position
  enddo
  ! reset the first_dump flag
  first_dump = .FALSE.

end subroutine dump_tile_diag_fields

! ============================================================================
subroutine dump_diag_field_with_sel(id, tiles, field, sel, time)
  integer :: id
  type(land_tile_list_type),   intent(in) :: tiles(:,:)
  type(tiled_diag_field_type), intent(in) :: field
  type(tile_selector_type)   , intent(in) :: sel
  type(time_type)            , intent(in) :: time ! current time
   
  ! ---- local vars
  integer :: i,j ! iterators
  integer :: is,ie,js,je,ks,ke ! array boundaries
  logical :: used ! value returned from send_data (ignored)
  real, allocatable :: buffer(:,:,:), weight(:,:,:), var(:,:,:)
  type(land_tile_enum_type)     :: ce, te
  type(land_tile_type), pointer :: tile
  
  ! calculate array boundaries
  is = lbound(tiles,1); ie = ubound(tiles,1)
  js = lbound(tiles,2); je = ubound(tiles,2)
  ks = field%offset   ; ke = field%offset + field%size - 1
  
  ! allocate and initialize temporary buffers
  allocate(buffer(is:ie,js:je,ks:ke), weight(is:ie,js:je,ks:ke))
  buffer(:,:,:) = 0.0
  weight(:,:,:) = 0.0
  
  ! accumulate data
  ce = first_elmt(tiles, is=is, js=js)
  te = tail_elmt (tiles)
  do while(ce /= te)
    tile => current_tile(ce)      ! get the pointer to current tile
    call get_elmt_indices(ce,i,j) ! get the indices of current tile
    ce = next_elmt(ce)           ! move to the next position
    
    if ( size(tile%diag%data) < ke )       cycle ! do nothing if there is no data in the buffer
    if ( .not.tile_is_selected(tile,sel) ) cycle ! do nothing if tile is not selected
    select case (field%opcode)
    case (OP_MEAN,OP_VAR,OP_STD)
       where(tile%diag%mask(ks:ke)) 
          buffer(i,j,:) = buffer(i,j,:) + tile%frac*tile%diag%data(ks:ke)
          weight(i,j,:) = weight(i,j,:) + tile%frac
       end where
    case (OP_SUM)
       where(tile%diag%mask(ks:ke)) 
          buffer(i,j,:) = buffer(i,j,:) + tile%diag%data(ks:ke)
          weight(i,j,:) = 1
       end where
    end select
  enddo

  ! normalize accumulated data
  where (weight>0) buffer=buffer/weight
  
  if (field%opcode == OP_VAR.or.field%opcode == OP_STD) then
     ! second loop to process the variance and standard deviation diagnostics.
     ! it may be possible to calc. var and std in one pass with weighted incremental 
     ! algorithm from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
     ! code her is more straightforward. buffer(:,:,:) already contains the mean,
     ! and weight(:,:,:) -- sum of  tile fractions
     allocate(var(is:ie,js:je,ks:ke))
     var(:,:,:) = 0.0
     ! the loop is somewhat different from the first, for no particular reason: 
     ! perhaps this way is better for performance?
     do j = js,je
     do i = is,ie
        ce = first_elmt(tiles(i,j))
        te = tail_elmt (tiles(i,j))
        do while(ce /= te)
           tile => current_tile(ce) ! get the pointer to current tile
           ce = next_elmt(ce)       ! move to the next position
    
           if ( size(tile%diag%data) < ke )       cycle ! do nothing if there is no data in the buffer
           if ( .not.tile_is_selected(tile,sel) ) cycle ! do nothing if tile is not selected
           where(tile%diag%mask(ks:ke))
              var(i,j,:) = var(i,j,:) + tile%frac*(tile%diag%data(ks:ke)-buffer(i,j,:))**2
           end where 
        enddo
     enddo
     enddo
     ! renormalize the variance or standard deviation. note that weight is 
     ! calculated in the first loop
     select case (field%opcode)
     case (OP_VAR)
         where (weight>0) buffer = var/weight
     case (OP_STD)
         where (weight>0) buffer = sqrt(var/weight)
     end select
     deallocate(var)
  endif
  
  ! send diag field
  used = send_data ( id, buffer, time, mask=weight>0 )   

  ! clean up temporary data
  deallocate(buffer,weight)

end subroutine dump_diag_field_with_sel

! ============================================================================
function register_cohort_diag_field(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, opt) result (id)

  integer :: id

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(*),     intent(in), optional :: opt ! tile aggregation operation

  integer :: i
  integer :: diag_ids(N_CHRT_FILTERS)
  type(cohort_diag_field_type), pointer :: new_fields(:)
  character(256) :: long_name_

  ! register tiled diag field for each of the samplings
  do i = 1,N_CHRT_FILTERS
     if(present(long_name)) then
         long_name_ = trim(long_name)//trim(chrt_filter_name(i))
         diag_ids(i) = register_tiled_diag_field( &
            module_name, trim(field_name)//trim(chrt_filter_suffix(i)), axes, init_time, &
            long_name_, units, missing_value, range, opt)
     else
         diag_ids(i) = register_tiled_diag_field( &
            module_name, trim(field_name)//trim(chrt_filter_suffix(i)), axes, init_time, &
            long_name, units, missing_value, range, opt)
     endif
  enddo
  id = 0
  if (any(diag_ids(:)>0)) then
     ! if there is not enough slots in the field table to add another one,
     ! allocate more space
     if(n_cfields>=size(cfields)) then
        allocate(new_fields(max(2*n_cfields,1)))
        new_fields(1:n_cfields) = cfields(1:n_cfields)
        deallocate(cfields)
        cfields => new_fields
     endif
     ! add the current field to the field table
     n_cfields = n_cfields+1
     id       = n_cfields     
     cfields(id)%ids = diag_ids
     id       = id+BASE_COHORT_FIELD_ID
  endif
end function register_cohort_diag_field


! ============================================================================
subroutine send_cohort_data_without_weight (id, buffer, cc, data, op)
  integer,                intent(in)    :: id        ! diag field ID
  type(diag_buff_type),   intent(inout) :: buffer    ! data buffer
  type(vegn_cohort_type), intent(in)    :: cc(:)     ! cohort array: used for filtering
  real,                   intent(in)    :: data(:)   ! data
  integer,                intent(in)    :: op        ! cohort aggregation operation

  real :: weight(size(data))
  weight(:) = 1.0
  call send_cohort_data_with_weight (id, buffer, cc, data, weight, op)
end subroutine send_cohort_data_without_weight

! ============================================================================
subroutine send_cohort_data_with_weight (id, buffer, cc, data, weight, op)
  integer,                intent(in)    :: id        ! diag field ID
  type(diag_buff_type),   intent(inout) :: buffer    ! data buffer
  type(vegn_cohort_type), intent(in)    :: cc(:)     ! cohort array: used for filtering
  real,                   intent(in)    :: data(:)      ! data
  real,                   intent(in)    :: weight(:) ! averaging weight
  integer,                intent(in)    :: op        ! cohort aggregation operation
  
  integer :: i
  real    :: value
  
  if (id<=0) return

  ! check data size
  if (size(cc)/=size(data)) call error_mesg(mod_name, &
     'size of cohort array is not equal to the size of data', FATAL)
  if (size(weight)/=size(data)) call error_mesg(mod_name,&
     'size of data is not equal to size of weight',FATAL)

  ! check that ID is in correct range
  if (id < BASE_COHORT_FIELD_ID.or.id>BASE_COHORT_FIELD_ID+n_cfields) &
     call error_mesg(mod_name,&
         'cohort diag field ID is out of range. Perhaps the field was not registred with register_cohort_diag_field?', &
         FATAL)

  i = id - BASE_COHORT_FIELD_ID
  ! TODO: generalize cohort subsampling
  if (cfields(i)%ids(1) > 0) then
     ! all cohorts
     value = aggregate(data,weight,op,mask=cc(:)%layer>0)
     call send_tile_data(cfields(i)%ids(1), value, buffer)
  endif
  if (cfields(i)%ids(2) > 0) then
     ! canopy cohorts
     value = aggregate(data,weight,op,mask=cc(:)%layer==1)
     call send_tile_data(cfields(i)%ids(2), value, buffer)
  endif
  if (cfields(i)%ids(3) > 0) then
     ! understory cohorts
     value = aggregate(data,weight,op,mask=cc(:)%layer>1)
     call send_tile_data(cfields(i)%ids(3), value, buffer)
  endif
end subroutine send_cohort_data_with_weight


! ============================================================================
function aggregate(data,weight,opcode,mask) result(ret)
  real :: ret
  real, intent(in) :: data(:),weight(:)
  integer, intent(in) :: opcode
  logical, intent(in) :: mask(:)

  real :: w
  integer :: i
  
  select case(opcode)
  case(OP_SUM)
     ret = sum(data*weight,mask=mask)
  case(OP_MEAN)
     w = sum(weight,mask=mask)
     ret = sum(data*weight,mask=mask)
     if (w/=0) ret = ret/w
  case(OP_MAX)
     ret = 0.0
     if (any(mask)) ret = maxval(data,mask)
  case(OP_MIN)
     ret = 0.0
     if (any(mask)) ret = minval(data,mask)
  case(OP_DOMINANT)
     ret = 0.0
     w=-HUGE(1.0)
     do i = 1,size(weight)
        if (mask(i).and.weight(i)>w) then 
           ret = data(i); w = weight(i)
        endif
     enddo
  case default
     call error_mesg(mod_name, 'unrecognized cohort data aggregation opcode', FATAL)
  end select
end function aggregate


end module land_tile_diag_mod
