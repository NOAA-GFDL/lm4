module tile_diag_buff_mod

implicit none
private

! ==== public interfaces =====================================================
public :: diag_buff_type
public :: new_diag_buff, delete_diag_buff, realloc_diag_buff

! storage for tile diagnostic data
type :: diag_buff_type
   real   , pointer :: data(:) => NULL()
   logical, pointer :: mask(:) => NULL()
end type diag_buff_type

! ==== module constants =====================================================
integer, parameter :: MIN_DIAG_BUFF_SIZE = 1

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
function new_diag_buff() result(buffer)
  type(diag_buff_type), pointer :: buffer
  
  integer :: m ! initial size of the buffer
  
  allocate(buffer)
  m = MIN_DIAG_BUFF_SIZE
  allocate(buffer%mask(m),buffer%data(m))
  ! initialize buffer content
  buffer%mask(:) = .FALSE.
  buffer%data(:) = 0.0

end function


! ============================================================================
subroutine delete_diag_buff(buffer)
  type(diag_buff_type), pointer :: buffer
  
  if(.not.associated(buffer)) return
  deallocate(buffer%mask,buffer%data)
  deallocate(buffer)
  
end subroutine


! ============================================================================
! reallocates buffer to have at least m elements
subroutine realloc_diag_buff(buffer, m)
  type(diag_buff_type), intent(inout) :: buffer
  integer             , intent(in)    :: m 

  real    , pointer :: new_data(:)
  logical , pointer :: new_mask(:)
  integer           :: n

  ! n is size of the original buffer; m is the current size of the buffer
  ! for all diagnostic fields
  n = size(buffer%data(:))
  ! do nothing if buffer is big enough
  if(n >= m) return

  allocate(new_data(m), new_mask(m))
  new_data(1:n) = buffer%data(1:n) ; new_data(n+1:m) = 0.0
  new_mask(1:n) = buffer%mask(1:n) ; new_mask(n+1:m) = .FALSE.
  deallocate(buffer%data, buffer%mask)

  buffer%data=>new_data
  buffer%mask=>new_mask

end subroutine


end module tile_diag_buff_mod
