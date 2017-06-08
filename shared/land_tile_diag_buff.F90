module tile_diag_buff_mod

implicit none
private

! ==== public interfaces =====================================================
public :: diag_buff_type
public :: init_diag_buff, realloc_diag_buff
! ==== end of public interfaces ==============================================

! storage for tile diagnostic data
type :: diag_buff_type
   real   , allocatable :: data(:)
   logical, allocatable :: mask(:)
end type diag_buff_type

! ==== module constants =====================================================
integer, parameter :: MIN_DIAG_BUFF_SIZE = 1

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine init_diag_buff(buffer)
  type(diag_buff_type), intent(inout) :: buffer

  allocate(buffer%data(MIN_DIAG_BUFF_SIZE), buffer%mask(MIN_DIAG_BUFF_SIZE))
  ! initialize buffer content
  buffer%mask(:) = .FALSE.
  buffer%data(:) = 0.0
end subroutine init_diag_buff

! ============================================================================
! reallocates buffer to have at least m elements
subroutine realloc_diag_buff(buffer, m)
  type(diag_buff_type), intent(inout) :: buffer
  integer             , intent(in)    :: m

  real    , allocatable :: new_data(:)
  logical , allocatable :: new_mask(:)
  integer               :: n

  ! n is size of the original buffer; m is the current size of the buffer
  ! for all diagnostic fields
  n = size(buffer%data(:))
  ! do nothing if buffer is big enough
  if(n >= m) return

  allocate(new_data(m), new_mask(m))
  new_data(1:n) = buffer%data(1:n) ; new_data(n+1:m) = 0.0
  new_mask(1:n) = buffer%mask(1:n) ; new_mask(n+1:m) = .FALSE.
  call move_alloc(new_data,buffer%data)
  call move_alloc(new_mask,buffer%mask)
end subroutine realloc_diag_buff

end module tile_diag_buff_mod
