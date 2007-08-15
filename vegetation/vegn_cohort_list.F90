module cohort_list_mod

use vegn_cohort_mod, only : vegn_cohort_type

implicit none
private
! ==== public interfaces =====================================================
public :: vegn_cohort_list_type
public :: vegn_cohort_enum_type

public :: cohort_list_init, cohort_list_end
public :: first_cohort, tail_cohort
public :: operator(==), operator(/=) ! comparison of enumerators
public :: next_cohort, prev_cohort   ! enumerator advance operations
public :: current_cohort ! returns pointer to cohort at a position
public :: insert_cohort  ! inserts a cohort at a given position, or appends it to a list
public :: remove_cohort  ! given a position, removes the cohort from the list
public :: erase_cohort   ! given a position, removes the cohort from the list and deletes it

! ==== end of public interfaces ==============================================
interface operator(==)
   module procedure enums_are_equal
end interface
interface operator(/=)
   module procedure enums_are_not_equal
end interface
interface insert_cohort
   module procedure insert_at_position, insert_in_list
end interface
interface remove_cohort
   module procedure remove_at_position
end interface
interface erase_cohort
   module procedure erase_at_position
end interface

! ==== types =================================================================
! vegn_cohort_list_type provides a container for the cohorts
type :: vegn_cohort_list_type
   private
   type(list_node_type), pointer :: head => NULL()
end type vegn_cohort_list_type

! vegn_cohort_enum_type provides a enumerator of cohorts -- a data structure
! that allows to walk through all cohorts in a container without bothering 
! with details of container implementation
type :: vegn_cohort_enum_type
   private
   type(list_node_type), pointer :: &
        node => NULL()                ! pointer to the current container node
end type vegn_cohort_enum_type

! private type -- used internally to implement lists
type :: list_node_type
   type(list_node_type), pointer :: prev => NULL()
   type(list_node_type), pointer :: next => NULL()
   type(vegn_cohort_type)   , pointer :: data => NULL()
end type list_node_type

contains

! ============================================================================
! cohort list constructor: initializes essential innards of cohort collection
! for future use. In current implementation, it is safe to call this function
! on a cohort list more then once.
subroutine cohort_list_init(list)
  type(vegn_cohort_list_type), intent(inout) :: list

  if (.not.associated(list%head)) then
     allocate(list%head)
     list%head%prev=>list%head
     list%head%next=>list%head
  endif
end subroutine cohort_list_init

! ============================================================================
! cohort list destructor: destroys the list of cohorts. NOTE that it also 
! deallocates all the cohorts that are still in the list.
subroutine cohort_list_end(list)
  type(vegn_cohort_list_type), intent(inout) :: list

  type(vegn_cohort_enum_type) :: ce

  if(associated(list%head)) then
     ce=first_cohort(list)
     do while(ce/=tail_cohort(list))
        call erase_at_position(ce)
     enddo
     deallocate(list%head)
  endif
end subroutine cohort_list_end

! ============================================================================
! returns the number of items currently stored in the list
function n_items_in_list(list) result (n)
  type(vegn_cohort_list_type), intent(in) :: list
  integer :: n

  type(list_node_type), pointer :: node

  n=0; 
  if(.not.associated(list%head)) return

  node => list%head%next
  do while ( .not.(associated(node,list%head)) )
     n = n+1
     node => node%next
  enddo
end function n_items_in_list

! ============================================================================
! inserts cohort in the list (more precisely, appends it to the end)
subroutine insert_in_list(cohort,list)
  type(vegn_cohort_type),            pointer :: cohort
  type(vegn_cohort_list_type), intent(inout) :: list

  call insert_at_position(cohort,tail_cohort(list))
end subroutine insert_in_list


! #### cohort container enumerator ###########################################

! ============================================================================
! return enumerator pointing to the first element of the list
function first_cohort(list) result(ce)
  type(vegn_cohort_enum_type) :: ce ! return value
  type(vegn_cohort_list_type), intent(in) :: list

  ce%node=>list%head%next 
end function first_cohort

! ============================================================================
! returns enumerator pointing to the tail of the list (the element just behind 
! end of the list)
function tail_cohort(list) result(ce)
  type(vegn_cohort_enum_type) :: ce ! return value
  type(vegn_cohort_list_type), intent(in) :: list

  ce%node=>list%head 
end function tail_cohort

! ============================================================================
! returns enumerator pointing to the next element of the container.
function next_cohort(pos) result(ce)
  type(vegn_cohort_enum_type) :: ce ! return value
  type(vegn_cohort_enum_type), intent(in) :: pos

  ce = pos
  ce%node => ce%node%next
end function next_cohort

! ============================================================================
! returns enumerator pointing to the previous element of the container.
function prev_cohort(pos) result(ce)
  type(vegn_cohort_enum_type) :: ce ! return value
  type(vegn_cohort_enum_type), intent(in) :: pos

  ce = pos
  ce%node => ce%node%prev
end function prev_cohort

! ============================================================================
! returns pointer to the element currently addressed by the enumerator
function current_cohort(ce) result(ptr)
  type(vegn_cohort_type), pointer :: ptr ! return value
  type(vegn_cohort_enum_type), intent(in) :: ce 

  ptr => ce%node%data
end function 

! ============================================================================
! returns TRUE if both enums refer to the same list node (and, therefore, cohort)
! or if both do not refer to anything. 
function enums_are_equal(pos1,pos2) result(ret)
  logical :: ret ! return value
  type(vegn_cohort_enum_type), intent(in) :: pos1,pos2
  
  if(associated(pos1%node)) then
     ret = associated(pos1%node,pos2%node)
  else
     ret = .not.associated(pos2%node)
  endif
end function enums_are_equal

! ============================================================================
! returns TRUE if two enumerators are not equal
function enums_are_not_equal(pos1,pos2) result(ret)
  logical :: ret ! return value
  type(vegn_cohort_enum_type), intent(in) :: pos1,pos2

  ret=.not.enums_are_equal(pos1,pos2)
end function enums_are_not_equal

! ============================================================================
! inserts cohort at the position indicated by enumerator
subroutine insert_at_position(cohort,ce)
  type(vegn_cohort_type),         pointer :: cohort
  type(vegn_cohort_enum_type), intent(in) :: ce

  ! local vars
  type(list_node_type), pointer :: node,n,p

  allocate(node)
  node%data=>cohort

  n=>ce%node  ; p=>n%prev

  node%next=>n ; node%prev=>p
  n%prev=>node ; p%next=>node
end subroutine insert_at_position

! ============================================================================
subroutine remove_at_position(ce)
  type(vegn_cohort_enum_type), intent(inout) :: ce

  type(list_node_type),pointer :: n,p
  n => ce%node%next
  p => ce%node%prev

  n%prev=>p ; p%next=>n
  deallocate(ce%node)
  ce%node=>n

end subroutine remove_at_position

! ============================================================================
subroutine erase_at_position(ce)
  type(vegn_cohort_enum_type), intent(inout) :: ce

  type(vegn_cohort_type), pointer :: cohort

  cohort=>current_cohort(ce)
  call remove_at_position(ce)
  deallocate(cohort)
end subroutine erase_at_position


end module cohort_list_mod
