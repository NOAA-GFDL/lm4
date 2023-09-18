module soilc_util_mod

use time_manager_mod, only : time_type
use land_constants_mod, only : N_C_TYPES, N_LITTER_POOLS, &
        c_longname, c_diagname,  l_longname, l_diagname
use tile_diag_base_mod, only : register_tiled_diag_field

implicit none; private

public :: register_soilc_diag_fields
public :: register_litter_diag_fields
public :: register_litter_soilc_diag_fields

contains

! ============================================================================
function replace_text (s,text,rep)  result(outs)
character(*), intent(in) :: s,text,rep
character(len(s)+100) :: outs     ! provide outs with extra 100 char len

integer      :: i, nt, nr

outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
do
   i = index(outs,text(:nt)) ; if (i == 0) exit
   outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
end do
end function replace_text

! ============================================================================
function register_soilc_diag_fields(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name) result (id)

  integer :: id(N_C_TYPES)

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(len=*), intent(in), optional :: op ! aggregation operation
  character(len=*), intent(in), optional :: standard_name

  integer :: i

  do i = 1, N_C_TYPES
     id(i) = register_tiled_diag_field(module_name, &
             trim(replace_text(field_name,'<ctype>',trim(c_diagname(i)))), &
             axes, init_time, &
             trim(replace_text(long_name,'<ctype>',trim(c_longname(i)))), &
             units, missing_value, range, op, standard_name)
  enddo
end function register_soilc_diag_fields

! ============================================================================
! registered an array of diag fields, one per litter pool
function register_litter_diag_fields(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name) result (id)

  integer :: id(N_LITTER_POOLS)

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(len=*), intent(in), optional :: op ! aggregation operation
  character(len=*), intent(in), optional :: standard_name

  integer :: i

  do i = 1, N_LITTER_POOLS
     id(i) = register_tiled_diag_field(module_name, &
             trim(replace_text(field_name,'<ltype>',trim(l_diagname(i)))), &
             axes, init_time, &
             trim(replace_text(long_name,'<ltype>',trim(l_longname(i)))), &
             units, missing_value, range, op, standard_name)
  enddo
end function register_litter_diag_fields

! ============================================================================
! registered a 2D array of diag fields, one per litter pool per carbon type
function register_litter_soilc_diag_fields(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name) result (id)

  integer :: id(N_LITTER_POOLS, N_C_TYPES)

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  character(len=*), intent(in), optional :: op ! aggregation operation
  character(len=*), intent(in), optional :: standard_name

  integer :: i, k
  character(128) :: name
  character(512) :: lname

  do i = 1, N_C_TYPES
     do k = 1, N_LITTER_POOLS
        name = replace_text(field_name,'<ctype>',trim(c_diagname(i)))
        name = replace_text(name,      '<ltype>',trim(l_diagname(k)))
        lname = replace_text(long_name,'<ctype>',trim(c_longname(i)))
        lname = replace_text(lname,    '<ltype>',trim(l_longname(k)))
        id(k,i) = register_tiled_diag_field(module_name, trim(name), axes, init_time, trim(lname), &
             units, missing_value, range, op, standard_name)
     enddo
  enddo
end function register_litter_soilc_diag_fields


end module