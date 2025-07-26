module river_tracers_mod

use mpp_mod, only : mpp_error
use fms_mod, only : stdlog, stdout, string, FATAL, NOTE

use field_manager_mod, only: fm_field_name_len, fm_string_len, &
   fm_type_name_len, fm_path_name_len, fm_dump_list, fm_get_length, &
   fm_get_current_list, fm_loop_over_list, fm_change_list
use fm_util_mod, only : fm_util_get_real, fm_util_get_logical, fm_util_get_string
use tracer_manager_mod, only : NO_TRACER
use table_printer_mod


implicit none
private

!--- version information ---------------------------------------------
character(len=*), parameter :: module_name = 'river_tracers_mod'
#include "../shared/version_variable.inc"

!--- public interface ------------------------------------------------
public :: river_tracers_init
public :: num_river_tracers
public :: river_tracer_index
public :: river_tracer_names
!--- end of public interface -----------------------------------------

!--- tracer-related constants, types, and data
character(*), parameter :: trtable='/land_mod/river_tracer' ! name of the field manager tracer table
integer :: num_species  ! index of last tracer in zero-based table "trdata"
! In river modules, three "tracers" are always defined and hardcoded
! to occupy slots 0, 1, and 2 of the trdata table: h2o (total water) is 0,
! heat (called "het") is 1, and ice is 3.

integer, public, parameter :: num_phys = 2 ! number of "physical" tracers: currently they are ice and heat content

type tracer_data_type
  character(fm_field_name_len) :: &
      name        = '', & ! name of the tracer
      units       = '', & ! units of the tracer
      flux_units  = '', & ! units of associated flux
      store_units = ''    ! units of associated storage
  character(fm_string_len)     :: longname = '' ! longname of the species
  ! tracer removal parameters:
  logical :: do_removal = .true.
  real :: &
      t_ref  = 298.0, &
      vf_ref = 0.0,   &
      q10    = 1.0,   &
      kinv   = 1.0
end type

! In river code, three "tracers" are always defined and hardcoded to occupy
! slots 0, 1, and 2 of the trdata table, in this specific order: h2o (total
! water), ice, and heat (called "het"). While they cannot be turned off through
! omission in the field table, their parameters can be set there. The other river
! tracers are optional.
!
! The tracer 0 (total water) is treated separately in the river physics code,
! while other tracers are assigned an index in the river%*_c arrays. For example,
! river%storage_c(nx,ny,nt) represents storage of tracers in the rivers, "tr"
! being the tracer index (in the range [1:num_species]). Ice and heat are
! sometimes called "physical traces", in contrast to other substances that may be
! present in the rivers.
!
! Size of "trdata" array can be larger than the total number of river tracers, for
! the convenience of allocation code.
type(tracer_data_type), allocatable, public, protected :: trdata(:) ! common tracer data

contains

! ===-----------------------------------------------------------------------=== !
! initialize river tracers
subroutine river_tracers_init()

 integer :: i, m, n
 character(fm_field_name_len) :: name ! name of the river tracer
 character(fm_type_name_len)  :: typ  ! type of the river tracer

 ! number of river tracers in the field table (can be 0)
 m = fm_get_length(trtable)

 ! dump river tracer table
 if(.not.fm_dump_list(trtable, recursive=.TRUE.)) &
    call mpp_error(NOTE, 'river_tracers_mod: Cannot dump field list "'//trtable//'"')

 ! allocating more space than absolutely necessary, in case water, and "physical
 ! tracers" (ice and heat) are not present in the user-supplied tracer table
 allocate(trdata(0:m+num_phys))
 ! initialize some parameters of the pre-defined species (water and "physical" tracers)
 trdata(0)%name = 'h2o'; trdata(0)%longname = 'h2o mass'
 trdata(0)%units = 'kg'; trdata(0)%flux_units = 'kg/m2/s'; trdata(0)%store_units = 'kg/m2'

 trdata(1)%name = 'ice'; trdata(1)%longname = 'ice mass'
 trdata(1)%units = 'kg/kg'; trdata(1)%flux_units = 'kg/m2/s'; trdata(1)%store_units = 'kg/m2'

 trdata(2)%name = 'het'; trdata(2)%longname = 'sensible heat content'
 trdata(2)%units = 'K'; trdata(2)%flux_units = 'W/m2'; trdata(2)%store_units = 'J/m2'

 ! read generic parameters of the tracers
 do while (fm_loop_over_list(trtable, name, typ, n))
    ! look for the tracer already in the table
    do i = 0,ubound(trdata,1)
       if (trim(trdata(i)%name)==trim(name)) exit ! found existing slot for this tracer
    enddo
    ! if tracer not found, look for an empty slot in the table
    if (i>=ubound(trdata,1)) then
       do i = 0, ubound(trdata,1)
          if (trim(trdata(i)%name)=='') exit ! found an empty slot
       enddo
    endif
    call read_river_tracer_data(name,trdata(i))
 enddo
 ! finally, calculate the actual number of tracers
 do num_species = ubound(trdata,1),0,-1
    if (trdata(num_species)%name/='') exit ! from loop
 enddo

 ! TODO: read specific tracer parameters. Different tracers might have different parameter sets.

 call print_river_tracer_data(stdout())
 call print_river_tracer_data(stdlog())

end subroutine river_tracers_init

! ===-----------------------------------------------------------------------=== !
subroutine river_tracer_names(tr,name,long_name,units,flux_units,store_units)
  integer, intent(in) :: tr
  character(*), intent(out), optional :: name, long_name, units, flux_units, store_units

  if (tr<0.or.tr>num_species) call mpp_error( FATAL, &
     'river_tracers_mod: tracer index '//string(tr)//' is outside of range of river tracers')
  if(present(name))        name        = trdata(tr)%name
  if(present(long_name))   long_name   = trdata(tr)%longname
  if(present(units))       units       = trdata(tr)%units
  if(present(store_units)) store_units = trdata(tr)%store_units
  if(present(flux_units))  flux_units  = trdata(tr)%flux_units
end subroutine river_tracer_names

! ===-----------------------------------------------------------------------=== !
! reads the field_table entry for specific tracers and fills in
! generic tracer parameters
subroutine read_river_tracer_data(name,tr)
  character(*), intent(in) :: name
  type(tracer_data_type), intent(inout) :: tr

  ! ---- local vars
  character(fm_path_name_len)  :: listname
  character(fm_path_name_len)  :: current_list

  current_list = fm_get_current_list()
  if (current_list .eq. ' ') call mpp_error(FATAL, 'river_tracers_mod: Could not get the current list')
  listname = trtable//'/'//trim(name)
  if (.not.fm_change_list(listname)) call mpp_error(FATAL,'river_tracers_mod: Cannot change field manager list to "'//trim(listname)//'"')

  tr%name        = name
  tr%longname    = fm_util_get_string('long_name',   caller='river_tracers_mod', default_value=tr%longname,    scalar=.true.)
  tr%units       = fm_util_get_string('units',       caller='river_tracers_mod', default_value=tr%units,       scalar=.true.)
  tr%flux_units  = fm_util_get_string('flux_units',  caller='river_tracers_mod', default_value=tr%flux_units,  scalar=.true.)
  tr%store_units = fm_util_get_string('store_units', caller='river_tracers_mod', default_value=tr%store_units, scalar=.true.)
  ! tracer removal parameters:
  tr%do_removal  = fm_util_get_logical('do_removal', caller='river_tracers_mod', default_value=tr%do_removal,  scalar=.true.)
#define __PARSE__(v) tr%v = fm_util_get_real(#v, caller='river_tracers_mod', default_value=tr%v, scalar=.true.)
  __PARSE__(t_ref)
  __PARSE__(vf_ref)
  __PARSE__(q10)
  __PARSE__(kinv)
#undef __PARSE__

  if (.not.fm_change_list(current_list)) call mpp_error(FATAL,'river_tracers_mod: Cannot change field manager list to "'//trim(listname)//'"')
end subroutine read_river_tracer_data

! ===-----------------------------------------------------------------------=== !
! prints a table of tracer data to specified output unit
subroutine print_river_tracer_data(unit)
  integer, intent(in) :: unit

  type(table_printer_type) :: table

  call init_with_headers(table, trdata(:)%name)
  call add_row(table, 'longname', trdata(:)%longname)
  call add_row(table, 'units',    trdata(:)%units)
  call add_row(table, 'flux_units', trdata(:)%flux_units)
  call add_row(table, 'store_units', trdata(:)%store_units)
  call add_row(table, 'do_removal', trdata(:)%do_removal)
  call add_row(table, 't_ref', trdata(:)%t_ref)
  call add_row(table, 'vf_ref', trdata(:)%vf_ref)
  call add_row(table, 'q10', trdata(:)%q10)
  call add_row(table, 'kinv', trdata(:)%kinv)

  call print(table,unit)
end subroutine print_river_tracer_data

! ===-----------------------------------------------------------------------=== !
integer function num_river_tracers()
   num_river_tracers = num_species
end function num_river_tracers

! ===-----------------------------------------------------------------------=== !
function river_tracer_index(name) result(tr)
   integer :: tr
   character(*), intent(in) :: name

   integer :: i

   tr = NO_TRACER
   do i = 1, num_species
      if (name==trdata(i)%name) then
         tr = i;
         exit
      endif
   enddo
end function river_tracer_index

end module river_tracers_mod