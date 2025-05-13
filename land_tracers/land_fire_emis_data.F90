module land_fire_emis_data_mod

use fms_mod, only: stdout, stdlog, error_mesg, string, NOTE, WARNING, FATAL
use field_manager_mod, only: fm_field_name_len, &
     fm_type_name_len, MODEL_ATMOS, MODEL_LAND, parse
use tracer_manager_mod, only: NO_TRACER, get_number_tracers, get_tracer_names, &
     get_tracer_index, query_method
use gex_mod, only : gex_get_index

use vegn_data_mod, only: nspecies, spdata
use land_data_mod, only: log_version
use table_printer_mod

implicit none
private

! ---- public interfaces
public :: init_fire_emis_data
public :: fire_emis_type
public :: n_fire_tr  ! number of fire tracers
public :: tr_gex_frp ! index of fire radiative power in GEX array of fields
public :: frdata


! ---- module constants
character(*), parameter :: module_name = 'land_fire_emis_data_mod'
character(*), parameter :: data_error_header = 'FIRE TRACER DATA FATAL ERROR'
#include "../shared/version_variable.inc"

! ---- structure holding fire emission data
type fire_emis_type
  character(fm_field_name_len) :: name = ''  ! name of the tracer
  integer :: tr_atm  = NO_TRACER ! index of this tracer in atmos tracer array
  integer :: tr_gex  = NO_TRACER ! index of this tracer in GEX tracer array
  real    :: fire_mw = 1.0       ! molecular weight of this fire tracers
  real, allocatable  :: efactors(:) ! emission factors for fire emissions of the
                                 ! tracer species
end type

! module variables
logical :: module_is_initialized =.FALSE.

integer, protected :: n_fire_tr  =  0 ! number of fire tracers
integer, protected :: tr_gex_frp = -1 ! index of fire radiative power in GEX array of fields
type(fire_emis_type), allocatable, protected :: frdata(:) ! fire emissions data
contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Read fire emission data from field tables and initialize data structures.
subroutine init_fire_emis_data()

  integer :: i, sp, nsp, tr
  integer :: total_errors
  character(fm_field_name_len) :: name ! name of the vegn tracer
  character(fm_type_name_len)  :: typ  ! type of the vegn tracer
  integer :: nt_atmos
  character(len=500) :: method
  character(len=500) :: parameters
  real    :: value ! temporary storage for parsing input
  real, allocatable :: value1(:) ! temporary storage solely for tracer data table output
  type(table_printer_type) :: table

  if (module_is_initialized) return ! do nothing further

  call log_version(version, module_name, __FILE__)

  ! get number of atmos_tracers
  call get_number_tracers (MODEL_ATMOS, num_tracers=nt_atmos)

  n_fire_tr = 0
  ! see if any of the atmos_tracers have bb_emis is land:lm4
  do tr = 1, nt_atmos
     call get_tracer_names (MODEL_ATMOS, tr, name = name)
     if(query_method('emissions2dbb', MODEL_ATMOS, tr, method, parameters)) then
        if (trim(method)=='land:lm4') then
           n_fire_tr=n_fire_tr+1
        endif
     endif
  enddo

  if (n_fire_tr .eq. 0) then
     call error_mesg(module_name, 'No interactive fire emission tracers found in ATMOS tracer table', NOTE)
     module_is_initialized = .TRUE.
     return
  endif

  allocate(frdata(1:n_fire_tr))

  i = 0; total_errors = 0
  ! register the frdata info
  do tr = 1, nt_atmos
     call get_tracer_names (MODEL_ATMOS, tr, name = name)
     method = ''; parameters = ''
     if(query_method('emissions2dbb', MODEL_ATMOS, tr, method, parameters)) then
        if (trim(method)=='land:lm4') then
           i = i + 1
           frdata(i)%name = trim(name)
           frdata(i)%tr_atm = get_tracer_index(MODEL_ATMOS,name)
           frdata(i)%tr_gex = gex_get_index( MODEL_LAND,MODEL_ATMOS, 'fire_emis_'//frdata(i)%name)
           if (frdata(i)%tr_gex.le.0) then
              total_errors = total_errors + 1
              call error_mesg(data_error_header, &
                    'Tracer "fire_emis_'//trim(frdata(i)%name)//'" not found in gex_lnd2atm', WARNING)
           endif
           if ( parse(parameters, 'mw', value) > 0 ) then
              frdata(i)%fire_mw = value
           endif
!            allocate(frdata(i)%efactors(0:nspecies-1))
           allocate(frdata(i)%efactors(nspecies))
           frdata(i)%efactors(:) = -1.0
           do sp = 0, nspecies-1
              if (trim(spdata(sp)%name)=='default') cycle ! skip emission for fake "default" species that should never appear in model's vegetation
              nsp = sp+1
              if ( parse(parameters, 'ef_'//trim(spdata(sp)%name), value) > 0 ) then
                 frdata(i)%efactors(nsp) = value
              else
                 total_errors = total_errors + 1
                 call error_mesg(data_error_header, &
                       trim(frdata(i)%name)//' fire emission factor for species "'//trim(spdata(sp)%name)//'" not found in the field table', WARNING)
              endif
           enddo
        endif
     endif
  enddo

  tr_gex_frp = gex_get_index( MODEL_LAND,MODEL_ATMOS, 'frp')
  if (tr_gex_frp.le.0) then
     total_errors = total_errors + 1
     call error_mesg(data_error_header, &
           'Tracer "frp" not found in gex_lnd2atm', NOTE)
  endif

  ! log tracer information
  call init_with_headers(table, frdata(:)%name)
  call add_row(table, 'atm.tr.number',   frdata(:)%tr_atm)
  call add_row(table, 'GEX.tr.number',   frdata(:)%tr_gex)
  call add_row(table, 'fire_mw',         frdata(:)%fire_mw)

  allocate(value1(n_fire_tr))
  do sp = 0, nspecies-1
     nsp = sp+1
     do tr = 1, n_fire_tr
        value1(tr) = frdata(tr)%efactors(nsp)
     enddo
     call add_row(table, 'ef_'//trim(spdata(sp)%name),value1(:))
  enddo
  call print(table,stdlog(),transposed=.TRUE.)
  call print(table,stdout(),transposed=.TRUE.)
  deallocate(value1)

  if (total_errors > 0) then
     call error_mesg(module_name, trim(string(total_errors))//' errors found in species parameters tables, look for "'//&
                                  data_error_header//'" in this output', FATAL)
  endif

  module_is_initialized = .TRUE.
end subroutine init_fire_emis_data

end module land_fire_emis_data_mod