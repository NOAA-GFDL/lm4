module land_fire_emis_data_mod

use mpp_mod, only: stdout, stdlog
use field_manager_mod, only: fm_field_name_len, &
     fm_type_name_len, MODEL_ATMOS, MODEL_LAND, parse
use tracer_manager_mod, only: NO_TRACER, get_number_tracers, get_tracer_names, &
     get_tracer_index, query_method

use vegn_data_mod, only: nspecies, spdata
use land_data_mod, only: log_version
use table_printer_mod

implicit none
private

! ---- public interfaces
public :: init_fire_emis_data
public :: fire_emis_type
public :: n_fire_tr ! number of fire tracers
public :: frdata


! ---- module constants
character(len=*), parameter :: module_name = 'glac_tile_mod'
#include "../shared/version_variable.inc"

! ---- structure holding fire emission data
type fire_emis_type
  character(fm_field_name_len) :: name = ''  ! name of the tracer
  integer :: tr_atm  = NO_TRACER ! index of this tracer in atmos tracer array
  real    :: fire_mw = 1.0                            ! molecular weights of fire tracers
  real    :: efactors(6) = (/ 93., 127., 127., 88., 63., 63. /)  ! emission factors for fire emissions of the tracer species
end type

! module variables
logical :: module_is_initialized =.FALSE.

integer, protected :: n_fire_tr ! number of fire tracers
type(fire_emis_type), allocatable, protected :: frdata(:) ! fire emissions data

contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Read fire emission data from field tables and initialize data structures.
subroutine init_fire_emis_data()

  integer :: i, nsp, tr
  integer :: trind
  character(fm_field_name_len) :: name ! name of the vegn tracer
  character(fm_type_name_len)  :: typ  ! type of the vegn tracer
  integer :: nt_atmos
  character(len=500) :: method
  character(len=500) :: parameters
  real    :: value ! temporary storage for parsing input
  type(table_printer_type) :: table

  call log_version(version, module_name, __FILE__)

  ! get number of atmos_tracers
  call get_number_tracers (MODEL_ATMOS, num_tracers=nt_atmos)

  n_fire_tr = 0
  ! see if any of the atmos_tracers have bb_emis is land:lm4
  do tr = 1, nt_atmos
     call get_tracer_names (MODEL_ATMOS, tr, name = name)
     trind = get_tracer_index(MODEL_ATMOS,name)
     if(query_method('emissions2dbb', MODEL_ATMOS, trind, method, parameters)) then
        if (trim(method)=='land:lm4') then
           n_fire_tr=n_fire_tr+1
        endif
     endif
  enddo

  allocate(frdata(1:n_fire_tr))

  i = 0
  ! register the frdata info
  do tr = 1, nt_atmos
     call get_tracer_names (MODEL_ATMOS, tr, name = name)
     trind = get_tracer_index(MODEL_ATMOS,name)
     method = ''; parameters = ''
     if(query_method('emissions2dbb', MODEL_ATMOS, trind, method, parameters)) then
        if (trim(method)=='land:lm4') then
           i = i + 1
           frdata(i)%name = trim(name)
           frdata(i)%tr_atm = get_tracer_index(MODEL_ATMOS,name)
           if ( parse(parameters, 'mw', value) > 0 ) then
              frdata(i)%fire_mw = value
           endif
           do nsp = 0, nspecies-1
              if ( parse(parameters, 'ef_'//trim(spdata(nsp)%name), value) > 0 ) then
                 frdata(i)%efactors(nsp+1) = value
              endif
           enddo
        endif
     endif
  enddo

  ! log tracer information
  call init_with_headers(table, frdata(:)%name)
  call add_row(table, 'atm.tr.number',   frdata(:)%tr_atm)
  call add_row(table, 'fire_mw',  frdata(:)%fire_mw)
  do nsp = 0, nspecies-1
     call add_row(table, 'ef_'//trim(spdata(nsp)%name),frdata(:)%efactors(nsp+1))
  enddo
!   call print(table,stdlog(),transposed=.TRUE.)
  call print(table,stdlog())
  call print(table,stdout())

  module_is_initialized = .TRUE.
end subroutine init_fire_emis_data

end module land_fire_emis_data_mod