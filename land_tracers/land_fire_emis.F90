module land_fire_emis_mod

#include "../shared/debug.inc"

! This stuff is boilerplate for making/reading namelists.
use mpp_mod, only: mpp_pe, mpp_root_pe
use mpp_mod, only: input_nml_file

use constants_mod,   only: PI
use time_manager_mod, only : time_type, get_date, days_in_month, operator(-), &
                             time_type_to_real
use fms_mod, only : check_nml_error, error_mesg, stdlog, stdout, &
      lowercase, WARNING, FATAL, NOTE
use fms2_io_mod, only: close_file, FmsNetcdfFile_t, open_file
use sphum_mod, only : qscomp
use diag_manager_mod, only : register_diag_field, send_data
use field_manager_mod , only : MODEL_ATMOS, MODEL_LAND, parse
use land_constants_mod, only : seconds_per_year
use land_tile_mod, only : land_tile_type
!!! dsward_cpl added several use statements, many to accompany vegn_tracer table creation
use field_manager_mod, only: fm_field_name_len, fm_string_len, &
     fm_type_name_len, fm_path_name_len, fm_dump_list, fm_get_length, &
     fm_get_current_list, fm_loop_over_list, fm_change_list
use fm_util_mod, only : fm_util_get_real, fm_util_get_logical, fm_util_get_string, fm_util_get_real_array
use tracer_manager_mod, only : NO_TRACER, get_number_tracers, get_tracer_names, get_tracer_index, &
                               query_method
use table_printer_mod
use mpp_mod, only: stdout, stdlog, mpp_error
!use land_data_mod, only: MAX_FR_TR
use vegn_data_mod, only: nspecies, spdata
use constants_mod, only: AVOGNO
use fms_mod, only : uppercase
use land_data_mod, only : lnd, log_version
use land_tile_diag_mod, only : set_default_diag_filter, &
        register_tiled_diag_field, send_tile_data


implicit none
private

! ==== public interfaces =====================================================
public  ::  land_fire_emis_init, land_fire_emis_end
public :: land_fire_emis
public :: fire_emis_type
!!! dsward_cpl end

integer, protected, public    :: n_fire_tr  ! number of fire emission tracers

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_fire_emis'

!namelist /land_fire_emis_nml/ &
!    fire_emission_factor_file

logical         :: module_is_initialized =.FALSE.
integer, parameter :: MAX_FR_TR = 99
real            :: delta_time
real            :: dt_fast_yr      ! fast time step in years
integer :: id_fire_emis(MAX_FR_TR)

type(table_printer_type) :: table


!--- Fire emissions type
type fire_emis_type
  character(fm_field_name_len) :: name        = ''  ! name of the tracer
  integer :: tr_atm  = NO_TRACER ! index of this tracer in atmos tracer array
  real  :: fire_mw = 1.0                            ! molecular weights of fire tracers
  real  :: efactors(6) = (/ 93., 127., 127., 88., 63., 63. /)  ! emission factors for fire emissions of the tracer species
end type
type(fire_emis_type), allocatable :: frdata(:) ! fire emissions data

!!! dsward_cpl end



contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#####################################################################
! initialize vegetation tracers
subroutine land_fire_emis_init(id_ug,frdata)
!####
  integer,intent(in) :: id_ug !<Unstructured axis id.
  type(fire_emis_type), allocatable, intent(inout) :: frdata(:)
 integer :: i, j, k, m, n, o, p, nsp, tr
 integer :: trind
 character(fm_field_name_len) :: name ! name of the vegn tracer
 character(fm_type_name_len)  :: typ  ! type of the vegn tracer
 integer :: nt_atmos
 character(len=32) :: ef_sp
 character(len=500)  :: method
 character(len=500) :: parameters
 real    :: value ! temporary storage for parsing input

 delta_time  = time_type_to_real(lnd%dt_fast)
 dt_fast_yr = delta_time/seconds_per_year

! get number of atmos_tracers
 call get_number_tracers (MODEL_ATMOS, num_tracers=nt_atmos)
 
  n_fire_tr = 0
! see if any of the atmos_tracers have bb_emis is lm4
  do tr = 1, nt_atmos
     call get_tracer_names (MODEL_ATMOS, tr, name = name)
     trind = get_tracer_index(MODEL_ATMOS,name)
     if(query_method('emissions2dbb', MODEL_ATMOS, trind, method, parameters)) then
        if (trim(method)=='lm4') then
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
        if (trim(method)=='lm4') then
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

 if (n_fire_tr .gt. MAX_FR_TR) call mpp_error(FATAL, 'Number of fire emission tracers defined exceeds the maximum of MAX_FR_TR -please increase MAX_FR_TR in vegn_data.F90')

  ! log tracer information
  call init_with_headers(table, frdata(:)%name)
  call add_row(table, 'atm.tr.number',   frdata(:)%tr_atm)
  call add_row(table, 'fire_mw',  frdata(:)%fire_mw)
  do nsp = 0, nspecies-1
     call add_row(table, 'ef_'//trim(spdata(nsp)%name),frdata(:)%efactors(nsp+1))
  enddo
  call print(table,stdlog())
  call print(table,stdout())

 
  do i = 1,n_fire_tr
     id_fire_emis(i) = register_tiled_diag_field( module_name, &
          trim(frdata(i)%name)//'_fire_emis', (/id_ug/), lnd%time, &
          'fire emissions', 'molecules/cm2/s', missing_value=-1.0)
  enddo
  module_is_initialized = .TRUE.
end subroutine land_fire_emis_init

subroutine land_fire_emis_end
   if (allocated(frdata)) then
      deallocate(frdata)
   else
      if (mpp_pe() == mpp_root_pe()) &
      call mpp_error(WARNING, 'frdata is not allocated; cannot deallocate')
   endif
end subroutine land_fire_emis_end


! ============================================================================
! Added this subroutine to compute fire emissions using factors read in from namelist
subroutine land_fire_emis(tile,frdata)
     
   type(land_tile_type), pointer :: tile
   type(fire_emis_type), allocatable, intent(inout) :: frdata(:)
     
   real  :: c_efactors(nspecies)  ! emission factors for fire emissions of C
   real  :: c_efact_def(nspecies)  ! emission factors for fire emissions of C
   real  :: sp_ave_ef ! emission factors averaged over cohorts
   real  :: c_ave_ef
     
!!! temporary:
   real  :: temp_csmoke_rate = 1.0
   real  :: csmoke_rate_daily = 0.0
     
   integer :: sp ! shorthand for cohort species 
   integer :: i,j,ln,k,l
  
     
   !!! Use C emission factors for these tracers
   !!! to compute dry matter burned from C lost.

   c_efactors(1:nspecies)=0.0
   c_efact_def=(/491.751, 464.989, 464.989, 489.416, 488.273, 488.273/)
   do i=1,n_fire_tr
      if (uppercase(frdata(i)%name(1:2))=='C') c_efactors(1:nspecies)=frdata(i)%efactors(1:nspecies)
   enddo

   do j = 1,n_fire_tr

      c_ave_ef=0.
      sp_ave_ef=0.

      if (associated(tile%vegn)) then
         associate(cc=>tile%vegn%cohorts)
             do i = 1, tile%vegn%n_cohorts

             sp = cc(i)%species+1   !!! added the plus one since efactors are indexed to one

             if (c_efactors(sp)==0.0) c_efactors(sp)=c_efact_def(sp)

             c_ave_ef=c_ave_ef+c_efactors(sp)
             sp_ave_ef=sp_ave_ef+frdata(j)%efactors(sp)

            end do
         end associate ! cc

      c_ave_ef=c_ave_ef/tile%vegn%n_cohorts
      sp_ave_ef=sp_ave_ef/tile%vegn%n_cohorts
      
      csmoke_rate_daily = tile%vegn%csmoke_rate * (1./(365.)) 
      if (csmoke_rate_daily < 0.0) csmoke_rate_daily = 0.0
      tile%vegn%fire_emis_land(j) = sp_ave_ef * &
                                    csmoke_rate_daily * &
                                    (1./(c_ave_ef * 1.E-3)) * &!! convert C to DM in grams
                                    1.E-4 * &                  !! m2_to_cm2
                                    (1./(24.*60.*60.)) * &     !! per_second
                                    (1./frdata(j)%fire_mw) * &
                                    AVOGNO

      endif      
   end do
  

   do i = 1,n_fire_tr
      call send_tile_data(id_fire_emis(i), tile%vegn%fire_emis_land(i),  tile%diag)
   enddo

end subroutine land_fire_emis 

end module land_fire_emis_mod 
