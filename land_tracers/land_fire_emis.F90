module land_fire_emis_mod

#include "../shared/debug.inc"

use fms_mod, only : input_nml_file, stdout, stdlog, check_nml_error, &
      error_mesg, stdlog, stdout, lowercase, uppercase, WARNING, FATAL, NOTE
use constants_mod, only: PI, AVOGNO
use time_manager_mod, only : time_type, time_type_to_real
use diag_manager_mod, only : register_diag_field, send_data
use field_manager_mod , only : MODEL_ATMOS, MODEL_LAND, parse
use land_constants_mod, only : seconds_per_year
use land_tile_mod, only : land_tile_type
use table_printer_mod
use vegn_data_mod, only: nspecies, spdata
use land_data_mod, only : lnd, log_version
use land_tile_diag_mod, only : set_default_diag_filter, &
        register_tiled_diag_field, send_tile_data
use land_fire_emis_data_mod, only : fire_emis_type, frdata, n_fire_tr

implicit none
private

! ==== public interfaces =====================================================
public  ::  land_fire_emis_init, land_fire_emis_end
public :: land_fire_emis
public :: fire_emis_type
!!! dsward_cpl end

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_fire_emis'
#include "../shared/version_variable.inc"

logical         :: module_is_initialized =.FALSE.
integer, allocatable :: id_fire_emis(:)


contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine land_fire_emis_init(id_ug)
  integer,intent(in) :: id_ug !<Unstructured axis id.
  integer :: i


  if (module_is_initialized) return

  ! currently the only action here is the registration of the diagnostic fields
  allocate(id_fire_emis(n_fire_tr))
  do i = 1,n_fire_tr
     id_fire_emis(i) = register_tiled_diag_field( module_name, &
          trim(frdata(i)%name)//'_fire_emis', (/id_ug/), lnd%time, &
          'fire emissions', 'molecules/cm2/s', missing_value=-1.0)
  enddo
  module_is_initialized = .TRUE.
end subroutine land_fire_emis_init

! Finish using the model: deallocate memory, etc.
subroutine land_fire_emis_end()
   if (allocated(id_fire_emis)) deallocate(id_fire_emis)

   module_is_initialized = .FALSE.
end subroutine land_fire_emis_end


! ============================================================================
! Added this subroutine to compute fire emissions using factors read in from namelist
subroutine land_fire_emis(tile)
   type(land_tile_type), intent(inout) :: tile

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
