module land_fire_emis_mod

#include "../shared/debug.inc"

use fms_mod, only : input_nml_file, stdout, stdlog, check_nml_error, &
      error_mesg, stdlog, stdout, lowercase, uppercase, WARNING, FATAL, NOTE
use constants_mod, only: PI, AVOGNO
use land_constants_mod, only : days_per_year
use time_manager_mod, only : time_type, time_type_to_real
use diag_manager_mod, only : register_diag_field, send_data
use field_manager_mod , only : MODEL_ATMOS, MODEL_LAND, parse
use land_constants_mod, only : seconds_per_year
use land_debug_mod, only : check_var_range
use land_tile_mod, only : land_tile_type
use vegn_tile_mod, only : vegn_tile_type
use table_printer_mod
use vegn_data_mod, only: nspecies, spdata
use land_data_mod, only : lnd, log_version
use land_tile_diag_mod, only : set_default_diag_filter, &
        register_tiled_diag_field, send_tile_data
use land_fire_emis_data_mod, only : fire_emis_type, frdata, n_fire_tr

implicit none
private

! ==== public interfaces =====================================================
public :: land_fire_emis_init, land_fire_emis_end
! public :: land_fire_emis
! public :: fire_emis_type
public :: update_fire_emissions
public :: diag_fire_emissions

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

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('land')

  ! currently the only action here is the registration of the diagnostic fields
  allocate(id_fire_emis(n_fire_tr))
  do i = 1,n_fire_tr
     id_fire_emis(i) = register_tiled_diag_field( module_name, &
          trim(frdata(i)%name)//'_fire_emis', (/id_ug/), lnd%time, &
          'fire emission of '//trim(frdata(i)%name), 'molecules/cm2/s', missing_value=-1.0)
  enddo
  module_is_initialized = .TRUE.
end subroutine land_fire_emis_init

! Finish using the model: deallocate memory, etc.
subroutine land_fire_emis_end()
   if (allocated(id_fire_emis)) deallocate(id_fire_emis)

   module_is_initialized = .FALSE.
end subroutine land_fire_emis_end

! given amount of burned carbon for each vegetation species and rate of carbon emission,
! calculate fire emissions for each of the fire tracers and ave it in the array in the
! vegetation tile
subroutine update_fire_emissions(vegn, burned_SP)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: burned_SP(0:nspecies-1) ! amount of burned biomass per species, kgC/(m2 s)

  integer :: tr ! fire tracer index
  integer :: sp ! species index

  real :: efact ! average emission factor for given tracer
  real :: c_ave_ef ! average C per dry mass (g C)/(kg DM)
  real :: w     ! averaging weight
  real :: s     ! sum of averaging weights
  real :: csmoke_rate_daily ! fire carbon emission, kg C/(m2 day)
!   integer :: i

  csmoke_rate_daily = max(0.0, vegn%csmoke_rate / days_per_year)
  call check_var_range(burned_SP,          0.0, HUGE(1.0), 'update_fire_emissions', 'burned_SP',         WARNING)

  do tr = 1,n_fire_tr
     ! calculate average [i.e. effective] emission factors based on the burned biomass
     efact = 0.0; c_ave_ef = 0.0; s = 0.0;
     do sp = 0,nspecies-1
        w     = max(0.0,burned_SP(sp)) ! averaging weight
        efact = efact + frdata(tr)%efactors(sp) * w
        c_ave_ef = c_ave_ef + spdata(sp)%c_per_dry_matter * w

        s     = s + w
     enddo

! old averaging, like Arman did
!      do i = 1,vegn%n_cohorts
!         w     = 1.0
!         sp    = vegn%cohorts(i)%species
!         efact = efact + frdata(tr)%efactors(sp) * w
!         c_ave_ef = c_ave_ef + spdata(sp)%c_per_dry_matter * w
!         s     = s + w
!      enddo

     if (s>0) then
        efact = efact/s
        c_ave_ef = c_ave_ef/s
     else
        ! This could be happening when we do not burn, or burn litter only: then
        ! we take the emission factors from the species of the first cohort (which
        ! always exists and is also the tallest).
        sp    = vegn%cohorts(1)%species
        efact = frdata(tr)%efactors(sp)
        c_ave_ef = spdata(sp)%c_per_dry_matter
     endif
     vegn%fire_emis_land(tr) = efact * &
        csmoke_rate_daily * &
        (1./(c_ave_ef * 1.E-3)) * & !! convert C to DM in grams
        1.E-4 * &                   !! m2_to_cm2
        (1./(24.*60.*60.)) * &        !! per_second
        (1./frdata(tr)%fire_mw) * &
        AVOGNO
    enddo
    call check_var_range(vegn%fire_emis_land(:),  0.0, HUGE(1.0), 'update_fire_emissions', 'fire_emis_land', WARNING)
end subroutine update_fire_emissions

! ============================================================================
! send fire emissions to diagnostics
subroutine diag_fire_emissions(tile)
   type(land_tile_type), intent(inout) :: tile

   integer :: tr
   if (associated(tile%vegn)) then
      do tr = 1,n_fire_tr
         call send_tile_data(id_fire_emis(tr), tile%vegn%fire_emis_land(tr),  tile%diag)
      enddo
   else
      do tr = 1,n_fire_tr
         call send_tile_data(id_fire_emis(tr), 0.0,  tile%diag)
      enddo
   endif
end subroutine

! ============================================================================
! Added this subroutine to compute fire emissions using factors read in from namelist
! subroutine land_fire_emis(tile)
!    type(land_tile_type), intent(inout) :: tile
!
!    real  :: sp_ave_ef ! emission factors averaged over cohorts
!    real  :: c_ave_ef
!
!    real  :: csmoke_rate_daily = 0.0
!
!    integer :: sp ! shorthand for cohort species
!    integer :: i,j
!
!
!    do j = 1,n_fire_tr
!       c_ave_ef=0.
!       sp_ave_ef=0.
!
!       if (associated(tile%vegn)) then
!          associate(cc=>tile%vegn%cohorts)
!             do i = 1, tile%vegn%n_cohorts
!                sp = cc(i)%species
!
!                c_ave_ef  = c_ave_ef  + spdata(sp)%c_per_dry_matter
!                sp_ave_ef = sp_ave_ef + frdata(j)%efactors(sp)
!             end do
!          end associate ! cc
!
!          c_ave_ef  = c_ave_ef/tile%vegn%n_cohorts
!          sp_ave_ef = sp_ave_ef/tile%vegn%n_cohorts
!
!          csmoke_rate_daily = tile%vegn%csmoke_rate * (1./(365.))
!          if (csmoke_rate_daily < 0.0) csmoke_rate_daily = 0.0
!          tile%vegn%fire_emis_land(j) = sp_ave_ef * &
!                                     csmoke_rate_daily * &
!                                     (1./(c_ave_ef * 1.E-3)) * &!! convert C to DM in grams
!                                     1.E-4 * &                  !! m2_to_cm2
!                                     (1./(24.*60.*60.)) * &     !! per_second
!                                     (1./frdata(j)%fire_mw) * &
!                                     AVOGNO
!       endif
!    end do
!
!    do i = 1,n_fire_tr
!       call send_tile_data(id_fire_emis(i), tile%vegn%fire_emis_land(i),  tile%diag)
!    enddo
! end subroutine land_fire_emis

end module land_fire_emis_mod
