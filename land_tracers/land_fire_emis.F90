module land_fire_emis_mod

#include "../shared/debug.inc"

use fms_mod, only : input_nml_file, stdout, stdlog, check_nml_error, &
      error_mesg, stdlog, stdout, lowercase, uppercase, WARNING, FATAL, NOTE
use field_manager_mod , only : MODEL_ATMOS, MODEL_LAND, parse
use constants_mod, only: PI, AVOGNO

use land_constants_mod, only : days_per_year
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, &
     get_tile_data, field_exists
use land_debug_mod, only : check_var_range
use land_tile_mod, only : land_tile_type
use vegn_tile_mod, only : vegn_tile_type
use vegn_accessors_mod, only : vegn_tile_exists, vegn_fire_emis_land_ptr
use vegn_data_mod, only: nspecies, spdata
use land_data_mod, only : lnd, log_version
use land_tile_diag_mod, only : set_default_diag_filter, &
        register_tiled_diag_field, send_tile_data
use land_fire_emis_data_mod, only : fire_emis_type, frdata, n_fire_tr

use table_printer_mod

implicit none
private

! ==== public interfaces =====================================================
public :: land_fire_emis_init, land_fire_emis_end
public :: save_fire_emis_restart
public :: update_fire_emissions
public :: diag_fire_emissions

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_fire_emis'
#include "../shared/version_variable.inc"

logical :: module_is_initialized =.FALSE.
integer, allocatable :: id_fire_emis(:)


contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ============================================================================
subroutine land_fire_emis_init(id_ug)
  integer,intent(in) :: id_ug !<Unstructured axis id.

  type(land_restart_type) :: restart
  logical :: restart_exists
  integer :: tr

  if (module_is_initialized) return

  call open_land_restart(restart,'INPUT/land_fire_emis.nc',restart_exists)
  if (restart_exists) then
     call error_mesg('land_fire_emis_init',&
          'reading NetCDF restarts "INPUT/land_fire_emis.res.nc"',&
          NOTE)
     do tr = 1, n_fire_tr
        if (field_exists(restart,trim(frdata(tr)%name)//'_fire_emis')) &
            call get_tile_data(restart,trim(frdata(tr)%name)//'_fire_emis',vegn_fire_emis_land_ptr,tr)
     enddo
     call free_land_restart(restart)
  else
     call error_mesg('land_fire_emis_init', 'cold-starting land fire emissions', NOTE)
  endif

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('land')

  ! register fire emission diagnostic fields
  allocate(id_fire_emis(n_fire_tr))
  do tr = 1,n_fire_tr
     id_fire_emis(tr) = register_tiled_diag_field( module_name, &
          trim(frdata(tr)%name)//'_fire_emis', (/id_ug/), lnd%time, &
          'fire emission of '//trim(frdata(tr)%name), 'molecules/cm2/s', missing_value=-1.0)
  enddo

  module_is_initialized = .TRUE.
end subroutine land_fire_emis_init

! ============================================================================
! Finish using the model: deallocate memory, etc.
subroutine land_fire_emis_end()
   if (allocated(id_fire_emis)) deallocate(id_fire_emis)

   module_is_initialized = .FALSE.
end subroutine land_fire_emis_end

! ============================================================================
subroutine save_fire_emis_restart(tile_dim_length,timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file (max number of tiles per grid cell)
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object
  integer :: tr

  if (n_fire_tr > 0) then
     call error_mesg('land_fire_emis_end','writing NetCDF restart',NOTE)

     ! create output file, including internal structure necessary for tile output
     filename = 'RESTART/'//trim(timestamp)//'land_fire_emis.nc'
     call init_land_restart(restart, filename, vegn_tile_exists, tile_dim_length)

     do tr = 1,n_fire_tr
        call add_tile_data(restart,trim(frdata(tr)%name)//'_fire_emis',vegn_fire_emis_land_ptr,tr,&
                           'fire emission of '//trim(frdata(tr)%name), 'molecules/cm2/s')
     enddo
     call save_land_restart(restart)
     call free_land_restart(restart)
  else
     call error_mesg('land_fire_emis_end','No fire tracers, NOT writing NetCDF restart',NOTE)
  endif
end subroutine save_fire_emis_restart

! ============================================================================
! given amount of burned carbon for each vegetation species and rate of carbon emission,
! calculate fire emissions for each of the fire tracers and ave it in the array in the
! vegetation tile
subroutine update_fire_emissions(vegn, burned_by_sp)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: burned_by_sp(0:nspecies-1) ! amount of burned biomass per species, kgC/(m2 s)

  integer :: tr ! fire tracer index
  integer :: sp ! species index

  real :: s     ! sum of averaging weights
  real :: csmoke_rate_daily ! fire carbon emission, kg C/(m2 day)
  real :: csmoke_by_sp(0:nspecies-1)

  call check_var_range(burned_by_sp, 0.0, HUGE(1.0), 'update_fire_emissions', 'burned_by_sp', WARNING)

  csmoke_rate_daily = max(0.0, vegn%csmoke_rate / days_per_year)
  ! Total burned biomass is typically lower than the total burned carbon because the
  ! latter includes burned litter. Here we assume that the species distribution in
  ! the litter is the same as in the burned biomass and scale the burned biomass to
  ! match the total rate of carbon emission.
  csmoke_by_sp(:) = max(burned_by_sp, 0.0)
  s = sum(csmoke_by_sp)
  if (s > 0) then
     csmoke_by_sp(:) = csmoke_by_sp(:) / s * csmoke_rate_daily
  else
     ! This could be happening when we do not burn, or burn litter only: then
     ! we take the emission factors from the species of the first cohort (which
     ! always exists and is also the tallest).
     csmoke_by_sp(:) = 0.0
     sp = vegn%cohorts(1)%species
     csmoke_by_sp(sp) = csmoke_rate_daily
  endif

  do tr = 1,n_fire_tr
     ! Use the emission factors for each species and with burned carbon rate for
     ! each species to calculate total emission.
     vegn%fire_emis_land(tr) = 0.0
     do sp = 0,nspecies-1
        vegn%fire_emis_land(tr) = vegn%fire_emis_land(tr) + &
             frdata(tr)%efactors(sp) * &
             csmoke_by_sp(sp) * &
             (1./(spdata(sp)%c_per_dry_matter * 1.E-3)) * & !! convert C to DM in grams
             1.E-4 * &                 !! m2_to_cm2
             (1./(24.*60.*60.)) * &    !! per_second
             (1./frdata(tr)%fire_mw) * &
             AVOGNO
     enddo
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
