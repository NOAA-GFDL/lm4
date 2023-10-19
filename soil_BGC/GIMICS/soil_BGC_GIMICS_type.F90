module soil_BGC_GIMICS_type_mod

use fms_mod, only: input_nml_file, check_nml_error, file_exist, close_file, &
        stdlog, mpp_pe, mpp_root_pe, error_mesg, FATAL, NOTE
use time_manager_mod, only: time_type_to_real
use constants_mod, only : PI,tfreeze

use land_constants_mod, only : N_LITTER_POOLS, seconds_per_year, &
        N_C_TYPES, C_FAST, C_SLOW, C_MIC, LITT_LEAF, LITT_CWOOD

use land_data_mod, only : log_version, lnd
use land_debug_mod, only : land_error_message

use tile_diag_buff_mod, only : diag_buff_type
use tile_diag_base_mod, only : set_default_diag_filter, &
        register_tiled_diag_field, send_tile_data, add_tiled_diag_field_alias, CMOR_NAME

use soil_BGC_type_mod, only : soil_BGC_t, deplete_pool
use soil_BGC_util_mod, only : register_soilc_diag_fields, register_litter_diag_fields, &
        register_litter_soilc_diag_fields
use soil_tile_mod, only : soil_tile_type, num_l, dz, soil_theta, soil_pClay
use vegn_data_mod, only : spdata
use vegn_tile_mod, only : vegn_tile_type
use vegn_cohort_mod, only : cohort_root_litter_profile

implicit none; private

! ---- public items
public :: soil_BGC_GIMICS_t
public :: new_soilc_GIMICS
public :: init_GIMICS_state
public :: read_soil_BGC_GIMICS_namelist, soil_BGC_diag_init_GIMICS
public :: save_equilibration_data ! logical flag triggering writing the data needed for equilibration of soil carbon

! ---- interfces
interface new_soilc_GIMICS
   module procedure soilc_GIMICS_ctor
   module procedure soilc_GIMICS_copy
end interface

! ---- constants
character(len=*), parameter :: module_name = 'soil_BCG_GIMICS_type_mod'
#include "../../shared/version_variable.inc"
real, parameter :: hours_per_year = seconds_per_year/3600.0

! ----  types

! GIMICS BGC pool
type GIMICS_BGC_pool
! concentrations of various carbon pools, [kgC/m3]
! slm: need initial values
    real :: metabolicLitterC  = 0.0
    real :: structuralLitterC = 0.0
    real :: protectedC        = 0.0
    real :: chemResistantC    = 0.0
    real :: availableC        = 0.0
    real :: microbesR         = 0.0
    real :: microbesK         = 0.0

! slm: are these prognostic or for diagnostics only?
    real :: DecompMrLm        = 0.0
    real :: DecompMrLs        = 0.0
    real :: DecompMrCa        = 0.0
    real :: DecompMkLm        = 0.0
    real :: DecompMkLs        = 0.0
    real :: DecompMkCa        = 0.0
    real :: OxidMrCc          = 0.0
    real :: OxidMkCc          = 0.0
    real :: Desorb            = 0.0
    real :: MrTau             = 0.0
    real :: MkTau             = 0.0
    real :: Resp              = 0.0
end type

!> @brief soil carbon data container for GIMICS soil carbon model
type, extends (soil_BGC_t) :: soil_BGC_GIMICS_t
  type(GIMICS_BGC_pool) :: litt (N_LITTER_POOLS) ! surface litter (leaf,coarse wood)
  ! slm: what is the dz that is associated with the surface litter? to get the total carbon, etc.
  type(GIMICS_BGC_pool), allocatable :: &
    rhiz(:),    & ! rhizosphere
    bulk(:)       ! bulk soil (i.e. soil that is not rhizosphere)
  real, allocatable :: fRhiz(:) ! fraction of rhizosphere in each layer, unitless, [0,1]
contains
  procedure :: merge => merge_GIMICS     ! merge another soil carbon tile into current one
  procedure :: total_C => total_C_GIMICS ! returns total C [kgC/m2]
  procedure :: total_N => total_N_GIMICS ! returns total N [kgN/m2]
  procedure :: rav_C => rav_C_GIMICS ! returns amounts of C [kgC/m2]
                                                       ! for legacy surface resistance calculations
  procedure :: get_DOC => get_zero_2D ! returns DOC, by type and by layer [slm: check with Minjin if we have (or can calculate) DOC in GIMICS]
  procedure :: get_DON => get_zero_2D ! returns DON, by type and by layer
  procedure :: get_nit => get_zero_1D ! returns nitrate by layer, kgN/m2
  procedure :: get_amm => get_zero_1D ! returns ammonium by layer, kgN/m2
  procedure :: get_littC => get_littC_GIMICS ! returns litter carbon, by litter pool, kgC/m2

  procedure :: add_soil_matter => add_soil_matter_GIMICS ! add new surface and sub-surface litter to soil carbon and nitrogen
  procedure :: add_root_litter => add_root_litter_GIMICS ! add new root litter to soil carbon and nitrogen
  procedure :: add_root_exudates => add_root_exudates_GIMICS ! add root exudates to soil carbon
  procedure :: burn_litter_frac => burn_litter_frac_GIMICS  ! burn a fraction of sfc litter and retuen amounts of burned carbon and nitrogen
  procedure :: tracer_leaching => tracer_leaching_GIMICS

  procedure :: deposit_N => deposit_N_GIMICS
  procedure :: active_root_N_uptake => active_root_N_uptake_GIMICS
  procedure :: myc_miner_N_uptake => myc_miner_N_uptake_GIMICS
  procedure :: myc_scavenger_N_uptake => myc_scavenger_N_uptake_GIMICS

  procedure :: spend_intermediate_pools => spend_intermediate_pools_GIMICS
  procedure :: dsdt => dsdt_GIMICS
  procedure :: step3 => step3_GIMICS
  procedure :: redistribute_peat_carbon => redistribute_peat_carbon_GIMICS
end type

! ---- module data

real :: delta_time ! fast (physical) time step [s]
real :: dt_fast_yr ! fast (physical) time step [yr] (year is defined as 365 days)
real :: dt_fast_hr ! fast (physical) time step [hr]
real :: dz_litt=0.1  ! slm: temporary, surface litter thickness [m]. To be replaced with dynamic thickness, based ol litter density


! namelist
real :: Vmod_Mr_Lm = 10.0    ! Modifies Vmax for Lm fluxes into Mr (unitless)
real :: Vmod_Mr_Ls = 2.0     ! Modifies Vmax for Ls fluxes into Mr (unitless)
real :: Vmod_Mr_Ca = 10.0    ! Modifies Vmax for Ca fluxes into Mr (unitless)
real :: Vmod_Mk_Lm = 3.0     ! Modifies Vmax for Lm fluxes into Mk (unitless)
real :: Vmod_Mk_Ls = 3.0     ! Modifies Vmax for Ls fluxes into Mk (unitless)
real :: Vmod_Mk_Ca = 2.0     ! Modifies Vmax for Ca fluxes into Mk (unitless)
real :: Vslope     = 0.063   ! Regression coefficient (ln(mgC/mgM/hr)/Celsius) (Eq 1 in Wieder et al., 2015)
real :: Vint       = 5.47    ! Regression intercept (ln(mgC/mgM/hr)) (Eq 1 in Wieder et al., 2015)
real :: aV         = 8e-6            ! Tuning coefficient (unitless) (Eq 1 in Wieder et al., 2015)

real :: Kmod_Mr_Lm = 0.125   ! Modifies Km for Lm fluxes into Mr (unitless)
real :: Kmod_Mr_Ls = 0.5     ! Modifies Km for Ls fluxes into Mr (unitless)
real :: Kmod_Mr_Ca = 0.25    ! Modifies Km for Ca fluxes into Mr (unitless)
real :: Kmod_Mk_Lm = 0.5     ! Modifies Km for Lm fluxes into Mk (unitless)
real :: Kmod_Mk_Ls = 0.25    ! Modifies Km for Ls fluxes into Mk (unitless)
real :: Kmod_Mk_Ca = 0.167   ! Modifies Km for Ca fluxes into Mk (unitless)
real :: Kslope_Lm = 0.017    ! Regression coefficient (ln(mgC/cm3)/Celsius) (Eq 2 in Wieder et al., 2015)
real :: Kslope_Ls = 0.027    ! Regression coefficient (ln(mgC/cm3)/Celsius) (Eq 2 in Wieder et al., 2015)
real :: Kslope_Ca = 0.017    ! Regression coefficient (ln(mgC/cm3)/Celsius) (Eq 2 in Wieder et al., 2015)
real :: Kint      = 3.19     ! Regression intercept (ln(mgC/cm3)) (Eq 2 in Wieder et al., 2015)
real :: aK        = 10.0     ! Tuning coefficient (unitless) (Eq 2 in Wieder et al., 2015)

real :: fI_Lm = 0.38464225   ! Partitioning of litter inputs to Lm (unitless)

real  :: eLm_Mr = 0.55       ! Microbial growth efficiency for fluxes from Lm to Mr (mg/mg)
real  :: eLs_Mr = 0.25       ! Microbial growth efficiency for fluxes from Ls to Mr (mg/mg)
real  :: eCa_Mr = 0.55       ! Microbial growth efficiency for fluxes from Ca to Mr (mg/mg)
real  :: eLm_Mk = 0.75       ! Microbial growth efficiency for fluxes from Lm to Mk (mg/mg)
real  :: eLs_Mk = 0.35       ! Microbial growth efficiency for fluxes from Ls to Mk (mg/mg)
real  :: eCa_Mk = 0.75       ! Microbial growth efficiency for fluxes from Ca to Mk (mg/mg)

real  :: Kmod_oxid_Mr = 4.0  ! Further modifies Km for oxidation of Cc
real  :: Kmod_oxid_Mk = 4.0  ! Further modifies Km for oxidation of Cc

real  :: min_anaerobic_resp_factor = 0.05
real  :: min_dry_resp_factor = 0.05
real  :: gas_diffusion_exp = 2.5 ! Exponent for gas diffusion power law dependence on theta
                             ! See Meslin et al 2010, SSAJ
real  :: substrate_diffusion_exp = 3.0  ! Exponent for theta dependence at low theta.
                             ! See Davison et al DAMM model paper

real  :: r_rhiz = 0.001      ! Radius of rhizosphere around fine root (m)

real  :: init_Mr = 1e-15 ! initial (cold-start) value of microbesR, kg/m3
real  :: init_Mk = 1e-15 ! initial (cold-start) value of microbesR, kg/m3

logical, protected :: save_equilibration_data = .FALSE. !< if TRUE, information for
                         !! soil BGC equilibration acceleration is saved to disk

namelist /soil_BGC_GIMICS_nml/ &
    Vmod_Mr_Lm, Vmod_Mr_Ls, Vmod_Mr_Ca, Vmod_Mk_Lm, Vmod_Mk_Ls, Vmod_Mk_Ca, Vslope, Vint, aV, &
    Kmod_Mr_Lm, Kmod_Mr_Ls, Kmod_Mr_Ca, Kmod_Mk_Lm, Kmod_Mk_Ls, Kmod_Mk_Ca, Kslope_Lm, Kslope_Ls, Kslope_Ca, Kint, aK, &
    fI_Lm, eLm_Mr, eLs_Mr, eCa_Mr, eLm_Mk, eLs_Mk, eCa_Mk, Kmod_oxid_Mr, Kmod_oxid_Mk, &
    min_anaerobic_resp_factor, min_dry_resp_factor, gas_diffusion_exp, substrate_diffusion_exp, &
! -----
    init_Mr, init_Mk, r_rhiz, &
    save_equilibration_data

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
!> @brief read namelist and set up few constants
subroutine read_soil_BGC_GIMICS_namelist()
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, __FILE__)
  read (input_nml_file, nml=soil_BGC_GIMICS_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_BGC_GIMICS_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=soil_BGC_GIMICS_nml)
  endif

  delta_time = time_type_to_real(lnd%dt_fast) ! store in a module variable for convenience
  dt_fast_yr = delta_time/seconds_per_year
  dt_fast_hr = delta_time/3600.0
end subroutine

! ============================================================================
!> @brief Register diagnostic fields
subroutine soil_BGC_diag_init_GIMICS(id_ug, id_zfull)
  integer,intent(in)  :: id_ug    !< Unstructured axis id
  integer,intent(in)  :: id_zfull !< Vertical (depth) axis id

  character(*), parameter :: diag_mod_name = 'soil'
  integer :: axes(2)

  ! define array of axis indices
  axes = [ id_ug,id_zfull ]

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  ! diag field registration goes here
end subroutine

! ============================================================================
!> @brief Create new (empty) soil carbon representation
!! @return Pointer to new soil carbon data structure
function soilc_GIMICS_ctor(soil) result(ptr)
  class(soil_BGC_GIMICS_t), pointer :: ptr
  type(soil_tile_type), intent(in) :: soil !< soil tile data

  allocate(ptr)
  allocate( &
      ptr%rhiz  (num_l), &
      ptr%bulk  (num_l), &
      ptr%fRhiz (num_l)  )
  ptr%fRhiz(:) = 0.0 ! is this reasonable?
end function

! ============================================================================
!> @brief Create a copy of existing soil carbon representation
!! @return Pointer to new soil carbon data structure
function soilc_GIMICS_copy(soilc) result(ptr)
  type(soil_BGC_GIMICS_t), pointer :: ptr
  type(soil_BGC_GIMICS_t), intent(in) :: soilc !< soil carbon data to copy

  allocate(ptr)
  ptr = soilc
end function

! ============================================================================
!> @brief Set initial (cold-start) values to the variables in GIMICS soil BGC data
! perhaps this can be done in constructor? not if the state depends on some other
! state variables not available in constructor, e.g. soil type or soil moisture
subroutine init_GIMICS_state(soilc)
  type(soil_BGC_GIMICS_t), intent(inout) :: soilc !< soil BGC data to initialize

  integer :: k

  do k = 1, num_l
     soilc%fRhiz(k) = 0.0
     soilc%bulk(k)%microbesR = init_Mr
     soilc%bulk(k)%microbesK = init_Mk
     ! the rest of the fields remain at their initial values of zero
  enddo
  do k = 1, N_LITTER_POOLS
     soilc%litt(k)%microbesR = init_Mr
     soilc%litt(k)%microbesK = init_Mk
     ! slm: set litter thickness here
     ! the rest of the fields remain at their initial values of zero
  enddo
end subroutine

! ============================================================================
!> @brief merge s1 into current soil carbon type s2, with given weights
subroutine merge_GIMICS(s2,w2,s1,w1)
  class(soil_BGC_GIMICS_t), intent(inout) :: s2    !< current soil carbon state
  class(soil_BGC_t)       , intent(in)    :: s1    !< soil carbon state to be merged into current
  real                    , intent(in)    :: w2,w1 !< merging weights

  real :: x1, x2 ! normalized relative weights
  real :: f1, f2 ! fractions of rhizosphere or bulk pools, for weight calculations
  real :: y1, y2 ! normalized relative weights for rhizosphere or bulk merges
  integer :: k

  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  select type(s1)
  type is (soil_BGC_GIMICS_t)
     ! merge surface litter pools
     ! slm: this is incorrect, we must take into account difference in surface litter thickness in pools
     do k = 1, N_LITTER_POOLS
        call merge_pools_GIMICS(s2%litt(k),x2,s1%litt(k),x1)
     enddo
     ! merge soil pools, layer by layer
     do k = 1,size(s2%rhiz)
        ! rhizosphere pools
        f1 = s1%fRhiz(k)         ; f2 = s2%fRhiz(k)
        y1 = x1*f1/(x1*f1+x2*f2) ; y2 = 1.0 - y1
        call merge_pools_GIMICS(s2%rhiz(k),y2,s1%rhiz(k),y1)
        ! bulk pools
        f1 = 1.0 - s1%fRhiz(k)   ; f2 = 1.0 - s2%fRhiz(k)
        y1 = x1*f1/(x1*f1+x2*f2) ; y2 = 1.0 - y1
        call merge_pools_GIMICS(s2%bulk(k),y2,s1%bulk(k),y1)
        ! update the rhizosphere fraction
        s2%fRhiz(k) = x1*s1%fRhiz(k) + x2*s2%fRhiz(k)
     enddo

  class default
     call land_error_message('merge_GIMICS: attempt to merge incompatible soil carbon types', FATAL)
  end select
end subroutine

! ============================================================================
!> @brief Merge GIMICS pool p1 into pool p2, with given weights
subroutine merge_pools_GIMICS(p2,w2,p1,w1)
  type(GIMICS_BGC_pool), intent(inout) :: p2
  type(GIMICS_BGC_pool), intent(in)    :: p1
  real,                  intent(in)    :: w2, w1

  real :: x1,x2 ! normalized weights for merging

  ! normalize wights
  x1 = w1/(w1+w2); x2 = 1.0-x1

#define __MERGE__(var) p2%var = x2*p2%var + x1*p1%var
  __MERGE__(metabolicLitterC)
  __MERGE__(structuralLitterC)
  __MERGE__(protectedC)
  __MERGE__(chemResistantC)
  __MERGE__(availableC)
  __MERGE__(microbesR)
  __MERGE__(microbesK)

  __MERGE__(DecompMrLm)
  __MERGE__(DecompMrLs)
  __MERGE__(DecompMrCa)
  __MERGE__(DecompMkLm)
  __MERGE__(DecompMkLs)
  __MERGE__(DecompMkCa)
  __MERGE__(OxidMrCc)
  __MERGE__(OxidMkCc)
  __MERGE__(Desorb)
  __MERGE__(MrTau)
  __MERGE__(MkTau)
  __MERGE__(Resp)
#undef __MERGE__
end subroutine


! ============================================================================
!> @brief Change the rhizosphere fraction in the soil
subroutine set_fRhiz(soilc,fRhiz)
  class(soil_BGC_GIMICS_t), intent(inout) :: soilc !< soil carbon data structure
  real, intent(in) :: fRhiz(:) !< new fractions of rhizosphere, by layer. Unitless, [0,1]

  integer :: k
  real :: wr,wb ! weights for rhizosphere and bulk

  do k = 1, size(soilc%fRhiz)
     if (fRhiz(k) < soilc%fRhiz(k)) then
        ! part of rhizosphere becomes bulk soil
        wr = soilc%fRhiz(k) - fRhiz(k)
        wb = 1.0 - soilc%fRhiz(k)
        call merge_pools_GIMICS(soilc%bulk(k),wb,soilc%rhiz(k),wr)
        soilc%fRhiz(k) = fRhiz(k)
     else if (fRhiz(k) > soilc%fRhiz(k)) then
        ! part of bulk soil becomes rhizosphere
        wb = fRhiz(k) - soilc%fRhiz(k)
        wr = soilc%fRhiz(k)
        call merge_pools_GIMICS(soilc%rhiz(k),wr,soilc%bulk(k),wb)
        soilc%fRhiz(k) = fRhiz(k)
     else
        ! do nothing, rhizosphere fraction did not change (or one of fRhiz is a NaN)
     endif
  enddo
end subroutine

! ============================================================================
!> @brief Given soil carbon state, return total soil C
!! @return Total soil carbon, kgC/m2
real function total_C_GIMICS(soilc) result(tot_C)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc !< soil carbon data structure

  integer :: k

  ! slm: what is dz associated with the surface litter?
  tot_C = 0.0
  do k = 1, N_LITTER_POOLS
     tot_C = tot_C + tot_pool_C(soilc%litt(k)) * dz_litt
  enddo

  do k = 1,num_l
     tot_C = tot_C + &
           ( tot_pool_C(soilc%rhiz(k)) * soilc%fRhiz(k)     &
           + tot_pool_C(soilc%bulk(k)) * (1-soilc%fRhiz(k)) &
           ) * dz(k)
  enddo
end function

! ============================================================================
!> @brief Given soil BGC pool, calculate total volumetric density of carbon, kgC/m3
!! @return Total soil carbon in the pool, kgC/m3
real function tot_pool_C(pool)
  type(GIMICS_BGC_pool), intent(in) :: pool
  tot_pool_C = pool%metabolicLitterC + pool%structuralLitterC &
             + pool%protectedC + pool%chemResistantC &
             + pool%availableC + pool%microbesR + pool%microbesK
end function

! ============================================================================
!> @brief Given soil carbon state, return total soil nitrogen
!! @return total soil nitrogen, kgN/m2
real function total_N_GIMICS(soilc) result(tot_N)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc ! soil carbon data structure
  tot_N = 0.0
end function

! ============================================================================
!> @brief Given soil carbon state, return carbon amount of litter relevant for surface
!! resistance calculations in legacy treatment of soil surface resistance
subroutine rav_C_GIMICS(soilc, fast_C,slow_C,dmic_C)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc !< soil carbon data structure
  real, intent(out) :: &
     fast_C,    & !< fast litter carbon, [kgC/m2]
     slow_C,    & !< slow litter carbon, [kgC/m2]
     dmic_C       !< mass of microbes in litter, [kgC/m2]
! following CORPSE example, we only report the leaf litter pool values
  associate (pool=>soilc%litt(LITT_LEAF))
     fast_C = pool%metabolicLitterC  * dz_litt
     slow_C = pool%structuralLitterC * dz_litt
     dmic_C = (pool%microbesR + pool%microbesR) * dz_litt
  end associate
end subroutine

! ============================================================================
!> @brief Given soil carbon state, returns total litter C per litter pool
subroutine get_littC_GIMICS(soilc, values)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc     !< soil carbon data structure
  real, intent(out)                     :: values(:) !< total C in litter, by litter pool, kgC/m2

  integer :: i
  do i = 1, N_LITTER_POOLS
     values(i) = tot_pool_C(soilc%litt(i))*dz_litt
  enddo
end subroutine

! ============================================================================
! Update the state of the soil BGC pools due to the soil microbiology and other
! processes, and accumulate heterotrophic respiration
subroutine dsdt_GIMICS(soilc, soil, vegn, diag, soilt, theta)
  class(soil_BGC_GIMICS_t)  , intent(inout) :: soilc
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! average soil temperature, deg K, [unused]
  real                , intent(in)    :: theta ! average soil moisture [unused]

  real, dimension(num_l) :: decomp_T, decomp_theta
  real, dimension(num_l) :: rhiz_frac
  real :: clay_frac ! fraction of clay, unitless in interval [0,1]. Should it be by-layer?
  integer :: k

  decomp_theta = soil_theta(soil)
  decomp_T     = soil%T(1:num_l) - tfreeze

  !  First surface litter is decomposed
  do k = 1,N_LITTER_POOLS
     call update_pool_GIMICS(soilc%litt(k), decomp_T(1), decomp_theta(1), fClay=0.0, is_sfc_litter=.TRUE.)
     ! accumulate loss of C to atmosphere [kgC/m2/year]
     vegn%rh=vegn%rh + soilc%litt(k)%Resp*dz_litt*hours_per_year
!      do i = 1, N_C_TYPES
!         call send_tile_data(id_litter_rsoil_C(k,i), litter_C_loss_rate(i), diag)
!         call send_tile_data(id_litter_rsoil_N(k,i), litter_N_loss_rate(i), diag)
!      enddo
     ! for budget check
!      vegn%fsc_out     = vegn%fsc_out     + litter_C_loss_rate(C_FAST)*dt_fast_yr
!      vegn%ssc_out     = vegn%ssc_out     + litter_C_loss_rate(C_SLOW)*dt_fast_yr
!      vegn%deadmic_out = vegn%deadmic_out + litter_C_loss_rate(C_MIC) *dt_fast_yr
  enddo

  ! Next we have to go through layers and decompose the soil carbon pools
  call rhizosphere_frac(vegn, rhiz_frac)
  call set_fRhiz(soilc,rhiz_frac)
  clay_frac = soil_pClay(soil)/100.0
  do k=1,num_l
     call update_pool_GIMICS(soilc%rhiz(k), decomp_T(k), decomp_theta(k), clay_frac, is_sfc_litter=.FALSE.)
     call update_pool_GIMICS(soilc%bulk(k), decomp_T(k), decomp_theta(k), clay_frac, is_sfc_litter=.FALSE.)
     ! accumulate loss of C to atmosphere [kgC/m2/year]
     vegn%rh = vegn%rh + soilc%rhiz(k)%Resp*dz(k)*hours_per_year*soilc%fRhiz(k)      &
                       + soilc%bulk(k)%Resp*dz(k)*hours_per_year*(1-soilc%fRhiz(k))
  enddo
!   do i = 1, N_C_TYPES
!      if (id_rsoil_C(i)>0) call send_tile_data(id_rsoil_C(i), C_loss_rate(:,i)/dz(1:num_l), diag)
!      if (id_rsoil_N(i)>0) call send_tile_data(id_rsoil_N(i), N_loss_rate(:,i)/dz(1:num_l), diag)
!   enddo

  ! ---- diagnostic section
!   call send_tile_data(id_rsoil_C(C_FAST), fast_C_loss(:)/(dz(1:num_l)*dt_fast_yr), diag)
!   call send_tile_data(id_rsoil_C(C_SLOW), slow_C_loss(:)/(dz(1:num_l)*dt_fast_yr), diag)
!   call send_tile_data(id_rsoil, vegn%rh, diag)

  ! TODO: arithmetic averaging of A does not seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that
!   if (id_asoil>0) call send_tile_data(id_asoil, sum(A(:))/size(A(:)), diag)
!   call send_tile_data(id_rh, vegn%rh/seconds_per_year, diag)
end subroutine

! ============================================================================
subroutine step3_GIMICS(soilc, diag)
  class(soil_BGC_GIMICS_t),   intent(inout) :: soilc
  type(diag_buff_type), intent(inout) :: diag

!   integer :: i, k
!
!   associate (soil=>soilc) ! to avoid renaming
!   call send_tile_data(id_fsc, sum(soil%fast_soil_C(:))+sum(soil%litter_SIMPLE_C(C_FAST,:)), diag)
!   call send_tile_data(id_ssc, sum(soil%slow_soil_C(:))+sum(soil%litter_SIMPLE_C(C_SLOW,:)), diag)
!   call send_tile_data(id_soil_C(C_FAST), soil%fast_soil_C(:)/dz(1:num_l), diag)
!   call send_tile_data(id_soil_C(C_SLOW), soil%slow_soil_C(:)/dz(1:num_l), diag)
!   call send_tile_data(id_total_soil_C, sum(soil%fast_soil_C(:))+sum(soil%slow_soil_C(:))+sum(soil%litter_SIMPLE_C(:,:)), diag)
!   do k = 1, N_LITTER_POOLS
!      if (id_litter_total_C(k)>0) call send_tile_data(id_litter_total_C(k), sum(soil%litter_SIMPLE_C(:,k)), diag)
!      do i = 1, N_C_TYPES
!         call send_tile_data(id_litter_C(k,i), soil%litter_SIMPLE_C(i,k), diag)
!      enddo
!   enddo
!
!   ! --- CMOR vars
!   if (id_csoilfast>0)   call send_tile_data(id_csoilfast,   sum(soil%fast_soil_C(:)), diag)
!   if (id_csoilmedium>0) call send_tile_data(id_csoilmedium, sum(soil%slow_soil_C(:)), diag)
!   call send_tile_data(id_csoilslow, 0.0, diag)
!   if (id_csoil>0)       call send_tile_data(id_csoil, sum(soil%fast_soil_C(:))+sum(soil%slow_soil_C(:)), diag)
!   if (id_cSoilLevels>0) call send_tile_data(id_cSoilLevels, soil%fast_soil_C(:)+soil%slow_soil_C(:), diag)
!   if (id_cLitter>0)     call send_tile_data(id_cLitter, sum(soil%litter_SIMPLE_C(:,:)), diag)
!   if (id_cLitterCwd>0)  call send_tile_data(id_cLitterCwd, sum(soil%litter_SIMPLE_C(:,LITT_CWOOD)), diag)
!   if (id_cLitterLeaf>0) call send_tile_data(id_cLitterLeaf, sum(soil%litter_SIMPLE_C(:,LITT_LEAF)), diag)
!   ! --- end of CMOR vars
!   end associate

end subroutine

! ============================================================================
!> @brief Update soil carbon pool
subroutine update_pool_GIMICS(pool, T, theta, fClay, is_sfc_litter)
  type(GIMICS_BGC_pool),intent(inout) :: pool
  real,    intent(in) :: T         !< Temperature [degC]
  real,    intent(in) :: theta     !< volumetric water content slm: [per unit soil volume, or per unit pore volume?]
  real,    intent(in) :: fClay     !< clay fraction, unitless, within [0,1] interval
  logical, intent(in) :: is_sfc_litter !< TRUE is the pool is surface litter: protected C is always zero in this case

  real:: Vmax_Mr_Lm, Vmax_Mr_Ls, Vmax_Mr_Ca, Vmax_Mk_Lm, Vmax_Mk_Ls, Vmax_Mk_Ca, &
         Km_Mr_Lm,   Km_Mr_Ls,   Km_Mr_Ca,   Km_Mk_Lm,   Km_Mk_Ls,   Km_Mk_Ca,   &
         fMrTau_Cp, fMkTau_Cp, fMrTau_Cc, fMkTau_Cc

  Vmax_Mr_Lm = theta_func(theta,1-theta,substrate_diffusion_exp,gas_diffusion_exp,min_anaerobic_resp_factor, min_dry_resp_factor) * exp(Vslope*T+Vint) * aV * Vmod_Mr_Lm ! mgC/mgM/h
  Vmax_Mr_Ls = theta_func(theta,1-theta,substrate_diffusion_exp,gas_diffusion_exp,min_anaerobic_resp_factor, min_dry_resp_factor) * exp(Vslope*T+Vint) * aV * Vmod_Mr_Ls
  Vmax_Mr_Ca = theta_func(theta,1-theta,substrate_diffusion_exp,gas_diffusion_exp,min_anaerobic_resp_factor, min_dry_resp_factor) * exp(Vslope*T+Vint) * aV * Vmod_Mr_Ca

  Vmax_Mk_Lm = theta_func(theta,1-theta,substrate_diffusion_exp,gas_diffusion_exp,min_anaerobic_resp_factor, min_dry_resp_factor) * exp(Vslope*T+Vint) * aV * Vmod_Mk_Lm
  Vmax_Mk_Ls = theta_func(theta,1-theta,substrate_diffusion_exp,gas_diffusion_exp,min_anaerobic_resp_factor, min_dry_resp_factor) * exp(Vslope*T+Vint) * aV * Vmod_Mk_Ls
  Vmax_Mk_Ca = theta_func(theta,1-theta,substrate_diffusion_exp,gas_diffusion_exp,min_anaerobic_resp_factor, min_dry_resp_factor) * exp(Vslope*T+Vint) * aV * Vmod_Mk_Ca

  Km_Mr_Lm = exp(Kslope_Lm*T+Kint) * aK * Kmod_Mr_Lm ! kgC/m3
  Km_Mr_Ls = exp(Kslope_Ls*T+Kint) * aK * Kmod_Mr_Ls
  Km_Mr_Ca = exp(Kslope_Ca*T+Kint) * aK * Kmod_Mr_Ca / (2.0*exp(-2.0*sqrt(fClay)))

  Km_Mk_Lm = exp(Kslope_Lm*T+Kint) * aK * Kmod_Mk_Lm
  Km_Mk_Ls = exp(Kslope_Ls*T+Kint) * aK * Kmod_Mk_Ls
  Km_Mk_Ca = exp(Kslope_Ca*T+Kint) * aK * Kmod_Mk_Ca / (2.0*exp(-2.0*sqrt(fClay)))



  pool%DecompMrLm = Vmax_Mr_Lm * pool%microbesR * (pool%metabolicLitterC/(Km_Mr_Lm+pool%metabolicLitterC)) ! kgC/m3/h
  pool%DecompMrLs = Vmax_Mr_Ls * pool%microbesR * (pool%structuralLitterC/(Km_Mr_Ls+pool%structuralLitterC))
  pool%DecompMrCa = Vmax_Mr_Ca * pool%microbesR * (pool%availableC/(Km_Mr_Ca+pool%availableC))

  pool%DecompMkLm = Vmax_Mk_Lm * pool%microbesK * (pool%metabolicLitterC/(Km_Mk_Lm+pool%metabolicLitterC))
  pool%DecompMkLs = Vmax_Mk_Ls * pool%microbesK * (pool%structuralLitterC/(Km_Mk_Ls+pool%structuralLitterC))
  pool%DecompMkCa = Vmax_Mk_Ca * pool%microbesK * (pool%availableC/(Km_Mk_Ca+pool%availableC))



  pool%Resp = (1-eLm_Mr)*pool%DecompMrLm + (1-eLs_Mr)*pool%DecompMrLs + (1-eCa_Mr)*pool%DecompMrCa + & ! kgC/m3/h
              (1-eLm_Mk)*pool%DecompMkLm + (1-eLs_Mk)*pool%DecompMkLs + (1-eCa_Mk)*pool%DecompMkCa



  pool%OxidMrCc = Vmax_Mr_Ls * pool%microbesR * (pool%chemResistantC/(Kmod_oxid_Mr*Km_Mr_Ls+pool%chemResistantC)) ! kgC/m3/h
  pool%OxidMkCc = Vmax_Mk_Ls * pool%microbesK * (pool%chemResistantC/(Kmod_oxid_Mk*Km_Mk_Ls+pool%chemResistantC))


  pool%MrTau = 5.2e-4 * exp(0.3*fI_Lm) * pool%microbesR ! kgC/m3/h
  pool%MkTau = 2.4e-4 * exp(0.1*fI_Lm) * pool%microbesK


  if (is_sfc_litter) then
     fMrTau_Cp = 0.0
     fMkTau_Cp = 0.0
  else
     fMrTau_Cp = 0.3*exp(1.3*fClay) ! 0.3646 unitless
     fMkTau_Cp = 0.2*exp(0.8*fClay) ! 0.2255
  endif
  fMrTau_Cc = 0.1*exp(-3*fI_Lm)  ! 0.0315
  fMkTau_Cc = 0.3*exp(-3*fI_Lm)  ! 0.0946
  !print *, fMrTau_Cp,fMkTau_Cp,fMrTau_Cc,fMkTau_Cc



  !Km_Mk_Ca=(5.2e-4 * exp(0.3*fI_Lm) * sqrt(NPP_bulk/100.0))
  !print *, Km_Mk_Ca
  !Km_Mk_Ca=(2.4e-4 * exp(0.1*fI_Lm) * sqrt(NPP_bulk/100.0))
  !print *, Km_Mk_Ca

  !Km_Mk_Ca=(1.5e-5*exp(-1.5*fClay))
  !print *, Km_Mk_Ca
  !Km_Mk_Ca = exp(Kslope_Ca*T+Kint) * aK * Kmod_Mk_Ca / (2.0*exp(-2.0*sqrt(fClay)))


  pool%metabolicLitterC  = pool%metabolicLitterC  - (pool%DecompMrLm+pool%DecompMkLm)*dt_fast_hr ! kgC/m3
  pool%structuralLitterC = pool%structuralLitterC - (pool%DecompMrLs+pool%DecompMkLs)*dt_fast_hr

  if (is_sfc_litter) then
     pool%Desorb = 0.0
     pool%protectedC = 0.0
  else
     pool%Desorb = (1.5e-5*exp(-1.5*fClay))*pool%protectedC ! kgC/m3/h
     pool%protectedC = pool%protectedC + (fMrTau_Cp*pool%MrTau + fMkTau_Cp*pool%MkTau - pool%Desorb)*dt_fast_hr ! kgC/m3
  endif

  pool%chemResistantC = pool%chemResistantC + (fMrTau_Cc*pool%MrTau + fMkTau_Cc*pool%MkTau - pool%OxidMrCc - pool%OxidMkCc)*dt_fast_hr ! gC/m3

  pool%availableC = pool%availableC + ((1-fMrTau_Cp-fMrTau_Cc)*pool%MrTau + (1-fMkTau_Cp-fMkTau_Cc)*pool%MkTau +& ! kgC/m3
                    pool%Desorb + pool%OxidMrCc + pool%OxidMkCc - pool%DecompMrCa - pool%DecompMkCa)*dt_fast_hr

  pool%microbesR = pool%microbesR + (eLm_Mr*pool%DecompMrLm + eLs_Mr*pool%DecompMrLs + eCa_Mr*pool%DecompMrCa - pool%MrTau)*dt_fast_hr ! gC/m3

  pool%microbesK = pool%microbesK + (eLm_Mk*pool%DecompMkLm + eLs_Mk*pool%DecompMkLs + eCa_Mk*pool%DecompMkCa - pool%MkTau)*dt_fast_hr
end subroutine

! ============================================================================
! note that in resp_denitrif the dependence on soil moisture is subtly different
real function theta_func(water_filled_porosity,air_filled_porosity,substrate_diffusion_exp,gas_diffusion_exp,min_anaerobic_resp_factor, min_dry_resp_factor)
  real, intent(in) :: water_filled_porosity ! fraction of pores filled with water
  real, intent(in) :: air_filled_porosity ! fraction of pores filled with water
  real, intent(in) :: substrate_diffusion_exp,gas_diffusion_exp,min_dry_resp_factor,min_anaerobic_resp_factor
  real::theta_resp_max,aerobic_max

  theta_resp_max=substrate_diffusion_exp/(gas_diffusion_exp*(1.0+substrate_diffusion_exp/gas_diffusion_exp))
  aerobic_max=theta_resp_max**substrate_diffusion_exp*(1.0-theta_resp_max)**gas_diffusion_exp

  ! Functional dependence on soil moisture, normalized so max is 1
  theta_func=(max(water_filled_porosity,0.0)**substrate_diffusion_exp)*(max(air_filled_porosity,0.0)**gas_diffusion_exp)/aerobic_max
  !theta_func=max(theta_func, 0.05)
  ! On the wet side of the function, make sure it does not go below min_anaerobic_resp_factor
  if(water_filled_porosity>theta_resp_max) theta_func=max(theta_func, min_anaerobic_resp_factor)
  ! On the dry side of the function, make sure it does not go below min_dry_resp_factor
  if(water_filled_porosity<theta_resp_max) theta_func=max(theta_func, min_dry_resp_factor)
end function theta_func

! ============================================================================
!> @brief Add exudates to the soil
subroutine add_root_exudates_GIMICS(soilc, exudateC, exudateN, ammonium, nitrate)
  class(soil_BGC_GIMICS_t), intent(inout)  :: soilC !< soil carbon data structure
  real,intent(in)           :: exudateC(:) !< (num_l) amount of C in exudate, kgC/m2 per layer
  real,intent(in), optional :: exudateN(:) !< (num_l) amount of N in exudate, kgN/m2 per layer
  real,intent(in), optional :: ammonium(:) !< (num_l) amount of ammonium in exudate, kgN/m2(?) per layer
  real,intent(in), optional :: nitrate (:) !< (num_l) amount of  nitrate in exudate, kgN/m2(?) per layer

!   real, dimension(size(soilc%org_matter)) :: NH4,NO3
  integer :: k

!   NH4(:)=0.0
!   NO3(:)=0.0
!   if(present(ammonium)) NH4=ammonium
!   if(present(nitrate))  NO3=nitrate

  do k=1,num_l
     soilC%rhiz(k)%metabolicLitterC = soilc%rhiz(k)%metabolicLitterC + exudateC(k)/dz(k) ! kgC/m3 slm: need factor in rhizosphere fraction
  enddo
end subroutine


! ============================================================================
!> @brief Add new root litter to soil carbon and nitrogen
subroutine add_root_litter_GIMICS(soilC, vegn, litterC, litterN)
  class(soil_BGC_GIMICS_t) , intent(inout) :: soilC !< soil carbon state
  type(vegn_tile_type)  , intent(in)    :: vegn !< vegetation state (for rhizosphere fraction calculation)
  real                  , intent(in)    :: litterC(:,:) !< new litter carbon content (num_l,N_C_TYPES), kgC/m2 of soil layer
  real                  , intent(in)    :: litterN(:,:) !< new litter nitrogen content kgN/m2 of soil layer

  call land_error_message('add_root_litter_GIMICS called -- this should never happen', FATAL)
end subroutine

! ============================================================================
subroutine add_soil_matter_GIMICS(soilc, vegn, &
        leaf_litter_C, wood_litter_C, root_litter_C, &
        leaf_litter_N, wood_litter_N, root_litter_N  )
  class(soil_BGC_GIMICS_t), intent(inout) :: soilc
  type(vegn_tile_type),  intent(inout) :: vegn
  real, intent(in), optional :: leaf_litter_C(:)   ! (N_C_TYPES) kgC/m2
  real, intent(in), optional :: wood_litter_C(:)   ! (N_C_TYPES) kgC/m2
  real, intent(in), optional :: root_litter_C(:,:) ! (num_l,N_C_TYPES) kgC/m2
  real, intent(in), optional :: leaf_litter_N(:)   ! (N_C_TYPES) kgN/m2
  real, intent(in), optional :: wood_litter_N(:)   ! (N_C_TYPES) kgN/m2
  real, intent(in), optional :: root_litter_N(:,:) ! (num_l,N_C_TYPES) kgN/m2

  integer :: k
  real :: leaf_litt_C(N_C_TYPES), leaf_litt_N(N_C_TYPES)
  real :: wood_litt_C(N_C_TYPES), wood_litt_N(N_C_TYPES)
  real :: root_litt_C(size(soilc%bulk),N_C_TYPES), &
          root_litt_N(size(soilc%bulk),N_C_TYPES)

  ! define values for optional arguments that may not be present
  if (present(leaf_litter_C)) then
     leaf_litt_C(:) = leaf_litter_C(:)
  else
     leaf_litt_C(:) = 0.0
  endif
  if (present(wood_litter_C)) then
     wood_litt_C(:) = wood_litter_C(:)
  else
     wood_litt_C(:) = 0.0
  endif
  if (present(root_litter_C)) then
     root_litt_C(:,:) = root_litter_C(:,:)
  else
     root_litt_C(:,:) = 0.0
  endif
  if (present(leaf_litter_N)) then
     leaf_litt_N(:) = leaf_litter_N(:)
  else
     leaf_litt_N(:) = 0.0
  endif
  if (present(wood_litter_N)) then
     wood_litt_N(:) = wood_litter_N(:)
  else
     wood_litt_N(:) = 0.0
  endif
  if (present(root_litter_N)) then
     root_litt_N(:,:) = root_litter_N(:,:)
  else
     root_litt_N(:,:) = 0.0
  endif

  call add_matter_GIMICS1(soilc%litt(LITT_LEAF),  dz_litt, leaf_litt_C, leaf_litt_N)
  call add_matter_GIMICS1(soilc%litt(LITT_CWOOD), dz_litt, wood_litt_C, wood_litt_N)
  do k = 1,size(soilc%bulk)
     call add_matter_GIMICS2(soilc%bulk(k), soilc%rhiz(k), dz(k), soilc%fRhiz(k), root_litt_C(k,:), root_litt_N(k,:))
  enddo

  ! accumulate litterfall diagnostics: it is sent to diag and then reset at every time step
  vegn%litterfall_C(:,LITT_LEAF)  = vegn%litterfall_C(:,LITT_LEAF)  + leaf_litt_C(:)
  vegn%litterfall_C(:,LITT_CWOOD) = vegn%litterfall_C(:,LITT_CWOOD) + wood_litt_C(:)
end subroutine

! ============================================================================
! transfer fraction of intermediate pools (used for smoothing out contributions
! of spiky processes, e.g. harvesting) defined in vegetation tile data structure
! to soil BGC pools.
subroutine spend_intermediate_pools_GIMICS(soilc, vegn)
  class(soil_BGC_GIMICS_t), intent(inout) :: soilc
  type(vegn_tile_type),     intent(inout) :: vegn

  integer :: i,k
  real :: deltafast, deltaslow
  real :: profile(num_l), profile1(num_l), psum ! for deposition profile calculation
  real :: litterC(num_l,N_C_TYPES) ! soil litter C input by layer and type
  real :: delta_C(N_C_TYPES,N_LITTER_POOLS)

  ! ---- surface
  vegn%litter_rate_C = MAX(0.0, MIN(vegn%litter_rate_C, vegn%litter_buff_C/dt_fast_yr))
  delta_C = vegn%litter_rate_C*dt_fast_yr

  do i = 1,N_LITTER_POOLS
     call add_matter_GIMICS1(soilc%litt(i), dz_litt, C=delta_C(:,i))
  enddo
  vegn%litter_buff_C = vegn%litter_buff_C - delta_C
  ! for litterfall diagnostics
  vegn%litterfall_C(:,:) = vegn%litterfall_C(:,:) + delta_C(:,:)

  ! ---- underground
  deltafast = 0.0; call deplete_pool(vegn%fsc_pool_bg, vegn%fsc_rate_bg, deltafast)
  deltaslow = 0.0; call deplete_pool(vegn%ssc_pool_bg, vegn%ssc_rate_bg, deltaslow)

  ! vertical profile of litter is proportional to the average of litter profiles
  ! of all cohorts, weighted with biomasses of fine roots. This does not seem to
  ! be a very good assumption, since fine roots sometimes die (mass is zero),
  ! but profile should not be zero in this case.
  profile(:) = 0.0
  do i = 1,vegn%n_cohorts
     associate(cc=>vegn%cohorts(i))
     call cohort_root_litter_profile(cc,dz,profile1)
     profile(:) = profile(:) + profile1(:)*cc%br*cc%nindivs
     end associate
  enddo
  psum = sum(profile)
  if (psum>0) then
     profile(:) = profile(:)/psum
  else
     profile(:) = 0.0
     profile(1) = 1.0
  endif
  do k = 1,num_l
     call add_matter_GIMICS2(soilc%bulk(k), soilc%rhiz(k), dz(k), soilc%fRhiz(k), C=[deltafast,deltaslow,0.0] * profile(k))
  enddo
end subroutine

! add carbon (and later nitrogen) to GIMICS soil BGC pool
subroutine add_matter_GIMICS1(pool, dz, C, N)
  type(GIMICS_BGC_pool), intent(inout) :: pool ! BGC pool to update
  real, intent(in)  :: dz                      ! layer thickness, m
  real, intent(in), optional :: C (N_C_TYPES)  ! (fast,slow,[dead]microbial), kgC/m2
  real, intent(in), optional :: N (N_C_TYPES)  ! (fast,slow,[dead]microbial), kgN/m2

  if (present(C)) then
     pool%metabolicLitterC  = pool%metabolicLitterC  + C(C_FAST)/dz + C(C_MIC)/dz ! kgC/m3
     pool%structuralLitterC = pool%structuralLitterC + C(C_SLOW)/dz
  endif
!   if (present(N)) then
!      ....
!   endif
end subroutine

! add carbon (and later nitrogen) to GIMICS soil BGC pool, distributing it between bulk soil and rhizosphere
subroutine add_matter_GIMICS2(bulk, rhiz, dz, rhiz_frac, C, N)
  type(GIMICS_BGC_pool), intent(inout) :: bulk, rhiz ! bulk soil and rhizosphere BGC pools, respectively
  real, intent(in) :: dz        ! layer thickness, m
  real, intent(in) :: rhiz_frac ! fraction of soil in rhizosphere
  real, intent(in), optional :: C (N_C_TYPES)  ! (fast,slow,[dead]microbial), kgC/m2
  real, intent(in), optional :: N (N_C_TYPES)  ! (fast,slow,[dead]microbial), kgN/m2

  if (present(C)) then
     bulk%metabolicLitterC  = bulk%metabolicLitterC  + (C(C_FAST)+ C(C_MIC))*(1-rhiz_frac)/dz ! kgC/m3
     bulk%structuralLitterC = bulk%structuralLitterC + C(C_SLOW)*(1-rhiz_frac)/dz

     rhiz%metabolicLitterC  = rhiz%metabolicLitterC  + (C(C_FAST)+ C(C_MIC))*(1-rhiz_frac)/dz ! kgC/m3
     rhiz%structuralLitterC = rhiz%structuralLitterC + C(C_SLOW)*(1-rhiz_frac)/dz
  endif
!   if (present(N)) then
!      ....
!   endif
end subroutine

! ============================================================================
subroutine redistribute_peat_carbon_GIMICS(soilc)
  class(soil_BGC_GIMICS_t), intent(inout) :: soilc

  call error_mesg('redistribute_peat_carbon_GIMICS','not implemented; should it be?', FATAL)
end subroutine

! ============================================================================
! Deposition of nitrogen
subroutine deposit_N_GIMICS(soilc, NH4, NO3, N_org)
  class(soil_BGC_GIMICS_t), intent(inout) :: soilc
  real, intent(in) :: NH4, NO3, N_org ! amounts of NH4, NO3, and organic nitrogen to deposit, kg N/m2
  ! do nothing now: nitrogen deposition is ignored
end subroutine

! ============================================================================
! Nitrogen uptake from the rhizosphere by roots (active transport across root-soil interface)
subroutine active_root_N_uptake_GIMICS(soilc, vegn, N_uptake, dt, update_pools)
  class(soil_BGC_GIMICS_t), intent(inout) :: soilc
  type(vegn_tile_type), intent(in)    :: vegn
  real,    intent(out) :: N_uptake(:) ! Nitrogen uptake, kg N per individual
  real,    intent(in)  :: dt ! in years
  logical, intent(in)  :: update_pools

  N_uptake       = 0.0
end subroutine

! ============================================================================
! Uptake of mineral N by mycorrhizal "scavengers" -- Should correspond to Arbuscular mycorrhizae
subroutine myc_scavenger_N_uptake_GIMICS(soilc,vegn,N_uptake_cohorts,myc_efficiency,dt,update_pools)
  class(soil_BGC_GIMICS_t),  intent(inout) :: soilc
  type(vegn_tile_type), intent(in) :: vegn
  real,intent(out) :: N_uptake_cohorts(:) ! Units: kgN/m2 per individual
  real, intent(in) :: dt  ! dt in years
  logical, intent(in) :: update_pools
  real, intent(out) :: myc_efficiency ! units: kgN/kg myc biomass C. Should give N uptake efficiency even when myc biomass is zero

  N_uptake_cohorts = 0.0
  myc_efficiency   = 0.0
end subroutine

! ============================================================================
! Uptake of mineral N by mycorrhizal "miners" -- Should correspond to Ecto mycorrhizae
subroutine myc_miner_N_uptake_GIMICS(soilc,soil,vegn,N_uptake_cohorts,C_uptake_cohorts,total_CO2prod,myc_efficiency,dt,update_pools)
  class(soil_BGC_GIMICS_t), intent(inout) :: soilc
  type(soil_tile_type), intent(in) :: soil
  type(vegn_tile_type), intent(in) :: vegn
  real,    intent(out) :: N_uptake_cohorts(:), C_uptake_cohorts(:)  ! Units kg/m2 of per individual
  real,    intent(out) :: total_CO2prod ! Units of kgC/m2 (not per individual)
  real,    intent(in)  :: dt  ! dt in years
  logical, intent(in)  :: update_pools
  real,    intent(out) :: myc_efficiency  ! units: kgN/kg myc biomass C. Should give N uptake efficiency even when myc biomass is zero
end subroutine

! ============================================================================
!> @brief Burn given fraction of litter and return amounts of burned C and N
subroutine burn_litter_frac_GIMICS(soilc, frac, burned_C, burned_N)
  class(soil_BGC_GIMICS_t), intent(inout) :: soilc !< soil carbon data structure
  real, intent(in)  :: frac(:)            !< (N_LITTER_POOLS) fraction of litter to burn [0,1]
  real, intent(out) :: burned_C, burned_N !< amounts of burned carbon and nitrogen

  integer :: i

  burned_C = 0.0; burned_N = 0.0
  do i = 1,N_LITTER_POOLS
     associate (pool=>soilc%litt(i))
     burned_C = burned_C + tot_pool_C(pool)*frac(i)
     pool%metabolicLitterC  = frac(i) * pool%metabolicLitterC
     pool%structuralLitterC = frac(i) * pool%structuralLitterC
     pool%protectedC        = frac(i) * pool%protectedC
     pool%chemResistantC    = frac(i) * pool%chemResistantC
     pool%availableC        = frac(i) * pool%availableC
     pool%microbesR         = frac(i) * pool%microbesR
     pool%microbesK         = frac(i) * pool%microbesK
     end associate
  enddo
end subroutine

! ============================================================================
subroutine tracer_leaching_GIMICS(soilc, diag, &
         wl, flow, div, &
         div_hlsp_DOC,div_hlsp_DON,&
         div_hlsp_NO3,div_hlsp_NH4,&
         ! output
         total_DOC_div, total_DON_div, total_NO3_div, total_NH4_div)

  class(soil_BGC_GIMICS_t), intent(inout) :: soilC
  type(diag_buff_type), intent(inout) :: diag
  !!xz check the unit of flow!!For CH's code, it should be kg/year or kg/delta_time unit.!!! I assume here the unit is mm/yr
  real, intent(in) :: flow(:), div(:), wl(:) ! flow (into layer) and wl in units of mm, downward is >0  !!!xz check the unit of dz (should be m in this subroutine), flow (shoul be mm)
  real, intent(in) :: div_hlsp_DOC(:,:) ! dim(N_C_TYPES, num_l) [kg C/m^2/s] net divergence loss from tile calculated in hlsp_hydrology
  real, intent(in) :: div_hlsp_DON(:,:) ! dim(N_C_TYPES, num_l) [kg N/m^2/s] net divergence
  real, intent(in) :: div_hlsp_NO3(:),div_hlsp_NH4(:) ! dim(num_l) [kg N/m^2/s] net divergence loss from tile calculated in hlsp_hydrology

  real, intent(out) :: total_DOC_div, total_DON_div, total_NO3_div, total_NH4_div

  total_DOC_div = 0.0
  total_DON_div = 0.0
  total_NO3_div = 0.0
  total_NH4_div = 0.0
end subroutine

! ============================================================================
subroutine get_zero_2D(soilC, values)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc ! soil carbon data structure
  real,                intent(out) :: values(:,:) ! (N_C_TYPES, num_l) ! [kg C/m^2] dissolved organic carbon

  values(:,:) = 0.0
end subroutine

! ============================================================================
!> @brief return zeros in 1D array
!! This function is used to retrieve substances and values that are not present in
!! soil carbon models, e.g. nitrate or ammonium when the nitrogen dynamics is not
!! implemented
subroutine get_zero_1D(soilC, values)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc !< soil carbon data structure (unused)
  real,                intent(out) :: values(:) !< returned values
  values(:) = 0.0
end subroutine

! ============================================================================
! slm: it is the same as in CORPSE. Should we move it to the vegetation modules,
! or keep it here to be able to make changes independent from CORPSE implementation?

!> @brief Calculate volumetric fraction of rhizosphere in each layer
subroutine rhizosphere_frac(vegn, rFrac)
  type(vegn_tile_type), intent(in)  :: vegn !< vegetation state
  real                , intent(out) :: rFrac(:)!< volumentric fraction of rhizosphere

!   rFrac = rhiz_frac
! slm: perhaps we should have a possibility to use constant rhizosphere fraction?
  real :: rhiz_vol(num_l)  ! volume of rhizosphere in each layer, m3/m3
  integer :: i

  ! first calculate the volume of rhizosphere
  rhiz_vol(:) = 0.0
  do i = 1,vegn%n_cohorts
     associate(cc=>vegn%cohorts(i),sp=>spdata(vegn%cohorts(i)%species))
     rhiz_vol(:) = rhiz_vol(:) + &
         PI*((r_rhiz+sp%root_r)**2-sp%root_r**2)*cc%root_length(1:num_l)*cc%nindivs
     end associate
  enddo
  rFrac(1:num_l) = max(0.0, min(1.0,rhiz_vol(:)))
end subroutine rhizosphere_frac

end module
