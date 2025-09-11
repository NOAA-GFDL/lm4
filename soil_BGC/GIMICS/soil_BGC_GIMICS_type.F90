module soil_BGC_GIMICS_type_mod

#include "../../shared/debug.inc"


use fms_mod, only: input_nml_file, check_nml_error, &
        stdlog, mpp_pe, mpp_root_pe, error_mesg, FATAL, NOTE, string, lowercase
use time_manager_mod, only: time_type, time_type_to_real
use constants_mod, only : PI,tfreeze

use land_constants_mod, only : N_LITTER_POOLS, l_diagname, c_diagname, c_longname, &
        seconds_per_year, N_C_TYPES, C_FAST, C_SLOW, C_MIC, LITT_LEAF, LITT_CWOOD, &
        MAX_SOIL_LEV

use land_data_mod, only : log_version, lnd
use land_debug_mod, only : land_error_message, check_var_range, check_conservation, &
        carbon_cons_tol, is_watch_point, is_watch_cell

use tile_diag_buff_mod, only : diag_buff_type
use tile_diag_base_mod, only : set_default_diag_filter, &
        register_tiled_diag_field, send_tile_data, add_tiled_diag_field_alias, CMOR_NAME, &
        CMOR_1M_DEPTH

use soil_BGC_type_mod, only : soil_BGC_t, deplete_pool, tracer_advection
use soil_BGC_util_mod, only : register_soilc_diag_fields, register_litter_diag_fields, &
        register_litter_soilc_diag_fields
use soil_tile_mod, only : soil_tile_type, gw_option, GW_TILED, &
         num_l, dz, zhalf, zfull, soil_theta, soil_porosity, soil_pClay, SOIL_ICE_POROSITY, soil_water_ice_porosity
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

! ---- interfaces
interface new_soilc_GIMICS
   module procedure soilc_GIMICS_ctor
   module procedure soilc_GIMICS_copy
end interface

! ---- constants
character(len=*), parameter :: module_name = 'soil_BCG_GIMICS_type_mod'
#include "../../shared/version_variable.inc"
real, parameter :: hours_per_year = seconds_per_year/3600.0

! ----  types

! GIMICS soil BGC pool
type GIMICS_BGC_pool
! concentrations of various carbon pools, [kgC/m3]
    real :: metabolicLitterC  = 0.0
    real :: structuralLitterC = 0.0
    real :: protectedC        = 0.0
    real :: chemResistantC    = 0.0
    real :: availableC        = 0.0
    real :: microbesR         = 0.0
    real :: microbesK         = 0.0
    real :: DOC               = 0.0

! tendencies, for diagnostics, [kgC/m3/hr]
    real :: DecompMrLm        = 0.0
    real :: DecompMrLs        = 0.0
    real :: DecompMrCa        = 0.0
    real :: DecompMrDOC       = 0.0
    real :: DecompMkLm        = 0.0
    real :: DecompMkLs        = 0.0
    real :: DecompMkCa        = 0.0
    real :: DecompMkDOC       = 0.0
    real :: OxidMrCc          = 0.0
    real :: OxidMkCc          = 0.0
    real :: Desorb            = 0.0
    real :: MrTau             = 0.0
    real :: MkTau             = 0.0
    real :: Resp              = 0.0
    real :: thetaF            = 0.0 ! soil moisture factor for decomposition [for diagnostics]
    ! since deposition of matter can come from multiple processes, for the diagnostics of tendencies
    ! we accumulate the deposition over time step, then send the data and reset accumulators to zero
    real :: InputStrC = 0.0 ! deposition of structural litter C during time step, kgC/m3/hr
    real :: InputMtbC = 0.0 ! deposition of metabolic litter C during time step, kgC/m3/hr
    real :: InputExdC = 0.0 ! deposition of exudate C during time step, kgC/m3/hr
end type

! GIMICS BGC surface litter pool: it is the same as the soil pool data structure,
! except adds variable "dz" to track the evolution of surface litter thickness.
type, extends(GIMICS_BGC_pool) :: GIMICS_BGC_litt
    real :: dz = 0.0 ! litter thickness, [m]
end type

!> @brief soil carbon data container for GIMICS soil carbon model
type, extends (soil_BGC_t) :: soil_BGC_GIMICS_t
  type(GIMICS_BGC_litt) :: litt (N_LITTER_POOLS) ! surface litter (leaf,coarse wood)
  type(GIMICS_BGC_pool), allocatable :: &
    rhiz(:),    & ! rhizosphere
    bulk(:)       ! bulk soil (i.e. soil that is not rhizosphere)
  real, allocatable :: fRhiz(:) ! fraction of rhizosphere in each layer, unitless, [0,1]

  real :: neg_litt_C(N_C_TYPES) = 0.0 ! cumulative value of negative C litter input to soil

contains
  procedure :: merge => merge_GIMICS     ! merge another soil carbon tile into current one
  procedure :: total_C => total_C_GIMICS ! returns total C [kgC/m2]
  procedure :: total_N => total_N_GIMICS ! returns total N [kgN/m2]
  procedure :: total_soil_C => total_soil_C_GIMICS ! returns total C in soil, excluding surface litter [kgC/m2]
  procedure :: total_soil_N => total_soil_N_GIMICS ! returns total N in soil, excluding surface litter [kgN/m2]
  procedure :: total_soil_C_to_depth => total_soil_C_to_depth_GIMICS ! returns total C in soil from the
                                         ! surface to the specified depth [kgC/m2]

  procedure :: rav_C => rav_C_GIMICS ! returns amounts of C [kgC/m2]
                                                       ! for legacy surface resistance calculations
  procedure :: get_DOC => get_DOC_GIMICS ! returns DOC, by type and by layer
  procedure :: get_DON => get_zero_2D ! returns DON, by type and by layer
  procedure :: get_layer_C => totC_by_layer_GIMICS ! returns total soil carbon by layer, kgC/m2
  procedure :: get_nit => get_zero_1D ! returns nitrate by layer, kgN/m2
  procedure :: get_amm => get_zero_1D ! returns ammonium by layer, kgN/m2
  procedure :: get_littC => get_littC_GIMICS ! returns litter carbon, by litter pool, kgC/m2

  procedure :: add_soil_matter => add_soil_matter_GIMICS ! add new surface and sub-surface litter to soil carbon and nitrogen
  procedure :: add_root_litter => add_root_litter_GIMICS ! add new root litter to soil carbon and nitrogen
  procedure :: add_root_exudates => add_root_exudates_GIMICS ! add root exudates to soil carbon
  procedure :: burn_litter_frac => burn_litter_frac_GIMICS  ! burn a fraction of sfc litter and return amounts of burned carbon and nitrogen
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


! namelist
real :: litt_theta_mod = 0.1 ! Modifies litt theta (unitless)
real :: theta_thres = 0.9    ! Modifies Vmax for Lm fluxes into Mr (unitless)
real :: thetaFmin_dry = 0.0    ! Modifies Vmax for Lm fluxes into Mr (unitless)
real :: thetaFmin_wet = 0.0    ! Modifies Vmax for Lm fluxes into Mr (unitless)

real :: Vmod_Mr_Lm = 10.0    ! Modifies Vmax for Lm fluxes into Mr (unitless)
real :: Vmod_Mr_Ls = 2.0     ! Modifies Vmax for Ls fluxes into Mr (unitless)
real :: Vmod_Mr_Ca = 10.0    ! Modifies Vmax for Ca fluxes into Mr (unitless)
real :: Vmod_Mr_DOC = 10.0    ! Modifies Vmax for Ca fluxes into Mr (unitless)

real :: Vmod_Mk_Lm = 3.0     ! Modifies Vmax for Lm fluxes into Mk (unitless)
real :: Vmod_Mk_Ls = 3.0     ! Modifies Vmax for Ls fluxes into Mk (unitless)
real :: Vmod_Mk_Ca = 2.0     ! Modifies Vmax for Ca fluxes into Mk (unitless)
real :: Vmod_Mk_DOC = 2.0     ! Modifies Vmax for Ca fluxes into Mk (unitless)

real :: Vslope     = 0.063   ! Regression coefficient (ln(mgC/mgM/hr)/Celsius) (Eq 1 in Wieder et al., 2015)
real :: Vint       = 5.47    ! Regression intercept (ln(mgC/mgM/hr)) (Eq 1 in Wieder et al., 2015)
real :: aV         = 8e-6    ! Tuning coefficient (unitless) (Eq 1 in Wieder et al., 2015)

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

! ml: need to fix this fraction in consistent with inputs from vegetation module
!     this fraction should be different for leaf litter, coarse wood litter, rhizosphere, and bulk soil respectively
!     This fraction is used when calculating microbe turnover rates and fractions of dead microbe inputs to Cc pools
real :: fI_Lm = 0.38464225   ! Partitioning of litter inputs to Lm (unitless)

real :: eLm_Mr = 0.55       ! Microbial growth efficiency for fluxes from Lm to Mr (mg/mg)
real :: eLs_Mr = 0.25       ! Microbial growth efficiency for fluxes from Ls to Mr (mg/mg)
real :: eCa_Mr = 0.55       ! Microbial growth efficiency for fluxes from Ca to Mr (mg/mg)
real :: eLm_Mk = 0.75       ! Microbial growth efficiency for fluxes from Lm to Mk (mg/mg)
real :: eLs_Mk = 0.35       ! Microbial growth efficiency for fluxes from Ls to Mk (mg/mg)
real :: eCa_Mk = 0.75       ! Microbial growth efficiency for fluxes from Ca to Mk (mg/mg)
real :: e_slope = 0.0       ! Microbial growth efficiency for fluxes from Ca to Mk (mg/mg)

real :: Kmod_oxid_Mr = 4.0  ! Further modifies Km for oxidation of Cc
real :: Kmod_oxid_Mk = 4.0  ! Further modifies Km for oxidation of Cc

real :: w_Lm = 0.1          ! Fluxes from Lm decomposition to DOC (mg/mg)
real :: w_Ls = 0.1          ! Fluxes from Ls decomposition to DOC (mg/mg)
real :: w_Ca = 0.1          ! Fluxes from Ca decomposition to DOC (mg/mg)

real :: fMrTau_DOC = 0.5    ! Partition between Ca and DOC for Mr turnover Fluxes excluding those to Cc and Cp
real :: fMkTau_DOC = 0.5    ! Partition between Ca and DOC for Mk turnover Fluxes excluding those to Cc and Cp

real :: fMrTau_Cp_a1  = 0.3      ! modifications from Zhang et al. (2019), from 0.3 to 0.13
real :: fMrTau_Cp_a2  = 0.2      ! modifications from Zhang et al. (2019), from 0.2 to 0.02
real :: fMrTau_Cc_a3  = -3.0     ! modifications from Zhang et al. (2019), from -3 to -2.61
real :: fMrTau_Cc_a4  = 0.1      ! modifications from Zhang et al. (2019), from 0.1 to 1.06
real :: fMrTau_Cc_a5  = 0.3     ! modifications from Zhang et al. (2019), from 0.3 to 8.93
real :: fMrTau_Cc_a3_litt  = -3.0     ! modifications from Zhang et al. (2019), from -3 to -2.61
real :: fMrTau_Cc_a4_litt  = 0.1      ! modifications from Zhang et al. (2019), from 0.1 to 1.06
real :: fMrTau_Cc_a5_litt  = 0.3     ! modifications from Zhang et al. (2019), from 0.3 to 8.93

real :: Desorb_kd  = 1.0      ! modifications from Zhang et al. (2019)
real :: Desorb_kdp = 0.0
real :: Desorb_clay = -1.5

logical :: highT_limit = .FALSE.    ! if FALSE, no limitation on decomposition for high temperature
logical :: lowT_limit  = .FALSE.    ! if FALSE, no limitation on decomposition for low temperature
logical :: DOC_cycling = .FALSE.    ! if FALSE, no doc cycling
logical :: Km_theta    = .FALSE.    ! if FALSE, theta function (soil water dependence) does not affect Km

real :: min_anaerobic_resp_factor = 0.05
real :: theta_func_orchidee_min = 0.25
real :: theta_func_orchidee_max = 0.8
real :: min_dry_resp_factor = 0.05
real :: gas_diffusion_exp = 2.5 ! Exponent for gas diffusion power law dependence on theta
                            ! See Meslin et al 2010, SSAJ
real :: substrate_diffusion_exp = 3.0  ! Exponent for theta dependence at low theta.
                            ! See Davison et al DAMM model paper


real :: tau_calib = 1.0   ! Microbial turnover rate calibration factor
real :: tau_beta = 1.66   ! Microbial turnover rate increases with growing microbial biomass density
                          ! Density-dependence exponent

real :: cw_r_cw = 30.0   !< coarse wood radius [cm]
real :: cw_z_cw = 5.0    !< coarse wood thickness accessible by microbes [cm]
real :: lf_f_cw = 1.0    !< leaf litter fraction accessible by microbes, set to 1.0 for coarse wood

real :: cw_r_lf = 1.0    !< coarse wood radius [cm], set to 1.0 for leaf
real :: cw_z_lf = 1.0    !< coarse wood thickness accessible by microbes [cm] , set to 1.0 leaf
real :: lf_f_lf = 0.1    !< leaf litter fraction accessible by microbes

real :: r_rhiz = 0.001      ! Radius of rhizosphere around fine root [m]
real :: litt_density = 22.0 ! C density of surface litter layer [kg/m3]
                            ! 22.0 roughly from Gaudinsky et al 2000, like in CORPSE
real :: min_litt_dz = 0.001 ! minimum litter thickness, [m]
real :: const_litt_dz = -9999.0 ! constant litter thickness, [m]
                            ! if set to value above zero, litter thickness is not updated
                            ! based on density and carbon mass; instead, the model uses
                            ! a constant value specified in this parameter

! Coefficients of turbation, per layer, at the bottom of each layer, [m2/yr]
! slm: eventually these coefficients should be calculated from environmental conditions,
!      e.g. ice content for bio/crio-turbation, or activity of bugs in the soil at
!      different layers, or water table depth, etc.
integer :: ii
! ml: these coefficients are for bioturbation between litter and the first soil layer, and also between soil layers
real :: K_turb(MAX_SOIL_LEV) = (/ &
    1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, &
    0.8e-4, 0.6e-4, 0.4e-4, 0.2e-4, 0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    &
    (0.0,ii=1,80) /)
! ml:  I assume bioturbation coefficient between litter and the first soil layer the same as the one for the 1st and 2nd soil layers
!      We could test sensitivity
real :: K_sfc_turb = 1.0e-4 ! coefficient of exchange between surface litter and soil, [m2/yr]
! ml: for diffusion of DOC, the coefficient should be much larger like 38.8e-4
real :: K_diff(MAX_SOIL_LEV) = (/ &
    38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, &
    38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4, 38.8e-4,    &
    (0.0,ii=1,80) /)
logical :: do_microbe_turb = .TRUE. ! if true, microbes are transported by turbation

real :: init_Mr = 1e-15 ! initial (cold-start) value of microbesR, [kg/m3]
real :: init_Mk = 1e-15 ! initial (cold-start) value of microbesR, [kg/m3]
real :: init_litt_dz = 1e-4 ! initial (cold-start) surface litter thickness, [m]

logical, protected :: save_equilibration_data = .FALSE. !< if TRUE, information for
                         !! soil BGC equilibration acceleration is saved to disk

character(32) :: theta_func_soil = 'CORPSE' ! or 'Yan2018', 'NONE'
character(32) :: theta_func_litt = 'CORPSE' ! or 'ORCHIDEE', 'Yan2018', 'NONE'

namelist /soil_BGC_GIMICS_nml/ &
    litt_theta_mod, theta_thres, thetaFmin_dry, thetaFmin_wet, Vmod_Mr_Lm, Vmod_Mr_Ls, Vmod_Mr_Ca, Vmod_Mr_DOC, Vmod_Mk_Lm, Vmod_Mk_Ls, Vmod_Mk_Ca, Vmod_Mk_DOC, Vslope, Vint, aV, &
    Kmod_Mr_Lm, Kmod_Mr_Ls, Kmod_Mr_Ca, Kmod_Mk_Lm, Kmod_Mk_Ls, Kmod_Mk_Ca, Kslope_Lm, Kslope_Ls, Kslope_Ca, Kint, aK, &
    fI_Lm, eLm_Mr, eLs_Mr, eCa_Mr, eLm_Mk, eLs_Mk, eCa_Mk, e_slope, Kmod_oxid_Mr, Kmod_oxid_Mk, &
    w_Lm, w_Ls, w_Ca, &
    theta_func_litt, theta_func_soil, &
    highT_limit,  lowT_limit, DOC_cycling, Km_theta, &
    min_anaerobic_resp_factor, min_dry_resp_factor, gas_diffusion_exp, substrate_diffusion_exp, theta_func_orchidee_min, theta_func_orchidee_max, &
    tau_calib, tau_beta, cw_r_cw, cw_z_cw, lf_f_cw, cw_r_lf, cw_z_lf, lf_f_lf, &
    fMrTau_DOC, fMkTau_DOC, fMrTau_Cp_a1, fMrTau_Cp_a2, fMrTau_Cc_a3, fMrTau_Cc_a4, fMrTau_Cc_a5, fMrTau_Cc_a3_litt, fMrTau_Cc_a4_litt, fMrTau_Cc_a5_litt, Desorb_kd, Desorb_kdp, Desorb_clay, &
! -----
    init_Mr, init_Mk, init_litt_dz, r_rhiz, litt_density, min_litt_dz, const_litt_dz, &
    K_turb, K_sfc_turb, K_diff, do_microbe_turb, &
    save_equilibration_data

! diag field IDs
integer :: id_total_soil_C
integer :: id_fRhiz, &
   id_sturb_metabolicC, id_sturb_structuralC, id_sturb_chemResistantC, id_sturb_availableC, id_sturb_DOC,&
   id_sturb_microbesR, id_sturb_microbesK, &
   id_negative_litter_C(N_C_TYPES), id_tot_negative_litter_C
! diag fields for rhizosphere, bulk soil, and total
integer, dimension(3) :: id_soilC, id_metabolicC, id_structuralC, id_protectedC, &
   id_chemResistantC, id_availableC, id_microbesR, id_microbesK, id_DOC, id_DecompMrLm, &
   id_DecompMrLs, id_DecompMrCa, id_DecompMkLm, id_DecompMkLs, id_DecompMkCa, id_DecompMrDOC, id_DecompMkDOC, &
   id_OxidMrCc, id_OxidMkCc, id_MrTau, id_MkTau, id_Resp, id_Desorb, id_thetaF, &
   ! input rates
   id_InputStrC, id_InputMtbC, id_InputExdC
! diag fields for column-integrated soil carbon pools:
integer :: id_clmn_metabolicC, id_clmn_structuralC, id_clmn_protectedC, &
   id_clmn_chemResistantC, id_clmn_availableC, id_clmn_microbesR, id_clmn_microbesK, &
   id_clmn_DOC, id_clmn_InputStrC, id_clmn_InputMtbC, id_clmn_InputExdC, &
   id_clmn_Resp, id_clmn_Desorb, id_clmn_Decomp, id_clmn_Oxid, id_clmn_Turnover

integer, dimension(N_LITTER_POOLS) :: id_litt_total_C, id_litt_dz, id_litt_thetaF, &
   id_litt_metabolicC, id_litt_structuralC, id_litt_chemResistantC, id_litt_availableC, &
   id_litt_microbesR, id_litt_microbesK, id_litt_DOC, id_litt_allC, &
   id_litt_DecompMrLm, id_litt_DecompMrLs, id_litt_DecompMrCa, id_litt_DecompMrDOC, &
   id_litt_DecompMkLm, id_litt_DecompMkLs, id_litt_DecompMkCa, id_litt_DecompMkDOC, &
   id_litt_OxidMrCc, id_litt_OxidMkCc, id_litt_MrTau, id_litt_MkTau, id_litt_Resp, &
   ! aggregated terms
   id_litt_Decomp, id_litt_Oxid, id_litt_Turnover, &
   ! input rates for surface litter pools
   id_litt_InputStrC, id_litt_InputMtbC, &
   ! turbation tendencies in surface litter pools
   id_lturb_metabolicC, id_lturb_structuralC, id_lturb_chemResistantC, id_lturb_availableC, id_lturb_DOC,&
   id_lturb_microbesR, id_lturb_microbesK

integer :: id_surf_DOC_loss, id_total_DOC_div_loss, id_sadvec_DOC, id_ladvec_DOC(N_LITTER_POOLS)

! CMIP/CMOR diag fields
integer :: id_rh, id_cSoil, id_cSoilLevels, id_cLitter, id_cLitterCwd, id_cLitterLeaf, &
   id_cSoilAbove1m, id_theta, id_theta_ice

! variables for CMOR/CMIP diagnostic calculations
real, allocatable :: mrs1m_weight(:) ! weights for mrs1m averaging

integer, parameter :: &
     THETA_F_NONE     = 0, &
     THETA_F_ORCHIDEE = 1, &
     THETA_F_YAN2018  = 2, &
     THETA_F_CORPSE   = 3, &
     THETA_F_A        = 4

integer :: theta_func_litt_option = -1 ! integer option corresponding to theta_func_litt namelist parameter
integer :: theta_func_soil_option = -1 ! integer option corresponding to theta_func_soil namelist parameter

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
!> @brief read namelist and set up few constants
subroutine read_soil_BGC_GIMICS_namelist()
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  character(256) :: msg   ! error message

  call log_version(version, module_name, __FILE__)
  read (input_nml_file, nml=soil_BGC_GIMICS_nml, iostat=io, iomsg=msg)
  ierr = check_nml_error(io, 'soil_BGC_GIMICS_nml :: '//msg)
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=soil_BGC_GIMICS_nml)
  endif

  ! parse options of BGC dynamics dependence on moisture
  select case(trim(lowercase(theta_func_litt)))
  case ('orchidee')
     theta_func_litt_option = THETA_F_ORCHIDEE
  case ('yan2018')
     theta_func_litt_option = THETA_F_YAN2018
  case ('corpse')
     theta_func_litt_option = THETA_F_CORPSE
  case ('none')
     theta_func_litt_option = THETA_F_NONE
  case default
     call error_mesg('read_soil_BGC_GIMICS_namelist',&
         'value "'//trim(theta_func_litt)//'" of theta_func_litt is incorrect : use "ORCHIDEE", "Yan2018", "CORPSE", or "none"', FATAL)
  end select

  select case(trim(lowercase(theta_func_soil)))
! slm: should we have ORCHIDEE option for soil too?
!   case ('orchidee')
!      theta_func_soil_option = THETA_F_ORCHIDEE
  case ('yan2018')
     theta_func_soil_option = THETA_F_YAN2018
  case ('corpse')
     theta_func_soil_option = THETA_F_CORPSE
  case ('none')
     theta_func_soil_option = THETA_F_NONE
  case ('afunc')
     theta_func_soil_option = THETA_F_A
  case default
     call error_mesg('read_soil_BGC_GIMICS_namelist',&
         'value "'//trim(theta_func_soil)//'" of theta_func_soil is incorrect : use "Yan2018", "CORPSE", "afunc", or "none"', FATAL)
  end select

  ! store time step in different units in module variables, for convenience
  delta_time = time_type_to_real(lnd%dt_fast) ! [s]
  dt_fast_yr = delta_time/seconds_per_year    ! [yr]
  dt_fast_hr = delta_time/3600.0              ! [hr]
end subroutine

! ============================================================================
!> @brief Register diagnostic fields
subroutine soil_BGC_diag_init_GIMICS(id_ug, id_zfull)
  integer,intent(in)  :: id_ug    !< Unstructured axis id
  integer,intent(in)  :: id_zfull !< Vertical (depth) axis id

  character(*), parameter :: diag_mod_name = 'soil_BGC_GIMICS'
  integer :: axes(2) ! IDs of diagnostic axes
  integer :: k

  ! define array of axis IDs
  axes = [ id_ug,id_zfull ]

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  ! diag field registration goes here
  id_total_soil_C = register_tiled_diag_field ( diag_mod_name, 'tot_soil_C', axes(1:1),  &
       lnd%time, 'total carbon, including soil and litter pools', 'kg C/m2', missing_value=-100.0 )

  id_soilC = register_3_diag_fields ( diag_mod_name, 'soilC', axes(:),  &
       lnd%time, 'Volumetric density of all soil carbon', 'kg C/m3', missing_value=-100.0 )

  id_metabolicC = register_3_diag_fields ( diag_mod_name, 'metabolicC', axes(:),  &
       lnd%time, 'Volumetric density of metabolic C', 'kg C/m3', missing_value=-100.0 )
  id_structuralC = register_3_diag_fields ( diag_mod_name, 'structuralC', axes(:),  &
       lnd%time, 'Volumetric density of structural C', 'kg C/m3', missing_value=-100.0 )
  id_protectedC = register_3_diag_fields ( diag_mod_name, 'protectedC', axes(:),  &
       lnd%time, 'Volumetric density of protected C', 'kg C/m3', missing_value=-100.0 )
  id_chemResistantC = register_3_diag_fields ( diag_mod_name, 'chemResistantC', axes(:),  &
       lnd%time, 'Volumetric density of chemically resistant C', 'kg C/m3', missing_value=-100.0 )
  id_availableC = register_3_diag_fields ( diag_mod_name, 'availableC', axes(:),  &
       lnd%time, 'Volumetric density of available C', 'kg C/m3', missing_value=-100.0 )
  id_microbesR = register_3_diag_fields ( diag_mod_name, 'microbesR', axes(:),  &
       lnd%time, 'Volumetric density of copiotrophic (R) microbes', 'kg C/m3', missing_value=-100.0 )
  id_microbesK = register_3_diag_fields ( diag_mod_name, 'microbesK', axes(:),  &
       lnd%time, 'Volumetric density of oligotrophic (K) microbes', 'kg C/m3', missing_value=-100.0 )
  id_DOC       = register_3_diag_fields ( diag_mod_name, 'DOC', axes(:),  &
       lnd%time, 'Volumetric density of DOC', 'kg C/m3', missing_value=-100.0 )

  id_clmn_metabolicC = register_tiled_diag_field ( diag_mod_name, 'clmn_metabolicC', axes(1:1),  &
       lnd%time, 'Column-integrated metabolic C in soil', 'kg C/m2', missing_value=-100.0 )
  id_clmn_structuralC = register_tiled_diag_field ( diag_mod_name, 'clmn_structuralC', axes(1:1),  &
       lnd%time, 'Column-integrated structural C in soil', 'kg C/m2', missing_value=-100.0 )
  id_clmn_protectedC = register_tiled_diag_field ( diag_mod_name, 'clmn_protectedC', axes(1:1),  &
       lnd%time, 'Column-integrated protected C in soil', 'kg C/m2', missing_value=-100.0 )
  id_clmn_chemResistantC = register_tiled_diag_field ( diag_mod_name, 'clmn_chemResistantC', axes(1:1),  &
       lnd%time, 'Column-integrated chemically resistant C in soil', 'kg C/m2', missing_value=-100.0 )
  id_clmn_availableC = register_tiled_diag_field ( diag_mod_name, 'clmn_availableC', axes(1:1),  &
       lnd%time, 'Column-integrated available C in soil', 'kg C/m2', missing_value=-100.0 )
  id_clmn_microbesR = register_tiled_diag_field ( diag_mod_name, 'clmn_microbesR', axes(1:1),  &
       lnd%time, 'Column-integrated copiotrophic (R) microbes in soil', 'kg C/m2', missing_value=-100.0 )
  id_clmn_microbesK = register_tiled_diag_field ( diag_mod_name, 'clmn_microbesk', axes(1:1),  &
       lnd%time, 'Column-integrated oligotrophic (K) microbes in soil', 'kg C/m2', missing_value=-100.0 )
  id_clmn_DOC = register_tiled_diag_field ( diag_mod_name, 'clmn_DOC', axes(1:1),  &
       lnd%time, 'Column-integrated DOC in soil', 'kg C/m2', missing_value=-100.0 )

  id_DecompMrLm = register_3_diag_fields ( diag_mod_name, 'DecompMrLm', axes(:),  &
       lnd%time, 'Rate of metabolic C decomposition by R microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_DecompMrLs = register_3_diag_fields ( diag_mod_name, 'DecompMrLs', axes(:),  &
       lnd%time, 'Rate of structural C decomposition by R microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_DecompMrCa = register_3_diag_fields ( diag_mod_name, 'DecompMrCa', axes(:),  &
       lnd%time, 'Rate of available C decomposition by R microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_DecompMrDOC = register_3_diag_fields ( diag_mod_name, 'DecompMrDOC', axes(:),  &
       lnd%time, 'Rate of DOC decomposition by R microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_OxidMrCc = register_3_diag_fields ( diag_mod_name, 'OxidMrCc', axes(:),  &
       lnd%time, 'Rate of chemically resistant C oxidation by R microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_MrTau = register_3_diag_fields ( diag_mod_name, 'MrTau', axes(:),  &
       lnd%time, 'Rate of R microbes overturning', 'kg C/m3/h', missing_value=-100.0 )

  id_DecompMkLm = register_3_diag_fields ( diag_mod_name, 'DecompMkLm', axes(:),  &
       lnd%time, 'Rate of metabolic C decomposition by K microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_DecompMkLs = register_3_diag_fields ( diag_mod_name, 'DecompMkLs', axes(:),  &
       lnd%time, 'Rate of structural C decomposition by K microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_DecompMkCa = register_3_diag_fields ( diag_mod_name, 'DecompMkCa', axes(:),  &
       lnd%time, 'Rate of available C decomposition by K microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_DecompMkDOC = register_3_diag_fields ( diag_mod_name, 'DecompMkDOC', axes(:),  &
       lnd%time, 'Rate of DOC decomposition by K microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_OxidMkCc = register_3_diag_fields ( diag_mod_name, 'OxidMkCc', axes(:),  &
       lnd%time, 'Rate of chemically resistant C oxidation by K microbes', 'kg C/m3/h', missing_value=-100.0 )
  id_MkTau = register_3_diag_fields ( diag_mod_name, 'MkTau', axes(:),  &
       lnd%time, 'Rate of K microbes overturning', 'kg C/m3/h', missing_value=-100.0 )

  id_InputStrC = register_3_diag_fields ( diag_mod_name, 'InputStrC', axes(:),  &
       lnd%time, 'Rate of input to structural C', 'kg C/m3/h', missing_value=-100.0 )
  id_InputMtbC = register_3_diag_fields ( diag_mod_name, 'InputMtbC', axes(:),  &
       lnd%time, 'Rate of input to metabolic C', 'kg C/m3/h', missing_value=-100.0 )
  id_InputExdC = register_3_diag_fields ( diag_mod_name, 'InputExdC', axes(:),  &
       lnd%time, 'Rate of C exudate input', 'kg C/m3/h', missing_value=-100.0 )

  id_clmn_InputStrC = register_tiled_diag_field ( diag_mod_name, 'clmn_InputStrC', axes(1:1),  &
       lnd%time, 'Column-integrated rate of input of structural C to soil', 'kg C/m2/h', missing_value=-100.0 )
  id_clmn_InputMtbC = register_tiled_diag_field ( diag_mod_name, 'clmn_InputMtbC', axes(1:1),  &
       lnd%time, 'Column-integrated rate of input of metabolic C to soil', 'kg C/m2/h', missing_value=-100.0 )
  id_clmn_InputExdC = register_tiled_diag_field ( diag_mod_name, 'clmn_InputExdC', axes(1:1),  &
       lnd%time, 'Column-integrated rate of exudate C input to soil', 'kg C/m2/h', missing_value=-100.0 )

  id_Resp = register_3_diag_fields ( diag_mod_name, 'Resp', axes(:),  &
       lnd%time, 'Rate of respiration', 'kg C/m3/h', missing_value=-100.0 )
  id_Desorb = register_3_diag_fields ( diag_mod_name, 'Desorb', axes(:),  &
       lnd%time, 'Rate of desorption', 'kg C/m3/h', missing_value=-100.0 )

  id_clmn_Resp = register_tiled_diag_field ( diag_mod_name, 'clmn_Resp', axes(1:1),  &
       lnd%time, 'Column-integrated rate of respiration in soil', 'kg C/m2/h', missing_value=-100.0 )
  id_clmn_Desorb = register_tiled_diag_field ( diag_mod_name, 'clmn_Desorb', axes(1:1),  &
       lnd%time, 'Column-integrated rate of desorption in soil', 'kg C/m2/h', missing_value=-100.0 )
  id_clmn_Decomp = register_tiled_diag_field ( diag_mod_name, 'clmn_Decomp', axes(1:1),  &
       lnd%time, 'Column-integrated rate of decomposition in soil', 'kg C/m2/h', missing_value=-100.0 )
  id_clmn_Oxid = register_tiled_diag_field ( diag_mod_name, 'clmn_Oxid', axes(1:1),  &
       lnd%time, 'Column-integrated rate of oxidation in soil', 'kg C/m2/h', missing_value=-100.0 )
  id_clmn_Turnover = register_tiled_diag_field ( diag_mod_name, 'clmn_Turnover', axes(1:1),  &
       lnd%time, 'Column-integrated rate of turnover in soil', 'kg C/m2', missing_value=-100.0 )

  id_fRhiz = register_tiled_diag_field ( diag_mod_name, 'fRhiz', axes(:),  &
       lnd%time, 'Volumetric fraction of rhizosphere', 'm3/m3', missing_value=-100.0 )
  id_thetaF = register_3_diag_fields ( diag_mod_name, 'ThetaFunc', axes(:),  &
       lnd%time, 'Soil moisture related factor for decomposition', '-', missing_value=-100.0 )

  id_theta = register_tiled_diag_field ( diag_mod_name, 'Theta', axes(:),  &
       lnd%time, 'Water-filled porosity (fraction of pores filled with water)', '-', missing_value=-100.0 )
  id_theta_ice = register_tiled_diag_field ( diag_mod_name, 'Theta_ice', axes(:),  &
       lnd%time, 'Water-filled porosity (fraction of pores filled with water)', '-', missing_value=-100.0 )

  ! litter fields
  id_litt_total_C(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_C', &
       axes(1:1),  lnd%time, '<ltype> litter total carbon', 'kg C/m2', missing_value=-100.0 )
  id_litt_dz(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_dz', &
       axes(1:1),  lnd%time, '<ltype> litter thickness', 'm', missing_value=-100.0 )

  id_litt_metabolicC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_metabolicC', axes(1:1),  &
       lnd%time, 'Volumetric density of metabolic C in <ltype> litter', 'kg C/m3', missing_value=-100.0 )
  id_litt_structuralC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_structuralC', axes(1:1),  &
       lnd%time, 'Volumetric density of structural C in <ltype> litter', 'kg C/m3', missing_value=-100.0 )
  id_litt_chemResistantC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_chemResistantC', axes(1:1),  &
       lnd%time, 'Volumetric density of chemically resistant C  in <ltype> litter', 'kg C/m3', missing_value=-100.0 )
  id_litt_availableC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_availableC', axes(1:1),  &
       lnd%time, 'Volumetric density of available C in <ltype> litter', 'kg C/m3', missing_value=-100.0 )
  id_litt_microbesR(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_microbesR', axes(1:1),  &
       lnd%time, 'Volumetric density of copiotrophic (R) microbes in <ltype> litter', 'kg C/m3', missing_value=-100.0 )
  id_litt_microbesK(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_microbesK', axes(1:1),  &
       lnd%time, 'Volumetric density of oligotrophic (K) microbes in <ltype> litter', 'kg C/m3', missing_value=-100.0 )
  id_litt_DOC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DOC', axes(1:1),  &
       lnd%time, 'Volumetric density of DOC in <ltype> litter', 'kg C/m3', missing_value=-100.0 )
  id_litt_allC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_allC', axes(1:1),  &
       lnd%time, 'Volumetric density of all C in <ltype> litter', 'kg C/m3', missing_value=-100.0 )

  id_litt_DecompMrLm(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DecompMrLm', axes(1:1),  &
       lnd%time, 'Rate of metabolic C decomposition by R microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_DecompMrLs(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DecompMrLs', axes(1:1),  &
       lnd%time, 'Rate of structural C decomposition by R microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_DecompMrCa(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DecompMrCa', axes(1:1),  &
       lnd%time, 'Rate of available C decomposition by R microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_DecompMrDOC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DecompMrDOC', axes(1:1),  &
       lnd%time, 'Rate of DOC decomposition by R microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_OxidMrCc(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_OxidMrCc', axes(1:1),  &
       lnd%time, 'Rate of chemically resistant C oxidation by R microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_MrTau(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_MrTau', axes(1:1),  &
       lnd%time, 'Rate of R microbes overturning in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )

  id_litt_DecompMkLm(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DecompMkLm', axes(1:1),  &
       lnd%time, 'Rate of metabolic C decomposition by K microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_DecompMkLs(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DecompMkLs', axes(1:1),  &
       lnd%time, 'Rate of structural C decomposition by K microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_DecompMkCa(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DecompMkCa', axes(1:1),  &
       lnd%time, 'Rate of available C decomposition by K microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_DecompMkDOC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DecompMkDOC', axes(1:1),  &
       lnd%time, 'Rate of DOC decomposition by K microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_OxidMkCc(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_OxidMkCc', axes(1:1),  &
       lnd%time, 'Rate of chemically resistant C oxidation by K microbes in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_MkTau(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_MkTau', axes(1:1),  &
       lnd%time, 'Rate of K microbes overturning in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )

  id_litt_Decomp(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_Decomp', axes(1:1),  &
       lnd%time, 'Rate of decomposition in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_Oxid(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_Oxid', axes(1:1),  &
       lnd%time, 'Rate of oxidation in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_Turnover(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_Turnover', axes(1:1),  &
       lnd%time, 'Rate of turnover in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )

  id_litt_InputStrC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_InputStrC', axes(1:1),  &
       lnd%time, 'Rate of input to structural C in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_InputMtbC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_InputMtbC', axes(1:1),  &
       lnd%time, 'Rate of input to metabolic C in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )

  id_litt_Resp(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_Resp', axes(1:1),  &
       lnd%time, 'Rate of respiration in <ltype> litter', 'kg C/m3/h', missing_value=-100.0 )
  id_litt_thetaF(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_ThetaFunc', axes(1:1),  &
       lnd%time, 'Soil moisture related factor for decomposition of <ltype> litter', '-', missing_value=-100.0 )

  ! tendencies due to turbation in soil
  id_sturb_metabolicC = register_tiled_diag_field( diag_mod_name, 'metabolicC_turb', axes(:), &
       lnd%time, 'Tendency of metabolic C due to turbation', 'kg C/(m3 yr)', missing_value = -1e20)
  id_sturb_structuralC = register_tiled_diag_field( diag_mod_name, 'structuralC_turb', axes(:), &
       lnd%time, 'Tendency of structural C due to turbation', 'kg C/(m3 yr)', missing_value = -1e20)
  id_sturb_chemResistantC = register_tiled_diag_field( diag_mod_name, 'chemResistantC_turb', axes(:), &
       lnd%time, 'Tendency of chemically resistant C due to turbation', 'kg C/(m3 yr)', missing_value = -1e20)
  id_sturb_availableC = register_tiled_diag_field( diag_mod_name, 'availableC_turb', axes(:), &
       lnd%time, 'Tendency of available C due to turbation', 'kg C/(m3 yr)', missing_value = -1e20)
  id_sturb_DOC = register_tiled_diag_field( diag_mod_name, 'DOC_turb', axes(:), &
       lnd%time, 'Tendency of DOC due to turbation', 'kg C/(m3 yr)', missing_value = -1e20)
  id_sturb_microbesR = register_tiled_diag_field( diag_mod_name, 'microbesR_turb', axes(:), &
       lnd%time, 'Tendency of R microbes due to turbation', 'kg C/(m3 yr)', missing_value = -1e20)
  id_sturb_microbesK = register_tiled_diag_field( diag_mod_name, 'microbesK_turb', axes(:), &
       lnd%time, 'Tendency of K microbes due to turbation', 'kg C/(m3 yr)', missing_value = -1e20)

  ! turbation exchange between surface litter and soil
  id_lturb_metabolicC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_metabolicC_turb', axes(1:1), &
       lnd%time, '<ltype> litter tendency of metabolic C due to turbation', 'kg C/(m3 yr)', missing_value=-100.0 )
  id_lturb_structuralC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_structuralC_turb', axes(1:1), &
       lnd%time, '<ltype> litter tendency of structural C due to turbation', 'kg C/(m3 yr)', missing_value=-100.0 )
  id_lturb_chemResistantC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_chemResistantC_turb', axes(1:1), &
       lnd%time, '<ltype> litter tendency of chemically resistant C due to turbation', 'kg C/(m3 yr)', missing_value=-100.0 )
  id_lturb_availableC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_availableC_turb', axes(1:1), &
       lnd%time, '<ltype> litter tendency of available C due to turbation', 'kg C/(m3 yr)', missing_value=-100.0 )
  id_lturb_DOC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DOC_turb', axes(1:1), &
       lnd%time, '<ltype> litter tendency of DOC due to turbation', 'kg C/(m3 yr)', missing_value=-100.0 )
  id_lturb_microbesR(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_microbesR_turb', axes(1:1), &
       lnd%time, '<ltype> litter tendency of R microbes due to turbation', 'kg C/(m3 yr)', missing_value=-100.0 )
  id_lturb_microbesK(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_microbesK_turb', axes(1:1), &
       lnd%time, '<ltype> litter tendency of K microbes due to turbation', 'kg C/(m3 yr)', missing_value=-100.0 )

  do k = 1, N_C_TYPES
     id_negative_litter_C(k) = register_tiled_diag_field(diag_mod_name, trim(c_diagname(k))//'_negative_litter_C', &
             axes(1:1), lnd%time, 'Cumulative negative '//trim(c_longname(k))//' carbon litter input', &
             'kg C/m2', missing_value = +1e20)
  enddo
  id_tot_negative_litter_C = register_tiled_diag_field(diag_mod_name, 'total_negative_litter_C', axes(1:1), &
       lnd%time, 'Total cumulative negative carbon litter input', 'kg C/m2', missing_value = +1e20)

  ! DOC-related fields
  id_total_DOC_div_loss = register_tiled_diag_field ( diag_mod_name, 'tot_DOC_div', axes(1:1), &
       lnd%time, 'total rate of DOC divergence loss', 'kg C/(m2 yr)', missing_value=-100.0)
  id_surf_DOC_loss = register_tiled_diag_field ( diag_mod_name, 'surf_DOC_loss', axes(1:1), &
       lnd%time, 'loss of top layer DOC to surface runoff due to efflux', 'kg C/(m2 yr)', &
       missing_value=-100.0)
  id_sadvec_DOC = register_tiled_diag_field ( diag_mod_name, 'DOC_advec', axes, &
       lnd%time, 'Tendency of DOC due to all advective processes',  'kg C/(m3 yr)', missing_value=-100.0)
  id_ladvec_DOC(:) = register_litter_diag_fields ( diag_mod_name, '<ltype>litt_DOC_advec', axes(1:1), &
       lnd%time, 'Tendency of <ltype> litter DOC due to all advective processes',  'kg C/(m3 yr)', missing_value=-100.0)

  ! CMOR fields
  ! set the default sub-sampling filter for the CMOR fields below
  call set_default_diag_filter('land')

  ! set up weights for vertical averaging
  allocate(mrs1m_weight(num_l))
  do k = 1,num_l
     ! zhalf(k) is the depth of layer k top; 1m - zhalf(k) is therefore the distance
     ! from layer top to 1m, and (1m - zhalf(k))/dz(k) is the fraction of layer that
     ! is within [0, 1m] range
     ! If zhalf(k) > 1m, then the fraction is zero (1m-zhalf(k)<0)
     ! If zhalf(k) + dz(k) < 1m (layer bottom is above 1m), then the fraction is 1
     !   because (1m-zhalf(k)>dz(k))
     mrs1m_weight(k) = min(1.0,max(0.0,(CMOR_1M_DEPTH-zhalf(k))/dz(k)))
  enddo

  id_rh = register_tiled_diag_field ( CMOR_NAME, 'rh', [ id_ug ], &
       lnd%time, 'Heterotrophic Respiration', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_heterotrophic_respiration', &
       fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_rh, CMOR_NAME, 'rhLut', axes(1:1),  &
       lnd%time, 'Soil Heterotrophic Respiration On Land Use Tile', 'kg m-2 s-1', &
       standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_heterotrophic_respiration', &
       fill_missing=.FALSE., missing_value=-100.0)
  id_cSoil = register_tiled_diag_field ( CMOR_NAME, 'cSoil', axes(1:1),  &
       lnd%time, 'Carbon in Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_mass_content_of_carbon', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_cSoil, CMOR_NAME, 'cSoilLut', axes(1:1),  &
       lnd%time, 'Carbon  In Soil Pool On Land Use Tiles', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_mass_content_of_carbon', fill_missing=.FALSE.)
  id_cSoilLevels = register_tiled_diag_field ( CMOR_NAME, 'cSoilLevels', axes(:),  lnd%time, &
       'Carbon mass in each model soil level (summed over all soil carbon pools in that level)', &
       'kg m-2', missing_value=-100.0, standard_name='soil_mass_content_of_carbon', &
       fill_missing=.TRUE.)
  id_cSoilAbove1m = register_tiled_diag_field ( CMOR_NAME, 'cSoilAbove1m', axes(1:1),  &
       lnd%time, 'Carbon mass in soil pool above 1m depth', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_mass_content_of_carbon', fill_missing=.TRUE.)
  id_cLitter = register_tiled_diag_field ( CMOR_NAME, 'cLitter', axes(1:1), &
       lnd%time, 'Carbon Mass in Litter Pool', 'kg m-2', &
       missing_value=-100.0, standard_name='litter_mass_content_of_carbon', &
       fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_cLitter, CMOR_NAME, 'cLitterLut', axes(1:1),  &
       lnd%time, 'carbon in above and belowground litter pools on land use tiles', &
       'kg m-2', missing_value=-100.0, &
       standard_name='litter_mass_content_of_carbon', fill_missing=.FALSE.)
  id_cLitterCwd = register_tiled_diag_field ( CMOR_NAME, 'cLitterCwd', axes(1:1), &
       lnd%time, 'Carbon Mass in Coarse Woody Debris', 'kg m-2', &
       missing_value=-100.0, standard_name='wood_debris_mass_content_of_carbon', &
       fill_missing=.TRUE.)
  id_cLitterLeaf = register_tiled_diag_field ( CMOR_NAME, 'cLitterLeaf', axes(1:1), &
       lnd%time, 'Carbon Mass in Leaf Debris', 'kg m-2', &
       missing_value=-100.0, standard_name='leaf_debris_mass_content_of_carbon', &
       fill_missing=.TRUE.)

end subroutine

! ============================================================================
!> register an array of three diag fields: for rhizosphere, bulk, and average
function register_3_diag_fields(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name) result (id)

  integer :: id(3)

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

  id(1) = register_tiled_diag_field(module_name, &
          'rhiz_'//field_name, &
          axes, init_time, &
          trim(long_name)//' in rhizosphere', &
          units, missing_value, range, op, standard_name)
  id(2) = register_tiled_diag_field(module_name, &
          'bulk_'//field_name, &
          axes, init_time, &
          trim(long_name)//' in bulk soil', &
          units, missing_value, range, op, standard_name)
  id(3) = register_tiled_diag_field(module_name, &
          field_name, &
          axes, init_time, &
          long_name, &
          units, missing_value, range, op, standard_name)
end function register_3_diag_fields

! ============================================================================
!> @brief Send rhizosphere and bulk soil data to diag, separately and an average
subroutine send_3_tile_data(id, rhiz, bulk, fRhiz, diag)
  integer, intent(in) :: id(3)
  real, intent(in) :: rhiz(:)  !< rhizosphere values
  real, intent(in) :: bulk(:)  !< bulk soil values
  real, intent(in) :: fRhiz(:) !< fraction of rhizosphere
  type(diag_buff_type), intent(inout) :: diag !< diagnostic buffer

  if (id(1)>0) call send_tile_data(id(1), rhiz(:), diag)
  if (id(2)>0) call send_tile_data(id(2), bulk(:), diag)
  ! slm: are we creating temp array here? any way to avoid doing it?
  if (id(3)>0) call send_tile_data(id(3), rhiz(:)*fRhiz(:) + bulk(:)*(1-fRhiz(:)), diag)
end subroutine send_3_tile_data

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
  ptr%fRhiz(:) = 0.0 ! slm: is this reasonable for initial rhizosphere fraction?
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
     ! the rest of the BGC pool fields remain at their initial values of zero.
     ! set initial surface litter thickness
     if (const_litt_dz > 0) then
        soilc%litt(k)%dz = const_litt_dz
     else
        soilc%litt(k)%dz = init_litt_dz
     endif

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
  real :: dz2

  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  select type(s1)
  type is (soil_BGC_GIMICS_t)
     ! merge surface litter pools
     ! for each state variable C of the pools, we must conserve mass, so that
     ! after merge
     ! C2'*dz2' = x1*C1*dz1 + x2*C2*dz2
     ! therefore
     ! C2' = x1*dz1/dz2'*C1 + x2*dz2/dz2'*C2
     ! dz2' can be set to dz2, or to any convenient value and then updated based on
     ! resulting BGC mass and density
     do k = 1, N_LITTER_POOLS
        ! set dz2' to the max of two litter thicknesses to avoid dividing by zero
        dz2 = max(s1%litt(k)%dz,s2%litt(k)%dz)
        if (dz2 > 0) then
           call combine_GIMICS_pools(s2%litt(k), x2*s2%litt(k)%dz/dz2, &
                                   s1%litt(k), x1*s1%litt(k)%dz/dz2  )
           s2%litt(k)%dz = dz2
           call update_litter_thickness(s2%litt(k))
        else
           ! Do nothing, since the mass is zero in both of the input pools.
           ! This assumes that if litter thickness is zero then the mass of every
           ! BGC component in litter is also zero.
        endif
     enddo
     ! Merge soil pools, layer by layer
     ! For soil layers thickness is the same in both BGC pools that are merged,
     ! so the averaging weights do not include thickness
     do k = 1,size(s2%rhiz)
        ! rhizosphere pools
        f1 = s1%fRhiz(k)         ; f2 = s2%fRhiz(k)
        if (f1>0.or.f2>0) then

	!   if (.not.(x1*f1+x2*f2 .gt. 0)) then
	!      __DEBUG5__(x1,f1,x2,f2,x1*f1+x2*f2)
	!      call land_error_message('denominator is 0 in rhiz tile merge',FATAL)
	!   endif

        !   y1 = x1*f1/(x1*f1+x2*f2) ; y2 = 1.0 - y1

	   if (x1*f1+x2*f2 > 0) then
              y1 = x1*f1/(x1*f1+x2*f2)
           else
              y1 = 0.0
           endif
           y2 = 1.0 - y1

           call combine_GIMICS_pools(s2%rhiz(k),y2,s1%rhiz(k),y1)
        endif
        ! bulk pools
        f1 = 1.0 - s1%fRhiz(k)   ; f2 = 1.0 - s2%fRhiz(k)
        if (f1>0.or.f2>0) then

	 !  if (.not.(x1*f1+x2*f2 .gt. 0)) then
	 !     call land_error_message('denominator is 0 in bulk tile merge',FATAL)
	 !  endif

         !  y1 = x1*f1/(x1*f1+x2*f2) ; y2 = 1.0 - y1

	   if (x1*f1+x2*f2 > 0) then
              y1 = x1*f1/(x1*f1+x2*f2)
           else
              y1 = 0.0
           endif
           y2 = 1.0 - y1

           call combine_GIMICS_pools(s2%bulk(k),y2,s1%bulk(k),y1)
        endif
        ! update the rhizosphere fraction
        s2%fRhiz(k) = x1*s1%fRhiz(k) + x2*s2%fRhiz(k)
     enddo

     s2%neg_litt_C(:)  = s1%neg_litt_C(:)*x1 + s2%neg_litt_C(:)*x2
  class default
     call land_error_message('merge_GIMICS: attempt to merge incompatible soil carbon types', FATAL)
  end select
end subroutine

! ============================================================================
!> Combines two GIMICS pools p1 and p2, with given weights, and puts the
!! result in p2.
!!
!! Currently, for each state variable X of the pool representing carbon concentration,
!! the resulting value is calculated as linear combination
!!  X_2' = w_1*X_1 + w_2*x_2
subroutine combine_GIMICS_pools(p2,w2,p1,w1)
  class(GIMICS_BGC_pool), intent(inout) :: p2
  class(GIMICS_BGC_pool), intent(in)    :: p1
  real,                   intent(in)    :: w2, w1

#define __MERGE__(var) p2%var = w2*p2%var + w1*p1%var
  __MERGE__(metabolicLitterC)
  __MERGE__(structuralLitterC)
  __MERGE__(protectedC)
  __MERGE__(chemResistantC)
  __MERGE__(availableC)
  __MERGE__(microbesR)
  __MERGE__(microbesK)
  __MERGE__(DOC)

  __MERGE__(DecompMrLm)
  __MERGE__(DecompMrLs)
  __MERGE__(DecompMrCa)
  __MERGE__(DecompMrDOC)
  __MERGE__(DecompMkLm)
  __MERGE__(DecompMkLs)
  __MERGE__(DecompMkCa)
  __MERGE__(DecompMkDOC)
  __MERGE__(OxidMrCc)
  __MERGE__(OxidMkCc)
  __MERGE__(Desorb)
  __MERGE__(MrTau)
  __MERGE__(MkTau)
  __MERGE__(Resp)
  __MERGE__(InputStrC)
  __MERGE__(InputMtbC)
  __MERGE__(InputExdC)
#undef __MERGE__
end subroutine

! ============================================================================
!> @brief Update thickness of a surface litter pool
subroutine update_litter_thickness(pool)
  type(GIMICS_BGC_litt), intent(inout) :: pool !< litter BGC pool

  real :: dz_new

  if (const_litt_dz > 0) then
     ! use constant litter thickness
     dz_new = const_litt_dz
  else
     dz_new = max(C_amount(pool)/litt_density, min_litt_dz)
  endif

  if (dz_new > 0) then
     ! change concentrations to keep the total amounts constant
     call scale_pool(pool,pool%dz/dz_new)
     pool%dz = dz_new
  else
     ! zero litter thickness will cause issues, for example in the litter-to-soil
     ! turbation, because it will lead to infinite changes in concentration
     ! even for finite changes in litter amount; most likely in other places too.
     call land_error_message('litter thickness is not positive', FATAL)
  endif
end subroutine update_litter_thickness


! ============================================================================
!> @brief scale all the values in the pool with specified factor
subroutine scale_pool(pool,f)
  class(GIMICS_BGC_litt), intent(inout) :: pool !< pool to update
  real, intent(in) :: f !< scaling factor

#define __SCALE__(var) pool%var = pool%var*f
! concentrations of various carbon pools, [kgC/m3]
  __SCALE__(metabolicLitterC)
  __SCALE__(structuralLitterC)
  __SCALE__(protectedC)
  __SCALE__(chemResistantC)
  __SCALE__(availableC)
  __SCALE__(microbesR)
  __SCALE__(microbesK)
  __SCALE__(DOC)

! tendencies, for diagnostics, [kgC/m3/hr]
  __SCALE__(DecompMrLm)
  __SCALE__(DecompMrLs)
  __SCALE__(DecompMrCa)
  __SCALE__(DecompMkLm)
  __SCALE__(DecompMkLs)
  __SCALE__(DecompMkCa)
  __SCALE__(DecompMrDOC)
  __SCALE__(DecompMkDOC)
  __SCALE__(OxidMrCc)
  __SCALE__(OxidMkCc)
  __SCALE__(Desorb)
  __SCALE__(MrTau)
  __SCALE__(MkTau)
  __SCALE__(Resp)
#undef __SCALE__
end subroutine

! ============================================================================
!> @brief Change the rhizosphere fraction in the soil
subroutine set_fRhiz(soilc,fRhiz)
  class(soil_BGC_GIMICS_t), intent(inout) :: soilc !< soil carbon data structure
  real, intent(in) :: fRhiz(:) !< new fractions of rhizosphere, by layer. Unitless, [0,1]

  integer :: k
  real :: wr,wb ! weights for rhizosphere and bulk
  real :: w     ! sum of the wr and wb, to normalize the weights

  do k = 1, size(soilc%fRhiz)
     if (fRhiz(k) < soilc%fRhiz(k)) then
        ! part of rhizosphere becomes bulk soil
        wr = soilc%fRhiz(k) - fRhiz(k)
        wb = 1.0 - soilc%fRhiz(k)
        w  = wr+wb
        call combine_GIMICS_pools(soilc%bulk(k),wb/w,soilc%rhiz(k),wr/w)
        soilc%fRhiz(k) = fRhiz(k)
     else if (fRhiz(k) > soilc%fRhiz(k)) then
        ! part of bulk soil becomes rhizosphere
        wb = fRhiz(k) - soilc%fRhiz(k)
        wr = soilc%fRhiz(k)
        w  = wr+wb
        call combine_GIMICS_pools(soilc%rhiz(k),wr/w,soilc%bulk(k),wb/w)
        soilc%fRhiz(k) = fRhiz(k)
     else
        ! do nothing, rhizosphere fraction did not change (or one of fRhiz is a NaN)
     endif
  enddo
end subroutine

! ============================================================================
!> @brief Given soil carbon state, return total C, including carbon in soil and
!! surface litter pools
!! @return Total carbon in soil and surface litter, kgC/m2
real function total_C_GIMICS(soilc) result(answer)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc !< soil carbon data structure

  integer :: k

  answer = sum(soilc%neg_litt_C)
  do k = 1, N_LITTER_POOLS
     answer = answer + C_amount(soilc%litt(k))
  enddo

  answer = answer + soilc%total_soil_C()
end function

! ============================================================================
! > @brief Given the fraction of rhizosphere in each layer and two arrays of
!! the concentrations (for rhizosphere and bulk soil), calculate total
!! amount in the entire soil
real function total_amount(fRhiz, rhiz, bulk) result(answer)
  real, intent(in) :: fRhiz(:) !< fraction of rhizosphere
  real, intent(in) :: rhiz(:)  !< concentration in the rhizosphere
  real, intent(in) :: bulk(:)  !< concentration in the bulk soil

  integer :: k

  answer = 0.0
  do k = 1,num_l
     answer = answer + &
           ( rhiz(k) * fRhiz(k)     &
           + bulk(k) * (1-fRhiz(k)) &
           ) * dz(k)
  enddo
end function

! ============================================================================
!> @brief Given soil carbon state, return total soil C
!! @return Total carbon in soil (not including surface litter), kgC/m2
real function total_soil_C_GIMICS(soilc) result(answer)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc !< soil carbon data structure

  integer :: k

  answer = 0.0
  do k = 1,num_l
     answer = answer + &
           ( C_density(soilc%rhiz(k)) * soilc%fRhiz(k)     &
           + C_density(soilc%bulk(k)) * (1-soilc%fRhiz(k)) &
           ) * dz(k)
  enddo
end function

! ============================================================================
!> @brief Given soil carbon state, and a depth, return total soil C in the layer
!! from the surface to the specified depth
!! @return total soil carbon in the depth range [0,arg], kgC/m2
real function total_soil_C_to_depth_GIMICS(soilc, arg) result(answer)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc !< soil carbon data structure
  real, intent(in) :: arg !< depth over which to calculate the total

  integer :: k ! layer counter
  real :: z    ! depth to the top of the current layer
  real :: dz1  ! thickness of the current layer that is within the interval [0,arg]

  answer = 0.0; z = 0.0
  do k = 1, num_l
     if (z.ge.arg) exit ! from loop
     dz1 = max(min(arg-z,dz(k)),0.0)
     answer = answer + &
           ( C_density(soilc%rhiz(k)) * soilc%fRhiz(k)     &
           + C_density(soilc%bulk(k)) * (1-soilc%fRhiz(k)) &
           ) * dz1
     z = z+dz(k)
  enddo
end function

! ============================================================================
!> @brief Given soil carbon state, return total soil N
!! @return Total nitrogen in soil (not including surface litter), kgN/m2
real function total_soil_N_GIMICS(soilc) result(answer)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc !< soil carbon data structure

  answer = 0.0
end function

! ============================================================================
!> @brief Given soil carbon state, return total soil C by layer, kgC/m2
subroutine totC_by_layer_GIMICS(soilc, values)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc !< soil carbon data structure
  real,                     intent(out) :: values(:) ! (num_l)

  integer :: k

  values(:) = 0.0
  do k = 1, min(size(values),num_l)
     values(k) = ( C_density(soilc%rhiz(k)) * soilc%fRhiz(k)     &
                 + C_density(soilc%bulk(k)) * (1-soilc%fRhiz(k)) &
                 ) * dz(k)
  enddo
end subroutine

! ============================================================================
!> @brief Given a BGC pool, calculate total volumetric density of carbon, kgC/m3
!! @return Volumetric density of carbon in the pool, kgC/m3
real function C_density(pool) result(answer)
  class(GIMICS_BGC_pool), intent(in) :: pool
  answer = pool%metabolicLitterC + pool%structuralLitterC &
         + pool%protectedC + pool%chemResistantC &
         + pool%availableC + pool%microbesR + pool%microbesK + pool%DOC
end function

! ============================================================================
!> @brief Given surface litter BGC pool, calculate total amount of carbon, kgC/m2
!! @return Total carbon in the pool, kgC/m2
real function C_amount(pool) result(answer)
  class(GIMICS_BGC_litt), intent(in) :: pool

  answer = C_density(pool) * pool%dz
end function

! ============================================================================
!> @brief Given soil carbon state, return total soil nitrogen
!! @return total soil nitrogen, kgN/m2
real function total_N_GIMICS(soilc) result(answer)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc ! soil carbon data structure
  answer = 0.0
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
     fast_C = pool%metabolicLitterC  * pool%dz ! slm: check the definition of fast/slow pools
     slow_C = pool%structuralLitterC * pool%dz
     dmic_C = (pool%microbesR + pool%microbesK) * pool%dz
  end associate
end subroutine

! ============================================================================
!> @brief Given soil carbon state, returns total litter C per litter pool
subroutine get_littC_GIMICS(soilc, values)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc     !< soil carbon data structure
  real, intent(out)                     :: values(:) !< total C in litter, by litter pool, kgC/m2

  integer :: i
  do i = 1, N_LITTER_POOLS
     values(i) = C_amount(soilc%litt(i))
  enddo
end subroutine

! ============================================================================
! @brief Given soil carbon state, returns amount of dissolved organic carbon (DOC) by
! type and by layer. In GIMICS, there is only one type of DOC, so first element of the
! first dimension is non-zero in the output values.
subroutine get_DOC_GIMICS(soilc, values)
  class(soil_BGC_GIMICS_t), intent(in)  :: soilc !< soil carbon data structure
  real, intent(out) :: values(:,:) !< dissolved organic carbon (N_C_TYPES, num_l) [kg C/m^2]

  integer :: k

  values(:,:)=0.0
  do k=1,num_l
     values(1,k)=( &
        soilc%rhiz(k)%DOC *    soilc%fRhiz(k) + &
        soilc%bulk(k)%DOC * (1-soilc%fRhiz(k))  &
        ) * dz(k)
  end do
end subroutine

! ============================================================================
! Update the state of the soil BGC pools due to the soil microbiology and other
! processes, and accumulate heterotrophic respiration
subroutine dsdt_GIMICS(soilc, soil, vegn, diag, soilt, theta)
  class(soil_BGC_GIMICS_t)  , intent(inout) :: soilc
  type(vegn_tile_type), intent(inout) :: vegn  !< vegetation data structure
  type(soil_tile_type), intent(inout) :: soil  !< soil data structure
  type(diag_buff_type), intent(inout) :: diag  !< diagnostic buffer
  real                , intent(in)    :: soilt !< average soil temperature, deg K, [unused]
  real                , intent(in)    :: theta !< average soil moisture [unused]

  real, dimension(num_l) :: decomp_T, decomp_theta, decomp_porosity, decomp_moist, &
                            decomp_theta_ice, decomp_water_ice_porosity
  real, dimension(num_l) :: rhiz_frac
  real :: clay_frac ! fraction of clay, unitless in interval [0,1]. Should it be by-layer?
  integer :: k

  vegn%rh = 0.0 ! reset rate of heterotrophic respiration

  decomp_theta = soil_theta(soil)
  decomp_porosity = soil_porosity(soil)
  decomp_moist = decomp_theta * decomp_porosity
  decomp_T     = soil%T(1:num_l) - tfreeze
  decomp_theta_ice = soil_ice_porosity(soil)
  decomp_water_ice_porosity = soil_water_ice_porosity(soil)

  !  First surface litter is decomposed
  do k = 1,N_LITTER_POOLS
     call check_GIMICS_pool(soilc%litt(k), trim(l_diagname(k))//'litt before update')
  enddo

!  do k = 1,N_LITTER_POOLS
     call update_GIMICS_pool(soilc%litt(1), decomp_T(1), decomp_theta(1), decomp_porosity(1), decomp_moist(1), fClay=0.0, cw_r=cw_r_lf, cw_z=cw_z_lf, lf_f=lf_f_lf, is_sfc_litter=.TRUE.)
     call update_GIMICS_pool(soilc%litt(2), decomp_T(1), decomp_theta(1), decomp_porosity(1), decomp_moist(1), fClay=0.0, cw_r=cw_r_cw, cw_z=cw_z_cw, lf_f=lf_f_cw, is_sfc_litter=.TRUE.)

     ! accumulate loss of C to atmosphere [kgC/m2/year]
     vegn%rh=vegn%rh + soilc%litt(1)%Resp*soilc%litt(1)%dz*hours_per_year
     vegn%rh=vegn%rh + soilc%litt(2)%Resp*soilc%litt(2)%dz*hours_per_year

!      do i = 1, N_C_TYPES
!         call send_tile_data(id_litter_rsoil_C(k,i), litter_C_loss_rate(i), diag)
!         call send_tile_data(id_litter_rsoil_N(k,i), litter_N_loss_rate(i), diag)
!      enddo
     ! for budget check
!      vegn%fsc_out     = vegn%fsc_out     + litter_C_loss_rate(C_FAST)*dt_fast_yr
!      vegn%ssc_out     = vegn%ssc_out     + litter_C_loss_rate(C_SLOW)*dt_fast_yr
!      vegn%deadmic_out = vegn%deadmic_out + litter_C_loss_rate(C_MIC) *dt_fast_yr

     call update_litter_thickness(soilc%litt(1))
     call update_litter_thickness(soilc%litt(2))
!  enddo



  do k = 1,N_LITTER_POOLS
     call check_GIMICS_pool(soilc%litt(k), trim(l_diagname(k))//'litt after update')
  enddo

  ! Next we have to go through layers and decompose the soil carbon pools
  call rhizosphere_frac(vegn, rhiz_frac)
  call set_fRhiz(soilc,rhiz_frac)
  clay_frac = soil_pClay(soil)/100.0
  do k=1,num_l
     call check_GIMICS_pool(soilc%rhiz(k), 'rhiz('//string(k)//') before update')
     call check_GIMICS_pool(soilc%bulk(k), 'bulk('//string(k)//') before update')

     call update_GIMICS_pool(soilc%rhiz(k), decomp_T(k), decomp_theta(k), decomp_porosity(k), decomp_moist(k), clay_frac, cw_r=1.0, cw_z=1.0, lf_f=1.0, is_sfc_litter=.FALSE.)
     call update_GIMICS_pool(soilc%bulk(k), decomp_T(k), decomp_theta(k), decomp_porosity(k), decomp_moist(k), clay_frac, cw_r=1.0, cw_z=1.0, lf_f=1.0, is_sfc_litter=.FALSE.)
     ! accumulate loss of C to atmosphere [kgC/m2/year]
     vegn%rh = vegn%rh + soilc%rhiz(k)%Resp*dz(k)*hours_per_year*soilc%fRhiz(k)      &
                       + soilc%bulk(k)%Resp*dz(k)*hours_per_year*(1-soilc%fRhiz(k))

     call check_GIMICS_pool(soilc%rhiz(k), 'rhiz('//string(k)//') after update')
     call check_GIMICS_pool(soilc%bulk(k), 'bulk('//string(k)//') after update')
  enddo

  ! calculate tendencies due to turbation
  ! slm: Minjin seems to apply turbation only to the four components of the soil carbon.
  !      For some reason, protectedC and microbes are not included?
  ! ml: microbes should be also moved via turbation, but not protected C.
  call turbation(soilc%litt(:)%metabolicLitterC, soilc%rhiz(:)%metabolicLitterC, soilc%bulk(:)%metabolicLitterC, &
                 soilc%fRhiz, soilc%litt(:)%dz, K_turb, decomp_theta, decomp_porosity, decomp_moist, id_sturb_metabolicC, id_lturb_metabolicC, diag, 'metabolicC')
  call turbation(soilc%litt(:)%structuralLitterC, soilc%rhiz(:)%structuralLitterC, soilc%bulk(:)%structuralLitterC, &
                 soilc%fRhiz, soilc%litt(:)%dz, K_turb, decomp_theta, decomp_porosity, decomp_moist, id_sturb_structuralC, id_lturb_structuralC, diag, 'structuralC')
  call turbation(soilc%litt(:)%chemResistantC, soilc%rhiz(:)%chemResistantC, soilc%bulk(:)%chemResistantC, &
                 soilc%fRhiz, soilc%litt(:)%dz, K_turb, decomp_theta, decomp_porosity, decomp_moist, id_sturb_chemResistantC, id_lturb_chemResistantC, diag, 'chemResistantC')
  call turbation(soilc%litt(:)%availableC, soilc%rhiz(:)%availableC, soilc%bulk(:)%availableC, &
                 soilc%fRhiz, soilc%litt(:)%dz, K_turb, decomp_theta, decomp_porosity, decomp_moist, id_sturb_availableC, id_lturb_availableC, diag, 'availableC')
  call turbation2(soilc%litt(:)%DOC, soilc%rhiz(:)%DOC, soilc%bulk(:)%DOC, &
                 soilc%fRhiz, soilc%litt(:)%dz, K_diff, decomp_theta, decomp_theta_ice, decomp_water_ice_porosity, decomp_porosity, decomp_moist, id_sturb_DOC, id_lturb_DOC, diag, 'DOC')
  if (do_microbe_turb) then
     call turbation(soilc%litt(:)%microbesR, soilc%rhiz(:)%microbesR, soilc%bulk(:)%microbesR, &
                    soilc%fRhiz, soilc%litt(:)%dz, K_turb, decomp_theta, decomp_porosity, decomp_moist, id_sturb_microbesR, id_lturb_microbesR, diag, 'microbesR', &
                    allow_flux_to_sfc_litter=.TRUE.)
     call turbation(soilc%litt(:)%microbesK, soilc%rhiz(:)%microbesK, soilc%bulk(:)%microbesK, &
                    soilc%fRhiz, soilc%litt(:)%dz, K_turb, decomp_theta, decomp_porosity, decomp_moist, id_sturb_microbesk, id_lturb_microbesK, diag, 'microbesK', &
                    allow_flux_to_sfc_litter=.TRUE.)
  ! perhaps it would be useful to have "else" statement here to send zeros to the
  ! diagnostics of microbe tendencies due to turbation
  endif


  do k = 1,N_LITTER_POOLS
     call check_GIMICS_pool(soilc%litt(k), trim(l_diagname(k))//'litt after turbation')
  enddo
  do k=1,num_l
     call check_GIMICS_pool(soilc%rhiz(k), 'rhiz('//string(k)//') after turbation')
     call check_GIMICS_pool(soilc%bulk(k), 'bulk('//string(k)//') after turbation')
  enddo

  ! slm: TODO: calculate horizontal exchange between rhizosphere and soil

  do k = 1,N_LITTER_POOLS
     call update_litter_thickness(soilc%litt(k))
  enddo
  do k = 1,N_LITTER_POOLS
     call check_GIMICS_pool(soilc%litt(k), trim(l_diagname(k))//'litt after thickness update')
  enddo

  call send_tile_data(id_rh, vegn%rh/seconds_per_year, diag)

  do k = 1, N_C_TYPES
     if (id_negative_litter_C(k)>0) call send_tile_data(id_negative_litter_C(k),soilc%neg_litt_C(k),diag)
  enddo
  if (id_tot_negative_litter_C>0) call send_tile_data(id_tot_negative_litter_C,sum(soilc%neg_litt_C),diag)

  call send_tile_data(id_theta, decomp_theta, diag)
  call send_tile_data(id_theta_ice, decomp_theta_ice, diag)

end subroutine dsdt_GIMICS

! ============================================================================
!> checks that values of GIMICS pool state variables are within the reasonable range
subroutine check_GIMICS_pool(pool, tag)
  class(GIMICS_BGC_pool), intent(inout) :: pool
  character(*), intent(in) :: tag

  call check_var_range(pool%metabolicLitterC , 0.0, HUGE(1.0), tag, 'metabolicLitterC',  FATAL)
  call check_var_range(pool%structuralLitterC, 0.0, HUGE(1.0), tag, 'structuralLitterC', FATAL)
  call check_var_range(pool%protectedC       , 0.0, HUGE(1.0), tag, 'protectedC',        FATAL)
  call check_var_range(pool%chemResistantC   , 0.0, HUGE(1.0), tag, 'chemResistantC',    FATAL)
  call check_var_range(pool%availableC       , 0.0, HUGE(1.0), tag, 'availableC',        FATAL)
  call check_var_range(pool%microbesR        , 0.0, HUGE(1.0), tag, 'microbesR',         FATAL)
  call check_var_range(pool%microbesK        , 0.0, HUGE(1.0), tag, 'microbesK',         FATAL)
  call check_var_range(pool%DOC              , 0.0, HUGE(1.0), tag, 'DOC',               FATAL)

  select type(pool)
  type is (GIMICS_BGC_litt)
     call check_var_range(pool%dz            , 0.0, HUGE(1.0), tag, 'dz',                FATAL)
  type is (GIMICS_BGC_pool)
     ! do nothing: soil layers dz is fixed
  end select

end subroutine check_GIMICS_pool

! ============================================================================
!> prints debug information for GIMICS pool
subroutine debug_GIMICS_pool(pool, tag)
  class(GIMICS_BGC_pool), intent(in) :: pool
  character(*),           intent(in) :: tag

  write (*,'(a16,":")',advance='NO') trim(tag)
  call dpri('mtbC',pool%metabolicLitterC)
  call dpri('strC',pool%structuralLitterC)
  call dpri('protC',pool%protectedC)
  call dpri('chemRC',pool%chemResistantC)
  call dpri('avlC',pool%availableC)
  call dpri('mR',pool%microbesR)
  call dpri('mK',pool%microbesK)
  call dpri('DOC',pool%DOC)
  select type(pool)
  type is (GIMICS_BGC_litt)
     call dpri('dz',pool%dz)
  type is (GIMICS_BGC_pool)
     ! do nothing
  end select
  write(*,*)
end subroutine debug_GIMICS_pool

! ============================================================================
!> @brief Update a soil carbon pools by crio/bio turbation processes in the soil
subroutine turbation(litt, rhiz, bulk, fRhiz, dz_litt, K_turb, theta, porosity, moist, id_turb_tend, id_litt_tend, diag, tag, &
    allow_flux_to_sfc_litter)
  real, intent(inout) :: litt(N_LITTER_POOLS)  !< concentration in litter(s), [kg/m3]
  real, intent(inout) :: rhiz(:)  !< concentration in rhizosphere, [kg/m3]
  real, intent(inout) :: bulk(:)  !< concentration in bulk soil, [kg/m3]
  real, intent(in)    :: fRhiz(:) !< fraction of rhizosphere, [m3/m3]
  real, intent(in)    :: dz_litt(N_LITTER_POOLS) !< litter thickness, [m]
  real, intent(in)    :: K_turb(:)
  real, intent(in)    :: theta(:)         !< volume of water per volume of air [m3/m3]
  real, intent(in)    :: porosity(:)      !< volume of air per volume of soil [m3/m3]
  real, intent(in)    :: moist(:)         !< volume of water per volume of soil [m3/m3]
  integer, intent(in) :: id_turb_tend !< diagnostic id for turbation tendency field
  integer, intent(in) :: id_litt_tend(N_LITTER_POOLS) !< diagnostic ids for turbation tendencies in litter
  type(diag_buff_type), intent(inout) :: diag !< diagnostic buffer
  character(*), intent(in) :: tag !< textual tag for error messages
  logical, intent(in), optional :: allow_flux_to_sfc_litter !< if TRUE, turbation flux
                      !! from soil to surface litter is allowed; otherwise only flux from
                      !! fitter to soil is allowed. Default if FALSE

  real, dimension(size(rhiz)) :: &
     c,     & ! average concentration in layer, [kg/m3]
     tend,  & ! tendency due to turbation, [kg/(m3 yr)]
     turb_limit, moist_op  !


  real :: f   ! proportionality factor for negative tendency application, unitless
  integer :: k
  real :: d   ! distance between centers of litter and soil layers, [m]
  real :: sfcFlux(N_LITTER_POOLS) ! flux from each litter pool to soil, [kg/(m2 yr)]
  logical :: allow_flux_to_litt ! allow flux from soil to litter


  allow_flux_to_litt = .FALSE.
  if (present (allow_flux_to_sfc_litter)) allow_flux_to_litt = allow_flux_to_sfc_litter



!do k = 1,size(rhiz)
!        moist_op(k) = 0.65 * porosity(k)
!
!        if (moist(k) .lt. moist_op(k)) then
!            turb_limit(k) = 1.0
!        else
!            turb_limit(k) = ((porosity(k) - moist(k))/(porosity(k) - moist_op(k)))**0.75
!        endif
!enddo



do k = 1,size(rhiz)

        if (theta(k) .lt. theta_thres) then
            turb_limit(k) = 1.0
        else
            turb_limit(k) = 0.0
        endif

enddo



  ! calculate average concentration in soil
  do k = 1,size(rhiz)
     c(k) = rhiz(k)*fRhiz(k) + bulk(k)*(1-fRhiz(k))
  enddo
  ! calculate fluxes from litter to soil
  do k = 1,N_LITTER_POOLS
     ! slm: it is questionable if we should use dz_litt to calculate gradient: presumably
     !      diffusion between litter and soil should not decrease as the litter gets
     !      thicker, but this formulation would make it so, because we divide by dz_litt.
     !      E.g., difusion from 20cm-thick liter would be almost 4 times slower than from
     !      5cm. Perhaps we should impose some maximum on d?
     d    = dz_litt(k)/2 + zfull(1) ! distance between centers of litter and top soil layer, [m]
     sfcFlux(k) = K_sfc_turb*(litt(k)-c(1))/d ! flux from litter to soil, [kg/(m2 yr)]
     if (.not.allow_flux_to_litt) then
        sfcFlux(k) = max(sfcFlux(k),0.0) ! disallow fluxes from soil to litter
     endif
     if (is_watch_point()) then
        write(*,'(a20,"(",a3,"):")',advance='NO') trim(tag),trim(l_diagname(k))
        __DEBUG5__(litt(k),c(1),dz_litt(k),d,sfcFlux(k)*dt_fast_yr)
     endif
     sfcFlux(k) = min(sfcFlux(k),litt(k)*dz_litt(k)/dt_fast_yr) ! to avoid depleting litter below zero
  enddo

  ! diffusion in soil
  call diffusion(c, K_turb, turb_limit, sum(sfcFlux), tend)

  !do k = 1,size(rhiz)
  !tend(k) = tend(k)*turb_limit(k)
  !enddo

  ! update soil concentrations
  do k = 1, size(rhiz)
     ! We apply turbation tendency differently depending on its sign: if a bug
     ! goes through the soil and consumed matter (negative tendency), it
     ! presumably does so proportionally to the concentrations in each of the soil
     ! pieces it encounters (rhiz or bulk). On the other hand, when it deposits
     ! carbon (positive tendency), it presumably drops the same concentration in
     ! rhiz or bulk.
     !
     ! The true mechanism of the exchange is probably much more complicated, as
     ! the consumption/dropping happens simultaneously; it also could be very
     ! different for cryoturbation
     if (tend(k).ge.0) then
        ! positive tendency: apply the same concentration increase to rhizosphere and
        ! bulk soil
        rhiz(k) = rhiz(k) + tend(k)*dt_fast_yr
        bulk(k) = bulk(k) + tend(k)*dt_fast_yr
     else
        ! negative tendency: reduce rhizosphere and bulk soil concentration proportionally
        f = (c(k)+tend(k)*dt_fast_yr)/c(k)
        rhiz(k) = rhiz(k)*f
        bulk(k) = bulk(k)*f
     endif
  enddo
  ! update litter concentrations
  do k = 1, N_LITTER_POOLS
     litt(k) = litt(k) - sfcFlux(k)/dz_litt(k)*dt_fast_yr
  enddo

  ! send soil tendencies to diagnostics
  call send_tile_data(id_turb_tend, tend, diag)
  ! send soil tendencies to diagnostics
  do k = 1,N_LITTER_POOLS
     call send_tile_data(id_litt_tend(k), -sfcFlux(k)/dz_litt(k), diag)
  enddo
  ! slm: possibly accumulate tendency for equilibrium concentrations

  ! Detect the situation when diffusion tendency leads to negative concentrations.
  call check_var_range(rhiz, 0.0, HUGE(1.0), 'after turbation of '//trim(tag), 'rhiz', FATAL)
  call check_var_range(bulk, 0.0, HUGE(1.0), 'after turbation of '//trim(tag), 'bulk', FATAL)
  do k = 1,N_LITTER_POOLS
     call check_var_range(litt(k), 0.0, HUGE(1.0), 'after turbation of '//trim(tag), trim(l_diagname(k))//'litt', FATAL)
  enddo
end subroutine



!DOC diffusion is not limited by soil water saturation
! ============================================================================
!> @brief Update DOC diffusion processes in the soil
subroutine turbation2(litt, rhiz, bulk, fRhiz, dz_litt, K_turb, theta, theta_ice, water_ice_porosity, porosity, moist, id_turb_tend, id_litt_tend, diag, tag, &
    allow_flux_to_sfc_litter)
  real, intent(inout) :: litt(N_LITTER_POOLS)  !< concentration in litter(s), [kg/m3]
  real, intent(inout) :: rhiz(:)  !< concentration in rhizosphere, [kg/m3]
  real, intent(inout) :: bulk(:)  !< concentration in bulk soil, [kg/m3]
  real, intent(in)    :: fRhiz(:) !< fraction of rhizosphere, [m3/m3]
  real, intent(in)    :: dz_litt(N_LITTER_POOLS) !< litter thickness, [m]
  real, intent(in)    :: K_turb(:)
  real, intent(in)    :: theta(:)         !< volume of water per volume of air [m3/m3]
  real, intent(in)    :: theta_ice(:)     !< volume of ice per volume of air [m3/m3]
  real, intent(in)    :: water_ice_porosity(:)     !< volume of water and ice per volume of air [m3/m3]
  real, intent(in)    :: porosity(:)      !< volume of air per volume of soil [m3/m3]
  real, intent(in)    :: moist(:)         !< volume of water per volume of soil [m3/m3]
  integer, intent(in) :: id_turb_tend !< diagnostic id for turbation tendency field
  integer, intent(in) :: id_litt_tend(N_LITTER_POOLS) !< diagnostic ids for turbation tendencies in litter
  type(diag_buff_type), intent(inout) :: diag !< diagnostic buffer
  character(*), intent(in) :: tag !< textual tag for error messages
  logical, intent(in), optional :: allow_flux_to_sfc_litter !< if TRUE, turbation flux
                      !! from soil to surface litter is allowed; otherwise only flux from
                      !! fitter to soil is allowed. Default if FALSE

  real, dimension(size(rhiz)) :: &
     c,     & ! average concentration in water, [kg/m3]
     cc,     & ! average concentration in layer, [kg/m3]
     tend,  & ! tendency due to turbation, [kg/(m3 yr)]
     turb_limit, moist_op, &  !
     h2o, h2o2


  real :: f   ! proportionality factor for negative tendency application, unitless
  integer :: k
  real :: d   ! distance between centers of litter and soil layers, [m]
  real :: sfcFlux(N_LITTER_POOLS) ! flux from each litter pool to soil, [kg/(m2 yr)]
  logical :: allow_flux_to_litt ! allow flux from soil to litter

  allow_flux_to_litt = .FALSE.
  if (present (allow_flux_to_sfc_litter)) allow_flux_to_litt = allow_flux_to_sfc_litter

  do k = 1,size(rhiz)
     turb_limit(k) = 1.0
  enddo

  ! calculate average concentration in soil
  do k = 1,size(rhiz)
     cc(k) = rhiz(k)*fRhiz(k) + bulk(k)*(1-fRhiz(k))

     h2o(k) = theta(k) * porosity(k)
     h2o2(k) = water_ice_porosity(k) * porosity(k)

     if (h2o(k) .gt. 0.0) then
        c(k) = rhiz(k)/h2o(k)*fRhiz(k) + bulk(k)/h2o(k)*(1-fRhiz(k))
     else
        c(k)=0
     endif
  enddo
  ! calculate fluxes from litter to soil
  do k = 1,N_LITTER_POOLS
     ! slm: it is questionable if we should use dz_litt to calculate gradient: presumably
     !      diffusion between litter and soil should not decrease as the litter gets
     !      thicker, but this formulation would make it so, because we divide by dz_litt.
     !      E.g., difusion from 20cm-thick liter would be almost 4 times slower than from
     !      5cm. Perhaps we should impose some maximum on d?
     d    = dz_litt(k)/2 + zfull(1) ! distance between centers of litter and top soil layer, [m]

     !if ((h2o2(1) .gt. 0.0) .and. (h2o(1) .gt. 0.0)) then
     if (h2o(1) .gt. 0.0) then
        sfcFlux(k) = (   2*K_turb(1)/1.5 * (h2o(1)*h2o(1)/ (h2o(1)*dz_litt(k)+h2o(1)*dz(1)) )   )*(litt(k)/h2o(1)-c(1)) ! flux from litter to soil, [kg/(m2 yr)]
        !(   2*K_turb(1)/1.5 * (h2o(1)*h2o(1)/ (h2o(1)*dz_litt(k)+h2o(1)*dz(1)) )   )

        if (.not.allow_flux_to_litt) then
           sfcFlux(k) = max(sfcFlux(k),0.0) ! disallow fluxes from soil to litter
        endif
        if (is_watch_point()) then
           write(*,'(a20,"(",a3,"):")',advance='NO') trim(tag),trim(l_diagname(k))
           __DEBUG5__(litt(k)/h2o(1),c(1),dz_litt(k),d,sfcFlux(k)*dt_fast_yr)
           __DEBUG5__(h2o(1),h2o2(1),theta(1),theta_ice(1),water_ice_porosity(1))
        endif
        sfcFlux(k) = min(sfcFlux(k),litt(k)*dz_litt(k)/dt_fast_yr) ! to avoid depleting litter below zero
     else
        sfcFlux(k) = 0.0
     endif
  enddo

  ! diffusion in soil
  call diffusion2(c, h2o, h2o2, theta_ice, K_turb, turb_limit, sum(sfcFlux), tend)

  !do k = 1,size(rhiz)
  !tend(k) = tend(k)*turb_limit(k)
  !enddo

  ! update soil concentrations
  do k = 1, size(rhiz)
     ! We apply turbation tendency differently depending on its sign: if a bug
     ! goes through the soil and consumed matter (negative tendency), it
     ! presumably does so proportionally to the concentrations in each of the soil
     ! pieces it encounters (rhiz or bulk). On the other hand, when it deposits
     ! carbon (positive tendency), it presumably drops the same concentration in
     ! rhiz or bulk.
     !
     ! The true mechanism of the exchange is probably much more complicated, as
     ! the consumption/dropping happens simultaneously; it also could be very
     ! different for cryoturbation
     if (tend(k).ge.0) then
        ! positive tendency: apply the same concentration increase to rhizosphere and
        ! bulk soil
        rhiz(k) = rhiz(k) + tend(k)*dt_fast_yr
        bulk(k) = bulk(k) + tend(k)*dt_fast_yr
     else
        ! negative tendency: reduce rhizosphere and bulk soil concentration proportionally
        f = (cc(k)+tend(k)*dt_fast_yr)/cc(k)
        rhiz(k) = rhiz(k)*f
        bulk(k) = bulk(k)*f
     endif
  enddo
  ! update litter concentrations
  do k = 1, N_LITTER_POOLS
     litt(k) = litt(k) - sfcFlux(k)/dz_litt(k)*dt_fast_yr
  enddo

  ! send soil tendencies to diagnostics
  call send_tile_data(id_turb_tend, tend, diag)
  ! send soil tendencies to diagnostics
  do k = 1,N_LITTER_POOLS
     call send_tile_data(id_litt_tend(k), -sfcFlux(k)/dz_litt(k), diag)
  enddo
  ! slm: possibly accumulate tendency for equilibrium concentrations

  ! Detect the situation when diffusion tendency leads to negative concentrations.
  call check_var_range(rhiz, 0.0, HUGE(1.0), 'after turbation of '//trim(tag), 'rhiz', FATAL)
  call check_var_range(bulk, 0.0, HUGE(1.0), 'after turbation of '//trim(tag), 'bulk', FATAL)
  do k = 1,N_LITTER_POOLS
     call check_var_range(litt(k), 0.0, HUGE(1.0), 'after turbation of '//trim(tag), trim(l_diagname(k))//'litt', FATAL)
  enddo
end subroutine

! ============================================================================
!> @brief Calculate tendency due to vertical bioturbation
subroutine diffusion(C,D,turb_limit,F0,tend)
  real, intent(in)  :: C(:) !< transported quantity, by layer, [kg/m3]
  real, intent(in)  :: D(:) !< coefficients of diffusion (at the layer's bottom), [m2/yr]
  real, intent(in)  :: turb_limit(:) !<
  real, intent(in)  :: F0   !< flux into the soil at the soil surface, [kg/(m2 yr)]
  real, intent(out) :: tend(:) !< tendencies due to diffusion [kg/(m3 yr)]

  integer :: k
  real    :: flux(0:num_l) ! flux at the lower boundary of the soil layer [kg/(m2 yr)],
      ! positive downward. flux(0) is on the top of the soil, flux(num_l) -- at the soil bottom
  flux(0) = F0
  do k = 1, num_l-1
     flux(k) = D(k)*(C(k)-C(k+1))/(zfull(k+1)-zfull(k))*turb_limit(k)
  enddo
  flux(num_l) = 0.0
  do k = 1,num_l
     tend(k) = ((flux(k-1)-flux(k))/dz(k))
  enddo
end subroutine

! ============================================================================
!> @brief Calculate tendency due to vertical diffusion of DOC
subroutine diffusion2(C,h2o,h2o2,theta_ice, D,turb_limit,F0,tend)
  real, intent(in)  :: C(:) !< DOC per water, [kg/m3]
  real, intent(in)  :: theta_ice(:) !< soil ice per air, by layer, [m3/m3]
  real, intent(in)  :: h2o(:) !<soil water content, [m3/m3]
  real, intent(in)  :: h2o2(:) !<soil water and ice content, [m3/m3]
  real, intent(in)  :: D(:) !< coefficients of diffusion (at the layer's bottom), [m2/yr]
  real, intent(in)  :: turb_limit(:) !<
  real, intent(in)  :: F0   !< flux into the soil at the soil surface, [kg/(m2 yr)]
  real, intent(out) :: tend(:) !< tendencies due to diffusion [kg/(m3 yr)]

  integer :: k
  real    :: flux(0:num_l) ! flux at the lower boundary of the soil layer [kg/(m2 yr)],
      ! positive downward. flux(0) is on the top of the soil, flux(num_l) -- at the soil bottom
  flux(0) = F0

  do k = 1, num_l-1
     !flux(k) = D(k)*(C(k)-C(k+1))/(zfull(k+1)-zfull(k))*turb_limit(k)

     !if ((h2o2(k) .gt. 0.0) .and.  (h2o(k) .gt. 0.0)) then
     if (h2o(k) .gt. 0.0) then

     flux(k) = (   2*D(k)/1.5 * (h2o(k)*h2o(k+1)/ (h2o(k+1)*dz(k)+h2o(k)*dz(k+1)) )   )*(C(k)-C(k+1))*turb_limit(k)

     else
     flux(k) = 0.0

     !(   2*D(k)/1.5 * (h2o(k)*h2o(k+1)/ (h2o(k+1)*dz(k)+h2o(k)*dz(k+1)) )   )
     endif

  enddo

  flux(num_l) = 0.0
  do k = 1,num_l
     tend(k) = ((flux(k-1)-flux(k))/dz(k))
  enddo
end subroutine




! ============================================================================
subroutine step3_GIMICS(soilc, diag)
  class(soil_BGC_GIMICS_t),   intent(inout) :: soilc
  type(diag_buff_type),       intent(inout) :: diag

  integer :: k
  real :: s
  real :: layer_C(num_l) ! total carbon amount per layer, kgC/m2
  real :: rhiz(num_l), bulk(num_l) ! total carbon density in bulk soil and rhizosphere, kgC/m3

  if (id_total_soil_C>0) call send_tile_data(id_total_soil_C, soilc%total_C(), diag)

  if (id_fRhiz > 0) call send_tile_data(id_fRhiz, soilc%fRhiz, diag)

  ! NOTE that IDs in calls to send_3_tile_data are arrays, so protecting them
  !      with "if" statements would be more involved

  ! vertical distribution of carbon pools
  call send_3_tile_data(id_metabolicC,     soilc%rhiz(:)%metabolicLitterC,  soilc%bulk(:)%metabolicLitterC,  soilc%fRhiz(:), diag)
  call send_3_tile_data(id_structuralC,    soilc%rhiz(:)%structuralLitterC, soilc%bulk(:)%structuralLitterC, soilc%fRhiz(:), diag)
  call send_3_tile_data(id_protectedC,     soilc%rhiz(:)%protectedC,        soilc%bulk(:)%protectedC,        soilc%fRhiz(:), diag)
  call send_3_tile_data(id_chemResistantC, soilc%rhiz(:)%chemResistantC,    soilc%bulk(:)%chemResistantC,    soilc%fRhiz(:), diag)
  call send_3_tile_data(id_availableC,     soilc%rhiz(:)%availableC,        soilc%bulk(:)%availableC,        soilc%fRhiz(:), diag)
  call send_3_tile_data(id_microbesR,      soilc%rhiz(:)%microbesR,         soilc%bulk(:)%microbesR,         soilc%fRhiz(:), diag)
  call send_3_tile_data(id_microbesK,      soilc%rhiz(:)%microbesK,         soilc%bulk(:)%microbesK,         soilc%fRhiz(:), diag)
  call send_3_tile_data(id_DOC,            soilc%rhiz(:)%DOC,               soilc%bulk(:)%DOC,               soilc%fRhiz(:), diag)

  ! totals of carbon pools
  if (id_clmn_metabolicC > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%metabolicLitterC, soilc%bulk(:)%metabolicLitterC)
     call send_tile_data(id_clmn_metabolicC, s, diag)
  endif
  if (id_clmn_structuralC > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%structuralLitterC, soilc%bulk(:)%structuralLitterC)
     call send_tile_data(id_clmn_structuralC, s, diag)
  endif
  if (id_clmn_protectedC > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%protectedC, soilc%bulk(:)%protectedC)
     call send_tile_data(id_clmn_protectedC, s, diag)
  endif
  if (id_clmn_chemResistantC > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%chemResistantC, soilc%bulk(:)%chemResistantC)
     call send_tile_data(id_clmn_chemResistantC, s, diag)
  endif
  if (id_clmn_availableC > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%availableC, soilc%bulk(:)%availableC)
     call send_tile_data(id_clmn_availableC, s, diag)
  endif
  if (id_clmn_microbesR > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%microbesR, soilc%bulk(:)%microbesR)
     call send_tile_data(id_clmn_microbesR, s, diag)
  endif
  if (id_clmn_microbesK > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%microbesK, soilc%bulk(:)%microbesK)
     call send_tile_data(id_clmn_microbesK, s, diag)
  endif
  if (id_clmn_DOC > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%DOC, soilc%bulk(:)%DOC)
     call send_tile_data(id_clmn_DOC, s, diag)
  endif

  ! decomposition rates
  call send_3_tile_data(id_DecompMrLm,     soilc%rhiz(:)%DecompMrLm,        soilc%bulk(:)%DecompMrLm,        soilc%fRhiz(:), diag)
  call send_3_tile_data(id_DecompMrLs,     soilc%rhiz(:)%DecompMrLs,        soilc%bulk(:)%DecompMrLs,        soilc%fRhiz(:), diag)
  call send_3_tile_data(id_DecompMrCa,     soilc%rhiz(:)%DecompMrCa,        soilc%bulk(:)%DecompMrCa,        soilc%fRhiz(:), diag)
  call send_3_tile_data(id_DecompMrDOC,    soilc%rhiz(:)%DecompMrDOC,       soilc%bulk(:)%DecompMrDOC,       soilc%fRhiz(:), diag)
  call send_3_tile_data(id_OxidMrCc,       soilc%rhiz(:)%OxidMrCc,          soilc%bulk(:)%OxidMrCc,          soilc%fRhiz(:), diag)
  call send_3_tile_data(id_MrTau,          soilc%rhiz(:)%MrTau,             soilc%bulk(:)%MrTau,             soilc%fRhiz(:), diag)

  call send_3_tile_data(id_DecompMkLm,     soilc%rhiz(:)%DecompMkLm,        soilc%bulk(:)%DecompMkLm,        soilc%fRhiz(:), diag)
  call send_3_tile_data(id_DecompMkLs,     soilc%rhiz(:)%DecompMkLs,        soilc%bulk(:)%DecompMkLs,        soilc%fRhiz(:), diag)
  call send_3_tile_data(id_DecompMkCa,     soilc%rhiz(:)%DecompMkCa,        soilc%bulk(:)%DecompMkCa,        soilc%fRhiz(:), diag)
  call send_3_tile_data(id_DecompMkDOC,    soilc%rhiz(:)%DecompMkDOC,       soilc%bulk(:)%DecompMkDOC,       soilc%fRhiz(:), diag)
  call send_3_tile_data(id_OxidMkCc,       soilc%rhiz(:)%OxidMkCc,          soilc%bulk(:)%OxidMkCc,          soilc%fRhiz(:), diag)
  call send_3_tile_data(id_MkTau,          soilc%rhiz(:)%MkTau,             soilc%bulk(:)%MkTau,             soilc%fRhiz(:), diag)

  call send_3_tile_data(id_Resp,           soilc%rhiz(:)%Resp,              soilc%bulk(:)%Resp,              soilc%fRhiz(:), diag)
  call send_3_tile_data(id_Desorb,         soilc%rhiz(:)%Desorb,            soilc%bulk(:)%Desorb,            soilc%fRhiz(:), diag)
  call send_3_tile_data(id_thetaF,         soilc%rhiz(:)%thetaF,            soilc%bulk(:)%thetaF,            soilc%fRhiz(:), diag)

  call send_3_tile_data(id_InputStrC,      soilc%rhiz(:)%InputStrC,         soilc%bulk(:)%InputStrC,         soilc%fRhiz(:), diag)
  call send_3_tile_data(id_InputMtbC,      soilc%rhiz(:)%InputMtbC,         soilc%bulk(:)%InputMtbC,         soilc%fRhiz(:), diag)
  call send_3_tile_data(id_InputExdC,      soilc%rhiz(:)%InputExdC,         soilc%bulk(:)%InputExdC,         soilc%fRhiz(:), diag)

  if (id_clmn_Resp > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%Resp, soilc%bulk(:)%Resp)
     call send_tile_data(id_clmn_Resp, s, diag)
  endif
  if (id_clmn_Desorb > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%Desorb, soilc%bulk(:)%Desorb)
     call send_tile_data(id_clmn_Desorb, s, diag)
  endif

  if (id_clmn_InputStrC > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%InputStrC, soilc%bulk(:)%InputStrC)
     call send_tile_data(id_clmn_InputStrC, s, diag)
  endif
  if (id_clmn_InputMtbC > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%InputMtbC, soilc%bulk(:)%InputMtbC)
     call send_tile_data(id_clmn_InputMtbC, s, diag)
  endif
  if (id_clmn_InputExdC > 0) then
     s = total_amount(soilc%fRhiz(:), soilc%rhiz(:)%InputExdC, soilc%bulk(:)%InputExdC)
     call send_tile_data(id_clmn_InputExdC, s, diag)
  endif

  if (id_clmn_Decomp > 0) then
     s = total_amount(soilc%fRhiz(:), &
        soilc%rhiz(:)%DecompMrLm  + soilc%rhiz(:)%DecompMkLm + &
        soilc%rhiz(:)%DecompMrLs  + soilc%rhiz(:)%DecompMkLs + &
        soilc%rhiz(:)%DecompMrCa  + soilc%rhiz(:)%DecompMkCa + &
        soilc%rhiz(:)%DecompMrDOC + soilc%rhiz(:)%DecompMkDOC, &

        soilc%bulk(:)%DecompMrLm  + soilc%bulk(:)%DecompMkLm + &
        soilc%bulk(:)%DecompMrLs  + soilc%bulk(:)%DecompMkLs + &
        soilc%bulk(:)%DecompMrCa  + soilc%bulk(:)%DecompMkCa + &
        soilc%bulk(:)%DecompMrDOC + soilc%bulk(:)%DecompMkDOC  )
     call send_tile_data(id_clmn_Decomp, s, diag)
  endif
  if (id_clmn_Oxid > 0) then
     s = total_amount(soilc%fRhiz(:), &
        soilc%rhiz(:)%OxidMrCc + soilc%rhiz(:)%OxidMkCc, &
        soilc%bulk(:)%OxidMrCc + soilc%bulk(:)%OxidMkCc  )
     call send_tile_data(id_clmn_Oxid, s, diag)
  endif
  if (id_clmn_Turnover > 0) then
     s = total_amount(soilc%fRhiz(:), &
        soilc%rhiz(:)%MrTau + soilc%rhiz(:)%MkTau, &
        soilc%bulk(:)%MrTau + soilc%bulk(:)%MkTau  )
     call send_tile_data(id_clmn_Turnover, s, diag)
  endif

  ! reset input accumulators for the next time step
  soilc%rhiz(:)%InputStrC = 0.0; soilc%bulk(:)%InputStrC = 0.0
  soilc%rhiz(:)%InputMtbC = 0.0; soilc%bulk(:)%InputMtbC = 0.0
  soilc%rhiz(:)%InputExdC = 0.0; soilc%bulk(:)%InputExdC = 0.0

  do k = 1, N_LITTER_POOLS
     if (id_litt_total_C(k)>0) call send_tile_data(id_litt_total_C(k), C_amount(soilc%litt(k)),  diag)
     if (id_litt_allC(k)>0)    call send_tile_data(id_litt_allC(k),    C_density(soilc%litt(k)), diag)
     call send_tile_data(id_litt_dz(k),             soilc%litt(k)%dz,                diag)

     ! surface litter carbon pools
     call send_tile_data(id_litt_metabolicC(k),     soilc%litt(k)%metabolicLitterC,  diag)
     call send_tile_data(id_litt_structuralC(k),    soilc%litt(k)%structuralLitterC, diag)
     call send_tile_data(id_litt_chemResistantC(k), soilc%litt(k)%chemResistantC,    diag)
     call send_tile_data(id_litt_availableC(k),     soilc%litt(k)%availableC,        diag)
     call send_tile_data(id_litt_microbesR(k),      soilc%litt(k)%microbesR,         diag)
     call send_tile_data(id_litt_microbesK(k),      soilc%litt(k)%microbesK,         diag)
     call send_tile_data(id_litt_DOC(k),            soilc%litt(k)%DOC,               diag)

     ! decomposition rates
     call send_tile_data(id_litt_DecompMrLm(k), soilc%litt(k)%DecompMrLm, diag)
     call send_tile_data(id_litt_DecompMrLs(k), soilc%litt(k)%DecompMrLs, diag)
     call send_tile_data(id_litt_DecompMrCa(k), soilc%litt(k)%DecompMrCa, diag)
     call send_tile_data(id_litt_DecompMrDOC(k),soilc%litt(k)%DecompMrDOC, diag)
     call send_tile_data(id_litt_OxidMrCc(k),   soilc%litt(k)%OxidMrCc,   diag)
     call send_tile_data(id_litt_MrTau(k),      soilc%litt(k)%MrTau,      diag)

     call send_tile_data(id_litt_DecompMkLm(k), soilc%litt(k)%DecompMkLm, diag)
     call send_tile_data(id_litt_DecompMkLs(k), soilc%litt(k)%DecompMkLs, diag)
     call send_tile_data(id_litt_DecompMkCa(k), soilc%litt(k)%DecompMkCa, diag)
     call send_tile_data(id_litt_DecompMkDOC(k),soilc%litt(k)%DecompMkDOC, diag)
     call send_tile_data(id_litt_OxidMkCc(k),   soilc%litt(k)%OxidMkCc,   diag)
     call send_tile_data(id_litt_MkTau(k),      soilc%litt(k)%MkTau,      diag)

     if (id_litt_Decomp(k) > 0) then
        call send_tile_data(id_litt_Decomp(k), &
           soilc%litt(k)%DecompMrLm  + soilc%litt(k)%DecompMkLm + &
           soilc%litt(k)%DecompMrLs  + soilc%litt(k)%DecompMkLs + &
           soilc%litt(k)%DecompMrCa  + soilc%litt(k)%DecompMkCa + &
           soilc%litt(k)%DecompMrDOC + soilc%litt(k)%DecompMkDOC, diag)
     endif
     if (id_litt_Oxid(k) > 0) then
        call send_tile_data(id_litt_Oxid(k), &
           soilc%litt(k)%OxidMrCc + soilc%litt(k)%OxidMkCc, diag)
     endif
     if (id_litt_Turnover(k) > 0) then
        call send_tile_data(id_litt_Turnover(k), &
           soilc%litt(k)%MrTau + soilc%litt(k)%MkTau, diag)
     endif
     call send_tile_data(id_litt_Resp(k),       soilc%litt(k)%Resp,       diag)
     call send_tile_data(id_litt_thetaF(k),     soilc%litt(k)%thetaF,     diag)

     call send_tile_data(id_litt_InputStrC(k),  soilc%litt(k)%InputStrC,  diag)
     call send_tile_data(id_litt_InputMtbC(k),  soilc%litt(k)%InputMtbC,  diag)
     ! reset input accumulators for the next timestep
     soilc%litt(k)%InputStrC = 0.0; soilc%litt(k)%InputMtbC = 0.0
  enddo

  do k = 1,num_l
     rhiz(k) = C_density(soilc%rhiz(k))
     bulk(k) = C_density(soilc%bulk(k))
     layer_C(k) = (rhiz(k)*soilc%fRhiz(k) + bulk(k)*(1-soilc%fRhiz(k))) * dz(k)
  enddo
  call send_3_tile_data(id_soilC, rhiz(:), bulk(:), soilc%fRhiz(:), diag)

!   call send_tile_data(id_fsc, sum(soil%fast_soil_C(:))+sum(soil%litter_SIMPLE_C(C_FAST,:)), diag)
!   call send_tile_data(id_ssc, sum(soil%slow_soil_C(:))+sum(soil%litter_SIMPLE_C(C_SLOW,:)), diag)
!   call send_tile_data(id_soil_C(C_FAST), soil%fast_soil_C(:)/dz(1:num_l), diag)
!   call send_tile_data(id_soil_C(C_SLOW), soil%slow_soil_C(:)/dz(1:num_l), diag)

  ! --- CMOR vars
  if (id_cSoilLevels  > 0) call send_tile_data(id_cSoilLevels, layer_C, diag)
  if (id_cSoilAbove1m > 0) call send_tile_data(id_cSoilAbove1m, sum(layer_C(:)*mrs1m_weight(:)), diag)

  if (id_cSoil>0) call send_tile_data(id_cSoil, soilc%total_soil_C(), diag)
! slm: in GIMICS, what is fast, medium, and slow carbon?
!   if (id_csoilfast>0)   call send_tile_data(id_csoilfast,   sum(soil%fast_soil_C(:)), diag)
!   if (id_csoilmedium>0) call send_tile_data(id_csoilmedium, sum(soil%slow_soil_C(:)), diag)
!   call send_tile_data(id_csoilslow, 0.0, diag)

  if (id_cLitter>0) then
     s = 0
     do k = 1, N_LITTER_POOLS
        s = s+C_amount(soilc%litt(k))
     enddo
     call send_tile_data(id_cLitter, s, diag)
  endif
  if (id_cLitterCwd>0)  call send_tile_data(id_cLitterCwd,  C_amount(soilc%litt(LITT_CWOOD)), diag)
  if (id_cLitterLeaf>0) call send_tile_data(id_cLitterLeaf, C_amount(soilc%litt(LITT_LEAF)),  diag)
  ! --- end of CMOR vars

end subroutine step3_GIMICS


! ============================================================================
!> @brief Update soil carbon pool
subroutine update_GIMICS_pool(pool, T, theta, porosity, moist, fClay, cw_r, cw_z, lf_f,is_sfc_litter)
  class(GIMICS_BGC_pool), intent(inout) :: pool
  real,    intent(in) :: T             !< Temperature [degC]
  real,    intent(in) :: theta         !< volume of water per volume of air [m3/m3]
  real,    intent(in) :: porosity      !< volume of air per volume of soil [m3/m3]
  real,    intent(in) :: moist         !< volume of water per volume of soil [m3/m3]
  real,    intent(in) :: fClay         !< clay fraction, unitless, within [0,1] interval
  real,    intent(in) :: cw_r          !< coarse wood radius [cm]
  real,    intent(in) :: cw_z          !< coarse wood thickness accessible by microbes [cm]
  real,    intent(in) :: lf_f          !< leaf litter fraction accessible by microbes
  logical, intent(in) :: is_sfc_litter !< TRUE is the pool is surface litter: protected C is always zero in this case

  real :: Vmax_Mr_Lm, Vmax_Mr_Ls, Vmax_Mr_Ca, Vmax_Mr_DOC, Vmax_Mk_Lm, Vmax_Mk_Ls, Vmax_Mk_Ca, Vmax_Mk_DOC, &
          Km_Mr_Lm,   Km_Mr_Ls,   Km_Mr_Ca,   Km_Mr_DOC, Km_Mk_Lm,   Km_Mk_Ls,   Km_Mk_Ca,   Km_Mk_DOC, &
          fMrTau_Cp, fMkTau_Cp, fMrTau_Cc, fMkTau_Cc

  real :: Vmax_base

  real :: moist_op, yan_theta_a

  moist_op = 0.65 * porosity

  if (fClay .le. 0.016) then
     yan_theta_a = 0.0

  else if (fClay .gt. 0.37) then
     yan_theta_a = 1.0

  else
     yan_theta_a = 2.8 * fClay - 0.046
  endif

  if (is_sfc_litter) then
     select case(theta_func_litt_option)
     case (THETA_F_ORCHIDEE)
        ! ORCHIDEE moisture function
        pool%thetaF = theta_func_orchidee(moist,theta_func_orchidee_min,theta_func_orchidee_max)
     case (THETA_F_YAN2018)
        ! Generalized, mechanistic soil moisture function from Yan et al. (2018)
        if (moist .lt. moist_op) then
            pool%thetaF = ((0.1 + moist_op)/(0.1 + moist)) * ((moist/moist_op)**(1+yan_theta_a*2))
        else
            pool%thetaF = ((porosity - moist)/(porosity - moist_op))**0.75
        endif
     case (THETA_F_CORPSE)
        ! CORPSE moisture function
        pool%thetaF = theta_func(theta*litt_theta_mod,1-theta*litt_theta_mod,substrate_diffusion_exp,gas_diffusion_exp,min_anaerobic_resp_factor, min_dry_resp_factor)
     case (THETA_F_NONE)
        pool%thetaF = 1.0
     case default
        call land_error_message('update_GIMICS_pool: incorrect value of theta_func_litt_option', FATAL)
     end select
  else
     select case(theta_func_soil_option)
     ! slm: should we have ORCHIDEE option for soil too?


     case (THETA_F_YAN2018)
     !Generalized, mechanistic soil moisture function from Yan et al. (2018)
        if (moist .lt. moist_op) then
            pool%thetaF = max(ThetaFmin_dry,((0.1 + moist_op)/(0.1 + moist)) * ((moist/moist_op)**(1+yan_theta_a*2)))
        else
            pool%thetaF = max(ThetaFmin_wet,((porosity - moist)/(porosity - moist_op))**0.75)
        endif


     case (THETA_F_CORPSE)
        ! CORPSE moisture function
        pool%thetaF = theta_func(theta*litt_theta_mod,1-theta*litt_theta_mod,substrate_diffusion_exp,gas_diffusion_exp,min_anaerobic_resp_factor, min_dry_resp_factor)



     case (THETA_F_A)

     if(theta <= 0.3) then
     pool%thetaF = 0.2;
     else if(theta <= 0.6) then
     pool%thetaF = 0.2+0.8*(theta-0.3)/0.3;
     else
     pool%thetaF = exp(2.3*(0.6-theta));
     endif


     case (THETA_F_NONE)
        pool%thetaF = 1.0
     case default
        call land_error_message('update_GIMICS_pool: incorrect value of theta_func_soil_option', FATAL)
     end select
  endif

  if (highT_limit) then

    if (lowT_limit) then
    Vmax_base =  pool%thetaF * max( exp(Vint)*aV, exp(Vslope*T+Vint)*aV) * (1/(1+exp(0.4*(T-45))))
    else
    Vmax_base =  pool%thetaF * exp(Vslope*T+Vint)*aV * (1/(1+exp(0.4*(T-45))))
    endif
  else

    if (lowT_limit) then
    Vmax_base =  pool%thetaF * max( exp(Vint)*aV, exp(Vslope*T+Vint)*aV)
    else
    Vmax_base =  pool%thetaF * exp(Vslope*T+Vint)*aV
    endif
  endif

  Vmax_Mr_Lm = Vmax_base * Vmod_Mr_Lm ! mgC/mgM/h
  Vmax_Mr_Ls = Vmax_base * Vmod_Mr_Ls
  Vmax_Mr_Ca = Vmax_base * Vmod_Mr_Ca
  Vmax_Mr_DOC = Vmax_base * Vmod_Mr_DOC

  Vmax_Mk_Lm = Vmax_base * Vmod_Mk_Lm
  Vmax_Mk_Ls = Vmax_base * Vmod_Mk_Ls
  Vmax_Mk_Ca = Vmax_base * Vmod_Mk_Ca
  Vmax_Mk_DOC = Vmax_base * Vmod_Mk_DOC

  if (Km_theta) then
     Km_Mr_Lm = exp(Kslope_Lm*T+Kint) * aK * pool%thetaF * Kmod_Mr_Lm ! kgC/m3
     Km_Mr_Ls = exp(Kslope_Ls*T+Kint) * aK * pool%thetaF * Kmod_Mr_Ls
     Km_Mr_Ca = exp(Kslope_Ca*T+Kint) * aK * pool%thetaF * Kmod_Mr_Ca
     Km_Mr_DOC = exp(Kslope_Ca*T+Kint) * aK * pool%thetaF * Kmod_Mr_Ca / (2.0*exp(-2.0*sqrt(fClay)))

     Km_Mk_Lm = exp(Kslope_Lm*T+Kint) * aK * pool%thetaF * Kmod_Mk_Lm
     Km_Mk_Ls = exp(Kslope_Ls*T+Kint) * aK * pool%thetaF * Kmod_Mk_Ls
     Km_Mk_Ca = exp(Kslope_Ca*T+Kint) * aK * pool%thetaF * Kmod_Mk_Ca
     Km_Mk_DOC = exp(Kslope_Ca*T+Kint) * aK * pool%thetaF * Kmod_Mk_Ca / (2.0*exp(-2.0*sqrt(fClay)))
  else
     Km_Mr_Lm = exp(Kslope_Lm*T+Kint) * aK * Kmod_Mr_Lm ! kgC/m3
     Km_Mr_Ls = exp(Kslope_Ls*T+Kint) * aK * Kmod_Mr_Ls
     Km_Mr_Ca = exp(Kslope_Ca*T+Kint) * aK * Kmod_Mr_Ca
     Km_Mr_DOC = exp(Kslope_Ca*T+Kint) * aK * Kmod_Mr_Ca / (2.0*exp(-2.0*sqrt(fClay)))

     Km_Mk_Lm = exp(Kslope_Lm*T+Kint) * aK * Kmod_Mk_Lm
     Km_Mk_Ls = exp(Kslope_Ls*T+Kint) * aK * Kmod_Mk_Ls
     Km_Mk_Ca = exp(Kslope_Ca*T+Kint) * aK * Kmod_Mk_Ca
     Km_Mk_DOC = exp(Kslope_Ca*T+Kint) * aK * Kmod_Mk_Ca / (2.0*exp(-2.0*sqrt(fClay)))
  endif

  if (is_sfc_litter) then
     ! litter accessible by microbes for decomposition
     !(pool%metabolicLitterC  * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)
     !(pool%structuralLitterC * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)

     pool%DecompMrLm  = Vmax_Mr_Lm * pool%microbesR * ((pool%metabolicLitterC  * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)/(Km_Mr_Lm+(pool%metabolicLitterC  * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)))
     pool%DecompMrLs  = Vmax_Mr_Ls * pool%microbesR * ((pool%structuralLitterC * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)/(Km_Mr_Ls+(pool%structuralLitterC * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)))
     pool%DecompMrCa  = Vmax_Mr_Ca * pool%microbesR * (pool%availableC/(Km_Mr_Ca+pool%availableC))
     pool%DecompMrDOC = Vmax_Mr_DOC * pool%microbesR * (pool%DOC/(Km_Mr_DOC+pool%DOC))
     pool%OxidMrCc    = Vmax_Mr_Ls * pool%microbesR * (pool%chemResistantC/(Kmod_oxid_Mr*Km_Mr_Ls+pool%chemResistantC)) ! kgC/m3/h

     pool%DecompMkLm  = Vmax_Mk_Lm * pool%microbesK * ((pool%metabolicLitterC  * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)/(Km_Mk_Lm+(pool%metabolicLitterC  * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)))
     pool%DecompMkLs  = Vmax_Mk_Ls * pool%microbesK * ((pool%structuralLitterC * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)/(Km_Mk_Ls+(pool%structuralLitterC * (cw_z * (2 * cw_r - cw_z) / cw_r**2) * lf_f)))
     pool%DecompMkCa  = Vmax_Mk_Ca * pool%microbesK * (pool%availableC/(Km_Mk_Ca+pool%availableC))
     pool%DecompMkDOC = Vmax_Mk_DOC * pool%microbesK * (pool%DOC/(Km_Mk_DOC+pool%DOC))
     pool%OxidMkCc    = Vmax_Mk_Ls * pool%microbesK * (pool%chemResistantC/(Kmod_oxid_Mk*Km_Mk_Ls+pool%chemResistantC))
  else
     pool%DecompMrLm  = Vmax_Mr_Lm * pool%microbesR * (pool%metabolicLitterC/(Km_Mr_Lm+pool%metabolicLitterC)) ! kgC/m3/h
     pool%DecompMrLs  = Vmax_Mr_Ls * pool%microbesR * (pool%structuralLitterC/(Km_Mr_Ls+pool%structuralLitterC))
     pool%DecompMrCa  = Vmax_Mr_Ca * pool%microbesR * (pool%availableC/(Km_Mr_Ca+pool%availableC))
     pool%DecompMrDOC = Vmax_Mr_DOC * pool%microbesR * (pool%DOC/(Km_Mr_DOC+pool%DOC))
     pool%OxidMrCc    = Vmax_Mr_Ls * pool%microbesR * (pool%chemResistantC/(Kmod_oxid_Mr*Km_Mr_Ls+pool%chemResistantC)) ! kgC/m3/h

     pool%DecompMkLm  = Vmax_Mk_Lm * pool%microbesK * (pool%metabolicLitterC/(Km_Mk_Lm+pool%metabolicLitterC))
     pool%DecompMkLs  = Vmax_Mk_Ls * pool%microbesK * (pool%structuralLitterC/(Km_Mk_Ls+pool%structuralLitterC))
     pool%DecompMkCa  = Vmax_Mk_Ca * pool%microbesK * (pool%availableC/(Km_Mk_Ca+pool%availableC))
     pool%DecompMkDOC = Vmax_Mk_DOC * pool%microbesK * (pool%DOC/(Km_Mk_DOC+pool%DOC))
     pool%OxidMkCc    = Vmax_Mk_Ls * pool%microbesK * (pool%chemResistantC/(Kmod_oxid_Mk*Km_Mk_Ls+pool%chemResistantC))
  endif


 !if (is_sfc_litter) then

  pool%MrTau = tau_calib * 5.2e-4 * exp(0.3*fI_Lm) * pool%microbesR**tau_beta ! kgC/m3/h
  pool%MkTau = tau_calib * 2.4e-4 * exp(0.1*fI_Lm) * pool%microbesK**tau_beta

  if (is_sfc_litter) then
     fMrTau_Cp = 0.0
     fMkTau_Cp = 0.0

     fMrTau_Cc = min(1.0-fMrTau_Cp, fMrTau_Cc_a4_litt*exp(fMrTau_Cc_a3_litt*fI_Lm))  ! 0.0315
     fMkTau_Cc = min(1.0-fMkTau_Cp, fMrTau_Cc_a5_litt*exp(fMrTau_Cc_a3_litt*fI_Lm))  ! 0.0946
  else
     fMrTau_Cp = min(1.0, fMrTau_Cp_a1*exp(1.3*fClay)) ! 0.3646 unitless
     fMkTau_Cp = min(1.0, fMrTau_Cp_a2*exp(0.8*fClay)) ! 0.2255

     fMrTau_Cc = min(1.0-fMrTau_Cp, fMrTau_Cc_a4*exp(fMrTau_Cc_a3*fI_Lm))  ! 0.0315
     fMkTau_Cc = min(1.0-fMkTau_Cp, fMrTau_Cc_a5*exp(fMrTau_Cc_a3*fI_Lm))  ! 0.0946
  endif

  !print *, fMrTau_Cp,fMkTau_Cp,fMrTau_Cc,fMkTau_Cc

  pool%metabolicLitterC  = pool%metabolicLitterC  - (pool%DecompMrLm+pool%DecompMkLm)*dt_fast_hr ! kgC/m3
  pool%structuralLitterC = pool%structuralLitterC - (pool%DecompMrLs+pool%DecompMkLs)*dt_fast_hr


  pool%chemResistantC    = pool%chemResistantC    - (pool%OxidMrCc+pool%OxidMkCc)*dt_fast_hr     &
                                                  + (fMrTau_Cc*pool%MrTau + fMkTau_Cc*pool%MkTau)*dt_fast_hr ! kgC/m3

  if (is_sfc_litter) then
     pool%Desorb     = 0.0
     pool%protectedC = 0.0
  else
     pool%Desorb = (Desorb_kd * 1.5e-5*exp(Desorb_clay*fClay) * exp(Desorb_kdp * pool%protectedC) )*pool%protectedC ! kgC/m3/h
     pool%protectedC = pool%protectedC + (fMrTau_Cp*pool%MrTau + fMkTau_Cp*pool%MkTau - pool%Desorb)*dt_fast_hr ! kgC/m3
  endif


  if (DOC_cycling) then
     pool%Resp = (1-max(0.01,eLm_Mr-e_slope*T))*pool%DecompMrLm*(1-w_Lm) + (1-max(0.01,eLs_Mr-e_slope*T))*pool%DecompMrLs*(1-w_Ls) + (1-max(0.01,eCa_Mr-e_slope*T))*pool%DecompMrCa*(1-w_Ca) +(1-max(0.01,eCa_Mr-e_slope*T))*pool%DecompMrDOC + & ! kgC/m3/h
                 (1-max(0.01,eLm_Mk-e_slope*T))*pool%DecompMkLm*(1-w_Lm) + (1-max(0.01,eLs_Mk-e_slope*T))*pool%DecompMkLs*(1-w_Ls) + (1-max(0.01,eCa_Mk-e_slope*T))*pool%DecompMkCa*(1-w_Ca) +(1-max(0.01,eCa_Mk-e_slope*T))*pool%DecompMkDOC

     pool%availableC     = pool%availableC        - (pool%DecompMrCa+pool%DecompMkCa)*dt_fast_hr &
                                                  + (pool%OxidMrCc+pool%OxidMkCc)*dt_fast_hr     &
                                                  + ((1-fMrTau_Cp-fMrTau_Cc)*(1-fMrTau_DOC)*pool%MrTau + (1-fMkTau_Cp-fMkTau_Cc)*(1-fMkTau_DOC)*pool%MkTau)*dt_fast_hr ! kgC/m3

     pool%DOC            = pool%DOC               - (pool%DecompMrDOC+pool%DecompMkDOC)*dt_fast_hr &
                                                  + pool%Desorb*dt_fast_hr &
                                                  + ((1-fMrTau_Cp-fMrTau_Cc)*(fMrTau_DOC)*pool%MrTau + (1-fMkTau_Cp-fMkTau_Cc)*(fMkTau_DOC)*pool%MkTau)*dt_fast_hr &
                                                  + ((pool%DecompMrLm+pool%DecompMkLm)*w_Lm + (pool%DecompMrLs+pool%DecompMkLs)*w_Ls + (pool%DecompMrCa+pool%DecompMkCa)*w_Ca)*dt_fast_hr

     pool%microbesR = pool%microbesR + (max(0.01,eLm_Mr-e_slope*T)*pool%DecompMrLm*(1-w_Lm) + max(0.01,eLs_Mr-e_slope*T)*pool%DecompMrLs*(1-w_Ls) + max(0.01,eCa_Mr-e_slope*T)*pool%DecompMrCa*(1-w_Ca) + max(0.01,eCa_Mr-e_slope*T)*pool%DecompMrDOC - pool%MrTau)*dt_fast_hr ! kgC/m3
     pool%microbesK = pool%microbesK + (max(0.01,eLm_Mk-e_slope*T)*pool%DecompMkLm*(1-w_Lm) + max(0.01,eLs_Mk-e_slope*T)*pool%DecompMkLs*(1-w_Ls) + max(0.01,eCa_Mk-e_slope*T)*pool%DecompMkCa*(1-w_Ca) + max(0.01,eCa_Mk-e_slope*T)*pool%DecompMkDOC - pool%MkTau)*dt_fast_hr
  else
     pool%Resp = (1-max(0.01,eLm_Mr-e_slope*T))*pool%DecompMrLm + (1-max(0.01,eLs_Mr-e_slope*T))*pool%DecompMrLs + (1-max(0.01,eCa_Mr-e_slope*T))*pool%DecompMrCa + & ! kgC/m3/h
                 (1-max(0.01,eLm_Mk-e_slope*T))*pool%DecompMkLm + (1-max(0.01,eLs_Mk-e_slope*T))*pool%DecompMkLs + (1-max(0.01,eCa_Mk-e_slope*T))*pool%DecompMkCa

     pool%availableC     = pool%availableC        - (pool%DecompMrCa+pool%DecompMkCa)*dt_fast_hr &
                                                  + (pool%OxidMrCc+pool%OxidMkCc)*dt_fast_hr     &
                                                  + ((1-fMrTau_Cp-fMrTau_Cc)*pool%MrTau + (1-fMkTau_Cp-fMkTau_Cc)*pool%MkTau)*dt_fast_hr & ! kgC/m3
                                                  + pool%Desorb*dt_fast_hr
     pool%DOC               = 0.0

     pool%microbesR = pool%microbesR + (max(0.01,eLm_Mr-e_slope*T)*pool%DecompMrLm + max(0.01,eLs_Mr-e_slope*T)*pool%DecompMrLs + max(0.01,eCa_Mr-e_slope*T)*pool%DecompMrCa - pool%MrTau)*dt_fast_hr ! kgC/m3
     pool%microbesK = pool%microbesK + (max(0.01,eLm_Mk-e_slope*T)*pool%DecompMkLm + max(0.01,eLs_Mk-e_slope*T)*pool%DecompMkLs + max(0.01,eCa_Mk-e_slope*T)*pool%DecompMkCa - pool%MkTau)*dt_fast_hr

  endif

end subroutine update_GIMICS_pool

! ============================================================================
! note that in resp_denitrif the dependence on soil moisture is subtly different
real function theta_func ( water_filled_porosity, air_filled_porosity, &
                           substrate_diffusion_exp, gas_diffusion_exp, &
                           min_anaerobic_resp_factor, min_dry_resp_factor )
  real, intent(in) :: water_filled_porosity ! fraction of pores filled with water
  real, intent(in) :: air_filled_porosity ! fraction of pores filled with water
  real, intent(in) :: substrate_diffusion_exp,gas_diffusion_exp,min_dry_resp_factor,min_anaerobic_resp_factor

  real :: theta_resp_max, aerobic_max

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
real function theta_func_orchidee (moist,theta_func_orchidee_min,theta_func_orchidee_max)

  real, intent(in) :: moist ! soil moisture, fraction of soil filled with water
  real, intent(in) :: theta_func_orchidee_min, theta_func_orchidee_max

  theta_func_orchidee=max(theta_func_orchidee_min,min(theta_func_orchidee_max,(1.1*(moist**2))+(2.4*moist)+0.29))

end function theta_func_orchidee


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
  real    :: deltaBulk, deltaRhiz

!   NH4(:)=0.0
!   NO3(:)=0.0
!   if(present(ammonium)) NH4=ammonium
!   if(present(nitrate))  NO3=nitrate

  do k=1,num_l
     if (soilc%fRhiz(k)>1.0e-10) then
         deltaRhiz = exudateC(k)/(dz(k)*soilc%fRhiz(k)) ! kgC/m3 of rhizosphere
         deltaBulk = 0.0
         ! slm: should we add some protection from very small rhizosphere fractions that may
         !      result in huge per-volume input?
     else ! rhizosphere does not exist, add exudates to bulk
         deltaRhiz = 0.0
         deltaBulk = exudateC(k)/dz(k) ! kgC/m3 of bulk soil == kgC/m3 of all soil, since there is no rhizosphere
     endif
     soilC%rhiz(k)%metabolicLitterC = soilc%rhiz(k)%metabolicLitterC + deltaRhiz
     soilC%bulk(k)%metabolicLitterC = soilc%bulk(k)%metabolicLitterC + deltaBulk

     ! save exudate inputs for diagnostics
     soilC%rhiz(k)%InputExdC = soilc%rhiz(k)%InputExdC + deltaRhiz/dt_fast_hr
     soilC%bulk(k)%InputExdC = soilc%bulk(k)%InputExdC + deltaBulk/dt_fast_hr
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

  if (is_watch_point()) then
     write(*,'(a25)',advance='NO') 'add_soil_matter_GIMICS:'
     __DEBUG2__(leaf_litt_C, wood_litt_C)
  endif
! + sanity check
  do k = 1,N_LITTER_POOLS
     call check_GIMICS_pool(soilc%litt(k), trim(l_diagname(k))//'litt before add matter')
  enddo
  do k=1,num_l
     call check_GIMICS_pool(soilc%rhiz(k), 'rhiz('//string(k)//') before add matter')
     call check_GIMICS_pool(soilc%bulk(k), 'bulk('//string(k)//') before add matter')
  enddo
! - sanity check

  ! sometimes litter input becomes negative (because vegn carbon_gain could be negative)
  ! To handle this situation an conserve carbon, we store negatives in a pool, and then
  ! fill it up later with positive inputs
  call borrow_to_negatives(wood_litt_C,soilc%neg_litt_C) ! borrow from wood litter first
  call borrow_to_negatives(leaf_litt_C,soilc%neg_litt_C) ! and from leaf litter second

  call add_matter_GIMICS1(soilc%litt(LITT_LEAF),  leaf_litt_C, leaf_litt_N)
  call add_matter_GIMICS1(soilc%litt(LITT_CWOOD), wood_litt_C, wood_litt_N)
  do k = 1,size(soilc%bulk)
     call add_matter_GIMICS2(soilc%bulk(k), soilc%rhiz(k), dz(k), soilc%fRhiz(k), root_litt_C(k,:), root_litt_N(k,:))
  enddo

! + sanity check
  do k = 1,N_LITTER_POOLS
     call check_GIMICS_pool(soilc%litt(k), trim(l_diagname(k))//'litt after add matter')
  enddo
  do k=1,num_l
     call check_GIMICS_pool(soilc%rhiz(k), 'rhiz('//string(k)//') af add matter')
     call check_GIMICS_pool(soilc%bulk(k), 'bulk('//string(k)//') before add matter')
  enddo
! - sanity check

  ! accumulate litterfall diagnostics: it is sent to diag and then reset at every time step
  vegn%litterfall_C(:,LITT_LEAF)  = vegn%litterfall_C(:,LITT_LEAF)  + leaf_litt_C(:)
  vegn%litterfall_C(:,LITT_CWOOD) = vegn%litterfall_C(:,LITT_CWOOD) + wood_litt_C(:)

contains
  ! given litter and amount of negative litter from previous time step, attempts to borrow
  ! positive carbon to reduce the amount of negativs
  subroutine borrow_to_negatives(litt, negatives)
    real, intent(inout) :: litt(:), negatives(:)

    litt      = litt + negatives
    negatives = min(litt,0.0)
    litt      = max(litt,0.0)
  end subroutine borrow_to_negatives

end subroutine add_soil_matter_GIMICS

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
     call add_matter_GIMICS1(soilc%litt(i), C=delta_C(:,i))
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

! ============================================================================
! add carbon (and later nitrogen) to GIMICS soil BGC pool
subroutine add_matter_GIMICS1(pool, C, N)
  type(GIMICS_BGC_litt), intent(inout) :: pool ! BGC pool to update
  real, intent(in), optional :: C (N_C_TYPES)  ! (fast,slow,[dead]microbial), kgC/m2
  real, intent(in), optional :: N (N_C_TYPES)  ! (fast,slow,[dead]microbial), kgN/m2

  real :: deltaMtb, deltaStr
  if (present(C)) then
     deltaMtb = (C(C_FAST) + C(C_MIC))/pool%dz ! kgC/m3
     deltaStr = C(C_SLOW)/pool%dz ! kgC/m3

     pool%metabolicLitterC  = pool%metabolicLitterC  + deltaMtb
     pool%structuralLitterC = pool%structuralLitterC + deltaStr
     ! update the deposition accumulators, for diagnostics
     pool%InputMtbC   = pool%InputMtbC + deltaMtb/dt_fast_hr
     pool%InputStrC   = pool%InputStrC + deltaStr/dt_fast_hr
  endif
!   if (present(N)) then
!      ....
!   endif
  call update_litter_thickness(pool)
end subroutine

! ============================================================================
subroutine debug_pool(pool, tag)
  class(GIMICS_BGC_pool), intent(in) :: pool
  character(*),           intent(in) :: tag

  write (*,'(a16,":")',advance='NO') trim(tag)
  call dpri('metabolicC',pool%metabolicLitterC)
  call dpri('structuralC',pool%structuralLitterC)
  call dpri('protectedC',pool%protectedC)
  call dpri('chemResistC',pool%chemResistantC)
  call dpri('availableC',pool%availableC)
  call dpri('microbesR',pool%microbesR)
  call dpri('microbesK',pool%microbesK)
  call dpri('DOC',pool%DOC)
  select type(pool)
  type is (GIMICS_BGC_litt)
     call dpri('dz',pool%dz)
  type is (GIMICS_BGC_pool)
     ! do nothing
  end select
  write(*,*)
end subroutine debug_pool

! ============================================================================
! add carbon (and later nitrogen) to GIMICS soil BGC pool, distributing it between bulk soil and rhizosphere
subroutine add_matter_GIMICS2(bulk, rhiz, dz, rhiz_frac, C, N)
  type(GIMICS_BGC_pool), intent(inout) :: bulk, rhiz ! bulk soil and rhizosphere BGC pools, respectively
  real, intent(in) :: dz        ! layer thickness, m
  real, intent(in) :: rhiz_frac ! fraction of soil in rhizosphere, currently unused because
                                ! we add the same concentrations to bulk soil and rhizosphere
  real, intent(in), optional :: C (N_C_TYPES)  ! (fast,slow,[dead]microbial), kgC/m2
  real, intent(in), optional :: N (N_C_TYPES)  ! (fast,slow,[dead]microbial), kgN/m2

  real :: deltaMtb, deltaStr ! increments of metabolic and structural C, respectively

  if (present(C)) then
     deltaMtb = (C(C_FAST)+ C(C_MIC))/dz ! kgC/m3
     deltaStr = C(C_SLOW)/dz             ! kgC/m3

     ! bulk soil
     bulk%metabolicLitterC  = bulk%metabolicLitterC  + deltaMtb
     bulk%structuralLitterC = bulk%structuralLitterC + deltaStr
     ! update the deposition accumulators, for diagnostics
     bulk%InputMtbC   = bulk%InputMtbC + deltaMtb/dt_fast_hr
     bulk%InputStrC   = bulk%InputStrC + deltaStr/dt_fast_hr

     ! rhizosphere
     rhiz%metabolicLitterC  = rhiz%metabolicLitterC  + deltaMtb
     rhiz%structuralLitterC = rhiz%structuralLitterC + deltaStr
     ! update the deposition accumulators, for diagnostics
     rhiz%InputMtbC  = rhiz%InputMtbC + deltaMtb/dt_fast_hr
     rhiz%InputStrC  = rhiz%InputStrC + deltaStr/dt_fast_hr
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

  integer :: k
  real :: f

  burned_C = 0.0; burned_N = 0.0
  do k = 1,N_LITTER_POOLS
     associate (pool=>soilc%litt(k))
     burned_C = burned_C + C_amount(pool)*frac(k)
     f = 1.0-frac(k)
     pool%metabolicLitterC  = f * pool%metabolicLitterC
     pool%structuralLitterC = f * pool%structuralLitterC
     pool%protectedC        = f * pool%protectedC
     pool%chemResistantC    = f * pool%chemResistantC
     pool%availableC        = f * pool%availableC
     pool%microbesR         = f * pool%microbesR
     pool%microbesK         = f * pool%microbesK
     pool%DOC               = f * pool%DOC

     call update_litter_thickness(pool)
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
  real, intent(in) :: flow(:), div(:), wl(:) ! flow (into layer) and wl in units of mm, downward is >0  !!!xz check the unit of dz (should be m in this subroutine), flow (should be mm)
  real, intent(in) :: div_hlsp_DOC(:,:) ! dim(N_C_TYPES, num_l) [kg C/m^2/s] net divergence loss from tile calculated in hlsp_hydrology
  real, intent(in) :: div_hlsp_DON(:,:) ! dim(N_C_TYPES, num_l) [kg N/m^2/s] net divergence
  real, intent(in) :: div_hlsp_NO3(:)   ! dim(num_l) [kg N/m^2/s] net divergence loss from tile calculated in hlsp_hydrology
  real, intent(in) :: div_hlsp_NH4(:)   ! dim(num_l) [kg N/m^2/s] net divergence loss from tile calculated in hlsp_hydrology

  real, intent(out) :: total_DOC_div, total_DON_div, total_NO3_div, total_NH4_div

  ! ---- local vars
  real :: surf_DOC_loss ! [kg C/m^2] loss from top layer to surface runoff loss from tile calculated in hlsp_hydrology
  real :: DOC0(0:num_l)  ! [kg C/m^2] initial amount of DOC per layer; layer 0 is surface litter
  real :: litt_DOC0(N_LITTER_POOLS) ! [kg C/m3] initial concentration of DOC per litter pool, for diagnostics of tendencies
  real :: DOC(0:num_l)  ! [kg C/m^2] amount of DOC per layer
  real :: DOC1(0:num_l)  ! [kg C/m^2] amount of DOC per layer
  real :: d_DOC(0:num_l)
  real :: div_loss(0:num_l)!xz

  ! soil flow-related variables expanded to include litter layer
  real :: flow_with_litter(0:num_l)
  real :: div_with_litter (0:num_l)
  real :: dz_with_litter  (0:num_l)

  real :: dz0 ! total surface litter thickness for concentration update in surface litter
  real :: scale ! scaling factor for concentrations update in soil layer

  real :: mass0, mass1 ! for mass conservation checks
  real :: carbon0, carbon1

  integer :: k

  real, parameter :: minwl    = 0.1 ! [mm]

  ! for C conservation checks across this entire subroutine
  carbon0 = soilc%total_C()

  total_DOC_div = 0.0
  total_DON_div = 0.0
  total_NO3_div = 0.0
  total_NH4_div = 0.0

!!!!!!!xz note: please make sure the unit of flow is ????
  flow_with_litter(0)=0.0
  flow_with_litter(1:num_l)=flow(1:num_l)  !mm
  flow_with_litter = flow_with_litter/1000 !xz change the div unit from mm to m

  div_with_litter(0)=0.0
  div_with_litter(1:num_l)=div(:)*delta_time ! div is in mm/s
  div_with_litter=div_with_litter/1000 !xz change the div unit from mm to m

  dz_with_litter(0)=sum(soilc%litt(:)%dz) ! slm: total litter thickness is the sum of all
                                          ! litter thicknesses. Is this reasonable?
  dz_with_litter(1:num_l) = dz(1:num_l) !!xz assume the unit of dz is m

  surf_DOC_loss = 0.0

  ! calculate amount of DOC in surface litter, kgC/m2
  DOC(0) = 0.0
  do k = 1, size(soilc%litt)
     DOC(0) = DOC(0) + soilc%litt(k)%DOC * soilc%litt(k)%dz
  enddo
  ! calculate amount of DOC in soil layers, kgC/m2
  do k = 1, num_l
     DOC(k)=(soilc%rhiz(k)%DOC*soilc%fRhiz(k) + soilc%bulk(k)%DOC*(1-soilc%fRhiz(k)))*dz(k)
  enddo
  ! store initial values of DOC, for scaling after update
  DOC0(:) = DOC(:)

  ! advect DOC with water flow
  call check_var_range(DOC(:), 0.0, HUGE(1.0), 'tracer_leaching_GIMICS', 'DOC before advection',FATAL)
  mass0 = sum(DOC(:))
  call tracer_advection(DOC(:),flow_with_litter(:),div_with_litter(:),dz_with_litter,d_DOC(:),div_loss(:),wl(:))!xz
  mass1 = sum(DOC(:))
  call check_conservation('tracer_leaching_GIMICS','DOC mass in advection', mass0,mass1, carbon_cons_tol )
  call check_var_range(DOC(:), 0.0, HUGE(1.0), 'tracer_leaching_GIMICS', 'DOC after advection',FATAL)

  if (gw_option == GW_TILED) then ! reset div_loss(1:num_l) according to values calculated in hlsp_hydrology
     div_loss(1:num_l) = div_hlsp_DOC(1,:)*delta_time ! slm: using only the first carbon type in GIMICS
     if (flow(1) < 0 .and. wl(1) > minwl) then  ! Add loss from top layer to runoff
        surf_DOC_loss = -DOC(1) * flow(1) / wl(1)
        surf_DOC_loss = min(surf_DOC_loss, DOC(1))
     end if
     div_loss(1) = min(div_loss(1), DOC(1) - surf_DOC_loss)
     do k=2,num_l
        div_loss(k) = min(div_loss(k), DOC(k))
     end do
     ! Note: if these limits are imposed, there will be an imbalance between inter-tile fluxes
     ! that will be effectively rectified by subtracting from the flux to stream. In rare
     ! situations, that could lead to a negative stream DOC flux.
  end if
  DOC(:)=DOC(:)-div_loss(:)
  ! slm: I do not understand why the calculations below operate on first soil layer,
  ! rather than surface litter.
  DOC(1)=DOC(1)-surf_DOC_loss !!xz This line does not exist in CH's code ; consider to add similar line to Nitrogen part
  ! Xin says this line was a mistake

  ! sum up the total loss due to water flow divergence
  total_DOC_div = surf_DOC_loss + sum(div_loss(1:))

  ! update DOC concentrations in the surface litter
  litt_DOC0(:) = soilc%litt(:)%DOC ! save old values by litter type, for tendency diagnostics
  if(DOC0(0)>0) then
     ! there was some DOC in surface litter initially: scale concentrations in pools
     ! proportionally
     soilc%litt(:)%DOC = soilc%litt(:)%DOC * DOC(0)/DOC0(0)
  else
     ! there were no DOC in surface litter: assign the same DOC concentration to all
     ! surface litter pools
     dz0 = sum(soilc%litt(:)%dz)
     if (dz0 > 0) then
        ! there is some surface litter: assign concentration
        soilc%litt(:)%DOC = DOC(0)/dz0
     else
        ! No surface litter: put surface DOC into top soil layer.
        ! An alternative would be to add it to surface loss (?)
        DOC(1) = DOC(1) + DOC(0); DOC(0) = 0.0
        soilc%litt(:)%DOC = 0.0
     endif
  endif

  ! update DOC concentrations in soil layers
  ! slm: this treatment distributes resulting DOC between bulk and rhizosphere
  ! proportionally to the values before advection (that is, scales them up or down).
  ! It seems reasonable at the first glance, but does this treatment describe what
  ! happens in advection correctly?
  !   There should be some exchange between rhizosphere and bulk soil due to water flow,
  ! in addition to diffusion that should also exist.
  !   Also, numerically, what if the initial value of DOC is positive but very small?
  ! Is it possible we get crazy numbers in this case?
  do k = 1,num_l
     if (DOC0(k)>0) then
        scale = DOC(k)/DOC0(k)
        soilc%rhiz(k)%DOC = scale * soilc%rhiz(k)%DOC
        soilc%bulk(k)%DOC = scale * soilc%bulk(k)%DOC
     else
        soilc%rhiz(k)%DOC = DOC(k)/dz(k)
        soilc%bulk(k)%DOC = DOC(k)/dz(k)
     endif
  enddo

  if (is_watch_point()) then
     write(*,*)'#### tracer_leaching_GIMICS'
     __DEBUG4__(total_DOC_div,surf_DOC_loss,div_loss(0),sum(div_loss(1:)))
     __DEBUG1__(soilc%total_C())
     ! calculate total DOC in surface litter, kgC/m2
     DOC1(0) = 0.0
     do k = 1, size(soilc%litt)
        DOC1(0) = DOC1(0) + soilc%litt(k)%DOC * soilc%litt(k)%dz
     enddo
     ! calculate total DOC in soil layers, kgC/m2
     do k=1,num_l
        DOC1(k)=(soilc%rhiz(k)%DOC*soilc%fRhiz(k) + soilc%bulk(k)%DOC*(1-soilc%fRhiz(k)))*dz(k)
     end do
     do k=0,num_l
        write(*,'(i2.2)', advance='NO') k
        call dpri('dz',dz_with_litter(k))
        call dpri('DOC0',DOC0(k))
        call dpri('DOC',DOC(k))
        call dpri('DOC1',DOC1(k))
        call dpri('diff',DOC1(k)-DOC(k))
        call dpri('d_DOC',d_DOC(k))
        call dpri('div_loss',div_loss(k))
        call dpri('flow',flow_with_litter(k))
        if (k>0) then
           call dpri('fRhiz',soilc%fRhiz(k))
        endif
        write(*,*)
     enddo
     __DEBUG3__(sum(DOC0),sum(DOC),sum(DOC)-sum(DOC0))
     __DEBUG2__(sum(DOC1),sum(DOC)-sum(DOC1))
  endif

  carbon1 = soilc%total_C() + total_DOC_div
  call check_conservation('tracer_leaching_GIMICS','carbon', carbon0, carbon1, carbon_cons_tol )

  ! diagnostics
  call send_tile_data(id_total_DOC_div_loss,total_DOC_div/dt_fast_yr, diag)
  ! it appears that surf_DOC_loss will be zero except in GW_TILED (hydroblocks) hydrology
  call send_tile_data(id_surf_DOC_loss, surf_DOC_loss/dt_fast_yr, diag)
  do k = 1, N_LITTER_POOLS
     call send_tile_data(id_ladvec_DOC(k), (soilc%litt(k)%DOC-litt_DOC0(k))/dt_fast_yr, diag)
  enddo
  if (id_sadvec_DOC>0) then
     call send_tile_data(id_sadvec_DOC, (DOC(1:num_l)-DOC0(1:num_l))/(dz(:)*dt_fast_yr), diag)
  endif
end subroutine tracer_leaching_GIMICS

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
! slm: it is the same as in CORPSE. Should we move it to one of the vegetation modules,
! or keep it here to be able to make changes independent from CORPSE implementation?

!> @brief Calculate volumetric fraction of rhizosphere in each layer
subroutine rhizosphere_frac(vegn, rFrac)
  type(vegn_tile_type), intent(in)  :: vegn !< vegetation state
  real                , intent(out) :: rFrac(:)!< volumetric fraction of rhizosphere

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
