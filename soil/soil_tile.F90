#include <fms_platform.h>

module soil_tile_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : file_exist, check_nml_error, &
     close_file, stdlog, read_data, error_mesg, FATAL
use constants_mod, only : &
     pi, tfreeze, rvgas, grav, dens_h2o, hlf, epsln
use land_constants_mod, only : NBANDS
use land_data_mod, only : log_version
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_SOIL, register_tile_selector
use land_io_mod, only : print_netcdf_error
use soil_carbon_mod, only : soil_carbon_option, SOILC_CORPSE, SOILC_CORPSE_N, &
    soil_pool, combine_pools, init_soil_pool, poolTotals, N_C_TYPES

implicit none
private

! ==== netcdf declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

! ==== public interfaces =====================================================
public :: soil_pars_type
public :: soil_tile_type

public :: new_soil_tile, delete_soil_tile
public :: soil_tiles_can_be_merged, merge_soil_tiles
public :: soil_is_selected
public :: get_soil_tile_tag
public :: soil_tile_stock_pe
public :: soil_tile_heat, soil_tile_carbon, soil_tile_nitrogen

public :: read_soil_data_namelist

public :: soil_ice_porosity
public :: soil_ave_ice_porosity

public :: soil_data_radiation
public :: soil_roughness
public :: soil_data_thermodynamics
public :: soil_data_hydraulic_properties
public :: soil_data_psi_for_rh
public :: soil_data_gw_hydraulics
public :: soil_data_gw_hydraulics_ar5
public :: soil_data_vwc_for_init_only
public :: soil_data_init_derive_subsurf_pars
public :: soil_data_init_derive_subsurf_pars_ar5
public :: soil_data_init_derive_subsurf_pars_tiled
public :: soil_ave_temp  ! calculate average soil temeperature
public :: soil_ave_theta0! calculate average soil moisture, pcm based on available water, zeta input
public :: soil_ave_theta1! calculate average soil moisture, ens based on all water
public :: soil_ave_wetness ! calculate average soil wetness
public :: soil_theta     ! returns array of soil moisture, for all layers
public :: soil_psi_stress ! return soil-water-stress index

! public data
public :: max_lev ! max number of soil layers (max dimension of arrays)
public :: num_l ! actual number of soil layers
public :: dz    ! layer thicknesses (m)
public :: zhalf ! depths of layer boundaries (m)
public :: zfull ! depths of layer centers (m)

public :: g_iso, g_vol, g_geo, g_RT
public :: num_storage_pts
public :: gw_zeta_s, gw_flux_table, gw_area_table
public :: gw_scale_length, gw_scale_relief, gw_scale_soil_depth, slope_exp
public :: num_zeta_pts, num_tau_pts
public :: log_tau, log_zeta_s, log_rho_table, gw_scale_perm, aspect
public :: use_alpha, z_ref, k0_macro_z, k0_macro_x, use_tau_fix
public :: retro_a0n1

public :: psi_wilt ! wilting water potential, m
! =====end of public interfaces ==============================================
interface new_soil_tile
   module procedure soil_tile_ctor
   module procedure soil_tile_copy_ctor
end interface

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'soil_tile_mod'
#include "../shared/version_variable.inc"

integer, parameter :: max_lev          = 100
integer, parameter, public :: n_dim_soil_types = 14      ! max size of lookup table
real,    parameter :: small            = 1.e-4
real,    parameter :: t_ref            = 293
real,    parameter :: g_RT             = grav / (rvgas*t_ref)
real,    parameter :: sigma_max        = 2.2
real,    parameter :: K_rel_min        = 1.e-12
real,    parameter, public :: initval  = 1.e36 ! For initializing variables

! from the modis brdf/albedo product user's guide:
real, parameter :: g_iso  = 1.
real, parameter :: g_vol  = 0.189184
real, parameter :: g_geo  = -1.377622
real, parameter :: g0_iso = 1.0
real, parameter :: g1_iso = 0.0
real, parameter :: g2_iso = 0.0
real, parameter :: g0_vol = -0.007574
real, parameter :: g1_vol = -0.070987
real, parameter :: g2_vol =  0.307588
real, parameter :: g0_geo = -1.284909
real, parameter :: g1_geo = -0.166314
real, parameter :: g2_geo =  0.041840

! geohydrology option selector (not used by lm2)
integer, parameter, public ::   &
     GW_LM2            = 0, &
     GW_LINEAR         = 1, &
     GW_HILL_AR5       = 2, &
     GW_HILL           = 3, &
     GW_TILED          = 4

! litter pool constants
integer, parameter, public :: &
     N_LITTER_POOLS    = 3, &
     LEAF              = 1, & ! leaf litter
     CWOOD             = 2, & ! coarse wood litter
     FWOOD             = 3    ! fine wood litter

character(16), parameter, public :: l_shortname(N_LITTER_POOLS) = [ &
     'leaf            ', &
     'coarsewood      ', &
     'finewood        '  ]

character(16), parameter, public :: l_longname(N_LITTER_POOLS) = [ &
     'leaf            ', &
     'coarse wood     ', &
     'fine wood       '  ]

! ==== types =================================================================
type :: soil_pars_type
  real vwc_wilt
  real vwc_fc
  real vwc_sat
  real vlc_min
  real awc_lm2
  real k_sat_ref
  real psi_sat_ref
  real chb
  ! Many of these thermal and hydraulic "parameters" will become prognostic
  ! with the peat model.
  real heat_capacity_dry
  real thermal_cond_dry
  real thermal_cond_sat
  real thermal_cond_exp
  real thermal_cond_scale
  real thermal_cond_weight
  real refl_dry_dir(NBANDS)
  real refl_dry_dif(NBANDS)
  real refl_sat_dir(NBANDS)
  real refl_sat_dif(NBANDS)
  real f_iso_dry(NBANDS)
  real f_vol_dry(NBANDS)
  real f_geo_dry(NBANDS)
  real f_iso_sat(NBANDS)
  real f_vol_sat(NBANDS)
  real f_geo_sat(NBANDS)
  real emis_dry
  real emis_sat
  real z0_momentum
  real tfreeze
  ! Note freezing point depression that differs among layers or tiles requires additional treatment
  ! for energy conservation.

  ! Used primarily with LM2-3.1 soil hydrology:
  real tau_groundwater     ! groundwater residence time (1/s)
  real hillslope_length    ! hillslope length for integrated-tile model (m)
  real hillslope_relief    ! hillslope max elevation for integrated-tile model (m)
  real hillslope_a         ! parameter of convergence/divergence of the hilslopes
  real hillslope_n         ! exponent in z(x) relationship for hillslopes
  real hillslope_zeta_bar  ! hillslope area-normalized elevation
  real zeta                ! soil depth normalized by hillslope elevation
  real tau                 ! effective residence time (1/s)
  integer storage_index
  real rsa_exp          ! riparian source-area exponent
  ! Used with full tiled Hillslope Model Hydrology:
  real k_sat_gw         ! saturated hydraulic conductivity at infinite depth, for bedrock (mm/s)
  real k_sat_sfc        ! saturated hydraulic conductivity at surface (mm/s)
  real soil_e_depth     ! soil depth scale for taking on properties of bedrock (m)
  real microtopo        ! microtopographic roughness length for surface water (m)
  real tile_hlsp_length ! horizontal length of tile parallel to hillslope (m) (proportional to
                        ! tile area)
  real tile_hlsp_slope  ! vertical slope of tile (-)
  real tile_hlsp_elev   ! elevation of center of tile above streambed at hillslope bottom (m)
  real tile_hlsp_hpos   ! horizontal position of tile center along hillslope (m)
  real tile_hlsp_width  ! width of tile perpendicular to hillslope, normalized to strm width (-)
                        ! (proportional to tile area)
!  real transm_bedrock    ! bedrock / inf-depth-wat-tab transmissivity (m^2/s, or vol wat/m/s)
  real disturb_scale    ! characteristic horizontal disturbance lengthscale within hillslope (m)
                        ! This will need to be set in transitions; otherwise, defaults to 1/2
                        ! tile_hlsp_width.

  real Qmax             ! Maximum carbon sorption capacity (kgC/m3 soil)
end type soil_pars_type


type :: soil_tile_type
   integer :: tag ! kind of the soil
   ! For hillslope model. Note, these indices establish behavior of tiles and can change
   ! during run so must be saved with restart.
   integer hidx_j         ! id of vertical tile position within hillslope
   integer hidx_k         ! id of hillslope parent topo type
       ! (Note multiple instances of each of these can exist, via, e.g. land use
       ! disturbance. So these indices function similarly to "tag".)

   type(soil_pars_type) :: pars

   real, allocatable ::  &
       wl(:)           , & ! liquid water, kg/m2
       ws(:)           , & ! solid water, kg/m2
       T(:)            , & ! temperature, degK
       groundwater(:)  , &
       groundwater_T(:), &
       w_fc(:)         , &
       w_wilt(:)       , &
       d_trans(:)      , &
       alpha(:)        , &
       k_macro_z(:)    , & ! Vertical macroporosity [mm/s]
       k_macro_x(:)    , & ! Horizontal macroporosity [mm/s]
       vwc_max(:)
   real :: Eg_part_ref
   real :: z0_scalar
   real :: geothermal_heat_flux
   ! data that were local to soil.f90
   real, allocatable ::  &
       heat_capacity_dry(:), &
       e(:), f(:),       &
       gw_flux_norm(:),  &
       gw_area_norm(:)
   ! added to avoid recalculation of soil hydraulics in case of Darcy uptake
   real          :: uptake_T
   ! These two are needed for tiled hillslope hydrology
   real, allocatable :: psi(:) ! soil water potential [m]
   real, allocatable :: hyd_cond_horz(:) ! soil hydraulic conductivity for inter-tile transfers [mm/s]
   ! flux variables for tiled hillslope hydrology
   real, allocatable :: div_hlsp(:) ! net groundwater divergence flux from tile to hillslope
                                     ! or stream [mm/s]
   real, allocatable :: div_hlsp_heat(:) ! net heat divergence flux associated with groundwater
                                     ! (relative to tfreeze) [W/m^2]

   ! soil carbon
   ! CENTURY-style values
   real, allocatable :: &
       fast_soil_C(:), & ! fast soil carbon pool, (kg C/m2), per layer
       slow_soil_C(:)    ! slow soil carbon pool, (kg C/m2), per layer
   ! values for CORPSE
   type(soil_pool) :: litter(N_LITTER_POOLS) ! Surface litter pools, just one layer
   type(soil_pool), allocatable :: soil_organic_matter(:) ! Soil carbon in soil layers, using soil_carbon_mod soil carbon pool type
   integer, allocatable   :: is_peat(:)             ! Keeps track of whether soil layer is peat, for redistribution
   real                   :: fast_DOC_leached !Carbon that has been leached out of the column
   real                   :: slow_DOC_leached !Carbon that has been leached out of the column
   real                   :: deadmic_DOC_leached !Carbon that has been leached out of the column
   real                   :: passive_N_uptake = 0.0 ! N uptake by water flux into roots
   real                   :: myc_scav_N_uptake = 0.0 ! N uptake by "scavenger" mycorrhizae (mostly corresponding to Arbuscular mycorrhizae)
   real                   :: myc_mine_N_uptake = 0.0 ! N uptake by "miner" mycorrhizae (corresponding to Ectomycorrhizae)
   real                   :: symbiotic_N_fixation = 0.0 ! N fixation by symbiotic microbes
   real                   :: active_root_N_uptake = 0.0 ! Mineral N uptake from rhizosphere by active transport across root membrane
   ! values for the diagnostic of carbon budget and soil carbon acceleration
   real, allocatable :: &
       asoil_in(:), &
       fsc_in(:), &
       ssc_in(:)

   ! For storing DOC fluxes in tiled model
   real, allocatable :: div_hlsp_DOC(:,:) ! dimension (n_c_types, num_l) [kg C/m^2/s] net flux of carbon pools
                                      ! out of tile
   real, allocatable :: div_hlsp_DON(:,:) ! dimension (n_c_types, num_l) [kg N/m^2/s] net flux of nitrogen pools
                                     ! out of tile
   real, allocatable :: div_hlsp_NO3(:)  ! dimension (num_l) [kg N/m^2/s] net flux of nitrate out of tile
   real, allocatable :: div_hlsp_NH4(:)  ! dimension (num_l) [kg N/m^2/s] net flux of ammonium out of tile
end type soil_tile_type

! ==== module data ===========================================================
integer, public :: gw_option

real, public :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

!---- namelist ---------------------------------------------------------------
character(256), public :: soil_type_file = 'INPUT/ground_type.nc'
real    :: psi_wilt              = -150.0  ! matric head at wilting
real, public :: comp             = 0.001  ! m^-1, dThdPsi at saturation
real    :: K_min                 = 0.     ! absolute lower limit on hydraulic cond
                                          ! used only when use_alt[2]_soil_hydraulics
real    :: K_max_matrix          = 1.e10
real    :: DThDP_max             = 1.e10
real    :: psi_min               = -1.e5  ! value beyond which psi(theta) is
                                          ! linearly extrapolated
real    :: k_over_B              = 2         ! reset to 0 for MCM
real    :: rate_fc               = 0.1/86400 ! 0.1 mm/d drainage rate at FC
real    :: sfc_heat_factor       = 1
real    :: z_sfc_layer           = 0.0
real    :: sub_layer_tc_fac      = 1.0
real    :: z_sub_layer_min       = 0.0
real    :: z_sub_layer_max       = 0.0
real    :: freeze_factor         = 1.0
real    :: aspect                = 1.0
real    :: zeta_mult             = 1.0  ! multiplier for root depth scale in stress index
integer :: num_l                 = 18        ! number of soil levels
real    :: dz(max_lev)           = (/ &
    0.02, 0.04, 0.04, 0.05, 0.05, 0.1, 0.1, 0.2, 0.2, &
    0.2,   0.4,  0.4,  0.4,  0.4, 0.4,  1.,  1.,  1., &
    0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)
                                              ! thickness (m) of model layers,
                                              ! from top down
logical :: use_lm2_awc           = .false.
logical :: lm2                   = .false.
logical :: use_alt_psi_for_rh    = .false.
logical :: use_tau_fix           = .false.
logical :: use_sat_fix           = .false.
logical :: use_comp_for_ic       = .false.
logical :: use_comp_for_push     = .false.
logical :: limit_hi_psi          = .false.
logical :: use_fluid_ice         = .false.
logical :: use_alt3_soil_hydraulics = .false.
logical :: limit_DThDP           = .false.
logical :: retro_a0n1            = .false.
! ---- remainder are used only for cold start ---------
character(32), public :: soil_to_use     = 'single-tile'
       ! 'multi-tile' for multiple soil types per grid cell, a tile per type
       ! 'single-tile' for geographically varying soil with single type per
       !     model grid cell [default]
       ! 'uniform' for global constant soil, e.g., to reproduce MCM
       ! 'tiled_by_hlsp': input data provided where soil type is a function of the hillslope type
                          ! (requires use of the hillslope model)
character(32) :: geohydrology_to_use = 'hill_ar5'
       ! 'lm2' for lm2 hydrology -- must be consistent with lm2 flag
       ! 'linear'
       ! 'hill' -- for integrated-tile hillslope geohydrology
       ! 'hill_ar5' -- integrated-tile hillslope geohydrology reproducing AR5 simulations
       ! 'tiled' -- for full hillslope-model geohydrology with hillslopes discretized into tiles
logical :: use_single_geo        = .false.   ! .true. for global gw res time,
                                             ! e.g., to recover MCM
logical :: use_alpha             = .true.    ! for vertical change in soil properties
integer, public :: soil_index_constant = 9   ! index of global constant soil,
                                             ! used when use_single_soil
real    :: gw_res_time           = 60.*86400 ! mean groundwater residence time,
                                             ! used when use_single_geo
real    :: rsa_exp_global        = 1.5
real    :: gw_scale_length       = 1.0
real    :: gw_scale_relief       = 1.0
real    :: gw_scale_soil_depth   = 1.0
real    :: slope_exp             = 0.0
real    :: gw_scale_perm         = 1.0
real    :: k0_macro_z            = 0.0
real    :: k0_macro_x            = 0.0
real    :: log_rho_max           = 2.0
real    :: z_ref                 = 0.0       ! depth where [psi/k]_sat = [psi/k]_sat_ref
real    :: geothermal_heat_flux_constant = 0.0  ! true continental average is ~0.065 W/m2
real    :: Dpsi_min_const        = -1.e16

real, dimension(n_dim_soil_types) :: &
  dat_w_sat=&
  (/ 0.380, 0.445, 0.448, 0.412, 0.414, 0.446, 0.424, 0.445, 0.445, 0.0, 0.0, 0.0, 0.0, 0.0/),&
  dat_awc_lm2=&
  (/ 0.063, 0.132, 0.109, 0.098, 0.086, 0.120, 0.101, 0.445, 0.150, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_k_sat_ref=&
  (/ 0.021, .0036, .0018, .0087, .0061, .0026, .0051, .0036, .0036, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_psi_sat_ref=&
  (/ -.059, -0.28, -0.27, -0.13, -0.13, -0.27, -0.16, -0.28, -0.28, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_chb=&
  (/   3.5,   6.4,  11.0,   4.8,   6.3,   8.4,   6.3,   6.4,   6.4, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
!  dat_heat_capacity_ref =&
!  (/ 1.8e6, 2.0e6, 2.6e6, 1.9e6, 2.2e6, 2.3e6, 2.1e6, 3.0e6,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
! previous (ref) values were based on typical water contents
! following dry values are based on w_min=(1-w_sat) w_org=0
! except for peat, where            w_org-(1-w_sat) w_min=0
! microscopic rho*c for w_min is 2650*733 and for w_org is 1300*1926
! (brutsaert 1982 evaporation into the atmosphere p 146)
! ignored air
  dat_heat_capacity_dry =&
  (/ 1.2e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.4e6,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
!  dat_thermal_cond_ref =&
!  (/   1.5,   0.8,  1.35,  1.15, 1.475, 1.075, 1.217,  0.39, 2.e-7, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
! previous (ref) values were based on typical water contents
! following dry and sat values and interpolating exponents are based on
! computations after deVries. i computed C M and F functions for
! unfrozen soil in
! spreadsheet Research\LaD2\soil thermal conductivity. the dry and
! sat values come right out of those computations. the exponent was
! done by eye. curves look like typical literature curves.
! TEMP: still need to treat freezing, maybe import deVries model into code.
  dat_thermal_cond_dry =&
  (/  0.14,  0.21,  0.20,  .175, 0.170, 0.205, 0.183,  0.05, 2.e-7, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_thermal_cond_sat =&
  (/  2.30,  1.50,  1.50,  1.90, 1.900, 1.500, 1.767,  0.50, 2.e-7, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_thermal_cond_scale =&
  (/  15.0,  0.50,   10.,  2.74,  12.2,  2.24,  4.22,   1.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_thermal_cond_exp =&
  (/   3.0,   5.0,   6.0,   4.0,   4.5,   5.5, 4.667,   1.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_thermal_cond_weight =&
  (/  0.20,  0.70,   0.7,  0.45, 0.450, 0.700, 0.533,   1.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_emis_dry=&
  (/ 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_emis_sat=&
  (/ 0.980, 0.975, 0.970, .9775, 0.975, .9725, 0.975, 0.975,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_z0_momentum=&
  (/  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01, 0.045, 0.0, 0.0, 0.0, 0.0, 0.0 /),&
  dat_tf_depr=&
  (/  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  !Coarse  Medium   Fine    CM     CF     MF    CMF    Peat    MCM
integer, parameter :: peat_soil_type = 8
real :: dat_refl_dry_dir(n_dim_soil_types,NBANDS); data dat_refl_dry_dir &
   / 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0, 5*0.0,      & ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0, 5*0.0   /     ! NIR
real :: dat_refl_dry_dif(n_dim_soil_types,NBANDS); data dat_refl_dry_dif &
   / 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0, 5*0.0,      & ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0, 5*0.0   /     ! NIR
real :: dat_refl_sat_dir(n_dim_soil_types,NBANDS); data dat_refl_sat_dir &
   / 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0, 5*0.0,      & ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0, 5*0.0   /     ! NIR
real :: dat_refl_sat_dif(n_dim_soil_types,NBANDS); data dat_refl_sat_dif &
   / 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0, 5*0.0,      & ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0, 5*0.0   /     ! NIR
  !Coarse  Medium   Fine    CM     CF     MF    CMF    Peat    MCM
integer, dimension(n_dim_soil_types), public :: &
  input_cover_types=&
  (/ 1,     2,     3,     4,     5,     6,     7,     8,     100,  &
     101,   102,   103,   104,   105 /)
character(len=4), dimension(n_dim_soil_types) :: &
  tile_names=&
  (/'c   ','m   ','f   ','cm  ','cf  ','mf  ','cmf ','peat','mcm ', &
    'A   ','B   ','C   ','D   ','E   '/)
real, dimension(n_dim_soil_types) :: clay = &  ! Clay percentage (for calculating Qmax)
  (/  5.0,  15.0,  60.0,  10.0,  32.5,  37.5,  26.67,  0.0,   30.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

real :: peat_soil_e_depth = -1  ! If positive (and GW_TILED or GW_HILL), override soil_e_depth for peat
real :: peat_kx0 = -1           ! If non-negative (and GW_TILED or GW_HILL), override macroporosity for peat
namelist /soil_data_nml/ psi_wilt, &
     soil_to_use, soil_type_file, tile_names, input_cover_types, &
     comp, K_min, K_max_matrix, DThDP_max, psi_min, k_over_B, &
     rate_fc, sfc_heat_factor, z_sfc_layer, &
     sub_layer_tc_fac, z_sub_layer_min, z_sub_layer_max, freeze_factor, &
     aspect, zeta_mult, num_l,  dz,                      &
     use_lm2_awc,    lm2, use_alt_psi_for_rh, &
     use_tau_fix, use_sat_fix, use_comp_for_ic, use_comp_for_push, limit_hi_psi, use_fluid_ice, &
     use_alt3_soil_hydraulics, &
     limit_DThDP, &
     use_single_geo, geohydrology_to_use, &
     use_alpha, &
     retro_a0n1, &
     soil_index_constant,         &
     gw_res_time,            rsa_exp_global,      &
     gw_scale_length, gw_scale_relief, gw_scale_soil_depth, &
     slope_exp, &
     gw_scale_perm, z_ref, geothermal_heat_flux_constant, &
     Dpsi_min_const, &
     k0_macro_z, k0_macro_x, log_rho_max, &
     dat_w_sat,               dat_awc_lm2,     &
     dat_k_sat_ref,            &
     dat_psi_sat_ref,               dat_chb,          &
     dat_heat_capacity_dry,       dat_thermal_cond_dry,   &
     dat_thermal_cond_sat,        dat_thermal_cond_exp,   &
     dat_thermal_cond_scale,        dat_thermal_cond_weight,   &
     dat_refl_dry_dir,            dat_refl_sat_dir,              &
     dat_refl_dry_dif,            dat_refl_sat_dif,              &
     dat_emis_dry,              dat_emis_sat,                &
     dat_z0_momentum,           dat_tf_depr,      clay,           &
     peat_soil_e_depth,         peat_kx0
!---- end of namelist --------------------------------------------------------

real    :: gw_hillslope_length   = 1000.
real    :: gw_hillslope_relief   =  100.
real    :: gw_hillslope_a        = 0.
real    :: gw_hillslope_n        = 1.
real    :: gw_hillslope_zeta_bar =    0.5
real    :: gw_soil_e_depth       =    4.
real    :: gw_perm               = 3.e-14   ! nominal permeability, m^2
real :: gw_flux_table(26,31), gw_area_table(26,31)
real, allocatable :: log_rho_table(:,:,:,:,:)
real, allocatable :: log_deficit_list(:)
real, allocatable :: log_zeta_s(:)
real, allocatable :: log_tau(:)

real, dimension(31 ) :: gw_zeta_s       = &
  (/ 1.0000000e-5, 1.5848932e-5, 2.5118864e-5, 3.9810717e-5, 6.3095737e-5, &
     1.0000000e-4, 1.5848932e-4, 2.5118864e-4, 3.9810717e-4, 6.3095737e-4, &
     1.0000000e-3, 1.5848932e-3, 2.5118864e-3, 3.9810717e-3, 6.3095737e-3, &
     1.0000000e-2, 1.5848932e-2, 2.5118864e-2, 3.9810717e-2, 6.3095737e-2, &
     1.0000000e-1, 1.5848932e-1, 2.5118864e-1, 3.9810717e-1, 6.3095737e-1, &
     1.0000000e+0, 1.5848932e+0, 2.5118864e+0, 3.9810717e+0, 6.3095737e+0, &
     1.0000000e+1 /)

real, dimension(26) :: gw_storage_norm = &
  (/ 0.,      0.04000, 0.08000, 0.12000, 0.16000, 0.20000, &
     0.24000, 0.28000, 0.32000, 0.36000, 0.40000, 0.44000, &
     0.48000, 0.52000, 0.56000, 0.60000, 0.64000, 0.68000, &
     0.72000, 0.76000, 0.80000, 0.84000, 0.88000, 0.92000, &
     0.96000, 1.00000   /)
real, dimension(26) :: gw_flux_norm_zeta_s_04 = &
  (/ 0.0e000, 7.04e-6, 1.14e-5, 1.85e-5, 3.01e-5, 4.89e-5, &
     6.95e-5, 8.10e-5, 9.26e-5, 1.42e-4, 2.93e-4, 6.14e-4, &
     1.25e-3, 2.47e-3, 4.76e-3, 8.98e-3, 1.66e-2, 3.02e-2, &
     5.41e-2, 9.56e-2, 1.67e-1, 2.88e-1, 4.92e-1, 8.36e-1, &
     1.53e+0, 1.00e+1   /)
real, dimension(26) :: gw_area_norm_zeta_s_04 = &
  (/ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
     0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
     0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
     0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
     3.48e-1, 1.00000  /)

integer :: num_sfc_layers, sub_layer_min, sub_layer_max
integer :: num_storage_pts, num_zeta_pts, num_tau_pts

real    :: zfull (max_lev)    ! depth of the soil layer centers, m
real    :: zhalf (max_lev+1)  ! depth of the soil layer interfaces, m

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_soil_data_namelist(soil_single_geo, soil_gw_option )
  logical, intent(out) :: soil_single_geo
  integer, intent(out) :: soil_gw_option
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i, rcode, ncid, varid, dimids(3)

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=soil_data_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=soil_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'soil_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  unit=stdlog()
  write(unit, nml=soil_data_nml)

  ! register selector for all soil tiles
  call register_tile_selector('soil', long_name='soil',&
       tag = SEL_SOIL, idata1 = 0, area_depends_on_time=.FALSE. )
  ! register selectors for tile-specific diagnostics
  do i=1, n_dim_soil_types
     call register_tile_selector(tile_names(i), long_name='',&
          tag = SEL_SOIL, idata1 = i, area_depends_on_time=.FALSE. )
  enddo
  num_sfc_layers = 0
  sub_layer_min = 0
  sub_layer_max = 0

  zhalf(1) = 0
  do i = 1, num_l
    zhalf(i+1) = zhalf(i) + dz(i)
    zfull(i)   = 0.5*(zhalf(i+1) + zhalf(i))

    if (zhalf(i)   < z_sub_layer_min+1.e-4) sub_layer_min = i
    if (zhalf(i+1) < z_sfc_layer+1.e-4) num_sfc_layers = i
    if (zhalf(i+1) < z_sub_layer_max+1.e-4) sub_layer_max = i
  enddo

!!$  write (*,*) 'min/max index of layers whose thermal cond is scaled:',sub_layer_min,sub_layer_max

  if (trim(geohydrology_to_use)=='lm2') then
     gw_option = GW_LM2
  else if (trim(geohydrology_to_use)=='linear') then
     gw_option = GW_LINEAR
  else if (trim(geohydrology_to_use)=='hill_ar5') then
     gw_option = GW_HILL_AR5
  else if (trim(geohydrology_to_use)=='hill') then
     gw_option = GW_HILL
  else if (trim(geohydrology_to_use)=='tiled') then
     gw_option = GW_TILED
  else
     call error_mesg('read_soil_data_namelist',&
        'option geohydrology_to_use="'//trim(geohydrology_to_use)//'" is invalid, use '// &
        ' "lm2", "linear", "hill", "hill_ar5", or "tiled"', &
        FATAL)
  endif

  if ((gw_option==GW_LM2).neqv.lm2) then
     call error_mesg('read_soil_data_namelist',&
        'geohydrlogy/LM2 option conflict: geohydrology_to_use must be consistent with the LM2 flag', &
        FATAL)
  endif

  if (gw_option==GW_HILL_AR5.and..not.use_single_geo) then
     num_storage_pts = 26
     call read_data('INPUT/geohydrology_table.nc', 'gw_flux_norm', &
                 gw_flux_table, no_domain=.true.)
     call read_data('INPUT/geohydrology_table.nc', 'gw_area_norm', &
                 gw_area_table, no_domain=.true.)
  else if (gw_option==GW_HILL) then
     __NF_ASRT__(nf_open('INPUT/geohydrology_table_2a2n.nc',NF_NOWRITE,ncid))
     __NF_ASRT__(nf_inq_varid(ncid,'log_rho_a0n1',varid))
     __NF_ASRT__(nf_inq_vardimid(ncid,varid,dimids))
     __NF_ASRT__(nf_inq_dimlen(ncid,dimids(1),num_storage_pts))
     __NF_ASRT__(nf_inq_dimlen(ncid,dimids(2),num_tau_pts))
     __NF_ASRT__(nf_inq_dimlen(ncid,dimids(3),num_zeta_pts))
     __NF_ASRT__(nf_close(ncid))

     allocate (log_rho_table(num_storage_pts, num_tau_pts, num_zeta_pts, 2, 2))
     allocate (log_deficit_list(num_storage_pts))
     allocate (log_tau(num_tau_pts))
     allocate (log_zeta_s(num_zeta_pts))

     if (.not.retro_a0n1) then
         call read_data('INPUT/geohydrology_table_2a2n.nc', 'log_rho_a0n1', &
                 log_rho_table(:,:,:,1,1), no_domain=.true.)
     else
         call read_data('INPUT/geohydrology_table_2a2n.nc', 'retro_log_rho_a0n1', &
                 log_rho_table(:,:,:,1,1), no_domain=.true.)
     endif
     call read_data('INPUT/geohydrology_table_2a2n.nc', 'log_rho_a0n2', &
             log_rho_table(:,:,:,1,2), no_domain=.true.)
     call read_data('INPUT/geohydrology_table_2a2n.nc', 'log_rho_a1n1', &
             log_rho_table(:,:,:,2,1), no_domain=.true.)
     call read_data('INPUT/geohydrology_table_2a2n.nc', 'log_rho_a1n2', &
             log_rho_table(:,:,:,2,2), no_domain=.true.)
     call read_data('INPUT/geohydrology_table_2a2n.nc', 'log_deficit', &
             log_deficit_list, no_domain=.true.)
     call read_data('INPUT/geohydrology_table_2a2n.nc', 'log_tau', &
             log_tau, no_domain=.true.)
     call read_data('INPUT/geohydrology_table_2a2n.nc', 'log_zeta_s', &
             log_zeta_s, no_domain=.true.)
  endif


  ! set up output arguments
  soil_single_geo = use_single_geo
  soil_gw_option  = gw_option

end subroutine read_soil_data_namelist


! ============================================================================
function soil_tile_ctor(tag, hidx_j, hidx_k) result(ptr)
  type(soil_tile_type), pointer :: ptr ! return value
  integer, intent(in)  :: tag ! kind of tile
  integer, intent(in)  :: hidx_j, hidx_k ! hillslope indices

  integer :: i

  allocate(ptr)
  ptr%tag = tag
  ptr%hidx_j = hidx_j
  ptr%hidx_k = hidx_k
  ! allocate storage for tile data
  allocate( ptr%wl                (num_l),  &
            ptr%ws                (num_l),  &
            ptr%T                 (num_l),  &
            ptr%groundwater       (num_l),  &
            ptr%groundwater_T     (num_l),  &
            ptr%w_fc              (num_l),  &
            ptr%w_wilt            (num_l),  &
            ptr%d_trans           (num_l),  &
            ptr%alpha             (num_l),  &
            ptr%k_macro_z         (num_l),  &
            ptr%k_macro_x         (num_l),  &
            ptr%vwc_max           (num_l),  &
            ptr%heat_capacity_dry (num_l),  &
            ptr%e                 (num_l),  &
            ptr%f                 (num_l),  &
            ptr%psi               (num_l),  &
            ptr%gw_flux_norm      (num_storage_pts),  &
            ptr%gw_area_norm      (num_storage_pts),  &
            ptr%hyd_cond_horz     (num_l),  &
            ptr%div_hlsp          (num_l),  &
            ptr%div_hlsp_heat     (num_l),  &
            ptr%fast_soil_C       (num_l),  &
            ptr%slow_soil_C       (num_l),  &
            ptr%fsc_in            (num_l),  &
            ptr%ssc_in            (num_l),  &
            ptr%asoil_in          (num_l),  &
            ptr%is_peat           (num_l),  &
            ptr%soil_organic_matter(num_l),  &
            ptr%div_hlsp_DOC      (N_C_TYPES, num_l), &
            ptr%div_hlsp_DON      (N_C_TYPES, num_l), &
            ptr%div_hlsp_NO3   (num_l) , &
            ptr%div_hlsp_NH4   (num_l)         )

  ! Initialize to catch use before appropriate
  !ptr%psi(:) = initval
  ptr%hyd_cond_horz(:) = initval
  ptr%div_hlsp(:)      = initval
  ptr%div_hlsp_heat(:) = initval
  ptr%div_hlsp_DOC(:,:) = initval
  ptr%div_hlsp_DON(:,:) = initval
  ptr%div_hlsp_NO3(:) = initval
  ptr%div_hlsp_NH4(:) = initval

  call soil_data_init_0d(ptr)
  do i=1,num_l
     call init_soil_pool(ptr%soil_organic_matter(i), Qmax=ptr%pars%Qmax)
  enddo
  do i = 1,N_LITTER_POOLS
     call init_soil_pool(ptr%litter(i), protectionRate=0.0, Qmax=0.0, max_cohorts=1)
  enddo
end function soil_tile_ctor


! ============================================================================
function soil_tile_copy_ctor(soil) result(ptr)
  type(soil_tile_type), pointer :: ptr ! return value
  type(soil_tile_type), intent(in) :: soil ! tile to copy

  allocate(ptr)
  ptr = soil ! copy all non-pointer members
  ! no need to allocate storage for allocatable components of the type, because
  ! F2003 takes care of that, and also takes care of copying data
end function soil_tile_copy_ctor


! ============================================================================
subroutine delete_soil_tile(ptr)
  type(soil_tile_type), pointer :: ptr

  ! no need to deallocate components of soil_tile, because F2003 takes care of
  ! allocatable components deallocation when soil_tile is deallocated
  deallocate(ptr)
end subroutine delete_soil_tile


! ============================================================================
subroutine soil_data_init_0d(soil)
  type(soil_tile_type), intent(inout) :: soil

  integer :: k, l
  real    :: comp_local
  real    :: z ! depth at top of current layer

  k = soil%tag

  soil%pars%vwc_sat           = dat_w_sat            (k)
  soil%pars%awc_lm2           = dat_awc_lm2          (k)
  soil%pars%k_sat_ref         = dat_k_sat_ref        (k)
  soil%pars%psi_sat_ref       = dat_psi_sat_ref      (k)
  soil%pars%chb               = dat_chb              (k)
  soil%pars%heat_capacity_dry = dat_heat_capacity_dry(k)
  soil%pars%thermal_cond_dry  = dat_thermal_cond_dry (k)
  soil%pars%thermal_cond_sat  = dat_thermal_cond_sat (k)
  soil%pars%thermal_cond_exp  = dat_thermal_cond_exp (k)
  soil%pars%thermal_cond_scale  = dat_thermal_cond_scale (k)
  soil%pars%thermal_cond_weight  = dat_thermal_cond_weight (k)
  soil%pars%refl_dry_dir      = dat_refl_dry_dir     (k,:)
  soil%pars%refl_dry_dif      = dat_refl_dry_dif     (k,:)
  soil%pars%refl_sat_dir      = dat_refl_sat_dir     (k,:)
  soil%pars%refl_sat_dif      = dat_refl_sat_dif     (k,:)
  soil%pars%emis_dry          = dat_emis_dry         (k)
  soil%pars%emis_sat          = dat_emis_sat         (k)
  soil%pars%z0_momentum       = dat_z0_momentum      (k)
  soil%pars%tfreeze           = tfreeze - dat_tf_depr(k)
  soil%pars%rsa_exp           = rsa_exp_global
  soil%pars%tau_groundwater   = gw_res_time
  soil%pars%hillslope_length  = gw_hillslope_length*gw_scale_length
  soil%pars%hillslope_zeta_bar= gw_hillslope_zeta_bar
  soil%pars%hillslope_relief  = gw_hillslope_relief*gw_scale_relief
  soil%pars%soil_e_depth      = gw_soil_e_depth*gw_scale_soil_depth
  soil%pars%hillslope_a       = gw_hillslope_a
  soil%pars%hillslope_n       = gw_hillslope_n
  soil%pars%k_sat_gw          = gw_perm*gw_scale_perm*9.8e9  ! m^2 to kg/(m2 s)
  soil%pars%storage_index     = 1
  soil%alpha                  = 1.0
  soil%fast_soil_C(:)         = 0.0
  soil%slow_soil_C(:)         = 0.0
  soil%asoil_in(:)            = 0.0
  soil%is_peat(:)             = 0
  soil%fsc_in(:)              = 0.0
  soil%ssc_in(:)              = 0.0

  soil%fast_DOC_leached      = 0.0
  soil%slow_DOC_leached      = 0.0
  soil%deadmic_DOC_leached   = 0.0

  comp_local = 0.0
  if (use_comp_for_push) comp_local = comp
  do l = 1, num_l
    z=zfull(l)
    soil%k_macro_z(l) = k0_macro_z &
                        * exp(-z/soil%pars%soil_e_depth)
    soil%k_macro_x(l) = k0_macro_x &
                        * exp(-z/soil%pars%soil_e_depth)
    soil%vwc_max(l)   = soil%pars%vwc_sat + comp_local*zfull(l)
  enddo

  select case(gw_option)
  case(GW_HILL_AR5)
     if (use_single_geo) then
        soil%gw_flux_norm = gw_flux_norm_zeta_s_04
        soil%gw_area_norm = gw_area_norm_zeta_s_04
     endif
     soil%pars%zeta = soil%pars%soil_e_depth / soil%pars%hillslope_relief
  case(GW_HILL, GW_TILED)
     if (use_single_geo) &
         call soil_data_init_derive_subsurf_pars ( soil )
  case default
     ! do nothing
  end select

  ! ---- derived constant soil parameters
  ! w_fc (field capacity) set to w at which hydraulic conductivity equals
  ! a nominal drainage rate "rate_fc"
  ! w_wilt set to w at which psi is psi_wilt
  if (use_lm2_awc) then
     soil%w_wilt(:) = 0.15
     soil%w_fc  (:) = 0.15 + soil%pars%awc_lm2
  else
     soil%w_wilt(:) = soil%pars%vwc_sat &
          *(soil%pars%psi_sat_ref/(psi_wilt*soil%alpha(:)))**(1/soil%pars%chb)
     soil%w_fc  (:) = soil%pars%vwc_sat &
          *(rate_fc/(soil%pars%k_sat_ref*soil%alpha(:)**2))**(1/(3+2*soil%pars%chb))
  endif

  soil%pars%vwc_wilt = soil%w_wilt(1)
  soil%pars%vwc_fc   = soil%w_fc  (1)

  soil%pars%vlc_min = soil%pars%vwc_sat*K_rel_min**(1/(3+2*soil%pars%chb))

  soil%z0_scalar = soil%pars%z0_momentum * exp(-k_over_B)

  soil%geothermal_heat_flux = geothermal_heat_flux_constant

  ! Init hlsp variables to initval to flag use before appropriate initialization
  !soil%pars%soil_e_depth = initval ! ZMS Let this keep gridcell value for now.
  !Also needed if .not. gw_option /= GW_TILE
  soil%pars%microtopo = initval
  soil%pars%tile_hlsp_length = initval
  soil%pars%tile_hlsp_slope = initval
  soil%pars%tile_hlsp_elev = initval
  soil%pars%tile_hlsp_hpos = initval
  soil%pars%tile_hlsp_width = initval
!  soil%pars%transm_bedrock = initval

  !Qmax in mgC/kg soil from Mayes et al 2012, converted to g/m3 using solid density of 2650 kg/m3
  if(clay(k) .le. 0) then
      soil%pars%Qmax= 0
  else
      soil%pars%Qmax = max(0.0,10**(.4833*log10(clay(k))+2.3282)*(1.0-dat_w_sat(k))*2650*1e-6)
  endif
end subroutine soil_data_init_0d

! ============================================================================
subroutine soil_data_init_derive_subsurf_pars ( soil )
  type(soil_tile_type), intent(inout) :: soil
  integer i, l, code, m_zeta, m_tau, m_a, m_n
  real :: alpha_inf_sq, alpha_sfc_sq
  real :: single_log_zeta_s, single_log_tau, frac_zeta, frac_tau

  real :: k_macro_x_local

  k_macro_x_local = k0_macro_x
  ! Override hydraulic properties for peat.
  if (soil%tag == peat_soil_type) then
     if (peat_soil_e_depth > 0.) soil%pars%soil_e_depth = peat_soil_e_depth
     if (peat_kx0 >= 0.) then
        do l=1,num_l
           soil%k_macro_x(l) = peat_kx0 &
                             * exp(-zfull(l)/soil%pars%soil_e_depth)
        end do
        k_macro_x_local = peat_kx0
     end if
  end if

  if (use_alpha) then
      soil%pars%k_sat_sfc = (soil%pars%k_sat_ref-soil%pars%k_sat_gw) &
                * exp(z_ref/soil%pars%soil_e_depth) + soil%pars%k_sat_gw
      alpha_inf_sq = soil%pars%k_sat_gw / soil%pars%k_sat_ref
      alpha_sfc_sq = soil%pars%k_sat_sfc / soil%pars%k_sat_ref
      do l = 1, num_l
        soil%alpha(l) = sqrt(alpha_inf_sq+(alpha_sfc_sq-alpha_inf_sq) &
                    *exp(-zfull(l)/soil%pars%soil_e_depth))
      enddo
  else
      soil%pars%k_sat_sfc = soil%pars%k_sat_ref
      soil%alpha = 1.0
  endif

  soil%pars%tau =    &
    (soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length) &
     / ((soil%pars%k_sat_sfc+k_macro_x_local)*soil%pars%soil_e_depth)

  soil%d_trans(1) = 1.
  do l = 1, num_l-1
     soil%d_trans(l+1) = exp(-zhalf(l+1)/soil%pars%soil_e_depth)
     soil%d_trans(l)   = soil%d_trans(l) - soil%d_trans(l+1)
  enddo
  soil%d_trans = soil%d_trans &
       * (soil%pars%k_sat_sfc + k_macro_x_local ) &
       * soil%pars%soil_e_depth
  soil%d_trans(num_l) = soil%d_trans(num_l) &
       + soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length

  if (soil%pars%hillslope_relief.le.0.) soil%pars%hillslope_relief = 1.e-10
  soil%pars%zeta = soil%pars%soil_e_depth &
                                        / soil%pars%hillslope_relief
  single_log_zeta_s = log10(soil%pars%zeta)
  single_log_zeta_s = max(single_log_zeta_s, log_zeta_s(1))
  single_log_zeta_s = min(single_log_zeta_s, log_zeta_s(num_zeta_pts))
  m_zeta = num_zeta_pts / 2
  code = 0
  do while (code.eq.0)
     if (single_log_zeta_s .lt. log_zeta_s(m_zeta)) then
        m_zeta = m_zeta - 1
     else if (single_log_zeta_s .gt. log_zeta_s(m_zeta+1)) then
        m_zeta = m_zeta + 1
     else
        code = 1
     endif
  enddo
  frac_zeta = (single_log_zeta_s - log_zeta_s(m_zeta)) &
     / (log_zeta_s(m_zeta+1) - log_zeta_s(m_zeta))
  single_log_tau = log10( soil%pars%tau )
  single_log_tau = max(single_log_tau, log_tau(1))
  single_log_tau = min(single_log_tau, log_tau(num_tau_pts))
  m_tau = num_tau_pts / 2
  code = 0
  do while (code.eq.0)
     if (single_log_tau .lt. log_tau(m_tau)) then
        m_tau = m_tau - 1
     else if (single_log_tau .gt. log_tau(m_tau+1)) then
        m_tau = m_tau + 1
     else
        code = 1
     endif
  enddo
  frac_tau = (single_log_tau - log_tau(m_tau)) &
     / (log_tau(m_tau+1) - log_tau(m_tau))
  m_a = 1+int(soil%pars%hillslope_a+0.1)
  m_n = int(soil%pars%hillslope_n+0.1)
  if (abs(real(m_a-1)-soil%pars%hillslope_a).gt.1.e-6 .or. &
      abs(real(m_n)-soil%pars%hillslope_n).gt.1.e-6 .or. &
      (m_a.ne.1 .and. m_a.ne.2) .or. (m_n.ne.1 .and. m_n.ne.2)) then
      write (*,*) 'a,n,m_a,m_n=', soil%pars%hillslope_a, &
          soil%pars%hillslope_n, m_a, m_n
      call error_mesg('soil_data_init_derive_subsurf_pars', &
          'hillslope a or n value is unacceptable', FATAL)
  endif
  do i = 1, num_storage_pts
    soil%gw_flux_norm(i) = min ( log_rho_max,                      &
       (1.-frac_zeta)*(1.-frac_tau)*log_rho_table(i , m_tau  , m_zeta  , m_a, m_n) &
     + (   frac_zeta)*(1.-frac_tau)*log_rho_table(i , m_tau  , m_zeta+1, m_a, m_n) &
     + (1.-frac_zeta)*(   frac_tau)*log_rho_table(i , m_tau+1, m_zeta  , m_a, m_n) &
     + (   frac_zeta)*(   frac_tau)*log_rho_table(i , m_tau+1, m_zeta+1, m_a, m_n) )
  enddo

end subroutine soil_data_init_derive_subsurf_pars


! ============================================================================
subroutine soil_data_init_derive_subsurf_pars_tiled ( soil, use_geohydrodata)
  type(soil_tile_type), intent(inout) :: soil
  logical, intent(in) :: use_geohydrodata ! use gcell data from geohydrology.nc to set soil properties
  ! else values will have been set by hillslope.nc in hlsp_init.
  integer i, l, code, m_zeta, m_tau
  real :: alpha_inf_sq, alpha_sfc_sq
  real :: single_log_zeta_s, single_log_tau, frac_zeta, frac_tau
  real :: k_macro_x_local

  k_macro_x_local = k0_macro_x
  ! Override hydraulic properties for peat.
  if (soil%tag == peat_soil_type) then
     if (peat_soil_e_depth > 0.) soil%pars%soil_e_depth = peat_soil_e_depth
     if (peat_kx0 >= 0.) then
        do l=1,num_l
           soil%k_macro_x(l) = peat_kx0 &
                             * exp(-zfull(l)/soil%pars%soil_e_depth)
        end do
        k_macro_x_local = peat_kx0
     end if
  end if

  if (use_alpha) then
      soil%pars%k_sat_sfc = (soil%pars%k_sat_ref-soil%pars%k_sat_gw) &
                * exp(z_ref/soil%pars%soil_e_depth) + soil%pars%k_sat_gw
      alpha_inf_sq = soil%pars%k_sat_gw / soil%pars%k_sat_ref
      alpha_sfc_sq = soil%pars%k_sat_sfc / soil%pars%k_sat_ref
      do l = 1, num_l
        soil%alpha(l) = sqrt(alpha_inf_sq+(alpha_sfc_sq-alpha_inf_sq) &
                    *exp(-zfull(l)/soil%pars%soil_e_depth))
      enddo
  else
      soil%pars%k_sat_sfc = soil%pars%k_sat_ref
      soil%alpha = 1.0
  endif

  if (use_geohydrodata) then
     soil%pars%tau =    &
       (soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length) &
        / ((soil%pars%k_sat_sfc+k_macro_x_local)*soil%pars%soil_e_depth)

!     ZMS: d_trans not needed
!     z_local = 0.
!     soil%d_trans(1) = 1.
!     do l = 1, num_l-1
!       z_local = z_local + dz(l)
!       soil%d_trans(l+1) = exp(-z_local/soil%pars%soil_e_depth)
!       soil%d_trans(l)   = soil%d_trans(l) - soil%d_trans(l+1)
!     enddo
!     soil%d_trans = soil%d_trans &
!          * (soil%pars%k_sat_sfc + k_macro_constant ) &
!          * soil%pars%soil_e_depth
!     soil%d_trans(num_l) = soil%d_trans(num_l) &
!          + soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length

     if (soil%pars%hillslope_relief.le.0.) soil%pars%hillslope_relief = 1.e-10
     soil%pars%zeta = soil%pars%soil_e_depth &
                                           / soil%pars%hillslope_relief
!     single_log_zeta_s = log10(soil%pars%zeta)
!     single_log_zeta_s = max(single_log_zeta_s, log_zeta_s(1))
!     single_log_zeta_s = min(single_log_zeta_s, log_zeta_s(num_zeta_pts))
!     m_zeta = num_zeta_pts / 2
!     code = 0
!     do while (code.eq.0)
!       if (single_log_zeta_s .lt. log_zeta_s(m_zeta)) then
!           m_zeta = m_zeta - 1
!         else if (single_log_zeta_s .gt. log_zeta_s(m_zeta+1)) then
!           m_zeta = m_zeta + 1
!         else
!           code = 1
!       endif
!     enddo
!     frac_zeta = (single_log_zeta_s - log_zeta_s(m_zeta)) &
!        / (log_zeta_s(m_zeta+1) - log_zeta_s(m_zeta))
!     single_log_tau = log10( soil%pars%tau )
!     single_log_tau = max(single_log_tau, log_tau(1))
!     single_log_tau = min(single_log_tau, log_tau(num_tau_pts))
!     m_tau = num_tau_pts / 2
!     code = 0
!     do while (code.eq.0)
!       if (single_log_tau .lt. log_tau(m_tau)) then
!           m_tau = m_tau - 1
!         else if (single_log_tau .gt. log_tau(m_tau+1)) then
!           m_tau = m_tau + 1
!         else
!           code = 1
!       endif
!     enddo
!     frac_tau = (single_log_tau - log_tau(m_tau)) &
!        / (log_tau(m_tau+1) - log_tau(m_tau))
!     do i = 1, num_storage_pts
!       soil%gw_flux_norm(i) = min ( log_rho_max,                      &
!          (1.-frac_zeta)*(1.-frac_tau)*log_rho_table(i , m_tau  , m_zeta  ) &
!        + (   frac_zeta)*(1.-frac_tau)*log_rho_table(i , m_tau  , m_zeta+1) &
!        + (1.-frac_zeta)*(   frac_tau)*log_rho_table(i , m_tau+1, m_zeta  ) &
!        + (   frac_zeta)*(   frac_tau)*log_rho_table(i , m_tau+1, m_zeta+1) )
!     enddo
  else
     soil%pars%tau =    &
       (soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length) &
        / ((soil%pars%k_sat_sfc+k_macro_x_local)*soil%pars%soil_e_depth)
        ! ZMS Check this.
     soil%pars%zeta = soil%pars%soil_e_depth &
                                                 / soil%pars%hillslope_relief
     ! ... Add more?
  end if ! use_geohydrodata

end subroutine soil_data_init_derive_subsurf_pars_tiled



! ============================================================================
subroutine soil_data_init_derive_subsurf_pars_ar5 ( soil )
  type(soil_tile_type), intent(inout) :: soil
  integer i, code, m_zeta
  real :: frac_zeta
  if (soil%pars%hillslope_relief.le.0.) &
     soil%pars%hillslope_relief =      &
        soil%pars%soil_e_depth / 10.**log_zeta_s(num_zeta_pts)
  soil%pars%zeta = soil%pars%soil_e_depth &
     / soil%pars%hillslope_relief
  soil%pars%zeta = max(soil%pars%zeta, gw_zeta_s(1))
  soil%pars%zeta = min(soil%pars%zeta, gw_zeta_s(31))
  m_zeta = 31 / 2
  code = 0
  do while (code.eq.0)
     if (soil%pars%zeta .lt. gw_zeta_s(m_zeta)) then
        m_zeta = m_zeta - 1
     else if (soil%pars%zeta .gt. gw_zeta_s(m_zeta+1)) then
        m_zeta = m_zeta + 1
     else
        code = 1
     endif
  enddo
  frac_zeta = (soil%pars%zeta - gw_zeta_s(m_zeta)) &
                / (gw_zeta_s(m_zeta+1) - gw_zeta_s(m_zeta))
  do i = 1, num_storage_pts
     soil%gw_flux_norm(i) = gw_flux_table(i,m_zeta) &
          + frac_zeta*(gw_flux_table(i,m_zeta+1)-gw_flux_table(i,m_zeta))
     soil%gw_area_norm(i) = gw_area_table(i,m_zeta) &
          + frac_zeta*(gw_area_table(i,m_zeta+1)-gw_area_table(i,m_zeta))
   enddo
end subroutine soil_data_init_derive_subsurf_pars_ar5


! =============================================================================
function soil_tiles_can_be_merged(soil1,soil2) result(response)
  logical :: response
  type(soil_tile_type), intent(in) :: soil1,soil2

  response = (soil1%tag==soil2%tag .and. soil1%hidx_j==soil2%hidx_j &
                .and. soil1%hidx_k==soil2%hidx_k)
! ZMS: If duplication of entire hillslopes is allowed to occur with land cover change (hillslope_topo_subdiv in
! hillslope_mod), then this may need to be modified. Could allow merging if mod(hidx_k,max_num_topo_hlsps) is
! the same for both, though this would deplete tiles in one hillslope, so it may be better addressed
! with a specific function that merges entire hillslopes.
end function


! =============================================================================
subroutine merge_soil_tiles(s1,w1,s2,w2)
  type(soil_tile_type), intent(in) :: s1
  type(soil_tile_type), intent(inout) :: s2
  real                , intent(in) :: w1,w2

  ! ---- local vars
  real    :: x1, x2 ! normalized relative weights
  real    :: gw, HEAT1, HEAT2 ! temporaries for groundwater and heat
  integer :: i

  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  ! combine state variables
  do i = 1,num_l
     ! calculate heat content at this level for both source tiles
     HEAT1 = &
          (s1%heat_capacity_dry(i)*dz(i)+clw*s1%Wl(i)+csw*s1%Ws(i))* &
          (s1%T(i)-tfreeze)
     HEAT2 = &
          (s2%heat_capacity_dry(i)*dz(i)+clw*s2%Wl(i)+csw*s2%Ws(i))* &
          (s2%T(i)-tfreeze)
     ! merge the amounts of water
     s2%Wl(i) = x1*s1%Wl(i) + x2*s2%Wl(i)
     s2%Ws(i) = x1*s1%Ws(i) + x2*s2%Ws(i)
     ! if the dry heat capacity of merged soil is to be changed, do it here
     ! ...
     ! calculate the merged temperature based on heat content
     s2%T(i) = tfreeze + (x1*HEAT1+x2*HEAT2)/ &
          (s2%heat_capacity_dry(i)*dz(i)+clw*s2%Wl(i)+csw*s2%Ws(i))

     ! calculate combined groundwater content
     gw = s1%groundwater(i)*x1 + s2%groundwater(i)*x2
     ! calculate combined groundwater temperature
     if (gw/=0) then
        s2%groundwater_T(i) = ( &
             s1%groundwater(i)*x1*(s1%groundwater_T(i)-tfreeze) + &
             s2%groundwater(i)*x2*(s2%groundwater_T(i)-tfreeze)   &
             ) / gw + tfreeze
     else
        s2%groundwater_T(i) = &
             s1%groundwater_T(i)*x1 + s2%groundwater_T(i)*x2
     endif
     s2%groundwater(i) = gw

  enddo
  s2%uptake_T    = s1%uptake_T*x1 + s2%uptake_T*x2
  ! merge soil carbon
  s2%fast_soil_C(:) = s1%fast_soil_C(:)*x1 + s2%fast_soil_C(:)*x2
  s2%slow_soil_C(:) = s1%slow_soil_C(:)*x1 + s2%slow_soil_C(:)*x2
  do i=1,num_l
    call combine_pools(s1%soil_organic_matter(i),s2%soil_organic_matter(i),w1,w2)
  enddo
  !is_peat is 1 or 0, so multiplying is like an AND operation
  s2%is_peat(:) = s1%is_peat(:) * s2%is_peat(:)
  do i = 1, N_LITTER_POOLS
     call combine_pools(s1%litter(i),s2%litter(i),w1,w2)
  enddo
  s2%asoil_in(:)    = s1%asoil_in(:)*x1 + s2%asoil_in(:)*x2
  s2%fsc_in(:)      = s1%fsc_in(:)*x1 + s2%fsc_in(:)*x2
  s2%ssc_in(:)      = s1%ssc_in(:)*x1 + s2%ssc_in(:)*x2

  s2%fast_DOC_leached=s1%fast_DOC_leached*x1 + s2%fast_DOC_leached*x2
  s2%slow_DOC_leached=s1%slow_DOC_leached*x1 + s2%slow_DOC_leached*x2
  s2%deadmic_DOC_leached=s1%deadmic_DOC_leached*x1 + s2%deadmic_DOC_leached*x2
end subroutine merge_soil_tiles

! =============================================================================
! returns true if tile fits the specified selector
function soil_is_selected(soil, sel)
  logical soil_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(soil_tile_type),      intent(in) :: soil

  soil_is_selected = (sel%idata1==0).or.(sel%idata1==soil%tag)
end function


! ============================================================================
! returns tag of the tile
function get_soil_tile_tag(soil) result(tag)
  integer :: tag
  type(soil_tile_type), intent(in) :: soil

  tag = soil%tag
end function



! ============================================================================
! compute average soil temperature with a given depth scale
function soil_ave_temp(soil, depth) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)                 :: depth ! averaging depth

  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  integer :: k

  A = 0 ; N = 0
  do k = 1, num_l
     w = dz(k) * exp(-zfull(k)/depth)
     A = A + soil%T(k) * w
     N = N + w
     if (zhalf(k+1).gt.depth) exit
  enddo
  A = A/N
end function soil_ave_temp


! ============================================================================
! compute average soil moisture with a given depth scale
function soil_ave_theta0(soil, zeta) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)                 :: zeta ! root depth scale
  real :: zeta2  !  adjusted length scale for soil-water-stress
  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  integer :: k

  A = 0 ; N = 0
  zeta2 = zeta*zeta_mult
  do k = 1, num_l
     w = exp(-zhalf(k)/zeta2)-exp(-zhalf(k+1)/zeta2)
     A = A + max(soil%wl(k)/(dens_h2o*dz(k))-soil%w_wilt(k),0.0)/&
          (soil%w_fc(k)-soil%w_wilt(k)) * w
     N = N + w
  enddo
  A = A/N
end function soil_ave_theta0


! ============================================================================
function soil_ave_theta1(soil, depth) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)                 :: depth ! averaging depth

  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  integer :: k

  A = 0 ; N = 0
  do k = 1, num_l
     w = dz(k) * exp(-zfull(k)/depth)
     A = A +min(max(soil%wl(k)/(dens_h2o*dz(k)),0.0)/&
          (soil%pars%vwc_sat),1.0) * w
     N = N + w
     if (zhalf(k+1).gt.depth) exit
  enddo
  A = A/N
end function soil_ave_theta1


! ============================================================================
! returns soil surface "wetness" -- fraction of the pores filled with water
subroutine soil_ave_wetness(soil, depth, SW, SI)
  type(soil_tile_type), intent(in) :: soil
  real                , intent(in) :: depth ! averaging depth
  real                , intent(out):: SW ! water volume / saturation
  real                , intent(out):: SI ! ice volume / saturation
  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  real    :: z ! current depth, m
  integer :: k

  SW = 0 ; SI = 0; N = 0 ; z = 0
  do k = 1, num_l
     w = dz(k) * exp(-(z+dz(k)/2)/depth)
     SW = SW + max(soil%wl(k)/(dens_h2o*dz(k)*soil%pars%vwc_sat),0.0) * w
     SI = SI + max(soil%ws(k)/(dens_h2o*dz(k)*soil%pars%vwc_sat),0.0) * w
     N = N + w
     z = z + dz(k)
  enddo
  SW = SW/N ; SI = SI/N
end subroutine soil_ave_wetness


! ============================================================================
! returns array of soil moisture
function soil_theta(soil) result (theta1)
  type(soil_tile_type), intent(in) :: soil
  real :: theta1(num_l)

  theta1(:) = min(max(soil%wl(:)/(dens_h2o*dz(1:num_l)),0.0)/(soil%pars%vwc_sat),1.0)
end function soil_theta


! ============================================================================
! Like soil_theta1 but for ice-filled porosity.
function soil_ice_porosity(soil) result(ice_porosity)
    type(soil_tile_type), intent(in) :: soil
    real :: ice_porosity(num_l)

    integer :: k

    do k=1, num_l
        ice_porosity(k) = min(max(soil%ws(k)/(dens_h2o*dz(k)),0.0)/(soil%pars%vwc_sat),1.0)
    enddo
end function soil_ice_porosity


! ============================================================================
function soil_ave_ice_porosity(soil,depth) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)                 :: depth ! averaging depth

  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  real    :: z ! current depth, m
  integer :: k

  A = 0 ; N = 0 ; z = 0
  do k = 1, num_l
     w = dz(k) * exp(-(z+dz(k)/2)/depth)
     A = A +min(max(soil%ws(k)/(dens_h2o*dz(k)),0.0)/&
          (soil%pars%vwc_sat),1.0) * w
     N = N + w
     z = z + dz(k)
     if (z.gt.depth) exit
  enddo
  A = A/N
end function soil_ave_ice_porosity


! ============================================================================
function soil_psi_stress(soil, zeta) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real,                 intent(in) :: zeta  ! root-mass depth scale
  real :: zeta2  !  adjusted length scale for soil-water-stress
  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  integer :: k

  A = 0 ; N = 0
  zeta2 = zeta*zeta_mult
  do k = 1, num_l
     w = exp(-zhalf(k)/zeta2)-exp(-zhalf(k+1)/zeta2)
     A = A + w * soil%psi(k)/psi_wilt
     N = N + w
  enddo
  A = A/N

end function soil_psi_stress

! ============================================================================
! compute bare-soil albedo, bare-soil emissivity, bare-soil roughness
! for scalar transport, and beta function
subroutine soil_data_radiation ( soil, cosz, use_brdf, soil_alb_dir, soil_alb_dif, soil_emis )
  type(soil_tile_type), intent(in)  :: soil
  real,                 intent(in)  :: cosz
  logical,              intent(in)  :: use_brdf
  real,                 intent(out) :: soil_alb_dir(NBANDS), soil_alb_dif(NBANDS), soil_emis
  ! ---- local vars
  real :: soil_sfc_vlc, blend, dry_value(NBANDS), sat_value(NBANDS)
  real :: zenith_angle, zsq, zcu

  soil_sfc_vlc  = soil%wl(1)/(dens_h2o*dz(1))
  blend         = max(0., min(1., soil_sfc_vlc/soil%pars%vwc_sat))
  if (use_brdf) then
     zenith_angle = acos(cosz)
     zsq = zenith_angle*zenith_angle
     zcu = zenith_angle*zsq
     dry_value =  soil%pars%f_iso_dry*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                + soil%pars%f_vol_dry*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                + soil%pars%f_geo_dry*(g0_geo+g1_geo*zsq+g2_geo*zcu)
     sat_value =  soil%pars%f_iso_sat*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                + soil%pars%f_vol_sat*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                + soil%pars%f_geo_sat*(g0_geo+g1_geo*zsq+g2_geo*zcu)
  else
     dry_value = soil%pars%refl_dry_dir
     sat_value = soil%pars%refl_sat_dir
  endif
  soil_alb_dir  = dry_value              + blend*(sat_value             -dry_value)
  soil_alb_dif  = soil%pars%refl_dry_dif + blend*(soil%pars%refl_sat_dif-soil%pars%refl_dry_dif)
  soil_emis     = soil%pars%emis_dry     + blend*(soil%pars%emis_sat    -soil%pars%emis_dry    )
end subroutine soil_data_radiation


! ============================================================================
! compute bare-soil albedo, bare-soil emissivity, bare-soil roughness
! for scalar transport, and beta function
subroutine soil_roughness ( soil, soil_z0s, soil_z0m )
  type(soil_tile_type), intent(in)  :: soil
  real,                 intent(out) :: soil_z0s, soil_z0m

  soil_z0s = soil%z0_scalar
  soil_z0m = soil%pars%z0_momentum
end subroutine soil_roughness

! ============================================================================
! compute soil thermodynamic properties.
subroutine soil_data_thermodynamics ( soil, vlc, vsc, &
                                      soil_E_max, thermal_cond)
  type(soil_tile_type), intent(inout) :: soil
  real,                 intent(in)  :: vlc(:)
  real,                 intent(in)  :: vsc(:)
  real,                 intent(out) :: soil_E_max
  real,                 intent(out) :: thermal_cond(:)
  real s, w, a, n, f

  integer l

  ! assign some index of water availability for snow-free soil

  soil_E_max = (soil%pars%k_sat_ref*soil%alpha(1)**2) &
               * (-soil%pars%psi_sat_ref/soil%alpha(1)) &
               * ((4.+soil%pars%chb)*vlc(1)/ &
                ((3.+soil%pars%chb)*soil%pars%vwc_sat))**(3.+soil%pars%chb) &
                / ((1.+3./soil%pars%chb)*dz(1))

     w = soil%pars%thermal_cond_weight
     a = soil%pars%thermal_cond_scale
     n = soil%pars%thermal_cond_exp
  do l = 1, num_sfc_layers
     soil%heat_capacity_dry(l) = sfc_heat_factor*soil%pars%heat_capacity_dry
     s = (vlc(l)+vsc(l))/soil%pars%vwc_sat
     thermal_cond(l)      = sfc_heat_factor * &
          ( soil%pars%thermal_cond_dry+ &
            (soil%pars%thermal_cond_sat-soil%pars%thermal_cond_dry) &
            *(w*s +(1-w)*(1+a**n)*(s**n)/(1+(a*s)**n))    )
     f = 1.
     if (vlc(l)+vsc(l).gt.0.) f = 1.+(freeze_factor-1.)*vsc(l)/(vlc(l)+vsc(l))
     thermal_cond(l) = f * thermal_cond(l)
  enddo
  do l = num_sfc_layers+1, num_l
     soil%heat_capacity_dry(l) = soil%pars%heat_capacity_dry
     s = (vlc(l)+vsc(l))/soil%pars%vwc_sat
     thermal_cond(l)  = &
          ( soil%pars%thermal_cond_dry+ &
            (soil%pars%thermal_cond_sat-soil%pars%thermal_cond_dry) &
            *(w*s +(1-w)*(1+a**n)*(s**n)/(1+(a*s)**n))    )
     f = 1.
     if (vlc(l)+vsc(l).gt.0.) f = 1.+(freeze_factor-1.)*vsc(l)/(vlc(l)+vsc(l))
     thermal_cond(l) = f * thermal_cond(l)
  enddo

  ! this is an additional factor intended for tuning annual T range in
  ! high latitudes. presumably other locations are insensitive to this
  ! global parameter, since they don't have freeze/thaw. this really is just a fudge.
  do l = sub_layer_min, sub_layer_max
     thermal_cond(l) = sub_layer_tc_fac * thermal_cond(l)
  enddo

end subroutine soil_data_thermodynamics


! ============================================================================
! compute soil hydraulic properties: wrapper to get all parameters for
! richards equation for full column. note that psi_for_rh is not passed.
subroutine soil_data_hydraulic_properties (soil, vlc, vsc, &
                    psi, DThDP, K_z, K_x, DKDP, DPsi_min, DPsi_max  )
  type(soil_tile_type),        intent(inout) :: soil
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: &
      psi, DThDP, K_z, K_x, DKDP
  real,                        intent(out) :: &
      DPsi_min, DPsi_max
  ! ---- local vars ----------------------------------------------------------
  real ::  psi_for_rh

  if (use_alt3_soil_hydraulics) then
      call soil_data_hydraulics_alt3 (soil, vlc, vsc, &
                    psi, DThDP, K_z, K_x, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )
  else
      call soil_data_hydraulics (soil, vlc, vsc, &
                    psi, DThDP, K_z, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )
      K_x = K_z
  endif

  ! For use in hillslope model
  soil%hyd_cond_horz(1:num_l) = K_x(1:num_l)
  where (soil%ws > 0. .or. soil%wl <= 0.) soil%hyd_cond_horz = epsln
  if (all(DThDP==0.)) soil%hyd_cond_horz(:) = epsln ! Will be "stiff"

end subroutine soil_data_hydraulic_properties


! ============================================================================
! wrapper to compute psi and the potentially different psi used for surface
! relative humidity computation.
! use of this wrapper allows us to avoid passing all the other properties
! back to soil module. (they're computed anyway.)
subroutine soil_data_psi_for_rh (soil, vlc, vsc, psi, psi_for_rh  )
  type(soil_tile_type),        intent(in) :: soil
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: psi
  real,                        intent(out) :: psi_for_rh
  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l) ::  DThDP, K_z, K_x, DKDP
  real :: DPsi_min, DPsi_max

  if (use_alt3_soil_hydraulics) then
      call soil_data_hydraulics_alt3 (soil, vlc, vsc, &
                    psi, DThDP, K_z, K_x, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )
  else
      call soil_data_hydraulics (soil, vlc, vsc, &
                    psi, DThDP, K_z, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )
  endif

end subroutine soil_data_psi_for_rh


! ============================================================================
! compute soil hydraulic properties.
subroutine soil_data_hydraulics (soil, vlc, vsc, &
                    psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )
  type(soil_tile_type),        intent(in) :: soil
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: &
      psi, DThDP, hyd_cond, DKDP
  real,                        intent(out) :: &
      DPsi_min, DPsi_max, psi_for_rh
  ! ---- local vars ----------------------------------------------------------
  integer l
  real :: vlc_loc, k_sat, alt_psi_for_rh
  real :: alpha_sq, f_psi
  logical flag

  ! ---- T-dependence of hydraulic properties --------------------------------
  ! k_sat   = soil%pars%k_sat0   !  * mu(t0)/mu(t), where mu is dynamic viscosity
  ! psi_sat = soil%pars%psi_sat0 !  * exp(c*(psi-psi0)), where c~+/-(?)0.0068
                     ! better approach would be to adopt air entrapment model
                     ! or at least to scale against surface tension model
  ! ---- water and ice dependence of hydraulic properties --------------------
  ! ---- (T-dependence can be added later)
  hyd_cond=1;DThDP=1;psi=1
  DKDP=0
  flag = .false.
  do l = 1, size(vlc)
     alpha_sq = soil%alpha(l)**2
     hyd_cond(l) = (soil%pars%k_sat_ref*alpha_sq)*  &
                (vlc(l)/soil%pars%vwc_sat)**(3+2*soil%pars%chb)
     if (hyd_cond(l).lt.1.e-12*soil%pars%k_sat_ref*alpha_sq) then
        vlc_loc     = soil%pars%vwc_sat*(1.e-12)**(1./(3+2*soil%pars%chb))
        hyd_cond(l) = 1.e-12*soil%pars%k_sat_ref*alpha_sq
        if (l.eq.1) flag = .true.
        if (vsc(l).eq.0.) then
           DThDP   (l) = -vlc_loc     &
                     *(vlc_loc/soil%pars%vwc_sat)**soil%pars%chb &
                     /(soil%pars%psi_sat_ref*soil%pars%chb/soil%alpha(l))
           psi     (l) = (soil%pars%psi_sat_ref/soil%alpha(l)) &
                     *(soil%pars%vwc_sat/vlc_loc)**soil%pars%chb &
                     + (vlc(l)-vlc_loc)/DThDP   (l)
           if (l.eq.1.and.vlc(1).gt.0.) then
              alt_psi_for_rh = &
                    (soil%pars%psi_sat_ref/soil%alpha(l)) &
                    *(soil%pars%vwc_sat/vlc(1)   )**soil%pars%chb
           else if (l.eq.1.and.vlc(1).le.0.) then
              alt_psi_for_rh = -1.e10
           endif
        else
           psi     (l) = ((soil%pars%psi_sat_ref/soil%alpha(l)) / 2.2) &
                   *(soil%pars%vwc_sat/vlc_loc   )**soil%pars%chb
           DThDP   (l) = 0.
           if (l.eq.1) alt_psi_for_rh = -1.e10
        endif
     else
        if (vsc(l).eq.0.) then
           if (vlc(l).le.soil%pars%vwc_sat) then
              psi     (l) = (soil%pars%psi_sat_ref/soil%alpha(l)) &
                      *(soil%pars%vwc_sat/vlc(l))**soil%pars%chb
              DKDP    (l) = -(2+3/soil%pars%chb)*hyd_cond(l) &
                                             /psi(l)
              DThDP   (l) = -vlc(l)/(psi(l)*soil%pars%chb)
           else
              psi(l) = soil%pars%psi_sat_ref/soil%alpha(l) &
                     + (vlc(l)-soil%pars%vwc_sat)/comp
              f_psi = min(max(1.-psi(l)/soil%pars%psi_sat_ref/soil%alpha(l),0.),1.)
              DThDP(l) = comp
              hyd_cond(l) = soil%pars%k_sat_ref*alpha_sq &
                               + f_psi * soil%k_macro_z(l)
           endif
        else
           psi     (l) = ((soil%pars%psi_sat_ref/soil%alpha(l)) / 2.2) &
                     *(soil%pars%vwc_sat/vlc(l))**soil%pars%chb
           DThDP   (l) = 0.
        endif
     endif
  enddo

  if (use_alt_psi_for_rh .and. flag) then
     psi_for_rh = alt_psi_for_rh
  else
     psi_for_rh = psi(1)
  endif

  if (DThDP(1).ne.0.) then
     DPsi_min =            -vlc(1) /DThDP(1)
     DPsi_max = (soil%pars%vwc_sat-vlc(1))/DThDP(1)
  else
     Dpsi_min = Dpsi_min_const
     DPsi_max = -psi(1)
  endif

end subroutine soil_data_hydraulics

! ============================================================================
! compute soil hydraulic properties.
subroutine soil_data_hydraulics_alt3 (soil, vlc, vsc, &
                    psi, DThDP, K_z, K_x, DKDP, DPsi_min, DPsi_max, &
                    psi_for_rh  )
  type(soil_tile_type),        intent(in) :: soil
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: &
      psi, DThDP, K_z, K_x, DKDP
  real,                        intent(out) :: &
      DPsi_min, DPsi_max, psi_for_rh
  ! ---- local vars ----------------------------------------------------------
  integer l
  real :: Ksat, Psat, Xsat, B, Xl, Xs, f_psi, X_min, Xl_eff
  real :: psimax ! [m] maximum feasible physical psi, above which comp increases to prevent
                 ! excessive pressures and smooth numerics
  real :: Xmax   ! [-] theta associated with psimax
  real, parameter :: supercomp = 0.01 ! [1/m] comp to use above psimax

  K_z=1;K_x=1;DThDP=1;psi=1
  DKDP=0
  do l = 1, size(vlc)
    Xl = vlc(l)
    Xs = vsc(l)
    Xl_eff = Xl
    if (use_fluid_ice) Xl_eff = Xl + Xs
    Ksat = soil%pars%k_sat_ref    * soil%alpha(l)**2
    Psat = soil%pars%psi_sat_ref  / soil%alpha(l)
    if (Xs.gt.0.) Psat = Psat / 2.2
    Xsat = soil%pars%vwc_sat
    B    = soil%pars%chb
    X_min = Xsat/(psi_min/Psat)**(1/B)
    if (Xl_eff.lt.X_min) then
        DThDP(l) = -X_min/(psi_min*B)
        if (limit_DThDP) DThDP(l) = min(DThDP(l), DThDP_max)
        psi(l) = psi_min + (Xl_eff-X_min)/DThDP(l)
        if (Xs.gt.0. .and. .not.use_fluid_ice) DThDP(l) = 0.
    else if (Xl_eff.le.Xsat) then
        psi(l) = Psat*(Xsat/Xl_eff)**B
        DThDP(l) = -Xl_eff/(psi(l)*B)
        if (limit_DThDP) DThDP(l) = min(DThDP(l), DThDP_max)
        if (Xs.gt.0. .and. .not.use_fluid_ice) DThDP(l) = 0.
    else
       if (comp.gt.0.) then
          psi(l) = Psat + (Xl_eff-Xsat)/comp
          if (Xs.le.0. .or. use_fluid_ice) then
              DThDP(l) = comp
          else
              DThDP(l) = 0.
          endif
          ! ZMS Adding code for stabilization of numerics
          psimax = soil%pars%hillslope_relief + zhalf(l+1)
          if (gw_option == GW_TILED) psimax = psimax - soil%pars%tile_hlsp_elev
          if (limit_hi_psi .and. psi(l) > psimax) then
             Xmax = comp*(psimax - Psat) + Xsat
             psi(l) = psimax + (Xsat - Xmax)/supercomp
             if (Xs.le.0. .or. use_fluid_ice) then
                DThDP(l) = supercomp
             else
                DThDP(l) = 0.
             end if
          endif
          ! end of ZMS addition
       else
           psi(l) = Psat
           DThDP(l) = -Xsat/(Psat*B)
       endif
    endif
    if (use_fluid_ice) then
        K_z(l) = Ksat*((Xl+Xs)/Xsat)**(3+2*B)
    else
        K_z(l) = Ksat*((  Xl )/Xsat)**(3+2*B)
    endif
    K_z(l) = min(max(K_min, K_z(l)), Ksat)
    K_z(l) = min(K_z(l), K_max_matrix)
    K_x(l) = K_z(l)
    f_psi = min(max(1.-psi(l)/Psat,0.),1.)
    K_z(l) = K_z(l) + f_psi * soil%k_macro_z(l)
    K_x(l) = K_x(l) + f_psi * soil%k_macro_x(l)
  enddo

  psi_for_rh = psi(1)

  if (DThDP(1).ne.0.) then
      Xl_eff = vlc(1)
      if (use_fluid_ice) Xl_eff = vlc(1) + vsc(1)
      DPsi_min =      -Xl_eff /DThDP(1)
      DPsi_max = (Xsat-Xl_eff)/DThDP(1)
      ! NOTE THAT TOTAL WATER CONTENT SHOULD BE SUBTRACTED TO GET
      ! DPSI_MAX, BUT IF VSC>0, THEN DThDp IS ZERO, AND WE WOULD
      ! BE IN THE OTHER BRANCH (BELOW).
  else
      Dpsi_min = Dpsi_min_const
      DPsi_max = -psi(1)
  endif

end subroutine soil_data_hydraulics_alt3

! ============================================================================
subroutine soil_data_gw_hydraulics_ar5(soil, storage_normalized, &
                                                     gw_flux, sat_frac )
  type(soil_tile_type), intent(inout)  :: soil
  real,                 intent(in)  :: storage_normalized
  real,                 intent(out) :: gw_flux
  real,                 intent(out) :: sat_frac

  integer :: code, m
  real :: recharge_normalized, frac

  ! storage_normalized is the fraction of soil above drainage base elevation
  ! that is below the water table
  code = 0
  m = soil%pars%storage_index
  do while (code.eq.0)
     if (storage_normalized .lt. gw_storage_norm(m)) then
        m = m - 1
     else if (storage_normalized .gt. gw_storage_norm(m+1)) then
        m = m + 1
     else
        code = 1
     endif
  enddo
  if (m.lt.1.or.m.gt.num_storage_pts-1) then
     write(*,*) '!!! *** m=',m, ' is outside the table in soil_data_gw_hydraulics_ar5 *** !!!'
     write(*,*) 'num_storage_pts=',num_storage_pts
     write(*,*) 'storage_normalized=',storage_normalized
     write(*,*) 'interval bounds:',gw_storage_norm(m),gw_storage_norm(m+1)
  endif
  frac = (storage_normalized-gw_storage_norm(m)) &
           /(gw_storage_norm(m+1)-gw_storage_norm(m))
  sat_frac = soil%gw_area_norm(m) &
               + frac*(soil%gw_area_norm(m+1)-soil%gw_area_norm(m))
  recharge_normalized = soil%gw_flux_norm(m) &
               + frac*(soil%gw_flux_norm(m+1)-soil%gw_flux_norm(m))
  gw_flux = recharge_normalized * soil%pars%k_sat_ref * soil%pars%soil_e_depth &
                * soil%pars%hillslope_relief &
                   / (soil%pars%hillslope_length * soil%pars%hillslope_length)
  soil%pars%storage_index = m

end subroutine soil_data_gw_hydraulics_ar5

! ============================================================================
subroutine soil_data_gw_hydraulics(soil, deficit_normalized, &
                                                     gw_flux, sat_frac )
  type(soil_tile_type), intent(inout)  :: soil
  real,                 intent(in)  :: deficit_normalized
  real,                 intent(out) :: gw_flux
  real,                 intent(out) :: sat_frac

  integer :: code, m, a_integer, n_integer
  real :: recharge_normalized, frac, log_deficit_normalized, xi_s
  real :: k_macro_x_local

  ! deficit_normalized is the fraction of soil above drainage base elevation
  ! that is NOT saturated/contributing to horizontal flow

  if (soil%tag == peat_soil_type .and. peat_kx0 >= 0.) then
     k_macro_x_local = peat_kx0
  else
     k_macro_x_local = k0_macro_x
  end if

  IF (deficit_normalized .LE. 0.) THEN
     m = 0
  ELSE IF (deficit_normalized .GE. 1.) THEN
     m = num_storage_pts
  ELSE
     log_deficit_normalized = log10(deficit_normalized)
     code = 0
     m = soil%pars%storage_index
     do while (code.eq.0.and.m.gt.0.and.m.lt.num_storage_pts)
        if (log_deficit_normalized .lt. log_deficit_list(m)) then
           m = m - 1
        else if (log_deficit_normalized .gt. log_deficit_list(m+1)) then
           m = m + 1
        else
           code = 1
        endif
     enddo
  ENDIF

  IF (m.eq.0) THEN
     recharge_normalized = 10.**soil%gw_flux_norm(1)
  ELSE IF (m.lt.num_storage_pts) then
     frac = (log_deficit_normalized-log_deficit_list(m)) &
              /(log_deficit_list(m+1)-log_deficit_list(m))
     recharge_normalized = 10.**(soil%gw_flux_norm(m) &
                  + frac*(soil%gw_flux_norm(m+1)-soil%gw_flux_norm(m)))
  ELSE
     recharge_normalized = 0.
  ENDIF

  if (use_tau_fix) then
     gw_flux = recharge_normalized * &
        ((soil%pars%k_sat_sfc+k_macro_x_local) * soil%pars%soil_e_depth + &
                  soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length) &
        * soil%pars%hillslope_relief &
        / (soil%pars%hillslope_length * soil%pars%hillslope_length)
  else
     gw_flux = recharge_normalized * &
        ((soil%pars%k_sat_sfc+k_macro_x_local) * soil%pars%soil_e_depth) &
        * soil%pars%hillslope_relief &
        / (soil%pars%hillslope_length * soil%pars%hillslope_length)
  endif

  a_integer = int(soil%pars%hillslope_a + 0.1)
  n_integer = int(soil%pars%hillslope_n + 0.1)

  if (a_integer.eq.0.and.n_integer.eq.1) then
     if (recharge_normalized.gt.1.) then
         sat_frac = 1. - 1./recharge_normalized
     else
         sat_frac = 0.
     endif
     if (use_sat_fix) then
         if (recharge_normalized.gt.1.+soil%pars%tau) then
             sat_frac = 1. - (1.+soil%pars%tau)/recharge_normalized
         else
             sat_frac = 0.
         endif
     endif
  endif

  if (a_integer.eq.1.and.n_integer.eq.1) then
     xi_s = -1.-(1.+soil%pars%tau)/recharge_normalized &
     + sqrt((1.+soil%pars%tau)*(1.+soil%pars%tau) &
             +4.*recharge_normalized*recharge_normalized) &
               /recharge_normalized
     xi_s = max(xi_s, 0.)
     sat_frac = (xi_s+0.5*xi_s*xi_s)/1.5
  endif

  if (a_integer.eq.0.and.n_integer.eq.2) then
     sat_frac = 1./(1.+2.*(1.+soil%pars%tau)/recharge_normalized)
  endif

  if (a_integer.eq.1.and.n_integer.eq.2) then
     xi_s = -(1.+2.*(1.+soil%pars%tau)/recharge_normalized)
     xi_s = xi_s + &
       sqrt(xi_s*xi_s+6.*(0.5+2.*(1.+soil%pars%tau)/recharge_normalized))
       xi_s = xi_s / (1.+4.*(1.+soil%pars%tau)/recharge_normalized)
     sat_frac = (xi_s+0.5*xi_s*xi_s)/1.5
  endif

  soil%pars%storage_index = min(max(1,m),num_storage_pts-1)

end subroutine soil_data_gw_hydraulics

! ============================================================================
subroutine soil_data_vwc_for_init_only (soil, psi, vwc)
  type(soil_tile_type), intent(in) :: soil
  real,                 intent(in) :: psi(:)
  real,                 intent(out):: vwc(:)
  vwc = soil%pars%vwc_sat* &
    ((soil%pars%psi_sat_ref/soil%alpha)/min(psi,(soil%pars%psi_sat_ref/soil%alpha))) &
     ** (1./ soil%pars%chb)
  if (use_comp_for_ic) then
      vwc = vwc + comp*max(psi-soil%pars%psi_sat_ref/soil%alpha,0.)
  endif
end subroutine soil_data_vwc_for_init_only

! ============================================================================
subroutine soil_tile_stock_pe (soil, twd_liq, twd_sol  )
  type(soil_tile_type),  intent(in)    :: soil
  real,                  intent(out)   :: twd_liq, twd_sol
  integer n

  twd_liq = 0.
  twd_sol = 0.
  do n=1, size(soil%wl)
    twd_liq = twd_liq + soil%wl(n) + soil%groundwater(n)
    twd_sol = twd_sol + soil%ws(n)
  enddo

end subroutine soil_tile_stock_pe


! ============================================================================
! returns soil tile heat content, J/m2
function soil_tile_heat (soil) result(heat) ; real heat
  type(soil_tile_type),  intent(in)  :: soil

  integer :: i

  heat = 0
  do i = 1, num_l
     heat = heat + &
          (soil%heat_capacity_dry(i)*dz(i)+clw*soil%Wl(i)+csw*soil%Ws(i))&
                           *(soil%T(i)-tfreeze) + &
          clw*soil%groundwater(i)*(soil%groundwater_T(i)-tfreeze) - &
          hlf*soil%ws(i)
  enddo
end function soil_tile_heat

! ============================================================================
! returns soil tile carbon content, kg C/m2
function soil_tile_carbon (soil); real soil_tile_carbon
  type(soil_tile_type),  intent(in)  :: soil

  real    :: temp
  integer :: i

  select case (soil_carbon_option)
  case (SOILC_CORPSE,SOILC_CORPSE_N)
     soil_tile_carbon = 0.0
     do i=1,num_l
	call poolTotals(soil%soil_organic_matter(i),totalCarbon=temp)
        soil_tile_carbon=soil_tile_carbon+temp
     enddo
     do i = 1,N_LITTER_POOLS
        call poolTotals(soil%litter(i),totalCarbon=temp)
        soil_tile_carbon=soil_tile_carbon+temp
     enddo
  case default
     soil_tile_carbon = sum(soil%fast_soil_C(:))+sum(soil%slow_soil_C(:))
  end select
end function soil_tile_carbon

! ============================================================================
! returns soil tile nitrogen content, kg N/m2
! this should return zero if soil_carbon_option is not SOILC_CORPSE_N
function soil_tile_nitrogen (soil); real soil_tile_nitrogen
  type(soil_tile_type),  intent(in)  :: soil

  real    :: temp
  integer :: i

  select case (soil_carbon_option)
  case (SOILC_CORPSE,SOILC_CORPSE_N)
     soil_tile_nitrogen = 0.0
     do i=1,num_l
        call poolTotals(soil%soil_organic_matter(i),totalNitrogen=temp)
        soil_tile_nitrogen=soil_tile_nitrogen+temp
     enddo
     do i = 1,N_LITTER_POOLS
        call poolTotals(soil%litter(i),totalNitrogen=temp)
        soil_tile_nitrogen=soil_tile_nitrogen+temp
     enddo
  case default
     soil_tile_nitrogen = 0.0
  end select
end function soil_tile_nitrogen

end module soil_tile_mod
