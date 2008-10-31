#include <fms_platform.h>

module soil_tile_mod

use fms_mod, only : &
     write_version_number, file_exist, open_namelist_file, check_nml_error, &
     close_file, stdlog
use constants_mod, only : &
     pi, tfreeze, rvgas, grav, dens_h2o
use land_constants_mod, only : &
     NBANDS
use land_io_mod, only : &
     init_cover_field
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_SOIL, register_tile_selector
use vegn_tile_mod, only : vegn_tile_type, vegn_data_uptake_frac_max

implicit none
private

! ==== public interfaces =====================================================
public :: soil_pars_type
public :: soil_prog_type
public :: soil_tile_type

public :: new_soil_tile, delete_soil_tile
public :: soil_tiles_can_be_merged, merge_soil_tiles
public :: soil_is_selected
public :: get_soil_tile_tag

public :: read_soil_data_namelist
public :: soil_cover_cold_start

public :: soil_data_beta
public :: soil_data_radiation
public :: soil_data_diffusion
public :: soil_data_thermodynamics
public :: soil_data_hydraulics
public :: soil_data_w_sat
public :: soil_ave_temp  ! calculate average soil temeperature
public :: soil_ave_theta ! calculate average soil moisture

public :: max_lev
! =====end of public interfaces ==============================================
interface new_soil_tile
   module procedure soil_tile_ctor
   module procedure soil_tile_copy_ctor
end interface

! ==== module constants ======================================================
character(len=*), parameter   :: &
     version     = '$Id: soil_tile.F90,v 16.0 2008/07/30 22:30:01 fms Exp $', &
     tagname     = '$Name: perth_2008_10 $', &
     module_name = 'soil_tile_mod'

integer, parameter :: max_lev          = 30 
integer, parameter :: n_dim_soil_types = 9       ! size of lookup table
real,    parameter :: psi_wilt         = -150.0  ! matric head at wilting
real,    parameter :: small            = 1.e-4
real,    parameter :: t_ref            = 293
real,    parameter :: g_RT             = grav / (rvgas*t_ref)
real,    parameter :: comp             = 0.001  ! m^-1

! ==== types =================================================================
type :: soil_pars_type
  real w_sat
  real awc_lm2
  real k_sat_ref
  real psi_sat_ref
  real chb
  real alpha              ! *** REPLACE LATER BY alpha(layer)
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
  real emis_dry
  real emis_sat
  real z0_momentum
  real tau_groundwater
  real rsa_exp         ! riparian source-area exponent
  real tfreeze
end type soil_pars_type


type :: soil_prog_type
  real wl
  real ws
  real T
  real groundwater
  real groundwater_T
end type soil_prog_type


type :: soil_tile_type
   integer :: tag ! kind of the soil
   type(soil_pars_type)               :: pars
   type(soil_prog_type), pointer :: prog(:)
   real,                 pointer :: w_fc(:)
   real,                 pointer :: w_wilt(:)
   real,                 pointer :: uptake_frac_max(:)
   real :: Eg_part_ref
   real :: z0_scalar
   ! data that were local to soil.f90
   real,                 pointer :: uptake_frac(:)
   real,                 pointer :: heat_capacity_dry(:)
   real,                 pointer :: e(:),f(:)
end type soil_tile_type

! ==== module data ===========================================================
real, public :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

!---- namelist ---------------------------------------------------------------
real    :: k_over_B              = 2         ! reset to 0 for MCM
real    :: rate_fc               = 0.1/86400 ! 0.1 mm/d drainage rate at FC
real    :: sfc_heat_factor       = 1
real    :: z_sfc_layer           = 0.0
integer :: num_l                 = 18        ! number of soil levels
real    :: dz(max_lev)           = (/ &
    0.02, 0.04, 0.04, 0.05, 0.05, 0.1, 0.1, 0.2, 0.2, &
    0.2,   0.4,  0.4,  0.4,  0.4, 0.4,  1.,  1.,  1., &
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)
                                              ! thickness (m) of model layers,
                                              ! from top down
logical :: use_lm2_awc           = .false.
logical :: lm2                   = .false.
! ---- remainder are used only for cold start ---------
character(len=16):: soil_to_use     = 'single-tile'
       ! 'multi-tile' for multiple soil types per grid cell, a tile per type
       ! 'single-tile' for geographically varying soil with single type per
       !     model grid cell [default]
       ! 'uniform' for global constant soil, e.g., to reproduce MCM

logical :: use_mcm_albedo        = .false.   ! .true. for CLIMAP albedo inputs
logical :: use_single_geo        = .false.   ! .true. for global gw res time,
                                             ! e.g., to recover MCM
integer :: soil_index_constant   = 9         ! index of global constant soil,
                                             ! used when use_single_soil
real    :: gw_res_time           = 60.*86400 ! mean groundwater residence time,
                                             ! used when use_single_geo
real    :: rsa_exp_global        = 1.5
  real, dimension(n_dim_soil_types) :: &
  dat_w_sat=&
  (/ 0.380, 0.445, 0.448, 0.412, 0.414, 0.446, 0.424, 0.445, 0.445   /),&
  dat_awc_lm2=&
  (/ 0.063, 0.132, 0.109, 0.098, 0.086, 0.120, 0.101, 0.445, 0.150   /),&
  dat_k_sat_ref=&
  (/ 0.021, .0036, .0018, .0087, .0061, .0026, .0051, .0036, .0036   /),&
  dat_psi_sat_ref=&
  (/ -.059, -0.28, -0.27, -0.13, -0.13, -0.27, -0.16, -0.28, -0.28   /),&
  dat_chb=&
  (/   3.5,   6.4,  11.0,   4.8,   6.3,   8.4,   6.3,   6.4,   6.4   /),&
!  dat_heat_capacity_ref =&
!  (/ 1.8e6, 2.0e6, 2.6e6, 1.9e6, 2.2e6, 2.3e6, 2.1e6, 3.0e6,   1.0   /),&
! previous (ref) values were based on typical water contents
! following dry values are based on w_min=(1-w_sat) w_org=0
! except for peat, where            w_org-(1-w_sat) w_min=0
! microscopic rho*c for w_min is 2650*733 and for w_org is 1300*1926
! (brutsaert 1982 evaporation into the atmosphere p 146)
! ignored air
  dat_heat_capacity_dry =&
  (/ 1.2e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.4e6,   1.0   /),&
!  dat_thermal_cond_ref =&
!  (/   1.5,   0.8,  1.35,  1.15, 1.475, 1.075, 1.217,  0.39, 2.e-7   /),&
! previous (ref) values were based on typical water contents
! following dry and sat values and interpolating exponents are based on
! computations after deVries. i computed C M and F functions for
! unfrozen soil in
! spreadsheet Research\LaD2\soil thermal conductivity. the dry and
! sat values come right out of those computations. the exponent was
! done by eye. curves look like typical literature curves.
! TEMP: still need to treat freezing, maybe import deVries model into code.
  dat_thermal_cond_dry =&
  (/  0.14,  0.21,  0.20,  .175, 0.170, 0.205, 0.183,  0.05, 2.e-7   /),&
  dat_thermal_cond_sat =&
  (/  2.30,  1.50,  1.50,  1.90, 1.900, 1.500, 1.767,  0.50, 2.e-7   /),&
  dat_thermal_cond_scale =&
  (/  15.0,  0.50,   10.,  2.74,  12.2,  2.24,  4.22,   1.0,   1.0   /),&
  dat_thermal_cond_exp =&
  (/   3.0,   5.0,   6.0,   4.0,   4.5,   5.5, 4.667,   1.0,   1.0   /),&
  dat_thermal_cond_weight =&
  (/  0.20,  0.70,   0.7,  0.45, 0.450, 0.700, 0.533,   1.0,   1.0   /),&
  dat_emis_dry=&
  (/ 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950,   1.0   /),&
  dat_emis_sat=&
  (/ 0.980, 0.975, 0.970, .9775, 0.975, .9725, 0.975, 0.975,   1.0   /),&
  dat_z0_momentum=&
  (/  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01, 0.045   /),&
  dat_tf_depr=&
  (/  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,   0.0   /)
  !Coarse  Medium   Fine    CM     CF     MF    CMF    Peat    MCM
real, dimension(n_dim_soil_types,NBANDS) :: &
  dat_refl_dry_dir=&
  (/ 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0     ,& ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0   /),& ! NIR
  dat_refl_dry_dif=&
  (/ 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0     ,& ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0   /),& ! NIR
  dat_refl_sat_dir=&
  (/ 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0     ,& ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0   /),& ! NIR
  dat_refl_sat_dif=&
  (/ 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0     ,& ! visible
     0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 0.333, 999.0   /)   ! NIR
  !Coarse  Medium   Fine    CM     CF     MF    CMF    Peat    MCM
integer, dimension(n_dim_soil_types) :: &
  input_cover_types=&
  (/ 1,     2,     3,     4,     5,     6,     7,     8,     100   /)
character(len=4), dimension(n_dim_soil_types) :: &
  tile_names=&
  (/'c   ','m   ','f   ','cm  ','cf  ','mf  ','cmf ','peat','mcm ' /)

namelist /soil_data_nml/ &
     soil_to_use, tile_names, input_cover_types, &
     k_over_B,             &
     rate_fc, sfc_heat_factor, z_sfc_layer,          &
     num_l,                   dz,                      &
     use_lm2_awc,    lm2, &
     use_mcm_albedo,            &
     use_single_geo,            soil_index_constant,         &
     gw_res_time,            rsa_exp_global,      &
     dat_w_sat,               dat_awc_lm2,     &
     dat_k_sat_ref,            &
     dat_psi_sat_ref,               dat_chb,          &
     dat_heat_capacity_dry,       dat_thermal_cond_dry,   &
     dat_thermal_cond_sat,        dat_thermal_cond_exp,   &
     dat_thermal_cond_scale,        dat_thermal_cond_weight,   &
     dat_refl_dry_dir,            dat_refl_sat_dir,              &
     dat_refl_dry_dif,            dat_refl_sat_dif,              &
     dat_emis_dry,              dat_emis_sat,                &
     dat_z0_momentum,           dat_tf_depr
!---- end of namelist --------------------------------------------------------

integer :: num_sfc_layers

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_soil_data_namelist(soil_num_l, soil_dz, soil_single_geo)
  integer, intent(out) :: soil_num_l
  real,    intent(out) :: soil_dz(:)
  logical, intent(out) :: soil_single_geo
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  real    :: z

  call write_version_number(version, tagname)
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
  write (stdlog(), nml=soil_data_nml)
  
  ! register selector for all soil tiles
  call register_tile_selector('soil', long_name='soil',&
       tag = SEL_SOIL, idata1 = 0 )
  ! register selectors for tile-specific diagnostics
  do i=1, n_dim_soil_types
     call register_tile_selector(tile_names(i), long_name='',&
          tag = SEL_SOIL, idata1 = i )
  enddo
  z = 0
  num_sfc_layers = 0

  do i = 1, num_l
    z = z + dz(i)
    if (z < z_sfc_layer+1.e-4) num_sfc_layers = i
  enddo

  ! set up output arguments
  soil_num_l      = num_l
  soil_dz         = dz
  soil_single_geo = use_single_geo

end subroutine 


! ============================================================================
function soil_tile_ctor(tag) result(ptr)
  type(soil_tile_type), pointer :: ptr ! return value
  integer, intent(in)  :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = tag
  ! allocate storage for tile data
  allocate( ptr%prog(num_l))
  allocate( ptr%w_fc              (num_l),  &
            ptr%w_wilt            (num_l),  &
            ptr%uptake_frac_max   (num_l),  &
            ptr%uptake_frac       (num_l),  &
            ptr%heat_capacity_dry (num_l),  &
            ptr%e                 (num_l),  &
            ptr%f                 (num_l)   )

  call soil_data_init_0d(ptr)
end function soil_tile_ctor


! ============================================================================
function soil_tile_copy_ctor(soil) result(ptr)
  type(soil_tile_type), pointer :: ptr ! return value
  type(soil_tile_type), intent(in) :: soil ! tile to copy

  allocate(ptr)
  ptr = soil ! copy all non-pointer members
  ! allocate storage for tile data
  allocate( ptr%prog(num_l))
  allocate( ptr%w_fc              (num_l),  &
            ptr%w_wilt            (num_l),  &
            ptr%uptake_frac_max   (num_l),  &
            ptr%uptake_frac       (num_l),  &
            ptr%heat_capacity_dry (num_l),  &
            ptr%e                 (num_l),  &
            ptr%f                 (num_l)   )
  ! copy all pointer members
  ptr%prog(:) = soil%prog(:)
  ptr%w_fc(:) = soil%w_fc(:)
  ptr%w_wilt(:) = soil%w_wilt(:)
  ptr%uptake_frac_max(:) = soil%uptake_frac_max(:)
  ptr%uptake_frac(:) = soil%uptake_frac(:)
  ptr%heat_capacity_dry(:) = soil%heat_capacity_dry(:)
  ptr%e(:) = soil%e(:)
  ptr%f(:) = soil%f(:)
end function soil_tile_copy_ctor


! ============================================================================
subroutine delete_soil_tile(ptr)
  type(soil_tile_type), pointer :: ptr

  deallocate(ptr%prog)
  deallocate(ptr%w_fc, ptr%w_wilt, ptr%uptake_frac_max, ptr%uptake_frac,&
             ptr%heat_capacity_dry, ptr%e, ptr%f)
  deallocate(ptr)
end subroutine delete_soil_tile


! ============================================================================
subroutine soil_data_init_0d(soil)
  type(soil_tile_type), intent(inout) :: soil
  
!  real tau_groundwater
!  real rsa_exp         ! riparian source-area exponent
  integer :: k 
  k = soil%tag

  soil%pars%w_sat             = dat_w_sat            (k)
  soil%pars%awc_lm2           = dat_awc_lm2          (k)
  soil%pars%k_sat_ref         = dat_k_sat_ref        (k)
  soil%pars%psi_sat_ref       = dat_psi_sat_ref      (k)
  soil%pars%chb               = dat_chb              (k)
  soil%pars%alpha             = 1
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
!   soil_type_tfreeze(:,:,iii)       = tfreeze - dat_tf_depr  (jjj)

  soil%pars%rsa_exp           = rsa_exp_global
  soil%pars%tau_groundwater   = gw_res_time

  ! ---- derived constant soil parameters
  ! w_fc (field capacity) set to w at which hydraulic conductivity equals
  ! a nominal drainage rate "rate_fc"
  ! w_wilt set to w at which psi is psi_wilt
  if (use_lm2_awc) then
     soil%w_wilt(:) = 0.15
     soil%w_fc  (:) = 0.15 + soil%pars%awc_lm2
  else
     soil%w_wilt(:) = soil%pars%w_sat &
          *(soil%pars%psi_sat_ref/(psi_wilt*soil%pars%alpha))**(1/soil%pars%chb)
     soil%w_fc  (:) = soil%pars%w_sat &
          *(rate_fc/(soil%pars%k_sat_ref*soil%pars%alpha**2))**(1/(3+2*soil%pars%chb))
  endif

  ! below made use of phi_e from parlange via entekhabi
  soil%Eg_part_ref = (-4*soil%w_fc(1)**2*soil%pars%k_sat_ref*soil%pars%psi_sat_ref*soil%pars%chb &
       /(pi*soil%pars%w_sat)) * (soil%w_fc(1)/soil%pars%w_sat)**(2+soil%pars%chb)   &
       *(2*pi/(3*soil%pars%chb**2*(1+3/soil%pars%chb)*(1+4/soil%pars%chb)))/2

  soil%z0_scalar = soil%pars%z0_momentum * exp(-k_over_B)
end subroutine 

! ============================================================================
function soil_cover_cold_start(land_mask, glonb, glatb) result (soil_frac)
! creates and initializes a field of fractional soil coverage
! should be called for global grid; othervise the parth that fills
! missing data points may fail to find any good data, or do it in nproc-
! dependent way
  logical, intent(in) :: land_mask(:,:)    ! global land mask
  real,    intent(in) :: glonb(:,:), glatb(:,:)! boundaries of the global grid cells
  real,    pointer    :: soil_frac (:,:,:) ! output-global map of soil fractional coverage

  allocate( soil_frac(size(land_mask,1),size(land_mask,2),n_dim_soil_types))

  call init_cover_field(soil_to_use, 'INPUT/ground_type.nc', 'cover','frac', &
       glonb, glatb, soil_index_constant, input_cover_types, soil_frac)
  
end function 


! =============================================================================
function soil_tiles_can_be_merged(soil1,soil2) result(response)
  logical :: response
  type(soil_tile_type), intent(in) :: soil1,soil2

  response = (soil1%tag==soil2%tag)
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
          (s1%heat_capacity_dry(i)+clw*s1%prog(i)%Wl+csw*s1%prog(i)%Ws)* &
          (s1%prog(i)%T-tfreeze)
     HEAT2 = &
          (s2%heat_capacity_dry(i)+clw*s2%prog(i)%Wl+csw*s2%prog(i)%Ws)* &
          (s2%prog(i)%T-tfreeze)
     ! merge the amounts of water
     s2%prog(i)%Wl = x1*s1%prog(i)%Wl + x2*s2%prog(i)%Wl
     s2%prog(i)%Ws = x1*s1%prog(i)%Ws + x2*s2%prog(i)%Ws
     ! if the dry heat capacity of merged soil is to be changed, do it here
     ! ...
     ! calculate the merged temperature based on heat content
     s2%prog(i)%T = tfreeze + (x1*HEAT1+x2*HEAT2)/ &
          (s2%heat_capacity_dry(i)+clw*s2%prog(i)%Wl+csw*s2%prog(i)%Ws)

     ! calculate combined groundwater content
     gw = s1%prog(i)%groundwater*x1 + s2%prog(i)%groundwater*x2
     ! calculate combined groundwater temperature
     if (gw/=0) then
        s2%prog(i)%groundwater_T = ( &
             s1%prog(i)%groundwater*x1*(s1%prog(i)%groundwater_T-tfreeze) + &
             s2%prog(i)%groundwater*x2*(s2%prog(i)%groundwater_T-tfreeze)   &
             ) / gw + tfreeze
     else
        s2%prog(i)%groundwater_T = &
             s1%prog(i)%groundwater_T*x1 + s2%prog(i)%groundwater_T*x2
     endif
     s2%prog(i)%groundwater = gw
  enddo
  
end subroutine

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
! compute beta function
! after Manabe (1969), but distributed vertically.
subroutine soil_data_beta ( soil, vegn, soil_beta, soil_water_supply )
  type(soil_tile_type), intent(inout)  :: soil
  type(vegn_tile_type), intent(in)     :: vegn
  real, intent(out) :: soil_beta
  real, intent(out) :: soil_water_supply ! supply of water to roots per unit
             ! active root biomass, kg/m2

  ! ---- local vars
  integer :: l
  real    :: vlc ! volumetric fraction of water in the layer

  call vegn_data_uptake_frac_max (vegn, dz(1:num_l), soil%uptake_frac_max )

  soil%uptake_frac = 0
  do l = 1, num_l
     vlc = max(0.0, soil%prog(l)%wl / (dens_h2o * dz(l)))
     soil%uptake_frac(l) = soil%uptake_frac_max(l) &
          * max(0.0, min(1.0,(vlc-soil%w_wilt(l))/&
               (0.75*(soil%w_fc(l)-soil%w_wilt(l)))))
  enddo
  soil_beta = sum(soil%uptake_frac)
  do l = 1, num_l
     if (soil_beta /= 0) then
          soil%uptake_frac(l) = soil%uptake_frac(l) / soil_beta
     else
          soil%uptake_frac(l) = soil%uptake_frac_max(l)
     endif
  enddo
  if (lm2) soil%uptake_frac = soil%uptake_frac_max

  ! calculate water supply to vegetation per unit active root biomass 
  soil_water_supply = 0
  do l = 1, num_l
     soil_water_supply = soil_water_supply + &
          soil%uptake_frac_max(l)*max(0.0,soil%prog(l)%wl-soil%w_wilt(l)*dens_h2o*dz(l))
  enddo

end subroutine soil_data_beta


! ============================================================================
! compute average soil temperature with a given depth scale
function soil_ave_temp(soil, depth) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)                 :: depth ! averaging depth

  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  real    :: z ! current depth, m
  integer :: k

  A = 0 ; N = 0 ; z = 0
  do k = 1, num_l
     w = dz(k) * exp(-(z+dz(k)/2)/depth)
     A = A + soil%prog(k)%T * w
     N = N + w
     z = z + dz(k)
  enddo
  A = A/N
end function soil_ave_temp


! ============================================================================
! compute average soil moisture with a given depth scale
function soil_ave_theta(soil, depth) result (A) ; real :: A
  type(soil_tile_type), intent(in) :: soil
  real, intent(in)                 :: depth ! averaging depth

  real    :: w ! averaging weight
  real    :: N ! normalizing factor for averaging
  real    :: z ! current depth, m
  integer :: k

  A = 0 ; N = 0 ; z = 0
  do k = 1, num_l
     w = dz(k) * exp(-(z+dz(k)/2)/depth)
     A = A + max(soil%prog(k)%wl/(dens_h2o*dz(k))-soil%w_wilt(k),0.0)/&
          (soil%w_fc(k)-soil%w_wilt(k)) * w
     N = N + w
     z = z + dz(k)
  enddo
  A = A/N
end function soil_ave_theta


! ============================================================================
! compute bare-soil albedo, bare-soil emissivity, bare-soil roughness
! for scalar transport, and beta function
subroutine soil_data_radiation ( soil, soil_alb_dir, soil_alb_dif, soil_emis )
  type(soil_tile_type), intent(in)  :: soil
  real,                 intent(out) :: soil_alb_dir(NBANDS), soil_alb_dif(NBANDS), soil_emis
  ! ---- local vars
  real :: soil_sfc_vlc, blend

  soil_sfc_vlc  = soil%prog(1)%wl/(dens_h2o*dz(1))
  blend         = soil_sfc_vlc/soil%pars%w_sat ! shouldn't it be restricted to [0,1] interval?
  soil_alb_dir  = soil%pars%refl_dry_dir + blend*(soil%pars%refl_sat_dir-soil%pars%refl_dry_dir)
  soil_alb_dif  = soil%pars%refl_dry_dif + blend*(soil%pars%refl_sat_dif-soil%pars%refl_dry_dif)
  soil_emis     = soil%pars%emis_dry     + blend*(soil%pars%emis_sat    -soil%pars%emis_dry    )
end subroutine soil_data_radiation


! ============================================================================
! compute bare-soil albedo, bare-soil emissivity, bare-soil roughness
! for scalar transport, and beta function
subroutine soil_data_diffusion ( soil, soil_z0s, soil_z0m )
  type(soil_tile_type), intent(in)  :: soil
  real,                 intent(out) :: soil_z0s, soil_z0m

  soil_z0s = soil%z0_scalar
  soil_z0m = soil%pars%z0_momentum
end subroutine soil_data_diffusion

! ============================================================================
! compute soil thermodynmamic properties.
subroutine soil_data_thermodynamics ( soil, vlc, vsc, &  
                                             soil_E_max, soil_rh, &
                                             thermal_cond)
  type(soil_tile_type), intent(inout) :: soil
  real,                 intent(in)  :: vlc(:)
  real,                 intent(in)  :: vsc(:)
  real,                 intent(out) :: soil_E_max
  real,                 intent(out) :: soil_rh
  real,                 intent(out) :: thermal_cond(:)
  real s, w, a, n

  integer l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  ! assign some index of water availability for snow-free soil
  soil_E_max = soil%Eg_part_ref / ( max(small, soil%w_fc(1) - vlc(1)) )  ! NEEDS T adj
  if (vlc(1)+vsc(1)>0) then
     soil_rh = exp ( ((soil%pars%psi_sat_ref/soil%pars%alpha) &
                    *(soil%pars%w_sat/(vlc(1)+vsc(1)))**soil%pars%chb)*g_RT )
  else
     soil_rh = 0
  endif
!!$  soil_rh = max(0.,min(1.,1.*(vlc(1)+vsc(1))/soil%pars%w_sat)) !********* TEMP FIX
!!$  soil_rh = max(0.,min(1.,(vlc(1)+vsc(1)-soil%w_wilt(1))/(soil%pars%w_sat-soil%w_wilt(1)))) !********* TEMP FIX

     w = soil%pars%thermal_cond_weight
     a = soil%pars%thermal_cond_scale
     n = soil%pars%thermal_cond_exp
  do l = 1, num_sfc_layers
     soil%heat_capacity_dry(l) = sfc_heat_factor*soil%pars%heat_capacity_dry
     s = (vlc(l)+vsc(l))/soil%pars%w_sat
     thermal_cond(l)      = sfc_heat_factor * &
          ( soil%pars%thermal_cond_dry+ &
            (soil%pars%thermal_cond_sat-soil%pars%thermal_cond_dry) &
            *(w*s +(1-w)*(1+a**n)*(s**n)/(1+(a*s)**n))    )
  enddo
  do l = num_sfc_layers+1, num_l
     soil%heat_capacity_dry(l) = soil%pars%heat_capacity_dry
     s = (vlc(l)+vsc(l))/soil%pars%w_sat
     thermal_cond(l)  = &
          ( soil%pars%thermal_cond_dry+ &
            (soil%pars%thermal_cond_sat-soil%pars%thermal_cond_dry) &
            *(w*s +(1-w)*(1+a**n)*(s**n)/(1+(a*s)**n))    )
  enddo
  
end subroutine soil_data_thermodynamics


! ============================================================================
! compute soil hydraulic properties.
subroutine soil_data_hydraulics (soil, vlc, vsc, &
                    psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max, tau_gw, &
                    soil_w_fc  )
  type(soil_tile_type),        intent(in) :: soil
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: &
      psi, DThDP, hyd_cond, DKDP, soil_w_fc
  real,                        intent(out) :: &
                     DPsi_min, DPsi_max, tau_gw
  ! ---- local vars ----------------------------------------------------------
  integer l
  real,  dimension(num_l) :: vlc_loc
  
  ! ---- T-dependence of hydraulic properties --------------------------------
  ! k_sat   = soil%pars%k_sat0   !  * mu(t0)/mu(t), where mu is dynamic viscosity
  ! psi_sat = soil%pars%psi_sat0 !  * exp(c*(psi-psi0)), where c~+/-(?)0.0068
                     ! better approach would be to adopt air entrapment model
                     ! or at least to scale against surface tension model


  ! ---- water and ice dependence of hydraulic properties --------------------
  ! ---- (T-dependence can be added later)
  hyd_cond=1;DThDP=1;psi=1
  do l = 1, num_l
    hyd_cond(l) = (soil%pars%k_sat_ref*soil%pars%alpha**2)*  &
                ! * mu(T)/mu(t_ref), where mu is dynamic viscosity
               (vlc(l)/soil%pars%w_sat)**(3+2*soil%pars%chb)
    if (hyd_cond(l).lt.1.e-12*soil%pars%k_sat_ref) then
      vlc_loc (l) = soil%pars%w_sat*(1.e-12)**(1./(3+2*soil%pars%chb))
      hyd_cond(l) = 1.e-12*soil%pars%k_sat_ref
      if (vsc(l).eq.0.) then
        DThDP   (l) = -vlc_loc(l)  &
                         *(vlc_loc(l)/soil%pars%w_sat)**soil%pars%chb &
                 /(soil%pars%psi_sat_ref*soil%pars%chb)
        psi     (l) = (soil%pars%psi_sat_ref/soil%pars%alpha) &
            *(soil%pars%w_sat/vlc_loc(l))**soil%pars%chb &
            + (vlc(l)-vlc_loc (l))/DThDP   (l)
        DKDP    (l) = 0.
      else
        psi     (l) = ((soil%pars%psi_sat_ref/soil%pars%alpha) / 2.2) &
            *(soil%pars%w_sat/vlc_loc(l))**soil%pars%chb
        DKDP    (l) = 0.
        DThDP   (l) = 0.
      endif
    else
      if (vsc(l).eq.0.) then
        if (vlc(l).le.soil%pars%w_sat) then
          psi     (l) = (soil%pars%psi_sat_ref/soil%pars%alpha) &
             *(soil%pars%w_sat/vlc(l))**soil%pars%chb
          DKDP    (l) = -(2+3/soil%pars%chb)*hyd_cond(l) &
                                                 /psi(l)
          DThDP   (l) = -vlc(l)/(psi(l)*soil%pars%chb)
        else
          psi(l) = soil%pars%psi_sat_ref &
             + (vlc(l)-soil%pars%w_sat)/comp
          DThDP(l) = comp
          hyd_cond(l) = soil%pars%k_sat_ref
          DKDP(l) = 0.
        endif
      else
        psi     (l) = ((soil%pars%psi_sat_ref/soil%pars%alpha) / 2.2) &
         *(soil%pars%w_sat/vlc(l))**soil%pars%chb
        DKDP    (l) = 0.
        DThDP   (l) = 0.
      endif
    endif
  enddo

  if (DThDP(1).ne.0.) then
    DPsi_min =            -vlc(1) /DThDP(1)
    DPsi_max = (soil%pars%w_sat-vlc(1))/DThDP(1)
  else
    Dpsi_min = -1.e16
    DPsi_max = -psi(1)
  endif

  soil_w_fc = soil%w_fc
  
  ! ---- groundwater ---------------------------------------------------------
  tau_gw = soil%pars%tau_groundwater
 
end subroutine soil_data_hydraulics


! ============================================================================
  subroutine soil_data_w_sat (soil, soil_w_sat  )
  type(soil_tile_type),  intent(in)  :: soil
  real,                  intent(out) :: soil_w_sat(:)

  soil_w_sat(1:num_l) = soil%pars%w_sat

end subroutine soil_data_w_sat

end module soil_tile_mod
