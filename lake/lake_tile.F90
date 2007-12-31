#include <fms_platform.h>

module lake_tile_mod

use fms_mod, only : &
     write_version_number, file_exist, open_namelist_file, check_nml_error, &
     close_file, stdlog
use constants_mod, only : &
     pi, tfreeze
use land_constants_mod, only : &
     NBANDS
use land_io_mod, only : &
     init_cover_field
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_LAKE, register_tile_selector

implicit none
private

! ==== public interfaces =====================================================
public :: lake_pars_type
public :: lake_prog_type
public :: lake_tile_type

public :: new_lake_tile, delete_lake_tile
public :: lake_can_be_merged
public :: lake_is_selected
public :: get_lake_tag

public :: read_lake_data_namelist
public :: lake_cover_cold_start

public :: lake_data_radiation
public :: lake_data_diffusion
public :: lake_data_thermodynamics
public :: lake_data_hydraulics

public :: max_lev
! =====end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter   :: &
     version = '', &
     tagname = '', &
     module_name = 'lake_tile_mod'

integer, parameter :: max_lev          = 30
integer, parameter :: n_dim_lake_types = 1  ! size of lookup table
real,    parameter :: psi_wilt         = -150.  ! matric head at wilting
real,    parameter :: comp             = 0.001  ! m^-1

! ==== types =================================================================
type :: lake_pars_type
  real w_sat
  real awc_lm2
  real k_sat_ref
  real psi_sat_ref
  real chb
  real alpha              ! *** REPLACE LATER BY alpha(layer)
  real heat_capacity_ref
  real thermal_cond_ref
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
end type lake_pars_type

type :: lake_prog_type
  real wl
  real ws
  real T
  real groundwater
  real groundwater_T
end type lake_prog_type

type :: lake_tile_type
   integer :: tag ! kind of the lake
   type(lake_prog_type), _ALLOCATABLE :: prog(:)
   type(lake_pars_type)               :: pars
   real,                 _ALLOCATABLE :: w_fc(:)
   real,                 _ALLOCATABLE :: w_wilt(:)
   real :: Eg_part_ref
   real :: z0_scalar
   real, _ALLOCATABLE :: e(:),f(:)
   real, _ALLOCATABLE :: heat_capacity_dry(:)
end type lake_tile_type

! ==== module data ===========================================================

!---- namelist ---------------------------------------------------------------
real    :: k_over_B              = 0.25      ! reset to 0 for MCM
real    :: rate_fc               = 0.1/86400 ! 0.1 mm/d drainage rate at FC
real    :: sfc_heat_factor       = 1
integer :: z_sfc_layer           = 0
integer :: num_l                 = 18           ! number of lake levels
real    :: dz(max_lev)           = (/ &
    0.02, 0.04, 0.04, 0.05, 0.05, 0.1, 0.1, 0.2, 0.2, &
    0.2,   0.4,  0.4,  0.4,  0.4, 0.4,  1.,  1.,  1., &
    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)
                                              ! thickness (m) of model layers,
                                              ! from top down
logical :: use_lm2_awc           = .false.
  integer :: n_map_1st_lake_type = 10
! ---- remainder are used only for cold start ---------
character(len=16):: lake_to_use     = 'single-tile'
       ! 'multi-tile' for tiled soil [default]
       ! 'single-tile' for geographically varying soil with single type per
       !     model grid cell
       ! 'uniform' for global constant soil, e.g., to reproduce MCM
logical :: use_single_lake       = .false.   ! true for single global lake,
                                             ! e.g., to recover MCM
logical :: use_mcm_albedo        = .false.   ! .true. for CLIMAP albedo inputs
logical :: use_single_geo        = .false.   ! .true. for global gw res time,
                                             ! e.g., to recover MCM
integer :: lake_index_constant   = 1         ! index of global constant lake,
                                             ! used when use_single_lake
real    :: gw_res_time           = 60.*86400 ! mean groundwater residence time,
                                             ! used when use_single_geo
real    :: rsa_exp_global        = 1.5
real, dimension(n_dim_lake_types) :: &
  dat_w_sat             =(/ 1.000   /),&
  dat_awc_lm2           =(/ 1.000   /),&
  dat_k_sat_ref         =(/ 0.021   /),&
  dat_psi_sat_ref       =(/ -.059   /),&
  dat_chb               =(/   3.5   /),&
  dat_heat_capacity_ref =(/ 8.4e7   /),&
  dat_thermal_cond_ref  =(/ 8.4e7   /),&
  dat_emis_dry          =(/ 0.950   /),&
  dat_emis_sat          =(/ 0.980   /),&
  dat_z0_momentum       =(/ 1.4e-4  /),&
  dat_tf_depr           =(/  0.00   /)
real, dimension(n_dim_lake_types, NBANDS) :: &
  dat_refl_dry_dif      =(/ 0.060   ,   & ! visible  
                            0.060   /), & ! NIR
  dat_refl_dry_dir      =(/ 0.060   ,   & ! visible
                            0.060   /), & ! NIR
  dat_refl_sat_dir      =(/ 0.060   ,   & ! visible
                            0.060   /), & ! NIR
  dat_refl_sat_dif      =(/ 0.060   ,   & ! visible
                            0.060   /)    ! NIR
integer, dimension(n_dim_lake_types) :: &
  input_cover_types     =(/ 10 /)
character(len=4), dimension(n_dim_lake_types) :: &
  tile_names            =(/ 'lake' /)

namelist /lake_data_nml/ &
     lake_to_use,input_cover_types, tile_names, &
     k_over_B,             &
     rate_fc, sfc_heat_factor, z_sfc_layer,          &
     num_l,                   dz,                      &
     use_lm2_awc,    n_map_1st_lake_type, &
     use_single_lake,           use_mcm_albedo,            &
     use_single_geo,            lake_index_constant,         &
     gw_res_time,            rsa_exp_global,      &
     dat_w_sat,               dat_awc_lm2,     &
     dat_k_sat_ref,            &
     dat_psi_sat_ref,               dat_chb,          &
     dat_heat_capacity_ref,         dat_thermal_cond_ref,   &
     dat_refl_dry_dir,            dat_refl_sat_dir,              &
     dat_refl_dry_dif,            dat_refl_sat_dif,              &
     dat_emis_dry,              dat_emis_sat,                &
     dat_z0_momentum,           dat_tf_depr

!---- end of namelist --------------------------------------------------------

integer :: num_sfc_layers

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_lake_data_namelist(lake_n_lev, lake_dz)
  integer, intent(out) :: lake_n_lev
  real,    intent(out) :: lake_dz(:)
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
        read (unit, nml=lake_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'lake_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  write (stdlog(), nml=lake_data_nml)

  ! initialize global module data here

  ! register selectors for tile-specific diagnostics
  do i=1, n_dim_lake_types
     call register_tile_selector(tile_names(i), long_name='',&
          tag = SEL_LAKE, idata1 = i )
  enddo

  ! initialize num_sfc_layers
  z = 0
  num_sfc_layers = 0
  do i = 1, num_l
     z = z + dz(i)
     if (z < z_sfc_layer+1.e-4) num_sfc_layers = i
  enddo

  ! set up output arguments
  lake_n_lev = num_l
  lake_dz    = dz
end subroutine 


! ============================================================================
function new_lake_tile(tag) result(ptr)
  type(lake_tile_type), pointer :: ptr ! return value
  integer, intent(in)           :: tag ! kind of lake

  allocate(ptr)
  ptr%tag = tag
  ! allocate storage for tile data
  ! allocate storage for tile data
  allocate(ptr%prog   (num_l),  &
           ptr%w_fc   (num_l),  &
           ptr%w_wilt (num_l),  &
           ptr%heat_capacity_dry(num_l),  &
           ptr%e      (num_l),  &
           ptr%f      (num_l)   )

  call init_lake_data_0d(ptr)

end function new_lake_tile


! ============================================================================
subroutine delete_lake_tile(ptr)
  type(lake_tile_type), pointer :: ptr

  deallocate(ptr%prog, ptr%w_fc, ptr%w_wilt,ptr%heat_capacity_dry, ptr%e, ptr%f)

  deallocate(ptr)
end subroutine delete_lake_tile


subroutine init_lake_data_0d(lake)
  type(lake_tile_type), intent(inout) :: lake

!  real tau_groundwater
!  real rsa_exp         ! riparian source-area exponent

  integer :: k
  k = lake%tag

  lake%pars%w_sat             = dat_w_sat            (k)
  lake%pars%awc_lm2           = dat_awc_lm2          (k)
  lake%pars%k_sat_ref         = dat_k_sat_ref        (k)
  lake%pars%psi_sat_ref       = dat_psi_sat_ref      (k)
  lake%pars%chb               = dat_chb              (k)
  lake%pars%alpha             = 1
  lake%pars%heat_capacity_ref = dat_heat_capacity_ref(k)
  lake%pars%thermal_cond_ref  = dat_thermal_cond_ref (k)
  lake%pars%refl_dry_dir      = dat_refl_dry_dir     (k,:)
  lake%pars%refl_dry_dif      = dat_refl_dry_dif     (k,:)
  lake%pars%refl_sat_dir      = dat_refl_sat_dir     (k,:)
  lake%pars%refl_sat_dif      = dat_refl_sat_dif     (k,:)
  lake%pars%emis_dry          = dat_emis_dry         (k)
  lake%pars%emis_sat          = dat_emis_sat         (k)
  lake%pars%z0_momentum       = dat_z0_momentum      (k)
  lake%pars%tfreeze           = tfreeze - dat_tf_depr(k)

!   lake_type_tfreeze(:,:,iii)       = tfreeze - dat_tf_depr  (k)

  lake%pars%tau_groundwater   = 86400.*30.
  lake%pars%rsa_exp           = rsa_exp_global

  ! -------- derived constant lake parameters --------
  ! w_fc (field capacity) set to w at which hydraulic conductivity equals
  ! a nominal drainage rate "rate_fc"
  ! w_wilt set to w at which psi is psi_wilt
  if (use_lm2_awc) then
     lake%w_wilt(:) = 0.15
     lake%w_fc  (:) = 0.15 + lake%pars%awc_lm2
  else
     lake%w_wilt(:) = lake%pars%w_sat &
          *(lake%pars%psi_sat_ref/(psi_wilt*lake%pars%alpha))**(1/lake%pars%chb)
     lake%w_fc  (:) = lake%pars%w_sat &
          *(rate_fc/(lake%pars%k_sat_ref*lake%pars%alpha**2))**(1/(3+2*lake%pars%chb))
  endif

  ! below made use of phi_e from parlange via entekhabi
  lake%Eg_part_ref = (-4*lake%w_fc(1)**2*lake%pars%k_sat_ref*lake%pars%psi_sat_ref*lake%pars%chb &
       /(pi*lake%pars%w_sat)) * (lake%w_fc(1)/lake%pars%w_sat)**(2+lake%pars%chb)   &
       *(2*pi/(3*lake%pars%chb**2*(1+3/lake%pars%chb)*(1+4/lake%pars%chb)))/2

  lake%z0_scalar = lake%pars%z0_momentum * exp(-k_over_B)

end subroutine 


! ============================================================================
function lake_cover_cold_start(land_mask, glonb, glatb) result (lake_frac)
! creates and initializes a field of fractional lake coverage
! should be called for global grid; othervise the parth that fills
! missing data points may fail to find any good data, or do it in nproc-
! dependent way
  logical, intent(in) :: land_mask(:,:)    ! global land mask
  real,    intent(in) :: glonb(:,:), glatb(:,:)! boundaries of the global grid cells
  real,    pointer    :: lake_frac (:,:,:) ! output-global map of lake fractional coverage

  allocate( lake_frac(size(land_mask,1),size(land_mask,2),n_dim_lake_types))

  call init_cover_field(lake_to_use, 'INPUT/ground_type.nc', 'cover','frac', &
       glonb, glatb, lake_index_constant, input_cover_types, lake_frac)
  
end function 

! =============================================================================
function lake_can_be_merged(lake1,lake2)
  logical :: lake_can_be_merged
  type(lake_tile_type), intent(in) :: lake1,lake2

  lake_can_be_merged = (lake1%tag==lake2%tag)
end function

! =============================================================================
! returns true if tile fits the specified selector
function lake_is_selected(lake, sel)
  logical lake_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(lake_tile_type),      intent(in) :: lake

  lake_is_selected = (sel%idata1 == lake%tag)
end function


! ============================================================================
! retruns tag of the tile
function get_lake_tag(lake) result(tag)
  integer :: tag
  type(lake_tile_type), intent(in) :: lake
  
  tag = lake%tag
end function


! ============================================================================
! compute bare-lake albedos and bare-lake emissivity
subroutine lake_data_radiation ( lake, lake_alb_dir, lake_alb_dif, lake_emis )
  type(lake_tile_type), intent(in)  :: lake
  real,                 intent(out) :: lake_alb_dir(:), lake_alb_dif(:), lake_emis

  ! ---- local vars
  real :: lake_sfc_vlc, blend

  ! ---- radiation properties
  lake_sfc_vlc = lake%prog(1)%wl/dz(1)
  blend        = lake_sfc_vlc/lake%pars%w_sat
  lake_alb_dir = lake%pars%refl_dry_dir + blend*(lake%pars%refl_sat_dir-lake%pars%refl_dry_dir)
  lake_alb_dif = lake%pars%refl_dry_dif + blend*(lake%pars%refl_sat_dif-lake%pars%refl_dry_dif)
  lake_emis = lake%pars%emis_dry   + blend*(lake%pars%emis_sat-lake%pars%emis_dry  )
end subroutine

! ============================================================================
! compute bare-lake roughness
subroutine lake_data_diffusion ( lake,lake_z0s, lake_z0m )
  type(lake_tile_type), intent(in)  :: lake
  real,                 intent(out) :: lake_z0s, lake_z0m

  ! ---- surface roughness
  lake_z0s = lake%z0_scalar
  lake_z0m = lake%pars%z0_momentum
end subroutine

! ============================================================================
! compute lake thermodynamic properties.
subroutine lake_data_thermodynamics ( lake_pars, vlc_sfc, vsc_sfc, &
     lake_rh, heat_capacity_dry, thermal_cond)
  type(lake_pars_type), intent(in)  :: lake_pars
  real,                 intent(in)  :: vlc_sfc
  real,                 intent(in)  :: vsc_sfc
  real,                 intent(out) :: lake_rh
  real,                 intent(out) :: heat_capacity_dry(:)
  real,                 intent(out) :: thermal_cond(:)

  ! ---- local vars
  integer l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  ! assign some index of water availability for snow-free lake
  lake_rh = 1

  ! these will eventually be functions of water content (or psi) and T.
  do l = 1, num_sfc_layers
     heat_capacity_dry(l) = sfc_heat_factor*lake_pars%heat_capacity_ref
     thermal_cond(l)  = sfc_heat_factor*lake_pars%thermal_cond_ref
  enddo
  do l = num_sfc_layers+1, num_l
     heat_capacity_dry(l) = lake_pars%heat_capacity_ref
     thermal_cond(l)  = lake_pars%thermal_cond_ref
  enddo

end subroutine

! ============================================================================
! compute lake hydraulic properties.
subroutine lake_data_hydraulics (lake, vlc, vsc, &
                    psi, DThDP, hyd_cond, DKDP, DPsi_min, DPsi_max, tau_gw, &
                    lake_w_fc  )
  type(lake_tile_type),        intent(in) :: lake
  real,                        intent(in),  dimension(:) :: vlc, vsc
  real,                        intent(out), dimension(:) :: &
      psi, DThDP, hyd_cond, DKDP, lake_w_fc
  real,                        intent(out) :: &
                     DPsi_min, DPsi_max, tau_gw
  ! ---- local vars ----------------------------------------------------------
  integer l
  real,  dimension(num_l) :: vlc_loc
  
  ! ---- T-dependence of hydraulic properties --------------------------------
  ! k_sat   = lake%pars%k_sat0   !  * mu(t0)/mu(t), where mu is dynamic viscosity
  ! psi_sat = lake%pars%psi_sat0 !  * exp(c*(psi-psi0)), where c~+/-(?)0.0068
                     ! better approach would be to adopt air entrapment model
                     ! or at least to scale against surface tension model


  ! ---- water and ice dependence of hydraulic properties --------------------
  ! ---- (T-dependence can be added later)
  hyd_cond=1;DThDP=1
  do l = 1, num_l
    hyd_cond(l) = (lake%pars%k_sat_ref*lake%pars%alpha**2)*  &
                  ! * mu(T)/mu(t_ref), where mu is dynamic viscosity
                 (vlc(l)/lake%pars%w_sat)**(3+2*lake%pars%chb)
    if (hyd_cond(l).lt.1.e-12*lake%pars%k_sat_ref) then
      vlc_loc (l) = lake%pars%w_sat*(1.e-12)**(1./(3+2*lake%pars%chb))
      hyd_cond(l) = 1.e-12*lake%pars%k_sat_ref
      if (vsc(l).eq.0.) then
        DThDP   (l) = -vlc_loc(l)  &
                         *(vlc_loc(l)/lake%pars%w_sat)**lake%pars%chb &
                 /(lake%pars%psi_sat_ref*lake%pars%chb)
        psi     (l) = (lake%pars%psi_sat_ref/lake%pars%alpha) &
            *(lake%pars%w_sat/vlc_loc(l))**lake%pars%chb &
            + (vlc(l)-vlc_loc (l))/DThDP   (l)
        DKDP    (l) = 0.
      else
        psi     (l) = ((lake%pars%psi_sat_ref/lake%pars%alpha) / 2.2) &
            *(lake%pars%w_sat/vlc_loc(l))**lake%pars%chb
        DKDP    (l) = 0.
        DThDP   (l) = 0.
      endif
    else
      if (vsc(l).eq.0.) then
         if (vlc(l).le.lake%pars%w_sat) then
           psi     (l) = (lake%pars%psi_sat_ref/lake%pars%alpha) &
              *(lake%pars%w_sat/vlc(l))**lake%pars%chb
           DKDP    (l) = -(2+3/lake%pars%chb)*hyd_cond(l) &
                                                  /psi(l)
           DThDP   (l) = -vlc(l)/(psi(l)*lake%pars%chb)
         else
           psi(l) = lake%pars%psi_sat_ref &
              + (vlc(l)-lake%pars%w_sat)/comp
           DThDP(l) = comp
           hyd_cond(l) = lake%pars%k_sat_ref
           DKDP(l) = 0.
         endif
       else
         psi     (l) = ((lake%pars%psi_sat_ref/lake%pars%alpha) / 2.2) &
                 *(lake%pars%w_sat/vlc(l))**lake%pars%chb
         DKDP    (l) = 0.
         DThDP   (l) = 0.
       endif
    endif
  enddo

  if (DThDP(1).ne.0.) then
    DPsi_min =            -vlc(1) /DThDP(1)
    DPsi_max = (lake%pars%w_sat-vlc(1))/DThDP(1)
  else
    Dpsi_min = -1.e16
    DPsi_max = -psi(1)
  endif

  lake_w_fc = lake%w_fc
  
  ! ---- groundwater ---------------------------------------------------------
  tau_gw = lake%pars%tau_groundwater
 
end subroutine lake_data_hydraulics


end module lake_tile_mod
