module vegn_tile_mod

use fms_mod, only : &
     write_version_number, stdlog, error_mesg, FATAL

use land_constants_mod, only : NBANDS
use land_io_mod, only : &
     init_cover_field
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_VEGN

use vegn_data_mod, only : NSPECIES, MSPECIES, NCMPT, C2B, n_dim_vegn_types, &
    read_vegn_data_namelist, spdata, &
    vegn_to_use,  input_cover_types, &
    mcv_min, mcv_lai, &
    use_bucket, use_mcm_masking, vegn_index_constant, &
    critical_root_density, &
    agf_bs, BSEED

use vegn_cohort_mod, only : vegn_cohort_type, vegn_phys_prog_type, &
     height_from_biomass, lai_from_biomass, update_bio_living_fraction, &
     cohort_data_uptake_frac_max, update_biomass_pools
use cohort_list_mod, only : &
     vegn_cohort_list_type, vegn_cohort_enum_type, &
     cohort_list_init, cohort_list_end, first_cohort, tail_cohort, next_cohort, &
     operator(/=), current_cohort

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_tile_type

public :: new_vegn_tile, delete_vegn_tile
public :: vegn_can_be_merged
public :: vegn_is_selected
public :: get_vegn_tag

public :: read_vegn_data_namelist
public :: vegn_cover_cold_start

public :: vegn_data_uptake_frac_max
public :: vegn_data_rs_min
public :: vegn_seed_supply
public :: vegn_seed_demand

public :: vegn_add_bliving
public :: update_derived_vegn_data  ! given state variables, calculate derived values
! =====end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter   :: &
     version = '$Id: vegn_tile.F90,v 15.0.2.1 2007/09/16 22:14:24 slm Exp $', &
     tagname = '$Name: omsk_2007_10 $', &
     module_name = 'vegn_tile_mod'

! ==== types =================================================================
type :: vegn_tile_type
   integer :: tag ! kind of the tile

   type(vegn_cohort_list_type) :: cohorts

   real :: age=0 ! tile age

   real :: npp=0 ! net primary productivity
   real :: nep=0 ! net ecosystem productivity

   real :: fast_soil_C=0  ! fast soil carbon pool, (kg C/m2)
   real :: slow_soil_C=0  ! slow soil carbon pool, (kg C/m2)

   ! fields for smoothing out the contribution of the spike-type processes (e.g. 
   ! harvesting) to the soil carbon pools over some period of time
   real :: fsc_pool=0, fsc_rate=0 ! for fast soil carbon
   real :: ssc_pool=0, ssc_rate=0 ! for slow soil carbon

   real :: csmoke_pool=0 ! carbon lost through fires, kg C/m2 
   real :: csmoke_rate=0 ! rate of release of the above to atmosphere, kg C/(m2 yr)

   ! values for the diagnostic of carbon budget and soil carbon acceleration
   real :: asoil_in=0
   real :: ssc_in=0, ssc_out=0
   real :: fsc_in=0, fsc_out=0
   real :: veg_in=0, veg_out=0

   real :: disturbance_rate(0:1) ! 1/year
   real :: lambda   ! cumulative drought months per year
   real :: fuel     ! fuel over dry months
   real :: litter   ! litter flux

   ! monthly accumulated/averaged values
   real :: theta_av ! realtive soil_moisture availablity not soil moisture
   real :: tsoil_av ! bulk soil temperature
   real :: tc_av    ! leaf temperature
   real :: precip_av! precipitation

   ! accumulation counters for long-term averages (monthly and annual). Having
   ! these counters in the tile is a bit stupid, since the values are the same for
   ! each tile, but it simplifies the current code, and they are going away when we
   ! switch to exponential averaging in any case.
   integer :: n_accum ! number of accumulated values for monthly averages
   integer :: nmn_acm ! number of accumulated values for annual averages
   ! annual-mean values
   real :: t_ann    ! annual mean T, degK
   real :: t_cold   ! average temperture of the coldest month, degK
   real :: p_ann    ! annual mean precip
   real :: ncm      ! number of cold months
   ! annual accumulated values
   real :: t_ann_acm  ! accumulated annual temperature for t_ann
   real :: t_cold_acm ! temperature of the coldest month in current year
   real :: p_ann_acm  ! accumulated annual precipitation for p_ann
   real :: ncm_acm    ! accumulated number of cold months


   ! it's probably possible to get rid of the fields below
   real :: rh=0 ! soil carbon lost to the atmosphere
   real :: total_biomass !
   real :: area_disturbed_by_treefall
   real :: area_disturbed_by_fire
   real :: total_disturbance_rate
end type vegn_tile_type

! ==== module data ===========================================================


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
function new_vegn_tile(tag) result(ptr)
  type(vegn_tile_type), pointer :: ptr ! return value
  integer, intent(in)  :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = tag

  call cohort_list_init(ptr%cohorts)

end function new_vegn_tile

! ============================================================================
subroutine delete_vegn_tile(vegn)
  type(vegn_tile_type), pointer :: vegn

  call cohort_list_end(vegn%cohorts)
  deallocate(vegn)
end subroutine delete_vegn_tile

! =============================================================================
function vegn_can_be_merged(vegn1,vegn2)
  logical :: vegn_can_be_merged
  type(vegn_tile_type), intent(in) :: vegn1,vegn2

  vegn_can_be_merged = (vegn1%tag==vegn2%tag)
end function

! ============================================================================
! given a vegetation tile with the state variables set up, calculate derived
! parameters to get a consistent state
subroutine update_derived_vegn_data(vegn)
  type(vegn_tile_type), intent(inout) :: vegn
  
  ! ---- local vars 
  type(vegn_cohort_enum_type) :: ci,ce ! cohort enumerators
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  integer :: sp ! shorthand for the vegetation species
  real :: cmc_max
  
  ! given that the cohort state variables are initialized, fill in
  ! the intermediate variables
  ci=first_cohort(vegn%cohorts) ; ce=tail_cohort (vegn%cohorts)
  do while(ci/=ce)
    cc=>current_cohort(ci) ; ci=next_cohort(ci)
    
    sp = cc%species
    ! set the physiology type according to species
    cc%pt     = spdata(sp)%pt
    ! calculate total biomass, calculate height
    cc%b      = cc%bliving + cc%bwood
    cc%height = height_from_biomass(cc%b);
    ! update fractions of the living biomass
    call update_bio_living_fraction(cc)
    cc%bs     = cc%bsw + cc%bwood;   
    cc%bstem  = agf_bs*cc%bs;
    cc%babove = cc%bl + agf_bs*cc%bs; 

    if(sp<NSPECIES) then ! LM3V species
       ! calculate the leaf area index based on the biomass of leaves
       cc%lai = lai_from_biomass(cc%bl, sp)
       ! calculate the root density as the total biomass below ground, in
       ! biomass (not carbon!) units
       cc%root_density = (cc%br + (cc%bsw+cc%bwood+cc%blv)*(1-agf_bs))*C2B
    else
       cc%height        = spdata(sp)%dat_height
       cc%lai           = spdata(sp)%dat_lai
       cc%root_density  = spdata(sp)%dat_root_density
    endif
    cc%sai = 0.035*cc%height
    cmc_max  = cc%lai*spdata(sp)%cmc_lai;
    if (cc%prog%Wl>0.and.cmc_max > 0) then
       cc%frac_moist = (cc%prog%Wl/cmc_max)**spdata(sp)%cmc_pow;
    else 
       cc%frac_moist = 0
    endif
    cc%leaf_size     = spdata(sp)%leaf_size
    cc%root_zeta     = spdata(sp)%dat_root_zeta
    cc%rs_min        = spdata(sp)%dat_rs_min
    cc%leaf_refl     = spdata(sp)%leaf_refl
    cc%leaf_tran     = spdata(sp)%leaf_tran
    cc%leaf_emis     = spdata(sp)%leaf_emis
    cc%snow_crit     = spdata(sp)%dat_snow_crit
  
    ! putting this initialization within the cohort loop is probably incorrect 
    ! in case of multiple-cohort vegetation, however for a single cohort it works
    cc%W_max   = spdata(sp)%cmc_lai*cc%lai
    cc%mcv_dry = max(mcv_min, mcv_lai*cc%lai)
  enddo
    
end subroutine

! ============================================================================
! vegn_data_uptake_frac_max probably should stay here, since it is called 
! by soil, and exposing internals of the vegetation tile to the other modules
! is not good. By internals I mean the cohort list
subroutine vegn_data_uptake_frac_max(vegn,dz,uptake_frac_max)
  type(vegn_tile_type), intent(in)  :: vegn
  real,                 intent(in)  :: dz(:)
  real,                 intent(out) :: uptake_frac_max(:)

  type(vegn_cohort_type), pointer :: cohort
  
  ! get the pointer to the first (and, currently, the only) cohort
  ! of the tile
  cohort => current_cohort(first_cohort(vegn%cohorts))
  call cohort_data_uptake_frac_max(cohort, dz, uptake_frac_max)

end subroutine

! ============================================================================
function vegn_data_rs_min ( vegn )
  real :: vegn_data_rs_min
  type(vegn_tile_type), intent(in)  :: vegn
  type(vegn_cohort_type), pointer :: cohort
  
  ! get the pointer to the first (and, currently, the only) cohort
  ! of the tile
  cohort => current_cohort(first_cohort(vegn%cohorts))
  vegn_data_rs_min = cohort%rs_min
end function


! ============================================================================
function vegn_seed_supply ( vegn )
  real :: vegn_seed_supply
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars 
  type(vegn_cohort_enum_type) :: ci,ce ! cohort enumerators
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  real :: vegn_bliving
  
  vegn_bliving = 0
  ci=first_cohort(vegn%cohorts) ; ce=tail_cohort (vegn%cohorts)
  do while(ci/=ce)
     cc=>current_cohort(ci) ; ci=next_cohort(ci)
     vegn_bliving = vegn_bliving + cc%bliving
  enddo
  vegn_seed_supply = MAX (vegn_bliving-BSEED, 0.0)
  
end function 

! ============================================================================
function vegn_seed_demand ( vegn )
  real :: vegn_seed_demand
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars 
  type(vegn_cohort_enum_type) :: ci,ce ! cohort enumerators
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort

  vegn_seed_demand = 0
  ci=first_cohort(vegn%cohorts) ; ce=tail_cohort (vegn%cohorts)
  do while(ci/=ce)
     cc=>current_cohort(ci) ; ci=next_cohort(ci)
     if(cc%bliving<BSEED.and.vegn%t_ann>253.16.and.vegn%p_ann>1E-6) then
        vegn_seed_demand = vegn_seed_demand + BSEED
     endif
  enddo
end function 

! ============================================================================
subroutine vegn_add_bliving ( vegn, delta )
  type(vegn_tile_type), intent(inout) :: vegn
  real :: delta ! increment of bliving

  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  cc=>current_cohort(first_cohort(vegn%cohorts))
  cc%bliving = cc%bliving + delta

  if (cc%bliving < 0)then
     call error_mesg('vegn_add_bliving','resulting bliving is less then 0', FATAL)
  endif
  call update_biomass_pools(cc)
end subroutine 

! ============================================================================
function vegn_cover_cold_start(land_mask, glonb, glatb) result (vegn_frac)
! creates and initializes a field of fractional vegn coverage
! should be called for global grid; otherwise the part that fills
! missing data points may fail to find any good data, or do it in nproc-
! dependent way
  logical, intent(in) :: land_mask(:,:)    ! global land mask
  real,    intent(in) :: glonb(:), glatb(:)! boundaries of the global grid cells
  real,    pointer    :: vegn_frac (:,:,:) ! output-global map of vegn fractional coverage

!  allocate( vegn_frac(size(land_mask,1),size(land_mask,2),n_dim_vegn_types))
  allocate( vegn_frac(size(land_mask,1),size(land_mask,2),MSPECIES))

  call init_cover_field(vegn_to_use, 'INPUT/cover_type.nc', 'cover','frac', &
       glonb, glatb, vegn_index_constant, input_cover_types, vegn_frac)
  
end function 

! =============================================================================
! returns true if tile fits the specified selector
function vegn_is_selected(vegn, sel)
  logical vegn_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(vegn_tile_type),      intent(in) :: vegn

  vegn_is_selected = (sel%idata1 == vegn%tag)
end function


! ============================================================================
! returns tag of the tile
function get_vegn_tag(vegn) result(tag)
  integer :: tag
  type(vegn_tile_type), intent(in) :: vegn
  
  tag = vegn%tag
end function

end module vegn_tile_mod
