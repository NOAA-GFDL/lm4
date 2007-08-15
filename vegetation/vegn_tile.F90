module vegn_tile_mod

use fms_mod, only : &
     write_version_number, stdlog

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
    init_pars_from_restart, &
    agf_bs

use vegn_cohort_mod, only : vegn_cohort_type, vegn_phys_prog_type, &
     height_from_biomass, lai_from_biomass, update_bio_living_fraction, &
     cohort_data_uptake_frac_max
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

public :: init_derived_vegn_data  ! given state variables, calculate derived values
! =====end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter   :: &
     version = '', &
     tagname = '', &
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

   ! values for the diagnostic of carbon budget and soil carbon acceleration
   real :: asoil_in=0
   real :: ssc_in=0, ssc_out=0
   real :: fsc_in=0, fsc_out=0
   real :: veg_in=0, veg_out=0

   ! it's probably possible to get rid of the fields below
   real :: rh=0 ! soil carbon lost to the atmosphere
   real :: total_biomass ! 
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
subroutine init_derived_vegn_data(vegn)
  type(vegn_tile_type), intent(inout) :: vegn
  
  ! ---- local vars 
  integer :: k ! shorthand for the vegetation tile tag
  type(vegn_cohort_enum_type) :: cc,ce ! cohort enumerators
  type(vegn_cohort_type), pointer :: cohort ! pointer to the cohort at hand
  real :: cmc_max
  
  k = vegn%tag

  ! given that the cohort state variables are initialized, fill in
  ! the intermediate variables
  cc=first_cohort(vegn%cohorts)
  ce=tail_cohort (vegn%cohorts)
  do while(cc/=ce)
    cohort=>current_cohort(cc) 
    cc=next_cohort(cc) ! advance the enumerator
    
    cohort%pt     = spdata(cohort%species)%pt
    cohort%b      = cohort%bliving + cohort%bwood
    call update_bio_living_fraction(cohort)
    cohort%bs     = cohort%bsw + cohort%bwood;   
    cohort%bstem  = agf_bs*cohort%bs;
    cohort%babove = cohort%bl + agf_bs*cohort%bs; 

    if(init_pars_from_restart) then
       ! the height has been read from restart
       ! calculate the leaf area index based on the biomass of leaves
       cohort%lai = lai_from_biomass(cohort%bl, cohort%species)
       ! calculate the root density as the total biomass below ground, in
       ! biomass (not carbon!) units
       cohort%root_density = (cohort%br + &
            (cohort%bsw+cohort%bwood+cohort%blv)*(1-agf_bs))*C2B
    else
       cohort%height        = spdata(k)%dat_height
       cohort%lai           = spdata(k)%dat_lai
       cohort%root_density  = spdata(k)%dat_root_density
    endif
    cohort%sai = 0.035*cohort%height
    cmc_max  = cohort%lai*spdata(cohort%species)%cmc_lai;
    if (cohort%prog%Wl>0.and.cmc_max > 0) then
       cohort%frac_moist = (cohort%prog%Wl/cmc_max)**spdata(cohort%species)%cmc_pow;
    else 
       cohort%frac_moist = 0
    endif
    cohort%leaf_size     = spdata(k)%leaf_size
    cohort%root_zeta     = spdata(k)%dat_root_zeta
    cohort%rs_min        = spdata(k)%dat_rs_min
    cohort%leaf_refl     = spdata(k)%leaf_refl
    cohort%leaf_tran     = spdata(k)%leaf_tran
    cohort%leaf_emis     = spdata(k)%leaf_emis
    cohort%snow_crit     = spdata(k)%dat_snow_crit
  
    ! putting this initialization within the cohort loop is probably incorrect 
    ! in case of multiple-cohort vegetation, however for a single cohort it works
    cohort%W_max   = spdata(k)%cmc_lai*cohort%lai
    cohort%mcv_dry = max(mcv_min, mcv_lai*cohort%lai)
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
