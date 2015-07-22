module vegn_tile_mod

use fms_mod, only : &
     write_version_number, stdlog, error_mesg, FATAL
use constants_mod, only : &
     tfreeze, hlf

use land_constants_mod, only : NBANDS
use land_io_mod, only : &
     init_cover_field
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_VEGN

use vegn_data_mod, only : &
     MSPECIES, spdata, &
     read_vegn_data_namelist, &
     vegn_to_use,  input_cover_types, &
     mcv_min, mcv_lai, &
     vegn_index_constant, &
     BSEED, LU_NTRL, LU_SCND, N_HARV_POOLS, &
     LU_SEL_TAG, SP_SEL_TAG, NG_SEL_TAG, &
     FORM_GRASS, &
     scnd_biomass_bins

use vegn_cohort_mod, only : vegn_cohort_type, update_biomass_pools

use soil_tile_mod, only : max_lev

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_tile_type

public :: new_vegn_tile, delete_vegn_tile
public :: vegn_tiles_can_be_merged, merge_vegn_tiles
public :: vegn_is_selected
public :: get_vegn_tile_tag
public :: vegn_tile_stock_pe
public :: vegn_tile_carbon ! returns total carbon per tile [kgC/m2]
public :: vegn_tile_bwood  ! returns total woody biomass of tile [kgC/m2]
public :: vegn_tile_heat   ! returns heat content of the vegetation [J/m2]
public :: vegn_tile_LAI    ! returns total LAI of vegetation [m2/m2]
public :: vegn_tile_SAI    ! returns total SAI of vegetation [m2/m2]

public :: read_vegn_data_namelist
public :: vegn_cover_cold_start

public :: vegn_seed_supply
public :: vegn_seed_demand

public :: vegn_tran_priority ! returns transition priority for land use 

public :: vegn_add_bliving
! =====end of public interfaces ==============================================

interface new_vegn_tile
   module procedure vegn_tile_ctor
   module procedure vegn_tile_copy_ctor
end interface


! ==== module constants ======================================================
character(len=*), parameter   :: &
     version = '$Id$', & 
     tagname = '$Name$', &
     module_name = 'vegn_tile_mod'

! ==== types =================================================================
type :: vegn_tile_type
   integer :: tag ! kind of the tile
   integer :: landuse = LU_NTRL

   integer :: n_cohorts = 0
   type(vegn_cohort_type), pointer :: cohorts(:)=>NULL()
   
   ! buffers holding amounts of stuff dropped by dead trees, etc
   real :: drop_wl=0, drop_ws=0 ! buffer accumulating amount of dropped water, kg/m2
   real :: drop_hl=0, drop_hs=0 ! buffer accumulating heat of dropped water, J/m2

   real :: age=0.0 ! tile age

   ! fields for smoothing out the contribution of the spike-type processes (e.g. 
   ! harvesting) to the soil carbon pools over some period of time
   real :: fsc_pool_ag=0.0, fsc_rate_ag=0.0 ! for fast soil carbon above ground
   real :: ssc_pool_ag=0.0, ssc_rate_ag=0.0 ! for slow soil carbon above ground
   real :: fsc_pool_bg=0.0, fsc_rate_bg=0.0 ! for fast soil carbon below ground
   real :: ssc_pool_bg=0.0, ssc_rate_bg=0.0 ! for slow soil carbon below ground
   
   real :: leaflitter_buffer_ag=0.0,      coarsewoodlitter_buffer_ag=0.0
   real :: leaflitter_buffer_rate_ag=0.0, coarsewoodlitter_buffer_rate_ag=0.0

   real :: csmoke_pool=0.0 ! carbon lost through fires, kg C/m2 
   real :: csmoke_rate=0.0 ! rate of release of the above to atmosphere, kg C/(m2 yr)

   real :: harv_pool(N_HARV_POOLS) = 0.0 ! pools of harvested carbon, kg C/m2
   real :: harv_rate(N_HARV_POOLS) = 0.0 ! rates of spending (release to the atmosphere), kg C/(m2 yr)

   ! uptake-related variables
   real :: root_distance(max_lev) ! characteristic half-distance between fine roots, m
   
   ! values for the diagnostic of carbon budget and soil carbon acceleration
   real :: ssc_out=0.0
   real :: fsc_out=0.0
   real :: deadmic_out=0.0
   real :: veg_in=0.0, veg_out=0.0

   real :: disturbance_rate(0:1) = 0 ! 1/year
   real :: lambda = 0.0 ! cumulative drought months per year
   real :: fuel   = 0.0 ! fuel over dry months
   real :: litter = 0.0 ! litter flux

   ! monthly accumulated/averaged values
   real :: theta_av_phen = 0.0 ! relative soil_moisture availability not soil moisture
   real :: theta_av_fire = 0.0
   real :: psist_av = 0.0 ! soil water stress index
   real :: tsoil_av = 0.0 ! bulk soil temperature
   real :: tc_av    = 0.0 ! leaf temperature
   real :: precip_av= 0.0 ! precipitation

   ! accumulation counters for long-term averages (monthly and annual). Having
   ! these counters in the tile is a bit stupid, since the values are the same for
   ! each tile, but it simplifies the current code, and they are going away when we
   ! switch to exponential averaging in any case.
   integer :: n_accum = 0 ! number of accumulated values for monthly averages
   integer :: nmn_acm = 0 ! number of accumulated values for annual averages
   ! annual-mean values
   real :: t_ann  = 0.0 ! annual mean T, degK
   real :: t_cold = 0.0 ! average temperature of the coldest month, degK
   real :: p_ann  = 0.0 ! annual mean precip
   real :: ncm    = 0.0 ! number of cold months
   ! annual accumulated values
   real :: t_ann_acm  = 0.0 ! accumulated annual temperature for t_ann
   real :: t_cold_acm = 0.0 ! temperature of the coldest month in current year
   real :: p_ann_acm  = 0.0 ! accumulated annual precipitation for p_ann
   real :: ncm_acm    = 0.0 ! accumulated number of cold months

   ! averaged quantities for PPA phenology
   real :: tc_daily = 0.0
   real :: tc_pheno = 0.0 ! smoothed canopy air temperature for phenology
   real :: daily_t_max = -HUGE(1.0) ! accumulator for daily max, used in GDD calculations
   real :: daily_t_min =  HUGE(1.0) ! accumulator for daily min, used in GDD calculations

   ! it's probably possible to get rid of the fields below
   real :: npp=0.0 ! net primary productivity
   real :: nep=0.0 ! net ecosystem productivity
   real :: rh =0.0 ! soil carbon lost to the atmosphere
end type vegn_tile_type

! ==== module data ===========================================================
real, public :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
function vegn_tile_ctor(tag) result(ptr)
  type(vegn_tile_type), pointer :: ptr ! return value
  integer, intent(in)  :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = tag
end function vegn_tile_ctor

! ============================================================================
function vegn_tile_copy_ctor(vegn) result(ptr)
  type(vegn_tile_type), pointer :: ptr ! return value
  type(vegn_tile_type), intent(in) :: vegn ! return value

  allocate(ptr)
  ! copy all non-pointer members
  ptr=vegn
  ! copy pointer members (cohorts)
  allocate(ptr%cohorts(ptr%n_cohorts))
  ptr%cohorts(:) = vegn%cohorts(1:ptr%n_cohorts)
end function vegn_tile_copy_ctor

! ============================================================================
subroutine delete_vegn_tile(vegn)
  type(vegn_tile_type), pointer :: vegn

  deallocate(vegn%cohorts)
  deallocate(vegn)
end subroutine delete_vegn_tile

! =============================================================================
function vegn_tiles_can_be_merged(vegn1,vegn2) result(response)
  logical :: response
  type(vegn_tile_type), intent(in) :: vegn1,vegn2

  real    :: b1, b2 
  integer :: i, i1, i2

  if (vegn1%landuse /= vegn2%landuse) then
     response = .false. ! different land use types can't be merged
  else if (vegn1%landuse == LU_SCND) then ! secondary vegetation tiles
     ! get tile wood biomasses
     b1 = vegn_tile_bwood(vegn1)
     b2 = vegn_tile_bwood(vegn2)
     ! find biomass bins where each the tiles belongs to
     i1 = 0 ; i2 = 0
     do i = 1, size(scnd_biomass_bins(:))
        if (b1>scnd_biomass_bins(i)) i1 = i
        if (b2>scnd_biomass_bins(i)) i2 = i
     enddo
     ! tiles can be merged only if biomasses belong to the same bin
     response = (i1 == i2)
  else
     response = .true. ! non-secondary tiles of the same land use type can always be merged
  endif
end function vegn_tiles_can_be_merged


! ============================================================================
! TODO: update merge_vegn_tiles for PPA
subroutine merge_vegn_tiles(t1,w1,t2,w2)
  type(vegn_tile_type), intent(in) :: t1
  type(vegn_tile_type), intent(inout) :: t2
  real, intent(in) :: w1, w2 ! relative weights
  
  ! ---- local vars
  real :: x1, x2 ! normalized relative weights
  real :: HEAT1, HEAT2 ! heat stored in respective canopies
  type(vegn_cohort_type), pointer :: c1, c2
  
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  ! TODO: reformulate merging cohorts when merging tiles in PPA mode
  !       probably combine all cohorts from two tiles + relayer + merge cohorts

  ! the following assumes that there is one, and only one, cohort per tile 
  c1 => t1%cohorts(1)
  c2 => t2%cohorts(1)
  ! define macro for merging cohort values
#define __MERGE__(field) c2%field = x1*c1%field + x2*c2%field
  HEAT1 = (clw*c1%Wl + csw*c1%Ws + c1%mcv_dry)*(c1%Tv-tfreeze)
  HEAT2 = (clw*c2%Wl + csw*c2%Ws + c2%mcv_dry)*(c2%Tv-tfreeze)
  __MERGE__(Wl)
  __MERGE__(Ws)

  __MERGE__(bl)      ! biomass of leaves, kg C/m2
  __MERGE__(blv)     ! biomass of virtual leaves (labile store), kg C/m2
  __MERGE__(br)      ! biomass of fine roots, kg C/m2
  __MERGE__(bsw)     ! biomass of sapwood, kg C/m2
  __MERGE__(bwood)   ! biomass of heartwood, kg C/m2
  __MERGE__(bliving) ! leaves, fine roots, and sapwood biomass
  __MERGE__(nsc)     ! non-structural carbon, kgC/indiv
  __MERGE__(bseed)   ! future progeny, kgC/indiv
  __MERGE__(age)     ! age of individual

  __MERGE__(carbon_gain) ! carbon gain during a day, kg C/m2
  __MERGE__(carbon_loss) ! carbon loss during a day, kg C/m2 [diag only]
  __MERGE__(bwood_gain)  ! heartwood gain during a day, kg C/m2

  ! should we do update_derived_vegn_data here? to get mcv_dry, etc
  call update_biomass_pools(c2)

  ! calculate the resulting dry heat capacity
  c2%mcv_dry = max(mcv_min,mcv_lai*c2%lai)
  ! update canopy temperature -- just merge it based on area weights if the heat 
  ! capacities are zero, or merge it based on the heat content if the heat contents
  ! are non-zero
  if(HEAT1==0.and.HEAT2==0) then
     __MERGE__(Tv)
  else
     c2%Tv = (HEAT1*x1+HEAT2*x2) / &
          (clw*c2%Wl + csw*c2%Ws + c2%mcv_dry) + tfreeze
  endif

#undef  __MERGE__
! re-define macro for tile values
#define __MERGE__(field) t2%field = x1*t1%field + x2*t2%field
  __MERGE__(drop_wl)
  __MERGE__(drop_ws)
  __MERGE__(drop_hl)
  __MERGE__(drop_hs)

  __MERGE__(age);
  
  __MERGE__(fsc_pool_ag); __MERGE__(fsc_rate_ag)
  __MERGE__(ssc_pool_ag); __MERGE__(ssc_rate_ag)
  __MERGE__(fsc_pool_bg); __MERGE__(fsc_rate_bg)
  __MERGE__(ssc_pool_bg); __MERGE__(ssc_rate_bg)
  
  __MERGE__(leaflitter_buffer_ag)
  __MERGE__(leaflitter_buffer_rate_ag)
  __MERGE__(coarsewoodlitter_buffer_ag)
  __MERGE__(coarsewoodlitter_buffer_rate_ag)

  __MERGE__(csmoke_pool)
  __MERGE__(csmoke_rate)

  __MERGE__(harv_pool)
  __MERGE__(harv_rate)

  ! do we need to merge these?
  __MERGE__(ssc_out)
  __MERGE__(fsc_out)
  __MERGE__(deadmic_out)
  __MERGE__(veg_in); __MERGE__(veg_out)
  
  ! or these?
  __MERGE__(disturbance_rate)
  __MERGE__(lambda)     ! cumulative drought months per year
  __MERGE__(fuel)       ! fuel over dry months
  __MERGE__(litter)     ! litter flux

  ! monthly accumulated/averaged values
  __MERGE__(theta_av_phen)   ! relative soil_moisture availability not soil moisture
  __MERGE__(theta_av_fire)
  __MERGE__(psist_av)   ! water potential divided by permanent wilting potential
  __MERGE__(tsoil_av)   ! bulk soil temperature
  __MERGE__(tc_av)      ! leaf temperature
  __MERGE__(precip_av)  ! precipitation

  ! annual-mean values
  __MERGE__(t_ann)      ! annual mean T, degK
  __MERGE__(t_cold)     ! average temperature of the coldest month, degK
  __MERGE__(p_ann)      ! annual mean precip
  __MERGE__(ncm)        ! number of cold months

  ! annual accumulated values
  __MERGE__(t_ann_acm)  ! accumulated annual temperature for t_ann
  __MERGE__(t_cold_acm) ! temperature of the coldest month in current year
  __MERGE__(p_ann_acm)  ! accumulated annual precipitation for p_ann
  __MERGE__(ncm_acm)    ! accumulated number of cold months

#undef __MERGE__

end subroutine merge_vegn_tiles


! ============================================================================
function vegn_seed_supply ( vegn )
  real :: vegn_seed_supply
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars 
  real :: vegn_bliving
  integer :: i
  
  vegn_bliving = 0
  do i = 1,vegn%n_cohorts
     vegn_bliving = vegn_bliving + vegn%cohorts(i)%bliving
  enddo
  vegn_seed_supply = MAX (vegn_bliving-BSEED, 0.0)
  
end function vegn_seed_supply

! ============================================================================
! TODO: do vegn_seed_demand and vegn_seed_supply make sense for PPA? or even the entire seed transport method?
function vegn_seed_demand ( vegn )
  real :: vegn_seed_demand
  type(vegn_tile_type), intent(in) :: vegn

  integer :: i

  vegn_seed_demand = 0
  do i = 1,vegn%n_cohorts
     if(vegn%cohorts(i)%bliving<BSEED.and.vegn%t_ann>253.16.and.vegn%p_ann>1E-6) then
        vegn_seed_demand = vegn_seed_demand + BSEED
     endif
  enddo
end function vegn_seed_demand

! ============================================================================
subroutine vegn_add_bliving ( vegn, delta )
  type(vegn_tile_type), intent(inout) :: vegn
  real :: delta ! increment of bliving

  vegn%cohorts(1)%bliving = vegn%cohorts(1)%bliving + delta

  if (vegn%cohorts(1)%bliving < 0)then
     call error_mesg('vegn_add_bliving','resulting bliving is less then 0', FATAL)
  endif
  call update_biomass_pools(vegn%cohorts(1))
end subroutine vegn_add_bliving





! ============================================================================
! given a vegetation patch, destination kind of transition, and "transition 
! intensity" value, this function returns a fraction of tile that will
! participate in transition.
!
! this function must be contiguous, monotonic, its value must be within
! interval [0,1]
!
! this function is used to determine what part of each tile is to be converted
! to another land use kind; the equation is solved to get "transition intensity" 
! tau for which total area is equal to requested. Tau is, therefore, a dummy
! parameter, and only relative values of the priority functions for tiles 
! participating in transition have any meaning. For most transitions the priority 
! function is just equal to tau: therefore there is no preference, and all tiles
! contribute equally to converted area. For secondary vegetation harvesting, 
! however, priority also depends on wood biomass, and therefore tiles
! with high wood biomass are harvested first.
function vegn_tran_priority(vegn, dst_kind, tau) result(pri)
  real :: pri
  type(vegn_tile_type), intent(in) :: vegn
  integer             , intent(in) :: dst_kind
  real                , intent(in) :: tau

  real :: vegn_bwood
  integer :: i

  if (vegn%landuse==LU_SCND.and.dst_kind==LU_SCND) then ! secondary biomass harvesting
     vegn_bwood = vegn_tile_bwood(vegn)
     pri = max(min(tau+vegn_bwood,1.0),0.0)
  else
     pri = max(min(tau,1.0),0.0)
  endif
end function vegn_tran_priority


! ============================================================================
function vegn_cover_cold_start(land_mask, lonb, latb) result (vegn_frac)
! creates and initializes a field of fractional vegn coverage
  logical, intent(in) :: land_mask(:,:)    ! land mask
  real,    intent(in) :: lonb(:,:), latb(:,:)! boundaries of the grid cells
  real,    pointer    :: vegn_frac (:,:,:) ! output: map of vegn fractional coverage

  allocate( vegn_frac(size(land_mask,1),size(land_mask,2),MSPECIES))

  call init_cover_field(vegn_to_use, 'INPUT/cover_type.nc', 'cover','frac', &
       lonb, latb, vegn_index_constant, input_cover_types, vegn_frac)
  
end function vegn_cover_cold_start

! =============================================================================
! returns true if tile fits the specified selector
function vegn_is_selected(vegn, sel)
  logical vegn_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(vegn_tile_type),      intent(in) :: vegn

  integer :: sp ! shorthand for vegetation species
  
  select case (sel%idata1)
  case (LU_SEL_TAG)
     vegn_is_selected = (sel%idata2 == vegn%landuse)
  case (SP_SEL_TAG)
     if (.not.associated(vegn%cohorts)) then
        vegn_is_selected = .FALSE.
     else
        vegn_is_selected = (sel%idata2 == vegn%cohorts(1)%species)
     endif
  case (NG_SEL_TAG)
     if (.not.associated(vegn%cohorts)) then
        vegn_is_selected = .FALSE.
     else
        sp = vegn%cohorts(1)%species
        vegn_is_selected = &
             (spdata(sp)%form==FORM_GRASS).and.&
             ((vegn%landuse==LU_NTRL).or.(vegn%landuse==LU_SCND))
     endif
  case default
     vegn_is_selected = .FALSE.
  end select  
     
end function vegn_is_selected


! ============================================================================
! returns tag of the tile
function get_vegn_tile_tag(vegn) result(tag)
  integer :: tag
  type(vegn_tile_type), intent(in) :: vegn
  
  tag = vegn%tag
end function get_vegn_tile_tag

! ============================================================================
! returns total wood biomass per tile 
function vegn_tile_bwood(vegn) result(bwood)
  real :: bwood
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars
  integer :: i

  bwood = 0
  do i = 1,vegn%n_cohorts
     bwood = bwood + vegn%cohorts(i)%bwood*vegn%cohorts(i)%nindivs
  enddo
end function vegn_tile_bwood

! ============================================================================
! returns total leaf area index
function vegn_tile_LAI(vegn) result(LAI) ; real LAI
  type(vegn_tile_type), intent(in) :: vegn

  integer :: i

  LAI = 0
  do i = 1,vegn%n_cohorts
     LAI = LAI + vegn%cohorts(i)%lai*vegn%cohorts(i)%nindivs
  enddo
end function vegn_tile_LAI

! ============================================================================
! returns total stem area index
function vegn_tile_SAI(vegn) result(SAI) ; real SAI
  type(vegn_tile_type), intent(in) :: vegn

  integer :: i

  SAI = 0
  do i = 1,vegn%n_cohorts
     SAI = SAI + vegn%cohorts(i)%sai*vegn%cohorts(i)%nindivs
  enddo
end function vegn_tile_SAI

! ============================================================================
subroutine vegn_tile_stock_pe (vegn, twd_liq, twd_sol  )
  type(vegn_tile_type),  intent(in)    :: vegn
  real,                  intent(out)   :: twd_liq, twd_sol
  integer i
  
  twd_liq = vegn%drop_wl ; twd_sol = vegn%drop_ws
  do i=1, vegn%n_cohorts
    twd_liq = twd_liq + vegn%cohorts(i)%wl * vegn%cohorts(i)%nindivs
    twd_sol = twd_sol + vegn%cohorts(i)%ws * vegn%cohorts(i)%nindivs
  enddo
end subroutine vegn_tile_stock_pe


! ============================================================================
! returns total carbon in the tile, kg C/m2
function vegn_tile_carbon(vegn) result(carbon) ; real carbon
  type(vegn_tile_type), intent(in)  :: vegn

  integer :: i

  carbon = 0
  do i = 1,vegn%n_cohorts
     carbon = carbon + &
         (vegn%cohorts(i)%bl  + vegn%cohorts(i)%blv + &
          vegn%cohorts(i)%br  + vegn%cohorts(i)%bwood + &
          vegn%cohorts(i)%bsw + vegn%cohorts(i)%bseed + &
          vegn%cohorts(i)%nsc + &
          vegn%cohorts(i)%carbon_gain + vegn%cohorts(i)%bwood_gain &
         )*vegn%cohorts(i)%nindivs
  enddo
  carbon = carbon + sum(vegn%harv_pool) + &
           vegn%fsc_pool_ag + vegn%ssc_pool_ag + &
           vegn%fsc_pool_bg + vegn%ssc_pool_bg + vegn%csmoke_pool
end function vegn_tile_carbon


! ============================================================================
! returns heat content of the vegetation, J/m2
function vegn_tile_heat (vegn) result(heat) ; real heat
  type(vegn_tile_type), intent(in)  :: vegn

  integer :: i

  heat = vegn%drop_hl+vegn%drop_hs
  do i = 1, vegn%n_cohorts
     heat = heat + &
            ( (clw*vegn%cohorts(i)%Wl + &
               csw*vegn%cohorts(i)%Ws + &
               vegn%cohorts(i)%mcv_dry)*(vegn%cohorts(i)%Tv-tfreeze) &
              - hlf*vegn%cohorts(i)%Ws &
            )*vegn%cohorts(i)%nindivs
  enddo
end function vegn_tile_heat

end module vegn_tile_mod
