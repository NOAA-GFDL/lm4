!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Land Model 4 (LM4).
!*
!* LM4 is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* LM4 is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with LM4.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
module vegn_cohort_mod

#include "../shared/debug.inc"

use constants_mod, only: PI

use land_constants_mod, only: NBANDS, mol_h2o, mol_air
use vegn_data_mod, only : spdata, &
   use_bucket, critical_root_density, &
   tg_c4_thresh, tg_c3_thresh, T_cold_tropical, &
   phen_ev1, phen_ev2, cmc_eps, sai_cover, N_limits_live_biomass, &
   SP_C4GRASS, SP_C3GRASS, SP_TEMPDEC, SP_TROPICAL, SP_EVERGR, &
   LEAF_OFF, LU_CROP, PHEN_EVERGREEN, PHEN_DECIDUOUS, FORM_GRASS, &
   ALLOM_EW, ALLOM_EW1, ALLOM_HML, PT_C3, PT_C4, &
   do_ppa, DBH_merge_rel, DBH_merge_abs, NSC_merge_rel, &
   snow_masking_option, permafrost_depth_thresh, permafrost_freq_thresh, &
   SNOW_MASKING_NONE, SNOW_MASKING_LM3, SNOW_MASKING_MCM, SNOW_MASKING_HEIGHT
use soil_tile_mod, only : soil_tile_type, max_lev, num_l, dz
use soil_carbon_mod, only : soil_carbon_option,SOILC_CORPSE_N
use fms_mod, only : error_mesg, FATAL
use land_debug_mod, only : is_watch_point

implicit none
private
! ==== public interfaces =====================================================
public :: vegn_cohort_type

! operations defined for cohorts
public :: plant_C, plant_N
public :: get_vegn_wet_frac
public :: vegn_data_cover
public :: cohort_uptake_profile
public :: cohort_root_litter_profile
public :: cohort_root_exudate_profile

public :: btotal ! returns cohort total biomass
public :: update_species ! updates cohort physiology, phenology type, and species
public :: leaf_area_from_biomass ! given leaf biomass, calculates leaf area
public :: height_from_biomass    ! given total biomass, calculated tree height
public :: update_bio_living_fraction
public :: update_biomass_pools
public :: update_cohort_root_properties
public :: init_cohort_allometry_ppa
public :: init_cohort_hydraulics
public :: cohorts_can_be_merged
public :: cohort_makes_seeds ! returns true if cohort can allocate to seeds
public :: cohort_can_reproduce ! returns true if cohort can reproduce
! ==== end of public interfaces ==============================================

! ==== types =================================================================
! vegn_cohort_type describes the data that belong to a vegetation cohort
type :: vegn_cohort_type
  real :: Wl ! liquid intercepted water, kg/individual
  real :: Ws ! solid intercepted water (i.e. snow), kg/individual
  real :: Tv ! canopy temperature, K

! ---- biological prognostic variables
! Currently bio prognostic variable is defined as anything that is saved in
! the restart; clearly some vars there are, strictly speaking, diagnostic,
! but saved for reproducibility (to avoid recalculation). exceptions :
! npp_previous_day is left outside, since it is obviously auxiliary; height
! is left outside
  integer :: species = 0   ! vegetation species
  real    :: bl      = 0.0 ! biomass of leaves, kg C/individual
  real    :: blv     = 0.0 ! biomass of virtual leaves (labile store), kg C/individual
  real    :: br      = 0.0 ! biomass of fine roots, kg C/individual
  real    :: bsw     = 0.0 ! biomass of sapwood, kg C/individual
  real    :: bwood   = 0.0 ! biomass of heartwood, kg C/individual
  real    :: bseed   = 0.0 ! biomass put aside for future progeny, kg C/individual
  real    :: nsc     = 0.0 ! non-structural carbon, kg C/individual
  real    :: bl_wilt = 0.0 ! biomass of leaves in wilted pool, kg C/individual

  real    :: bliving = 0.0 ! leaves, fine roots, and sapwood biomass
  integer :: status  = 0   ! growth status of plant
  real    :: leaf_age= 0.0 ! age of leaf in days since budburst

  ! Nitrogen vegetation pools
  real :: stored_N = 0.0
  real :: leaf_N = 0.0
  real :: wood_N = 0.0
  real :: sapwood_N = 0.0
  real :: root_N = 0.0
  real :: seed_N = 0.0
  real :: total_N = 0.0 ! sum of plant nitrogen pools (does not includes fixers and mycorrahize)

! ---- physical parameters
  real    :: height       = 0.0 ! vegetation height, m
  real    :: zbot         = 0.0 ! height of bottom of the canopy, m
  real    :: lai          = 0.0 ! leaf area index, m2/m2
  real    :: sai          = 0.0 ! stem area index, m2/m2
  real    :: leaf_size    = 0.0 ! leaf dimension, m
  real    :: root_density = 0.0 ! total biomass below ground, kg B/m2 -- used in bucket formulation of water uptake
  real    :: root_zeta    = 0.0 ! e-folding depth of root biomass, m -- used in many places
  real    :: rs_min       = 0.0 ! min. stomatal conductance -- used in VEGN_PHOT_SIMPLE only
  real    :: leaf_refl(NBANDS) = 0.0 ! reflectance of leaf, per band
  real    :: leaf_tran(NBANDS) = 0.0 ! transmittance of leaf, per band
  real    :: leaf_emis    = 0.0 ! emissivity of leaf
  real    :: snow_crit    = 0.0 ! later parameterize this as snow_mask_fac*height

! ---- PPA-related variables
  real    :: age          = 0.0 ! age of cohort, years
  real    :: dbh          = 0.0 ! diameter at breast height, m
  real    :: crownarea    = 1.0 ! crown area, m2/individual
  real    :: leafarea     = 0.0 ! total area of leaves, m2/individual
  real    :: nindivs      = 1.0 ! density of vegetation, individuals/m2
  integer :: layer        = 1   ! the layer of this cohort (numbered from top)
  integer :: firstlayer   = 0   ! 0 = never been in the first layer; 1 = at least one year in first layer
  real    :: layerfrac    = 0.0 ! fraction of layer area occupied by this cohort, m2 of cohort per m2 of ground
  real    :: bl_max       = 0.0 ! Max. leaf biomass, kg C/individual
  real    :: br_max       = 0.0 ! Max. fine root biomass, kg C/individual
  real    :: topyear      = 0.0 ! the years that a plant in top layer
  ! isa 20170707 - required to decompose stem of grasses
  real    :: bsw_max      = 0.0 ! Max. sapwood biomass, kg C/individual

  ! TODO: figure out how to do starvation mortality without hard-coded assumption
  !       of the mortality call time step
  real    :: BM_ys        = 0.0 ! bwood + bsw at the end the previous year, for starvation
                                ! mortality calculation.
  real    :: DBH_ys

! Adam Wolf
  real    :: psi_r  = 0.0 ! psi of root (root-stem interface), m of water
  real    :: psi_x  = 0.0 ! psi of xylem (stem-leaf interface), Pa
  real    :: psi_l  = 0.0 ! psi of leaf (leaf-substomatal cavity interface), Pa
  real    :: Kxa    = 0.0 ! conductivity of stem kg/m2 swa /s /(Pa/m m height)
  real    :: Kla    = 0.0 ! conductivity of leaf kg/m2 leaf /s /Pa
  real    :: Kxi    = 0.0 ! conductivity of stem kg/indiv /s /(Pa/m Pa)
  real    :: Kli    = 0.0 ! conductivity of leaf kg/indiv /s /Pa

  real    :: w_scale = 1.0 ! water stress reduction of stomatal conductance, unitless,
                          ! for diagnostics only

! ---- uptake-related variables
  real    :: br_profile(max_lev)  = 0.0 ! normalized vertical distribution of root biomass
  real    :: root_length(max_lev) = 0.0 ! individual root length per unit depth, m of root/m
  real    :: K_r = 0.0 ! root membrane permeability per unit area, kg/(m3 s)
  real    :: r_r = 0.0 ! radius of fine roots, m
  real    :: uptake_frac(max_lev) = 0.0 ! normalized vertical distribution of uptake

! ---- auxiliary variables
  real    :: Wl_max  = 0.0 ! maximum liquid water content of canopy, kg/individual
  real    :: Ws_max  = 0.0 ! maximum solid water content of canopy, kg/individual
  real    :: mcv_dry = 0.0 ! heat capacity of dry canopy J/(K individual)
  real    :: cover

  real :: An_op = 0.0 ! mol C/(m2 of leaf per year)
  real :: An_cl = 0.0 ! mol C/(m2 of leaf per year)

  real :: carbon_gain = 0.0 ! carbon gain since last growth, kg C/individual
  real :: carbon_loss = 0.0 ! carbon loss since last growth, kg C/individual
  real :: bwood_gain  = 0.0 ! wood gain since last growth, kg C/individual

  real :: nitrogen_gain = 0.0 ! nitrogen gain since last growth, kg N/individual
  real :: nitrogen_loss = 0.0 ! nitrogen loss since last growth, kg N/individual

  ! used in fast time scale calculations
  real :: npp_previous_day     = 0.0
  real :: npp_previous_day_tmp = 0.0

  ! new allocation fractions, Jan2 03
  real :: Pl = 0.0          ! fraction of living biomass in leaves
  real :: Pr = 0.0          ! fraction of living biomass in fine roots
  real :: Psw= 0.0          ! fraction of living biomass in sapwood
  real :: Psw_alphasw = 0.0 ! fraction of sapwood times
                            ! retirement rate of sapwood into wood
  real :: extinct = 0.0     ! light extinction coefficient in the canopy for photosynthesis calculations

  ! ens introduce for growth respiration
  real :: growth_previous_day     = 0.0 ! kgC/individual, pool of growth respiration
  real :: growth_previous_day_tmp = 0.0 ! kgC/individual per year, rate of release of
                                        ! growth respiration to the atmosphere
  !02/08/17
  real :: brsw = 0.0

  ! for phenology
  real :: gdd = 0.0
  real :: theta_av_phen = 0.0 ! relative soil moisture availability for phenology
  real :: psist_av_phen = 0.0 ! water stress for phenology

  ! for Light Saber
  real :: An_newleaf_daily = 0.0

  ! BNS: Maximum leaf biomass, to be used in the context of nitrogen limitation as suggested by Elena
  ! This will be either fixed or calculated as a function of nitrogen uptake or availability
  real :: nitrogen_stress = 0.0

  real :: scav_C = 0.0, scav_N = 0.0 ! Biomass of "scavenger" type mycorrhizae (corresponding to Arbuscular mycorrhizae)
  real :: mine_C = 0.0, mine_N = 0.0 ! Biomass of "miner" type mycorrhizae
  real :: nfix_C = 0.0, nfix_N = 0.0 ! Biomass of symbiotic N fixing microbes

  ! C and N reservoirs of symbiotic microbes (used for growth and transfers to plants)
  real :: scav_C_reservoir = 0.0, scav_N_reservoir = 0.0
  real :: mine_C_reservoir = 0.0, mine_N_reservoir = 0.0
  real :: nfix_C_reservoir = 0.0, nfix_N_reservoir = 0.0

  ! Allows cohorts to remember previous marginal gains so shifts in N allocation
  ! can have time scales
  real :: scav_mgain_smoothed = 0.0
  real :: mine_mgain_smoothed = 0.0
  real :: nfix_mgain_smoothed = 0.0
  real :: exud_mgain_smoothed = 0.0

  real :: max_monthly_scav_alloc = 0.0
  real :: max_monthly_mine_alloc = 0.0
  real :: max_monthly_Nfix_alloc = 0.0
  real :: scav_alloc_accum = 0.0
  real :: mine_alloc_accum = 0.0
  real :: Nfix_alloc_accum = 0.0
  real :: max_scav_alloc = 0.0
  real :: max_mine_alloc = 0.0
  real :: max_Nfix_alloc = 0.0
  real :: nitrogen_stress_smoothed = 1.0
end type vegn_cohort_type

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
! returns biomass of an individual of given cohort, kgC
pure elemental real function plant_C(cc)
  type(vegn_cohort_type), intent(in) :: cc
  plant_C = &
          cc%bl  + cc%blv + &
          cc%br  + cc%bwood + &
          cc%bsw + cc%bseed + &
          cc%nsc + &
          cc%carbon_gain + cc%bwood_gain + &
          cc%growth_previous_day + &
          ! Mycorrhizal and N fixer biomass added by B. Sulman
          cc%scav_C + cc%mine_C + cc%nfix_C + &
          cc%scav_C_reservoir + cc%mine_C_reservoir + cc%nfix_C_reservoir
end function plant_C

! ============================================================================
! returns nitrogen of an individual of given cohort, kgN
pure elemental real function plant_N(cc)
  type(vegn_cohort_type), intent(in) :: cc
  plant_N = &
         cc%root_N + cc%sapwood_N + &
         cc%wood_N + cc%leaf_N + cc%seed_N + &
         cc%stored_N + &
         ! Symbionts are counted as part of veg, not part of soil
         cc%nfix_N + cc%mine_N + cc%scav_N + &
         cc%scav_N_reservoir + cc%mine_N_reservoir + cc%nfix_N_reservoir
end function plant_N

! ============================================================================
! calculates functional dependence of wet canopy function f = x**p and its
! derivative, but approximates it with linear function in the vicinity
! of zero, to make sure that derivative does not become infinite
subroutine wet_frac(w, w_max, p, eps, f, DfDw)
  real, intent(in) :: &
       w, &     ! water content
       w_max, & ! maximum storage capacity
       p, &     ! exponent of the function
       eps      ! neighborhood of zero where we approximate x**p with linear function
  real, intent(out) :: &
       f, & ! function value
       DfDw ! function derivative w.r.t. water content

  if ( w > w_max ) then
     f = 1; DfDw = 0;
  else if ( w < 0 ) then
     f = 0; DfDw = 0;
  else
     if ( w/w_max <= eps ) then
        f = eps**(p-1)*w/w_max; DfDw = eps**(p-1)/w_max
     else
        f = (w/w_max)**p;   DfDw = p/w_max*(w/w_max)**(p-1)
     endif
  endif
end subroutine wet_frac

! ============================================================================
! given amount of water and snow, returns combined fraction of covered canopy
subroutine get_vegn_wet_frac (cohort, &
     fw, DfwDwl, DfwDws, &
     fs, DfsDwl, DfsDws  )
  type(vegn_cohort_type), intent(in)  :: cohort
  real, intent(out), optional :: &
       fw, DfwDwl, DfwDws, & ! water-covered fraction of canopy and its derivatives
       fs, DfsDwl, DfsDws    ! snow-covered fraction of canopy and its derivatives

  ! ---- local vars
  integer :: sp  ! shorthand for current cohort species
  real    :: fw0 ! total water-covered fraction (without overlap)
  real    :: &   ! local variable for fractions and derivatives
       fwL, DfwDwlL, DfwDwsL, & ! water-covered fraction of canopy and its derivatives
       fsL, DfsDwlL, DfsDwsL    ! snow-covered fraction of canopy and its derivatives

  sp = cohort%species

  ! snow-covered fraction
  if(cohort%Ws_max > 0) then
     call wet_frac(cohort%Ws, cohort%Ws_max, spdata(sp)%csc_pow, cmc_eps, fsL,  DfsDwsL)
  else
     fsL = 0.0; DfsDwsL=0.0
  endif
  DfsDwlL = 0

  ! wet fraction
  if(cohort%Wl_max > 0) then
     call wet_frac(cohort%Wl, cohort%Wl_max, spdata(sp)%cmc_pow, cmc_eps, fw0, DfwDwlL)
  else
     fw0 = 0.0; DfwDwlL=0.0
  endif
  ! take into account overlap by snow
  fwL     = fw0*(1-fsL)
  DfwDwlL = DfwDwlL*(1-fsL)
  DfwDwsL = -fw0*DfsDwsL

  ! assign result to output parameters, if present
  if (present(fw))     fw = fwL
  if (present(DfwDwl)) DfwDwl = DfwDwlL
  if (present(DfwDws)) DfwDws = DfwDwsL
  if (present(fs))     fs = fsL
  if (present(DfsDwl)) DfsDwl = DfsDwlL
  if (present(DfsDws)) DfsDws = DfsDwsL

end subroutine get_vegn_wet_frac


! ============================================================================
! given cohort and snow depth, calculates the cover (1-fraction of gaps), and
! snow factor for radiation
subroutine vegn_data_cover ( cohort, snow_depth, vegn_cover, &
                                         vegn_cover_snow_factor )
  type(vegn_cohort_type), intent(inout)  :: cohort
  real, intent(in)  :: snow_depth
  real, intent(out), optional :: vegn_cover
  real, intent(out), optional :: vegn_cover_snow_factor

  real :: f

  if(sai_cover) then
     cohort%cover = 1 - exp(-max(cohort%lai, cohort%sai))
  else
     cohort%cover = 1 - exp(-cohort%lai)
  endif

  select case(snow_masking_option)
  case (SNOW_MASKING_MCM)
     f = 1 - min(1.0, 0.5*sqrt(max(snow_depth,0.)/cohort%snow_crit))
  case (SNOW_MASKING_LM3, SNOW_MASKING_HEIGHT)
     f = cohort%snow_crit / (max(snow_depth,0.0) + cohort%snow_crit)
  case (SNOW_MASKING_NONE)
     f = 1.0
  end select
  cohort%cover = cohort%cover * f

  if (present(vegn_cover)) vegn_cover = cohort%cover
  if (present(vegn_cover_snow_factor)) vegn_cover_snow_factor = f
end subroutine vegn_data_cover


! ============================================================================
! returns properties of the fine roots
subroutine update_cohort_root_properties(soil, cohort)
  type(soil_tile_type),   intent(in)     :: soil
  type(vegn_cohort_type), intent(inout)  :: cohort
  ! note that in LM3, when the density of individuals per m2 is imposed to be 1,
  ! "per unit depth" measures are the same as volumetric density measures.

  integer :: l
  real :: factor, z
  real :: vbr ! density of fine roots biomass per unit depth, kg C/m

  associate(sp => spdata(cohort%species))
  cohort%br_profile(:) = 0.0
  z = 0
  do l = 1, num_l
     if (z+dz(l)/2>permafrost_depth_thresh.and.soil%frozen_freq(l)>permafrost_freq_thresh) exit ! from loop
     ! so that the rest of profile remains zero.

     cohort%br_profile(l) = exp(-z/cohort%root_zeta) - exp(-(z+dz(l))/cohort%root_zeta)
     z = z + dz(l)
  enddo

  factor = 1.0/sum(cohort%br_profile)
  cohort%br_profile(:) = cohort%br_profile(:)*factor
  do l = 1, num_l
     ! calculate the vertical fine root biomass density [kgC/m] for current layer
     ! NOTE: sum(brv*dz) must be equal to cohort%br, which is achieved by normalizing
     ! factor
     vbr = cohort%br * cohort%br_profile(l)*factor/dz(l)
     ! calculate fine root length per unit depth
     cohort%root_length(l) = vbr*sp%srl
  enddo

  cohort%K_r = sp%root_perm
  cohort%r_r = sp%root_r
  end associate
end subroutine update_cohort_root_properties


! ============================================================================
! calculates vertical distribution of active roots: given layer thicknesses,
! returns fraction of active roots per level
subroutine cohort_uptake_profile(cohort, dz, uptake_frac_max, vegn_uptake_term)
  type(vegn_cohort_type), intent(in)  :: cohort
  real, intent(in)  :: dz(:)
  real, intent(out) :: uptake_frac_max(:)
  real, intent(out) :: vegn_uptake_term(:)

  real, parameter :: res_scaler = mol_air/mol_h2o  ! scaling factor for water supply
  ! NOTE: there is an inconsistency there between the
  ! units of stomatal conductance [mol/(m2 s)], and the units of humidity deficit [kg/kg],
  ! in the calculations of water demand. Since the uptake options other than LINEAR cannot
  ! use res_scaler, in this code the units of humidity deficit are converted to mol/mol,
  ! and the additional factor is introduced in res_scaler to ensure that the LINEAR uptake
  ! gives the same results.

  integer :: l
  real    :: z, sum_rf

  if (use_bucket) then
     uptake_frac_max(1) = dz(1)
     z = dz(1)
     do l = 2, size(dz)
        if (cohort%root_density*(exp(-z/cohort%root_zeta)-&
            exp(-(z+dz(l))/cohort%root_zeta))/dz(l) > critical_root_density) &
        then
           uptake_frac_max(l) = dz(l)
        else
           uptake_frac_max(l) = 0
        endif
        z = z + dz(l)
     enddo
  else
     ! linear scaling, LM3V
     uptake_frac_max = cohort%br_profile(1:size(uptake_frac_max))
  endif

  sum_rf = sum(uptake_frac_max)
  if(sum_rf>0) &
       uptake_frac_max(:) = uptake_frac_max(:)/sum_rf

  if (cohort%br <= 0) then
     vegn_uptake_term(:) = 0.0
  else
     vegn_uptake_term(:) = uptake_frac_max(:) * &
          res_scaler * spdata(cohort%species)%dfr * cohort%br
  endif

end subroutine cohort_uptake_profile


! ============================================================================
! returns normalizes vertical profile for root litter
subroutine cohort_root_litter_profile(cohort, dz, profile)
  type(vegn_cohort_type), intent(in) :: cohort
  real, intent(in)  :: dz(:) ! layer thickneses, m
  real, intent(out) :: profile(:)

  profile = cohort%br_profile(1:size(profile))
end subroutine cohort_root_litter_profile

! ============================================================================
! returns normalizes vertical profile for root exudates
subroutine cohort_root_exudate_profile(cohort, dz, profile)
  type(vegn_cohort_type), intent(in) :: cohort
  real, intent(in)  :: dz(:) ! layer thickneses, m
  real, intent(out) :: profile(:)

  call cohort_root_litter_profile(cohort,dz,profile)
end subroutine cohort_root_exudate_profile

! ============================================================================
function btotal(c)
  real :: btotal ! returned value
  type(vegn_cohort_type), intent(in) :: c

  btotal = c%bliving+c%bwood
end function btotal


! ============================================================================
function c3c4(c, temp, precip) result (pt)
  integer :: pt
  type(vegn_cohort_type), intent(in) :: c
  real,              intent(in) :: temp   ! temperature, degK
  real,              intent(in) :: precip ! precipitation, ???

  real :: pc4

  ! Rule based on analysis of ED global output; equations from JPC, 2/02
  if(btotal(c) < tg_c4_thresh) then
    pc4=exp(-0.0421*(273.16+25.56-temp)-(0.000048*(273.16+25.5-temp)*precip));
  else
    pc4=0.0;
  endif

  if(pc4>0.5) then
    pt=PT_C4
  else
    pt=PT_C3
  endif

end function c3c4


! ============================================================================
! given current conditions, returns type of phenology.
function phenology_type(c, cm)
  integer :: phenology_type
  type(vegn_cohort_type), intent(in) :: c  ! cohort (not used)
  real, intent(in) :: cm ! number of cold months

  real :: pe  ! prob evergreen

  ! GCH, Rule based on analysis of ED global output; equations from JPC, 2/02
  ! GCH, Parameters updated 2/9/02 from JPC
  pe = 1.0/(1.0+((1.0/0.00144)*exp(-0.7491*cm)));

  if(pe>phen_ev1 .and. pe<phen_ev2) then
     phenology_type = PHEN_EVERGREEN ! its evergreen
  else
     phenology_type = PHEN_DECIDUOUS ! its deciduous
  endif
end function phenology_type


! ============================================================================
! given a cohort, climatology, and land use type, determines and updates
! physiology type, phenology type, and species of the cohort
subroutine update_species(c, t_ann, t_cold, p_ann, cm, landuse)
  type(vegn_cohort_type), intent(inout) :: c    ! cohort to update
  real,              intent(in) :: t_ann   ! annual-mean temperature, degK
  real,              intent(in) :: t_cold  ! average temperature of the coldest month, degK
  real,              intent(in) :: p_ann   ! annual-mean precipitation, mm/yr
  real,              intent(in) :: cm      ! number of cold months
  integer,           intent(in) :: landuse ! land use type

  integer :: spp
  integer :: pt    ! physiology type
  integer :: phent ! type of phenology

  pt    = c3c4(c,t_ann,p_ann)
  phent = phenology_type(c, cm)

  if(landuse == LU_CROP) phent = PHEN_DECIDUOUS ! crops cannot be evergreen

  if(pt==PT_C4) then
     spp=SP_C4GRASS;  ! c4 grass
  else if(phent==1) then
     spp=SP_EVERGR;   ! evergreen non-grass
  else if(btotal(c) < tg_c3_thresh) then
     spp=SP_C3GRASS;  ! c3 grass
  else if ( t_cold > T_cold_tropical ) then  ! ens,slm Jun 21 2003 to prohibit tropical forest in coastal cells
     spp=SP_TROPICAL; ! tropical deciduous non-grass
  else
     spp=SP_TEMPDEC;  ! temperate deciduous non-grass
  endif

  ! reset leaf age to zero if species are changed
  if (spp/=c%species) c%leaf_age = 0.0

  c%species = spp
end subroutine update_species


! ============================================================================
function leaf_area_from_biomass(bl,species,layer,firstlayer) result (area)
  real :: area ! returned value
  real,    intent(in) :: bl      ! biomass of leaves, kg C/individual
  integer, intent(in) :: species ! species
  integer, intent(in) :: layer, firstlayer

  if (layer > 1 .AND. firstlayer == 0) then
     area = bl/(spdata(species)%LMA_understory_factor*spdata(species)%LMA) ! half thickness for leaves in understory
  else
     area = bl/spdata(species)%LMA
  endif
end function leaf_area_from_biomass


! ============================================================================
function height_from_biomass(btotal) result (height)
    real :: height ! return value
    real, intent(in) :: btotal ! total biomass

    height = 24.19*(1.0-exp(-0.19*btotal))
end function height_from_biomass


! ============================================================================
! calculates fractions of living biomass in different compartments
subroutine update_bio_living_fraction(c)
  type(vegn_cohort_type), intent(inout) :: c

  real    :: D  ! inverse denominator
  integer :: sp ! species, for convenience

  sp     = c%species
  D      = 1/(1 + spdata(sp)%c1 + c%height*spdata(sp)%c2)

  c%Pl   = D
  c%Pr   = spdata(sp)%c1 * D
  c%Psw  = 1 - c%Pl - c%Pr
  c%Psw_alphasw = spdata(sp)%c3 * spdata(sp)%alpha_leaf * D

end subroutine update_bio_living_fraction


! ============================================================================
! redistribute living biomass pools in a given cohort, and update related
! properties (height, lai, sai)
subroutine update_biomass_pools(c)
  type(vegn_cohort_type), intent(inout) :: c

  real :: biomass_N_demand  ! Live biomass in excess of max (used in N limitation system) -- BNS
  real :: potential_stored_N ! N storage if there is no N-caused change in biomass allocation
  real :: available_N,N_demand,extra_C

  if (do_ppa) return ! in PPA mode, do nothing at all

  c%height  = height_from_biomass(btotal(c))
  c%total_N = c%stored_N+c%leaf_N+c%wood_N+c%root_N+c%sapwood_N

  call update_bio_living_fraction(c)


  ! Stress increases as stored N declines relative to total biomass N demand
  ! N stress is calculated based on "potential pools" without N limitation
  ! biomass_N_demand=(c%bliving*c%Pl/leaf_live_c2n + c%bliving*c%Pr/froot_live_c2n + c%bliving*c%Psw/wood_fast_c2n)
  ! Elena suggests using 2*(root N + leaf N) as storage target
  ! bliving is already increased after cgain
  biomass_N_demand=(c%bliving*c%Pl/spdata(c%species)%leaf_live_c2n + c%bliving*c%Pr/spdata(c%species)%froot_live_c2n)
  if(c%status == LEAF_OFF) then
    potential_stored_N = c%total_N - biomass_N_demand*2 - c%bwood/spdata(c%species)%wood_c2n-c%bliving*c%Psw/spdata(c%species)%sapwood_c2n
  else
    potential_stored_N = c%total_N - biomass_N_demand - c%bwood/spdata(c%species)%wood_c2n-c%bliving*c%Psw/spdata(c%species)%sapwood_c2n
  endif
  ! c%nitrogen_stress = biomass_N_demand/c%total_N

  ! Spring physical analogy -- restoring force proportional to distance from target (equal to demand*2.0)
  ! Leaving spring constant 1.0 for now
  ! Stress is normalized by N demand so it is an index that does not depend on total biomass

  ! if (c%total_N>0.0) then
  !  c%nitrogen_stress = -1.0 * ((potential_stored_N) - 2.0*biomass_N_demand)/abs(c%total_N)
    c%nitrogen_stress = (2.0*biomass_N_demand - potential_stored_N)/abs(biomass_N_demand)
  ! else
  !   c%nitrogen_stress = 0.0
  ! endif


  ! if(c%total_N > biomass_N_demand) then
  !   c%nitrogen_stress = (biomass_N_demand/(c%total_N-biomass_N_demand-c%bwood/wood_fast_c2n))**1 ! This is demand/storage, constrained to be >0
  !
  ! else
  !   c%nitrogen_stress = 5.0
  ! endif
  c%nitrogen_stress = min(c%nitrogen_stress,5.0)
  c%nitrogen_stress = max(c%nitrogen_stress,0.05)


  ! Move some biomass from sapwood to root based on N stress
  c%Pr=c%Pr+c%Psw*c%nitrogen_stress*spdata(c%species)%N_stress_root_factor
  c%Psw=c%Psw-c%Psw*c%nitrogen_stress*spdata(c%species)%N_stress_root_factor

  c%bsw = c%Psw*c%bliving;
  if(c%status == LEAF_OFF) then
     c%blv = c%Pl*c%bliving + c%Pr*c%bliving;
     c%bl  = 0;
     c%br  = 0;
    !  c%nitrogen_stress = 0
  else
     c%blv = 0;
     c%bl  = c%Pl*c%bliving;
     c%br  = c%Pr*c%bliving;
  endif


  if(N_limits_live_biomass) then

    if(c%total_N<0.0) then
      __DEBUG4__(c%leaf_N,c%wood_N,c%root_N,c%sapwood_N)
      __DEBUG3__(c%status,c%stored_N,c%total_N)
      call error_mesg('update_biomass_pools','totalN<0',FATAL)
    endif
      available_N = max(0.0,(c%total_N - c%bwood/spdata(c%species)%wood_c2n))
      N_demand = &
          (c%bl/spdata(c%species)%leaf_live_c2n+c%bsw/spdata(c%species)%sapwood_c2n+c%br/spdata(c%species)%froot_live_c2n)

      if(N_demand>available_N) then
        extra_C = (c%bl+c%br+c%bsw)*(1.0-available_N/N_demand)
        ! __DEBUG3__(available_N/N_demand,extra_C,(c%bl+c%br+c%bsw)*(available_N/N_demand))
        c%bl=c%bl*available_N/N_demand
        c%br=c%br*available_N/N_demand
        c%bsw=c%bsw*available_N/N_demand
        c%blv=c%blv+extra_C
        c%Pl=c%Pl*available_N/N_demand
        c%Pr=c%Pr*available_N/N_demand
        c%Psw=c%Psw*available_N/N_demand
        if(c%bsw<0) then
          __DEBUG3__(c%Pl,c%Pr,c%Psw)
          __DEBUG5__(N_demand,available_N,c%stored_N,c%total_N,c%bwood/spdata(c%species)%wood_c2n)
          call error_mesg('update_biomass_pools','bsw<0',FATAL)
        endif
      endif
    endif


  c%leaf_N=c%bl/spdata(c%species)%leaf_live_c2n
  c%wood_N=c%bwood/spdata(c%species)%wood_c2n
  c%sapwood_N=c%bsw/spdata(c%species)%sapwood_c2n
  c%root_N=c%br/spdata(c%species)%froot_live_c2n
  c%stored_N=c%total_N-(c%leaf_N+c%wood_N+c%root_N+c%sapwood_N)

  ! In rare cases wood N content growth could cause stored N to go below zero
  ! In that case, make sure stored N is zero by reducing wood N
  if (c%stored_N<0 .AND. N_limits_live_biomass) then
     c%wood_N=c%wood_N+c%stored_N
     c%stored_N=0.0
  endif

end subroutine update_biomass_pools


! ==============================================================================
subroutine init_cohort_allometry_ppa(cc, height, nsc_frac, nsn_frac)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: height
  real, intent(in) :: nsc_frac ! amount of nsc, as a fraction of bl_max
                               ! NOTE: typically nsc_frac is greater than 1
  real, intent(in) :: nsn_frac ! nonstructural N as fraction of bl_max N content

  real :: bw ! biomass of wood+sapwood, based on allometry and height
  real :: Dwood ! diameter of heartwood
  real :: BL_c, BL_u ! max leaf bimass for canopy and understory, respectively

  associate(sp=>spdata(cc%species))
  cc%height = height
  select case (sp%allomt)
  case (ALLOM_EW,ALLOM_EW1)
     cc%dbh = (height/sp%alphaHT)**(1.0/sp%thetaHT)
     bw = sp%rho_wood * sp%alphaBM * cc%dbh**sp%thetaBM
     Dwood = sqrt(max(0.0,cc%dbh**2 - 4/pi*sp%alphaCSASW * cc%dbh**sp%thetaCSASW))
     cc%bwood = sp%rho_wood * sp%alphaBM * cc%dbh**sp%thetaBM
  case (ALLOM_HML)
     cc%dbh = (sp%gammaHT/(sp%alphaHT/height - 1))**(1.0/sp%thetaHT)
     bw = sp%rho_wood * sp%alphaBM * cc%dbh**2 * cc%height
     ! isa and es 201701 - CSAsw different for ALLOM_HML
     !Dwood = sqrt(max(0.0,cc%dbh**2 - 4/pi*sp%alphaCSASW * cc%dbh**sp%thetaCSASW))
     Dwood = sqrt(max(0.0,cc%dbh**2 - 4/pi* sp%phiCSA * &
                  cc%DBH**(sp%thetaCA + sp%thetaHT) / (sp%gammaHT + cc%DBH** sp%thetaHT)))
     cc%bwood = sp%rho_wood * sp%alphaBM * Dwood**2 * cc%height
  end select
  if (sp%lifeform == FORM_GRASS) then ! isa 20170705 - grasses do not form heartwood
     cc%bwood = 0.0
  endif
  cc%bsw   = bw - cc%bwood
  ! update derived quantyties based on the allometry
  cc%crownarea = sp%alphaCA * cc%dbh**sp%thetaCA
  cc%bl      = 0.0
  cc%br      = 0.0
  !02/08/17 ens
  cc%brsw = sp%branch_wood_frac * ( cc%bsw + cc%bwood )
  cc%brsw = min(cc%brsw, 0.6*cc%bsw)

  BL_c = sp%LMA * sp%LAImax * cc%crownarea * (1.0-sp%internal_gap_frac)
  BL_u = BL_c * sp%understory_lai_factor
  if(sp%lifeform == FORM_GRASS) then
     cc%bl_max = BL_c
     cc%bsw_max = sp%rho_wood * sp%alphaBM * cc%height * cc%dbh**2
  else
     cc%bl_max = BL_u + min(cc%topyear/5.0,1.0)*(BL_c - BL_u)
     cc%bsw_max = 0.0
  endif
  cc%br_max  = sp%phiRL*cc%bl_max/(sp%LMA*sp%SRA)
  cc%nsc     = nsc_frac * cc%bl_max
  cc%blv     = 0.0
  cc%bseed   = 0.0
  cc%bliving = cc%br + cc%bl + cc%bsw + cc%blv

  if (soil_carbon_option==SOILC_CORPSE_N) then
    cc%stored_N   = nsn_frac * cc%bl_max/sp%leaf_live_c2n
    cc%seed_N     = 0.0
    cc%wood_N     = cc%bwood/sp%wood_c2n
    cc%sapwood_N  = cc%bsw/sp%sapwood_c2n
    cc%leaf_N     = (cc%bl+cc%blv)/sp%leaf_live_c2n
    cc%root_N     = cc%br/sp%froot_live_c2n
    cc%total_N    = cc%stored_N+cc%leaf_N+cc%wood_N+cc%root_N+cc%sapwood_N+cc%seed_N
  else
    cc%stored_N   = 0.0
    cc%seed_N     = 0.0
    cc%wood_N     = 0.0
    cc%sapwood_N  = 0.0
    cc%leaf_N     = 0.0
    cc%root_N     = 0.0
    cc%total_N    = 0.0
  endif

  ! Should these have nonzero initial values on initialization/reproduction?
  cc%scav_C = 0.0 ; cc%scav_N = 0.0
  cc%mine_C = 0.0 ; cc%mine_N = 0.0
  cc%nfix_C = 0.0 ; cc%nfix_N = 0.0

  cc%scav_C_reservoir = 0.0 ; cc%scav_N_reservoir = 0.0
  cc%mine_C_reservoir = 0.0 ; cc%mine_N_reservoir = 0.0
  cc%nfix_C_reservoir = 0.0 ; cc%nfix_N_reservoir = 0.0

  cc%growth_previous_day     = 0.0
  cc%growth_previous_day_tmp = 0.0
  cc%carbon_gain = 0.0
  cc%carbon_loss = 0.0
  cc%DBH_ys      = cc%DBH
  cc%BM_ys       = cc%bsw + cc%bwood
  cc%npp_previous_day     = 0.0
  cc%npp_previous_day_tmp = 0.0
  cc%leaf_age    = 0.0

  if (is_watch_point()) then
     write(*,*) '#### init_cohort_allometry_ppa output ####'
     __DEBUG2__(cc%height, cc%dbh)
     __DEBUG5__(cc%bl, cc%br, cc%bsw, cc%bwood, cc%nsc)
     __DEBUG5__(cc%leaf_N, cc%root_N, cc%sapwood_N, cc%wood_N, cc%stored_N)
     __DEBUG3__(cc%bseed, cc%seed_N, cc%blv)
     __DEBUG2__(cc%bl_max, cc%br_max)
  endif
  end associate
end subroutine init_cohort_allometry_ppa

! ==============================================================================
! adam wolf
! Stored in cc as K per tissue area.
!  spdata stores pressure units in m
!  Potential in units m to be consistent with Sergei's (conversion at namelist read)
subroutine init_cohort_hydraulics(cc, init_psi)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: init_psi

!     TODO: rootarea, leafarea, stemarea are calculated independently in init_cohort_hydraulics.
!     Is it consistent with the rest of the code? Is there any way to avoid this calculations, for consistency?
!     perhaps move Kri, Kxi, Kli initialization/update to update_derived_vegn_data?
  real :: rootarea, stemarea

  associate(sp=>spdata(cc%species))
  rootarea = cc%br * sp%srl * 2*PI * sp%root_r  ! (kg/indiv)(m/kg)(m2/m)

  ! isa and es 201701 - CSAsw different for ALLOM_HML
  select case (sp%allomt)
  case (ALLOM_EW,ALLOM_EW1)
     stemarea = sp%alphaCSASW * cc%DBH**sp%thetaCSASW
  case (ALLOM_HML)
     stemarea = sp%phiCSA * cc%DBH**(sp%thetaCA + sp%thetaHT) / (sp%gammaHT + cc%DBH** sp%thetaHT)
  end select
  ! isa 20170715
  if(sp%lifeform == FORM_GRASS) then
    stemarea = PI * cc%DBH * cc%DBH / 4.0 ! trunk cross-sectional area
  endif

  cc%Kxa = sp%Kxam
  cc%Kla = sp%Klam

  cc%Kxi = cc%Kxa * stemarea / cc%height
  cc%Kli = cc%Kla * cc%leafarea

  !write (*,*) ' Kli  ', cc%Kli

  cc%psi_r  = init_psi
  cc%psi_x  = init_psi
  cc%psi_l  = init_psi

  end associate
end subroutine init_cohort_hydraulics

! ============================================================================
function cohorts_can_be_merged(c1,c2); logical cohorts_can_be_merged
   type(vegn_cohort_type), intent(in) :: c1,c2

   logical :: sameSpecies, sameLayer, sameSize

   sameSpecies = c1%species == c2%species
   sameLayer   = (c1%layer == c2%layer) .and. (c1%firstlayer == c2%firstlayer)
   sameSize    = (abs(c1%DBH - c2%DBH) <= DBH_merge_rel*max(abs(c1%DBH),abs(c2%DBH))) .or.  &
                 (abs(c1%DBH - c2%DBH) <  DBH_merge_abs)
   if (spdata(c1%species)%lifeform == FORM_GRASS) then
      sameSize = sameSize .and. (abs(c1%nsc - c2%nsc) <= NSC_merge_rel*max(abs(c1%nsc), abs(c2%nsc)))
   endif

   cohorts_can_be_merged = &
        sameSpecies .and. sameLayer .and. sameSize
end function cohorts_can_be_merged

! ============================================================================
! returns TRUE if cohort can allocate carbon to seeds.
logical function cohort_makes_seeds(cc, G_WF) result(answer)
  type(vegn_cohort_type), intent(in) :: cc
  real, intent(in) :: G_WF ! amount of carbon spent on new sapwood growth and seeds

  answer = (cc%layer == 1.or.spdata(cc%species)%reproduces_in_understory) &
           .and. cc%age > spdata(cc%species)%maturalage
  if (soil_carbon_option==SOILC_CORPSE_N.AND.N_limits_live_biomass) then
     answer = answer .AND. &
        .NOT.(cc%nitrogen_stress > spdata(cc%species)%max_n_stress_for_seed_production &
              .OR. spdata(cc%species)%v_seed*G_WF/spdata(cc%species)%seed_c2n>0.1*cc%stored_N )
  endif
end function cohort_makes_seeds

! ============================================================================
logical function cohort_can_reproduce(cc)
  type(vegn_cohort_type), intent(in) :: cc

  cohort_can_reproduce = (cc%layer==1.or.spdata(cc%species)%reproduces_in_understory) &
                         .and. cc%age > spdata(cc%species)%maturalage
end function cohort_can_reproduce

end module vegn_cohort_mod
