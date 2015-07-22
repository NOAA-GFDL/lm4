module vegn_cohort_mod

use constants_mod, only: PI

use land_constants_mod, only: NBANDS, mol_h2o, mol_air
use vegn_data_mod, only : spdata, &
   use_mcm_masking, use_bucket, critical_root_density, &
   tg_c4_thresh, tg_c3_thresh, l_fract, fsc_liv, &
   phen_ev1, phen_ev2, cmc_eps
use vegn_data_mod, only : PT_C3, PT_C4, CMPT_ROOT, CMPT_LEAF, &
   SP_C4GRASS, SP_C3GRASS, SP_TEMPDEC, SP_TROPICAL, SP_EVERGR, &
   LEAF_OFF, LU_CROP, PHEN_EVERGREEN, PHEN_DECIDIOUS, &
   do_ppa
use soil_tile_mod, only : max_lev

implicit none
private
! ==== public interfaces =====================================================
public :: vegn_cohort_type

! operations defined for cohorts
public :: get_vegn_wet_frac
public :: vegn_data_cover
public :: cohort_uptake_profile
public :: cohort_root_properties
 
public :: btotal ! returns cohort total biomass
public :: c3c4   ! returns physiology type for given biomasses and conditions
public :: phenology_type ! returns type of phenology for given conditions
public :: update_species ! updates cohort physiology, phenology type, and species
public :: leaf_area_from_biomass ! given leaf biomass, calculates leaf area
public :: height_from_biomass    ! given total biomass, calculated tree height
public :: update_bio_living_fraction
public :: update_biomass_pools
public :: init_cohort_allometry_ppa
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id$', &
     tagname = '$Name$'

! ==== types =================================================================
type :: vegn_phys_prog_type
end type vegn_phys_prog_type

! vegn_cohort_type describes the data that belong to a vegetation cohort
type :: vegn_cohort_type
  real :: Wl ! liquid water content of canopy, kg/individual
  real :: Ws ! solid water (snow) content of canopy, kg/individual
  real :: Tv ! canopy temperature, K

! ---- biological prognostic variables
! Currently bio prognostic variable is defined as anything that's saved in
! the restart; clearly some vars there are, strictly speaking, diagnostic,
! but saved for reproducibility (to avoid recalculation). exceptions :
! npp_previous_day is left outside, since it's obviously auxiliary; height
! is left outside
  integer :: species = 0   ! vegetation species
  real    :: bl      = 0.0 ! biomass of leaves, kg C/individual
  real    :: blv     = 0.0 ! biomass of virtual leaves (labile store), kg C/individual
  real    :: br      = 0.0 ! biomass of fine roots, kg C/individual
  real    :: bsw     = 0.0 ! biomass of sapwood, kg C/individual
  real    :: bwood   = 0.0 ! biomass of heartwood, kg C/individual
  real    :: bseed   = 0.0 ! biomass put aside for future progeny, kg C/individual
  real    :: nsc     = 0.0 ! non-structural carbon, kg C/individual

  real    :: bliving = 0.0 ! leaves, fine roots, and sapwood biomass
  integer :: status  = 0   ! growth status of plant
  real    :: leaf_age= 0.0 ! age of leaf in days since budburst

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
  
  ! TODO: figure out how to do starvation mortality without hard-coded assumption
  !       of the mortality call time step
  real    :: BM_ys        = 0.0 ! bwood + bsw at the end the previous year, for starvation 
                                ! mortality calculation.
  real    :: DBH_ys

! ---- uptake-related variables
  real    :: root_length(max_lev) = 0.0 ! individual's root length per unit depth, m of root/m
  real    :: K_r = 0.0 ! root membrane permeability per unit area, kg/(m3 s)
  real    :: r_r = 0.0 ! radius of fine roots, m
  real    :: uptake_frac(max_lev) = 0.0 ! normalized vertical distribution of uptake

! ---- auxiliary variables 
  real    :: Wl_max  = 0.0 ! maximum liquid water content of canopy, kg/individual
  real    :: Ws_max  = 0.0 ! maximum solid water content of canopy, kg/individual
  real    :: mcv_dry = 0.0 ! heat capacity of dry canopy J/(K individual)
  real    :: cover

  integer :: pt = 0  ! physiology type
  integer :: phent = 0 ! phenology type

  real :: gpp  = 0.0 ! gross primary productivity kg C/timestep
  real :: npp  = 0.0 ! net primary productivity kg C/timestep

  real :: resp = 0.0 ! plant respiration
  real :: resl = 0.0 ! leaf respiration
  real :: resr = 0.0 ! root respiration
  real :: resg = 0.0 ! growth respiration

  real :: An_op = 0.0 ! mol C/(m2 of leaf per year)
  real :: An_cl = 0.0 ! mol C/(m2 of leaf per year)
  
  real :: carbon_gain = 0.0 ! carbon gain since last growth, kg C/individual
  real :: carbon_loss = 0.0 ! carbon loss since last growth, kg C/individual
  real :: bwood_gain  = 0.0 ! wood gain since last growth, kg C/individual

  ! used in fast time scale calculations
  real :: npp_previous_day     = 0.0
  real :: npp_previous_day_tmp = 0.0

  ! lena added this for storing previous size stomatal opening and lwnet 
  ! for computing canopy air T and q at the next step
  
  real :: gs = 0.0
  real :: gb = 0.0

  real :: ds = 0.0

  ! new allocation fractions, Jan2 03
  real :: Pl = 0.0          ! fraction of living biomass in leaves
  real :: Pr = 0.0          ! fraction of living biomass in fine roots
  real :: Psw= 0.0          ! fraction of living biomass in sapwood
  real :: Psw_alphasw = 0.0 ! fraction of sapwood times 
                            ! retirement rate of sapwood into wood
  real :: extinct = 0.0     ! light extinction coefficient in the canopy for photosynthesis calculations
  
! for phenology
  real :: gdd = 0.0

! in LM3V the cohort structure has a handy pointer to the tile it belongs to;
! so operations on cohort can update tile-level variables. In this code, it is
! probably impossible to have this pointer here: it needs to be of type
! "type(vegn_tile_type), pointer", which means that the vegn_cohort_mod needs to
! use vegn_tile_mod. But the vegn_tile_mod itself uses vegn_cohort_mod, and 
! this would create a circular dependency of modules, something that's 
! prohibited in FORTRAN.
!  type(vegn_tile_type), pointer :: cp
end type vegn_cohort_type

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
! calculates functional dependence of wet canopy function f = x**p and its 
! derivative, but approximates it with linear function in the vicinity 
! of zero, to make sure that derivative doesn't become infinite
subroutine wet_frac(w, w_max, p, eps, f, DfDw)
  real, intent(in) :: &
       w, &     ! water content
       w_max, & ! maximum storage capacity
       p, &     ! exponent of the function
       eps      ! neighborhood of zero where we approximate x**p with linear function
  real, intent(out) :: &
       f, & ! function value
       DfDw ! it's derivative

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
  
end subroutine


! ============================================================================
! given cohort and snow depth, calculates the cover (1-fraction of gaps), and
! snow factor for radiation
! TODO: we probably need to revise the snow correction, because currently it
! does not take into account the height of the cohort, and therefore will
! be the same for all cohorts (of the same species). snow_crit doesn't
! depend on height, it's just prescribed (per species)
subroutine vegn_data_cover ( cohort, snow_depth, vegn_cover, &
                                         vegn_cover_snow_factor )
  type(vegn_cohort_type), intent(inout)  :: cohort
  real, intent(in)  :: snow_depth
  real, intent(out), optional :: vegn_cover
  real, intent(out), optional :: vegn_cover_snow_factor

  cohort%cover = 1 - exp(-cohort%lai)
  if (use_mcm_masking) then
     if (present(vegn_cover_snow_factor)) vegn_cover_snow_factor =  &
           (1 - min(1., 0.5*sqrt(max(snow_depth,0.)/cohort%snow_crit)))
     cohort%cover = cohort%cover * &
           (1 - min(1., 0.5*sqrt(max(snow_depth,0.)/cohort%snow_crit)))
  else
     if (present(vegn_cover_snow_factor)) vegn_cover_snow_factor =  &
           cohort%snow_crit / &
          (max(snow_depth,0.0) + cohort%snow_crit)
     cohort%cover = cohort%cover * &
           cohort%snow_crit / &
          (max(snow_depth,0.0) + cohort%snow_crit)
  endif
  if (present(vegn_cover)) vegn_cover = cohort%cover
end subroutine vegn_data_cover


! ============================================================================
! returns properties of the fine roots
subroutine cohort_root_properties(cohort, dz, vrl, K_r, r_r)
  type(vegn_cohort_type), intent(in)  :: cohort
  real, intent(in)  :: dz(:)
  real, intent(out) :: &
       vrl(:), & ! volumetric fine root length, m/m3
       K_r,    & ! root membrane permeability per unit area, kg/(m3 s)
       r_r       ! radius of fine roots, m

  integer :: sp, l
  real :: factor, z
  real :: vbr ! volumetric biomass of fine roots, kg C/m3

  sp = cohort%species

  factor = 1.0/(1.0-exp(-sum(dz)/cohort%root_zeta))
  z = 0
  do l = 1, size(dz)
     ! calculate the volumetric fine root biomass density [kgC/m3] for current layer
     ! NOTE: sum(brv*dz) must be equal to cohort%br, which is achieved by normalizing
     ! factor
     vbr = cohort%br * &
          (exp(-z/cohort%root_zeta) - exp(-(z+dz(l))/cohort%root_zeta))*factor/dz(l)
     ! calculate the volumetric fine root length
     vrl(l) = vbr*spdata(sp)%srl

     z = z + dz(l)
  enddo

  K_r = spdata(sp)%root_perm
  r_r = spdata(sp)%root_r

end subroutine 


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
  ! in the calculations of water demand. Since the uptake options other than LINEAR can't 
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
     !linear scaling, LM3V
     z = 0
     do l = 1, size(dz)
        uptake_frac_max(l) = (exp(-z/cohort%root_zeta)    &
                - exp(-(z+dz(l))/cohort%root_zeta))
        uptake_frac_max(l) = &
                max( uptake_frac_max(l), 0.0)
        z = z + dz(l)
     enddo

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

end subroutine 


! ============================================================================
function btotal(c)
  real :: btotal ! returned value
  type(vegn_cohort_type), intent(in) :: c
  
  btotal = c%bliving+c%bwood
end function


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
  
end function


! ============================================================================
! given current conditions, returns type of phenology.
function phenology_type(c, cm)
  integer :: phenology_type
  type(vegn_cohort_type), intent(in) :: c  ! cohort (not used???)
  real, intent(in) :: cm ! number of cold months
   
  real :: pe  ! prob evergreen
   
  ! GCH, Rule based on analysis of ED global output; equations from JPC, 2/02
  ! GCH, Parameters updated 2/9/02 from JPC
  pe = 1.0/(1.0+((1.0/0.00144)*exp(-0.7491*cm)));
  
  if(pe>phen_ev1 .and. pe<phen_ev2) then
     phenology_type = PHEN_EVERGREEN ! its evergreen
  else
     phenology_type = PHEN_DECIDIOUS ! its deciduous
  endif
end function


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

  c%pt    = c3c4(c,t_ann,p_ann)
  c%phent = phenology_type(c, cm) 
  
  if(landuse == LU_CROP) c%phent = 0  ! crops can't be evergreen
  
  if(c%pt==PT_C4) then
     spp=SP_C4GRASS;  ! c4 grass
  else if(c%phent==1) then
     spp=SP_EVERGR;   ! evergreen non-grass
  else if(btotal(c) < tg_c3_thresh) then
     spp=SP_C3GRASS;  ! c3 grass
  else if ( t_cold > 278.16 ) then  ! ens,slm Jun 21 2003 to prohibit tropical forest in coastal cells
     spp=SP_TROPICAL; ! tropical deciduous non-grass
  else 
     spp=SP_TEMPDEC;  ! temperate deciduous non-grass
  endif

  ! reset leaf age to zero if species are changed
  if (spp/=c%species) c%leaf_age = 0.0

  c%species = spp
end subroutine


! ============================================================================
function leaf_area_from_biomass(bl,species,layer,firstlayer) result (area)
  real :: area ! returned value
  real,    intent(in) :: bl      ! biomass of leaves, kg C/individual
  integer, intent(in) :: species ! species
  integer, intent(in) :: layer, firstlayer

  if (layer > 1 .AND. firstlayer == 0) then
     area = bl/(0.5*spdata(species)%LMA) ! half thickness for leaves in understory
  else
     area = bl/spdata(species)%LMA    
  endif
end function 


! ============================================================================
function height_from_biomass(btotal) result (height)
    real :: height ! return value
    real, intent(in) :: btotal ! total biomass
    
    height = 24.19*(1.0-exp(-0.19*btotal))
end function


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
  c%Psw_alphasw = spdata(sp)%c3 * spdata(sp)%alpha(CMPT_LEAF) * D
  
end subroutine update_bio_living_fraction


! ============================================================================
! redistribute living biomass pools in a given cohort, and update related 
! properties (height, lai, sai)
subroutine update_biomass_pools(c)
  type(vegn_cohort_type), intent(inout) :: c

  if (do_ppa) return ! in PPA mode, do nothing at all

  c%height = height_from_biomass(btotal(c))
  call update_bio_living_fraction(c);
  c%bsw = c%Psw*c%bliving;
  if(c%status == LEAF_OFF) then
     c%blv = c%Pl*c%bliving + c%Pr*c%bliving;
     c%bl  = 0;
     c%br  = 0;
  else
     c%blv = 0;
     c%bl  = c%Pl*c%bliving;
     c%br  = c%Pr*c%bliving;
  endif
end subroutine 


! ============================================================================
! calculate tree height, DBH, height, and crown area by bwood and density 
! The allometry equations are from Ray Dybzinski et al. 2011 and Farrior et al. in review
!         HT = alphaHT * DBH ** (gamma-1)   ! DBH --> Height
!         CA = alphaCA * DBH ** gamma       ! DBH --> Crown Area
!         BM = alphaBM * DBH ** (gamma + 1) ! DBH --> tree biomass
subroutine init_cohort_allometry_ppa(cc)
  type(vegn_cohort_type), intent(inout) :: cc

  real    :: btot ! total biomass per individual, kg C

  btot = max(0.0001,cc%bwood+cc%bsw)
  associate(sp=>spdata(cc%species))
     cc%DBH        = (btot / sp%alphaBM) ** ( 1.0/sp%thetaBM )
!    cc%treeBM     = sp%alphaBM * cc%dbh ** sp%thetaBM
     cc%height     = sp%alphaHT * cc%dbh ** sp%thetaHT
     cc%crownarea  = sp%alphaCA * cc%dbh ** sp%thetaCA

     ! calculations of bl_max and br_max are here only for the sake of the
     ! diagnostics, because otherwise those fields are inherited from the 
     ! parent cohort and produce spike in the output, even though these spurious
     ! values are not used by the model
     cc%bl_max = sp%LMA   * sp%LAImax        * cc%crownarea
     cc%br_max = sp%phiRL * sp%LAImax/sp%SRA * cc%crownarea 
  end associate
end subroutine 

end module vegn_cohort_mod
