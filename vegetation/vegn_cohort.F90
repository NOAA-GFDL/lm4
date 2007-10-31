module vegn_cohort_mod

use land_constants_mod, only: NBANDS
use vegn_data_mod, only : spdata, &
   use_mcm_masking, use_bucket, critical_root_density, &
   tg_c4_thresh, tg_c3_thresh, l_fract, fsc_liv
use vegn_data_mod, only : PT_C3, PT_C4, CMPT_ROOT, CMPT_LEAF, &
   SP_C4GRASS, SP_C3GRASS, SP_TEMPDEC, SP_TROPICAL, SP_EVERGR, &
   LEAF_OFF, LU_CROP, PHEN_EVERGREEN, PHEN_DECIDIOUS
   
implicit none
private
! ==== public interfaces =====================================================
public :: vegn_phys_prog_type
public :: vegn_cohort_type

! operations defined for cohorts
!public :: new_cohort, delete_cohort
public :: vegn_data_heat_capacity
public :: vegn_data_intrcptn_cap
public :: get_vegn_wet_frac
public :: vegn_data_cover
public :: cohort_data_uptake_frac_max


public :: btotal ! returns cohort total biomass
public :: c3c4   ! returns physiology type for given biomasses and conditions
public :: phenology_type ! returns type of phenology for given conditions
public :: update_species ! updates cohort physiology, phenology type, and species
public :: height_from_biomass
public :: lai_from_biomass
public :: update_bio_living_fraction
public :: update_biomass_pools
! ==== end of public interfaces ==============================================

! ==== types =================================================================
type :: vegn_phys_prog_type
  real Wl
  real Ws
  real Tv
end type vegn_phys_prog_type

! vegn_cohort_type describes the data that belong to a vegetation cohort
type :: vegn_cohort_type
  type(vegn_phys_prog_type) :: prog

! ---- biological prognostic variables
! Currently bio prognostic variable is defined as anything that's saved in
! the restart; clearly some vars there are, strictly speaking, diagnostic,
! but saved for reproducibility (to avoid recalculation). exceptions :
! npp_previous_day is left outside, since it's obviously auxiliary; height
! is left outside
  integer :: species ! vegetation species
  real    :: bl      ! biomass of leaves, kg C/m2
  real    :: blv     ! biomass of virtual leaves (labile store), kg C/m2
  real    :: br      ! biomass of fine roots, kg C/m2
  real    :: bsw     ! biomass of sapwood, kg C/m2
  real    :: bwood   ! biomass of heartwood, kg C/m2

  real    :: bliving ! leaves, fine roots, and sapwood biomass
  integer :: status  ! growth status of plant

! ---- physical parameters
  real    :: height     ! vegetation height, m
  real    :: lai        ! leaf area index, m2/m2
  real    :: sai        ! stem area index, m2/m2
  real    :: leaf_size  ! leaf dimension, m
  real root_density
  real root_zeta
  real rs_min
  real    :: leaf_refl(NBANDS) ! reflectance of leaf, per band
  real    :: leaf_tran(NBANDS) ! transmittance of leaf, per band
  real    :: leaf_emis         ! emissivity of leaf
  real snow_crit  ! later parameterize this as snow_mask_fac*height

! ---- auxiliary variables 

  real    :: W_max   ! maximum water comtent of canopy, kg/(m2 of ground)
  real    :: mcv_dry ! heat capacity of dry canopy
  real    :: cover

  integer :: pt      ! physiology type
  integer :: phent

  real :: b       ! total biomass
  real :: babove  ! total above ground biomass
  real :: bs      ! structural biomass: stem + structural roots
  real :: bstem   ! stem biomass

  real :: gpp  = 0.0 ! gross primary productivity kg C/timestep
  real :: npp  = 0.0 ! net primary productivity kg C/timestep
  real :: npp2 = 0.0 ! temporarily stores eddy_npp
  real :: miami_npp ! stores miami-model npp

  real :: resp = 0.0 ! plant respiration
  real :: resl = 0.0 ! leaf respiration
  real :: resr = 0.0 ! root respiration
  real :: resg    ! growth respiration
  real :: md      ! plant tissue maintenance kg C/timestep

  real :: An_op ! mol C/(m2 of leaf per year)
  real :: An_cl ! mol C/(m2 of leaf per year)
  
  real :: carbon_gain = 0.0 ! carbon gain during the month
  real :: carbon_loss = 0.0 ! carbon loss during the month
  real :: bwood_gain  = 0.0 !

  ! used in fast time scale calculations
  real :: npp_previous_day
  real :: npp_previous_day_tmp = 0.0

  ! lena added this for storing previous size stomatal opening and lwnet 
  ! for computing canopy air T and q at the next step
  
  real :: gs;
  real :: gb;
  real :: frac_moist;
!moved to prog%Wl  real :: cmc;
!moved to prog%Tv  real :: tleaf ! temperature of leaves, degK

  real :: ds;

  ! new allocation fractions, Jan2 03
  real :: Pl          ! fraction of living biomass in leaves
  real :: Pr          ! fraction of living biomass in fine roots
  real :: Psw         ! fraction of living biomass in sapwood
  real :: Psw_alphasw ! fraction of sapwood times 
                      ! retirement rate of sapwood into wood
  real :: extinct ! light extinction coefficient in the canopy for photosynthesis calculations
  
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
subroutine vegn_data_heat_capacity ( cohort, mcv )
  type(vegn_cohort_type), intent(in)  :: cohort
  real                  , intent(out) :: mcv
  
  mcv = cohort%mcv_dry        ! later add term in vegn_pars%lai
end subroutine


! ============================================================================
subroutine vegn_data_intrcptn_cap ( cohort, vegn_W_max )
  type(vegn_cohort_type), intent(in)  :: cohort
  real                  , intent(out) :: vegn_W_max

  vegn_W_max = cohort%W_max  
end subroutine

! ============================================================================
! calculates functional dependence of wet canopy function f = x**p and its 
! derivative, but approximates it with linear function in the eps-vicinity 
! of zero, to make sure that derivative doesn't become infinite
subroutine wet_frac(w, w_max, p, eps, f, DfDw)
  real, intent(in) :: &
       w, &     ! water content
       w_max, & ! maximum storage capacity
       p, &     ! exponent of the function
       eps      ! neighbourhood of zero where we approximate x**p with linear function
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
  real, intent(out) :: &
       fw, DfwDwl, DfwDws, & ! water-covered fraction of canopy and its derivatives
       fs, DfsDwl, DfsDws    ! snow-covered fraction of canopy and its derivatives

  ! ---- local constants
  real, parameter :: eps = 0.01 ! value of w/w_max for transition to linear function

  ! ---- local vars
  integer :: sp  ! current cohort species
  real    :: fw0 ! total water-covered fraction (without overlap)

  if (cohort%W_max <= 0) then
     fw = 0.0; DfwDwl=0.0; DfwDws=0.0
     fs = 0.0; DfsDwl=0.0; DfsDws=0.0
     return
  endif

  sp = cohort%species

  ! snow-covered fraction
  call wet_frac(cohort%prog%Ws, cohort%W_max, spdata(sp)%cmc_pow, eps, fs,  DfsDws)
  DfsDwl = 0

  ! wet fraction
  call wet_frac(cohort%prog%Wl, cohort%W_max, spdata(sp)%cmc_pow, eps, fw0, DfwDwl)
  ! take into account overlap by snow
  fw     = fw0*(1-fs)
  DfwDwl = DfwDwl*(1-fs)
  DfwDws = -fw0*DfsDws

!!$  fs = max(cohort%prog%Ws, 0.0)/cohort%W_max; DfwDwl = 1.0/cohort%W_max
!!$  fw = max(cohort%prog%Wl, 0.0)/cohort%W_max; DfsDws = 1.0/cohort%W_max
!!$  ! shouldn't be there some check to make sure a_v<=1 and a_vs<=1?
end subroutine


! ============================================================================
subroutine vegn_data_cover ( cohort, snow_depth, vegn_cover, &
                                         vegn_cover_snow_factor )
  type(vegn_cohort_type), intent(inout)  :: cohort
  real, intent(in)  :: snow_depth
  real, intent(out) :: vegn_cover
  real, intent(out) :: vegn_cover_snow_factor

  cohort%cover = 1 - exp(-cohort%lai)
  if (use_mcm_masking) then
     vegn_cover_snow_factor =  &
           (1 - min(1., 0.5*sqrt(max(snow_depth,0.)/cohort%snow_crit)))
     cohort%cover = cohort%cover * &
           (1 - min(1., 0.5*sqrt(max(snow_depth,0.)/cohort%snow_crit)))
  else
     vegn_cover_snow_factor =  &
           cohort%snow_crit / &
          (max(snow_depth,0.0) + cohort%snow_crit)
     cohort%cover = cohort%cover * &
           cohort%snow_crit / &
          (max(snow_depth,0.0) + cohort%snow_crit)
  endif
  vegn_cover = cohort%cover
end subroutine vegn_data_cover


! ============================================================================
! calculates vertical distribution of active roots: given layer thicknesses,
! returns fraction of active roots per level
subroutine cohort_data_uptake_frac_max(cohort,dz,uptake_frac_max)
  type(vegn_cohort_type), intent(in)  :: cohort
  real, intent(in)  :: dz(:)
  real, intent(out) :: uptake_frac_max(:)

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
     z = 0
     do l = 1, size(dz)
        uptake_frac_max(l) = &
                cohort%root_density*(exp(-z/cohort%root_zeta)    &
                - exp(-(z+dz(l))/cohort%root_zeta))/dz(l)
        uptake_frac_max(l) = &
                max( uptake_frac_max(l)-critical_root_density, 0.0)*dz(l)
        z = z + dz(l)
     enddo
  endif
  
  sum_rf = sum(uptake_frac_max)
  if(sum_rf>0) &
       uptake_frac_max(:) = uptake_frac_max(:)/sum_rf
  
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
  real,              intent(in) :: temp   ! temperatire, degK
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
  
  if(pe>0.5 .and. pe<0.9) then
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
  
  c%species = spp
end subroutine


! ============================================================================
function height_from_biomass(btotal) result (height)
  real :: height ! returned value
  real, intent(in) :: btotal

  ! GCH, Function from JPC 2/9/02
!  height = 24.19*(1.0-exp(-0.19*(c%bliving+c%bwood)))
  height = 24.19*(1.0-exp(-0.19*btotal))
end function


! ============================================================================
function lai_from_biomass(bl,species) result (lai)
  real :: lai ! returned value
  real,    intent(in) :: bl      ! biomass of leaves, kg C/m2
  integer, intent(in) :: species ! species

  lai = bl*(spdata(species)%specific_leaf_area);   
end function


! ============================================================================
! calculates fractions of living biomass in differerent compartments
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

  c%b      = c%bliving + c%bwood;
  c%height = height_from_biomass(c%b);
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
  c%lai = lai_from_biomass(c%bl,c%species)
  c%sai = 0.035*c%height ! Federer and Lash,1978
end subroutine 


end module vegn_cohort_mod
