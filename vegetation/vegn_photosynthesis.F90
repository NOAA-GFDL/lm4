module vegn_photosynthesis_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif
use fms_mod, only: error_mesg, FATAL, WARNING, file_exist, close_file, check_nml_error, stdlog, &
      mpp_pe, mpp_root_pe, lowercase
use constants_mod,      only : TFREEZE, PI, rdgas, dens_h2o, grav
use sphum_mod,          only : qscomp

use land_constants_mod, only : Rugas, seconds_per_year, mol_h2o, mol_air, d608
use land_numerics_mod,  only : gammaU, gamma
use land_debug_mod,     only : is_watch_point, check_var_range
use land_data_mod,      only : log_version
use soil_tile_mod,      only : soil_tile_type, psi_wilt
use vegn_tile_mod,      only : vegn_tile_type
use vegn_data_mod,      only : PT_C4, FORM_GRASS, spdata, T_transp_min, &
                               ALLOM_EW, ALLOM_EW1, ALLOM_HML
use vegn_cohort_mod,    only : vegn_cohort_type, get_vegn_wet_frac
use uptake_mod,         only : darcy2d_uptake, darcy2d_uptake_solver
implicit none
private

! ==== public interfaces =====================================================
public :: vegn_photosynthesis_init
public :: vegn_photosynthesis
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_photosynthesis_mod'
#include "../shared/version_variable.inc"

! values for selector of CO2 option used for photosynthesis
integer, public, parameter :: &
    VEGN_PHOT_CO2_PRESCRIBED  = 1, &
    VEGN_PHOT_CO2_INTERACTIVE = 2

! values for internal vegetation photosynthesis option selector
integer, parameter :: &
    VEGN_PHOT_SIMPLE  = 1, &   ! zero photosynthesis
    VEGN_PHOT_LEUNING = 2      ! photosynthesis according to simplified Leuning model

! values for internal vegetation water stress option selector
integer, public, parameter :: & ! water limitation options
    WSTRESS_NONE       = 1, &  ! no water limitation
    WSTRESS_LM3        = 2, &  ! LM3-like
    WSTRESS_HYDRAULICS = 3     ! plant hydraulics formulation

real, parameter :: gs_lim = 0.25 ! maximum stomatal cond for Leuning, mol/(m2 s)
real, parameter :: b      = 0.01 ! minimum stomatal cond for Leuning, mol/(m2 s)

! ==== module variables ======================================================
integer :: vegn_phot_option     = -1 ! selector of the photosynthesis option
integer, public, protected :: vegn_phot_co2_option = -1 ! internal selector of co2 option for photosynthesis
integer :: water_stress_option  = -1 ! selector of the water stress option

character(32) :: photosynthesis_to_use = 'simple' ! or 'leuning'

logical       :: Kok_effect  = .FALSE. ! if TRUE, Kok effect is taken in photosynthesis
!real          :: light_kok   = 0.00004 !mol_of_quanta/(m^2s) PAR
!real          :: Inib_factor = 0.5

character(32) :: co2_to_use_for_photosynthesis = 'prescribed' ! or 'interactive'
   ! specifies what co2 concentration to use for photosynthesis calculations:
   ! 'prescribed'  : a prescribed value is used, equal to co2_for_photosynthesis
   !      specified below.
   ! 'interactive' : concentration of co2 in canopy air is used
real, public, protected :: co2_for_photosynthesis = 350.0e-6 ! concentration of co2 for
   ! photosynthesis calculations, mol/mol. Ignored if co2_to_use_for_photosynthesis is
   ! not 'prescribed'

real :: lai_eps = 1e-5 ! threshold for switching to linear approximation for Ag_l, m2/m2

character(32) :: water_stress_to_use = 'lm3' ! type of water stress formulation:
   ! 'lm3', 'plant-hydraulics', or 'none'
logical :: hydraulics_repair = .TRUE.

namelist /photosynthesis_nml/ &
    photosynthesis_to_use, &
    Kok_effect,  &
    co2_to_use_for_photosynthesis, co2_for_photosynthesis, &
    lai_eps, &
    water_stress_to_use, hydraulics_repair

contains


! ============================================================================
subroutine vegn_photosynthesis_init()

  integer :: unit, io, ierr
  
  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=photosynthesis_nml, iostat=io)
    ierr = check_nml_error(io, 'photosynthesis_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=photosynthesis_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'photosynthesis_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  unit=stdlog()
  if (mpp_pe() == mpp_root_pe()) then
     write(unit, nml=photosynthesis_nml)
  endif

  ! convert symbolic names of photosynthesis options into numeric IDs to
  ! speed up selection during run-time
  if (trim(lowercase(photosynthesis_to_use))=='simple') then
     vegn_phot_option = VEGN_PHOT_SIMPLE
  else if (trim(lowercase(photosynthesis_to_use))=='leuning') then
     vegn_phot_option = VEGN_PHOT_LEUNING
  else
     call error_mesg('vegn_photosynthesis_init',&
          'vegetation photosynthesis option photosynthesis_to_use="'//&
          trim(photosynthesis_to_use)//'" is invalid, use "simple" or "leuning"',&
          FATAL)
  endif

  ! convert symbolic names of photosynthesis CO2 options into numeric IDs to
  ! speed up selection during run-time
  if (trim(lowercase(co2_to_use_for_photosynthesis))=='prescribed') then
     vegn_phot_co2_option = VEGN_PHOT_CO2_PRESCRIBED
  else if (trim(lowercase(co2_to_use_for_photosynthesis))=='interactive') then
     vegn_phot_co2_option = VEGN_PHOT_CO2_INTERACTIVE
  else
     call error_mesg('vegn_photosynthesis_init',&
          'vegetation photosynthesis option co2_to_use_for_photosynthesis="'//&
          trim(co2_to_use_for_photosynthesis)//'" is invalid, use "prescribed" or "interactive"',&
          FATAL)
  endif

  ! parse the options for the water stress
  if (trim(lowercase(water_stress_to_use))=='none') then
     water_stress_option = WSTRESS_NONE
  else if (trim(lowercase(water_stress_to_use))=='lm3') then
     water_stress_option = WSTRESS_LM3
  else if (trim(lowercase(water_stress_to_use))=='plant-hydraulics') then
     water_stress_option = WSTRESS_HYDRAULICS
  else
     call error_mesg('vegn_photosynthesis_init',&
        'option water_stress_to_use="'//trim(water_stress_to_use)//'" is invalid, use '// &
        ' "lm3", "plant-hydraulics", or "none"', &
        FATAL)
  endif
end subroutine vegn_photosynthesis_init


! ============================================================================
! compute stomatal conductance, photosynthesis and respiration
! updates cohort%An_op and cohort%An_cl
subroutine vegn_photosynthesis ( soil, vegn, cohort, &
     PAR_dn, PAR_net, cana_T, cana_q, cana_co2, p_surf, drag_q, &
     soil_beta, soil_water_supply, con_v_v, &
! output:
     evap_demand, stomatal_cond, RHi, &
     lai_kok, Anlayer,lai_light)
  type(soil_tile_type),   intent(in)    :: soil
  type(vegn_tile_type),   intent(in)    :: vegn
  type(vegn_cohort_type), intent(inout) :: cohort
  real, intent(in)  :: PAR_dn   ! downward PAR at the top of the cohort, W/m2
  real, intent(in)  :: PAR_net  ! net PAR absorbed by the cohort, W/m2
  real, intent(in)  :: cana_T   ! temperature of canopy air, deg K
  real, intent(in)  :: cana_q   ! specific humidity in canopy air space, kg/kg
  real, intent(in)  :: cana_co2 ! co2 concentration in canopy air space, mol CO2/mol dry air
  real, intent(in)  :: p_surf   ! surface pressure
  real, intent(in)  :: drag_q   ! drag coefficient for specific humidity
  real, intent(in)  :: soil_beta
  real, intent(in)  :: soil_water_supply ! max supply of water to roots, kg/(indiv s)
  real, intent(in)  :: con_v_v  ! one-sided foliage-CAS conductance per unit ground area
  real, intent(out) :: evap_demand   ! transpiration water demand, kg/(indiv s)
  real, intent(out) :: stomatal_cond ! stomatal conductance, m/s
  real, intent(out) :: RHi      ! relative humidity inside leaf, at the point of vaporization
  real, intent(out) :: lai_kok  ! LAI value for light inhibition m2/m2
  real, intent(out) :: Anlayer
  real, intent(out) :: lai_light ! LAI at which Ag=Resp
  
  select case (vegn_phot_option)
  case(VEGN_PHOT_SIMPLE)
     ! beta non-unity only for "beta" models
     stomatal_cond = soil_beta / (cohort%rs_min  + (1-soil_beta)/drag_q)
     cohort%An_op  = 0
     cohort%An_cl  = 0
     RHi = 1.0
     evap_demand   = 0
  case(VEGN_PHOT_LEUNING)
     call vegn_photosynthesis_Leuning (soil, vegn, cohort, &
            PAR_dn, PAR_net, cana_T, cana_q, cana_co2, p_surf, &
            soil_water_supply, con_v_v, evap_demand, stomatal_cond, RHi, &
            lai_kok, Anlayer, lai_light )
  case default
     call error_mesg('vegn_stomatal_cond', &
          'invalid vegetation photosynthesis option', FATAL)
  end select

end subroutine vegn_photosynthesis


subroutine vegn_photosynthesis_Leuning (soil, vegn, cohort, &
     PAR_dn, PAR_net, cana_T, cana_q, cana_co2, p_surf, &
     soil_water_supply, con_v_v, &
     evap_demand, stomatal_cond, RHi, &
     lai_kok, Anlayer, lai_light)
  type(soil_tile_type),   intent(in)    :: soil
  type(vegn_tile_type),   intent(in)    :: vegn
  type(vegn_cohort_type), intent(inout) :: cohort
  real, intent(in)  :: PAR_dn   ! downward PAR at the top of the canopy, W/m2
  real, intent(in)  :: PAR_net  ! net PAR absorbed by the canopy, W/m2
  real, intent(in)  :: cana_T   ! temperature of canopy air, deg K
  real, intent(in)  :: cana_q   ! specific humidity in canopy air space, kg/kg
  real, intent(in)  :: cana_co2 ! co2 concentration in canopy air space, mol CO2/mol dry air
  real, intent(in)  :: p_surf   ! surface pressure
  real, intent(in)  :: soil_water_supply ! max supply of water to roots per unit
                                ! active root biomass per second, kg/(indiv s)
  real, intent(in)  :: con_v_v  ! one-sided foliage-CAS conductance per unit ground area
  real, intent(out) :: evap_demand ! transpiration water demand, kg/(indiv s)
  real, intent(out) :: stomatal_cond ! stomatal conductance, m/s
  real, intent(out) :: RHi      ! relative humidity inside leaf, at the point of vaporization
  real, intent(out) :: lai_kok  ! LAI value for light inhibition m2/m2
  real, intent(out) :: Anlayer
  real, intent(out) :: lai_light ! LAI at which Ag=Resp
  
  ! ---- local vars
  integer :: sp      ! shortcut for cohort%species
  real    :: fw, fs  ! wet and snow-covered fraction of leaves, dimensionless
  real    :: psyn    ! net photosynthesis, mol C/(m2 of leaves s)
  real    :: resp    ! leaf respiration, mol C/(m2 of leaves s)

  real    :: leaf_q  ! saturated specific humidity at leaf temperature, kg/kg
  real    :: ds      ! water vapor deficit, kg/kg
  real    :: w_scale ! water stress scaling factor, dimensionless

  real    :: rho     ! density of canopy air, kg/m3
  real    :: gs      ! max stomatal conductance per individual, m/s
  real    :: gb      ! aerodynamic conductance of canopy air per individual, m/s
  real    :: fdry    ! fraction of canopy not covered by intercepted water/snow
  real    :: evap_demand_old

  ! set the default values for outgoing parameters, overriden by the calculations
  ! in gs_leuning
  lai_kok   = 0.0
  Anlayer   = 0.0
  lai_light = 0.0
  
  if(cohort%lai <= 0) then
     ! no leaves means no photosynthesis and zero stomatal conductance, of course
     cohort%An_op  = 0
     cohort%An_cl  = 0
     stomatal_cond = 0
     evap_demand   = 0
     RHi           = 1.0
     ! TODO: call vegn_hydraulics?
     return
  endif

  sp = cohort%species ! shorthand, for convenience

  ! calculate humidity deficit, kg/kg
  call qscomp(cohort%Tv, p_surf, leaf_q)
  ds = max(leaf_q-cana_q,0.0)

  ! calculate non-water-limited stomatal conductance, photosynthesis, and respiration
  call gs_Leuning(PAR_dn, PAR_net, cohort%Tv, ds, cohort%lai, &
       cohort%leaf_age, p_surf, sp, cana_co2, cohort%extinct, &
       cohort%layer, &
       ! output:
       stomatal_cond, psyn, resp, &
       lai_kok,Anlayer, lai_light)

  ! scale down stomatal conductance and photosynthesis due to leaf wetness
  ! effect
  call get_vegn_wet_frac (cohort, fw=fw, fs=fs)
  stomatal_cond = stomatal_cond*(1-spdata(sp)%wet_leaf_dreg*fw);
  if (psyn > 0) &
     psyn = psyn*(1-spdata(sp)%wet_leaf_dreg*fw);
  ! limit stomatal conductance
  if (stomatal_cond > gs_lim) then
     if(psyn > 0) psyn = psyn*gs_lim/stomatal_cond
     stomatal_cond = gs_lim
  endif

  ! calculate transpiration water demand
  RHi = 1.0 ! max relative humidity of in-canopy air
  ! convert units of stomatal conductance from mol/(m2 s) to m/(s indiv) by
  ! multiplying it by a volume of a mole of gas and by leaf area of individual
  gs = stomatal_cond * Rugas * cohort%Tv / p_surf * cohort%leafarea
  ! mol_h2o/molair = 18.0e-3( molar mass of water, kg/mol_air) /wtmair/1000.0 ( molar mass of air, kg). where wtmair=2.896440E+01

  ! aerodynamic conductance per individual
  gb = con_v_v * cohort%crownarea ! foliage-CAS conductance per m2 => per individual

  if (psyn > 0) psyn = psyn*gb/(gs+gb) ! decrease photosynthesis due to CO2 supply towards the leaf by 1/gb

  rho = p_surf/(rdgas*cana_T*(1+d608*cana_q)) ! canopy air density
  !fdry = 1-fw-fs ! fraction of canopy that is dry
  fdry = 1
  ! transpiration demand, kg/(indiv s)
  evap_demand = rho * fdry * gs*gb/(gs+gb) * (leaf_q*RHi-cana_q)
  ! in LM3 old demand = gs_w*ds*mol_air/mol_h2o;
  ! ! calculate humidity deficit, kg/kg
  !  call qscomp(tl, p_surf, hl)
  ! ds = max(hl-ea,0.0)

  evap_demand_old = stomatal_cond * ds * 18.0e-3/2.896440E+01/1000.0 * cohort%leafarea

   ! ens test evap_demand = rho  * gs * (leaf_q*RHi-cana_q)
  ! scale down stomatal conductance and photosynthesis due to water stress
  select case (water_stress_option)
  case (WSTRESS_LM3)
     if (evap_demand>soil_water_supply) then
        w_scale = soil_water_supply/evap_demand
     else
        w_scale = 1.0
     endif


  case (WSTRESS_HYDRAULICS)
     call vegn_hydraulics(soil, vegn, cohort, p_surf, cana_T, cana_q, gb, gs, fdry, &
          w_scale, RHi )
  case (WSTRESS_NONE)
     w_scale = 1.0
  case default
     call error_mesg('vegn_stomatal_cond', 'invalid vegetation water stress option', FATAL)
  end select
!print *, 'w_scale ', w_scale, 'evap_demand', evap_demand, ' old ',evap_demand_old, 'water_supply', &
!soil_water_supply, ' gb red ', gb/(gs+gb), ' fdry ', fdry, ' rho ', rho

  stomatal_cond=w_scale*stomatal_cond
  if(psyn > 0) psyn = psyn*w_scale
  if(psyn < 0.and.stomatal_cond>b) stomatal_cond=b
  if (is_watch_point()) then
     __DEBUG4__(w_scale, stomatal_cond, psyn, b)
     __DEBUG2__(gs,stomatal_cond)
  endif
  ! store the calculated photosynthesis and photorespiration for future use
  ! in carbon_int
  cohort%An_op  = psyn * seconds_per_year
  cohort%An_cl  = resp * seconds_per_year
  ! convert units of stomatal conductance from mol/(m2 s) to m/s by
  ! multiplying it by a volume of a mole of gas
  stomatal_cond   = stomatal_cond * Rugas * cohort%Tv / p_surf
  ! convert stomatal conductance from units per unit area of leaf to the
  ! units per unit area of cohort
  stomatal_cond = stomatal_cond*cohort%lai
  ! store w_scale for diagnostics
  cohort%w_scale = w_scale

end subroutine vegn_photosynthesis_Leuning


! ============================================================================
! calculates non-water-stressed photosynthesis and stomatal conductance,
! vertically averaged over the cohort canopy
subroutine gs_Leuning(rad_top, rad_net, tl, ds, lai, leaf_age, &
                   p_surf, pft, ca, kappa, layer, &
                   gs, apot, acl, &
                   lai_kok, Anlayer, lai_light)
  real,    intent(in)    :: rad_top ! PAR dn on top of the canopy, w/m2
  real,    intent(in)    :: rad_net ! PAR net on top of the canopy, w/m2
  real,    intent(in)    :: tl   ! leaf temperature, degK
  real,    intent(in)    :: ds   ! humidity deficit, kg/kg
  real,    intent(in)    :: lai  ! leaf area index
  real,    intent(in)    :: leaf_age ! age of leaf since bud burst (deciduous), days
  real,    intent(in)    :: p_surf ! surface pressure, Pa
  integer, intent(in)    :: pft  ! species
  real,    intent(in)    :: ca   ! concentration of CO2 in the canopy air space, mol CO2/mol dry air
  real,    intent(in)    :: kappa! canopy extinction coefficient (move inside f(pft))
  integer, intent(in)    :: layer ! canopy layer of the cohort
  ! note that the output is per area of leaf; to get the quantities per area of
  ! land, multiply them by LAI
  real,    intent(out)   :: gs   ! stomatal conductance, mol/(m2 s)
  real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
  real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)
  !#### Modified by PPG 2017-11-03
  real,    intent(out)   :: lai_kok ! Lai at which kok effect is considered
  real,    intent(out)   :: Anlayer
  real,    intent(out)   :: lai_light ! Lai at which Ag=Resp
  ! ---- local vars
  ! photosynthesis
  real :: vm;
  real :: kc,ko; ! Michaelis-Menten constants for CO2 and O2, respectively
  real :: ci;
  real :: capgam; ! CO2 compensation point
  real :: f2,f3;
  real :: coef0,coef1;

  real :: Resp;

  ! conductance related
  real :: do1;

  ! miscellaneous
  real :: dum2;
  real, parameter :: light_crit = 0

  ! new average computations
  real :: lai_eq;
  real, parameter :: rad_phot = 0.0000046 ! PAR conversion factor of J -> mol of quanta
  real :: light_top;
  real :: par_net;
  real :: Ag;
  real :: An;
  real :: Ag_l;
  real :: Ag_rb;
  real :: anbar;
  real :: gsbar;
  real :: w_scale;
  real, parameter :: p_sea = 1.0e5 ! sea level pressure, Pa

  if (is_watch_point()) then
     write(*,*) '####### gs_leuning input #######'
     __DEBUG2__(rad_top, rad_net)
     __DEBUG1__(tl)
     __DEBUG1__(ds)
     __DEBUG1__(lai)
     __DEBUG1__(leaf_age)
     __DEBUG1__(p_surf)
     __DEBUG1__(pft)
     __DEBUG1__(ca)
     __DEBUG1__(kappa)
     write(*,*) '####### end of ### gs_leuning input #######'
  endif

  do1=0.09 ; ! kg/kg
  
  if (spdata(pft)%lifeform == FORM_GRASS) do1=0.15;


  ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
  ! empirical relationship from McCree is light=rn*0.0000046
  light_top = rad_top*rad_phot;
  par_net   = rad_net*rad_phot;

  ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986

  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;

  vm=spdata(pft)%Vmax*exp(3000.0*(1.0/288.2-1.0/tl));

  if (layer > 1) vm=vm*spdata(pft)%Vmax_understory_factor ! reduce vmax in the understory

  !decrease Vmax due to aging of temperate deciduous leaves
  !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
  if (spdata(pft)%leaf_age_tau>0 .and. leaf_age>spdata(pft)%leaf_age_onset) then
     vm=vm*exp(-(leaf_age-spdata(pft)%leaf_age_onset)/spdata(pft)%leaf_age_tau)
  endif

  capgam=0.5*kc/ko*0.21*0.209; ! Farquhar & Caemmerer 1982

  ! Find respiration for the whole canopy layer

  !Resp=spdata(pft)%gamma_resp*vm*lai !/layer
  ! Find respiration for the whole canopy layer
  if (light_top>spdata(pft)%light_kok .and. spdata(pft)%inib_factor>0.0) then
     lai_kok=min(log(light_top/spdata(pft)%light_kok)/kappa,lai)
  else
     lai_kok = 0.0
  endif
  if (Kok_effect) then
     ! modify vm for Vmax later and add a temperature function to it.
     Resp=(1-spdata(pft)%inib_factor)*spdata(pft)%gamma_resp*vm*lai_kok+spdata(pft)%gamma_resp*vm*(lai-lai_kok)
  else
     Resp=spdata(pft)%gamma_resp*vm*lai
  endif

  if (layer > 1) Resp = Resp*spdata(pft)%Resp_understory_factor ! reduce gamma_resp by 50% per Steve's experiments - need a reference

  Resp=Resp/((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));


  ! ignore the difference in concentrations of CO2 near
  !  the leaf and in the canopy air, rb=0.

  Ag_l=0.;
  Ag_rb=0.;
  Ag=0.;
  anbar=-Resp/lai;
  gsbar=b;
  lai_light = 0.0

  ! find the LAI level at which gross photosynthesis rates are equal
  ! only if PAR is positive
  if ( light_top > light_crit ) then
     if (spdata(pft)%pt==PT_C4) then ! C4 species
        coef0=(1+ds/do1)/spdata(pft)%m_cond;
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
        if (ci>capgam) then
           f2=vm;
           f3=18000.0*vm*ci;

           dum2=min(f2,f3)

           ! find LAI level at which rubisco limited rate is equal to light limited rate
           lai_eq = -log(dum2/(kappa*spdata(pft)%alpha_phot*light_top))/kappa;
           lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

           ! gross photosynthesis for light-limited part of the canopy
           if (lai>lai_eps) then
              Ag_l   = spdata(pft)%alpha_phot * par_net &
                   * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
           else
              ! approximation for very small LAI: needed because the general formula above
              ! produces division by zero if LAI is very close to zero
              Ag_l = spdata(pft)%alpha_phot * par_net &
                   * (lai-lai_eq)/lai
           endif
           ! gross photosynthesis for rubisco-limited part of the canopy
           Ag_rb  = dum2*lai_eq

           Ag=(Ag_l+Ag_rb)/ &
             ((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
           An=Ag-Resp;
           anbar=An/lai;

           if(anbar>0.0) then
               gsbar=anbar/(ci-capgam)/coef0;
           endif
        endif ! ci>capgam
     else ! C3 species
        coef0=(1+ds/do1)/spdata(pft)%m_cond;
        coef1=kc*(1.0+0.209/ko);
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
        f2=vm*(ci-capgam)/(ci+coef1);
        f3=vm/2.;
        dum2=min(f2,f3);
        if (ci>capgam) then
           ! find LAI level at which rubisco limited rate is equal to light limited rate
           lai_eq=-log(dum2*(ci+2.*capgam)/(ci-capgam)/ &
                       (spdata(pft)%alpha_phot*light_top*kappa))/kappa;
           lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

           ! gross photosynthesis for light-limited part of the canopy
           if (lai>lai_eps) then
              Ag_l   = spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * par_net &
                   * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
           else
              ! approximation for very small LAI: needed because the general formula above
              ! produces division by zero if LAI is very close to zero
              Ag_l = spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * par_net &
                   * (lai-lai_eq)/lai
           endif

           ! gross photosynthesis for rubisco-limited part of the canopy
           Ag_rb  = dum2*lai_eq

           Ag=(Ag_l+Ag_rb)/((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
           An=Ag-Resp;
           anbar=An/lai;

           if(anbar>0.0) then
               gsbar=anbar/(ci-capgam)/coef0;
           endif
        endif ! ci>capgam
     endif
  endif ! light is available for photosynthesis

  apot = anbar
  gs   = gsbar
  acl  = -Resp/lai

  if (is_watch_point()) then
     __DEBUG4__(gs, apot, acl, ds)
  endif
  !write (*,*) 'apot  ',apot,  'ag_rb ', ag_rb,'vm ', vm, 'ci_f ', (ci-capgam)/(ci+coef1)
end subroutine gs_Leuning


! ==============================================================================
! calculate the stomatal conductance reduction due to water stress
subroutine vegn_hydraulics(soil, vegn, cc, p_surf, cana_T, cana_q, gb, gs0, fdry, &
      w_scale, RHi )
  type(soil_tile_type),   intent(in)    :: soil
  type(vegn_tile_type),   intent(in)    :: vegn
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in)  :: cana_T   ! canopy air temperature, degK
  real, intent(in)  :: cana_q   ! canopy air specific humidity, kg/kg
  real, intent(in)  :: p_surf   ! surface pressure, N/m2 = kg/m/s2 = Pa
  real, intent(in)  :: gb       ! aerodynamic conductance of canopy air per individual, m/s
  real, intent(in)  :: gs0      ! non-water-limited stomatal conductance per individual, m/s
  real, intent(in)  :: fdry     ! fraction of canopy not covered by intercepted water/snow
  real, intent(out) :: w_scale  ! reduction of stomatal conductance due to water stress
  real, intent(out) :: RHi      ! relative humidity inside leaf, at the point of vaporization

  ! ---- local constants ----
  real, parameter :: Vm = 18.0e-6 ! m3/mol partial molar volume of liquid water

  ! ---- local vars ----
  real :: rho ! density of canopy air, kg/m3
  real :: gs, DgsDpl, DgsDTl  ! stomatal conductance per individual, m/s, and its derivatives
  real :: Et0, DetDpl, DetDTl ! transpiration, kg/(s indiv), and its derivatives
  real :: ur0_(size(soil%wl)), DurDpr_(size(soil%wl)) ! water flux from roots, kg/(s indiv), and its derivatives
  real :: ur0, DurDpr ! water flux from roots, kg/(s indiv), and its derivatives
  real :: ux0, DuxDpr, DuxDpx ! water flux in stem, kg/(s indiv), and its derivatives
  real :: ul0, DulDpx, DulDpl ! water flux in leaves, kg/(s indiv), and its derivatives
  real :: qsat, DqsatDT ! saturated spec. humidity and its derivative at the leaf temperature
  real :: gamma_r, ar, br, gamma_x, ax, bx, gamma_l ! variables for back substitution
  ! TODO: define CSAsw in LM3 case
  real :: CSAsw ! cross-section area of sapwood, m2
  real :: psi_r ! psi of root (root-stem interface), Pa
  real :: delta_pr, delta_px, delta_pl ! tendencies of water potentials, Pa
  integer :: n

  real, parameter :: m2pa = dens_h2o*grav ! unit conversion factor, pa per m of water
  real, parameter :: small_Et0 = 1e-10 ! small transpiration value used to get to
  ! to the region of non-zero derivatives w.r.t. psi_r

  associate(sp=>spdata(cc%species), vegn_T=>cc%Tv)
  if (is_watch_point()) then
     write(*,*)'######### vegn_hydraulics input ###########'
     __DEBUG3__(cana_T, cana_q, p_surf)
     __DEBUG3__(gb,gs0,fdry)
     __DEBUG3__(cc%psi_l, cc%psi_x, cc%psi_r)
     __DEBUG2__(sp%cl, sp%dl)
     __DEBUG2__(sp%cx, sp%dx)
     __DEBUG2__(cc%Kla, sp%Kxam)
     __DEBUG2__(cc%leafarea,cc%height)

     write(*,*)'######### end of vegn_hydraulics input ###########'
  endif

  ! calculate explicit estimate of the transpiration (per individual),
  ! and its derivatives w.r.t. leaf water potential and leaf temperature
  ! TODO: what is gong on with fw and fl in vegn_step_1? Are we double-counting?
  !       are we consistent with the photosynthesis?
  call qscomp(vegn_T, p_surf, qsat, DqsatDT)
  rho = p_surf/(rdgas*cana_T*(1+d608*cana_q)) ! canopy air density
  RHi = exp(cc%psi_l*Vm/(Rugas*vegn_T)) ! relative humidity inside the leaf
  gs = gs0 * exp(-(cc%psi_l/sp%dl)**sp%cl)
  DgsDpl = - gs0 * exp(-(cc%psi_l/sp%dl)**sp%cl)*sp%cl/sp%dl*(cc%psi_l/sp%dl)**(sp%cl-1)
  DgsDTl = 0.0

  Et0 = rho * fdry * gs*gb/(gs+gb) * (qsat*RHi-cana_q)

  if (is_watch_point()) then
     __DEBUG2__(RHi,cc%Tv)
     __DEBUG2__(gb,gs0)
     __DEBUG2__(gs,DgsDpl)
     __DEBUG2__(qsat,cana_q)
     __DEBUG1__(Et0)
!     __DEBUG1__(vegn%root_distance)
!     __DEBUG1__(soil%psi)
!     __DEBUG1__(soil%wl)
!     __DEBUG1__(soil%ws)
!     __DEBUG3__(cc%root_length, cc%K_r, cc%r_r)
  endif

  if (Et0<=0) then
     ! negative transpiration is not allowed
     Et0    = 0.0 ; DEtDpl = 0.0 ; DetDpl = 0.0
     ! set all water potentials to water potential of roots (which is the max
     ! of all water potentials)
     cc%psi_x = cc%psi_r
     cc%psi_l = cc%psi_r
  else
     DEtDpl = rho * fdry * gb/(gs+gb) * ( gs*qsat*Vm/(Rugas*cc%Tv)*RHi +      &
                                  DgsDpl*gb/(gs+gb)*(qsat*Rhi-cana_q)  )
     !dTl derivative is not used
     DEtDTl = rho * fdry * gs*gb/(gs+gb)*RHi*(DqsatDT - qsat*cc%psi_l*Vm/(Rugas*vegn_T**2))

     ! TODO: generalize for linear water uptake

     ! calculate water uptake by roots and its derivative w.r.t. water potential at
     ! root-stem interface
     psi_r = cc%psi_r
     call darcy2d_uptake ( soil, psi_r/m2pa, vegn%root_distance, cc%root_length, &
         cc%K_r, cc%r_r, ur0_, DurDpr_ )
     ur0 = sum(ur0_); DurDpr = sum(DurDpr_)/m2pa ! converting derivative from head (m) to Pascal
     if (ur0<epsilon(1.0).and.abs(DurDpr)<epsilon(1.0)) then
        ! this can happen if initial psi_r is higher than soil water potential,
        ! and uptake_one_way is true. Or there is no water in the soil. Or soil
        ! is frozen

        ! In this case, we solve for Darcy-flow uptake, find the root water potential
        ! to satisfy transpiration by the vegetation, use resulting psi as initial
        ! condition
        if (is_watch_point()) then
           write (*,*)'###### ur0 and DurDpr are 0; adjusting psi_r #####'
        endif
        call darcy2d_uptake_solver     (soil, max(Et0,small_Et0), vegn%root_distance, &
                cc%root_length, cc%K_r, cc%r_r, &
                ur0_, psi_r, n)
        psi_r = psi_r*m2pa

        ! shift initial condition for psi_r a bit below minimum soil psi, so that
        ! we do have some suction
   !     psi_r = minval(soil%psi) - 1.0   ! in m of water
   !     psi_r = max(psi_r,psi_wilt)*m2pa ! but we don't want to go below psi_wilt

        ! recalculate root water flow and its derivative, now that we believe the
        ! derivative is non-zero
        call darcy2d_uptake ( soil, psi_r/m2pa, vegn%root_distance, cc%root_length, &
            cc%K_r, cc%r_r, ur0_, DurDpr_ )
        ur0 = sum(ur0_); DurDpr = sum(DurDpr_)/m2pa ! converting derivative from head (m) to Pascal
        if (is_watch_point()) then
           write (*,*)'###### finished adjusting psi_r #####'
           __DEBUG1__(psi_r)
        endif
     endif

     ! calculate stem flow and its derivatives
     ! isa and es 201701 - CSAsw different for ALLOM_HML
     select case(sp%allomt)
     case (ALLOM_EW,ALLOM_EW1)
        CSAsw  = sp%alphaCSASW * cc%DBH**sp%thetaCSASW ! cross-section area of sapwood
     case (ALLOM_HML)
        CSAsw = sp%phiCSA * cc%DBH**(sp%thetaCA + sp%thetaHT) / (sp%gammaHT + cc%DBH** sp%thetaHT)
     end select
     ! isa 20170705 - grasses don't form heartwood
     if (sp%lifeform == FORM_GRASS) then
        CSAsw = PI * cc%DBH * cc%DBH / 4.0 ! trunk cross-sectional area = sapwood area
     endif

     cc%Kxi    =  sp%Kxam / cc%height * CSAsw
     ux0    = -cc%Kxi*sp%dx/sp%cx*gamma(1/sp%cx)*( &
                     gammaU((   psi_r/sp%dx)**sp%cx, 1/sp%cx) - &
                     gammaU((cc%psi_x/sp%dx)**sp%cx, 1/sp%cx))
     DuxDpr =  cc%Kxi*exp(-(   psi_r/sp%dx)**sp%cx)
     DuxDpx = -cc%Kxi*exp(-(cc%psi_x/sp%dx)**sp%cx)

   ! calculate leaf flow and its derivatives
     cc%Kli = cc%Kla * cc%leafarea
     ul0    = -cc%Kli*sp%dl/sp%cl*gamma(1/sp%cl)*(gammaU((cc%psi_x/sp%dl)**sp%cl, 1/sp%cl) &
                                             - gammaU((cc%psi_l/sp%dl)**sp%cl, 1/sp%cl))
     DulDpx =  cc%Kli*exp(-(cc%psi_x/sp%dl)**sp%cl)
     DulDpl = -cc%Kli*exp(-(cc%psi_l/sp%dl)**sp%cl)

     ! do forward elimination
     gamma_r = 1/(DuxDpr - DurDpr)
     ar = (ur0-ux0)*gamma_r
     br = -gamma_r*DuxDpx

     gamma_x = 1/(DulDpx - DuxDpx - br*DuxDpr)
     ax = (ux0-ul0+ar*DuxDpr)*gamma_x
     bx = -gamma_x*DulDpl

     if(DetDpl - DulDpl - bx*DulDpx/=0) then
        gamma_l = 1/(DetDpl - DulDpl - bx*DulDpx)
     else
        ! this happened at least once due to precision loss
        gamma_l = 0.0
     endif
     ! calculate tendencies
     delta_pl = gamma_l*(ul0 - Et0 + ax*DulDpx)
     delta_px = ax+bx*delta_pl
     delta_pr = ar+br*delta_px

     ! update water potentials
     cc%psi_r = min(psi_r    + delta_pr,0.0)
     cc%psi_x = min(cc%psi_x + delta_px,cc%psi_r)
     cc%psi_l = min(cc%psi_l + delta_pl,cc%psi_x)
  endif

  call check_var_range(cc%psi_l,-HUGE(1.0),0.0,'vegn_hydraulics','psi_l',WARNING)
  call check_var_range(cc%psi_x,-HUGE(1.0),0.0,'vegn_hydraulics','psi_x',WARNING)
  call check_var_range(cc%psi_r,-HUGE(1.0),0.0,'vegn_hydraulics','psi_r',WARNING)
  w_scale = exp(-(cc%psi_l/sp%dl)**sp%cl)
  if (is_watch_point()) then
!     __DEBUG1__(ur0_)
!     __DEBUG1__(DurDpr_)
     __DEBUG3__(cc%Kli, cc%Kxi, CSAsw)
     __DEBUG3__(sp%alphaCSASW, cc%DBH, sp%thetaCSASW)
     __DEBUG2__(ur0, DurDpr)
     __DEBUG3__(ux0, DuxDpr, DuxDpx)
     __DEBUG3__(ul0, DulDpx, DulDpl)
     __DEBUG3__(Et0, DetDpl, DetDtl)
     __DEBUG3__(gamma_r, ar, br)
     __DEBUG3__(gamma_x, ax, bx)
     __DEBUG1__(gamma_l)
     __DEBUG3__(delta_pl, delta_px, delta_pr)
     __DEBUG3__(cc%psi_l, cc%psi_x, cc%psi_r)
     __DEBUG1__(w_scale)
  endif

  end associate

end subroutine vegn_hydraulics

end module vegn_photosynthesis_mod
