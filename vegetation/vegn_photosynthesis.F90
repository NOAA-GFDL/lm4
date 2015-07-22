module vegn_photosynthesis_mod

#include "../shared/debug.inc"

use fms_mod,            only : write_version_number, error_mesg, FATAL
use constants_mod,      only : TFREEZE, PI, rdgas, dens_h2o, grav 
use sphum_mod,          only : qscomp

use land_constants_mod, only : Rugas, seconds_per_year, mol_h2o, mol_air, d608
use land_numerics_mod,  only : m44inv, gammaU
use land_debug_mod,     only : is_watch_point
use soil_tile_mod,      only : soil_tile_type
use vegn_data_mod,      only : PT_C4, FORM_GRASS, spdata, T_transp_min
use vegn_cohort_mod,    only : vegn_cohort_type, get_vegn_wet_frac

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_photosynthesis_init
public :: vegn_photosynthesis
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id$', &
   tagname = '$Name$', &
   module_name = 'vegn_photosynthesis'
! values for internal vegetation photosynthesis option selector
integer, parameter :: VEGN_PHOT_SIMPLE  = 1 ! zero photosynthesis
integer, parameter :: VEGN_PHOT_LEUNING = 2 ! photosynthesis according to simplified Leuning model

! values for internal vegetation water stress option selector
integer, public, parameter :: & ! water limitation options
   WSTRESS_NONE       = 1, &  ! no water limitation
   WSTRESS_LM3        = 2, &  ! LM3-like
   WSTRESS_HYDRAULICS = 3     ! plant hydraulics formulation

! ==== module variables ======================================================
integer :: vegn_phot_option    = -1 ! selector of the photosynthesis option
integer :: water_stress_option = -1 ! selector of the water stress option
logical :: hydraulics_repair    = .TRUE.

contains


! ============================================================================
subroutine vegn_photosynthesis_init(photosynthesis_to_use, water_stress_to_use, hydraulics_repair_in)
  character(*), intent(in) :: photosynthesis_to_use
  character(*), intent(in) :: water_stress_to_use
  logical,      intent(in) :: hydraulics_repair_in

  call write_version_number(version, tagname)

  ! convert symbolic names of photosynthesis options into numeric IDs to
  ! speed up selection during run-time
  if (trim(photosynthesis_to_use)=='simple') then
     vegn_phot_option = VEGN_PHOT_SIMPLE
  else if (trim(photosynthesis_to_use)=='leuning') then
     vegn_phot_option = VEGN_PHOT_LEUNING
  else
     call error_mesg('vegn_photosynthesis_init',&
          'vegetation photosynthesis option photosynthesis_to_use="'//&
          trim(photosynthesis_to_use)//'" is invalid, use "simple" or "leuning"',&
          FATAL)
  endif

  ! parse the options for the water stress
  if (trim(water_stress_to_use)=='none') then
     water_stress_option = WSTRESS_NONE
  else if (trim(water_stress_to_use)=='lm3') then
     water_stress_option = WSTRESS_LM3
  else if (trim(water_stress_to_use)=='plant-hydraulics') then
     water_stress_option = WSTRESS_HYDRAULICS
  else
     call error_mesg('vegn_photosynthesis_init',&
        'option water_stress_to_use="'//trim(water_stress_to_use)//'" is invalid, use '// &
        ' "lm3", "plant-hydraulics", or "none"', &
        FATAL)
  endif

  hydraulics_repair = hydraulics_repair_in
end subroutine vegn_photosynthesis_init


! ============================================================================
! compute stomatal conductance, photosynthesis and respiration
! updates cohort%An_op and cohort%An_cl
subroutine vegn_photosynthesis ( cohort, &
     PAR_dn, PAR_net, cana_T, cana_q, cana_co2, p_surf, drag_q, &
     soil_beta, soil_water_supply, con_v_v, &
     stomatal_cond )
  ! TODO: Add back evaporative demand calculations
  type(vegn_cohort_type), intent(inout) :: cohort
  real, intent(in)  :: PAR_dn   ! downward PAR at the top of the canopy, W/m2 
  real, intent(in)  :: PAR_net  ! net PAR absorbed by the canopy, W/m2
  real, intent(in)  :: cana_T   ! temperature of canopy air, deg K
  real, intent(in)  :: cana_q   ! specific humidity in canopy air space, kg/kg
  real, intent(in)  :: cana_co2 ! co2 concentration in canopy air space, mol CO2/mol dry air
  real, intent(in)  :: p_surf   ! surface pressure
  real, intent(in)  :: drag_q   ! drag coefficient for specific humidity
  real, intent(in)  :: soil_beta
  real, intent(in)  :: soil_water_supply ! max supply of water to roots per unit
                                ! active root biomass per second, kg/(indiv s)
  real, intent(in)  :: con_v_v  ! one-sided foliage-CAS conductance per unit ground area        
  real, intent(out) :: stomatal_cond ! stomatal conductance, m/s

  ! ---- local vars
  real    :: water_supply ! water supply per m2 of leaves
  real    :: fw, fs ! wet and snow-covered fraction of leaves
  real    :: psyn   ! net photosynthesis, mol C/(m2 of leaves s)
  real    :: resp   ! leaf respiration, mol C/(m2 of leaves s)
  real    :: stomatal_cond_limited

  select case (vegn_phot_option)

  case(VEGN_PHOT_SIMPLE)
     ! beta non-unity only for "beta" models
     stomatal_cond = soil_beta / (cohort%rs_min  + (1-soil_beta)/drag_q)
     cohort%An_op  = 0
     cohort%An_cl  = 0
     psyn = 0
     resp = 0

  case(VEGN_PHOT_LEUNING)
     if(cohort%lai > 0) then
        ! recalculate the water supply to mol H20 per m2 of leaf per second
        water_supply = soil_water_supply/(mol_h2o*cohort%leafarea)
      
        call get_vegn_wet_frac (cohort, fw=fw, fs=fs)
        call gs_Leuning(PAR_dn, PAR_net, cohort%Tv, cana_q, cohort%lai, &
             cohort%leaf_age, p_surf, water_supply, cohort%species, &
             cana_co2, cohort%extinct, fs+fw, cohort%layer, &
             ! output:
             stomatal_cond, psyn, resp )
        if (water_stress_option == WSTRESS_HYDRAULICS) then
           ! stomatal cond is wetness-modified per-leaf-area 
           stomatal_cond_limited = stomatal_cond
           call vegn_hydraulics(cohort, cana_T, cana_q, p_surf, &
                                con_v_v, stomatal_cond_limited)
           psyn = psyn*stomatal_cond_limited/stomatal_cond
           stomatal_cond = stomatal_cond_limited
        endif
        ! store the calculated photosynthesis and photorespiration for future use
        ! in carbon_int
        cohort%An_op  = psyn * seconds_per_year
        cohort%An_cl  = resp * seconds_per_year
        ! convert stomatal conductance from units per unit area of leaf to the 
        ! units per unit area of cohort
        stomatal_cond = stomatal_cond*cohort%lai
     else
        ! no leaves means no photosynthesis and no stomatal conductance either
        cohort%An_op  = 0
        cohort%An_cl  = 0
        stomatal_cond = 0
     endif

  case default
     call error_mesg('vegn_stomatal_cond', &
          'invalid vegetation photosynthesis option', FATAL)
  end select

end subroutine vegn_photosynthesis


! ============================================================================
subroutine gs_Leuning(rad_top, rad_net, tl, ea, lai, leaf_age, &
                   p_surf, ws, pft, ca, kappa, leaf_wet, layer, &
                   gs, apot, acl)
  real,    intent(in)    :: rad_top ! PAR dn on top of the canopy, w/m2
  real,    intent(in)    :: rad_net ! PAR net on top of the canopy, w/m2
  real,    intent(in)    :: tl   ! leaf temperature, degK
  real,    intent(in)    :: ea   ! specific humidity in the canopy air (?), kg/kg
  real,    intent(in)    :: lai  ! leaf area index
  real,    intent(in)    :: leaf_age ! age of leaf since budburst (deciduous), days
  real,    intent(in)    :: p_surf ! surface pressure, Pa
  real,    intent(in)    :: ws   ! water supply, mol H2O/(m2 of leaf s)
  integer, intent(in)    :: pft  ! species
  real,    intent(in)    :: ca   ! concentration of CO2 in the canopy air space, mol CO2/mol dry air
  real,    intent(in)    :: kappa! canopy extinction coefficient (move inside f(pft))
  real,    intent(in)    :: leaf_wet ! fraction of leaf that's wet or snow-covered
  integer, intent(in)    :: layer ! canopy layer of the cohort
  ! note that the output is per area of leaf; to get the quantities per area of
  ! land, multiply them by LAI
  real,    intent(out)   :: gs   ! stomatal conductance, m/s
  real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
  real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)

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
  real :: b;
  real :: ds;  ! humidity deficit, kg/kg
  real :: hl;  ! saturated specific humidity at the leaf temperature, kg/kg
  real :: do1;
  
  ! miscellaneous
  real :: dum2;
  real, parameter :: light_crit = 0;
  real, parameter :: gs_lim = 0.25;

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
  ! soil water stress
  real :: Ed,an_w,gs_w;

  if (is_watch_point()) then
     write(*,*) '####### gs_leuning input #######'
     __DEBUG2__(rad_top, rad_net)
     __DEBUG1__(tl)
     __DEBUG1__(ea)
     __DEBUG1__(lai)
     __DEBUG1__(leaf_age)
     __DEBUG1__(p_surf)
     __DEBUG1__(ws)
     __DEBUG1__(pft)
     __DEBUG1__(ca)
     __DEBUG1__(kappa)
     __DEBUG1__(leaf_wet)
     write(*,*) '####### end of ### gs_leuning input #######'
  endif

  b=0.01;
  do1=0.09 ; ! kg/kg
  if (spdata(pft)%form == FORM_GRASS) do1=0.15;


  ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
  ! empirical relationship from McCree is light=rn*0.0000046
  light_top = rad_top*rad_phot;
  par_net   = rad_net*rad_phot;
  
  ! calculate humidity deficit, kg/kg
  call qscomp(tl, p_surf, hl)
  ds = max(hl-ea,0.0)

  ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986

  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
  vm=spdata(pft)%Vmax*exp(3000.0*(1.0/288.2-1.0/tl));
  !decrease Vmax due to aging of temperate deciduous leaves 
  !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
  if (spdata(pft)%leaf_age_tau>0 .and. leaf_age>spdata(pft)%leaf_age_onset) then
     vm=vm*exp(-(leaf_age-spdata(pft)%leaf_age_onset)/spdata(pft)%leaf_age_tau)
  endif

  capgam=0.5*kc/ko*0.21*0.209; ! Farquhar & Caemmerer 1982

  ! Find respiration for the whole canopy layer
  
  Resp=spdata(pft)%gamma_resp*vm*lai/layer
  Resp=Resp/((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
  
  
  ! ignore the difference in concentrations of CO2 near
  !  the leaf and in the canopy air, rb=0.
 
  Ag_l=0.;
  Ag_rb=0.;
  Ag=0.;
  anbar=-Resp/lai;
  gsbar=b;

 
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
           Ag_l   = spdata(pft)%alpha_phot * par_net &
                * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
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
           Ag_l   = spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * par_net &
                * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
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
  
  an_w=anbar;
  if (an_w > 0.) then
     an_w=an_w*(1-spdata(pft)%wet_leaf_dreg*leaf_wet);
  endif
  
  gs_w=gsbar*(1-spdata(pft)%wet_leaf_dreg*leaf_wet);

  if (gs_w > gs_lim) then
      if(an_w > 0.) an_w = an_w*gs_lim/gs_w;
      gs_w = gs_lim;
  endif

  if (water_stress_option == WSTRESS_LM3) then
    ! find water availability
    ! diagnostic demand

    Ed=gs_w*ds*mol_air/mol_h2o; ! mol/m2/s
    ! the factor mol_air/mol_h2o makes units of gs_w and humidity deficit ds compatible:
    ! ds*mol_air/mol_h2o is the humidity deficit in [mol_h2o/mol_air]

    if (Ed>ws) then
       w_scale=ws/Ed;
       gs_w=w_scale*gs_w;
       if(an_w > 0.0) an_w = an_w*w_scale;
       if(an_w < 0.0.and.gs_w >b) gs_w=b;
       if (is_watch_point()) then
          write(*,*)'#### gs is water-limited'
          __DEBUG4__(w_scale, gs_w, an_w, b)
       endif
    endif
  endif
  ! convert units of stomatal conductance to m/s from mol/(m2 s) by
  ! multiplying it by a volume of a mole of gas
  gs   = gs_w * Rugas * Tl / p_surf
  apot = an_w
  acl  = -Resp/lai

  if (is_watch_point()) then
     __DEBUG3__(gs, apot, acl)
  endif
end subroutine gs_Leuning


! ==============================================================================
subroutine vegn_hydraulics(cc, cana_T, cana_q, p_surf, con_v_v, gs_w)
  type(vegn_cohort_type), intent(inout) :: cc
  real,    intent(in)    :: cana_T   ! temperature in the canopy air C
  real,    intent(in)    :: cana_q   ! specific humidity in the canopy air space (common to all cohorts), kg/kg
  real,    intent(in)    :: p_surf   ! surface pressure, N/m2 = kg/m/s2 = Pa
  real,    intent(in)    :: con_v_v  ! cohort-specific one-sided foliage-CAS conductance per unit ground area, m/s
  real,    intent(inout) :: gs_w     ! stomatal conductance, per leaf area, wetness already accounted m/s

  ! ---- local vars     
  real :: rho ! kg/m3
  real :: gb, gs ! kg/m2/s
  real :: rb, rs ! m2-s/kg  
  real :: vegn_T, vegn_q, surf_q, dq ! K, kg
  real :: Kra, Kxa, Kxam, Kla  ! K_a are units per tissue area kg / m2 tissue / s / Pa
  real :: Kri, Kxi, Kxim, Kli  ! K_i are units individual kg / s / Pa
  real :: psi_rs, psi_r, psi_x, psi_l, psi_l_min, psi_l_min_prior  ! state variables Pa = kg/m/s2
  real :: psi_bar ! psi used in stress calcs. Units m, for comparison to spdata in m
  real :: cx, cl  ! exponents to Weibull (unitless)
  real :: dx, dl  ! denominators to Weibull Pa
  real :: psi_xs, psi_ls, psi_tlp ! paramaters (s=star) Pa = kg/m/s2
  real :: uro, uxo, ulo, EAo, uco, u ! flows kg/indiv/s
  real :: dur_dPR, dux_dPX, dux_dPR, dul_dPL, dul_dPX, dEA_dPL, dEA_dqs, duc_dqs ! derivatives wrt state variables 
  real :: Kpi, dpsi ! diagnostics: whole-plant K and expected dPsi as E/K

  real :: rootarea, stemarea ! tissue areas, m2/indiv

  ! factors to put units of A into mmol/MPa instead of kg/Pa
!  real, parameter :: a1 = mmol_per_kg_h2o*1e6 ! (mmol/MPa) = (kg/Pa)(mmol/kg)(Pa/MPa)
!  real, parameter :: a2 = mmol_per_kg_h2o/1e3  ! (mmol/(mmol/mol)) = (kg/(kg/kg)(mmol/kg)(mmol/kg)(kg/mmol)(mol/mmol)
  
  real, parameter :: eps = 0.0001 ! convergence on flux, kg/s  ! Sperry 98 uses 0.0018kg/s
  real, parameter :: Pa_per_m = dens_h2o*grav ! Pascal per meter of water
  integer, parameter :: niter = 1 ! iterations before checking on convergence
  integer :: i ! iterator 
  real, parameter :: dfac = 1 ! for update, slows down convergence
  real :: crit ! convergence test
  real :: DY(4,1), DU(4,1), A(4,4), AINV(4,4) ! matrices for flows
  logical :: OK_FLAG ! req'd for matrix inversion

  associate(sp=>spdata(cc%species)) 

  ! calculate humidity deficit, kg/kg
  vegn_T = cc%Tv ! units K
  call qscomp(vegn_T, p_surf, vegn_q) !
  dq = max(vegn_q-cana_q,0.0)
  
  if(vegn_q <= cana_q .or. vegn_T <= T_transp_min) then
     ! vegetation is too cold, or air is condensing on leaf
     return
  endif

  rootarea = cc%br * sp%srl * 2*PI * sp%root_r  ! (kg/indiv)(m/kg)(m2/m)
  stemarea = sp%alphaCSASW * cc%DBH**sp%thetaCSASW

  ! MAKE INTERMEDIATE CALCULATIONS

  ! rdgas units: J/kg/deg = kg*m2/s2 / kg / deg    
  rho = p_surf/(rdgas*cana_T *(1+d608*cana_q)) ! units kg/m3

  ! conductance per individual kg (air)/s
  ! (kg/indiv/s) = (kg/m2leaf/s)(m2leaf/m2ground)(m2ground/indiv)
  gb = con_v_v*rho*cc%crownarea ! foliage-CAS conductance per m2 => per individual. 
                             ! conceptual model is all leaves put moisture into a volume, and that volume exchanges with pbl
                             ! in proportion to crown area.
  ! (kg/indiv/s) = (m/s)(kg/m3)(m2leaf/indiv)
  gs = gs_w*rho*cc%leafarea  ! stomatal conductance per m2 leaf => per individual
    
  if (hydraulics_repair) then
    ! take hydraulics from species data and cohort biomass allocation
    ! head (m) in denominator
    Kri = sp%Kram * rootarea / Pa_per_m 
    Kxi = sp%Kxam / cc%height * stemarea / Pa_per_m 
    Kli = sp%Klam * cc%leafarea / Pa_per_m 
  else
    ! take hydraulics from state variables
    Kri = cc%Kra * rootarea / Pa_per_m 
    Kxi = cc%Kxa / cc%height * stemarea / Pa_per_m 
    Kxim = sp%Kxam / cc%height * stemarea / Pa_per_m 
    Kli = cc%Kla * cc%leafarea / Pa_per_m 
    psi_xs = sp%dx*log(cc%Kxa/sp%Kxam)**(1./sp%cx) / Pa_per_m ! invert Weibull to find psi where Kxa is equal to K(psi)
  endif

  dx = sp%dx * Pa_per_m     
  dl = sp%dl * Pa_per_m 
  cx = sp%cx
  cl = sp%cl

  psi_rs = min(cc%psi_rs * Pa_per_m, 0.0)
  psi_l  = cc%psi_l * Pa_per_m

  ! Assume prior psi_l is plausible and compute initial ET demand, thus initial potentials 
  rb = 1./gb
  rs = 1./(gs*exp(-(psi_l/dl)**cl))    
  surf_q = vegn_q - dq*(rs/(rs+rb))

  EAo = gs*exp(-(psi_l/dl)**cl)*(vegn_q - surf_q)    
  psi_r = psi_rs - EAo/Kri
  psi_x = psi_r  - EAo/(Kxi*exp(-(psi_r /dx)**cx))
  psi_l = psi_x  - EAo/(Kli*exp(-(psi_x /dl)**cl))

  ! Diagnostics
  Kpi = 1 / ( 1/Kri + 1/(Kxi*exp(-(psi_r/dx)**cx)) + 1/(Kli*exp(-(psi_x/dl)**cl)) )
  dpsi = EAo/Kpi

  ! ENTER WHILE LOOP
  do i = 1,niter
 
    ! CALCULATE FLUXES AND DERIVATIVES
    ! Right now this is instantly repairable. Need to modify for retained damage
    ! Is it as simple as changing KRM to KR etc?

    ! LHS: FLUXES   
    ! note the premult by gamma; conventionally incomplete gammas are normalized by 1./gamma()
    ! Eq 1

    uro = Kri*(psi_rs-psi_r)
    
    ! Eq 2
    if (hydraulics_repair) then 
      uxo = -Kxi*dx/cx*gamma(1./cx)*(gammaU((psi_r/dx)**cx, 1./cx) - gammaU((psi_x/dx)**cx, 1./cx))
    else
      ! Piecewise integration over clamped PLC curve
      if (psi_r .le. psi_xs .and. psi_x .lt. psi_xs) then
        ! Only in curved range
        uxo = -Kxim*dx/cx*gamma(1./cx)*(gammaU((psi_r/dx)**cx, 1./cx) - gammaU((psi_x/dx)**cx, 1./cx))
      elseif (psi_r .gt. psi_xs .and. psi_x .ge. psi_xs) then
        ! Only in constant range
        uxo = -Kxi*(psi_x - psi_xs)
      else  
        ! In curved and constant range
      uxo = -Kxim*dx/cx*gamma(1./cx)*(gammaU((psi_r/dx)**cx, 1./cx) - gammaU((psi_xs/dx)**cx, 1./cx)) &
        - Kxi*(psi_x - psi_xs)
      endif
    endif

    ! Eq 3-5
    ulo = -Kli*dl/cl*gamma(1./cl)*(gammaU((psi_x/dl)**cl, 1./cl) - gammaU((psi_l/dl)**cl, 1./cl))
    EAo = gs*exp(-(psi_l/dl)**cl)*(vegn_q - surf_q)
    uco = gb*(surf_q-cana_q)
  
    !RHS: DERIVATIVES WRT UNKNOWNS
    !eq1
    dur_dPR =   Kri

    !eq2
    if (hydraulics_repair) then 
      dux_dPX =  Kxi*exp(-(psi_x/dx)**cx)
      dux_dPR = -Kxi*exp(-(psi_r/dx)**cx)
    else
      ! Piecewise integration over clamped PLC curve
      if (psi_r .lt. psi_xs) then
        dux_dPR = -Kxim*exp(-(psi_r/dx)**cx)
      else
        dux_dPR = -Kxi
      endif
      if (psi_x .lt. psi_xs) then
        dux_dPX =  Kxim*exp(-(psi_x/dx)**cx)
      else
        dux_dPX =  Kxi
      endif
    endif
  
    !eq3
    dul_dPL =  Kli*exp(-(psi_l/dl)**cl)
    dul_dPX = -Kli*exp(-(psi_x/dl)**cl)
  
    !eq4 
    dEA_dPL = -gs*(cl/dl)*(psi_l/dl)**(cl - 1.)*exp(-(psi_l/dl)**cl)*(vegn_q-surf_q)    
    dEA_dqs = -gs*exp(-(psi_l/dl)**cl)    
  
    !eq5 = 
    duc_dqs =  gb
  
    ! Solve update to Psi vector
    ! note F90 is column-major so the rows specified here need to be transposed to be actual rows
    A = transpose(reshape( &
        (/ dur_dPR-dux_dPR,        -dux_dPX, 0              , 0                   , &
           dux_dPR        , dux_dPX-dul_dPX,        -dul_dPL, 0                   , &
           0              , dul_dPX        , dul_dPL-dEA_dPL,        -dEA_dqs  , &
           0              , 0              , dEA_dPL        , dEA_dqs-duc_dqs /), shape(A)))
  
    DU = reshape((/ uxo - uro, &
            ulo - uxo, &
            EAo - ulo, &
            uco - EAo /), shape(DU))
   
    ! perhaps easier to do LU decomp
    ! see update land model fast
    ! remember to declare indx (size of DU), ierr
    ! call ludcmp(A,indx,ierr) -- decompses A in place
    ! DY = DU
    ! call lubksb(A,indx,DY) -- ! replaces DU with solution DY      

    call m44inv(A, AINV, OK_FLAG)
    ! recommend to implement by hand
    DY(1,1) = DU(1,1)*AINV(1,1) + DU(2,1)*AINV(1,2)
    DY(2,1) = DU(1,1)*AINV(2,1) + DU(2,1)*AINV(2,2) + DU(3,1)*AINV(2,3)
    DY(3,1) =                     DU(2,1)*AINV(3,2) + DU(3,1)*AINV(3,3) + DU(4,1)*AINV(3,4)
    DY(4,1) =                                         DU(3,1)*AINV(4,3) + DU(4,1)*AINV(4,4)

    !check for convergence
    crit = abs(DY(1,1)/psi_r) + abs(DY(2,1)/psi_x) + abs(DY(3,1)/psi_l) + abs(DY(4,1)/surf_q)

    __DEBUG2__(psi_l,DY(3,1))
    !update. dfac slows convergence
    psi_r  = psi_r  + DY(1,1)/dfac
    psi_x  = psi_x  + DY(2,1)/dfac
    psi_l  = psi_l  + DY(3,1)/dfac
    surf_q = surf_q + DY(4,1)/dfac
    __DEBUG1__(psi_l)
    
    ! irrelevant when niter == 1
    if (crit .le. eps) exit

  enddo ! solver loop

  ! POSTPROCESS SOLUTION
  ! all U equal at solution, but in a single iteration perhaps not. 
  ! check which term is longest to converge
  u = -Kli*dl/cl*gamma(1./cl)*(gammaU((psi_x/dl)**cl, 1./cl) - gammaU((psi_l/dl)**cl, 1./cl))
  if (u .gt. 0.) then
     if (psi_rs - psi_r .gt. 0.) Kri = u/(psi_rs - psi_r)
     if (psi_r - psi_x .gt. 0.) Kxi = u/(psi_r - psi_x)
     if (psi_x - psi_l .gt. 0.) Kli = u/(psi_x - psi_l)
  endif

  ! UPDATE STATE VARIABLES IN COHORT
  ! return all pressure units to m
  ! pressure in numerator
  cc%psi_r = psi_r / Pa_per_m
  cc%psi_x = psi_x / Pa_per_m
  cc%psi_l = psi_l / Pa_per_m
  ! pressure in denominator
  cc%Kxi = Kxi * Pa_per_m
  cc%Kxa = Kxi * cc%height / stemarea * Pa_per_m
  cc%Kli = Kli * Pa_per_m
  cc%Kla = Kli / cc%leafarea * Pa_per_m 

  ! Kra is defined here each timestep - be careful if hydraulics_repair is .true.
  ! This can lead to unexpected results even if do_virtual_hydraulics .eq. .true.
  ! Kra is used in all soil uptake solvers, while Kri is used to get psi_rs in soil_step_2
  if (hydraulics_repair) then
    cc%Kri = sp%Kram * rootarea ! stays in meters
    cc%Kra = sp%Kram
  else
    cc%Kri = Kri * Pa_per_m
    cc%Kra = cc%Kri / rootarea
  endif


  ! Update gs_w for output
  !gs_w = gs_w*exp(-psi_l/psi_ls)
  gs_w = gs_w*exp(-(psi_l/dl)**cl)    ! gs now same Weibull with same parms as Kl
 
  ! TODO: make sure this honors the minimum cuticular conductance


  ! CANOPY DIEBACK
  ! TODO: either completely disable or fix this code. As is, it doesn't even make
  ! a try at carbon conservation. Apparently, Adam has not tested it, perhaps
  ! because he says he only ran the code on one point (where I think wilting is 
  ! not likely)
  psi_tlp = sp%psi_tlp * Pa_per_m 
  psi_l_min = psi_l - (psi_x - psi_l) ! assume uniform distribution, with psi_l as average and psi_x as max
  if (psi_l_min .lt. psi_tlp) then
    cc%bl = cc%bl - cc%bl_wilt ! kill previously wilted
    cc%bl_wilt = cc%bl * (psi_l_min - psi_tlp)/(psi_l_min - psi_x) ! add additional wilt if need be. psi_x cannot get below psi_tlp

                    !  min  tlp L       X
                    !   |=======|=======|
                    !        |==========|
                    !   |wilt| unwilted |
  else    
    if (cc%bl_wilt .gt. 0.) then !restore previously wilted to the extent possible
      psi_l_min_prior = psi_x - (psi_l_min - psi_x)*(cc%bl + cc%bl_wilt)/cc%bl
      cc%bl_wilt = cc%bl_wilt * (psi_tlp - psi_l_min)/(psi_l_min_prior - psi_l_min)

                    ! minp tlp   min   L    X
                    !             |====|====|
                    !       |===============|
                    !  |wilt|recov| unwiltd |
    endif
  endif

  end associate

end subroutine vegn_hydraulics

end module vegn_photosynthesis_mod
