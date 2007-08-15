module vegn_photosynthesis_mod

use fms_mod,            only : write_version_number, error_mesg, FATAL
use constants_mod,      only : TFREEZE 
use sphum_mod,          only : qscomp

use land_constants_mod, only : BAND_VIS, Rugas
use land_debug_mod,     only : is_watch_point
use vegn_data_mod,      only : MSPECIES, PT_C4, spdata
use vegn_tile_mod,      only : vegn_tile_type
use vegn_cohort_mod,    only : vegn_cohort_type
use cohort_list_mod,    only : first_cohort, current_cohort

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_photosynthesis_init
public :: vegn_photosynthesis
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: ', &
   tagname = '$Name: ' ,&
   module_name = 'vegn_photosynthesis'
! values for internal vegetation radiation option selector
integer, parameter :: VEGN_PHOT_SIMPLE  = 1 ! zero photosynthesis
integer, parameter :: VEGN_PHOT_LEUNING = 2 ! photosynthesis according to simplified Leuning model

! ==== module variables ======================================================
integer :: vegn_phot_option = -1 ! selector of the photosynthesis option


contains


! ============================================================================
subroutine vegn_photosynthesis_init(photosynthesis_to_use)
  character(*), intent(in) :: photosynthesis_to_use

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

end subroutine vegn_photosynthesis_init


! ============================================================================
! compute stomatal conductance, photosynthesis and respiration
subroutine vegn_photosynthesis ( vegn, &
     PAR_dn, PAR_net, cana_q, cana_co2, p_surf, drag_q, &
     soil_beta, soil_water_supply,&
     stomatal_cond, psyn, resp )
  type(vegn_tile_type), intent(in) :: vegn
  real, intent(in)  :: PAR_dn   ! downward PAR at the top of the canopy, W/m2 
  real, intent(in)  :: PAR_net  ! net PAR absorbed by the canopy, W/m2
  real, intent(in)  :: cana_q   ! specific humidity in canopy air space, kg/kg
  real, intent(in)  :: cana_co2 ! co2 concentration in canopy air space, ?
  real, intent(in)  :: p_surf   ! surface pressure
  real, intent(in)  :: drag_q   ! drag coefficient for specific humidity
  real, intent(in)  :: soil_beta
  real, intent(in)  :: soil_water_supply ! max supply of water to roots per unit
                                ! active root biomass per second, kg/(m2 s)
  real, intent(out) :: stomatal_cond ! stomatal conductance, m/s(?)
  real, intent(out) :: psyn     ! net photosynthesis, mol C/(m2 s)
  real, intent(out) :: resp     ! leaf respiration, mol C/(m2 s)


  ! ---- local constants
  real, parameter :: res_scaler = 20.0    ! scaling factor for water supply
  real, parameter :: mol_H20    = 18.0e-3 ! molar mass of H20, kg

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cohort
  integer :: sp ! vegetation species
  real    :: leaf_abs ! leaf absorption coefficient for PAR
  real    :: kappa    ! canopy extinction coefficient
  real    :: water_supply ! water supply per m2 of leaves

  ! get the pointer to the first (and, currently, the only) cohort
  cohort => current_cohort(first_cohort(vegn%cohorts))

  select case (vegn_phot_option)

  case(VEGN_PHOT_SIMPLE)
     ! beta non-unity only for "beta" models
     stomatal_cond = soil_beta / (cohort%rs_min  + (1-soil_beta)/drag_q)
     psyn = 0
     resp = 0

  case(VEGN_PHOT_LEUNING)
     if(cohort%lai > 0) then
        ! assign species type to local var, purely for convenience 
        sp = cohort%species
        ! calculate leaf absorption coefficient
        leaf_abs = 1 - spdata(sp)%leaf_refl(BAND_VIS) - spdata(sp)%leaf_tran(BAND_VIS)
        ! calculate total water supply, taking into account fine root biomass and
        ! fine root diameter
        water_supply = soil_water_supply * res_scaler * spdata(sp)%dfr * cohort%br
        ! recalculate the water supply to mol H20 per m2 of leaf per second
        water_supply = water_supply/(mol_H20*cohort%lai)
        call gs_Leuning(PAR_dn, PAR_net, cohort%prog%Tv, cana_q, cohort%lai, &
             p_surf, water_supply, sp, cana_co2, &
             cohort%extinct, leaf_abs, stomatal_cond, psyn, resp, cohort%pt)
        ! convert stomatal conductance, photosynthesis and leaf respiration from units
        ! per unit area of leaf to the units per unit area of land
        stomatal_cond = stomatal_cond*cohort%lai
        psyn          = psyn         *cohort%lai
        resp          = resp         *cohort%lai
     else
        ! no leaves means no photosynthesis and no stomatal conductance either
        stomatal_cond = 0
        psyn          = 0
        resp          = 0
     endif

  case default
     call error_mesg('vegn_stomatal_cond', &
          'invalid vegetation photosynthesis option', FATAL)
  end select

end subroutine vegn_photosynthesis


! ============================================================================
subroutine gs_Leuning(rad_top, rad_net, tl, ea, lai, &
                   p_surf, ws, pft, ca, &
                   kappa, leaf_abs,  &
                   gs, apot, acl, pt)
  real,    intent(in)    :: rad_top ! PAR dn on top of the canopy, w/m2
  real,    intent(in)    :: rad_net ! PAR net on top of the canopy, w/m2
  real,    intent(in)    :: tl   ! leaf temperature, degK
  real,    intent(in)    :: ea   ! specific humidity in the canopy air (?), kg/kg
  real,    intent(in)    :: lai  ! leaf area index

  real,    intent(in)    :: p_surf ! surface pressure, Pa
  real,    intent(in)    :: ws   ! water supply, mol H20/(m2 of leaf s)
  integer, intent(in)    :: pft  ! species
  real,    intent(in)    :: ca   ! concentartion of CO2 in the canopy air space
  real,    intent(in)    :: kappa! canopy extinction coefficient (move inside f(pft))
  real,    intent(in)    :: leaf_abs ! leaf absorption (move inside f(pft))
  ! note that the output is per area of leaf; to get the quantities per area of
  ! land, multiply them by LAI
  real,    intent(out)   :: gs   ! stomatal conductance, m/s
  real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
  real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)
  integer, intent(in)    :: pt   ! physiology type (C3 or C4)

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
  real :: ds;
  real :: hl;  ! saturated specific humidity at the leaf temperature, kg/kg
  real :: do1;
  
  ! misceleneous
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
#define __DEBUG__(x) write(*,*) #x , x
     write(*,*) '####### gs_leuning input #######'
     __DEBUG__(rad_top) 
     __DEBUG__(rad_net)
     __DEBUG__(tl)
     __DEBUG__(ea)
     __DEBUG__(lai)
     __DEBUG__(p_surf)
     __DEBUG__(ws)
     __DEBUG__(pft)
     __DEBUG__(ca)
     __DEBUG__(kappa)
     __DEBUG__(leaf_abs)
     __DEBUG__(pt)
     write(*,*) '####### end of ### gs_leuning input #######'
  endif

  b=0.01;
  do1=0.09 ; ! kg/kg
  if (pft < 2) do1=0.15;


  ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
  ! empirical relationship from McCree is light=rn*0.0000046
  light_top = rad_top*rad_phot;
  par_net   = rad_net*rad_phot;
  
  ! calculate humidity deficit
  call qscomp(tl, p_surf, hl)
  ds = max(hl-ea,0.0)

  ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986

  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
  vm=spdata(pft)%Vmax*exp(3000.0*(1.0/288.2-1.0/tl));
  capgam=0.5*kc/ko*0.21*0.209; ! Farquhar & Caemmerer 1982

  ! Find respiration for the whole canopy layer
  
  Resp=spdata(pft)%gamma_resp*vm*lai;
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
     if (pt==PT_C4) then ! C4 species
        coef0=(1+ds/do1)/spdata(pft)%m_cond;
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
        if (ci>capgam) then
           ! if (C4)
           !   f1=alpha4*light*kappa*1.00;
           f2=vm;
           f3=18000.0*vm*ci;       
           ! /*if (f1<f2) an=f1;
           !  else an=f1;
           !  if (f3<an) an=f3;*/
       
           dum2=min(f2,f3)
           
           ! find LAI level at which rubisco limited rate is equal to light limited rate
           lai_eq=-log(dum2/(kappa*spdata(pft)%alpha_phot*light_top))/kappa;
           if (lai_eq > lai) then
              Ag_l  = 0.0; 
              Ag_rb = dum2*lai;
           else
              if (lai_eq > 0.01) then
                 ! solve for both gross photosynthesis
                 Ag_l = spdata(pft)%alpha_phot*par_net*exp(-kappa*lai_eq);
                 Ag_rb= dum2*lai_eq;        
              else     
                 Ag_rb= 0.0;
                 Ag_l = spdata(pft)%alpha_phot*par_net;
              endif
           endif
           
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
                       (spdata(pft)%alpha_phot*light_top*leaf_abs))/kappa;
           if (lai_eq > lai) then
              ! all rubisco or transport limited
              Ag_l=0;
              Ag_rb=dum2*lai;          
           else 
              if(lai_eq > 0.000001) then
                 ! solve for both rates
                 Ag_l=spdata(pft)%alpha_phot*(ci-capgam)/(ci+2.*capgam)*par_net*exp(-kappa*lai_eq);
                 Ag_rb=dum2*lai_eq;           
              else ! only light limited rate 
                 Ag_rb=0;
                 Ag_l=spdata(pft)%alpha_phot*(ci-capgam)/(ci+2.*capgam)*par_net;
              endif
           endif
          
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
  gs_w=gsbar;

  if (gs_w > gs_lim) then
      if(an_w > 0.) an_w = an_w*gs_lim/gs_w;
      gs_w = gs_lim;
  endif

#if 1
  ! find water availability
  ! diagnostic demand

  Ed=gs_w*ds;

  if (Ed>ws) then
     w_scale=ws/Ed;
     gs_w=w_scale*gs_w;
     if(an_w > 0.0) an_w = an_w*w_scale;
     if(an_w < 0.0.and.gs_w >b) gs_w=b;
  endif
  gs=gs_w;
  apot=an_w;
  acl=-Resp/lai;
  

#else
! no water limitation on stomata
   gs=gsbar;  
   apot=anbar; 
   acl=-Resp/lai; 
#endif 

   ! finally, convert units of stomatal conductance to m/s from mol/(m2 s) by
   ! multiplying it by a volume of a mole of gas
   gs = gs * Rugas * Tl / p_surf

   if (is_watch_point()) then
      __DEBUG__(gs)
      __DEBUG__(apot)
      __DEBUG__(acl)
   endif
end subroutine gs_Leuning

end module vegn_photosynthesis_mod
