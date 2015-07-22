! ============================================================================
! updates carbon pools and rates on the fast time scale
! ============================================================================
module vegn_dynamics_mod

#include "../shared/debug.inc"

use fms_mod, only: write_version_number
use time_manager_mod, only: time_type

use constants_mod, only : PI,tfreeze
use land_constants_mod, only : seconds_per_year, mol_C
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type
use vegn_data_mod, only : spdata, &
     CMPT_VLEAF, CMPT_SAPWOOD, CMPT_ROOT, CMPT_WOOD, CMPT_LEAF, LEAF_ON, LEAF_OFF, &
     fsc_liv, fsc_wood, K1, K2, soil_carbon_depth_scale, C2B, agf_bs, &
     l_fract, mcv_min, mcv_lai, do_ppa
use vegn_tile_mod, only: vegn_tile_type, vegn_tile_carbon
use soil_tile_mod, only: soil_tile_type, soil_ave_temp, soil_ave_theta, clw, csw
use vegn_cohort_mod, only : vegn_cohort_type, &
     update_biomass_pools, update_bio_living_fraction, update_species, &
     leaf_area_from_biomass

use land_debug_mod, only : is_watch_point
use land_numerics_mod, only : rank_descending

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_dynamics_init

public :: vegn_carbon_int_lm3   ! fast time-scale integrator of carbon balance
public :: vegn_carbon_int_ppa   ! fast time-scale integrator of carbon balance
public :: vegn_growth       ! slow time-scale redistributor of accumulated carbon
public :: vegn_daily_npp    ! updates values of daily-average npp
public :: vegn_phenology    !
public :: vegn_biogeography !

public :: vegn_reproduction_ppa ! reproduction for PPA case
public :: relayer_cohorts
public :: vegn_mergecohorts_ppa ! merge cohorts
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id$', &
   tagname = '$Name$' ,&
   module_name = 'vegn'
real, parameter :: GROWTH_RESP=0.333  ! fraction of npp lost as growth respiration


! ==== module data ===========================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

! diagnostic field IDs
integer :: id_npp, id_nep, id_gpp, id_fast_soil_C, id_slow_soil_C, id_rsoil, id_rsoil_fast
integer :: id_resp, id_resl, id_resr, id_resg, id_asoil
integer :: id_soilt, id_theta, id_litter


contains

! ============================================================================
subroutine vegn_dynamics_init(id_lon, id_lat, time, delta_time)
  integer        , intent(in) :: id_lon ! ID of land longitude (X) axis 
  integer        , intent(in) :: id_lat ! ID of land latitude (Y) axis
  type(time_type), intent(in) :: time       ! initial time for diagnostic fields
  real           , intent(in) :: delta_time ! fast time step, s

  call write_version_number(version, tagname)

  ! set up global variables
  dt_fast_yr = delta_time/seconds_per_year

  ! register diagnostic fields
  id_gpp = register_tiled_diag_field ( module_name, 'gpp',  &
       (/id_lon,id_lat/), time, 'gross primary productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_npp = register_tiled_diag_field ( module_name, 'npp',  &
       (/id_lon,id_lat/), time, 'net primary productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_nep = register_tiled_diag_field ( module_name, 'nep',  &
       (/id_lon,id_lat/), time, 'net ecosystem productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_litter = register_tiled_diag_field (module_name, 'litter', (/id_lon,id_lat/), &
       time, 'litter productivity', 'kg C/(m2 year)', missing_value=-100.0)
  id_fast_soil_C = register_tiled_diag_field ( module_name, 'fsc',  &
       (/id_lon,id_lat/), time, 'fast soil carbon', 'kg C/m2', &
       missing_value=-100.0 )
  id_slow_soil_C = register_tiled_diag_field ( module_name, 'ssc',  &
       (/id_lon,id_lat/), time, 'slow soil carbon', 'kg C/m2', &
       missing_value=-100.0 )
  id_resp = register_tiled_diag_field ( module_name, 'resp', (/id_lon,id_lat/), &
       time, 'respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_resl = register_tiled_diag_field ( module_name, 'resl', (/id_lon,id_lat/), &
       time, 'leaf respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_resr = register_tiled_diag_field ( module_name, 'resr', (/id_lon,id_lat/), &
       time, 'root respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_resg = register_tiled_diag_field ( module_name, 'resg', (/id_lon,id_lat/), &
       time, 'growth respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_rsoil = register_tiled_diag_field ( module_name, 'rsoil',  &
       (/id_lon,id_lat/), time, 'soil respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_fast = register_tiled_diag_field ( module_name, 'rsoil_fast',  &
       (/id_lon,id_lat/), time, 'fast soil carbon respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_asoil = register_tiled_diag_field ( module_name, 'asoil',  &
       (/id_lon,id_lat/), time, 'aerobic activity modifier', &
       missing_value=-100.0 )
  id_soilt = register_tiled_diag_field ( module_name, 'tsoil_av',  &
       (/id_lon,id_lat/), time, 'average soil temperature for carbon decomposition', 'degK', &
       missing_value=-100.0 )
  id_theta = register_tiled_diag_field ( module_name, 'theta',  &
       (/id_lon,id_lat/), time, 'average soil wetness for carbon decomposition', 'm3/m3', &
       missing_value=-100.0 )
end subroutine vegn_dynamics_init

! ============================================================================
subroutine vegn_carbon_int_lm3(vegn, soilt, theta, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: soilt ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  type(vegn_cohort_type), pointer :: cc
  type(vegn_cohort_type), pointer :: c(:) ! for debug only
  real :: resp, resl, resr, resg ! respiration terms accumulated for all cohorts 
  real :: md_alive, md_wood;
  real :: md ! plant tissue maintenance, kg C/timestep
  real :: gpp ! gross primary productivity per tile
  integer :: sp ! shorthand for current cohort specie
  integer :: i

  c=>vegn%cohorts(1:vegn%n_cohorts)
  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_lm3 input ####'
     __DEBUG2__(soilt,theta)
     __DEBUG1__(c%npp_previous_day)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%carbon_loss)
     __DEBUG1__(c%bwood_gain)
  endif

  !  update plant carbon
  vegn%npp = 0
  resp = 0 ; resl = 0 ; resr = 0 ; resg = 0 ; gpp = 0
  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)
     sp = cc%species

     call plant_respiration(cc,soilt);
   
     cc%gpp = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea;
     cc%npp = cc%gpp - cc%resp;
   
     if(cc%npp_previous_day > 0) then
        cc%resg = GROWTH_RESP*cc%npp_previous_day;
        cc%npp  = cc%npp  - GROWTH_RESP*cc%npp_previous_day;
        cc%resp = cc%resp + GROWTH_RESP*cc%npp_previous_day;
     else
        cc%resg = 0;
     endif
     ! accumulate npp for the current day
     cc%npp_previous_day_tmp = cc%npp_previous_day_tmp + cc%npp; 
     
     cc%carbon_gain = cc%carbon_gain + cc%npp*dt_fast_yr;
     
     ! check if leaves/roots are present and need to be accounted in maintenance
     if(cc%status == LEAF_ON) then
        md_alive = (cc%Pl * spdata(sp)%alpha(CMPT_LEAF) + &
                    cc%Pr * spdata(sp)%alpha(CMPT_ROOT))* &
              cc%bliving*dt_fast_yr;    
     else
        md_alive = 0
     endif
     
     ! compute branch and coarse wood losses for tree types
     md_wood =0;
     if (sp > 1) then
        md_wood = 0.6 *cc%bwood * spdata(sp)%alpha(CMPT_WOOD)*dt_fast_yr;
     endif
        
     md = md_alive + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain - md_wood;
!     if (cc%bwood_gain < 0.0) cc%bwood_gain=0.0; ! potential non-conservation ?
     cc%carbon_gain = cc%carbon_gain - md;
     cc%carbon_loss = cc%carbon_loss + md; ! used in diagnostics only

     ! add maintenance demand from leaf and root pools to fast soil carbon
     vegn%fast_soil_C = vegn%fast_soil_C + (   fsc_liv *md_alive +    fsc_wood *md_wood)*cc%nindivs;
     vegn%slow_soil_C = vegn%slow_soil_C + ((1-fsc_liv)*md_alive + (1-fsc_wood)*md_wood)*cc%nindivs;

     ! for budget tracking
!/*     cp->fsc_in+= data->fsc_liv*md_alive+data->fsc_wood*md_wood; */
!/*     cp->ssc_in+= (1.- data->fsc_liv)*md_alive+(1-data->fsc_wood)*md_wood; */
     vegn%fsc_in  = vegn%fsc_in + (     1 *md_alive+   0 *md_wood)*cc%nindivs;
     vegn%ssc_in  = vegn%ssc_in + ((1.- 1)*md_alive+(1-0)*md_wood)*cc%nindivs;

     vegn%veg_in  = vegn%veg_in  + cc%npp*cc%nindivs*dt_fast_yr;
     vegn%veg_out = vegn%veg_out + (md_alive+md_wood)*cc%nindivs;

     ! accumulate tile-level NPP and GPP
     vegn%npp = vegn%npp + cc%npp*cc%nindivs
     gpp = gpp + cc%gpp*cc%nindivs
     ! accumulate respiration terms for tile-level reporting
     resp = resp + cc%resp*cc%nindivs ; resl = resl + cc%resl*cc%nindivs
     resr = resr + cc%resr*cc%nindivs ; resg = resg + cc%resg*cc%nindivs
     ! update cohort age
     cc%age = cc%age + dt_fast_yr
  enddo
  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_lm3 output ####'
     __DEBUG1__(c%species)
     __DEBUG1__(c%bl)
     __DEBUG1__(c%br)
     __DEBUG1__(c%bsw)
     __DEBUG1__(c%bwood)
     __DEBUG1__(c%An_op)
     __DEBUG1__(c%An_cl)
     __DEBUG1__(c%lai)
     __DEBUG1__(c%npp)
     __DEBUG1__(c%gpp)
     __DEBUG1__(c%resp)
     __DEBUG1__(c%resl)
     __DEBUG1__(c%resr)
     __DEBUG1__(c%resg)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%carbon_loss)
     __DEBUG1__(c%bwood_gain)
  endif

  ! update soil carbon
  call Dsdt(vegn, diag, soilt, theta)

  ! NEP is equal to NNP minus soil respiration
  vegn%nep = vegn%npp - vegn%rh

  call update_soil_pools(vegn)
  vegn%age = vegn%age + dt_fast_yr;


  ! ---- diagnostic section
  call send_tile_data(id_gpp,gpp,diag)
  call send_tile_data(id_npp,vegn%npp,diag)
  call send_tile_data(id_nep,vegn%nep,diag)
  call send_tile_data(id_litter,vegn%litter,diag)
  call send_tile_data(id_resp, resp, diag)
  call send_tile_data(id_resl, resl, diag)
  call send_tile_data(id_resr, resr, diag)
  call send_tile_data(id_resg, resg, diag)
  call send_tile_data(id_soilt,soilt,diag)
  call send_tile_data(id_theta,theta,diag)
  
end subroutine vegn_carbon_int_lm3


! ============================================================================
! fast time step, Modified by Weng 2011-10-05
subroutine vegn_carbon_int_ppa (vegn, tsoil, theta, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: tsoil ! average temperature of soil for soil carbon decomposition, deg K
  real, intent(in) :: theta ! average soil wetness, unitless
  type(diag_buff_type), intent(inout) :: diag

  type(vegn_cohort_type), pointer :: cc
  type(vegn_cohort_type), pointer :: c(:) ! for debug only
  real :: resp, resl, resr, resg ! respiration terms accumulated for each tree 
  real :: md_alive, md_wood;
  real :: gpp ! gross primary productivity per tile
  real :: cgrowth ! carbon spent during this time step for growth, kg C/individual
  real :: deltaBL, deltaBR ! leaf and fine root carbon tendencies
  integer :: sp ! shorthand for current cohort specie
  integer :: i
  real :: cmass0, cmass1 ! for debug only

  c=>vegn%cohorts(1:vegn%n_cohorts)
  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_ppa input ####'
     __DEBUG2__(tsoil,theta)
     __DEBUG1__(c%bl)
     __DEBUG1__(c%br)
     __DEBUG1__(c%bsw)
     __DEBUG1__(c%bwood)
     __DEBUG1__(c%nsc)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%carbon_loss)
     __DEBUG1__(c%bwood_gain)
     cmass0 = vegn_tile_carbon(vegn) 
  endif
  ! update plant carbon for all cohorts
  vegn%npp = 0
  resp = 0 ; resl = 0 ; resr = 0 ; resg = 0 ; gpp = 0
  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)
     sp = cc%species

     ! that was in eddy_npp_PPA
     call plant_respiration(cc,tsoil)
     cc%gpp  = (cc%An_op - cc%An_cl)*mol_C*cc%leafarea
     cc%resp = MIN(cc%nsc * 0.005/dt_fast_yr, cc%resp) ! TODO: convert to sane form
     ! calculate the amount of carbon that is spent on growth
     if (cc%nsc > 0 .AND. cc%status == LEAF_ON) then
        cgrowth = 0.002 * cc%nsc/dt_fast_yr ! TODO: convert to sane form
     else
        cgrowth = 0.0
     endif
     cc%resg = GROWTH_RESP*cgrowth ! growth respiration
     cc%resp = cc%resp + cc%resg
     cc%npp  = cc%gpp - cc%resp;
     cc%nsc  = cc%nsc + (cc%npp - cgrowth)*dt_fast_yr

     ! the C that will be added to living biomass of a tree, Weng 2011-10-05
     cc%carbon_gain = cc%carbon_gain + cgrowth*dt_fast_yr 

!    Weng, 2012-01-31
     ! update leaf and fine root biomass due to overturning, and calculate 
     ! maintenance demand
     if(cc%status == LEAF_ON) then
        deltaBL = cc%bl * spdata(sp)%alpha(CMPT_LEAF) * dt_fast_yr
        deltaBR = cc%br * spdata(sp)%alpha(CMPT_ROOT) * dt_fast_yr
     else
        deltaBL = 0; deltaBR = 0
     endif
     md_alive = deltaBL + deltaBR
     cc%bl = cc%bl - deltaBL
     cc%br = cc%br - deltaBR

     cc%carbon_loss = cc%carbon_loss + md_alive  ! used in diagnostics only

     ! compute branch and coarse wood losses for tree types
     md_wood = 0.0
     if (sp > 1) then ! TODO: add input data selector for woody/grass species
        md_wood = 0.6 *cc%bwood * spdata(sp)%alpha(CMPT_WOOD)*dt_fast_yr
     endif
     cc%bwood_gain  = cc%bwood_gain - md_wood

!    It's not necessary to calculate wood C gain here with the new allocation scheme
!     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr ! - md_wood  !
!    In the new scheme, bwood_gain is driven by the turnover of sapwood, which is driven by NPP allocated to it in turn
!     if (cc%bwood_gain < 0.0) cc%bwood_gain=0.0; ! potential non-conservation ?

     ! add maintenance demand from leaf and root pools to fast soil carbon
     vegn%fast_soil_C = vegn%fast_soil_C + (   fsc_liv *md_alive +    fsc_wood *md_wood)*cc%nindivs;
     vegn%slow_soil_C = vegn%slow_soil_C + ((1-fsc_liv)*md_alive + (1-fsc_wood)*md_wood)*cc%nindivs;

     ! accumulate tile-level NPP and GPP
     vegn%npp = vegn%npp + cc%npp * cc%nindivs
     gpp      = gpp      + cc%gpp * cc%nindivs

     ! accumulate respiration terms for tile-level reporting
     resp = resp + cc%resp*cc%nindivs ; resl = resl + cc%resl*cc%nindivs
     resr = resr + cc%resr*cc%nindivs ; resg = resg + cc%resg*cc%nindivs
  enddo 

  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int_ppa output ####'
     __DEBUG1__(c%species)
     __DEBUG1__(c%bl)
     __DEBUG1__(c%br)
     __DEBUG1__(c%bsw)
     __DEBUG1__(c%bwood)
     __DEBUG1__(c%nsc)
     __DEBUG1__(c%An_op)
     __DEBUG1__(c%An_cl)
     __DEBUG1__(c%lai)
     __DEBUG1__(c%npp)
     __DEBUG1__(c%gpp)
     __DEBUG1__(c%resp)
     __DEBUG1__(c%resg)
     __DEBUG1__(c%carbon_gain)
     __DEBUG1__(c%bwood_gain)
     cmass1 = vegn_tile_carbon(vegn)
     __DEBUG3__(cmass1-cmass0-(gpp-resp)*dt_fast_yr,cmass1,cmass0)
     
     write(*,*)'#### end of vegn_carbon_int_ppa output ####'
  endif

  ! update soil carbon
  call Dsdt(vegn, diag, tsoil, theta)

  ! NEP is equal to NNP minus soil respiration
  vegn%nep = vegn%npp - vegn%rh

  call update_soil_pools(vegn)
  vegn%age = vegn%age + dt_fast_yr;

! ------ diagnostic section
  call send_tile_data(id_gpp,gpp,diag)
  call send_tile_data(id_npp,vegn%npp,diag)
  call send_tile_data(id_nep,vegn%nep,diag)
  call send_tile_data(id_litter,vegn%litter,diag)
  call send_tile_data(id_resp, resp, diag)
  call send_tile_data(id_resl, resl, diag)
  call send_tile_data(id_resr, resr, diag)
  call send_tile_data(id_resg, resg, diag)
  call send_tile_data(id_soilt,tsoil,diag)
  call send_tile_data(id_theta,theta,diag)
!  cc => vegn%cohorts(1)
!  call send_tile_data(id_nsc_c1,cc%nsc,diag)
  
end subroutine vegn_carbon_int_ppa



! ============================================================================
! updates cohort biomass pools, LAI, SAI, and height using accumulated 
! carbon_gain and bwood_gain
subroutine vegn_growth (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: cmass0,cmass1 ! for debug only
  integer :: i

  if (is_watch_point()) then
     cmass0 = vegn_tile_carbon(vegn)
  endif

  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)

     if (do_ppa) then
        call cohort_growth_ppa(cc)
     else
        cc%bwood   = cc%bwood   + cc%bwood_gain
        cc%bliving = cc%bliving + cc%carbon_gain
        
        if(cc%bliving < 0) then
           cc%bwood    = cc%bwood+cc%bliving
           cc%bliving  = 0
           if (cc%bwood < 0) &
                cc%bwood = 0 ! in principle, that's not conserving carbon
        endif
        
        call update_biomass_pools(cc)
        
        ! reset carbon acculmulation terms
        cc%carbon_gain = 0
        cc%carbon_loss = 0
        cc%bwood_gain  = 0
     endif

     if (cc%status == LEAF_ON) then
        cc%leaf_age = cc%leaf_age + 1.0
        ! limit the maximum leaf age by the leaf time span (reciprocal of leaf 
        ! turnover rate alpha) for given species. alpha is in 1/year, factor of
        ! 365 converts the result to days.
        if (spdata(cc%species)%alpha(CMPT_LEAF) > 0) &
             cc%leaf_age = min(cc%leaf_age,365.0/spdata(cc%species)%alpha(CMPT_LEAF))
     endif
  end do

  if (is_watch_point()) then
     cmass1 = vegn_tile_carbon(vegn)
     write(*,*)"############## vegn_growth #################"
     __DEBUG3__(cmass1-cmass0, cmass1, cmass0)
  endif
end subroutine vegn_growth


subroutine cohort_growth_ppa(c)
  type(vegn_cohort_type), intent(inout) :: c
  
  real :: CSAtot ! total cross section area, m2
  real :: CSAsw  ! Sapwood cross sectional area, m2
  real :: CSAwd  ! Heartwood cross sectional area, m2
  real :: DBHwd  ! diameter of heartwood at breast height, m
  real :: deltaBL, deltaBR ! tendencies of the leaf and root biomass, kgC/individual
  real :: deltaB ! sum of deltaBL and deltaBR
  real :: BSWmax ! max sapwood biomass, kg C/individual
  real :: deltaBSW ! tendency of sapwood biomass, kgC/individual
  real :: deltaBwood ! tendency of wood biomass, kgC/individual
  real :: b
  integer :: sp ! species

  sp = c%species
  
  ! calculate max biomasses of canopy and fine roots
  if(c%layer > 1)then
      c%MLmax  = spdata(sp)%LMA   * spdata(sp)%LAImax                * c%crownarea/4
      c%MFRmax = spdata(sp)%phiRL * spdata(sp)%LAImax/spdata(sp)%SRA * c%crownarea/4
  else
      c%MLmax  = spdata(sp)%LMA   * spdata(sp)%LAImax * c%crownarea
      c%MFRmax = spdata(sp)%phiRL * spdata(sp)%LAImax/spdata(sp)%SRA*c%crownarea 
  endif

  ! from biomass_allocation_ppa
  ! TODO: what if carbon_gain is not 0, but leaves are OFF (marginal case? or
  ! typical in lm3?)
  if (c%status == LEAF_ON) then
     ! for the first day of growing season: leaves and fine roots recover from blv
     if(c%blv > 0) then
        b = (c%blv + c%bl + c%br)
        c%bl  = c%MLmax/(c%MLmax+c%MFRmax)  * b
        c%br  = c%MFRmax/(c%MLmax+c%MFRmax) * b
        c%blv = 0
     endif
     
     ! for the following days in growing season: growth is driven by NPP
     deltaBL = max(0.0, c%MLmax  - c%bl)
     deltaBR = max(0.0, c%MFRmax - c%br)
     deltaB = deltaBL+deltaBR
     if (deltaB > c%carbon_gain) then
        deltaBL = deltaBL*c%carbon_gain/deltaB
        deltaBR = deltaBR*c%carbon_gain/deltaB
     endif
     deltaBSW = c%carbon_gain - deltaBL - deltaBR

     c%bl    = c%bl  + deltaBL
     c%br    = c%br  + deltaBR
     c%bsw   = c%bsw + deltaBSW

     ! recalculate DBH for upadated biomasses
     b = MAX(0.0001,(c%bwood+c%bsw) ) 
     c%DBH    = (b / spdata(sp)%alphaBM) ** (1.0 / (spdata(sp)%thetaBM))
  endif
  ! calculate DBH, BLmax, BRmax, BSWmax, bwood_gain using allometric relationships
  ! Weng 2012-01-31 update_bio_living_fraction
  CSAsw    = spdata(sp)%alphaCSASW * c%DBH ** spdata(sp)%thetaCSASW
  CSAtot   = PI *(c%DBH/2.0)**2
  CSAwd    = MAX(0.0,CSAtot - CSAsw)
  DBHwd    = 2*SQRT(CSAwd/PI)
  BSWmax   = spdata(sp)%alphaBM * (c%DBH ** spdata(sp)%thetaBM - DBHwd ** spdata(sp)%thetaBM)

  deltaBwood = max(c%bsw - BSWmax,0.0)

  c%bwood   = c%bwood + deltaBwood + c%bwood_gain
  c%bsw     = c%bsw   - deltaBwood

  ! from Ensheng's vegn_growth_ppa
  if(c%layer == 1 .AND.                  &
     c%age > spdata(sp)%maturalage .AND. &
     c%bl > 0.75*c%MLmax .AND.           &
     c%bseed < spdata(sp)%fecundity * c%crownarea)then
       c%bseed = c%bseed + 0.005 * c%bl ! convert it to sane form: use e-folding time to describe rate
       c%bl = c%bl * (1. - 0.005)
  endif
  ! reset carbon acculmulation terms
  c%carbon_gain = 0
  c%carbon_loss = 0
  c%bwood_gain  = 0
   
end subroutine cohort_growth_ppa

! ============================================================================
subroutine Dsdt(vegn, diag, soilt, theta)
  type(vegn_tile_type), intent(inout) :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! soil temperature, deg K 
  real                , intent(in)    :: theta

  real :: fast_C_loss
  real :: slow_C_loss
  real :: A  ! decomp rate reduction due to moisture and temperature
  
  A=A_function(soilt,theta);
  
  fast_C_loss = vegn%fast_soil_C*A*K1*dt_fast_yr;
  slow_C_loss = vegn%slow_soil_C*A*K2*dt_fast_yr;
  
  vegn%fast_soil_C = vegn%fast_soil_C - fast_C_loss;
  vegn%slow_soil_C = vegn%slow_soil_C - slow_C_loss;

  ! for budget check
  vegn%fsc_out = vegn%fsc_out + fast_C_loss;
  vegn%ssc_out = vegn%ssc_out + slow_C_loss;

  ! loss of C to atmosphere and leaching
  vegn%rh =   (fast_C_loss+slow_C_loss)/dt_fast_yr;
  ! vegn%rh_fast = fast_C_loss/dt_fast_yr;

  ! accumulate decomposition rate reduction for the soil carbon restart output
  vegn%asoil_in = vegn%asoil_in + A

  ! ---- diagnostic section
  call send_tile_data(id_fast_soil_C, vegn%fast_soil_C, diag)
  call send_tile_data(id_slow_soil_C, vegn%slow_soil_C, diag)
  call send_tile_data(id_rsoil_fast, fast_C_loss/dt_fast_yr, diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)
  call send_tile_data(id_asoil, A, diag)

end subroutine Dsdt


! ============================================================================
! The combined reduction in decomposition rate as a funciton of TEMP and MOIST
! Based on CENTURY Parton et al 1993 GBC 7(4):785-809 and Bolker's copy of
! CENTURY code
function A_function(soilt, theta) result(A)
  real :: A                 ! return value, resulting reduction in decomposition rate
  real, intent(in) :: soilt ! effective temperature for soil carbon decomposition
  real, intent(in) :: theta 

  real :: soil_temp; ! temperature of the soil, deg C
  real :: Td; ! rate multiplier due to temp
  real :: Wd; ! rate reduction due to mositure

  ! coefficeints and terms used in temperaturex term
  real :: Topt,Tmax,t1,t2,tshl,tshr;

  soil_temp = soilt-273.16;

  ! EFFECT OF TEMPERATURE
  ! from Bolker's century code
  Tmax=45.0;
  if (soil_temp > Tmax) soil_temp = Tmax;
  Topt=35.0;
  tshr=0.2; tshl=2.63;
  t1=(Tmax-soil_temp)/(Tmax-Topt);
  t2=exp((tshr/tshl)*(1.-t1**tshl));
  Td=t1**tshr*t2;

  if (soil_temp > -10) Td=Td+0.05;
  if (Td > 1.) Td=1.;

  ! EFFECT OF MOISTURE
  ! Linn and Doran, 1984, Soil Sci. Amer. J. 48:1267-1272
  ! This differs from the Century Wd
  ! was modified by slm/ens based on the figures from the above paper 
  !     (not the reported function)

  if(theta <= 0.3) then
     Wd = 0.2;
  else if(theta <= 0.6) then
     Wd = 0.2+0.8*(theta-0.3)/0.3;
  else 
     Wd = exp(2.3*(0.6-theta));
  endif

  A = (Td*Wd); ! the combined (multiplicative) effect of temp and water
               ! on decomposition rates
end function A_function


! ============================================================================
subroutine plant_respiration(cc, tsoil)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: tsoil
  
  real :: tf,tfs;
  real :: r_leaf, r_vleaf, r_stem, r_root
  
  integer :: sp ! shorthand for cohort species
  sp = cc%species

  tf = exp(3000.0*(1.0/288.16-1.0/cc%prog%Tv));
  tf = tf / ( &
            (1.0+exp(0.4*(5.0-cc%prog%Tv+273.16)))*&
            (1.0+exp(0.4*(cc%prog%Tv - 273.16-45.0)))&
            )

  tfs = exp(3000.0*(1.0/288.16-1.0/tsoil));
  tfs = tfs / ( &
              (1.0+exp(0.4*(5.0-tsoil+273.16)))* &
              (1.0+exp(0.4*(tsoil - 273.16-45.0)))&
              )

  r_leaf = -mol_C*cc%An_cl*cc%leafarea;
  r_vleaf = spdata(sp)%beta(CMPT_VLEAF)   * cc%blv*tf;
  r_stem  = spdata(sp)%beta(CMPT_SAPWOOD) * cc%bsw*tf;
  r_root  = spdata(sp)%beta(CMPT_ROOT)    * cc%br*tfs;

  cc%resp = r_leaf + r_vleaf + r_stem + r_root;
  cc%resl = r_leaf;
  cc%resr = r_root;
end subroutine plant_respiration


! ============================================================================
! calculates prev. day average NPP from accumualted values
subroutine vegn_daily_npp(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  integer :: n_fast_step;
  integer :: i
  type(vegn_cohort_type), pointer :: cc

  n_fast_step = 1.0/365.0/dt_fast_yr;
  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)
     vegn%cohorts(i)%npp_previous_day=vegn%cohorts(i)%npp_previous_day_tmp/n_fast_step;
     vegn%cohorts(i)%npp_previous_day_tmp=0.0
  enddo
end subroutine vegn_daily_npp


! =============================================================================
subroutine vegn_phenology(vegn, wilt)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: wilt ! ratio of wilting to saturated water content

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc
  real :: leaf_litter,root_litter;    
  real :: theta_crit; ! critical ratio of average soil water to sat. water
  integer :: i
  
  vegn%litter = 0

  do i = 1,vegn%n_cohorts   
     cc => vegn%cohorts(i)
     
     if(is_watch_point())then
        write(*,*)'####### vegn_phenology #######'
        __DEBUG4__(vegn%theta_av, wilt, spdata(cc%species)%cnst_crit_phen, spdata(cc%species)%fact_crit_phen)
        __DEBUG1__(cc%species)
        __DEBUG2__(vegn%tc_av,spdata(cc%species)%tc_crit)
     endif
     ! if drought-deciduous or cold-deciduous species
     ! temp=10 degrees C gives good growing season pattern        
     ! spp=0 is c3 grass,1 c3 grass,2 deciduous, 3 evergreen
     ! assumption is that no difference between drought and cold deciduous
     cc%status = LEAF_ON; ! set status to indicate no leaf drop
      
     if(cc%species < 4 )then! deciduous species
        ! actually either fact_crit_phen or cnst_crit_phen is zero, enforced
        ! by logic in the vegn_data.F90
        theta_crit = spdata(cc%species)%cnst_crit_phen &
              + wilt*spdata(cc%species)%fact_crit_phen
        theta_crit = max(0.0,min(1.0, theta_crit))
        if ( (vegn%theta_av < theta_crit) &
             .or.(vegn%tc_av < spdata(cc%species)%tc_crit) ) then
           cc%status = LEAF_OFF; ! set status to indicate leaf drop 
           cc%leaf_age = 0;
           
           leaf_litter = (1.0-l_fract)*cc%bl*cc%nindivs;
           root_litter = (1.0-l_fract)*cc%br*cc%nindivs;
           
           ! add to patch litter flux terms
           vegn%litter = vegn%litter + leaf_litter + root_litter;
           
           vegn%fast_soil_C = vegn%fast_soil_C +    fsc_liv *(leaf_litter+root_litter);
           vegn%slow_soil_C = vegn%slow_soil_C + (1-fsc_liv)*(leaf_litter+root_litter);
	     
           ! vegn%fsc_in+=data->fsc_liv*(leaf_litter+root_litter);
           ! vegn%ssc_in+=(1.0-data->fsc_liv)*(leaf_litter+root_litter);
           vegn%fsc_in  = vegn%fsc_in  + leaf_litter+root_litter;
           vegn%veg_out = vegn%veg_out + leaf_litter+root_litter;
           
           cc%blv = cc%blv + l_fract*(cc%bl+cc%br);
           cc%bl  = 0.0;
           cc%br  = 0.0;
           cc%lai = 0.0;
           
           ! update state
           cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw;
           call update_bio_living_fraction(cc);   
        endif
     endif
  enddo
end subroutine vegn_phenology


! =============================================================================
subroutine vegn_biogeography(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  integer :: i

  do i = 1, vegn%n_cohorts   
     call update_species(vegn%cohorts(i), vegn%t_ann, vegn%t_cold, &
          vegn%p_ann*seconds_per_year, vegn%ncm, vegn%landuse)
  enddo
end subroutine

! =============================================================================
! The stuff below comes from she_update.c -- it looks like it belongs here, 
! since it is essentially a part of the carbon integration (update_patch_fast
! is only called immediately after carbon_int in lm3v)
! =============================================================================


! =============================================================================
subroutine update_soil_pools(vegn)
  type(vegn_tile_type), intent(inout) :: vegn
  
  ! ---- local vars
  real :: delta;

  ! update fsc input rate so that intermediate fsc pool is never
  ! depleted below zero; on the other hand the pool can be only 
  ! depleted, never increased
  vegn%fsc_rate = MAX( 0.0, MIN(vegn%fsc_rate, vegn%fsc_pool/dt_fast_yr));
  delta = vegn%fsc_rate*dt_fast_yr;
  vegn%fast_soil_C = vegn%fast_soil_C + delta;
  vegn%fsc_pool    = vegn%fsc_pool    - delta;

  ! update ssc input rate so that intermediate ssc pool is never
  ! depleted below zero; on the other hand the pool can be only 
  ! depleted, never increased
  vegn%ssc_rate = MAX(0.0, MIN(vegn%ssc_rate, vegn%ssc_pool/dt_fast_yr));
  delta = vegn%ssc_rate*dt_fast_yr;
  vegn%slow_soil_C = vegn%slow_soil_C + delta;
  vegn%ssc_pool    = vegn%ssc_pool    - delta;
end subroutine update_soil_pools


! ============================================================================
function cohort_can_reproduce(cc); logical cohort_can_reproduce
  type(vegn_cohort_type), intent(in) :: cc
  
  cohort_can_reproduce = (cc%layer == 1 .and. cc%age > spdata(cc%species)%maturalage)
end function 


! ============================================================================
! the reproduction of each canopy cohort, yearly time step
! calculate the new cohorts added in this step and states:
! tree density, DBH, woddy and living biomass
subroutine vegn_reproduction_ppa (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

! ---- local vars
  type(vegn_cohort_type), pointer :: parent, cc ! parent and child cohort pointers
  type(vegn_cohort_type), pointer :: ccold(:)   ! pointer to old cohort array
  integer :: newcohorts ! number of new cohorts to be created
  integer :: i, k ! cohort indices

! Check if reproduction happens
  newcohorts = 0
  do i = 1, vegn%n_cohorts
     if (cohort_can_reproduce(vegn%cohorts(i))) &
        newcohorts=newcohorts + 1
  enddo

  if (newcohorts == 0) return ! do nothing if no cohorts are ready for reproduction

  ccold => vegn%cohorts ! keep old cohort information
  allocate(vegn%cohorts(vegn%n_cohorts+newcohorts))
  vegn%cohorts(1:vegn%n_cohorts) = ccold(1:vegn%n_cohorts) ! copy old cohort information
  deallocate (ccold)

!   set up new cohorts
  k = vegn%n_cohorts
  do i = 1,vegn%n_cohorts
    if (.not.cohort_can_reproduce(vegn%cohorts(i))) cycle ! nothing to do

    k = k+1 ! new cohort index
    
    vegn%cohorts(k) = vegn%cohorts(i) ! copy all values into the new cohort then change some of them below
    cc => vegn%cohorts(k) ; parent => vegn%cohorts(i)

    ! update child cohort parameters
    cc%nindivs = 1.0 * parent%nindivs  ! It should be dependent on the fecundity of species
    cc%age     = 0.2 * parent%age      ! 0.
    cc%bliving = 0.2 * parent%bliving  ! 0.02  ! these values need to be tested that if they are
    cc%bwood   = 0.2 * parent%bwood    ! 0.01  ! OK to initialize plant growth

    ! partition living pools
    call update_biomass_pools(cc)

    cc%gpp          = 0.
    cc%npp          = 0.  
    cc%resp         = 0.
    cc%resl         = 0.
    cc%resr         = 0.
    cc%resg         = 0.
    cc%carbon_gain  = 0.
    cc%carbon_loss  = 0.
    cc%bwood_gain   = 0.
    cc%npp_previous_day     = 0.
    cc%npp_previous_day_tmp = 0.
    cc%leaf_age = 0.0
    ! we assume that the newborn cohort is dry; since nindivs of the parent
    ! doesn't change we don't need to do anything with its Wl and Ws to 
    ! conserve water (since Wl and Ws are per individual)
    cc%prog%Wl = 0 ; cc%prog%Ws = 0
    ! TODO: make sure that energy is conserved in reproduction
    cc%prog%Tv = parent%prog%Tv

    ! update C of parent cohort
    parent%bliving = parent%bliving - cc%bliving * cc%nindivs/parent%nindivs
    parent%bwood   = parent%bwood   - cc%bwood   * cc%nindivs/parent%nindivs
    call update_biomass_pools(parent)
  enddo
  
  vegn%n_cohorts = k
  
end subroutine vegn_reproduction_ppa


function cohorts_can_be_merged(c1,c2); logical cohorts_can_be_merged
   type(vegn_cohort_type), intent(in) :: c1,c2
   cohorts_can_be_merged = (c1%layer == c2%layer) &
                     .and. (c1%species == c2%species)
end function 

! ============================================================================
subroutine merge_cohorts(c1,c2)
  type(vegn_cohort_type), intent(in) :: c1
  type(vegn_cohort_type), intent(inout) :: c2
  
  real :: x1, x2 ! normalized relative weights
  real :: HEAT1, HEAT2 ! heat stored in respective canopies

  x1 = c1%nindivs/(c1%nindivs+c2%nindivs)
  x2 = 1-x1
  ! update number of individuals in merged cohort
  c2%nindivs = c1%nindivs+c2%nindivs
#define __MERGE__(field) c2%field = x1*c1%field + x2*c2%field
  HEAT1 = (clw*c1%prog%Wl + csw*c1%prog%Ws + c1%mcv_dry)*(c1%prog%Tv-tfreeze)
  HEAT2 = (clw*c2%prog%Wl + csw*c2%prog%Ws + c2%mcv_dry)*(c2%prog%Tv-tfreeze)
  __MERGE__(prog%Wl)
  __MERGE__(prog%Ws)
 
  __MERGE__(bl)      ! biomass of leaves, kg C/indiv
  __MERGE__(blv)     ! biomass of virtual leaves (labile store), kg C/indiv
  __MERGE__(br)      ! biomass of fine roots, kg C/indiv
  __MERGE__(bsw)     ! biomass of sapwood, kg C/indiv
  __MERGE__(bwood)   ! biomass of heartwood, kg C/indiv
  __MERGE__(bliving) ! leaves, fine roots, and sapwood biomass, kgC/indiv
  __MERGE__(nsc)     ! non-structural carbon, kgC/indiv
  __MERGE__(bseed)   ! future progeny, kgC/indiv
  __MERGE__(age)     ! age of individual

  __MERGE__(carbon_gain) ! carbon gain during a day, kg C/indiv
  __MERGE__(carbon_loss) ! carbon loss during a day, kg C/indiv [diag only]
  __MERGE__(bwood_gain)  ! heartwood gain during a day, kg C/indiv

  call update_biomass_pools(c2)

  ! calculate the resulting dry heat capacity
  c2%leafarea = leaf_area_from_biomass(c2%bl, c2%species)
  c2%mcv_dry = max(mcv_min,mcv_lai*c2%leafarea)
  ! update canopy temperature -- just merge it based on area weights if the heat 
  ! capacities are zero, or merge it based on the heat content if the heat contents
  ! are non-zero
  if(HEAT1==0.and.HEAT2==0) then
     __MERGE__(prog%Tv)
  else
     c2%prog%Tv = (HEAT1*x1+HEAT2*x2) / &
          (clw*c2%prog%Wl + csw*c2%prog%Ws + c2%mcv_dry) + tfreeze
  endif
#undef  __MERGE__
  
end subroutine merge_cohorts

! ============================================================================
! Merge similar cohorts in a tile
subroutine vegn_mergecohorts_ppa(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

! ---- local vars
  type(vegn_cohort_type), pointer :: cc(:) ! array to hold new cohohrts
  logical :: merged(vegn%n_cohorts)        ! mask to skip cohorts that were already merged
  integer :: i,j,k

  allocate(cc(vegn%n_cohorts))

  merged(:)=.FALSE. ; k = 0
  do i = 1, vegn%n_cohorts 
     if(merged(i)) cycle ! skip cohorts that were alraedy merged
     k = k+1
     cc(k) = vegn%cohorts(i)
     ! try merging therest of the cohorts into current one
     do j = i+1, vegn%n_cohorts
        if (merged(j)) cycle ! skip cohorts that are already merged
        if (cohorts_can_be_merged(cc(k),vegn%cohorts(j))) then
           call merge_cohorts(vegn%cohorts(j),cc(k))
           merged(j) = .TRUE.
        endif
     enddo
  enddo
  ! at this point, k is the number of new cohorts
  vegn%n_cohorts = k
  deallocate(vegn%cohorts)
  vegn%cohorts=>cc

  ! note that the size of the vegn%cohorts may be larger than vegn%n_cohorts
  
end subroutine vegn_mergecohorts_ppa


! =============================================================================
! given an array of cohorts, create a new array with old cohorts re-arranged
! in layers according to their height and crown areas.
subroutine relayer_cohorts (vegn)
  type(vegn_tile_type), intent(inout) :: vegn ! input cohorts

  ! ---- local constants
  real, parameter :: tolerance = 1e-6 
  
  ! ---- local vars
  integer :: idx(vegn%n_cohorts) ! indices of cohorts in decreasing height order
  integer :: i ! new cohort index
  integer :: k ! old cohort index
  integer :: L ! layer index (top-down)
  integer :: N0,N1 ! initial and final number of cohorts 
  real    :: frac ! fraction of the layer covered so far by the canopies
  type(vegn_cohort_type), pointer :: cc(:),new(:)
  real    :: nindivs
  
  ! rank cohorts in descending order by height. For now, assume that they are 
  ! in order
  N0 = vegn%n_cohorts; cc=>vegn%cohorts
  call rank_descending(cc(1:N0)%height,idx)
  
  ! calculate max possible number of new cohorts : it is equal to the number of
  ! old cohorts, plus the number of layers -- since the number of full layers is 
  ! equal to the maximum number of times an input cohort can be split by a layer 
  ! boundary.
  N1 = vegn%n_cohorts + int(sum(cc(1:N0)%nindivs*cc(1:N0)%crownarea))
  allocate(new(N1))

  ! copy cohort information to the new cohorts, splitting the old cohorts that 
  ! stride the layer boundaries
  i = 1 ; k = 1 ; L = 1 ; frac = 0.0 ; nindivs = cc(idx(k))%nindivs
  do 
     new(i)         = cc(idx(k))
     new(i)%nindivs = min(nindivs,(1-frac)/cc(idx(k))%crownarea)
     new(i)%layer   = L
     frac = frac+new(i)%nindivs*new(i)%crownarea
     nindivs = nindivs - new(i)%nindivs
     
     if (abs(nindivs*cc(idx(k))%crownarea)<tolerance) then
       new(i)%nindivs = new(i)%nindivs + nindivs ! allocate the remainder of individuals to the last cohort
       if (k==N0) exit ! end of loop
       k = k+1 ; nindivs = cc(idx(k))%nindivs  ! go to the next input cohort
     endif
     
     if (abs(1-frac)<tolerance) then
       L = L+1 ; frac = 0.0              ! start new layer
     endif

     i = i+1
  enddo
  
  ! replace the array of cohorts
  deallocate(vegn%cohorts)
  vegn%cohorts => new ; vegn%n_cohorts = i
end subroutine relayer_cohorts

end module vegn_dynamics_mod
