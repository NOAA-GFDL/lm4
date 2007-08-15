! ============================================================================
! updates carbon pools and rates on the fast time scale
! ============================================================================
module vegn_carbon_int_mod

use fms_mod, only: write_version_number
use time_manager_mod, only: time_type

use land_constants_mod, only : seconds_per_year, mol_C
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type
use vegn_data_mod, only : spdata, &
     CMPT_VLEAF, CMPT_SAPWOOD, CMPT_ROOT, CMPT_WOOD, CMPT_LEAF, LEAF_ON, LEAF_OFF, &
     fsc_liv, fsc_wood, K1, K2
use vegn_tile_mod, only: vegn_tile_type
use cohort_list_mod, only : vegn_cohort_enum_type, first_cohort, tail_cohort, &
     current_cohort, next_cohort, operator(/=)
use vegn_cohort_mod, only : vegn_cohort_type, height_from_biomass, lai_from_biomass, &
     update_bio_living_fraction

implicit none
public

! ==== public interfaces =====================================================
public :: vegn_carbon_int_init
public :: vegn_carbon_int_end

public :: carbon_int ! fast time-scale integrator of carbon balance
public :: growth     ! slow time-scale redistributor of accumulated carbon
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: vegn_carbon_int.F90,v 15.0 2007/08/14 03:59:51 fms Exp $', &
   tagname = '$Name: omsk $' ,&
   module_name = 'vegn'
real, parameter :: GROWTH_RESP=0.333  ! fraction of npp lost as growth respiration


! ==== module data ===========================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

! diagnostic field IDs
integer :: id_npp, id_nep, id_rsoil, id_rsoil_fast


contains

! ============================================================================
subroutine vegn_carbon_int_init(id_lon, id_lat, time, delta_time)
  integer        , intent(in) :: id_lon ! ID of land longitude (X) axis 
  integer        , intent(in) :: id_lat ! ID of land latitude (Y) axis
  type(time_type), intent(in) :: time       ! initial time for diagnostic fields
  real           , intent(in) :: delta_time ! fast time step, s

  call write_version_number(version, tagname)

  ! set up global variables
  dt_fast_yr = delta_time/seconds_per_year

  ! register diagnostic fields
  id_npp = register_tiled_diag_field ( module_name, 'npp',  &
       (/id_lon,id_lat/), time, 'net primary productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_nep = register_tiled_diag_field ( module_name, 'nep',  &
       (/id_lon,id_lat/), time, 'net ecosystem productivity', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil = register_tiled_diag_field ( module_name, 'rsoil',  &
       (/id_lon,id_lat/), time, 'soil respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_fast = register_tiled_diag_field ( module_name, 'rsoil_fast',  &
       (/id_lon,id_lat/), time, 'fast soil carbon respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
end subroutine vegn_carbon_int_init


! ============================================================================
subroutine vegn_carbon_int_end
end subroutine vegn_carbon_int_end


! ============================================================================
subroutine carbon_int(vegn, diag, soilt, theta)
  type(vegn_tile_type), intent(inout) :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! temperature of soil for soil carbon decomposition
  real                , intent(in)    :: theta

  type(vegn_cohort_type), pointer :: cc
  type(vegn_cohort_enum_type) :: ci,ce ! cohort enumerators

  real :: md_alive, md_wood;
  integer :: sp ! copy of current cohort species, for convinience

  !  update plant carbon
  ci = first_cohort(vegn%cohorts);
  ce = tail_cohort (vegn%cohorts);
  do while(ci/=ce)   
     cc => current_cohort(ci)
     ci = next_cohort(ci) ! advance the enumerator for the next step

     sp = cc%species

     call eddy_npp(cc,soilt);
     ! npp2 is for diagnostics and comparison
     cc%npp2 = cc%miami_npp;  ! treat miami npp as above+below npp
     
     cc%carbon_gain = cc%carbon_gain + cc%npp*dt_fast_yr;
     
     ! check if leaves/roots are present and need to be accounted in maintanence
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
        
     cc%md = md_alive + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain - md_wood;
     if (cc%bwood_gain < 0.0) cc%bwood_gain=0.0; ! potential non-conservation ?
     cc%carbon_gain = cc%carbon_gain - cc%md;
     cc%carbon_loss = cc%carbon_loss + cc%md; ! used in diagnostics only

     ! add md from leaf and root pools to fast soil carbon
     vegn%fast_soil_C = vegn%fast_soil_C +    fsc_liv *md_alive +    fsc_wood *md_wood;
     vegn%slow_soil_C = vegn%slow_soil_C + (1-fsc_liv)*md_alive + (1-fsc_wood)*md_wood;

     ! for budget tracking
!/*     cp->fsc_in+= data->fsc_liv*md_alive+data->fsc_wood*md_wood; */
!/*     cp->ssc_in+= (1.- data->fsc_liv)*md_alive+(1-data->fsc_wood)*md_wood; */
     vegn%fsc_in  = vegn%fsc_in + 1*md_alive+0*md_wood;
     vegn%ssc_in  = vegn%ssc_in + (1.- 1)*md_alive+(1-0)*md_wood;

     vegn%veg_in  = vegn%veg_in  + cc%npp*dt_fast_yr;
     vegn%veg_out = vegn%veg_out + md_alive+md_wood;

  enddo

  ! update soil carbon
  call Dsdt(vegn, diag, soilt, theta)

end subroutine carbon_int


! ============================================================================
! updates cohort biomass pools, LAI, SAI, and height using accumulated 
! carbon_gain and bwood_gain
subroutine growth (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  type(vegn_cohort_enum_type)     :: ci,ce ! cohort enumerators

  ci = first_cohort(vegn%cohorts);
  ce = tail_cohort (vegn%cohorts);
  do while(ci/=ce)   
     cc => current_cohort(ci)
     ci = next_cohort(ci) ! advance the enumerator for the next step

     cc%bwood   = cc%bwood   + cc%bwood_gain
     cc%bliving = cc%bliving + cc%carbon_gain
     
     if(cc%bliving < 0) then
        cc%bwood    = cc%bwood+cc%bliving
        cc%bliving  = 0
        if (cc%bwood < 0) &
             cc%bwood = 0 ! in principle, that's not conserving carbon
     endif
     
     cc%b           = cc%bliving + cc%bwood
     cc%height = height_from_biomass(cc%bliving+cc%bwood)
     ! update biomass compartment fractions
     call update_bio_living_fraction(cc)
     ! update biomass compartments
     cc%bsw = cc%Psw*cc%bliving
     if(cc%status == LEAF_OFF) then
        ! leaves are off, no fine roots, comes out of bliving
        cc%bl = 0
        cc%br = 0
        cc%blv = cc%Pl*cc%bliving + cc%Pr*cc%bliving
     else
        ! leaves and fine roots are on
        cc%bl = cc%Pl*cc%bliving
        cc%br = cc%Pr*cc%bliving
        cc%blv = 0
     endif
     
     cc%lai = lai_from_biomass(cc%bl, cc%species)
     cc%sai = 0.035 * cc%height ! Federer and Lash,1978
     
     ! reset carbon acculmulation terms
     cc%carbon_gain = 0
     cc%carbon_loss = 0
     cc%bwood_gain  = 0
  end do
end subroutine growth


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
  call send_tile_data(id_rsoil_fast, fast_C_loss/dt_fast_yr, diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)

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
subroutine eddy_npp(cc, tsoil)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: tsoil

  call plant_respiration(cc,tsoil);

  cc%gpp = (cc%An_op - cc%An_cl)*mol_C*cc%lai;
  cc%npp = cc%gpp - cc%resp;

  if(cc%npp_previous_day > -0.00001/2500.0) then
     cc%resg = GROWTH_RESP*cc%npp_previous_day;
     cc%npp  = cc%npp  - GROWTH_RESP*cc%npp_previous_day;
     cc%resp = cc%resp + GROWTH_RESP*cc%npp_previous_day;
  else
     cc%resg = 0;
  endif

  ! update, accumulate
  cc%npp_previous_day_tmp = cc%npp_previous_day_tmp + cc%npp; 
end subroutine eddy_npp


! ============================================================================
subroutine plant_respiration(cc, tsoil)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: tsoil
  
  real :: tf,tfs;
  real :: r_leaf, r_vleaf, r_stem, r_root
  
  integer :: sp ! copy of cohort species
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

  r_leaf = -mol_C*cc%An_cl*cc%lai;
  r_vleaf = spdata(sp)%beta(CMPT_VLEAF)   * cc%blv*tf;
  r_stem  = spdata(sp)%beta(CMPT_SAPWOOD) * cc%bsw*tf;
  r_root  = spdata(sp)%beta(CMPT_ROOT)    * cc%br*tfs;
  
  cc%resp = r_leaf + r_vleaf + r_stem + r_root;
  cc%resl = r_leaf;
  cc%resr = r_root;
end subroutine plant_respiration


! ============================================================================
! calculates prev. day average NPP from accumualted values
subroutine daily_npp(vegn, dt_fast_yr)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: dt_fast_yr

  integer :: n_fast_step;
  type(vegn_cohort_type), pointer :: cc
  type(vegn_cohort_enum_type) :: ci,ce ! cohort enumerators

  n_fast_step = 86400/dt_fast_yr;
  ci = first_cohort(vegn%cohorts);
  ce = tail_cohort (vegn%cohorts);
  do while(ci/=ce)   
     cc => current_cohort(ci)
     ci = next_cohort(ci) ! advance the enumerator for the next step

     cc%npp_previous_day=cc%npp_previous_day_tmp/n_fast_step;
     cc%npp_previous_day_tmp=0.0
  enddo
end subroutine daily_npp


! =============================================================================
! The stuff below comes from she_update.c -- it looks like it belongs here, 
! since it is essentially a part of the carbon integration (update_patch_fast
! is only called immediately after carbon_int in lm3v)
! =============================================================================


! =============================================================================
subroutine update_patch_fast(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  call update_biomass(vegn);
!!$  call update_hite(vegn);
  call update_npp(vegn);
  call update_soil_pools(vegn);
  
  vegn%age = vegn%age + dt_fast_yr;
end subroutine update_patch_fast


! =============================================================================
subroutine update_biomass(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc
  type(vegn_cohort_enum_type) :: ci,ce ! cohort enumerators

  ci = first_cohort(vegn%cohorts);
  ce = tail_cohort (vegn%cohorts);
  vegn%total_biomass = 0
  do while(ci/=ce)   
     cc => current_cohort(ci)
     ci = next_cohort(ci) ! advance the enumerator for the next step

     vegn%total_biomass = vegn%total_biomass + cc%b
  end do
end subroutine update_biomass


! =============================================================================
!!$subroutine update_hite(vegn)
!!$  type(vegn_tile_type), intent(inout) :: vegn
!!$
!!$  ! ---- local vars
!!$  type(vegn_cohort_type), pointer :: cc
!!$  type(vegn_cohort_enum_type) :: ci,ce ! cohort enumerators
!!$
!!$  ci = first_cohort(vegn%cohorts);
!!$  ce = tail_cohort (vegn%cohorts);
!!$  vegn%height = 0
!!$  do while(ci/=ce)   
!!$     cc => current_cohort(ci)
!!$     ci = next_cohort(ci) ! advance the enumerator for the next step
!!$
!!$     cc%height = height_from_biomass(cc%b)
!!$     vegn%height = max(cc%height,vegn%height)
!!$  end do
!!$end subroutine update_hite
  

! =============================================================================
subroutine update_npp(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc
  type(vegn_cohort_enum_type) :: ci,ce ! cohort enumerators

  vegn%npp = 0.0;
  vegn%nep = 0.0;
  ! vegn%npp2 = 0.0;

  ci = first_cohort(vegn%cohorts);
  ce = tail_cohort (vegn%cohorts);
  do while(ci/=ce)   
     if (cc%height > 0) then
        vegn%npp = vegn%npp + cc%npp
        vegn%nep = vegn%npp + cc%npp
        ! vegn%npp2 = vegn%npp2 + cc%npp2
     end if
  end do
  vegn%nep = vegn%nep - vegn%rh;
end subroutine update_npp
  

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


end module vegn_carbon_int_mod
