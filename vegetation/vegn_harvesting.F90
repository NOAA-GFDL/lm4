module vegn_harvesting_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

#include "../shared/debug.inc"

use fms_mod, only : string, error_mesg, FATAL, NOTE, &
     mpp_pe, file_exist, close_file, &
     check_nml_error, stdlog, mpp_root_pe
use mpp_io_mod, only : axistype, mpp_get_atts, mpp_get_axis_data, &
     mpp_open, mpp_close, MPP_RDONLY, MPP_WRONLY, MPP_ASCII
use land_data_mod, only : log_version
use vegn_data_mod, only : &
     N_LU_TYPES, LU_PAST, LU_CROP, LU_NTRL, LU_SCND, &
     HARV_POOL_PAST, HARV_POOL_CROP, HARV_POOL_CLEARED, HARV_POOL_WOOD_FAST, &
     HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW, &
     spdata, agf_bs
use soil_tile_mod, only : LEAF, CWOOD, soil_tile_type
use soil_mod, only : add_root_litter
use vegn_tile_mod, only : vegn_tile_type, vegn_tile_lai
use vegn_cohort_mod, only : vegn_cohort_type, update_biomass_pools
use soil_carbon_mod, only: soil_carbon_option, &
     SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, SOILC_CORPSE_N, N_C_TYPES, add_litter

implicit none
private

! ==== public interface ======================================================
public :: vegn_harvesting_init
public :: vegn_harvesting_end

public :: vegn_harvesting

public :: vegn_graze_pasture
public :: vegn_harvest_cropland
public :: vegn_cut_forest

public :: do_harvesting
! ==== end of public interface ===============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_harvesting_mod'
#include "../shared/version_variable.inc"
real, parameter :: ONETHIRD = 1.0/3.0
integer, parameter :: DAILY = 1, ANNUAL = 2

! ==== module data ===========================================================

! ---- namelist variables ----------------------------------------------------
logical :: do_harvesting       = .TRUE.  ! if true, then harvesting of crops and pastures is done
real :: grazing_intensity      = 0.25    ! fraction of biomass removed each time by grazing
real :: grazing_residue        = 0.1     ! fraction of the grazed biomass transferred into soil pools
real :: frac_wood_wasted_harv  = 0.25    ! fraction of wood wasted while harvesting
real :: frac_wood_wasted_clear = 0.25    ! fraction of wood wasted while clearing land for pastures or crops
logical :: waste_below_ground_wood = .TRUE. ! If true, all the wood below ground (1-agf_bs fraction of bwood
        ! and bsw) is wasted. Old behavior assumed this to be FALSE.
real :: frac_wood_fast         = ONETHIRD ! fraction of wood consumed fast
real :: frac_wood_med          = ONETHIRD ! fraction of wood consumed with medium speed
real :: frac_wood_slow         = ONETHIRD ! fraction of wood consumed slowly
real :: crop_seed_density      = 0.1     ! biomass of seeds left after crop harvesting, kg/m2
real :: min_lai_for_grazing    = 0.0     ! no grazing if LAI lower than this threshold
character(16) :: grazing_frequency = 'annual' ! or 'daily'
namelist/harvesting_nml/ do_harvesting, grazing_intensity, grazing_residue, &
     frac_wood_wasted_harv, frac_wood_wasted_clear, waste_below_ground_wood, &
     frac_wood_fast, frac_wood_med, frac_wood_slow, &
     crop_seed_density, &
     grazing_frequency, min_lai_for_grazing

integer :: grazing_freq = -1 ! inidicator of grazing frequency (ANNUAL or DAILY)

contains ! ###################################################################

! ============================================================================
subroutine vegn_harvesting_init
  integer :: unit, ierr, io

  call log_version(version, module_name, &
  __FILE__)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=harvesting_nml, iostat=io)
  ierr = check_nml_error(io, 'harvesting_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=harvesting_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'harvesting_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=harvesting_nml)
  endif

  if (frac_wood_fast+frac_wood_med+frac_wood_slow/=1.0) then
     call error_mesg('vegn_harvesting_init', &
          'sum of frac_wood_fast, frac_wood_med, and frac_wood_slow must be 1.0',&
          FATAL)
  endif
  ! parse the grazing frequency parameter
  select case(grazing_frequency)
  case('annual')
     grazing_freq = ANNUAL
  case('daily')
     grazing_freq = DAILY
     ! scale grazing intensity for daily frequency
     grazing_intensity = grazing_intensity/365.0
  case default
     call error_mesg('vegn_harvesting_init','grazing_frequency must be "annual" or "daily"',FATAL)
  end select
end subroutine vegn_harvesting_init


! ============================================================================
subroutine vegn_harvesting_end
end subroutine vegn_harvesting_end


! ============================================================================
! harvest vegetation in a tile
subroutine vegn_harvesting(vegn,soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  if (.not.do_harvesting) &
       return ! do nothing if no harvesting requested

  select case(vegn%landuse)
  case(LU_PAST)  ! pasture
     call vegn_graze_pasture    (vegn,soil)
  case(LU_CROP)  ! crop
     call vegn_harvest_cropland (vegn)
  end select
end subroutine


! ============================================================================
subroutine vegn_graze_pasture(vegn,soil)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil

  ! ---- local vars
  real ::  bdead0, balive0, bleaf0, blv0, bfroot0 ! initial combined biomass pools
  real ::  bdead1, balive1, bleaf1, blv1, bfroot1 ! updated combined biomass pools
  integer :: i
  real :: carbon_lost
  real :: delta_leaf, delta_root, delta_wood
  real,dimension(n_C_types) :: leaflitter_C,woodlitter_C,bglitter_C,leaflitter_N,woodlitter_N,bglitter_N
  real :: wood_n2c

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1,vegn%n_cohorts
     associate (cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))

  if ( vegn_tile_LAI(vegn) .lt. min_lai_for_grazing ) return

     ! This makes sure biomass pools are correct before calculating changes
     ! in leaf biomass and such, just in case it wasn't called before
     call update_biomass_pools(cc);

         if(cc%bwood>0) then
           wood_n2c=cc%wood_N/cc%bwood
         else
           wood_n2c=0.0
         endif
     ! calculate total biomass pools for the patch
     balive0 =  cc%bl + cc%blv + cc%br
     bleaf0  =  cc%bliving*cc%Pl  !cc%bl + cc%blv
     bfroot0 =  cc%bliving*cc%Pr  !cc%br
     if(cc%bl+cc%br>0) then  ! Not leaf off/deciduous winter
       blv0 = cc%blv
     else
       blv0    =  cc%blv - cc%bliving*(cc%Pl+cc%Pr) ! Excess carbon (due to N limitation)
     endif
     bdead0  =  cc%bwood + cc%bsw
     ! only potential leaves are consumed
     carbon_lost=cc%bliving*cc%Pl*grazing_intensity
     vegn%harv_pool(HARV_POOL_PAST) = vegn%harv_pool(HARV_POOL_PAST) + &
          carbon_lost*(1-grazing_residue) ;
     cc%bliving = cc%bliving - carbon_lost;
     cc%stored_N = cc%stored_N - carbon_lost/sp%leaf_live_c2n

     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc);

     ! calculate new combined vegetation biomass pools
     balive1 =  cc%bl + cc%blv + cc%br
     bleaf1  =  cc%bliving*cc%Pl  !cc%bl + cc%blv
     bfroot1 =  cc%bliving*cc%Pr  !cc%br
     if(cc%bl+cc%br>0) then
       blv1 = cc%blv
     else
       blv1    =  cc%blv - cc%bliving*(cc%Pl+cc%Pr) ! Excess carbon (due to N limitation)
     endif
     bdead1  =  cc%bwood + cc%bsw


     ! update intermediate soil carbon pools
     select case(soil_carbon_option)
     case(SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + grazing_residue*( &
             sp%fsc_liv*(balive0-balive1)+sp%fsc_wood*(bdead0-bdead1))
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + grazing_residue*( &
             (1-sp%fsc_liv)*(balive0-balive1)+ (1-sp%fsc_wood)*(bdead0-bdead1))
     case(SOILC_CORPSE, SOILC_CORPSE_N)
        if(blv0 < blv1) then ! Some biomass was re-absorbed due to N limitation. Reduce litter.
          delta_leaf=bleaf0-bleaf1   - (blv1-blv0)*(bleaf0-bleaf1)/(bleaf0+bfroot0+bdead0-bleaf1-bfroot1-bdead1)
          delta_root=bfroot0-bfroot1 - (blv1-blv0)*(bfroot0-bfroot1)/(bleaf0+bfroot0+bdead0-bleaf1-bfroot1-bdead1)
          delta_wood=bdead0-bdead1   - (blv1-blv0)*(bdead0-bdead1)/(bleaf0+bfroot0+bdead0-bleaf1-bfroot1-bdead1)
        else   ! Virtual leaves decreased. Send to leaf litter
          delta_leaf=bleaf0+blv0-bleaf1-blv1
          delta_root=bfroot0-bfroot1
          delta_wood=bdead0-bdead1
        endif


       leaflitter_C=(/(delta_leaf)*sp%fsc_liv,(delta_leaf)*(1-sp%fsc_liv),0.0/)*grazing_residue
       woodlitter_C=(/(delta_wood)*sp%fsc_wood,(delta_wood)*(1-sp%fsc_wood),0.0/)*agf_bs*grazing_residue
       bglitter_C=(/(sp%fsc_froot*(delta_root) +(1-agf_bs)*sp%fsc_wood*(delta_wood)),&
                     (1.0-sp%fsc_froot)*(delta_root) +(1-agf_bs)*(1.0-sp%fsc_wood)*(delta_wood),0.0/)*grazing_residue

       if(soil_carbon_option == SOILC_CORPSE_N) then
         leaflitter_N=leaflitter_C/sp%leaf_live_c2n
         woodlitter_N=woodlitter_C/sp%leaf_live_c2n
         bglitter_N=(/grazing_residue*(sp%fsc_froot*(delta_root)/sp%froot_live_c2n +(1-agf_bs)*sp%fsc_wood*(delta_wood)*wood_n2c),&
                      grazing_residue*((1-sp%fsc_froot)*(delta_root)/sp%froot_live_c2n +  (1-agf_bs)*(1-sp%fsc_wood)*(delta_wood)*wood_n2c),&
                      0.0/)
       else
         leaflitter_N=(/0.0,0.0,0.0/)
         woodlitter_N=(/0.0,0.0,0.0/)
         bglitter_N=(/0.0,0.0,0.0/)
       endif

       if (grazing_freq==DAILY) then
         ! Put carbon directly in soil pools
         call add_litter(soil%litter(LEAF),leaflitter_C,leaflitter_N)
         call add_litter(soil%litter(CWOOD),woodlitter_C,woodlitter_N)
         call add_root_litter(soil,vegn,bglitter_C,bglitter_N)
       else
         vegn%litter_buff_C(:,LEAF) = vegn%litter_buff_C(:,LEAF) + &
              [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*(delta_leaf)*grazing_residue
         vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + &
              [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*agf_bs*(delta_wood)*grazing_residue

         vegn%fsc_pool_bg=vegn%fsc_pool_bg + bglitter_C(1)
         vegn%ssc_pool_bg = vegn%ssc_pool_bg + bglitter_C(2)


         vegn%litter_buff_N(:,LEAF) = vegn%litter_buff_N(:,LEAF) + &
            [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*(delta_leaf)*grazing_residue/sp%leaf_live_c2n
         vegn%litter_buff_N(:,CWOOD) = vegn%litter_buff_N(:,CWOOD) + &
            [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*agf_bs*(delta_wood)*grazing_residue/sp%wood_c2n


         vegn%fsn_pool_bg=vegn%fsn_pool_bg + bglitter_N(1)
         vegn%ssn_pool_bg = vegn%ssn_pool_bg + bglitter_N(2)

       endif
     case default
        call error_mesg('vegn_graze_pasture','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select
     end associate
  enddo
end subroutine vegn_graze_pasture


! ================================================================================
subroutine vegn_harvest_cropland(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real :: fraction_harvested;    ! fraction of biomass harvested this time
  real :: bdead, balive, btotal; ! combined biomass pools
  integer :: i

  balive = 0 ; bdead = 0
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i))
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
     end associate
  enddo
  btotal = balive+bdead;

  ! calculate harvested fraction: cut everything down to seed level
  fraction_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0);

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))
     ! use for harvest only aboveg round living biomass and waste the correspondent below living and wood
     vegn%harv_pool(HARV_POOL_CROP) = vegn%harv_pool(HARV_POOL_CROP) + &
          cc%bliving*(cc%Pl + cc%Psw*agf_bs)*fraction_harvested
     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + fraction_harvested*(sp%fsc_liv*cc%bliving*cc%Pr + &
             sp%fsc_wood*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)))
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + fraction_harvested*((1-sp%fsc_liv)*cc%bliving*cc%Pr + &
             (1-sp%fsc_wood)*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)))
     case (SOILC_CORPSE, SOILC_CORPSE_N)
        vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + &
               [sp%fsc_wood, 1-sp%fsc_wood, 0.0] * fraction_harvested*agf_bs*cc%bwood

        vegn%fsc_pool_bg = vegn%fsc_pool_bg + fraction_harvested*(&
               sp%fsc_froot*cc%bliving*cc%Pr + &
               (1-agf_bs)*sp%fsc_wood*(cc%bwood + cc%bliving*cc%Psw))
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + fraction_harvested*(&
               (1-sp%fsc_froot)*cc%bliving*cc%Pr + &
               (1-agf_bs)*(1-sp%fsc_wood)*(cc%bwood + cc%bliving*cc%Psw))

        if (soil_carbon_option == SOILC_CORPSE_N) then
           vegn%litter_buff_N(:,CWOOD) = vegn%litter_buff_N(:,CWOOD) + &
               [sp%fsc_wood, 1-sp%fsc_wood, 0.0] * fraction_harvested*agf_bs*cc%wood_N

           vegn%fsn_pool_bg = vegn%fsn_pool_bg + fraction_harvested*(&
                   sp%fsc_froot*cc%root_N + &
                   (1-agf_bs)*(sp%fsc_wood*cc%wood_N + sp%fsc_liv*cc%sapwood_N))
           vegn%ssn_pool_bg = vegn%ssn_pool_bg + fraction_harvested*(&
                   (1-sp%fsc_froot)*cc%root_N + &
                   (1-agf_bs)*(cc%wood_N*(1-sp%fsc_wood) + cc%sapwood_N*(1-sp%fsc_liv)))
        endif
     case default
        call error_mesg('vegn_harvest_cropland','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select

    !  cc%stored_N = cc%stored_N - cc%bliving*cc%Pl*(1-fraction_harvested)/spdata(sp)%leaf_live_c2n
    !  cc%stored_N = cc%stored_N - cc%bliving*cc%Pr*(1-fraction_harvested)/spdata(sp)%froot_live_c2n
      cc%stored_N = cc%stored_N*(1-fraction_harvested)
      cc%leaf_N=cc%leaf_N*(1-fraction_harvested)
      cc%root_N=cc%root_N*(1-fraction_harvested)
     cc%wood_N = cc%wood_N*(1-fraction_harvested)
     cc%sapwood_N = cc%sapwood_N*(1-fraction_harvested)
     ! Should stored N be lost or retained?
     ! redistribute leftover biomass between biomass pools

     cc%bliving = cc%bliving * (1-fraction_harvested);
     cc%bwood   = cc%bwood   * (1-fraction_harvested);

     call update_biomass_pools(cc);
     end associate
  enddo
end subroutine vegn_harvest_cropland


! ============================================================================
! for now cutting forest is the same as harvesting cropland --
! we basically cut down everything, leaving only seeds
subroutine vegn_cut_forest(vegn, new_landuse)
  type(vegn_tile_type), intent(inout) :: vegn
  integer, intent(in) :: new_landuse ! new land use type that gets assigned to
                                     ! the tile after the wood harvesting

  ! ---- local vars
  real :: frac_harvested;        ! fraction of biomass harvested this time
  real :: frac_wood_wasted       ! fraction of wood wasted during transition
  real :: frac_wood_wasted_ag    ! fraction of above-ground wood wasted during transition
  real :: wood_harvested         ! anount of harvested wood, kgC/m2
  real :: bdead, balive, bleaf, bfroot, btotal; ! combined biomass pools
  real :: delta
  integer :: i

  balive = 0 ; bdead = 0 ; bleaf = 0 ; bfroot = 0 ;
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     associate (cc=>vegn%cohorts(i))
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bleaf  = bleaf  + cc%bl + cc%blv
     bfroot = bfroot + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
     end associate
  enddo
  btotal = balive+bdead;

  ! calculate harvested fraction: cut everything down to seed level
  frac_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0);

  ! define fraction of wood wasted, based on the transition type
  if (new_landuse==LU_SCND) then
     frac_wood_wasted = frac_wood_wasted_harv
  else
     frac_wood_wasted = frac_wood_wasted_clear
  endif
  ! take into accont that all wood below ground is wasted; also the fraction
  ! of waste calculated above is lost from the above-ground part of the wood
  frac_wood_wasted_ag=frac_wood_wasted
  if (waste_below_ground_wood) then
     frac_wood_wasted = (1-agf_bs) + agf_bs*frac_wood_wasted
  endif

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     associate(cc => vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))

     ! calculate total amount of harvested wood, minus the wasted part
     wood_harvested = (cc%bwood+cc%bsw)*frac_harvested*(1-frac_wood_wasted)

     ! distribute harvested wood between pools
     if (new_landuse==LU_SCND) then
        ! this is harvesting, distribute between 3 different wood pools
        vegn%harv_pool(HARV_POOL_WOOD_FAST) = vegn%harv_pool(HARV_POOL_WOOD_FAST) &
             + wood_harvested*frac_wood_fast
        vegn%harv_pool(HARV_POOL_WOOD_MED) = vegn%harv_pool(HARV_POOL_WOOD_MED) &
             + wood_harvested*frac_wood_med
        vegn%harv_pool(HARV_POOL_WOOD_SLOW) = vegn%harv_pool(HARV_POOL_WOOD_SLOW) &
             + wood_harvested*frac_wood_slow
     else
        ! this is land clearance: everything goes into "cleared" pool
        vegn%harv_pool(HARV_POOL_CLEARED) = vegn%harv_pool(HARV_POOL_CLEARED) &
             + wood_harvested
     endif

     ! distribute wood and living biomass between fast and slow intermediate
     ! soil carbon pools according to fractions specified through the namelists
     delta = (cc%bwood+cc%bsw)*frac_harvested*frac_wood_wasted;
     if(delta<0) call error_mesg('vegn_cut_forest', &
          'harvested amount of dead biomass ('//string(delta)//' kgC/m2) is below zero', &
          FATAL)

     select case (soil_carbon_option)
     case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + delta*(1-sp%fsc_wood)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + delta*   sp%fsc_wood

        delta = balive * frac_harvested;
        if(delta<0) call error_mesg('vegn_cut_forest', &
             'harvested amount of live biomass ('//string(delta)//' kgC/m2) is below zero', &
             FATAL)
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + delta*(1-sp%fsc_liv) ;
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + delta*   sp%fsc_liv  ;
     case (SOILC_CORPSE, SOILC_CORPSE_N)
        delta = (cc%bwood+cc%bsw)*frac_harvested*agf_bs*frac_wood_wasted_ag;
        vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + &
            [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*delta

        delta = (cc%bl+cc%blv) * frac_harvested;
        if(delta<0) call error_mesg('vegn_cut_forest', &
                'harvested amount of live biomass ('//string(delta)//' kgC/m2) is below zero', &
                FATAL)

        vegn%litter_buff_C(:,LEAF) = vegn%litter_buff_C(:,LEAF) + &
            [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*delta

        vegn%ssc_pool_bg = vegn%ssc_pool_bg + cc%br*frac_harvested*(1-sp%fsc_froot)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + cc%br*frac_harvested*sp%fsc_froot

        if(waste_below_ground_wood) then
          vegn%ssc_pool_bg = vegn%ssc_pool_bg + (cc%bwood+cc%bsw)*frac_harvested*(1-agf_bs)*(1-sp%fsc_wood)
          vegn%fsc_pool_bg = vegn%fsc_pool_bg + (cc%bwood+cc%bsw)*frac_harvested*(1-agf_bs)*sp%fsc_wood
        endif

        if (soil_carbon_option == SOILC_CORPSE_N) then
            vegn%litter_buff_N(:,CWOOD) = vegn%litter_buff_N(:,CWOOD) + (&
                  [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*cc%wood_N +&
                  [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*(cc%sapwood_N+cc%stored_N)&
               )*frac_harvested*agf_bs*frac_wood_wasted_ag
            vegn%litter_buff_N(:,LEAF) = vegn%litter_buff_N(:,LEAF) + &
               [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*cc%leaf_N*frac_harvested
            vegn%ssn_pool_bg = vegn%ssn_pool_bg + cc%root_N*frac_harvested*(1-sp%fsc_froot)
            vegn%fsn_pool_bg = vegn%fsn_pool_bg + cc%root_N*frac_harvested*sp%fsc_froot

            if(waste_below_ground_wood) then
              vegn%ssn_pool_bg = vegn%ssn_pool_bg + (cc%wood_N*(1-sp%fsc_wood)+(cc%sapwood_N+cc%stored_N)*(1-sp%fsc_liv))*frac_harvested*(1-agf_bs)
              vegn%fsn_pool_bg = vegn%fsn_pool_bg + (cc%wood_N*sp%fsc_wood+(cc%sapwood_N+cc%stored_N)*sp%fsc_liv)*frac_harvested*(1-agf_bs)
            endif
        endif
     case default
        call error_mesg('vegn_phenology','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select

     cc%bliving = cc%bliving*(1-frac_harvested);
     cc%bwood   = cc%bwood*(1-frac_harvested);
     cc%leaf_N = cc%leaf_N*(1-frac_harvested)
     cc%root_N = cc%root_N*(1-frac_harvested)
     cc%wood_N = cc%wood_N*(1-frac_harvested)
     cc%sapwood_N = cc%sapwood_N*(1-frac_harvested)
     cc%stored_N = cc%stored_N*(1-frac_harvested)
     ! Should stored N be lost or retained?
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc);
     end associate
  enddo
end subroutine vegn_cut_forest

end module
