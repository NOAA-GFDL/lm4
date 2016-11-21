module vegn_harvesting_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use constants_mod, only : tfreeze
use fms_mod, only : string, error_mesg, FATAL, NOTE, &
     mpp_pe, file_exist, close_file, &
     check_nml_error, stdlog, mpp_root_pe
use mpp_io_mod, only : axistype, mpp_get_atts, mpp_get_axis_data, &
     mpp_open, mpp_close, MPP_RDONLY, MPP_WRONLY, MPP_ASCII
use land_data_mod, only : log_version
use vegn_data_mod, only : do_ppa, &
     N_LU_TYPES, LU_PAST, LU_CROP, LU_NTRL, LU_SCND, &
     HARV_POOL_PAST, HARV_POOL_CROP, HARV_POOL_CLEARED, HARV_POOL_WOOD_FAST, &
     HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW, &
     agf_bs, fsc_liv, fsc_froot, fsc_wood
use vegn_tile_mod, only : vegn_tile_type, vegn_relayer_cohorts_ppa
use vegn_disturbance_mod, only : kill_plants_ppa
use vegn_cohort_mod, only : &
     vegn_cohort_type, update_biomass_pools
use soil_tile_mod, only: num_l, LEAF, CWOOD
use soil_carbon_mod, only: soil_carbon_option, &
     SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, N_C_TYPES, C_CEL, C_LIG

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

! TODO: possibly move all definition of cpw,clw,csw in one place
real, parameter :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

! ==== module data ===========================================================

! ---- namelist variables ----------------------------------------------------
logical :: do_harvesting       = .TRUE.  ! if true, then harvesting of crops and pastures is done
real :: grazing_intensity      = 0.25    ! fraction of biomass removed each time by grazing
real :: grazing_residue        = 0.1     ! fraction of the grazed biomass transferred into soil pools
real :: frac_wood_wasted_harv  = 0.25    ! fraction of wood wasted while harvesting
real :: frac_wood_wasted_clear = 0.25    ! fraction of wood wasted while clearing land for pastures or crops
logical :: waste_below_ground_wood = .TRUE. ! If true, all the wood below ground (1-agf_bs fraction of bwood
        ! and bsw) is wasted. Old behavior assumed this to be FALSE. NOTE that in PPA this option has no
        ! effect; below-ground wood is always wasted.
real :: wood_harv_DBH          = 0.05    ! DBH above which trees are harvested, m
real :: frac_trampled          = 0.9     ! fraction of small trees that get trampled during harvesting
real :: frac_wood_fast         = ONETHIRD ! fraction of wood consumed fast
real :: frac_wood_med          = ONETHIRD ! fraction of wood consumed with medium speed
real :: frac_wood_slow         = ONETHIRD ! fraction of wood consumed slowly
real :: crop_seed_density      = 0.1     ! biomass of seeds left after crop harvesting, kg/m2
namelist/harvesting_nml/ do_harvesting, grazing_intensity, grazing_residue, &
     frac_wood_wasted_harv, frac_wood_wasted_clear, waste_below_ground_wood, &
     wood_harv_DBH, frac_trampled, &
     frac_wood_fast, frac_wood_med, frac_wood_slow, &
     crop_seed_density

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
     ierr = 1
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
end subroutine vegn_harvesting_init


! ============================================================================
subroutine vegn_harvesting_end
end subroutine vegn_harvesting_end


! ============================================================================
! harvest vegetation in a tile
subroutine vegn_harvesting(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  if (.not.do_harvesting) &
       return ! do nothing if no harvesting requested

  select case(vegn%landuse)
  case(LU_PAST)  ! pasture
     call vegn_graze_pasture    (vegn)
  case(LU_CROP)  ! crop
     call vegn_harvest_cropland (vegn)
  end select
end subroutine


! ============================================================================
subroutine vegn_graze_pasture(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  if (do_ppa) then
     call vegn_graze_pasture_ppa(vegn)
  else
     call vegn_graze_pasture_lm3(vegn)
  endif
end subroutine vegn_graze_pasture

! ============================================================================
subroutine vegn_harvest_cropland(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  if (do_ppa) then
     call vegn_harvest_crop_ppa(vegn)
  else
     call vegn_harvest_crop_lm3(vegn)
  endif
end subroutine vegn_harvest_cropland

! ============================================================================
subroutine vegn_cut_forest(vegn, new_landuse)
  type(vegn_tile_type), intent(inout) :: vegn
  integer, intent(in) :: new_landuse ! new land use type that gets assigned to
                                     ! the tile after the wood harvesting

  if (do_ppa) then
     call vegn_cut_forest_ppa(vegn, new_landuse)
  else
     call vegn_cut_forest_lm3(vegn, new_landuse)
  endif
end subroutine vegn_cut_forest

! ============================================================================
subroutine vegn_graze_pasture_lm3(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real ::  bdead0, balive0, bleaf0, bfroot0, btotal0 ! initial combined biomass pools
  real ::  bdead1, balive1, bleaf1, bfroot1, btotal1 ! updated combined biomass pools
  type(vegn_cohort_type), pointer :: cc ! shorthand for the current cohort
  integer :: i

  balive0 = 0 ; balive1 = 0
  bdead0  = 0 ; bdead1  = 0
  bleaf0  = 0 ; bleaf1  = 0
  bfroot0 = 0 ; bfroot1 = 0

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1,vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive0 = balive0 + cc%bl + cc%blv + cc%br
     bleaf0  = bleaf0  + cc%bl + cc%blv
     bfroot0 = bfroot0 + cc%br
     bdead0  = bdead0  + cc%bwood + cc%bsw
     ! only potential leaves are consumed
     vegn%harv_pool(HARV_POOL_PAST) = vegn%harv_pool(HARV_POOL_PAST) + &
          cc%bliving*cc%Pl*grazing_intensity*(1-grazing_residue)
     cc%bliving = cc%bliving - cc%bliving*cc%Pl*grazing_intensity

     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc)

     ! calculate new combined vegetation biomass pools
     balive1 = balive1 + cc%bl + cc%blv + cc%br
     bleaf1  = bleaf1  + cc%bl + cc%blv
     bfroot1 = bfroot1 + cc%br
     bdead1  = bdead1  + cc%bwood + cc%bsw
  enddo
  btotal0 = balive0 + bdead0
  btotal1 = balive1 + bdead1

  ! update intermediate soil carbon pools
  select case(soil_carbon_option)
  case(SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
     vegn%fsc_pool_bg = vegn%fsc_pool_bg + &
          (fsc_liv*(balive0-balive1)+fsc_wood*(bdead0-bdead1))*grazing_residue
     vegn%ssc_pool_bg = vegn%ssc_pool_bg + &
          ((1-fsc_liv)*(balive0-balive1)+ (1-fsc_wood)*(bdead0-bdead1))*grazing_residue
  case(SOILC_CORPSE)
     vegn%litter_buff_C(:,LEAF) = vegn%litter_buff_C(:,LEAF) + &
          [fsc_liv, 1-fsc_liv, 0.0]*(bleaf0-bleaf1)*grazing_residue
     vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + &
          [fsc_wood, 1-fsc_wood, 0.0]*agf_bs*(bdead0-bdead1)*grazing_residue
     vegn%fsc_pool_bg = vegn%fsc_pool_bg + grazing_residue*( &
          fsc_froot*(bfroot0-bfroot1) + &
          fsc_liv*(bleaf0-bleaf1) + &
          (1-agf_bs)*fsc_wood*(bdead0-bdead1))
     vegn%ssc_pool_bg = vegn%ssc_pool_bg + grazing_residue*(&
          (1-fsc_froot)*(bfroot0-bfroot1) + &
          (1-fsc_liv)*(bleaf0-bleaf1) + &
          (1-agf_bs)*(1-fsc_wood)*(bdead0-bdead1))
  case default
     call error_mesg('vegn_graze_pasture_lm3','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select
end subroutine vegn_graze_pasture_lm3

! ================================================================================
subroutine vegn_harvest_crop_lm3(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  real :: fraction_harvested;    ! fraction of biomass harvested this time
  real :: bdead, balive, btotal; ! combined biomass pools
  integer :: i

  balive = 0 ; bdead = 0
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
  enddo
  btotal = balive+bdead

  ! calculate harvested fraction: cut everything down to seed level
  fraction_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0)

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! use for harvest only above-ground living biomass and waste the correspondent below living and wood
     vegn%harv_pool(HARV_POOL_CROP) = vegn%harv_pool(HARV_POOL_CROP) + &
          cc%bliving*(cc%Pl + cc%Psw*agf_bs)*fraction_harvested
     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + fraction_harvested*(fsc_liv*cc%bliving*cc%Pr + &
             fsc_wood*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)))
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + fraction_harvested*((1-fsc_liv)*cc%bliving*cc%Pr + &
             (1-fsc_wood)*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)))
     case (SOILC_CORPSE)
        vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + &
               [fsc_wood, 1-fsc_wood, 0.0] * fraction_harvested*agf_bs*cc%bwood

        vegn%fsc_pool_bg = vegn%fsc_pool_bg + fraction_harvested*(&
               fsc_froot*cc%bliving*cc%Pr + &
               (1-agf_bs)*fsc_wood*(cc%bwood + cc%bliving*cc%Psw))
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + fraction_harvested*(&
               (1-fsc_froot)*cc%bliving*cc%Pr + &
               (1-agf_bs)*(1-fsc_wood)*(cc%bwood + cc%bliving*cc%Psw))
     case default
        call error_mesg('vegn_harvest_crop_lm3','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select
     cc%bliving = cc%bliving * (1-fraction_harvested)
     cc%bwood   = cc%bwood   * (1-fraction_harvested)
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc)
  enddo
end subroutine vegn_harvest_crop_lm3

! ============================================================================
! for now cutting forest is the same as harvesting cropland --
! we basically cut down everything, leaving only seeds
subroutine vegn_cut_forest_lm3(vegn, new_landuse)
  type(vegn_tile_type), intent(inout) :: vegn
  integer, intent(in) :: new_landuse ! new land use type that gets assigned to
                                     ! the tile after the wood harvesting

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  real :: frac_harvested         ! fraction of biomass harvested this time
  real :: frac_wood_wasted       ! fraction of wood wasted during transition
  real :: frac_wood_wasted_ag    ! fraction of above-ground wood wasted during transition
  real :: wood_harvested         ! amount of harvested wood, kgC/m2
  real :: bdead, balive, bleaf, bfroot, btotal; ! combined biomass pools
  real :: delta
  integer :: i

  balive = 0 ; bdead = 0 ; bleaf = 0 ; bfroot = 0
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bleaf  = bleaf  + cc%bl + cc%blv
     bfroot = bfroot + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
  enddo
  btotal = balive+bdead

  ! calculate harvested fraction: cut everything down to seed level
  frac_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0)

  ! define fraction of wood wasted, based on the transition type
  if (new_landuse==LU_SCND) then
     frac_wood_wasted = frac_wood_wasted_harv
  else
     frac_wood_wasted = frac_wood_wasted_clear
  endif
  ! take into account that all wood below ground is wasted; also the fraction
  ! of waste calculated above is lost from the above-ground part of the wood
  frac_wood_wasted_ag=frac_wood_wasted
  if (waste_below_ground_wood) then
     frac_wood_wasted = (1-agf_bs) + agf_bs*frac_wood_wasted
  endif

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)

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
     ! soil carbon pools according to fractions specified thorough the namelists
     delta = (cc%bwood+cc%bsw)*frac_harvested*frac_wood_wasted
     if(delta<0) call error_mesg('vegn_cut_forest_lm3', &
          'harvested amount of dead biomass ('//string(delta)//' kgC/m2) is below zero', &
          FATAL)

     select case (soil_carbon_option)
     case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + delta*(1-fsc_wood)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + delta*   fsc_wood

        delta = balive * frac_harvested
        if(delta<0) call error_mesg('vegn_cut_forest_lm3', &
             'harvested amount of live biomass ('//string(delta)//' kgC/m2) is below zero', &
             FATAL)
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + delta*(1-fsc_liv)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + delta*   fsc_liv
     case (SOILC_CORPSE)
        delta = (cc%bwood+cc%bsw)*frac_harvested*agf_bs*frac_wood_wasted_ag
        vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + &
            [fsc_wood, 1-fsc_wood, 0.0]*delta

        delta = (cc%bl+cc%blv) * frac_harvested
        if(delta<0) call error_mesg('vegn_cut_forest_lm3', &
                'harvested amount of live biomass ('//string(delta)//' kgC/m2) is below zero', &
                FATAL)

        vegn%litter_buff_C(:,LEAF) = vegn%litter_buff_C(:,LEAF) + &
            [fsc_liv, 1-fsc_liv, 0.0]*delta

        vegn%ssc_pool_bg = vegn%ssc_pool_bg + cc%br*frac_harvested*(1-fsc_froot)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + cc%br*frac_harvested*fsc_froot

        if(waste_below_ground_wood) then
          vegn%ssc_pool_bg = vegn%ssc_pool_bg + (cc%bwood+cc%bsw)*frac_harvested*(1-agf_bs)*(1-fsc_wood)
          vegn%fsc_pool_bg = vegn%fsc_pool_bg + (cc%bwood+cc%bsw)*frac_harvested*(1-agf_bs)*fsc_wood
        endif
     case default
        call error_mesg('vegn_cut_forest_lm3','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select

     cc%bliving = cc%bliving*(1-frac_harvested)
     cc%bwood   = cc%bwood*(1-frac_harvested)
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc)
  enddo
end subroutine vegn_cut_forest_lm3

! ============================================================================
subroutine vegn_graze_pasture_ppa(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  real    :: ag_litt ! accumulated above-ground litter, kg C/m2
  integer :: i

  ag_litt = 0.0
  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i))
     vegn%harv_pool(HARV_POOL_PAST) = vegn%harv_pool(HARV_POOL_PAST) + &
          cc%bl*grazing_intensity*(1-grazing_residue)*cc%nindivs
     ag_litt = ag_litt + cc%bl*grazing_intensity*grazing_residue*cc%nindivs
     cc%bl = cc%bl*(1-grazing_intensity)
     end associate
  enddo

  ! add carbon to intermediate pools
  select case (soil_carbon_option)
  case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
     ! add litter to intermediate carbon pools
     vegn%fsc_pool_ag = vegn%fsc_pool_ag + ag_litt*fsc_liv
     vegn%ssc_pool_ag = vegn%ssc_pool_ag + ag_litt*(1-fsc_liv)
     ! note that we assume microbial biomass is zero
  case (SOILC_CORPSE)
     vegn%litter_buff_C(:,LEAF) = vegn%litter_buff_C(:,LEAF) + ag_litt*[fsc_liv,1-fsc_liv,0.0]
  case default
     call error_mesg('vegn_graze_pasture_ppa','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select

end subroutine vegn_graze_pasture_ppa

! ============================================================================
! NOTE that the PPA harvest would not work properly if applied only once per year, because
! at the end of the year the leaves are down in NH, and therefore harvest amount will be
! unrealistically small. The same is true about the grazing.
subroutine vegn_harvest_crop_ppa(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  real, dimension(N_C_TYPES) :: &
     ag_litt(N_C_TYPES), & ! accumulated above-ground litter, kg C/m2
     bg_litt(N_C_TYPES)    ! accumulated below-ground litter, kg C/m2
  real :: ndead  ! number of individuals killed in cohort
  integer :: i

  ag_litt = 0.0; bg_litt = 0.0
  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i))
     ndead = cc%nindivs ! harvest everything
     ! take care of water and energy conservation
     vegn%drop_wl = vegn%drop_wl + cc%wl*ndead
     vegn%drop_ws = vegn%drop_ws + cc%ws*ndead
     vegn%drop_hl = vegn%drop_hl + clw*cc%wl*ndead*(cc%Tv-tfreeze)
     vegn%drop_hs = vegn%drop_hs + csw*cc%ws*ndead*(cc%Tv-tfreeze)

     ! add C to harvest and litter pools. This is slightly different from kill_plants_ppa
     ! because above-ground part of nsc and sapwood are getting harvested, rather than going
     ! into litter
     vegn%harv_pool(HARV_POOL_CROP) = vegn%harv_pool(HARV_POOL_CROP) + ndead*( &
         cc%bl + cc%bseed + cc%carbon_gain + cc%growth_previous_day + &
         cc%bsw*agf_bs + cc%nsc*agf_bs &
         )
     ag_litt(:) = ag_litt(:) + [fsc_wood, 1-fsc_wood, 0.0]*(cc%bwood+cc%bwood_gain)*agf_bs*ndead
     bg_litt(:) = bg_litt(:) + ndead * [ &
         fsc_froot    *cc%br + fsc_wood    *(cc%bsw+cc%bwood+cc%bwood_gain)*(1-agf_bs) + cc%nsc*(1-agf_bs), &
         (1-fsc_froot)*cc%br + (1-fsc_wood)*(cc%bsw+cc%bwood+cc%bwood_gain)*(1-agf_bs), &
         0.0 ]

     ! reduce the number of individuals in cohort
     cc%nindivs = cc%nindivs-ndead
     end associate
  enddo

  ! add carbon to intermediate pools
  select case (soil_carbon_option)
  case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
     ! add litter to intermediate carbon pools
     vegn%fsc_pool_ag = vegn%fsc_pool_ag + ag_litt(C_CEL)
     vegn%ssc_pool_ag = vegn%ssc_pool_ag + ag_litt(C_LIG)
     ! note that we assume microbial biomass is zero
  case (SOILC_CORPSE)
     vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + ag_litt(:)
  case default
     call error_mesg('vegn_harvest_crop_ppa','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select
  vegn%fsc_pool_bg = vegn%fsc_pool_bg + bg_litt(C_CEL)
  vegn%ssc_pool_bg = vegn%ssc_pool_bg + bg_litt(C_LIG)

  call vegn_relayer_cohorts_ppa(vegn)

end subroutine vegn_harvest_crop_ppa

! ============================================================================
subroutine vegn_cut_forest_ppa(vegn, new_landuse)
  type(vegn_tile_type), intent(inout) :: vegn
  integer, intent(in) :: new_landuse ! new land use type that gets assigned to
                                     ! the tile after the wood harvesting

  ! ---- local vars
  real :: frac_wood_wasted       ! fraction of wood wasted during transition
  real, dimension(N_C_TYPES) :: &
     leaf_litt(N_C_TYPES), & ! accumulated leaf litter, kg C/m2
     wood_harv(N_C_TYPES), & ! accumulated wood harvest, kg C/m2
     wood_litt(N_C_TYPES)    ! accumulated wood litter, kg C/m2
  real :: root_litt(num_l, N_C_TYPES) ! accumulated root litter per soil layer, kgC/m2
  real :: ndead  ! number of individuals killed in cohort
  real :: dbh_min
  integer :: i

  ! define fraction of wood wasted, based on the transition type
  if (new_landuse==LU_SCND) then
     ! this is wood haresting
     frac_wood_wasted = frac_wood_wasted_harv
     dbh_min = wood_harv_DBH
  else
     ! this is land clearance
     frac_wood_wasted = frac_wood_wasted_clear
     dbh_min = -HUGE(1.0) ! in clearance, everything is harvested
  endif

  leaf_litt=0.0; wood_harv=0.0; wood_litt=0.0; root_litt=0.0
  do i = 1, vegn%n_cohorts
     associate (cc=>vegn%cohorts(i))
     if (cc%dbh > dbh_min) then
        ! these trees are harvested
        ndead = cc%nindivs
        call kill_plants_ppa(cc,vegn,ndead,0.0, leaf_litt, wood_harv, root_litt)
     else
        ! these trees are too small to be harvested, so a part of them get trampled
        ! and goes to waste, the rest stays
        ndead = cc%nindivs * frac_trampled
        call kill_plants_ppa(cc,vegn,ndead,0.0, leaf_litt, wood_litt, root_litt)
     endif
     end associate
     ! note that below-ground wood all goes to litter, by construction of kill_plants_ppa
  enddo

  ! distribute harvested wood between pools
  if (new_landuse==LU_SCND) then
     ! this is harvesting, distribute between 3 different wood pools
     vegn%harv_pool(HARV_POOL_WOOD_FAST) = vegn%harv_pool(HARV_POOL_WOOD_FAST) &
          + sum(wood_harv)*frac_wood_fast*(1-frac_wood_wasted)
     vegn%harv_pool(HARV_POOL_WOOD_MED) = vegn%harv_pool(HARV_POOL_WOOD_MED) &
          + sum(wood_harv)*frac_wood_med*(1-frac_wood_wasted)
     vegn%harv_pool(HARV_POOL_WOOD_SLOW) = vegn%harv_pool(HARV_POOL_WOOD_SLOW) &
          + sum(wood_harv)*frac_wood_slow*(1-frac_wood_wasted)
  else
     ! this is land clearance: everything goes into "cleared" pool
     vegn%harv_pool(HARV_POOL_CLEARED) = vegn%harv_pool(HARV_POOL_CLEARED) &
          + sum(wood_harv)*(1-frac_wood_wasted)
  endif

  select case (soil_carbon_option)
  case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
     ! add litter to intermediate carbon pools
     vegn%fsc_pool_ag = vegn%fsc_pool_ag + &
         wood_harv(C_CEL)*frac_wood_wasted + wood_litt(C_CEL) + leaf_litt(C_CEL)
     vegn%ssc_pool_ag = vegn%ssc_pool_ag + &
         wood_harv(C_LIG)*frac_wood_wasted + wood_litt(C_LIG) + leaf_litt(C_LIG)
     ! note that we assume microbial biomass is zero
  case (SOILC_CORPSE)
     vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + &
        wood_litt(:) + wood_harv(:)*frac_wood_wasted
     vegn%litter_buff_C(:,LEAF) = vegn%litter_buff_C(:,LEAF) + &
        leaf_litt(:)
  case default
     call error_mesg('vegn_cut_forest_ppa','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select
  vegn%fsc_pool_bg = vegn%fsc_pool_bg + sum(root_litt(:,C_CEL))
  vegn%ssc_pool_bg = vegn%ssc_pool_bg + sum(root_litt(:,C_LIG))

  call vegn_relayer_cohorts_ppa(vegn)

end subroutine vegn_cut_forest_ppa

end module
