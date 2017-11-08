module vegn_harvesting_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

#include "../shared/debug.inc"

use constants_mod, only : tfreeze
use fms_mod, only : string, error_mesg, FATAL, NOTE, &
     mpp_pe, file_exist, close_file, &
     check_nml_error, stdlog, mpp_root_pe, lowercase
use diag_manager_mod, only : register_static_field, send_data

use land_io_mod, only : read_field
use land_debug_mod, only : check_conservation, carbon_cons_tol, nitrogen_cons_tol, &
     land_error_message, do_check_conservation
use land_utils_mod, only : check_conservation_1, check_conservation_2
use land_data_mod, only : log_version, lnd
use vegn_data_mod, only : do_ppa, &
     N_LU_TYPES, LU_PAST, LU_CROP, LU_NTRL, LU_SCND, LU_RANGE, &
     HARV_POOL_PAST, HARV_POOL_CROP, HARV_POOL_CLEARED, HARV_POOL_WOOD_FAST, &
     HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW, &
     nspecies, spdata, agf_bs, LEAF_OFF
use land_tile_mod, only : land_tile_type
use soil_tile_mod, only : num_l, LEAF, CWOOD, soil_tile_type, &
     soil_tile_carbon, soil_tile_nitrogen
use vegn_tile_mod, only : vegn_tile_type, vegn_relayer_cohorts_ppa, vegn_tile_lai, &
     vegn_tile_carbon, vegn_tile_nitrogen
use soil_util_mod, only : add_root_litter
use vegn_cohort_mod, only : vegn_cohort_type, update_biomass_pools
use vegn_util_mod, only : kill_plants_ppa, add_seedlings_ppa
use soil_carbon_mod, only: soil_carbon_option, add_litter, C_FAST, C_SLOW, C_MIC, &
     SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, SOILC_CORPSE_N, N_C_TYPES, &
     deadmic_slow_frac

implicit none
private

! ==== public interface ======================================================
public :: vegn_harvesting_init
public :: vegn_harvesting_end

public :: vegn_harvesting

public :: vegn_cut_forest
! ==== end of public interface ===============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_harvesting_mod'
#include "../shared/version_variable.inc"
real, parameter :: ONETHIRD = 1.0/3.0
integer, parameter :: DAILY = 1, ANNUAL = 2
integer, parameter :: CROP_SCHEDULE_LM3 = 1, CROP_SCHEDULE_PRESCRIBED = 2

! TODO: possibly move all definition of cpw,clw,csw in one place
real, parameter :: &
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

! ==== module data ===========================================================

! ---- namelist variables ----------------------------------------------------
logical, public, protected  :: do_harvesting       = .TRUE.  ! if true, then harvesting of crops and pastures is done

character(16) :: grazing_frequency = 'daily' ! or 'annual'
real :: grazing_intensity_past  = 3.65 ! fraction of leaf biomass removed by grazing annually, roughly 1% per day.
real :: grazing_intensity_range = 3.65 ! fraction of leaf biomass removed by grazing annually, roughly 1% per day.
  ! NOTE that for daily grazing, grazing_intensity/365 fraction of leaf biomass is removed
  ! every day. E.g. if the desired intensity is 1% of leaves per day, set grazing_intensity
  ! to 3.65
real :: grazing_residue        = 0.1     ! fraction of the grazed biomass transferred into soil pools
real :: min_lai_for_grazing_past  = 0.0     ! no grazing if LAI lower than this threshold
real :: min_lai_for_grazing_range = 0.0     ! no grazing if LAI lower than this threshold
  ! NOTE that in CORPSE mode regardless of the grazing frequency soil carbon input from
  ! grazing still goes to intermediate pools, and then it is transferred from
  ! these pools to soil/litter carbon pools with constant rates over the next year.
  ! In CENTURY mode grazing residue is deposited to soil directly in case of daily
  ! grazing frequency; still goes through intermediate pools in case of annual grazing.
real :: max_grazing_height_past  = 9999.0  ! m, no grazing of vegetation above this height.
real :: max_grazing_height_range = 3.0     ! m, no grazing of vegetation above this height.

real :: wood_harv_DBH          = 0.05    ! DBH above which trees are harvested in PPA, m
real :: frac_trampled          = 0.9     ! fraction of small trees that get trampled during harvesting in PPA
real :: frac_wood_wasted_harv  = 0.25    ! fraction of wood wasted while harvesting
real :: frac_wood_wasted_clear = 0.25    ! fraction of wood wasted while clearing land for pastures or crops
logical :: waste_below_ground_wood = .TRUE. ! If true, all the wood below ground (1-agf_bs fraction of bwood
        ! and bsw) is wasted. Old behavior assumed this to be FALSE. NOTE that in PPA this option has no
        ! effect; below-ground wood is always wasted.
real :: frac_wood_fast         = ONETHIRD ! fraction of wood consumed fast
real :: frac_wood_med          = ONETHIRD ! fraction of wood consumed with medium speed
real :: frac_wood_slow         = ONETHIRD ! fraction of wood consumed slowly

character(16) :: crop_schedule = 'lm3' ! or 'prescribed'
character(len=1024) :: crop_schedule_file  = '' ! input data set of crop planting and harvesting dates
real :: crop_seed_density      = 0.1   ! biomass of seeds left after crop harvesting, kg/m2
character(32) :: crop_species  = ''    ! name of the species used for crops

namelist/harvesting_nml/ do_harvesting, &
     ! pasture and rangeland grazing parameters
     grazing_frequency,  &
     grazing_residue,    &
     grazing_intensity_past, grazing_intensity_range, &
     max_grazing_height_past, max_grazing_height_range, &
     min_lai_for_grazing_past, min_lai_for_grazing_range, &
     ! wood harvesting and clearance parameters
     wood_harv_DBH, frac_trampled, &
     frac_wood_wasted_harv, frac_wood_wasted_clear, waste_below_ground_wood, &
     frac_wood_fast, frac_wood_med, frac_wood_slow, &
     ! crop harvesting and planting parameters
     crop_species, crop_seed_density, crop_schedule, crop_schedule_file

integer :: grazing_freq = -1 ! indicator of grazing frequency (ANNUAL or DAILY)
integer :: crop_schedule_option = -1 ! selected planting/harvesting schedule option
real, allocatable :: crop_planting_day(:) ! day of year when planting is done
real, allocatable :: crop_harvest_day(:)  ! day of year when harvesting is done
integer :: crop_species_idx = -1 ! index of crop species

integer :: id_crop_planting_day, id_crop_harvest_day

contains ! ###################################################################

! ============================================================================
subroutine vegn_harvesting_init(id_ug)
  integer, intent(in) :: id_ug ! id of the unstructured grid diagnostic axis

  integer :: unit, ierr, io, i
  logical :: used

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
  ! parse the grazing frequency parameter
  select case(lowercase(grazing_frequency))
  case('annual')
     grazing_freq = ANNUAL
  case('daily')
     grazing_freq = DAILY
     ! scale grazing intensity for daily frequency
     grazing_intensity_past  = grazing_intensity_past/365.0
     grazing_intensity_range = grazing_intensity_range/365.0
  case default
     call error_mesg('vegn_harvesting_init','grazing_frequency must be "annual" or "daily"',FATAL)
  end select

  allocate (crop_planting_day(lnd%ls:lnd%le),crop_harvest_day(lnd%ls:lnd%le))
  select case(trim(lowercase(crop_schedule)))
  case('lm3')
     crop_schedule_option = CROP_SCHEDULE_LM3
     crop_planting_day = 1.0
     crop_harvest_day  = 1.0
  case('prescribed')
     crop_schedule_option = CROP_SCHEDULE_PRESCRIBED
     ! read input data
     call read_field( crop_schedule_file, 'plantingdy', crop_planting_day, interp='nearest' )
     call read_field( crop_schedule_file, 'harvestdy',  crop_harvest_day,  interp='nearest' )
  case default
     call error_mesg('vegn_harvesting_init','crop_schedule must be "lm3" or "prescribed"',FATAL)
  end select

  id_crop_harvest_day = register_static_field ( 'vegn', 'crop_harvest_day', (/id_ug/), &
         'day of year when cropa are harvested', 'day', missing_value = -1.0 )
  id_crop_planting_day = register_static_field ( 'vegn', 'crop_planting_day', (/id_ug/), &
         'day of year when crops are planted', 'day', missing_value = -1.0 )

  if ( id_crop_harvest_day > 0 )  used = send_data ( id_crop_harvest_day, crop_harvest_day, lnd%time )
  if ( id_crop_planting_day > 0 ) used = send_data ( id_crop_planting_day, crop_planting_day, lnd%time )

  ! find crop species; it's OK if it is not found, perhaps we do not need it -- e.g.
  ! we are running potential vegetation, or in LM3 mode where it is not used (yet)
  ! For now, use single species for crop everywhere. This will obviously need to be changed
  ! when we switch to more sophisticated treatment of agriculture.
  do i = 0, nspecies-1
     if (trim(spdata(i)%name)==trim(crop_species)) then
        crop_species_idx = i
        exit ! from loop
     endif
  enddo
  if (crop_species_idx<0) then
     call error_mesg('vegn_harvesting_init','crop species "'//trim(crop_species)//'" not found in the list of species',NOTE)
  else
     call error_mesg('vegn_harvesting_init','crop species "'//trim(crop_species)//'" is species #'//string(crop_species_idx)//&
                     ' in the list of species', NOTE)
  endif
end subroutine vegn_harvesting_init


! ============================================================================
subroutine vegn_harvesting_end
   if (allocated(crop_harvest_day))  deallocate(crop_harvest_day)
   if (allocated(crop_planting_day)) deallocate(crop_planting_day)
end subroutine vegn_harvesting_end


! ============================================================================
! harvest vegetation in a tile
subroutine vegn_harvesting(tile, end_of_year, end_of_month, end_of_day, day_of_year, l)
  type(land_tile_type), intent(inout) :: tile
  logical, intent(in) :: end_of_year, end_of_month, end_of_day ! indicators of respective period boundaries
  integer, intent(in) :: day_of_year ! current day of year
  integer, intent(in) :: l ! index of current grid cell in unstructured grid

  if (.not.do_harvesting) return ! do nothing if no harvesting requested

  associate(vegn=>tile%vegn)
  select case(vegn%landuse)
  case(LU_PAST)  ! pasture
     if ((end_of_day  .and. grazing_freq==DAILY).or. &
         (end_of_year .and. grazing_freq==ANNUAL)) then
        call vegn_graze_pasture (tile)
     endif
  case(LU_RANGE)  ! rangeland
     if ((end_of_day  .and. grazing_freq==DAILY).or. &
         (end_of_year .and. grazing_freq==ANNUAL)) then
        call vegn_graze_rangeland (tile)
     endif
  case(LU_CROP)  ! crop
     select case(crop_schedule_option)
     case (CROP_SCHEDULE_LM3)
        if (end_of_year) call vegn_harvest_cropland (tile)
     case (CROP_SCHEDULE_PRESCRIBED)
        if (end_of_day.AND.day_of_year==nint(crop_harvest_day(l))) then
           call vegn_harvest_cropland (tile)
        endif
        if (end_of_day.AND.day_of_year==nint(crop_planting_day(l))) then
           call vegn_plant_crop (tile)
        endif
     end select ! crop_schedule_option
  end select
  end associate
end subroutine vegn_harvesting


! ============================================================================
subroutine vegn_graze_pasture(tile)
  type(land_tile_type), intent(inout) :: tile

  if (do_ppa) then
     call vegn_graze_pasture_ppa(tile, min_lai_for_grazing_past, grazing_intensity_past, max_grazing_height_past)
  else
     call vegn_graze_pasture_lm3(tile, min_lai_for_grazing_past, grazing_intensity_past)
  endif
end subroutine vegn_graze_pasture


! ============================================================================
subroutine vegn_graze_rangeland(tile)
  type(land_tile_type), intent(inout) :: tile

  if (do_ppa) then
     call vegn_graze_pasture_ppa(tile, min_lai_for_grazing_range, grazing_intensity_range, max_grazing_height_range)
  else
     call vegn_graze_pasture_lm3(tile, min_lai_for_grazing_range, grazing_intensity_range)
  endif
end subroutine vegn_graze_rangeland

! ============================================================================
subroutine vegn_harvest_cropland(tile)
  type(land_tile_type), intent(inout) :: tile

  if (do_ppa) then
     call vegn_harvest_crop_ppa(tile)
  else
     call vegn_harvest_crop_lm3(tile)
  endif
end subroutine vegn_harvest_cropland


! ============================================================================
subroutine vegn_plant_crop(tile)
  type(land_tile_type), intent(inout) :: tile

  if (do_ppa) then
     call vegn_plant_crop_ppa(tile)
  else
     ! do nothing at the moment -- later add turning phenology on
  endif
end subroutine vegn_plant_crop

! ============================================================================
subroutine vegn_cut_forest(tile, new_landuse)
  type(land_tile_type), intent(inout) :: tile
  integer, intent(in) :: new_landuse ! new land use type that gets assigned to
                                     ! the tile after the wood harvesting

  if (do_ppa) then
     call vegn_cut_forest_ppa(tile, new_landuse)
  else
     call vegn_cut_forest_lm3(tile, new_landuse)
  endif
end subroutine vegn_cut_forest

! ============================================================================
subroutine vegn_graze_pasture_lm3(tile, min_lai_for_grazing, grazing_intensity)
  type(land_tile_type), intent(inout) :: tile
  real, intent(in) :: min_lai_for_grazing
  real, intent(in) :: grazing_intensity

  ! ---- local vars
  real ::  bdead0, balive0, bleaf0, blv0, bfroot0 ! initial combined biomass pools
  real ::  bdead1, balive1, bleaf1, blv1, bfroot1 ! updated combined biomass pools
  integer :: i
  real :: carbon_lost
  real :: delta_leaf, delta_root, delta_wood
  real,dimension(N_C_TYPES) :: leaflitter_C,woodlitter_C,bglitter_C,leaflitter_N,woodlitter_N,bglitter_N
  real :: wood_n2c

  associate(vegn=>tile%vegn,soil=>tile%soil)
  ! do nothing if the LAI is less that the lower grazing limit
  if ( vegn_tile_LAI(vegn) <= min_lai_for_grazing ) return

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1,vegn%n_cohorts
     associate (cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))

     ! This makes sure biomass pools are correct before calculating changes
     ! in leaf biomass and such, just in case it wasn't called before
     call update_biomass_pools(cc);

     ! In multiple-cohort scenario, should we be adding over all cohorts?
     if(cc%bliving*cc%Pl/sp%LMA < min_lai_for_grazing) continue

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
        blv0 = cc%blv - cc%bliving*(cc%Pl+cc%Pr) ! Excess carbon (due to N limitation)
     endif
     bdead0  =  cc%bwood + cc%bsw
     ! only potential leaves are consumed
     carbon_lost=cc%bliving*cc%Pl*grazing_intensity
     vegn%harv_pool(HARV_POOL_PAST) = vegn%harv_pool(HARV_POOL_PAST) + &
          carbon_lost*(1-grazing_residue) ;
     cc%bliving = cc%bliving - carbon_lost;

     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc);

     ! calculate new combined vegetation biomass pools
     balive1 =  cc%bl + cc%blv + cc%br
     bleaf1  =  cc%bliving*cc%Pl  !cc%bl + cc%blv
     bfroot1 =  cc%bliving*cc%Pr  !cc%br
     if(cc%bl+cc%br>0) then
        blv1 = cc%blv
     else
        blv1 = cc%blv - cc%bliving*(cc%Pl+cc%Pr) ! Excess carbon (due to N limitation)
     endif
     bdead1  = cc%bwood + cc%bsw

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

        ! We're not removing belowground portion of what was grazed, so that needs to be clawed back from harvest pool
        vegn%harv_pool(HARV_POOL_PAST) = vegn%harv_pool(HARV_POOL_PAST) - (1.0-grazing_residue)*(delta_root+(1-agf_bs)*delta_wood)

        if(soil_carbon_option == SOILC_CORPSE_N) then
           leaflitter_N=leaflitter_C/sp%leaf_live_c2n
           woodlitter_N=woodlitter_C/sp%leaf_live_c2n
           bglitter_N=(/grazing_residue*(sp%fsc_froot*(delta_root)/sp%froot_live_c2n +(1-agf_bs)*sp%fsc_wood*(delta_wood)*wood_n2c),&
                        grazing_residue*((1-sp%fsc_froot)*(delta_root)/sp%froot_live_c2n +  (1-agf_bs)*(1-sp%fsc_wood)*(delta_wood)*wood_n2c),&
                        0.0/)

           cc%stored_N = cc%stored_N - delta_leaf/sp%leaf_live_c2n - delta_wood*wood_n2c - delta_root/sp%froot_live_c2n
           vegn%harv_pool_nitrogen(HARV_POOL_PAST) = vegn%harv_pool_nitrogen(HARV_POOL_PAST) + &
                delta_leaf/sp%leaf_live_c2n*(1-grazing_residue) + delta_wood*agf_bs*wood_n2c*(1-grazing_residue)
        else
           leaflitter_N = 0.0
           woodlitter_N = 0.0
           bglitter_N   = 0.0
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
  end associate ! vegn, soil
end subroutine vegn_graze_pasture_lm3


! ================================================================================
subroutine vegn_harvest_crop_lm3(tile)
  type(land_tile_type), intent(inout) :: tile

  ! ---- local vars
  real :: fraction_harvested;    ! fraction of biomass harvested this time
  real :: bdead, balive, btotal; ! combined biomass pools
  integer :: i

  associate (vegn=>tile%vegn)
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

           vegn%harv_pool_nitrogen(HARV_POOL_CROP) = vegn%harv_pool_nitrogen(HARV_POOL_CROP) + &
                (cc%bliving*cc%Pl/sp%leaf_live_c2n + agf_bs*cc%sapwood_N)*fraction_harvested

           ! Make sure stored nitrogen loss is sensible when leaves are off:
           ! Amounts harvested should be determined by potential leaf tissues, but this N would still be in stored pool so needs to be subtracted
           if (cc%status == LEAF_OFF) then
              vegn%harv_pool_nitrogen(HARV_POOL_CROP) = vegn%harv_pool_nitrogen(HARV_POOL_CROP) + &
                 (max(cc%stored_N,0.0)-(cc%bliving*cc%Pl/sp%leaf_live_c2n + cc%bliving*cc%Pr/sp%froot_live_c2n))*fraction_harvested
              cc%stored_N = cc%stored_N - &
                 (max(cc%stored_N,0.0)-(cc%bliving*cc%Pl/sp%leaf_live_c2n + cc%bliving*cc%Pr/sp%froot_live_c2n))*fraction_harvested
           else
              ! In this case stored N can just be lost because it doesn't contain N from harvested potential leaves and roots
              vegn%harv_pool_nitrogen(HARV_POOL_CROP) = vegn%harv_pool_nitrogen(HARV_POOL_CROP) + max(cc%stored_N,0.0)*fraction_harvested
              cc%stored_N = max(cc%stored_N,0.0)*(1-fraction_harvested)
           endif
           ! Subtract N lost from N pools. Leaf and root can get subtracted from stored and things will be rebalanced later.
           ! Some stored N also was lost
           cc%stored_N = cc%stored_N -  &
                      cc%bliving*fraction_harvested*(cc%Pl/sp%leaf_live_c2n + cc%Pr/sp%froot_live_c2n)
           cc%sapwood_N = cc%sapwood_N*(1-fraction_harvested)
           cc%wood_N = cc%wood_N*(1-fraction_harvested)
       endif

     case default
        call error_mesg('vegn_harvest_cropland','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select


     ! redistribute leftover biomass between biomass pools

     cc%bliving = cc%bliving * (1-fraction_harvested);
     cc%bwood   = cc%bwood   * (1-fraction_harvested);

     call update_biomass_pools(cc);
     end associate
  enddo
  end associate ! vegn
end subroutine vegn_harvest_crop_lm3


! ============================================================================
! for now cutting forest is the same as harvesting cropland --
! we basically cut down everything, leaving only seeds
subroutine vegn_cut_forest_lm3(tile, new_landuse)
  type(land_tile_type), intent(inout) :: tile
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

  associate (vegn=>tile%vegn)
  if (new_landuse==LU_RANGE) return ! do nothing for the conversion to rangeland

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

        ! distribute harvested wood between pools
        if (new_landuse==LU_SCND) then
           ! this is harvesting, distribute between 3 different wood pools
           vegn%harv_pool_nitrogen(HARV_POOL_WOOD_FAST) = vegn%harv_pool_nitrogen(HARV_POOL_WOOD_FAST) &
                + (cc%bwood/sp%wood_c2n+cc%bsw/sp%sapwood_c2n)*frac_harvested*(1-frac_wood_wasted)*frac_wood_fast
           vegn%harv_pool_nitrogen(HARV_POOL_WOOD_MED) = vegn%harv_pool_nitrogen(HARV_POOL_WOOD_MED) &
                + (cc%bwood/sp%wood_c2n+cc%bsw/sp%sapwood_c2n)*frac_harvested*(1-frac_wood_wasted)*frac_wood_med
           vegn%harv_pool_nitrogen(HARV_POOL_WOOD_SLOW) = vegn%harv_pool_nitrogen(HARV_POOL_WOOD_SLOW) &
                + (cc%bwood/sp%wood_c2n+cc%bsw/sp%sapwood_c2n)*frac_harvested*(1-frac_wood_wasted)*frac_wood_slow
        else
           ! this is land clearance: everything goes into "cleared" pool
           vegn%harv_pool_nitrogen(HARV_POOL_CLEARED) = vegn%harv_pool_nitrogen(HARV_POOL_CLEARED) &
                + (cc%bwood/sp%wood_c2n+cc%bsw/sp%sapwood_c2n)*frac_harvested*(1-frac_wood_wasted)
        endif

    case default
        call error_mesg('vegn_cut_forest','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
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
  end associate ! vegn
end subroutine vegn_cut_forest_lm3

! ============================================================================
subroutine vegn_graze_pasture_ppa(tile, min_lai_for_grazing, grazing_intensity, max_grazing_height)
  type(land_tile_type), intent(inout) :: tile
  real, intent(in) :: min_lai_for_grazing
  real, intent(in) :: grazing_intensity
  real, intent(in) :: max_grazing_height ! meters. Vegetation taller than this threshold is not grazed

  real    :: littC, littN ! litter from grazing, kg/m2
  real    :: LAI ! leaf area index of *grazed* vegetation. Does not include plants taller
                 ! than max grazing height.
  integer :: i
  real    :: c1, c2, n1, n2 ! for conservation checks

  associate (vegn=>tile%vegn)
  c1 = vegn_tile_carbon(vegn)
  n1 = vegn_tile_nitrogen(vegn)


  LAI = 0.0
  do i = 1,vegn%n_cohorts
     if (vegn%cohorts(i)%height < max_grazing_height) &
             LAI = LAI + vegn%cohorts(i)%lai*vegn%cohorts(i)%nindivs
  enddo
  if (LAI < min_lai_for_grazing) return

  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))
     if (cc%height >= max_grazing_height) continue ! do nothing to vegetation browsers cannot reach.

     vegn%harv_pool(HARV_POOL_PAST) = vegn%harv_pool(HARV_POOL_PAST) + &
          cc%bl*grazing_intensity*(1-grazing_residue)*cc%nindivs
     vegn%harv_pool_nitrogen(HARV_POOL_PAST) = vegn%harv_pool_nitrogen(HARV_POOL_PAST) + &
          cc%leaf_N*grazing_intensity*(1-grazing_residue)*cc%nindivs
     littC = cc%bl     * grazing_intensity*grazing_residue*cc%nindivs
     littN = cc%leaf_N * grazing_intensity*grazing_residue*cc%nindivs
     cc%bl     = cc%bl     * (1-grazing_intensity)
     cc%leaf_N = cc%leaf_N * (1-grazing_intensity)

     ! add carbon to intermediate pools
     select case (soil_carbon_option)
     case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
        ! add litter to intermediate carbon pools
        vegn%fsc_pool_ag = vegn%fsc_pool_ag + littC*sp%fsc_liv
        vegn%ssc_pool_ag = vegn%ssc_pool_ag + littC*(1-sp%fsc_liv)
        ! note that we assume microbial biomass is zero
     case (SOILC_CORPSE, SOILC_CORPSE_N)
        vegn%litter_buff_C(:,LEAF) = vegn%litter_buff_C(:,LEAF) + littC*[sp%fsc_liv,1-sp%fsc_liv,0.0]
        vegn%litter_buff_N(:,LEAF) = vegn%litter_buff_N(:,LEAF) + littN*[sp%fsc_liv,1-sp%fsc_liv,0.0]
     case default
        call error_mesg('vegn_graze_pasture_ppa','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select
     end associate
  enddo

  c2 = vegn_tile_carbon(vegn)
  call check_conservation('vegn_graze_pasture_ppa','carbon',c1,c2,carbon_cons_tol,FATAL)
  n2 = vegn_tile_nitrogen(vegn)
  call check_conservation('vegn_graze_pasture_ppa','nitrogen',n1,n2,nitrogen_cons_tol,FATAL)
  end associate ! vegn
end subroutine vegn_graze_pasture_ppa

! ============================================================================
! NOTE that the PPA harvest would not work properly if applied only once per year, because
! at the end of the year the leaves are down in NH, and therefore harvest amount will be
! unrealistically small. The same is true about the grazing.
subroutine vegn_harvest_crop_ppa(tile)
  type(land_tile_type), intent(inout) :: tile

  real, dimension(N_C_TYPES) :: &
     leaf_litt_C, leaf_litt_N, & ! accumulated leaf litter, kg/m2
     wood_litt_C, wood_litt_N    ! accumulated wood litter, kg/m2
  real, dimension(num_l, N_C_TYPES) :: &
     root_litt_C, root_litt_N    ! accumulated root litter per soil layer, kg/m2
  real :: ndead  ! number of individuals killed in cohort
  integer :: i
  real :: c1, c2, n1, n2 ! for conservation checks

  associate(vegn=>tile%vegn, soil=>tile%soil)
  c1 = vegn_tile_carbon(vegn)
  n1 = vegn_tile_nitrogen(vegn)

  leaf_litt_C=0.0; wood_litt_C=0.0; root_litt_C=0.0
  leaf_litt_N=0.0; wood_litt_N=0.0; root_litt_N=0.0
  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))
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
     vegn%harv_pool_nitrogen(HARV_POOL_CROP) = vegn%harv_pool_nitrogen(HARV_POOL_CROP) + ndead*( &
         cc%leaf_N + cc%seed_N + cc%sapwood_N*agf_bs + cc%stored_N*agf_bs &
         )
     ! subtract C and N harvest from cohort
     cc%bl = 0.0; cc%bseed=0.0; cc%carbon_gain = 0.0; cc%growth_previous_day=0.0
     cc%bsw = cc%bsw*(1-agf_bs); cc%brsw = cc%brsw*(1-agf_bs); cc%nsc=cc%nsc*(1-agf_bs)
     cc%leaf_N = 0.0;  cc%seed_N = 0.0
     cc%sapwood_N=cc%sapwood_N*(1-agf_bs); cc%stored_N = cc%stored_N*(1-agf_bs)

     call kill_plants_ppa(cc,vegn,ndead,0.0, leaf_litt_C, wood_litt_C, root_litt_C, &
                                             leaf_litt_N, wood_litt_N, root_litt_N  )
     end associate
  enddo

  ! add carbon to intermediate pools
  select case (soil_carbon_option)
  case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
     ! add litter to intermediate carbon pools
     vegn%fsc_pool_ag = vegn%fsc_pool_ag + wood_litt_C(C_FAST) + leaf_litt_C(C_FAST)
     vegn%ssc_pool_ag = vegn%ssc_pool_ag + wood_litt_C(C_SLOW) + leaf_litt_C(C_SLOW)
     ! note that we assume microbial biomass is zero
  case (SOILC_CORPSE,SOILC_CORPSE_N)
     vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + wood_litt_C(:)
     vegn%litter_buff_N(:,CWOOD) = vegn%litter_buff_N(:,CWOOD) + wood_litt_N(:)
     vegn%litter_buff_C(:,LEAF)  = vegn%litter_buff_C(:,LEAF)  + leaf_litt_C(:)
     vegn%litter_buff_N(:,LEAF)  = vegn%litter_buff_N(:,LEAF)  + leaf_litt_N(:)
  case default
     call error_mesg('vegn_harvest_crop_ppa','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select

  vegn%fsc_pool_bg = vegn%fsc_pool_bg + sum(root_litt_C(:,C_FAST))+sum(root_litt_C(:,C_MIC))
  vegn%fsn_pool_bg = vegn%fsn_pool_bg + sum(root_litt_N(:,C_FAST))+sum(root_litt_N(:,C_MIC))
  vegn%ssc_pool_bg = vegn%ssc_pool_bg + sum(root_litt_C(:,C_SLOW))
  vegn%ssn_pool_bg = vegn%ssn_pool_bg + sum(root_litt_N(:,C_SLOW))

  call vegn_relayer_cohorts_ppa(vegn)

  c2 = vegn_tile_carbon(vegn)
  call check_conservation('vegn_harvest_crop_ppa','carbon',c1,c2,carbon_cons_tol,FATAL)
  n2 = vegn_tile_nitrogen(vegn)
  call check_conservation('vegn_harvest_crop_ppa','nitrogen',n1,n2,nitrogen_cons_tol,FATAL)
  end associate ! vegn, soil
end subroutine vegn_harvest_crop_ppa

! ============================================================================
subroutine vegn_cut_forest_ppa(tile, new_landuse)
  type(land_tile_type), intent(inout) :: tile
  integer, intent(in) :: new_landuse ! new land use type that gets assigned to
                                     ! the tile after the wood harvesting

  ! ---- local vars
  real :: frac_wood_wasted       ! fraction of wood wasted during transition
  real, dimension(N_C_TYPES) :: &
     leaf_litt_C, leaf_litt_N, & ! accumulated leaf litter, kg/m2
     wood_harv_C, wood_harv_N, & ! accumulated wood harvest, kg/m2
     wood_litt_C, wood_litt_N    ! accumulated wood litter, kg/m2
  real, dimension(num_l, N_C_TYPES) :: &
     root_litt_C, root_litt_N    ! accumulated root litter per soil layer, kg/m2
  real :: ndead  ! number of individuals killed in cohort
  real :: dbh_min
  integer :: i
  real :: c1, c2, n1, n2 ! for conservation checks

  if (new_landuse==LU_RANGE) return ! do nothing for the conversion to rangeland

  associate(vegn=>tile%vegn, soil=>tile%soil)
  c1 = vegn_tile_carbon(vegn)
  n1 = vegn_tile_nitrogen(vegn)

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

  leaf_litt_C=0.0; wood_harv_C=0.0; wood_litt_C=0.0; root_litt_C=0.0
  leaf_litt_N=0.0; wood_harv_N=0.0; wood_litt_N=0.0; root_litt_N=0.0
  do i = 1, vegn%n_cohorts
     associate (cc=>vegn%cohorts(i))
     if (cc%dbh > dbh_min) then
        ! these trees are harvested
        ndead = cc%nindivs
        call kill_plants_ppa(cc,vegn,ndead,0.0, leaf_litt_C, wood_harv_C, root_litt_C, &
                                                leaf_litt_N, wood_harv_N, root_litt_N  )
     else
        ! these trees are too small to be harvested, so a part of them get trampled
        ! and goes to waste, the rest stays
        ndead = cc%nindivs * frac_trampled
        call kill_plants_ppa(cc,vegn,ndead,0.0, leaf_litt_C, wood_harv_C, root_litt_C, &
                                                leaf_litt_N, wood_harv_N, root_litt_N  )
     endif
     end associate
     ! note that below-ground wood all goes to litter, by construction of kill_plants_ppa
  enddo

  ! distribute harvested wood between pools
  if (new_landuse==LU_SCND) then
     ! this is harvesting, distribute between 3 different wood pools
     vegn%harv_pool(HARV_POOL_WOOD_FAST) = vegn%harv_pool(HARV_POOL_WOOD_FAST) &
          + sum(wood_harv_C)*frac_wood_fast*(1-frac_wood_wasted)
     vegn%harv_pool_nitrogen(HARV_POOL_WOOD_FAST) = vegn%harv_pool_nitrogen(HARV_POOL_WOOD_FAST) &
          + sum(wood_harv_N)*frac_wood_fast*(1-frac_wood_wasted)
     vegn%harv_pool(HARV_POOL_WOOD_MED) = vegn%harv_pool(HARV_POOL_WOOD_MED) &
          + sum(wood_harv_C)*frac_wood_med*(1-frac_wood_wasted)
     vegn%harv_pool_nitrogen(HARV_POOL_WOOD_MED) = vegn%harv_pool_nitrogen(HARV_POOL_WOOD_MED) &
          + sum(wood_harv_N)*frac_wood_med*(1-frac_wood_wasted)
     vegn%harv_pool(HARV_POOL_WOOD_SLOW) = vegn%harv_pool(HARV_POOL_WOOD_SLOW) &
          + sum(wood_harv_C)*frac_wood_slow*(1-frac_wood_wasted)
     vegn%harv_pool_nitrogen(HARV_POOL_WOOD_SLOW) = vegn%harv_pool_nitrogen(HARV_POOL_WOOD_SLOW) &
          + sum(wood_harv_N)*frac_wood_slow*(1-frac_wood_wasted)
  else
     ! this is land clearance: everything goes into "cleared" pool
     vegn%harv_pool(HARV_POOL_CLEARED) = vegn%harv_pool(HARV_POOL_CLEARED) &
          + sum(wood_harv_C)*(1-frac_wood_wasted)
     vegn%harv_pool_nitrogen(HARV_POOL_CLEARED) = vegn%harv_pool_nitrogen(HARV_POOL_CLEARED) &
          + sum(wood_harv_N)*(1-frac_wood_wasted)
  endif

  select case (soil_carbon_option)
  case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
     ! add litter to intermediate carbon pools
     vegn%fsc_pool_ag = vegn%fsc_pool_ag + &
         wood_harv_C(C_FAST)*frac_wood_wasted + wood_litt_C(C_FAST) + leaf_litt_C(C_FAST)
     vegn%ssc_pool_ag = vegn%ssc_pool_ag + &
         wood_harv_C(C_SLOW)*frac_wood_wasted + wood_litt_C(C_SLOW) + leaf_litt_C(C_SLOW)
     ! note that we assume microbial biomass is zero
  case (SOILC_CORPSE,SOILC_CORPSE_N)
     vegn%litter_buff_C(:,CWOOD) = vegn%litter_buff_C(:,CWOOD) + &
        wood_litt_C(:) + wood_harv_C(:)*frac_wood_wasted
     vegn%litter_buff_N(:,CWOOD) = vegn%litter_buff_N(:,CWOOD) + &
        wood_litt_N(:) + wood_harv_N(:)*frac_wood_wasted
     vegn%litter_buff_C(:,LEAF) = vegn%litter_buff_C(:,LEAF) + leaf_litt_C(:)
     vegn%litter_buff_N(:,LEAF) = vegn%litter_buff_N(:,LEAF) + leaf_litt_N(:)
  case default
     call error_mesg('vegn_cut_forest_ppa','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select
  vegn%fsc_pool_bg = vegn%fsc_pool_bg + sum(root_litt_C(:,C_FAST))+sum(root_litt_C(:,C_MIC))
  vegn%fsn_pool_bg = vegn%fsn_pool_bg + sum(root_litt_N(:,C_FAST))+sum(root_litt_N(:,C_MIC))
  vegn%ssc_pool_bg = vegn%ssc_pool_bg + sum(root_litt_C(:,C_SLOW))
  vegn%ssn_pool_bg = vegn%ssn_pool_bg + sum(root_litt_N(:,C_SLOW))

  call vegn_relayer_cohorts_ppa(vegn)

  c2 = vegn_tile_carbon(vegn)
  call check_conservation('vegn_cut_forest_ppa','carbon',c1,c2,carbon_cons_tol,FATAL)
  n2 = vegn_tile_nitrogen(vegn)
  call check_conservation('vegn_cut_forest_ppa','nitrogen',n1,n2,nitrogen_cons_tol,FATAL)
  end associate ! vegn,soil
end subroutine vegn_cut_forest_ppa


subroutine vegn_plant_crop_ppa(tile)
  type(land_tile_type), intent(inout) :: tile

  ! list of pools we borrow seeds from, highest priority first
  integer, parameter :: seed_source_pools(4) = &
     [ HARV_POOL_CROP, HARV_POOL_PAST, HARV_POOL_CLEARED, HARV_POOL_WOOD_FAST ]

  integer :: i, p
  real, dimension(0:nspecies-1) :: seedC, seedN ! seed biomass, kg/m2
  real :: deltaC, deltaN ! amount we borrow from harvest pools, kg/m2
  real :: c1,c2,n1,n2 ! for conservation checking

  associate (vegn=>tile%vegn, soil=>tile%soil)
  c1 = vegn_tile_carbon(vegn) + soil_tile_carbon(soil)
  n1 = vegn_tile_nitrogen(vegn) + soil_tile_nitrogen(soil)

  if (crop_species_idx<0) then
     call error_mesg('vegn_plant_crop_ppa','crop species "'//trim(crop_species)//'" not found in the list of species',FATAL)
  endif
  
  ! borrow biomass (crop_seed_density) from harvest pools, in order of preference
  seedC(:) = 0.0; seedN(:) = 0.0
  do i = 1, size(seed_source_pools)
     p = seed_source_pools(i)
     deltaC = max(crop_seed_density-seedC(crop_species_idx),0.0)
     deltaC = min(deltaC,vegn%harv_pool(p))
     if (vegn%harv_pool(p)>0) then
        deltaN = deltaC*vegn%harv_pool_nitrogen(p)/vegn%harv_pool(p)
     else
        deltaN = 0.0
     endif
     vegn%harv_pool(p)          = vegn%harv_pool(p)          - deltaC
     vegn%harv_pool_nitrogen(p) = vegn%harv_pool_nitrogen(p) - deltaN
     seedC(crop_species_idx) = seedC(crop_species_idx) + deltaC
     seedN(crop_species_idx) = seedN(crop_species_idx) + deltaN
     if (seedC(crop_species_idx) >= crop_seed_density) exit ! from loop
  enddo
  call add_seedlings_ppa(vegn,soil,seedC,seedN)

  c2 = vegn_tile_carbon(vegn) + soil_tile_carbon(soil)
  n2 = vegn_tile_nitrogen(vegn) + soil_tile_nitrogen(soil)
  call check_conservation('vegn_plant_crop_ppa','carbon',   c1,c2,carbon_cons_tol,   FATAL)
  call check_conservation('vegn_plant_crop_ppa','nitrogen', n1,n2,nitrogen_cons_tol, FATAL)
  end associate ! vegn,soil
end subroutine vegn_plant_crop_ppa

end module
