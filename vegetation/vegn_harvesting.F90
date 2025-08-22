module vegn_harvesting_mod

#include "../shared/debug.inc"

use constants_mod, only : tfreeze
use fms_mod, only : string, error_mesg, FATAL, NOTE, WARNING, &
     mpp_pe, check_nml_error, stdlog, mpp_root_pe, lowercase
use mpp_mod, only : mpp_sum, input_nml_file
use mpp_domains_mod, only : mpp_global_sum, BITWISE_EXACT_SUM, mpp_pass_UG_to_SG
use diag_manager_mod, only : register_static_field, send_data

use land_constants_mod, only : N_C_TYPES, C_FAST, C_SLOW, C_MIC, LITT_LEAF, LITT_CWOOD, &
     seconds_per_year
use land_io_mod, only : read_field
use land_debug_mod, only : string_from_time, land_error_message, check_conservation, &
     do_check_conservation, carbon_cons_tol, nitrogen_cons_tol, check_var_range
use land_utils_mod, only : check_conservation_1, check_conservation_2
use land_data_mod, only : log_version, lnd
use vegn_data_mod, only : do_ppa, is_cropland, &
     N_LU_TYPES, LU_PAST, LU_RAINF, LU_IRRIG, LU_NTRL, LU_SCND, LU_RANGE, &
     HARV_POOL_PAST, HARV_POOL_CROP, HARV_POOL_CLEARED, HARV_POOL_WOOD_FAST, &
     HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW, PT_C3, PT_C4, LEAF_OFF, &
     nspecies, spdata, agf_bs, NO_DATE, NO_CROP, &
     IRRIGATED_MAIZE, IRRIGATED_SOYBEAN, IRRIGATED_RICE, IRRIGATED_SPRING_WHEAT, IRRIGATED_WINTER_WHEAT, &
     RAINFED_MAIZE, RAINFED_SOYBEAN, RAINFED_RICE, RAINFED_SPRING_WHEAT, RAINFED_WINTER_WHEAT, &
     crop_name, landuse_name
use land_tile_mod, only : land_tile_type, land_tile_enum_type, land_tile_map, &
     first_elmt, loop_over_tiles, land_tile_nitrogen, land_tile_carbon
use soil_tile_mod, only : num_l, dz
use vegn_tile_mod, only : vegn_relayer_cohorts_ppa, vegn_mergecohorts_ppa, &
     vegn_tile_LAI, vegn_tile_type
use vegn_cohort_mod, only : update_biomass_pools, cohort_root_litter_profile
use vegn_util_mod, only : kill_plants_ppa, add_seedlings_ppa
use soil_BGC_SIMPLE_type_mod, only: soil_BGC_SIMPLE_t
use soil_BGC_CORPSE_type_mod, only: soil_BGC_CORPSE_t, do_CORPSE_nitrogen => do_nitrogen
use vegn_crop_mod, only: vegn_crop_init, compute_crop_calendars, vegn_crop_end, save_crop_restart
use fms2_io_mod, only: close_file, FmsNetcdfFile_t, open_file

implicit none
private

! ==== public interface ======================================================
public :: vegn_harvesting_init
public :: vegn_harvesting_end

public :: vegn_harvesting
public :: crop_seed_transport

public :: vegn_cut_forest
public :: save_harvesting_restart
! ==== end of public interface ===============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_harvesting_mod'
#include "../shared/version_variable.inc"
real, parameter :: ONETHIRD = 1.0/3.0
integer, parameter :: & ! grazing frequency options
  GRAZING_DAILY  = 1, &
  GRAZING_ANNUAL = 2
integer, parameter :: & ! crop scedule options
  CROP_SCHEDULE_LM3        = 1, &
  CROP_SCHEDULE_PRESCRIBED = 2, &
  CROP_SCHEDULE_COMPUTED   = 3
integer, parameter :: & ! crop type distribution options
  CROP_DISTR_LM3           = 1, &
  CROP_DISTR_LUH2_DOMINANT = 2, &
  CROP_DISTR_MIRCA2000     = 3

! TODO: possibly move all definition of cpw,clw,csw in one place
real, parameter :: &
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

! ==== module data ===========================================================

! ---- namelist variables ----------------------------------------------------
logical, public, protected  :: do_harvesting = .TRUE.  ! if true, then planting and harvesting of crops and pastures is done

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

character(16) :: crop_schedule = 'lm3' ! or 'prescribed' or 'computed'
character(1024) :: crop_schedule_file  = '' ! input data set of crop planting and harvesting dates
character(32) :: crop_distribution = 'lm3' ! or 'MIRCA2000' or 'luh2-dominant'
character(1024) :: luh2_state_file  = '' ! input data set of LU states (for C3/C4 crop distribution in case of 'luh2-dominant' crop species distribution)
character(32) :: c3_crop_species  = '' ! name of the species used for C3 crops
character(32) :: c4_crop_species  = '' ! name of the species used for C4 crops
character(32) :: maize_crop_species = ''
character(32) :: spring_wheat_crop_species = ''
character(32) :: winter_wheat_crop_species = ''
character(32) :: rice_crop_species  = ''
character(32) :: soybean_crop_species = ''
real :: crop_seed_density      = 0.1   ! biomass of seeds left after crop harvesting, kg/m2
real :: crop_seed_c2n          = 30    ! crop seed C:N ratio, used to calculate N demand for crop seed transport
logical, public, protected :: allow_weeds_on_crops = .FALSE. ! if TRUE, seeds transported
        ! from outside of cropland can start growing on croplands; if FALSE they are not
        ! allowed to germinate.
logical :: clear_crop_before_planting = .TRUE. ! if TRUE, all vegetation is removed from
        ! croplands right before planting; otherwise planting adds crops to existing
        ! (presumably small) vegetation
logical, public, protected :: clear_all_on_conversion_to_crop = .TRUE. ! if TRUE
        ! then all vegetation is removed in transition any -> crop; otherwise for some LU
        ! types (e.g. pastures) vegetation remains unchanged, resulting in crops contaminated
        ! by other species (including woody) that happened to grow there.
logical :: transport_crop_seeds = .TRUE. ! if true, seeds are transported horizontally
        ! to satisfy the demand

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
     crop_schedule, crop_schedule_file, &
     crop_distribution, luh2_state_file, &
     c3_crop_species, c4_crop_species, maize_crop_species, spring_wheat_crop_species, winter_wheat_crop_species, rice_crop_species, soybean_crop_species, &
     crop_seed_density, allow_weeds_on_crops, clear_crop_before_planting, clear_all_on_conversion_to_crop, &
     transport_crop_seeds, crop_seed_c2n

integer :: grazing_freq = -1 ! indicator of grazing frequency (GRAZING_ANNUAL or GRAZING_DAILY)
integer :: crop_schedule_option = -1 ! selected planting/harvesting schedule option
integer :: crop_distribution_option = -1
real, allocatable :: crop_planting_day(:) ! day of year when planting is done
real, allocatable :: crop_harvest_day(:)  ! day of year when harvesting is done
integer :: c3_crop_idx = -1, c4_crop_idx = -1, maize_crop_idx = -1 ! index of crop species
integer :: spring_wheat_crop_idx = -1, winter_wheat_crop_idx = -1, rice_crop_idx = -1, soybean_crop_idx = -1

integer :: id_crop_planting_day, id_crop_harvest_day

real :: tot_area_land ! global land area, m2 (for normalization in conservation checks)

contains ! ###################################################################

! ============================================================================
subroutine vegn_harvesting_init(id_ug)
  integer, intent(in) :: id_ug ! id of the unstructured grid diagnostic axis

  integer :: ierr, io, i
  logical :: used
  type(FmsNetcdfFile_t) :: fileobj
  logical :: exists

  call log_version(version, module_name, __FILE__)

  read (input_nml_file, nml=harvesting_nml, iostat=io)
  ierr = check_nml_error(io, 'harvesting_nml')
  if (mpp_pe() == mpp_root_pe()) then
     write(stdlog(), nml=harvesting_nml)
  endif

  if (frac_wood_fast+frac_wood_med+frac_wood_slow/=1.0) then
     call error_mesg('vegn_harvesting_init', &
          'sum of frac_wood_fast, frac_wood_med, and frac_wood_slow must be 1.0',&
          FATAL)
  endif
  ! parse the grazing frequency parameter
  select case(lowercase(grazing_frequency))
  case('annual')
     grazing_freq = GRAZING_ANNUAL
  case('daily')
     grazing_freq = GRAZING_DAILY
     ! scale grazing intensity for daily frequency
     grazing_intensity_past  = grazing_intensity_past/365.0
     grazing_intensity_range = grazing_intensity_range/365.0
  case default
     call error_mesg('vegn_harvesting_init','grazing_frequency must be "annual" or "daily"',FATAL)
  end select

  select case(trim(lowercase(crop_schedule)))
  case('lm3')
     allocate (crop_planting_day(lnd%ls:lnd%le),crop_harvest_day(lnd%ls:lnd%le))
     crop_schedule_option = CROP_SCHEDULE_LM3
     crop_planting_day = 1.0
     crop_harvest_day  = 1.0
  case('prescribed')
     crop_schedule_option = CROP_SCHEDULE_PRESCRIBED
     ! read input data
     exists = open_file(fileobj, crop_schedule_file, "read")
     if (.not. exists) then
        call error_mesg("vegn_harvesting_init", trim(crop_schedule_file)//" does not exist.", FATAL)
     endif
     allocate (crop_planting_day(lnd%ls:lnd%le),crop_harvest_day(lnd%ls:lnd%le))
     call read_field( fileobj, 'plantingdy', crop_planting_day, interp='nearest' )
     call read_field( fileobj, 'harvestdy',  crop_harvest_day,  interp='nearest' )
     call close_file(fileobj)
  case('computed')
     crop_schedule_option = CROP_SCHEDULE_COMPUTED
  case default
     call error_mesg('vegn_harvesting_init','crop_schedule must be "lm3", "prescribed", or "computed"',FATAL)
  end select

  ! initialize crop C3/C4 distribution option
  if(crop_schedule_option == CROP_SCHEDULE_COMPUTED) then
     crop_distribution_option = CROP_DISTR_MIRCA2000
  else
     id_crop_harvest_day = register_static_field ( 'vegn', 'crop_harvest_day', (/id_ug/), &
            'day of year when crops are harvested', 'day', missing_value = -1.0 )
     id_crop_planting_day = register_static_field ( 'vegn', 'crop_planting_day', (/id_ug/), &
            'day of year when crops are planted', 'day', missing_value = -1.0 )
     if ( id_crop_harvest_day > 0 )  used = send_data ( id_crop_harvest_day, crop_harvest_day, lnd%time )
     if ( id_crop_planting_day > 0 ) used = send_data ( id_crop_planting_day, crop_planting_day, lnd%time )

     select case (trim(lowercase(crop_distribution)))
     case ('lm3')
        crop_distribution_option = CROP_DISTR_LM3
     case ('luh2-dominant')
        crop_distribution_option = CROP_DISTR_LUH2_DOMINANT
        ! open input state file
        ! add variables to varset
        call error_mesg('vegn_harvesting_init','crop_distribution "luh2-dominant" is not implemented yet',FATAL)
     case default
        call error_mesg('vegn_harvesting_init','crop_distribution must be "lm3" or "luh2-dominant"',FATAL)
     end select
  endif

  ! find crop species; it is OK if it is not found, perhaps we do not need it -- e.g.
  ! we are running potential vegetation, or in LM3 mode where it is not used (yet)
  ! For now, use single species for crop everywhere. This will obviously need to be changed
  ! when we switch to more sophisticated treatment of agriculture.
  do i = 0, nspecies-1
     if (trim(spdata(i)%name)==trim(c3_crop_species)) c3_crop_idx = i
     if (trim(spdata(i)%name)==trim(c4_crop_species)) c4_crop_idx = i
     if (trim(spdata(i)%name)==trim(maize_crop_species))   maize_crop_idx   = i
     if (trim(spdata(i)%name)==trim(spring_wheat_crop_species)) spring_wheat_crop_idx = i
     if (trim(spdata(i)%name)==trim(winter_wheat_crop_species)) winter_wheat_crop_idx = i
     if (trim(spdata(i)%name)==trim(rice_crop_species))    rice_crop_idx    = i
     if (trim(spdata(i)%name)==trim(soybean_crop_species)) soybean_crop_idx = i
  enddo
  if (c3_crop_idx<0) then
     call error_mesg('vegn_harvesting_init','C3 crop species "'//trim(c3_crop_species)//'" not found in the list of species',NOTE)
  else
     call error_mesg('vegn_harvesting_init','C3 crop species "'//trim(c3_crop_species)//'" is #'//string(c3_crop_idx)//&
                     ' in the list of species', NOTE)
  endif
  if (c4_crop_idx<0) then
     call error_mesg('vegn_harvesting_init','C4 crop species "'//trim(c4_crop_species)//'" not found in the list of species',NOTE)
  else
     call error_mesg('vegn_harvesting_init','C4 crop species "'//trim(c4_crop_species)//'" is #'//string(c4_crop_idx)//&
                     ' in the list of species', NOTE)
  endif
  if (maize_crop_idx<0) then
     call error_mesg('vegn_harvesting_init','maize crop species "'//trim(maize_crop_species)//'" not found in the list of species',NOTE)
  else
     call error_mesg('vegn_harvesting_init','maize crop species "'//trim(maize_crop_species)//'" is #'//string(maize_crop_idx)//&
                     ' in the list of species', NOTE)
  endif
  if (spring_wheat_crop_idx<0) then
     call error_mesg('vegn_harvesting_init','spring wheat crop species "'//trim(spring_wheat_crop_species)//'" not found in the list of species',NOTE)
  else
     call error_mesg('vegn_harvesting_init','spring wheat crop species "'//trim(spring_wheat_crop_species)//'" is #'//string(spring_wheat_crop_idx)//&
                     ' in the list of species', NOTE)
  endif
  if (winter_wheat_crop_idx<0) then
     call error_mesg('vegn_harvesting_init','winter wheat crop species "'//trim(winter_wheat_crop_species)//'" not found in the list of species',NOTE)
  else
     call error_mesg('vegn_harvesting_init','winter wheat crop species "'//trim(winter_wheat_crop_species)//'" is #'//string(winter_wheat_crop_idx)//&
                     ' in the list of species', NOTE)
  endif
  if (rice_crop_idx<0) then
     call error_mesg('vegn_harvesting_init','rice crop species "'//trim(rice_crop_species)//'" not found in the list of species',NOTE)
  else
     call error_mesg('vegn_harvesting_init','rice crop species "'//trim(rice_crop_species)//'" is #'//string(rice_crop_idx)//&
                     ' in the list of species', NOTE)
  endif
  if (soybean_crop_idx<0) then
     call error_mesg('vegn_harvesting_init','soybean crop species "'//trim(soybean_crop_species)//'" not found in the list of species',NOTE)
  else
     call error_mesg('vegn_harvesting_init','soybean crop species "'//trim(soybean_crop_species)//'" is #'//string(soybean_crop_idx)//&
                     ' in the list of species', NOTE)
  endif

  ! calculate total land and soil areas
  tot_area_land = sum(lnd%ug_area)
  call mpp_sum(tot_area_land)
  if(crop_schedule_option == CROP_SCHEDULE_COMPUTED) call vegn_crop_init( id_ug )
end subroutine vegn_harvesting_init

! ============================================================================
subroutine vegn_harvesting_end
   if (allocated(crop_harvest_day))  deallocate(crop_harvest_day)
   if (allocated(crop_planting_day)) deallocate(crop_planting_day)
   if(crop_schedule_option == CROP_SCHEDULE_COMPUTED) call vegn_crop_end()
end subroutine vegn_harvesting_end

! ============================================================================
! harvest vegetation in a tile
subroutine vegn_harvesting(tile, end_of_year, end_of_month, end_of_day, day_of_year, l)
  type(land_tile_type), intent(inout) :: tile
  logical, intent(in) :: end_of_year, end_of_month, end_of_day ! indicators of respective period boundaries
  integer, intent(in) :: day_of_year ! current day of year
  integer, intent(in) :: l ! index of current grid cell in unstructured grid
  logical :: a_crop_is_active
  integer :: chosen_crp

  if (.not.do_harvesting) return ! do nothing if no harvesting requested
  if (crop_schedule_option == CROP_SCHEDULE_COMPUTED) then
     call compute_crop_calendars(tile%vegn, tile%diag, L)
  endif

  associate(vegn=>tile%vegn)
  select case(vegn%landuse)
  case(LU_PAST)  ! pasture
     if ((end_of_day  .and. grazing_freq==GRAZING_DAILY).or. &
         (end_of_year .and. grazing_freq==GRAZING_ANNUAL)) then
        call vegn_graze_pasture (tile)
     endif
  case(LU_RANGE)  ! rangeland
     if ((end_of_day  .and. grazing_freq==GRAZING_DAILY).or. &
         (end_of_year .and. grazing_freq==GRAZING_ANNUAL)) then
        call vegn_graze_rangeland (tile)
     endif
  case(LU_RAINF, LU_IRRIG)  ! crop
     select case(crop_schedule_option)
     case (CROP_SCHEDULE_LM3)
        if (end_of_year) then
           call vegn_harvest_cropland (tile,'AA')
           call vegn_plant_crop (tile)
        endif
     case (CROP_SCHEDULE_PRESCRIBED)
        if (end_of_day.AND.day_of_year==nint(crop_harvest_day(l))) then
           call vegn_harvest_cropland (tile,'BB')
        endif
        if (end_of_day.AND.day_of_year==nint(crop_planting_day(l))) then
           call vegn_plant_crop (tile)
        endif
     case (CROP_SCHEDULE_COMPUTED)
        if(end_of_year) then
           if(vegn%Crop%grass_is_active) then
              call vegn_harvest_cropland (tile,'CC')
              vegn%Crop%grass_is_active = .FALSE.
           endif
           if(vegn%Crop%chosen_calendars(1,1) == NO_DATE .AND. vegn%Crop%chosen_calendars(1,2) == NO_DATE) then
              ! There is no crop to plant so plant grass on the LM3 schedule.
              call vegn_plant_crop (tile, chosen_crop=NO_CROP)
              vegn%Crop%grass_is_active = .TRUE.
           endif
        endif
        if(end_of_day .AND. (day_of_year==vegn%Crop%chosen_calendars(1,1) .or. day_of_year==vegn%Crop%chosen_calendars(1,2))) then
           ! Today is either the main or second season planting date
           if(vegn%Crop%grass_is_active) then
              ! harvest grass before planting crop
              call vegn_harvest_cropland (tile,'DD')
              vegn%Crop%grass_is_active = .FALSE.
           endif
           a_crop_is_active = vegn%Crop%chosen_crop_is_active(1) .or. vegn%Crop%chosen_crop_is_active(2)
           if(day_of_year==vegn%Crop%chosen_calendars(1,1)) then
             ! Today is the 1st crop's planting date
             if(.not.a_crop_is_active) then
               chosen_crp = vegn%Crop%chosen_crop(1)
               call vegn_plant_crop (tile, chosen_crop=chosen_crp)
               vegn%Crop%chosen_crop_is_active(1) = .TRUE.
             endif
           else if(day_of_year==vegn%Crop%chosen_calendars(1,2)) then
             ! Today is the 2nd crop's planting date
             if(.not.a_crop_is_active) then
               chosen_crp = vegn%Crop%chosen_crop(2)
               call vegn_plant_crop (tile, chosen_crop=chosen_crp)
               vegn%Crop%chosen_crop_is_active(2) = .TRUE.
             endif
           endif
        endif
        if(end_of_day .AND. (day_of_year==vegn%Crop%chosen_calendars(2,1) .or. day_of_year==vegn%Crop%chosen_calendars(2,2))) then
           ! Today is either the main or second season harvest date. A crop is in the ground and ready to harvest.
           if(day_of_year==vegn%Crop%chosen_calendars(2,1)) then
             ! Today is the 1st crop's harvest date
             if(vegn%Crop%chosen_crop_is_active(1)) then
               chosen_crp = vegn%Crop%chosen_crop(1) ! debug_pjp
               call vegn_harvest_cropland (tile,'EE '//trim(crop_name(chosen_crp))) ! debug_pjp
               vegn%Crop%chosen_crop_is_active(1) = .FALSE.
             endif
           else if(day_of_year==vegn%Crop%chosen_calendars(2,2)) then
             ! Today is the 2nd crop's harvest date
             if(vegn%Crop%chosen_crop_is_active(2)) then
               call vegn_harvest_cropland (tile,'FF')
               vegn%Crop%chosen_crop_is_active(2) = .FALSE.
             endif
           endif
        endif
     end select ! crop_schedule_option
  end select ! vegn%landuse
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
subroutine vegn_harvest_cropland(tile, ctag) ! debug_pjp
  type(land_tile_type), intent(inout) :: tile
  character(len=*), intent(in) :: ctag ! debug_pjp

  if (do_ppa) then
     call vegn_harvest_crop_ppa(tile, ctag) ! debug_pjp
  else
     call vegn_harvest_crop_lm3(tile)
  endif
end subroutine vegn_harvest_cropland

! ============================================================================
subroutine vegn_plant_crop(tile, chosen_crop)
  type(land_tile_type), intent(inout) :: tile
  integer, optional, intent(in) :: chosen_crop

  if (do_ppa) then
     call vegn_plant_crop_ppa(tile, chosen_crop)
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
  integer :: i,k
  real :: carbon_lost
  real :: delta_leaf, delta_root, delta_wood
  real,dimension(N_C_TYPES) :: leaflitter_C,woodlitter_C,leaflitter_N,woodlitter_N
  real :: bglitter_C(num_l,N_C_TYPES) ! below-ground (root) C litter, by layer
  real :: bglitter_N(num_l,N_C_TYPES) ! below-ground (root) N litter, by layer
  real :: profile(num_l) ! normalized root litter profile: sum(profile) == 1.0
  real :: wood_n2c

  associate(vegn=>tile%vegn,soil=>tile%soil)
  ! do nothing if the LAI is less that the lower grazing limit
  if ( vegn_tile_LAI(vegn) <= min_lai_for_grazing ) return

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1,vegn%n_cohorts
     associate (cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))

     ! This makes sure biomass pools are correct before calculating changes
     ! in leaf biomass and such, just in case it was not called before
     call update_biomass_pools(cc)

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
     vegn%harv_pool_C(HARV_POOL_PAST) = vegn%harv_pool_C(HARV_POOL_PAST) + &
          carbon_lost*(1-grazing_residue)
     cc%bliving = cc%bliving - carbon_lost

     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc)

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
     select type(soilc=>tile%soilc)
     class is (soil_BGC_SIMPLE_t)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + grazing_residue*( &
             sp%fsc_liv*(balive0-balive1)+sp%fsc_wood*(bdead0-bdead1))
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + grazing_residue*( &
             (1-sp%fsc_liv)*(balive0-balive1)+ (1-sp%fsc_wood)*(bdead0-bdead1))
     class is (soil_BGC_CORPSE_t)
        if(blv0 < blv1) then ! Some biomass was re-absorbed due to N limitation. Reduce litter.
           delta_leaf=bleaf0-bleaf1   - (blv1-blv0)*(bleaf0-bleaf1)/(bleaf0+bfroot0+bdead0-bleaf1-bfroot1-bdead1)
           delta_root=bfroot0-bfroot1 - (blv1-blv0)*(bfroot0-bfroot1)/(bleaf0+bfroot0+bdead0-bleaf1-bfroot1-bdead1)
           delta_wood=bdead0-bdead1   - (blv1-blv0)*(bdead0-bdead1)/(bleaf0+bfroot0+bdead0-bleaf1-bfroot1-bdead1)
        else   ! Virtual leaves decreased. Send to leaf litter
           delta_leaf=bleaf0+blv0-bleaf1-blv1
           delta_root=bfroot0-bfroot1
           delta_wood=bdead0-bdead1
        endif

        leaflitter_C = [ delta_leaf*sp%fsc_liv,  delta_leaf*(1-sp%fsc_liv),  0.0 ] * grazing_residue
        woodlitter_C = [ delta_wood*sp%fsc_wood, delta_wood*(1-sp%fsc_wood), 0.0 ] * grazing_residue * agf_bs
        call cohort_root_litter_profile(cc,dz,profile)
        do k = 1,num_l
           bglitter_C(k,:) = profile(k) * grazing_residue * &
               [      sp%fsc_froot *delta_root + (1-agf_bs)*     sp%fsc_wood *delta_wood, &
                 (1.0-sp%fsc_froot)*delta_root + (1-agf_bs)*(1.0-sp%fsc_wood)*delta_wood, &
                 0.0  ]
        enddo
        ! We are not removing belowground portion of what was grazed, so that needs to be clawed back from harvest pool
        vegn%harv_pool_C(HARV_POOL_PAST) = vegn%harv_pool_C(HARV_POOL_PAST) - (1.0-grazing_residue)*(delta_root+(1-agf_bs)*delta_wood)

        if(do_CORPSE_nitrogen) then
           leaflitter_N=leaflitter_C/sp%leaf_live_c2n
           woodlitter_N=woodlitter_C/sp%leaf_live_c2n
           do k = 1,num_l
              bglitter_N(k,:) = profile(k) * grazing_residue * &
                   [    sp%fsc_froot *delta_root/sp%froot_live_c2n + (1-agf_bs)*   sp%fsc_wood *delta_wood*wood_n2c, &
                     (1-sp%fsc_froot)*delta_root/sp%froot_live_c2n + (1-agf_bs)*(1-sp%fsc_wood)*delta_wood*wood_n2c, &
                     0.0  ]
           enddo
           cc%stored_N = cc%stored_N - delta_leaf/sp%leaf_live_c2n - delta_wood*wood_n2c - delta_root/sp%froot_live_c2n
           vegn%harv_pool_N(HARV_POOL_PAST) = vegn%harv_pool_N(HARV_POOL_PAST) + &
                delta_leaf/sp%leaf_live_c2n*(1-grazing_residue) + delta_wood*agf_bs*wood_n2c*(1-grazing_residue)
        else
           leaflitter_N = 0.0
           woodlitter_N = 0.0
           bglitter_N   = 0.0
        endif

       if (grazing_freq==GRAZING_DAILY) then
          ! Put carbon directly in soil pools
          call soilc%add_soil_matter(vegn, &
             leaf_litter_C=leaflitter_C, leaf_litter_N=leaflitter_N, &
             wood_litter_C=woodlitter_C, wood_litter_N=woodlitter_N, &
             root_litter_C=bglitter_C,   root_litter_N=bglitter_N    )
!           call add_litter(soilc%litter_corpse(LITT_LEAF),leaflitter_C,leaflitter_N)
!           call add_litter(soilc%litter_corpse(LITT_CWOOD),woodlitter_C,woodlitter_N)
!           call soilc%add_root_litter(vegn,bglitter_C,bglitter_N)
       else
          vegn%litter_buff_C(:,LITT_LEAF) = vegn%litter_buff_C(:,LITT_LEAF) + &
               [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*(delta_leaf)*grazing_residue
          vegn%litter_buff_C(:,LITT_CWOOD) = vegn%litter_buff_C(:,LITT_CWOOD) + &
               [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*agf_bs*(delta_wood)*grazing_residue

          vegn%fsc_pool_bg = vegn%fsc_pool_bg + sum(bglitter_C(:,C_FAST))
          vegn%ssc_pool_bg = vegn%ssc_pool_bg + sum(bglitter_C(:,C_SLOW))


          vegn%litter_buff_N(:,LITT_LEAF) = vegn%litter_buff_N(:,LITT_LEAF) + &
             [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*(delta_leaf)*grazing_residue/sp%leaf_live_c2n
          vegn%litter_buff_N(:,LITT_CWOOD) = vegn%litter_buff_N(:,LITT_CWOOD) + &
             [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*agf_bs*(delta_wood)*grazing_residue/sp%wood_c2n

          vegn%fsn_pool_bg = vegn%fsn_pool_bg + sum(bglitter_N(:,C_FAST))
          vegn%ssn_pool_bg = vegn%ssn_pool_bg + sum(bglitter_N(:,C_SLOW))
       endif
     class default
        call error_mesg('vegn_graze_pasture_lm3','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select
     end associate
  enddo
  end associate ! vegn, soil
end subroutine vegn_graze_pasture_lm3

! ================================================================================
subroutine vegn_harvest_crop_lm3(tile)
  type(land_tile_type), intent(inout) :: tile

  ! ---- local vars
  real :: fraction_harvested    ! fraction of biomass harvested this time
  real :: bdead, balive, btotal ! combined biomass pools
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
  btotal = balive+bdead

  ! calculate harvested fraction: cut everything down to seed level
  fraction_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0)

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))
     ! use for harvest only above-ground living biomass and waste the correspondent below living and wood
     vegn%harv_pool_C(HARV_POOL_CROP) = vegn%harv_pool_C(HARV_POOL_CROP) + &
          cc%bliving*(cc%Pl + cc%Psw*agf_bs)*fraction_harvested
     select type(soilc => tile%soilc)
     class is (soil_BGC_SIMPLE_t)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + fraction_harvested*(sp%fsc_liv*cc%bliving*cc%Pr + &
             sp%fsc_wood*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)))
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + fraction_harvested*((1-sp%fsc_liv)*cc%bliving*cc%Pr + &
             (1-sp%fsc_wood)*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)))
     class is (soil_BGC_CORPSE_t)
        vegn%litter_buff_C(:,LITT_CWOOD) = vegn%litter_buff_C(:,LITT_CWOOD) + &
               [sp%fsc_wood, 1-sp%fsc_wood, 0.0] * fraction_harvested*agf_bs*cc%bwood

        vegn%fsc_pool_bg = vegn%fsc_pool_bg + fraction_harvested*(&
               sp%fsc_froot*cc%bliving*cc%Pr + &
               (1-agf_bs)*sp%fsc_wood*(cc%bwood + cc%bliving*cc%Psw))
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + fraction_harvested*(&
               (1-sp%fsc_froot)*cc%bliving*cc%Pr + &
               (1-agf_bs)*(1-sp%fsc_wood)*(cc%bwood + cc%bliving*cc%Psw))

        if (do_CORPSE_nitrogen) then
           vegn%litter_buff_N(:,LITT_CWOOD) = vegn%litter_buff_N(:,LITT_CWOOD) + &
               [sp%fsc_wood, 1-sp%fsc_wood, 0.0] * fraction_harvested*agf_bs*cc%wood_N

           vegn%fsn_pool_bg = vegn%fsn_pool_bg + fraction_harvested*(&
                   sp%fsc_froot*cc%root_N + &
                   (1-agf_bs)*(sp%fsc_wood*cc%wood_N + sp%fsc_liv*cc%sapwood_N))
           vegn%ssn_pool_bg = vegn%ssn_pool_bg + fraction_harvested*(&
                   (1-sp%fsc_froot)*cc%root_N + &
                   (1-agf_bs)*(cc%wood_N*(1-sp%fsc_wood) + cc%sapwood_N*(1-sp%fsc_liv)))

           vegn%harv_pool_N(HARV_POOL_CROP) = vegn%harv_pool_N(HARV_POOL_CROP) + &
                (cc%bliving*cc%Pl/sp%leaf_live_c2n + agf_bs*cc%sapwood_N)*fraction_harvested

           ! Make sure stored nitrogen loss is sensible when leaves are off:
           ! Amounts harvested should be determined by potential leaf tissues, but this N would still be in stored pool so needs to be subtracted
           if (cc%status == LEAF_OFF) then
              vegn%harv_pool_N(HARV_POOL_CROP) = vegn%harv_pool_N(HARV_POOL_CROP) + &
                 (max(cc%stored_N,0.0)-(cc%bliving*cc%Pl/sp%leaf_live_c2n + cc%bliving*cc%Pr/sp%froot_live_c2n))*fraction_harvested
              cc%stored_N = cc%stored_N - &
                 (max(cc%stored_N,0.0)-(cc%bliving*cc%Pl/sp%leaf_live_c2n + cc%bliving*cc%Pr/sp%froot_live_c2n))*fraction_harvested
           else
              ! In this case stored N can just be lost because it does not contain N from harvested potential leaves and roots
              vegn%harv_pool_N(HARV_POOL_CROP) = vegn%harv_pool_N(HARV_POOL_CROP) + max(cc%stored_N,0.0)*fraction_harvested
              cc%stored_N = max(cc%stored_N,0.0)*(1-fraction_harvested)
           endif
           ! Subtract N lost from N pools. Leaf and root can get subtracted from stored and things will be rebalanced later.
           ! Some stored N also was lost
           cc%stored_N = cc%stored_N -  &
                      cc%bliving*fraction_harvested*(cc%Pl/sp%leaf_live_c2n + cc%Pr/sp%froot_live_c2n)
           cc%sapwood_N = cc%sapwood_N*(1-fraction_harvested)
           cc%wood_N = cc%wood_N*(1-fraction_harvested)
       endif

     class default
        call error_mesg('vegn_harvest_crop_lm3','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select

     ! redistribute leftover biomass between biomass pools
     cc%bliving = cc%bliving * (1-fraction_harvested)
     cc%bwood   = cc%bwood   * (1-fraction_harvested)
     call update_biomass_pools(cc)
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
  real :: frac_harvested         ! fraction of biomass harvested this time
  real :: frac_wood_wasted       ! fraction of wood wasted during transition
  real :: frac_wood_wasted_ag    ! fraction of above-ground wood wasted during transition
  real :: wood_harvested         ! amount of harvested wood, kgC/m2
  real :: bdead, balive, bleaf, bfroot, btotal ! combined biomass pools
  real :: delta
  integer :: i

  if (new_landuse==LU_RANGE) return ! do nothing for the conversion to rangeland

  associate (vegn=>tile%vegn)

  balive = 0 ; bdead = 0 ; bleaf = 0 ; bfroot = 0
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
     associate(cc => vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))

     ! calculate total amount of harvested wood, minus the wasted part
     wood_harvested = (cc%bwood+cc%bsw)*frac_harvested*(1-frac_wood_wasted)

     ! distribute harvested wood between pools
     if (new_landuse==LU_SCND) then
        ! this is harvesting, distribute between 3 different wood pools
        vegn%harv_pool_C(HARV_POOL_WOOD_FAST) = vegn%harv_pool_C(HARV_POOL_WOOD_FAST) &
             + wood_harvested*frac_wood_fast
        vegn%harv_pool_C(HARV_POOL_WOOD_MED) = vegn%harv_pool_C(HARV_POOL_WOOD_MED) &
             + wood_harvested*frac_wood_med
        vegn%harv_pool_C(HARV_POOL_WOOD_SLOW) = vegn%harv_pool_C(HARV_POOL_WOOD_SLOW) &
             + wood_harvested*frac_wood_slow
        ! store harvested wood amount, for diagnostics
        vegn%amount_wood_harv_C = wood_harvested
     else
        ! this is land clearance: everything goes into "cleared" pool
        vegn%harv_pool_C(HARV_POOL_CLEARED) = vegn%harv_pool_C(HARV_POOL_CLEARED) &
             + wood_harvested
        ! store cleared wood amount, for diagnostics
        vegn%amount_wood_cleared_C = wood_harvested
     endif

     ! distribute wood and living biomass between fast and slow intermediate
     ! soil carbon pools according to fractions specified through the namelists
     delta = (cc%bwood+cc%bsw)*frac_harvested*frac_wood_wasted
     if(delta<0) call land_error_message('vegn_cut_forest_lm3: '// &
          'harvested amount of dead biomass ('//string(delta)//' kgC/m2) is below zero', &
          FATAL)

     select type (soilc => tile%soilc)
     class is (soil_BGC_SIMPLE_t)
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + delta*(1-sp%fsc_wood)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + delta*   sp%fsc_wood

        delta = balive * frac_harvested
        if(delta<0) call land_error_message('vegn_cut_forest_lm3: '// &
             'harvested amount of live biomass ('//string(delta)//' kgC/m2) is below zero', &
             FATAL)
        vegn%ssc_pool_bg = vegn%ssc_pool_bg + delta*(1-sp%fsc_liv)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + delta*   sp%fsc_liv
     class is (soil_BGC_CORPSE_t)
        delta = (cc%bwood+cc%bsw)*frac_harvested*agf_bs*frac_wood_wasted_ag
        vegn%litter_buff_C(:,LITT_CWOOD) = vegn%litter_buff_C(:,LITT_CWOOD) + &
            [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*delta

        delta = (cc%bl+cc%blv) * frac_harvested
        if(delta<0) call land_error_message('vegn_cut_forest_lm3: '// &
                'harvested amount of live biomass ('//string(delta)//' kgC/m2) is below zero', &
                FATAL)

        vegn%litter_buff_C(:,LITT_LEAF) = vegn%litter_buff_C(:,LITT_LEAF) + &
            [sp%fsc_liv, 1-sp%fsc_liv, 0.0]*delta

        vegn%ssc_pool_bg = vegn%ssc_pool_bg + cc%br*frac_harvested*(1-sp%fsc_froot)
        vegn%fsc_pool_bg = vegn%fsc_pool_bg + cc%br*frac_harvested*sp%fsc_froot

        if(waste_below_ground_wood) then
          vegn%ssc_pool_bg = vegn%ssc_pool_bg + (cc%bwood+cc%bsw)*frac_harvested*(1-agf_bs)*(1-sp%fsc_wood)
          vegn%fsc_pool_bg = vegn%fsc_pool_bg + (cc%bwood+cc%bsw)*frac_harvested*(1-agf_bs)*sp%fsc_wood
        endif

        if (do_CORPSE_nitrogen) then
            vegn%litter_buff_N(:,LITT_CWOOD) = vegn%litter_buff_N(:,LITT_CWOOD) + (&
                  [sp%fsc_wood, 1-sp%fsc_wood, 0.0]*cc%wood_N +&
                  [sp%fsc_liv,  1-sp%fsc_liv,  0.0]*(cc%sapwood_N+cc%stored_N)&
               )*frac_harvested*agf_bs*frac_wood_wasted_ag
            vegn%litter_buff_N(:,LITT_LEAF) = vegn%litter_buff_N(:,LITT_LEAF) + &
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
           vegn%harv_pool_N(HARV_POOL_WOOD_FAST) = vegn%harv_pool_N(HARV_POOL_WOOD_FAST) &
                + (cc%bwood/sp%wood_c2n+cc%bsw/sp%sapwood_c2n)*frac_harvested*(1-frac_wood_wasted)*frac_wood_fast
           vegn%harv_pool_N(HARV_POOL_WOOD_MED) = vegn%harv_pool_N(HARV_POOL_WOOD_MED) &
                + (cc%bwood/sp%wood_c2n+cc%bsw/sp%sapwood_c2n)*frac_harvested*(1-frac_wood_wasted)*frac_wood_med
           vegn%harv_pool_N(HARV_POOL_WOOD_SLOW) = vegn%harv_pool_N(HARV_POOL_WOOD_SLOW) &
                + (cc%bwood/sp%wood_c2n+cc%bsw/sp%sapwood_c2n)*frac_harvested*(1-frac_wood_wasted)*frac_wood_slow
        else
           ! this is land clearance: everything goes into "cleared" pool
           vegn%harv_pool_N(HARV_POOL_CLEARED) = vegn%harv_pool_N(HARV_POOL_CLEARED) &
                + (cc%bwood/sp%wood_c2n+cc%bsw/sp%sapwood_c2n)*frac_harvested*(1-frac_wood_wasted)
        endif

     class default
        call error_mesg('vegn_cut_forest_lm3','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
     end select

     cc%bliving = cc%bliving*(1-frac_harvested)
     cc%bwood   = cc%bwood*(1-frac_harvested)
     cc%leaf_N = cc%leaf_N*(1-frac_harvested)
     cc%root_N = cc%root_N*(1-frac_harvested)
     cc%wood_N = cc%wood_N*(1-frac_harvested)
     cc%sapwood_N = cc%sapwood_N*(1-frac_harvested)
     cc%stored_N = cc%stored_N*(1-frac_harvested)
     ! Should stored N be lost or retained?
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc)
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

  real :: littC, littN ! litter from grazing of individual cohorts, kg/m2
  real :: buffC(N_C_TYPES), buffN(N_C_TYPES) ! accumulators of litter, kg/m2
  real :: LAI ! leaf area index of *grazed* vegetation. Does not include plants taller
              ! than max grazing height.
  integer :: i
  ! variables for conservation checks
  real :: lmass0, fmass0, heat0, cmass0, nmass0

  call check_conservation_1(tile, lmass0,fmass0,cmass0,nmass0,heat0)

  associate (vegn=>tile%vegn,soil=>tile%soil)
  LAI = 0.0
  do i = 1,vegn%n_cohorts
     if (vegn%cohorts(i)%height < max_grazing_height) &
             LAI = LAI + vegn%cohorts(i)%leafarea*vegn%cohorts(i)%nindivs
  enddo
  if (LAI < min_lai_for_grazing) return

  buffC(:) = 0.0; buffN(:) = 0.0
  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))
     if (cc%height >= max_grazing_height) continue ! do nothing to vegetation browsers cannot reach.

     vegn%harv_pool_C(HARV_POOL_PAST) = vegn%harv_pool_C(HARV_POOL_PAST) + &
          cc%bl*grazing_intensity*(1-grazing_residue)*cc%nindivs
     vegn%harv_pool_N(HARV_POOL_PAST) = vegn%harv_pool_N(HARV_POOL_PAST) + &
          cc%leaf_N*grazing_intensity*(1-grazing_residue)*cc%nindivs
     littC = cc%bl     * grazing_intensity*grazing_residue*cc%nindivs
     littN = cc%leaf_N * grazing_intensity*grazing_residue*cc%nindivs
     cc%bl     = cc%bl     * (1-grazing_intensity)
     cc%leaf_N = cc%leaf_N * (1-grazing_intensity)

     ! accumulate litter input
     buffC(:) = buffC(:) + littC*[sp%fsc_liv,1-sp%fsc_liv,0.0]
     buffN(:) = buffN(:) + littN*[sp%fsc_liv,1-sp%fsc_liv,0.0]
     end associate ! cohorts and spdata
  enddo
  ! If grazing is daily, litter goes directly to the soil litter pools; otherwise (in
  ! case of annual grazing), it goes into intermediate buffers to be gradually transferred
  ! into the soil pools later.
  if (grazing_freq==GRAZING_DAILY) then
     call tile%soilc%add_soil_matter(vegn, leaf_litter_C=buffC, leaf_litter_N=buffN )
  else
     ! litter goes to intermediate pool directly
     vegn%litter_buff_C(:,LITT_LEAF) = vegn%litter_buff_C(:,LITT_LEAF) + buffC(:)
     vegn%litter_buff_N(:,LITT_LEAF) = vegn%litter_buff_N(:,LITT_LEAF) + buffN(:)
  endif
  end associate ! vegn

  call check_conservation_2(tile,'vegn_graze_pasture_ppa',lmass0,fmass0,cmass0,nmass0,heat0)
end subroutine vegn_graze_pasture_ppa

! ============================================================================
! NOTE that the PPA harvest would not work properly if applied only once per year, because
! at the end of the year the leaves are down in NH, and therefore harvest amount will be
! unrealistically small. The same is true about the grazing.
subroutine vegn_harvest_crop_ppa(tile, ctag) ! debug_pjp
character(len=*), intent(in) :: ctag ! debug_pjp
  type(land_tile_type), intent(inout) :: tile

  real, dimension(N_C_TYPES) :: &
     leaf_litt_C, leaf_litt_N, & ! accumulated leaf litter, kg/m2
     wood_litt_C, wood_litt_N    ! accumulated wood litter, kg/m2
  real, dimension(num_l, N_C_TYPES) :: &
     root_litt_C, root_litt_N    ! accumulated root litter per soil layer, kg/m2
  real :: ndead  ! number of individuals killed in cohort
  real :: dheat  ! heat residual in cohort merge
  integer :: i
  ! variables for conservation checks
  real :: lmass0, fmass0, heat0, cmass0, nmass0

  call check_conservation_1(tile, lmass0,fmass0,cmass0,nmass0,heat0)

  associate(vegn=>tile%vegn)
  leaf_litt_C=0.0; wood_litt_C=0.0; root_litt_C=0.0
  leaf_litt_N=0.0; wood_litt_N=0.0; root_litt_N=0.0
  do i = 1, vegn%n_cohorts
     associate(cc=>vegn%cohorts(i), sp=>spdata(vegn%cohorts(i)%species))
     ndead = cc%nindivs ! harvest everything

     ! add C to harvest and litter pools. This is slightly different from kill_plants_ppa
     ! because above-ground part of nsc and sapwood are getting harvested, rather than going
     ! into litter
     vegn%harv_pool_C(HARV_POOL_CROP) = vegn%harv_pool_C(HARV_POOL_CROP) + ndead*( &
         cc%bl + cc%bseed + cc%carbon_gain + cc%growth_previous_day + &
         cc%bsw*agf_bs + cc%nsc*agf_bs &
         )
     vegn%harv_pool_N(HARV_POOL_CROP) = vegn%harv_pool_N(HARV_POOL_CROP) + ndead*( &
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
  vegn%litter_buff_C(:,LITT_CWOOD) = vegn%litter_buff_C(:,LITT_CWOOD) + wood_litt_C(:)
  vegn%litter_buff_N(:,LITT_CWOOD) = vegn%litter_buff_N(:,LITT_CWOOD) + wood_litt_N(:)
  vegn%litter_buff_C(:,LITT_LEAF)  = vegn%litter_buff_C(:,LITT_LEAF)  + leaf_litt_C(:)
  vegn%litter_buff_N(:,LITT_LEAF)  = vegn%litter_buff_N(:,LITT_LEAF)  + leaf_litt_N(:)

  vegn%fsc_pool_bg = vegn%fsc_pool_bg + sum(root_litt_C(:,C_FAST))+sum(root_litt_C(:,C_MIC))
  vegn%fsn_pool_bg = vegn%fsn_pool_bg + sum(root_litt_N(:,C_FAST))+sum(root_litt_N(:,C_MIC))
  vegn%ssc_pool_bg = vegn%ssc_pool_bg + sum(root_litt_C(:,C_SLOW))
  vegn%ssn_pool_bg = vegn%ssn_pool_bg + sum(root_litt_N(:,C_SLOW))

  call vegn_relayer_cohorts_ppa(vegn)
  call vegn_mergecohorts_ppa(vegn, dheat)
  tile%e_res_2 = tile%e_res_2 - dheat
  end associate ! vegn

  call check_conservation_2(tile,trim(ctag)//' vegn_harvest_crop_ppa',lmass0,fmass0,cmass0,nmass0,heat0)
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
  real :: dbh_min ! minimum DBH of harvested trees
  real :: dheat  ! heat residual in cohort merge
  integer :: i
  ! variables for conservation checks
  real :: lmass0, fmass0, heat0, cmass0, nmass0

  if (new_landuse==LU_RANGE) return ! do nothing for the conversion to rangeland

  call check_conservation_1(tile, lmass0,fmass0,cmass0,nmass0,heat0)

  associate(vegn=>tile%vegn)
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
        call kill_plants_ppa(cc,vegn,ndead,0.0, leaf_litt_C, wood_litt_C, root_litt_C, &
                                                leaf_litt_N, wood_litt_N, root_litt_N  )
     endif
     end associate
     ! note that below-ground wood all goes to litter, by construction of kill_plants_ppa
  enddo

  ! distribute harvested wood between pools
  if (new_landuse==LU_SCND) then
     ! this is harvesting, distribute between 3 different wood pools
     vegn%harv_pool_C(HARV_POOL_WOOD_FAST) = vegn%harv_pool_C(HARV_POOL_WOOD_FAST) &
          + sum(wood_harv_C)*frac_wood_fast*(1-frac_wood_wasted)
     vegn%harv_pool_N(HARV_POOL_WOOD_FAST) = vegn%harv_pool_N(HARV_POOL_WOOD_FAST) &
          + sum(wood_harv_N)*frac_wood_fast*(1-frac_wood_wasted)
     vegn%harv_pool_C(HARV_POOL_WOOD_MED) = vegn%harv_pool_C(HARV_POOL_WOOD_MED) &
          + sum(wood_harv_C)*frac_wood_med*(1-frac_wood_wasted)
     vegn%harv_pool_N(HARV_POOL_WOOD_MED) = vegn%harv_pool_N(HARV_POOL_WOOD_MED) &
          + sum(wood_harv_N)*frac_wood_med*(1-frac_wood_wasted)
     vegn%harv_pool_C(HARV_POOL_WOOD_SLOW) = vegn%harv_pool_C(HARV_POOL_WOOD_SLOW) &
          + sum(wood_harv_C)*frac_wood_slow*(1-frac_wood_wasted)
     vegn%harv_pool_N(HARV_POOL_WOOD_SLOW) = vegn%harv_pool_N(HARV_POOL_WOOD_SLOW) &
          + sum(wood_harv_N)*frac_wood_slow*(1-frac_wood_wasted)
     ! store harvested wood amount, for diagnostics. We could send the diagnostics
     ! from here, but we are currently not merging the diag buffers when merging
     ! land tiles, and therefore the output would be incorrect if the tiles
     ! that are just harvested are merged after land use transitions and before
     ! dumping the diag (which is likely). Same note applies to amount_wood_cleared_*
     ! below
     vegn%amount_wood_harv_C = sum(wood_harv_C)*(1-frac_wood_wasted)
     vegn%amount_wood_harv_N = sum(wood_harv_N)*(1-frac_wood_wasted)
  else
     ! this is land clearance: everything goes into "cleared" pool
     vegn%harv_pool_C(HARV_POOL_CLEARED) = vegn%harv_pool_C(HARV_POOL_CLEARED) &
          + sum(wood_harv_C)*(1-frac_wood_wasted)
     vegn%harv_pool_N(HARV_POOL_CLEARED) = vegn%harv_pool_N(HARV_POOL_CLEARED) &
          + sum(wood_harv_N)*(1-frac_wood_wasted)
     ! cleared wood amount, for diagnostics
     vegn%amount_wood_cleared_C = sum(wood_harv_C)*(1-frac_wood_wasted)
     vegn%amount_wood_cleared_N = sum(wood_harv_N)*(1-frac_wood_wasted)
  endif

  vegn%litter_buff_C(:,LITT_CWOOD) = vegn%litter_buff_C(:,LITT_CWOOD) + &
     wood_litt_C(:) + wood_harv_C(:)*frac_wood_wasted
  vegn%litter_buff_N(:,LITT_CWOOD) = vegn%litter_buff_N(:,LITT_CWOOD) + &
     wood_litt_N(:) + wood_harv_N(:)*frac_wood_wasted
  vegn%litter_buff_C(:,LITT_LEAF) = vegn%litter_buff_C(:,LITT_LEAF) + leaf_litt_C(:)
  vegn%litter_buff_N(:,LITT_LEAF) = vegn%litter_buff_N(:,LITT_LEAF) + leaf_litt_N(:)

  vegn%fsc_pool_bg = vegn%fsc_pool_bg + sum(root_litt_C(:,C_FAST))+sum(root_litt_C(:,C_MIC))
  vegn%fsn_pool_bg = vegn%fsn_pool_bg + sum(root_litt_N(:,C_FAST))+sum(root_litt_N(:,C_MIC))
  vegn%ssc_pool_bg = vegn%ssc_pool_bg + sum(root_litt_C(:,C_SLOW))
  vegn%ssn_pool_bg = vegn%ssn_pool_bg + sum(root_litt_N(:,C_SLOW))

  call vegn_relayer_cohorts_ppa(vegn)
  call vegn_mergecohorts_ppa(vegn, dheat)
  tile%e_res_2 = tile%e_res_2 - dheat
  end associate ! vegn

  call check_conservation_2(tile,'vegn_cut_forest_ppa',lmass0,fmass0,cmass0,nmass0,heat0)
end subroutine vegn_cut_forest_ppa

! ============================================================================
! this function uses the same rule as LM3 does for biogeographic C3/C4,
! distribution except it disregards the biomass, that is returns the
! physiology type for grases that would be optimal for given annual T and P.
function biogeographic_physiology_type(temp, precip) result (pt)
  integer :: pt
  real, intent(in) :: temp   ! temperature, degK
  real, intent(in) :: precip ! precipitation, ???

  real :: pc4

  ! Rule based on analysis of ED global output; equations from JPC, 2/02
  pc4=exp(-0.0421*(273.16+25.56-temp)-(0.000048*(273.16+25.5-temp)*precip))

  if(pc4>0.5) then
    pt=PT_C4
  else
    pt=PT_C3
  endif
end function biogeographic_physiology_type

! ============================================================================
subroutine vegn_plant_crop_ppa(tile, chosen_crop)
  type(land_tile_type), intent(inout) :: tile
  integer, optional, intent(in) :: chosen_crop

  ! list of pools we borrow seeds from, highest priority first
  integer, parameter :: seed_source_pools(6) = &
     [ HARV_POOL_CROP, HARV_POOL_PAST, HARV_POOL_CLEARED, HARV_POOL_WOOD_FAST, HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW ]

  integer :: i, p, pt, crop_species_idx
  real, dimension(0:nspecies-1) :: seedC, seedN ! seed biomass, kg/m2
  real :: deltaC, deltaN ! amount we borrow from harvest pools, kg/m2
  ! variables for conservation checks
  real :: lmass0, fmass0, heat0, cmass0, nmass0

  call check_conservation_1(tile, lmass0,fmass0,cmass0,nmass0,heat0)

  ! prepare cropland for planting: right now just kill all vegetation; in the
  ! future we possibly need to add some soil carbon mixing by plows, perhaps
  ! other agricultural processes
  if (clear_crop_before_planting) call vegn_cut_forest_ppa(tile, tile%vegn%landuse)

  ! determine crop species: now using the same biogeography rules that LM3 was using
  ! to determine c3/c4 photosynthesis type
  associate (vegn=>tile%vegn, soil=>tile%soil)
  select case(crop_distribution_option)
  case (CROP_DISTR_LM3)
     pt = biogeographic_physiology_type(tile%vegn%t_ann, tile%vegn%p_ann*seconds_per_year)
     select case(pt)
     case (PT_C3)
        crop_species_idx = c3_crop_idx
        if (crop_species_idx<0) call land_error_message('vegn_plant_crop_ppa: C3 crop species "'//trim(c3_crop_species)//'" not found.', FATAL)
     case (PT_C4)
        crop_species_idx = c4_crop_idx
        if (crop_species_idx<0) call land_error_message('vegn_plant_crop_ppa: C4 crop species "'//trim(c4_crop_species)//'" not found.', FATAL)
     case default
        call land_error_message('vegn_plant_crop_ppa: unknown physiology type '//string(pt)//'; this should never happen.', FATAL)
     end select
  case (CROP_DISTR_MIRCA2000)
     if(.not.present(chosen_crop)) then
        call land_error_message('vegn_plant_crop_ppa: calling argument chosen_crop must be present when crop_schedule = computed', FATAL)
     endif
     select case(chosen_crop)
     case (NO_CROP)
        pt = biogeographic_physiology_type(tile%vegn%t_ann, tile%vegn%p_ann*seconds_per_year)
        select case(pt)
        case (PT_C3)
           crop_species_idx = c3_crop_idx
        case (PT_C4)
           crop_species_idx = c4_crop_idx
       case default
          call land_error_message('vegn_plant_crop_ppa: unknown physiology type '//string(pt)//'; this should never happen.', FATAL)
        end select
     case (IRRIGATED_MAIZE)
       crop_species_idx = maize_crop_idx
     case (IRRIGATED_SOYBEAN)
       crop_species_idx = soybean_crop_idx
     case (IRRIGATED_RICE)
       crop_species_idx = rice_crop_idx
     case (IRRIGATED_SPRING_WHEAT)
       crop_species_idx = spring_wheat_crop_idx
     case (IRRIGATED_WINTER_WHEAT)
       crop_species_idx = winter_wheat_crop_idx
     case (RAINFED_MAIZE)
       crop_species_idx = maize_crop_idx
     case (RAINFED_SOYBEAN)
       crop_species_idx = soybean_crop_idx
     case (RAINFED_RICE)
       crop_species_idx = rice_crop_idx
     case (RAINFED_SPRING_WHEAT)
       crop_species_idx = spring_wheat_crop_idx
     case (RAINFED_WINTER_WHEAT)
       crop_species_idx = winter_wheat_crop_idx
     case default
        call land_error_message('vegn_plant_crop_ppa: invalid crop type number='//string(chosen_crop)//'; this should never happen.', FATAL)
     end select
  case default
     call error_mesg('vegn_plant_crop_ppa','Unknown crop distribution option; this should never happen.', FATAL)
  end select

  ! borrow biomass (crop_seed_density) from harvest pools, in order of preference
  seedC(:) = 0.0; seedN(:) = 0.0
  do i = 1, size(seed_source_pools)
     p = seed_source_pools(i)
     deltaC = max(crop_seed_density-seedC(crop_species_idx),0.0)
     deltaC = min(deltaC,vegn%harv_pool_C(p))
     if (vegn%harv_pool_C(p)>0) then
        deltaN = deltaC*vegn%harv_pool_N(p)/vegn%harv_pool_C(p)
     else
        deltaN = 0.0
     endif
     vegn%harv_pool_C(p) = vegn%harv_pool_C(p) - deltaC
     vegn%harv_pool_N(p) = vegn%harv_pool_N(p) - deltaN
     seedC(crop_species_idx) = seedC(crop_species_idx) + deltaC
     seedN(crop_species_idx) = seedN(crop_species_idx) + deltaN
     if (seedC(crop_species_idx) >= crop_seed_density) exit ! from loop
  enddo
  call check_var_range(seedC(crop_species_idx),0.99*crop_seed_density,HUGE(1.0),'vegn_plant_crop_ppa','seedC',WARNING)

  call add_seedlings_ppa(vegn,soil,tile%soilc,seedC,seedN, prob_est = 1.0, prob_ger = 1.0)
  end associate ! vegn,soil

  call check_conservation_2(tile,'vegn_plant_crop_ppa', lmass0,fmass0,cmass0,nmass0,heat0)
end subroutine vegn_plant_crop_ppa

! ============================================================================
! transport crops horizontally to satisfy demand on the planting day
subroutine crop_seed_transport(day_of_year)
  integer :: day_of_year

  ! local vars
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  integer ::  l ! current point index
  real :: total_seed_supply_C, total_seed_supply_N
  real :: total_seed_demand_C, total_seed_demand_N
  real :: crop_seed_supply_C, crop_seed_supply_N, crop_seed_demand_C, crop_seed_demand_N
  real :: f_supply_C, f_supply_N ! fraction of the supply that gets spent
  real :: f_demand_C, f_demand_N ! fraction of the demand that gets satisfied

  ! supply/demand on unstructured grid
  real, dimension (lnd%ls:lnd%le) :: &
        crop_seed_supply_C_UG, crop_seed_supply_N_UG, &
        crop_seed_demand_C_UG, crop_seed_demand_N_UG

  real :: btot0, btot1 ! total carbon, for conservation check only
  real :: ntot0, ntot1 ! total nitrogen, for conservation check only

  if(.not.transport_crop_seeds) return

  ! + conservation check part 1
  if (do_check_conservation) then
     btot0 = 0.0; ntot0 = 0.0
     ce = first_elmt(land_tile_map,lnd%ls)
     do while (loop_over_tiles(ce,tile,l))
        btot0 = btot0 + lnd%ug_area(l) * tile%frac * land_tile_carbon(tile)
        ntot0 = ntot0 + lnd%ug_area(l) * tile%frac * land_tile_nitrogen(tile)
     end do
     ! this will likely not reproduce across PE count, but that is OK since it is only
     ! used for conservaton checks
     call mpp_sum(btot0); call mpp_sum(ntot0)
  end if
  ! - conservation check part 1

  crop_seed_supply_C_UG(:) = 0.0; crop_seed_demand_C_UG(:) = 0.0
  crop_seed_supply_N_UG(:) = 0.0; crop_seed_demand_N_UG(:) = 0.0
  do l = lnd%ls, lnd%le
     ce = first_elmt(land_tile_map(l))
     do while (loop_over_tiles(ce,tile))
        if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

        call crop_seed_supply(tile%vegn,crop_seed_supply_C,crop_seed_supply_N)
        crop_seed_supply_C_UG(l) = crop_seed_supply_C_UG(l) + crop_seed_supply_C*tile%frac*lnd%ug_area(l)
        crop_seed_supply_N_UG(l) = crop_seed_supply_N_UG(l) + crop_seed_supply_N*tile%frac*lnd%ug_area(l)
        call crop_seed_demand(tile%vegn,l,day_of_year,crop_seed_demand_C,crop_seed_demand_N)
        crop_seed_demand_C_UG(l) = crop_seed_demand_C_UG(l) + crop_seed_demand_C*tile%frac*lnd%ug_area(l)
        crop_seed_demand_N_UG(l) = crop_seed_demand_N_UG(l) + crop_seed_demand_N*tile%frac*lnd%ug_area(l)
     enddo
  enddo
  ! sum totals globally
  total_seed_demand_C = land_global_sum_UG(crop_seed_demand_C_UG)
  total_seed_demand_N = land_global_sum_UG(crop_seed_demand_N_UG)
  total_seed_supply_C = land_global_sum_UG(crop_seed_supply_C_UG)
  total_seed_supply_N = land_global_sum_UG(crop_seed_supply_N_UG)
  ! if either demand or supply are zeros we don't need (or can't) transport anything
  if (total_seed_demand_C==0)then
     return
  end if
  if (total_seed_supply_C==0)then
     call error_mesg('crop_seed_transport '//string_from_time(lnd%time), &
        'total seed C supply is zero, but demand is not:'//string(total_seed_demand_C), NOTE)
     return
  endif

  ! calculate the fraction of the supply that is going to be used
  f_supply_C = MIN(total_seed_demand_C/total_seed_supply_C, 1.0)
  if (total_seed_supply_N > 0) then
     f_supply_N = MIN(total_seed_demand_N/total_seed_supply_N, 1.0)
  else
     total_seed_supply_N = 0.0
     f_supply_N          = 1.0
  endif
  ! calculate the fraction of the demand that is going to be satisfied
  f_demand_C = MIN(total_seed_supply_C/total_seed_demand_C, 1.0)
  if (total_seed_demand_N > 0) then
     f_demand_N = MIN(total_seed_supply_N/total_seed_demand_N, 1.0)
  else
     total_seed_demand_N = 0.0
     f_demand_N          = 1.0
  endif
  ! note that either f_supply or f_demand is 1; the mass conservation law in the
  ! following calculations is satisfied since
  ! f_demand*total_seed_demand == f_supply*total_seed_supply
  call error_mesg('crop_seed_transport '//string_from_time(lnd%time), &
     'fraction of C demand satisfied='//string(f_demand_C)//' fraction of C supply used ='//string(f_supply_C), NOTE)
  call error_mesg('crop_seed_transport '//string_from_time(lnd%time), &
     'fraction of N demand satisfied='//string(f_demand_N)//' fraction of N supply used ='//string(f_supply_N), NOTE)

  ! redistribute part (or possibly all) of the supply to satisfy part (or possibly all)
  ! of the demand. This relies on the assumption that supply and demand did not change
  ! since the last calculation above
  ce = first_elmt(land_tile_map, lnd%ls)
  do while (loop_over_tiles(ce,tile,l))
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body
     call crop_seed_supply(tile%vegn,crop_seed_supply_C,crop_seed_supply_N)
     call crop_seed_demand(tile%vegn,l,day_of_year,crop_seed_demand_C,crop_seed_demand_N)
     tile%vegn%harv_pool_C(HARV_POOL_CROP) = tile%vegn%harv_pool_C(HARV_POOL_CROP) + &
                      f_demand_C*crop_seed_demand_C - f_supply_C*crop_seed_supply_C
     tile%vegn%harv_pool_N(HARV_POOL_CROP) = tile%vegn%harv_pool_N(HARV_POOL_CROP) + &
                      f_demand_N*crop_seed_demand_N - f_supply_N*crop_seed_supply_N
  enddo

  ! + conservation check part 2
  if (do_check_conservation) then
     btot1 = 0.0; ntot1 = 0.0
     ce = first_elmt(land_tile_map,lnd%ls)
     do while (loop_over_tiles(ce,tile,l))
        btot1 = btot1 + lnd%ug_area(l) * tile%frac * land_tile_carbon(tile)
        ntot1 = ntot1 + lnd%ug_area(l) * tile%frac * land_tile_nitrogen(tile)
     end do
     call mpp_sum(btot1) ; call mpp_sum(ntot1)
     if (mpp_pe()==mpp_root_pe()) then
        call check_conservation ('vegn_reproduction_ppa','total carbon', &
             btot0/tot_area_land, btot1/tot_area_land, carbon_cons_tol, severity=FATAL)
        call check_conservation ('vegn_reproduction_ppa','total nitrogen', &
             ntot0/tot_area_land, ntot1/tot_area_land, nitrogen_cons_tol, severity=FATAL)
     endif
  end if
  ! - conservation check part 2

end subroutine crop_seed_transport

! given an array on structural grid, calculates sum of all elements in a way that
! is supposed to reproduce across different PE counts and layouts
function land_global_sum_UG(a) result(s)
  real, intent(in) :: a(:) ! data to sum up
  real :: s ! resulting sum

  real :: a2D(lnd%is:lnd%ie,lnd%js:lnd%je) ! input field on structured grid

  a2D = 0.0
  call mpp_pass_UG_to_SG(lnd%ug_domain, a, a2D)
  s = mpp_global_sum(lnd%sg_domain, a2D, flags=BITWISE_EXACT_SUM)
end function

! ============================================================================
subroutine crop_seed_supply(vegn, crop_seed_supply_C, crop_seed_supply_N)
   type(vegn_tile_type), intent(in) :: vegn
   real, intent(out) :: crop_seed_supply_C, crop_seed_supply_N
   crop_seed_supply_C = MAX(vegn%harv_pool_C(HARV_POOL_CROP)-crop_seed_density,0.0)
   if (vegn%harv_pool_C(HARV_POOL_CROP)>0) then
      crop_seed_supply_N = crop_seed_supply_C*vegn%harv_pool_N(HARV_POOL_CROP)/vegn%harv_pool_C(HARV_POOL_CROP)
   else
      crop_seed_supply_N = 0.0
   endif
end subroutine crop_seed_supply

! ============================================================================
subroutine crop_seed_demand(vegn, l, day_of_year, crop_seed_demand_C, crop_seed_demand_N)
   type(vegn_tile_type), intent(in) :: vegn
   integer, intent(in) :: l ! index of grid cell
   integer, intent(in) :: day_of_year
   real, intent(out) :: crop_seed_demand_C, crop_seed_demand_N

   crop_seed_demand_C = 0.0; crop_seed_demand_N = 0.0
   if (is_cropland(vegn%landuse)) then
     if(crop_schedule_option == CROP_SCHEDULE_COMPUTED) then
        if(day_of_year==vegn%Crop%chosen_calendars(1,1) .or. day_of_year==vegn%Crop%chosen_calendars(1,2)) then
          crop_seed_demand_C = MAX(crop_seed_density               - vegn%harv_pool_C(HARV_POOL_CROP),0.0)
          crop_seed_demand_N = MAX(crop_seed_density/crop_seed_c2n - vegn%harv_pool_N(HARV_POOL_CROP),0.0)
        endif
     else
        if(day_of_year==nint(crop_planting_day(l))) then
          crop_seed_demand_C = MAX(crop_seed_density               - vegn%harv_pool_C(HARV_POOL_CROP),0.0)
          crop_seed_demand_N = MAX(crop_seed_density/crop_seed_c2n - vegn%harv_pool_N(HARV_POOL_CROP),0.0)
        endif
     endif
   endif
end subroutine

! ============================================================================
subroutine save_harvesting_restart(tile_dim_length,timestamp)
   integer, intent(in) :: tile_dim_length
   character(*), intent(in) :: timestamp

   if(crop_schedule_option == CROP_SCHEDULE_COMPUTED) call save_crop_restart(tile_dim_length,timestamp)
end subroutine save_harvesting_restart

end module
