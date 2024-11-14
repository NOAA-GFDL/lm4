 module vegn_crop_mod

! This module executes an algorithm that determines dates for planting and harvesting of four
! major agricultural crops as a function of local climate. These crops are: maize, soybean,
! wheat and rice. The algorithm determines the most favorable dates for planting and a range of
! dates over which climatic conditions are favorable, and the corresponding harvest dates.
! Up to four pairs of planting and harvest dates are possible for each crop type. Where climatic
! conditions allow it, two crops may be grown on the same land within the same year. This is most
! common in low latitudes. This yields "main" and "second" cropping seasons for each crop type.
! For each of these, planting and harvest dates may change depending on whether or not irrigation
! is available. Indeed, the possibility of cultivating a particular crop at all may depend on the
! availability of irrigation. This yields "rainfed" and "irrigated" dates. The rainfed and irrigated
! calendars often differ where climates are arid or semi-arid, but are most often identical where
! precipitation is adaquate. In all cases, planting and harvests dates are set the the value of
! "NO_DATE" when climatic conditions do not allow cultivation. For example, an arid location where
! irrigation is available may have the rainfed dates set to NO_DATE and have valid calendar dates
! for the irrigated calendar. These dates, taken together, are referred to as a crop calendar.
! vegn_tile.F90 has more inline documentation regarding the crop_calendars and associated data.

! There are two routines that execute the crop calendar algorithm, CCA_Wheat and CCA_Maize_Soybean_Rice.
! Wheat is subdivided into two broad catagories: spring wheat and winter wheat. Winter wheat requires
! a period of low temperatures early in its growth cycle whereas spring wheat does not. Therefore,
! winter wheat is planted in autumn and spring wheat in the spring. This fact greatly impacts the
! dates of planting and harvest, requiring them to be treated as separate crops. subroutine CCA_Wheat
! has a single calling argument telling it if it is dealing with spring or winter wheat.
! Therefore, CCA_Wheat is called twice, once for each catagory of wheat. Aside from the low
! temperature requirement, the algorithm for spring and winter wheat is similar in many ways,
! allowing both wheat types to be handled within a single routine.

! Wheat is rarely, if ever, double-cropped with itself. Therefore, subroutine CCA_Wheat returns
! crop calendars for only a main season of rainfed and irrigated wheat.
! The code fills in the second season dates with NO_DATE after the calls to CCA_Wheat.

! Subroutine CCA_Maize_Soybean_Rice executes the algorithm for Maize, Soybean and Rice. The algorithm
! is identical for these three crops, except for the numerical values of some key parameters, allowing
! all three crops to be handled by a single routine. Unlike CCA_Wheat, CCA_Maize_Soybean_Rice returns
! dates for main and second season crops instead of irrigated and rainfed. Therefore,
! CCA_Maize_Soybean_Rice is called twice, once each for the irrigated and rainfed calendars.
! The water source (irrigated or rainfed) is specified by an input argument to the routine.

! The returned values from these routines comprise a set of crop calendars.
! All dates are given as the day number of the year, 1 to 365.
! These dates are stored in array crop_calendars.
! See vegn_tile.F90 for documentation regarding array crop_calendars.

! subroutine crop_selection selects the crop or crops that will be cultivated on each vegetation
! tile which contains the crop landuse type. This routine chooses two crop calendars from the
! crop_calendars array to serve as the main and second season crops. These dates are stored in array
! chosen_calendars. The crop type is stored in a separate array, chosen_crop.
! See vegn_tile.F90 for documentation regarding arrays chosen_calendars and chosen_crop.
! These two arrays are then used by vegn_harvesting_mod to decide what and when to plant and harvest.

#include "../shared/debug.inc"
 use mpp_mod, only: input_nml_file, get_unit
 use fms_mod, only: error_mesg, NOTE, WARNING, FATAL, file_exist, check_nml_error, stdlog, stdout
 use time_manager_mod, only: time_type, set_date, get_date, operator(-), set_time, operator(+), length_of_year, operator(//), operator(<)
 use constants_mod, only: TFREEZE, SECONDS_PER_DAY, PI
 use land_tile_mod, only: land_tile_type, land_tile_enum_type, first_elmt, loop_over_tiles, land_tile_map
 use vegn_tile_mod, only: vegn_tile_type
 use vegn_data_mod, only: NO_DATE, NO_CROP, MAIZE, SOYBEAN, RICE, SPRING_WHEAT, WINTER_WHEAT, IRRIGATED, RAINFED, MAIN_SEASON, SECOND_SEASON, season_name, &
                          IDLE, ACTIVE_ON_LM3_SCHEDULE, crop_name, num_crop_types, num_crop_cal, num_crop_seasons, num_crop_water_sources, &
                          num_crop_periods, water_source_name, LU_CROP
 use land_data_mod, only: lnd
 use land_tile_io_mod, only: land_restart_type, init_land_restart, open_land_restart, save_land_restart, &
                             free_land_restart, add_restart_axis, add_tile_data, get_tile_data, field_exists, add_int_tile_data, get_int_tile_data
 use astronomy_mod, only: get_orbital_parameters, get_ref_date_of_ae
 use diag_manager_mod,  only: diag_axis_init, send_data, register_static_field, diag_field_add_attribute
 use land_tile_diag_mod,only: register_tiled_diag_field, send_tile_data, diag_buff_type, set_default_diag_filter
 use land_numerics_mod, only: ludcmp, lubksb
 use land_io_mod, only: init_cover_field, read_field
 use fms2_io_mod, only: close_file, FmsNetcdfFile_t, open_file

 implicit none
 private

 public :: vegn_crop_init, vegn_crop_end, save_crop_restart, compute_crop_calendars
 character(len=4), private, parameter :: module_name = 'crop'
 integer, parameter :: num_angles=3600, num_water=2
 integer, parameter :: days_in_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
 character(len=3), parameter :: month_name(12) = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)
 integer, parameter :: num_m = 4
 real, parameter :: aPTTtH_range_SW(2) = (/800.,967./)
 real, parameter :: aPTTtH_range_WW(2) = (/800.,837./)
 real, parameter :: unsuitable = 1000. ! suitablity index is set to unsuitable whenever it exeeds SI_crit. (Is this necessary?)
 real, parameter :: aPTT_interval = 200. ! Wheat suitablity is tested at intervals of aPTT_interval units of photo-thermal time.
 character(len=9), parameter :: cwater(num_water) = (/'irrigated','rainfed  '/)
 character(len=5), parameter :: cseason(num_crop_seasons) = (/' main'," 2'nd"/)
 character(len=16), parameter :: restart_file_name = 'crop.nc'

!-----------------------------------------------------------
! used by routines that compute_day_length
 real :: ecc, obliq, per
 real :: orb_angle(0:num_angles)
 type(time_type) :: period_time_type, autumnal_eq_ref
!-----------------------------------------------------------

 real, parameter, dimension(0:num_m) :: central_T_Maize_orig_units  = (/18.96, 21.80, 23.41, 23.39, 21.12/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_Maize            = (/30.87, 13.87,  8.87,  7.92, 14.10/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_Maize_orig_units  = (/ 3.41,  4.26,  4.46,  4.34,  3.82/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_Maize_orig_units = (/ 3.37,  4.48,  4.49,  5.60,  5.62/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_Maize             = (/.5618, .5842, .5819, .5563, .5188/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_Maize            = (/.001449, .002331, .002202, .001230, .000514/) ! fraction^2

 real, parameter, dimension(0:num_m) :: central_T_Soy_orig_units  = (/21.19, 23.58, 24.35, 23.06, 19.63/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_Soy            = (/21.39,  7.78,  4.41,  7.08, 18.95/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_Soy_orig_units  = (/ 4.45,  4.87,  4.65,  4.02,  3.24/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_Soy_orig_units = (/ 3.75,  5.76,  4.56,  3.06,  3.32/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_Soy             = (/.5837, .5917, .5724, .5356, .4926/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_Soy            = (/.001053, .001243, .000980, .000597, .000441/) ! fraction^2

 real, parameter, dimension(0:num_m) :: central_T_SW_orig_units  = (/12.86, 15.96, 18.35, 20.31, 21.11/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_SW            = (/24.11,  8.33,  5.54,  4.18, 10.62/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_SW_orig_units  = (/ 1.59,  1.94,  1.89,  1.70,  1.49/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_SW_orig_units = (/ 0.58,  0.75,  0.89,  1.02,  1.14/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_SW             = (/.4755, .5291, .5519, .5514, .5331/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_SW            = (/.014516, .017204, .010947, .005253, .002207/) ! fraction^2

 real, parameter, dimension(0:num_m) :: central_T_WW_orig_units  = (/14.47, 12.38, 17.31, 20.49, 22.43/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_WW            = (/17.93,  5.30,  2.57,  2.91,  4.52/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_WW_orig_units  = (/ 1.80,  1.90,  2.04,  2.12,  2.16/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_WW_orig_units = (/ 0.71,  0.48,  0.74,  1.00,  1.43/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_WW             = (/.3992, .4994, .5711, .5860, .5834/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_WW            = (/.002711, .009490, .004064, .002138, .001261/) ! fraction^2

 real, parameter, dimension(0:num_m) :: central_T_Rice_orig_units  = (/24.10, 25.06, 25.93, 26.16, 24.83/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_Rice            = (/33.18, 16.08,  8.35,  6.92, 12.34/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_Rice_orig_units  = (/ 7.19,  8.61,  8.85,  7.44,  5.23/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_Rice_orig_units = (/11.36, 12.17, 11.58, 13.69, 13.60/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_Rice             = (/.5373, .5415, .5347, .5185, .4977/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_Rice            = (/.001309, .001683, .001709, .001265, .000884/) ! fraction^2

 integer, parameter :: GP_Maize = 149, GP_Soy = 142,  GP_Rice = 137
 real,    parameter :: SI_crit_Maize = 26.0, SI_crit_Soy = 28.0, SI_crit_Rice = 20.0, SI_crit_SW = 22.0, SI_crit_WW = 30.0

 real, dimension(0:num_m) :: central_T_Maize, central_P_Maize, variance_P_Maize ! model units
 real, dimension(0:num_m) :: central_T_Soy,   central_P_Soy,   variance_P_Soy   ! model units
 real, dimension(0:num_m) :: central_T_SW,    central_P_SW,    variance_P_SW    ! model units
 real, dimension(0:num_m) :: central_T_WW,    central_P_WW,    variance_P_WW    ! model units
 real, dimension(0:num_m) :: central_T_Rice,  central_P_Rice,  variance_P_Rice  ! model units

! For all crop types, potential planting dates are tested at 5 day
! intervals starting with Jan 5 and ending Dec 31.
! This results in 73 potential planting dates (365/5 = 73)
! Harvest dates, however, can be any day of the year.
 integer, parameter :: num_test_days = 365/5, num_seasons_Rice = 2

 real, allocatable :: day_length(:,:)
 logical :: crop_mod_initialized = .FALSE.
 real, dimension(12,12) :: X_ludcmp
 integer :: indx_ludcmp(12)
 character(len=64) :: restart_fieldname(num_crop_seasons,num_crop_types,num_crop_water_sources)

 real :: t_mid_month(0:13)
! Note that t_mid_month is dimensioned (0:13) where t_mid_month(0) is negative because it is the middle of Dec
! of the previous year and t_mid_month(13) is > 365. because it is the middle of Jan of the following year.

 integer :: id_crop_calendars(2,num_crop_seasons,num_crop_water_sources), id_chosen_calendars(num_crop_seasons)
 integer :: id_T_ave, id_P_ave, id_potential_crop, id_status, id_chosen_crop

 real :: weight_climate=.10
 real :: max_planting_SI_SW = 9.75
 real :: Tbase_Wheat = 5.0 + TFREEZE
 integer :: length_of_vernalization_period = 40
 real :: max_T_for_vernalization = 8.0 + TFREEZE
 real :: min_planting_T_Wheat = 5.0 + TFREEZE
 real :: absolute_min_T_for_Wheat = -8.0 + TFREEZE
 character(len=9) :: water_source = 'rainfed  ' ! valid options are 'rainfed' and 'irrigated'
 real :: SI_crit(num_crop_types)
 integer :: GP(num_crop_types-2)
 real, dimension(0:num_m,num_crop_types) :: central_T, central_P, central_D, variance_T, variance_P, variance_D

 namelist / vegn_crop_nml / weight_climate, &
            max_planting_SI_SW, Tbase_Wheat, length_of_vernalization_period, &
            max_T_for_vernalization, min_planting_T_Wheat, absolute_min_T_for_Wheat, &
            water_source
 contains
!============================================================================
 subroutine compute_crop_calendars(vegn,diag,L)
 type(vegn_tile_type), intent(inout) :: vegn
 type(diag_buff_type), intent(inout) :: diag
 integer, intent(in) :: L ! index of grid cell which contains this tile
 integer :: n1, n2, second, minute, hour, day0, day1, month0, month1, year0, year1, iph, iseason, mth, iwater, ipref, pot_crop
 logical :: new_month
 real, dimension(12) :: rhs
 character(len=128) :: outname

 integer, dimension(2) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end
 real :: aPTTtH_range(2)
 character(len=2) :: Wheat_type

 if(.not.crop_mod_initialized) call error_mesg('compute_crop_calendars','vegn_crop_init has not been called', FATAL)
 call get_date(lnd%time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
 call get_date(lnd%time, year0,month0,day0,hour,minute,second)
 new_month = month0 /= month1
 if(new_month) then
    ! compute new climatological temperature and precipitation rates and new crop calendars once per month
    vegn%Crop%tc_av_climate(month1) = weight_climate*vegn%tc_av + (1-weight_climate)*vegn%Crop%tc_av_climate(month1)
    vegn%Crop%precip_av_climate(month1) = weight_climate*vegn%precip_av + (1-weight_climate)*vegn%Crop%precip_av_climate(month1)
    vegn%Crop%T_mid_mth = 4*vegn%Crop%tc_av_climate ! 4*tc_av_climate is the rhs. lubksb overwrites it with the solution.
    call lubksb(X_ludcmp, indx_ludcmp, vegn%Crop%T_mid_mth)
    vegn%Crop%P_mid_mth = 4*vegn%Crop%precip_av_climate ! 4*tcprecip_av_climate is the rhs. lubksb overwrites it with the solution.
    call lubksb(X_ludcmp, indx_ludcmp, vegn%Crop%P_mid_mth)
    crop_loop: do ipref=1,num_crop_types
      pot_crop = vegn%Crop%potential_crop(ipref)
      if(pot_crop == NO_CROP) exit crop_loop
      if(pot_crop == SPRING_WHEAT) then
        aPTTtH_range = aPTTtH_range_SW
        Wheat_type = 'SW'
      else if(pot_crop == WINTER_WHEAT) then
        aPTTtH_range = aPTTtH_range_WW
        Wheat_type = 'WW'
      endif
     if(pot_crop == SPRING_WHEAT .or. pot_crop == WINTER_WHEAT) then
        call CCA_Wheat(L, vegn, Wheat_type, central_T(:,pot_crop), variance_T(:,pot_crop), central_P(:,pot_crop), & ! intent(in)
                       variance_P(:,pot_crop), central_D(:,pot_crop), variance_D(:,pot_crop), & ! intent(in)
                       SI_crit(pot_crop), max_planting_SI_SW, Tbase_Wheat, aPTTtH_range, &      ! intent(in)
                       length_of_vernalization_period, max_T_for_vernalization, min_planting_T_Wheat, & ! intent(in)
                       pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)
        vegn%Crop%crop_calendars(:,1,MAIN_SEASON,ipref,IRRIGATED) = (/pday(1), hday(1)/)
        vegn%Crop%crop_calendars(:,1,MAIN_SEASON,ipref,RAINFED  ) = (/pday(2), hday(2)/)
        vegn%Crop%crop_calendars(:,2,MAIN_SEASON,ipref,IRRIGATED) = (/pday_beg(1), hday_beg(1)/)
        vegn%Crop%crop_calendars(:,2,MAIN_SEASON,ipref,RAINFED  ) = (/pday_beg(2), hday_beg(2)/)
        vegn%Crop%crop_calendars(:,3,MAIN_SEASON,ipref,IRRIGATED) = (/pday_end(1), hday_end(1)/)
        vegn%Crop%crop_calendars(:,3,MAIN_SEASON,ipref,RAINFED  ) = (/pday_end(2), hday_end(2)/)
        ! Two seasons of Wheat do not exist anywhere, and may not even be possible.
        vegn%Crop%crop_calendars(:,:,SECOND_SEASON,ipref,IRRIGATED) = NO_DATE
        vegn%Crop%crop_calendars(:,:,SECOND_SEASON,ipref,RAINFED  ) = NO_DATE
      else
        water_loop: do iwater=1,num_water
          call CCA_Maize_Soybean_Rice(L, vegn, trim(cwater(iwater)), GP(pot_crop), central_T(:,pot_crop), variance_T(:,pot_crop), & ! intent(in)
                        central_P(:,pot_crop), variance_P(:,pot_crop), central_D(:,pot_crop), variance_D(:,pot_crop), SI_crit(pot_crop), & ! intent(in)
                        pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)
          if(iwater==1) then
            vegn%Crop%crop_calendars(:,1,MAIN_SEASON,  ipref,IRRIGATED) = (/ pday(1), hday(1)/)
            vegn%Crop%crop_calendars(:,1,SECOND_SEASON,ipref,IRRIGATED) = (/ pday(2), hday(2)/)
            vegn%Crop%crop_calendars(:,2,MAIN_SEASON,  ipref,IRRIGATED) = (/ pday_beg(1), hday_beg(1)/)
            vegn%Crop%crop_calendars(:,2,SECOND_SEASON,ipref,IRRIGATED) = (/ pday_beg(2), hday_beg(2)/)
            vegn%Crop%crop_calendars(:,3,MAIN_SEASON,  ipref,IRRIGATED) = (/ pday_end(1), hday_end(1)/)
            vegn%Crop%crop_calendars(:,3,SECOND_SEASON,ipref,IRRIGATED) = (/ pday_end(2), hday_end(2)/)
          endif
          if(iwater==2) then
            vegn%Crop%crop_calendars(:,1,MAIN_SEASON,  ipref,RAINFED) = (/ pday(1), hday(1)/)
            vegn%Crop%crop_calendars(:,1,SECOND_SEASON,ipref,RAINFED) = (/ pday(2), hday(2)/)
            vegn%Crop%crop_calendars(:,2,MAIN_SEASON,  ipref,RAINFED) = (/ pday_beg(1), hday_beg(1)/)
            vegn%Crop%crop_calendars(:,2,SECOND_SEASON,ipref,RAINFED) = (/ pday_beg(2), hday_beg(2)/)
            vegn%Crop%crop_calendars(:,3,MAIN_SEASON,  ipref,RAINFED) = (/ pday_end(1), hday_end(1)/)
            vegn%Crop%crop_calendars(:,3,SECOND_SEASON,ipref,RAINFED) = (/ pday_end(2), hday_end(2)/)
          endif
        enddo water_loop
      endif
    enddo crop_loop
    if(vegn%Crop%status == IDLE .or. vegn%Crop%status == ACTIVE_ON_LM3_SCHEDULE) then
      call crop_selection(vegn)
    endif
 endif ! if(new_month)
 call send_tile_data(id_T_ave, vegn%Crop%tc_av_climate,diag)
 call send_tile_data(id_P_ave, vegn%Crop%precip_av_climate,diag)
 call send_tile_data(id_status, real(vegn%Crop%status),diag)
 call send_tile_data(id_chosen_calendars(1), real(vegn%Crop%chosen_calendars(:,1)),diag)
 call send_tile_data(id_chosen_calendars(2), real(vegn%Crop%chosen_calendars(:,2)),diag)
 call send_tile_data(id_chosen_crop, real(vegn%Crop%chosen_crop),diag)
 do iwater=1,num_crop_water_sources
   do iseason=1,num_crop_seasons
     do iph=1,2
       call send_tile_data(id_crop_calendars(iph,iseason,iwater), real(vegn%Crop%crop_calendars(iph,1,iseason,:,iwater)),diag)
     enddo
   enddo
 enddo

 end subroutine compute_crop_calendars
!======================================================================================================================================================
 subroutine crop_selection(vegn)
 type(vegn_tile_type), intent(inout) :: vegn
 integer :: iwater, ipref_1, ipref_2, iseason, iperiod, pday1, hday1, pday2, hday2, cbeg
 logical :: found_a_1st_crop, found_a_2nd_crop
 character(len=12) :: crop_name_tmp
 character(len=256) :: text

 if(trim(water_source) == 'rainfed') then
   iwater = RAINFED
 else if(trim(water_source) == 'irrigated') then
   iwater = IRRIGATED
 else
   call error_mesg('subroutine crop_selection', trim(water_source)//' is an invalid value of water_source', FATAL)
 endif

! Terminology:
! "CCA" refers to the crop calendar algorithm.
! "CSA" refers to the crop selection algorithm. CSA is executed within this subroutine.
! "main season" and "second season" refer to the two seasons of a single crop type from the CCA.
! "1st" and "2nd" refer to the first and second cropping periods selected by the CSA.
! "1st" and "2nd" do not have to be the same crop type and "2nd" does not have to be one of the second season crops from the CCA.

 found_a_1st_crop = .false. ! Will later be set to .true. if conditions are suitable for any of the crops.
 crop_loop_1: do ipref_1=1,num_crop_types
   if(vegn%Crop%crop_calendars(1,1,1,ipref_1,iwater) /= NO_DATE) then
     ! If there is an optimal planting date for the main season of this crop then the
     ! crop calendar algorithm has determined that the climate is suitable for this
     ! crop and this will be the 1st cropping period of the crop selection algorithm.
     vegn%Crop%chosen_calendars(:,1) = vegn%Crop%crop_calendars(:,1,1,ipref_1,iwater)
     vegn%Crop%chosen_crop(1) = vegn%Crop%potential_crop(ipref_1)
     found_a_1st_crop = .true.
     exit crop_loop_1
   endif
 enddo crop_loop_1
 if(.not.found_a_1st_crop) then
   ! The crop calendar algorithm has determined that the climate is
   ! not suitable for any of the crops in the crop priority list.
   vegn%Crop%chosen_calendars(:,:) = NO_DATE
   vegn%Crop%chosen_crop(:) = NO_CROP
   return
 endif

! Choose the second season of the same crop as the 2nd crop if a second season exists. 
 pday2 = vegn%Crop%crop_calendars(1,1,2,ipref_1,iwater)
 hday2 = vegn%Crop%crop_calendars(2,1,2,ipref_1,iwater)
 if(pday2 /= NO_DATE .and. hday2 /= NO_DATE) then
   ! A second crop of the same crop type can be grown
   vegn%Crop%chosen_calendars(:,2) = (/pday2,hday2/)
   vegn%Crop%chosen_crop(2) = vegn%Crop%potential_crop(ipref_1)
   return
 endif

 ! A 2nd crop of the same type as the 1st crop cannot be grown. Look for a different crop type.
 ! It is possible to grow a 2nd crop only if its growing period does not overlap with that of the 1st crop.
 found_a_2nd_crop = .false.
 pday1 = vegn%Crop%chosen_calendars(1,1)
 hday1 = vegn%Crop%chosen_calendars(2,1)
 crop_loop_2: do ipref_2=ipref_1+1,num_crop_types
   if(vegn%Crop%potential_crop(ipref_2) == NO_CROP) exit crop_loop_2
   crop_name_tmp = crop_name(vegn%Crop%potential_crop(ipref_2))
   season_loop: do iseason=1,num_crop_seasons
     period_loop: do iperiod=1,num_crop_periods
       pday2 = vegn%Crop%crop_calendars(1,iperiod,iseason,ipref_2,iwater)
       hday2 = vegn%Crop%crop_calendars(2,iperiod,iseason,ipref_2,iwater)
       if(no_overlap(pday1, hday1, pday2, hday2)) then
         ! This growing period does not overlap with the growing period of the main crop. Make this the 2nd growing period.
         vegn%Crop%chosen_calendars(1,2) = pday2
         vegn%Crop%chosen_calendars(2,2) = hday2
         vegn%Crop%chosen_crop(2) = vegn%Crop%potential_crop(ipref_2)
         found_a_2nd_crop = .true.
         return
       endif
     enddo period_loop
   enddo season_loop
 enddo crop_loop_2
 if(.not.found_a_2nd_crop) then
   vegn%Crop%chosen_calendars(1,2) = NO_DATE
   vegn%Crop%chosen_calendars(2,2) = NO_DATE
   vegn%Crop%chosen_crop(2) = NO_CROP
 endif
 end subroutine crop_selection
!======================================================================================================================================================
 logical function no_overlap(tbeg1, tend1, tbeg2, tend2)
! returns .true. if there is no overlap between the time intervals tbeg1 to tend1 and tbeg2 to tend2
 integer, intent(in) :: tbeg1, tend1, tbeg2, tend2
 logical :: within_t1(365)
 integer :: GP1, GP2, dd, doy

 if(any((/tbeg1,tend1,tbeg2,tend2/) == (/NO_DATE,NO_DATE,NO_DATE,NO_DATE/))) then
   no_overlap = .false.
   return
 endif
 if(tend1 > tbeg1) then
   GP1 = tend1 - tbeg1
 else
   GP1 = tend1 - tbeg1 + 365
 endif
 within_t1 = .false.
 do dd=tbeg1,tbeg1+GP1
   doy = modulo_no_zero(dd,365)
   within_t1(doy) = .true.
 enddo
 if(tend2 > tbeg2) then
   GP2 = tend2 - tbeg2
 else
   GP2 = tend2 - tbeg2 + 365
 endif
 no_overlap = .true.
 do dd=tbeg2,tbeg2+GP2
   doy = modulo_no_zero(dd,365)
   if(within_t1(doy)) then
     no_overlap = .false.
     return
   endif
 enddo
 end function no_overlap
!======================================================================================================================================================
 subroutine CCA_Maize_Soybean_Rice(L, vegn, water, GP, central_T, variance_T, central_P, variance_P, &
                     central_D, variance_D, SI_crit, pday, pday_beg, pday_end, hday, hday_beg, hday_end)
 integer, intent(in) :: L
 type(vegn_tile_type), intent(in) :: vegn
 character(len=*), intent(in) :: water
 integer, intent(in) :: GP
 real, intent(in) :: central_T(0:), variance_T(0:), central_P(0:), variance_P(0:), central_D(0:), variance_D(0:), SI_crit
 integer, dimension(num_crop_seasons), intent(out) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end
 integer :: k, km, kp, daybeg, mths_after, doy, k_at_SI_min, ktest, day, iseason, max_range_length
 real :: Temp, Prec, annual_SI_min, dlen
 real, dimension(num_test_days) :: TSI, DSI, PSI, SI
!---------------------------------------------------------------------
 if(trim(water) /= 'irrigated' .and. trim(water) /= 'rainfed') then
   call error_mesg('CCA_Maize_Soybean_Rice ERROR: '//trim(water), 'is not a valid value of water', FATAL)
 endif
 k_loop_1: do k=1,num_test_days ! Compute the suitability index at 5 day intervals, starting with Jan 1
   TSI(k) = 0.0
   PSI(k) = 0.0
   DSI(k) = 0.0
   daybeg = 5*k
   mths_loop: do mths_after=0,num_m
     day = daybeg + 30*mths_after
     doy = modulo_no_zero(day,365)
     km = doy/5
     kp = km+1
     Temp = interp_between_mid_mths(doy, vegn%Crop%T_mid_mth)
     TSI(k) = TSI(k) + (Temp - central_T(mths_after))**2/variance_T(mths_after)
     if(mths_after < 4) then
       ! Month 4 is not tested for precip or day length
       Prec = interp_between_mid_mths(doy, vegn%Crop%P_mid_mth)
       if(trim(water) == 'irrigated') Prec = max(Prec,central_P(mths_after))
       PSI(k) = PSI(k) + (Prec - central_P(mths_after))**2/variance_P(mths_after)
       dlen = .2*((doy-5*km)*day_length(kp,L) + (5*kp-doy)*day_length(km,L))
       DSI(k) = DSI(k) + (dlen - central_D(mths_after))**2/variance_D(mths_after)
     endif
   enddo mths_loop
   SI(k) = TSI(k) + DSI(k) + PSI(k)
 enddo k_loop_1

 annual_SI_min = HUGE(1.0)
 k_loop_2: do k=1,num_test_days ! Find the k index of the annual minimum of SI
   if(SI(k) < annual_SI_min) then
     annual_SI_min = SI(k)
     k_at_SI_min = k
   endif
 enddo k_loop_2

 if(annual_SI_min < SI_crit) then
   pday(1) = 5*k_at_SI_min
 else
   pday(:) = NO_DATE
   pday_beg(:) = NO_DATE
   pday_end(:) = NO_DATE
   hday(:) = NO_DATE
   hday_beg(:) = NO_DATE
   hday_end(:) = NO_DATE
   return
 endif

 ktest = modulo_no_zero(k_at_SI_min+36,num_test_days)
 if(SI(ktest) < SI_crit) then
   pday(2) = 5*ktest
   max_range_length = 179 - GP
   max_range_length = 5*(max_range_length/5) ! Round down to the nearest multiple of 5
   do iseason=1,num_crop_seasons
     call find_planting_date_range(pday(iseason), SI_crit, SI, max_range_length, pday_beg(iseason), pday_end(iseason))
   enddo
 else
   call find_planting_date_range(pday(1), SI_crit, SI, 360, pday_beg(1), pday_end(1))
   pday(2) = NO_DATE
   pday_beg(2) = NO_DATE
   pday_end(2) = NO_DATE
 endif

 do iseason=1,num_crop_seasons
   if(pday(iseason) /= NO_DATE) then
     hday_beg(iseason) = modulo_no_zero(pday_beg(iseason) + GP, 365)
     hday(iseason)     = modulo_no_zero(pday(iseason)     + GP, 365)
     hday_end(iseason) = modulo_no_zero(pday_end(iseason) + GP, 365)
   else
     hday_beg(iseason) = NO_DATE
     hday(iseason)     = NO_DATE
     hday_end(iseason) = NO_DATE
   endif
 enddo

 end subroutine CCA_Maize_Soybean_Rice
!======================================================================================================================================================
 subroutine CCA_Wheat(L, vegn, Wtype, central_T, variance_T, central_P, variance_P, central_D, variance_D, & ! intent(in)
                      SI_crit, max_planting_SI_SW, Tbase, aPTTtH_range, & ! intent(in)
                      length_of_vernalization_period, max_T_for_vernalization, min_planting_T, & ! intent(in)
                      pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)

! 1. Compute dates of accumulated photo-thermal time at intervals of 200 units, from zero to 800, for each
! candidate Optimal Planting Date (OPD) starting with Jan 5 and at five day intervals throughout the year.

! 2. If the accumulated photo-thermal time does not exceed 800 units starting from any date then the climate is deemed unsuitable for wheat.

! 3. Compute suitability index using climatic conditions at intervals of 200 units of accumulated photo-thermal time
! The suitability index is specific to the water source and variety:
! irrigated winter wheat, rainfed winter wheat, irrigated spring wheat, rainfed spring wheat

! 4. If the suitability index for the specific type of wheat being tested exceeds the critical
! value at all tested dates throughout the year then the climate is deemed unsuitable.

! 5. Reduce the candidate OPDs to those for which the suitability index is below the critical value.
! The corresponding harvest dates are the dates when the accumulated photo-thermal time reaches 837 units or the maximum, starting from the candidate OPD.

! 6. Reduce the candidate OPDs to those for which the temperature never drops below -7°C before the corresponding harvest date.

! 7. Reduce the candidate OPDs to those which are warmer than 5°C.

! 8. For winter wheat: Reduce the candidate OPDs to those for which temperature drops below 7°C for at least 40 days between the planting and harvest dates.
! For spring wheat: Reduce the candidate OPDs to those for which temperature remains above 5°C between the planting and harvest dates.

! 9. For winter wheat: The predicted OPD is the date of minimum suitability index among the remaining candidate OPDs.
! The predicted harvest date is the date when the accumulated photo-thermal time reaches 837 units or the maximum, starting from the candidate OPD.
! For spring wheat: The predicted OPD is the date of minimum suitability index among the remaining candidate OPDs if the minimum is between 3.25 and 9.0
! or, if the minimum is below 3.25, the date prior to the the date of the minimum when it reaches 3.25

!10. The range of suitable dates includes all contiguous dates having a suitability index below critial before and after the OPD.
!======================================================================================================================================================
 integer, intent(in) :: L ! index of grid cell which contains this tile
 type(vegn_tile_type), intent(inout) :: vegn
 character(len=2), intent(in) :: Wtype ! 'SW' or 'WW'
 real, intent(in) :: central_T(0:), variance_T(0:), central_P(0:), variance_P(0:), central_D(0:), variance_D(0:) ! At intervals of 200 aPTT units after planting. (0) is planting day.
 real, intent(in) :: SI_crit, max_planting_SI_SW, Tbase, aPTTtH_range(2)
 integer, intent(in) :: length_of_vernalization_period
 real, intent(in) :: max_T_for_vernalization, min_planting_T
 integer, dimension(num_water), intent(out) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end

 character(len=12) :: chwater(num_water)
 integer :: k, k2, k2m, k2p, km, kp, daybeg, crossing_point, kautumn, k_of_ann_SI_min, kk, kkp, k_of_ann_SI_max, doy, iwater
 real :: Temp, Prec, TSI_test, DSI_test, dlen
 real :: PSI_test(2), SI_test(2) ! first element for irrigated, second for rainfed
 real :: annual_SI_max, annual_SI_min
 real :: SI(num_test_days,2) ! Suitability Index. Computed at 5 day intervals from Jan 5 to Dec 31.
 integer :: crossing_day_400(num_test_days) ! Date at which accumulated photo-thermal time (aPTT) since planting reaches 400 units
 integer :: crossing_days(0:num_m)
 integer :: hday_list(num_test_days) ! Remember the values for each test date then choose the one that corresponds to the annual minimum suitability index.
 real :: aPTTtH_list(num_test_days)  ! Remember the values for each test date then choose the one that corresponds to the annual minimum suitability index.
 logical :: passes_other_criteria

 chwater(1) = 'irrigated '//Wtype
 chwater(2) = 'rainfed '//Wtype
 pday = NO_DATE
 hday = NO_DATE
 hday_list = NO_DATE
 aPTTtH_list = 0.0
 ! compute SI at 5 day intervals from Jan 5 to Dec 31. SI_test is used for this.
 ! SI_test is not loaded into SI unless it passes the suitability test.
 k_loop_1: do k=1,num_test_days
   daybeg = 5*k
! Steps 1, 2, 6 and the harvest day of Step 5 are all handled within subroutine days_of_aPTT_crossings
   call days_of_aPTT_crossings(L, daybeg, num_m, Tbase, aPTT_interval, aPTTtH_range, vegn%Crop%T_mid_mth, & ! intent(in)
                       hday_list(k), aPTTtH_list(k), crossing_days) ! intent(out)
   if(any(crossing_days(:) == (/NO_DATE,NO_DATE,NO_DATE,NO_DATE,NO_DATE/))) then
     SI(k,:) = unsuitable ! Steps 2, 6 and planting day of Step 5
     cycle k_loop_1 ! cycle k loop if aPTT never reaches 800 or if the temperature drops below -7°C before 800 units of aPTT is reached.
   endif
   crossing_day_400(k) = crossing_days(2)
   SI_test = 0.0
   do crossing_point=0,num_m ! Step 3
     doy = crossing_days(crossing_point)
     Temp = interp_between_mid_mths(doy, vegn%Crop%T_mid_mth)
     TSI_test = (Temp - central_T(crossing_point))**2/variance_T(crossing_point)
     Prec = interp_between_mid_mths(doy, vegn%Crop%P_mid_mth)
     PSI_test(1) = (max(Prec,central_P(crossing_point)) - central_P(crossing_point))**2/variance_P(crossing_point)
     PSI_test(2) = (Prec - central_P(crossing_point))**2/variance_P(crossing_point)
     km = doy/5
     kp = km+1 ! modulo_no_zero is not used here because kp is used to index day_length, which is dimensioned (0:num_test_days+1)
     dlen = .2*((doy-5*km)*day_length(kp,L) + (5*kp-doy)*day_length(km,L))
     DSI_test = (dlen - central_D(crossing_point))**2/variance_D(crossing_point)
     SI_test(1) = SI_test(1) + TSI_test + DSI_test + PSI_test(1)
     SI_test(2) = SI_test(2) + TSI_test + DSI_test + PSI_test(2)
     if(SI_test(1) > SI_crit) then
       SI(k,:) = unsuitable ! If SI_test(1) exceeds critical, then so does SI_test(2)
       cycle k_loop_1
     endif
   enddo ! do crossing_point=0,num_m
   SI(k,1) = SI_test(1) ! Conditions are suitable for planting irrigated Wheat on day of the year 5*k, provided it passes the tests in k_loop_2 and k_loop_3
   if(SI_test(2) > SI_crit) then
     SI(k,2) = unsuitable
   else
     SI(k,2) = SI_test(2) ! Conditions are suitable for planting rainfed Wheat on day of the year 5*k, provided it passes the tests in k_loop_2 and k_loop_3
   endif
 enddo k_loop_1

 water_loop: do iwater=1,num_water
   passes_other_criteria = .true.
   k_loop_2: do k=1,num_test_days ! Step 8: check that vernalization is possible for winter wheat and that temperature remains above 5C for spring wheat.
     if(SI(k,iwater)==unsuitable) cycle k_loop_2 ! Step 4 If the suitability index for the specific type of wheat being tested exceeds the critical value at
                                                 ! all tested dates throughout the year then the tests within k_loop_2 are not necessary and will be skipped.
     if(Wtype == 'SW') then
       ! If the temperature drops below 5°C during the growing period then flag it as unsuitable for planting.
       if(T_goes_below_5C_during_GP(5*k, hday_list(k), vegn%Crop%T_mid_mth)) then
         SI(k,iwater) = unsuitable ! Step 8
         passes_other_criteria = .false.
       endif
     endif
     if(Wtype == 'WW') then
       if(.not.vernalization_is_possible(5*k, crossing_day_400(k), vegn%Crop%T_mid_mth, length_of_vernalization_period, max_T_for_vernalization)) then
         SI(k,iwater) = unsuitable ! Step 8
         passes_other_criteria = .false.
       endif
     endif
   enddo k_loop_2

   k_loop_3: do k=1,num_test_days ! Do not plant when the temperature is below min_planting_T (default value is 5°C)
     if(SI(k,iwater) == unsuitable) cycle k_loop_3
     Temp = interp_between_mid_mths(5*k, vegn%Crop%T_mid_mth)
     if(Temp < min_planting_T) then
       SI(k,iwater) = unsuitable ! Step 7
       passes_other_criteria = .false.
     endif
   enddo k_loop_3

   annual_SI_min = unsuitable
   pday(iwater) = NO_DATE
   hday(iwater) = NO_DATE
   k_of_ann_SI_min = index_of_annual_SI_min(SI(:,iwater)) ! returns -1 if unsuitable all year
   k_of_ann_SI_max = index_of_annual_SI_max(SI(:,iwater)) ! never returns -1, returns the index of an unsuitable date if there are any.
   if(k_of_ann_SI_min > 0) then ! compute predicted planting and harvest days. Both remain zero if a planting date cannot be found.
     if(SI(k_of_ann_SI_min,iwater) /= unsuitable) then
       if(Wtype == 'SW') then
         annual_SI_min = SI(k_of_ann_SI_min,iwater)
         if(annual_SI_min > max_planting_SI_SW) then
           pday(iwater) = 5*k_of_ann_SI_min ! Step 9 If the annual minimum SI is above max_planting_SI_SW, then the predicted optimal planting date is the date of the minimum.
           hday(iwater) = hday_list(k_of_ann_SI_min) ! Step 9
         else
           k_loop_4: do k=k_of_ann_SI_min-1,k_of_ann_SI_min-72,-1
             ! Go back in time until SI reaches max_planting_SI_SW or until the end of the suitable period is reached.
             kk = modulo_no_zero(k, num_test_days)
             if(SI(kk,iwater) > max_planting_SI_SW .or. SI(kk,iwater) == unsuitable) then
               kkp = modulo_no_zero(k+1,num_test_days)
               pday(iwater) = 5*kkp ! Step 9 If the annual minimum SI is below 3.25 then the predicted optimal
                                    ! planting date is the date before the minimum when SI first drops below 3.25
               hday(iwater) = hday_list(kkp)
               exit k_loop_4
             endif
           enddo k_loop_4
         endif
       endif
       if(Wtype == 'WW') then
         annual_SI_min = SI(k_of_ann_SI_min,iwater)
         pday(iwater) = 5*k_of_ann_SI_min
         hday(iwater) = hday_list(k_of_ann_SI_min)
       endif
     endif
   endif

   annual_SI_max = SI(k_of_ann_SI_max,iwater)
   if(annual_SI_min > SI_crit) then ! Step 10 Find the range of dates over which conditions are suitable for planting
     ! SI is above SI_crit all year
     pday_beg(iwater) = NO_DATE
     pday_end(iwater) = NO_DATE
     hday_beg(iwater) = NO_DATE
     hday_end(iwater) = NO_DATE
   else if(annual_SI_max < SI_crit) then
     ! SI is below SI_crit all year
     pday_beg(iwater) = 5
     pday_end(iwater) = 365
     hday_beg(iwater) = 5
     hday_end(iwater) = 365
   else
     k_loop_5: do k=k_of_ann_SI_min-1,k_of_ann_SI_min-num_test_days+1,-1 ! Go back in time to find the first date where SI < SI_crit
       k2 = modulo_no_zero(k,num_test_days)
       if(SI(k2,iwater) > SI_crit) then
         k2p = modulo_no_zero(k2+1,num_test_days)
         pday_beg(iwater) = 5*k2p
         hday_beg(iwater) = hday_list(k2p)
         exit k_loop_5
       endif
     enddo k_loop_5
     k_loop_6: do k=k_of_ann_SI_min+1,k_of_ann_SI_min+num_test_days-1 ! Go forward in time to find the last date where SI < SI_crit
       k2 = modulo_no_zero(k,num_test_days)
       if(SI(k2,iwater) > SI_crit) then
         k2m = modulo_no_zero(k2-1,num_test_days)
         pday_end(iwater) = 5*k2m
         hday_end(iwater) = hday_list(k2m)
         exit k_loop_6
       endif
     enddo k_loop_6
   endif
 enddo water_loop

 end subroutine CCA_Wheat
!======================================================================================================================================================
 function vernalization_is_possible(pday, day_aPPT_400, T_mid_mth, length_of_vernalization_period, max_T_for_vernalization) result(It_is)
 integer, intent(in) :: pday, day_aPPT_400
 real, intent(in) :: T_mid_mth(12)
 integer, intent(in) :: length_of_vernalization_period
 real, intent(in) :: max_T_for_vernalization
 integer :: day, num_cold_days, day400
 logical :: It_is
 real :: Temp

 It_is = .false.
 if(day_aPPT_400 < pday) then
   day400 = day_aPPT_400 + 365
 else
   day400 = day_aPPT_400
 endif
 num_cold_days = 0
 day_loop: do day=pday,day400
   Temp = interp_between_mid_mths(day, T_mid_mth)
   if(Temp < max_T_for_vernalization) then
     num_cold_days = num_cold_days + 1
     if(num_cold_days > length_of_vernalization_period) then
       It_is = .true.
       exit day_loop
     endif
   endif
 enddo day_loop
 end function vernalization_is_possible
!======================================================================================================================================================
 function T_goes_below_5C_during_GP(pday, hday, T_mid_mth) result(It_does)
 integer, intent(in) :: pday, hday
 real, intent(in) :: T_mid_mth(12)
 integer :: day, hdayy
 logical :: It_does
 real :: Temp

 It_does = .false.
 if(hday < pday) then
   hdayy = hday + 365
 else
   hdayy = hday
 endif
 day_loop: do day=pday,hdayy
   Temp = interp_between_mid_mths(day, T_mid_mth)
   if(Temp < 5.0+TFREEZE) then
     It_does = .true.
     exit day_loop
   endif
 enddo day_loop
 end function T_goes_below_5C_during_GP
!======================================================================================================================================================
 function index_of_annual_SI_min(SI) result(k_of_ann_SI_min)
 real, intent(in) :: SI(num_test_days)
 integer :: k_of_ann_SI_min, k
 real :: ann_min

 k_of_ann_SI_min = -1
 ann_min = unsuitable
 do k=1,num_test_days
   if(SI(k) == unsuitable) cycle
   if(SI(k) < ann_min) then
     k_of_ann_SI_min = k
     ann_min = SI(k)
   endif
 enddo
 end function index_of_annual_SI_min
!======================================================================================================================================================
 function index_of_annual_SI_max(SI) result(k_of_ann_SI_max)
 real, intent(in) :: SI(num_test_days)
 integer :: k_of_ann_SI_max, k
 real :: ann_max

 ann_max = 0.0
 do k=1,num_test_days
   if(SI(k) > ann_max) then
     k_of_ann_SI_max = k
     ann_max = SI(k)
   endif
 enddo
 end function index_of_annual_SI_max
!======================================================================================================================================================
 subroutine days_of_aPTT_crossings(L, daybeg, num_m, Tbase, aPTT_interval, aPTTtH_range, T_mid_mth, & ! intent(in)
                                   harvestday, harvest_aPTT, crossing_days) ! intent(out)
 integer, intent(in) :: L, daybeg, num_m
 real, intent(in) :: Tbase, aPTT_interval, aPTTtH_range(2)
 real, intent(in) :: T_mid_mth(12)

 ! crossing_days(m) = day of year when aPTT reaches m*aPTT_interval
 ! crossing_days(m) = zero if m*aPTT_interval is never reached or if temperature drops below -7C before aPTT reaches a value of num_m*aPTT_interval
 ! The date being tested is not suitable for planting either Spring or Winter Wheat if any of crossing_days(:) returned is zero

 integer, intent(out) :: harvestday
 real, intent(out) :: harvest_aPTT
 integer, intent(out) :: crossing_days(0:num_m)

 integer :: dd, m, km, kp, day_of_aPTT_max, doy_today, doy_tomorrow
 real :: aPTT_target, T_today, dlen, aPTT_max, aPTT(365)

 aPTT_target = aPTT_interval
 aPTT_max = 0.0
 day_of_aPTT_max = 0
 aPTT(:) = 0.0
 m = 0
 crossing_days(m) = modulo_no_zero(daybeg+1,365)
 crossing_days(1:num_m) = NO_DATE
  day_loop: do dd=daybeg,daybeg+364
   doy_today    = modulo_no_zero(dd,365)
   doy_tomorrow = modulo_no_zero(dd+1,365)
   if(aPTT(doy_today) > aPTTtH_range(2)) exit day_loop
   T_today = interp_between_mid_mths(doy_today, T_mid_mth)
   if(T_today < absolute_min_T_for_Wheat) exit day_loop
   km = doy_today/5
   kp = km + 1
   dlen = .2*((doy_today-5*km)*day_length(kp,L) + (5*kp-doy_today)*day_length(km,L))
   aPTT(doy_tomorrow) = aPTT(doy_today) + dlen*max(0.0,T_today-Tbase)
   if(aPTT(doy_tomorrow) > aPTT_max) then
     aPTT_max = aPTT(doy_tomorrow)
     day_of_aPTT_max = doy_tomorrow
   endif
   if(aPTT(doy_today) <= aPTT_target .and. aPTT(doy_tomorrow) > aPTT_target .and. m < num_m) then
     m = m + 1
     crossing_days(m) = doy_tomorrow
     aPTT_target = aPTT_target + aPTT_interval
   endif
 enddo day_loop
 if(aPTT_max > aPTTtH_range(1)) then
   harvestday = day_of_aPTT_max
   harvest_aPTT = aPTT_max
 else
   harvestday = NO_DATE
   harvest_aPTT = 0.0
 endif
 end subroutine days_of_aPTT_crossings
!======================================================================================================================================================
 subroutine find_planting_date_range(pday, SI_crit, SI, max_range_length, pday_beg, pday_end)
 integer, intent(in) :: pday
 real, intent(in) :: SI_crit, SI(num_test_days)
 integer, intent(in) :: max_range_length
 integer, intent(out) :: pday_beg, pday_end
 character(len=256) :: mesg
 integer :: k, k2, k2m, k2p, k_at_pday, krange

 k_at_pday = pday/5
 if(k_at_pday < 1 .or. k_at_pday > num_test_days) then
   mesg = 'ERROR1 in subroutine find_planting_date_range: invalid value of pday. pday=      '
   write(mesg(76:81),'(i6)') pday
   call error_mesg('find_planting_date_range in vegn_crop_mod',trim(mesg), FATAL)
 endif
 if(SI(k_at_pday) > SI_crit) then
   mesg = 'ERROR2 in subroutine find_planting_date_range: suitability index at planting day exceeds SI_crit'
   call error_mesg('find_planting_date_range in vegn_crop_mod',trim(mesg), FATAL)
 endif

 krange = nint(0.1*max_range_length)
 k2 = modulo_no_zero(k_at_pday-krange,num_test_days)
 pday_beg = 5*k2 ! overwritten below if an unsuitable date is found before k_at_pday - krange is reached.
 k_loop_1: do k=k_at_pday-1,k_at_pday-krange,-1 ! Go back in time to find the first date where SI < SI_crit
   k2 = modulo_no_zero(k,num_test_days)
   if(SI(k2) > SI_crit) then
     k2p = modulo_no_zero(k2+1,num_test_days) ! move forward one step to get first date where SI < SI_crit
     pday_beg = 5*k2p
     exit k_loop_1
   endif
 enddo k_loop_1

 k2 = modulo_no_zero(k_at_pday+krange,num_test_days)
 pday_end = 5*k2 ! overwritten below if an unsuitable date is found before k_at_pday + krange is reached.
 k_loop_2: do k=k_at_pday+1,k_at_pday+krange ! Go forward in time to find the last date where SI < SI_crit
   k2 = modulo_no_zero(k,num_test_days)
   if(SI(k2) > SI_crit) then
     k2m = modulo_no_zero(k2-1,num_test_days) ! back up one step to get last date where SI < SI_crit
     pday_end = 5*k2m
     exit k_loop_2
   endif
 enddo k_loop_2
 end subroutine find_planting_date_range
!======================================================================================================================================================
 real function interp_between_mid_mths(day_of_year, midmonth_values)
 integer, intent(in) :: day_of_year
 real, intent(in) :: midmonth_values(12)
 real :: tmp(0:13), w0, w1, day_of_year_local
 integer :: mon

 day_of_year_local = modulo_no_zero(day_of_year,365)
 tmp( 0) = midmonth_values(12)
 tmp(13) = midmonth_values( 1)
 tmp(1:12) = midmonth_values
  do mon=12,0,-1
   if(day_of_year_local >= t_mid_month(mon)) exit
 enddo
! Note that t_mid_month is dimensioned (0:13) where t_mid_month(0) is negative because it is the middle of Dec
! of the previous year and t_mid_month(13) is > 365. because it is the middle of Jan of the following year.
 w0 = ( day_of_year_local - t_mid_month(mon)) / (t_mid_month(mon+1) - t_mid_month(mon))
 w1 = (t_mid_month(mon+1) - day_of_year_local) / (t_mid_month(mon+1) - t_mid_month(mon))
 interp_between_mid_mths = w1*tmp(mon) + w0*tmp(mon+1)
 end function interp_between_mid_mths
!======================================================================================================================================================
 subroutine read_crop_namelist
 integer :: outunit, io, ierr

  read(input_nml_file, nml=vegn_crop_nml, iostat=io)
  ierr = check_nml_error(io, 'vegn_crop_nml')
  outunit = stdlog()
  write(outunit, nml=vegn_crop_nml)
 end subroutine read_crop_namelist
!======================================================================================================================================================
 subroutine vegn_crop_init(id_ug)
 integer, intent(in) :: id_ug
 type(land_restart_type) :: restart
 type(FmsNetcdfFile_t) :: fileobj
 logical :: restart_exists, used, init_clim_exists
 real :: Dmm, Dm0, Dmp
 real, dimension(12,12) :: X
 real :: max_frac
 integer :: dummyi, L, m, k, k_of_max_frac, dom, doy ! L = index of grid cell which contains this tile
 type(land_tile_enum_type) :: ce
 type(land_tile_type), pointer :: tile
 character(len=3) :: month_name(12) = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)
 integer :: day_ae, month_ae, year_ae, hour_ae, minute_ae, second_ae, ierr, icrop, iwater, iseason
 real, allocatable, dimension(:,:) :: MIRCA_crop_frac, crop_frac_tmp, Tclim, Pclim
 integer, allocatable, dimension(:,:) :: potential_crop
 character(len=256) :: infile, text

 call read_crop_namelist
 do iwater=1,num_crop_water_sources
 do icrop=1,num_crop_types
 do iseason=1,num_crop_seasons
   restart_fieldname(iseason,icrop,iwater) = trim(water_source_name(iwater))//'_'//trim(crop_name(icrop))//'_'//trim(season_name(iseason))
 enddo
 enddo
 enddo
 infile = ''
! Read the restart data
 infile = 'INPUT/'//trim(restart_file_name)
 call open_land_restart(restart,trim(infile),restart_exists)
 if(restart_exists) then
   call error_mesg('vegn_crop_init', 'reading NetCDF restart', NOTE)
   call get_tile_data(restart, 'tc_av_climate', 'month', vegn_tc_av_climate_ptr)
   call get_tile_data(restart, 'precip_av_climate', 'month', vegn_precip_av_climate_ptr)
   call get_tile_data(restart, 'T_mid_mth', 'month', vegn_T_mid_mth_ptr)
   call get_tile_data(restart, 'P_mid_mth', 'month', vegn_P_mid_mth_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,1,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_1_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,1,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_1_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,2,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_2_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,2,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_2_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,3,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_3_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,3,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_3_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,4,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_4_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,4,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_4_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,5,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_5_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,5,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_5_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,1,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_1_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,1,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_1_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,2,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_2_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,2,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_2_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,3,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_3_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,3,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_3_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,4,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_4_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,4,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_4_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,5,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_5_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,5,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_5_2_ptr)
   call get_int_tile_data(restart, 'chosen_calendars', 'plant_harvest', 'crop_seasons', vegn_chosen_calendars_ptr)
   call get_int_tile_data(restart, 'potential_crop', 'crop_types',vegn_potential_crop_ptr)
   call get_int_tile_data(restart, 'chosen_crop', 'crop_seasons', vegn_chosen_crop_ptr)
   call get_int_tile_data(restart, 'status', vegn_status_ptr)
 else
   call error_mesg('vegn_crop_init', 'cold starting vegn_crop_mod', NOTE)
   allocate(Tclim(lnd%ls:lnd%le,12), Pclim(lnd%ls:lnd%le,12))
   init_clim_exists = open_file(fileobj, "INPUT/initial_climatology.nc", "read")
   if (.not. init_clim_exists) then
     call error_mesg('soil_init', 'INPUT/initial_climatology.nc does not exist', FATAL)
   endif
   call read_field(fileobj, 'tc_av_climate',     Tclim, 'nearest')
   call read_field(fileobj, 'precip_av_climate', Pclim, 'nearest')
   call close_file(fileobj)
   ce = first_elmt(land_tile_map, ls=lnd%ls)
   do while(loop_over_tiles(ce,tile,L,k))
     if(.not.associated(tile%vegn)) cycle
     tile%vegn%Crop%tc_av_climate = Tclim(L,:)
     tile%vegn%Crop%precip_av_climate = Pclim(L,:)
     tile%vegn%Crop%T_mid_mth = Tclim(L,:)
     tile%vegn%Crop%P_mid_mth = Pclim(L,:)
     tile%vegn%Crop%crop_calendars   = NO_DATE
     tile%vegn%Crop%chosen_calendars = NO_DATE
     tile%vegn%Crop%status = IDLE
   enddo
   deallocate(Tclim, Pclim)

! Read the MIRCA crop fractions. The crop with the largest fraction becomes the potential_crop.
! These crop fractions are the sum of irrigated and rainfed fractions from the MIRCA2000 data.
   allocate(MIRCA_crop_frac(lnd%ls:lnd%le,num_crop_types), crop_frac_tmp(lnd%ls:lnd%le,num_crop_types))
   allocate(potential_crop(lnd%ls:lnd%le,num_crop_types))
   dummyi = 0
   call init_cover_field('single-tile', 'INPUT/crop_frac.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, (/ 1,-1,-1,-1,-1/), crop_frac_tmp)
   MIRCA_crop_frac(:,1) = crop_frac_tmp(:,1) ! MAIZE
   call init_cover_field('single-tile', 'INPUT/crop_frac.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, (/-1, 2,-1,-1,-1/), crop_frac_tmp)
   MIRCA_crop_frac(:,2) = crop_frac_tmp(:,2) ! SOYBEAN
   call init_cover_field('single-tile', 'INPUT/crop_frac.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, (/-1,-1, 3,-1,-1/), crop_frac_tmp)
   MIRCA_crop_frac(:,3) = crop_frac_tmp(:,3) ! RICE
   call init_cover_field('single-tile', 'INPUT/crop_frac.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, (/-1,-1,-1, 4,-1/), crop_frac_tmp)
   MIRCA_crop_frac(:,4) = crop_frac_tmp(:,4) ! SPRING_WHEAT
   call init_cover_field('single-tile', 'INPUT/crop_frac.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, (/-1,-1,-1,-1, 5/), crop_frac_tmp)
   MIRCA_crop_frac(:,5) = crop_frac_tmp(:,5) ! WINTER_WHEAT

   L_loop_2: do L=lnd%ls,lnd%le
     m_loop: do m=1,num_crop_types
       max_frac = 0.0
       k_of_max_frac = 0
       do k=1,num_crop_types
         if(MIRCA_crop_frac(L,k) > max_frac) then
           max_frac = MIRCA_crop_frac(L,k)
           k_of_max_frac = k
         endif
       enddo
       if(k_of_max_frac == 0) then
         potential_crop(L,m) = NO_CROP
       else
         potential_crop(L,m) = k_of_max_frac
         MIRCA_crop_frac(L,k_of_max_frac) = 0.0
       endif
     enddo m_loop
   enddo L_loop_2

   ce = first_elmt(land_tile_map, ls=lnd%ls)
   do while(loop_over_tiles(ce,tile,L,k))
     if(.not.associated(tile%vegn)) cycle
     tile%vegn%Crop%potential_crop = potential_crop(L,:)
     call crop_selection(tile%vegn)
   enddo
   deallocate(MIRCA_crop_frac, crop_frac_tmp, potential_crop)
 endif
!------------------------------------------------------------------------------------------------------------------------------------------------------
 ce = first_elmt(land_tile_map, ls=lnd%ls)
!------------------------------------------------------------------------------------------------------------------------------------------------------
! Convert units of means and variances
 central_T_Maize  = central_T_Maize_orig_units + TFREEZE
 central_P_Maize  = central_P_Maize_orig_units/SECONDS_PER_DAY
 variance_P_Maize = variance_P_Maize_orig_units/SECONDS_PER_DAY**2
 central_T_Soy    = central_T_Soy_orig_units + TFREEZE
 central_P_Soy    = central_P_Soy_orig_units/SECONDS_PER_DAY
 variance_P_Soy   = variance_P_Soy_orig_units/SECONDS_PER_DAY**2
 central_T_SW     = central_T_SW_orig_units + TFREEZE
 central_P_SW     = central_P_SW_orig_units/SECONDS_PER_DAY
 variance_P_SW    = variance_P_SW_orig_units/SECONDS_PER_DAY**2
 central_T_WW     = central_T_WW_orig_units + TFREEZE
 central_P_WW     = central_P_WW_orig_units/SECONDS_PER_DAY
 variance_P_WW    = variance_P_WW_orig_units/SECONDS_PER_DAY**2
 central_T_Rice   = central_T_Rice_orig_units + TFREEZE
 central_P_Rice   = central_P_Rice_orig_units/SECONDS_PER_DAY
 variance_P_Rice  = variance_P_Rice_orig_units/SECONDS_PER_DAY**2
 do icrop=1,num_crop_types
   if(icrop == MAIZE) then
     SI_crit(icrop) = SI_crit_Maize
     GP(icrop) = GP_Maize
     central_T(:,icrop)  = central_T_Maize
     central_P(:,icrop)  = central_P_Maize
     central_D(:,icrop)  = central_D_Maize
     variance_T(:,icrop) = variance_T_Maize
     variance_P(:,icrop) = variance_P_Maize
     variance_D(:,icrop) = variance_D_Maize
   else if(icrop == SOYBEAN) then
     SI_crit(icrop) = SI_crit_Soy
     GP(icrop) = GP_Soy
     central_T(:,icrop)  = central_T_Soy
     central_P(:,icrop)  = central_P_Soy
     central_D(:,icrop)  = central_D_Soy
     variance_T(:,icrop) = variance_T_Soy
     variance_P(:,icrop) = variance_P_Soy
     variance_D(:,icrop) = variance_D_Soy
   else if(icrop == RICE) then
     SI_crit(icrop) = SI_crit_Rice
     GP(icrop) = GP_Rice
     central_T(:,icrop)  = central_T_Rice
     central_P(:,icrop)  = central_P_Rice
     central_D(:,icrop)  = central_D_Rice
     variance_T(:,icrop) = variance_T_Rice
     variance_P(:,icrop) = variance_P_Rice
     variance_D(:,icrop) = variance_D_Rice
   else if(icrop == SPRING_WHEAT) then
     SI_crit(icrop) = SI_crit_SW
     central_T(:,icrop)  = central_T_SW
     central_P(:,icrop)  = central_P_SW
     central_D(:,icrop)  = central_D_SW
     variance_T(:,icrop) = variance_T_SW
     variance_P(:,icrop) = variance_P_SW
     variance_D(:,icrop) = variance_D_SW
   else if(icrop == WINTER_WHEAT) then
     SI_crit(icrop) = SI_crit_WW
     central_T(:,icrop)  = central_T_WW
     central_P(:,icrop)  = central_P_WW
     central_D(:,icrop)  = central_D_WW
     variance_T(:,icrop) = variance_T_WW
     variance_P(:,icrop) = variance_P_WW
     variance_D(:,icrop) = variance_D_WW
   endif
 enddo
!------------------------------------------------------------------------------------------------------------------------------------------------------
! compute coefficients used to compute mid-month climatological values of temperature and precip rates from monthly means
 X = 0.0
 do m=2,11
   Dmm = days_in_month(m-1)
   Dm0 = days_in_month(m)
   Dmp = days_in_month(m+1)
   X(m,m-1) = Dm0/(Dmm+Dm0)
   X(m,m ) = Dmm/(Dmm+Dm0) + Dmp/(Dm0+Dmp) + 2.0
   X(m,m+1) = Dm0/(Dm0+Dmp)
 enddo
 Dmm = days_in_month(12)
 Dm0 = days_in_month(1)
 Dmp = days_in_month(2)
 X(1,12) = Dm0/(Dmm+Dm0)
 X(1, 1) = Dmm/(Dmm+Dm0) + Dmp/(Dm0+Dmp) + 2.0
 X(1, 2) = Dm0/(Dm0+Dmp)
 Dmm = days_in_month(11)
 Dm0 = days_in_month(12)
 Dmp = days_in_month(1)
 X(12,11) = Dm0/(Dmm+Dm0)
 X(12,12) = Dmm/(Dmm+Dm0) + Dmp/(Dm0+Dmp) + 2.0
 X(12, 1) = Dm0/(Dm0+Dmp)
 X_ludcmp = X
 call ludcmp(X_ludcmp, indx_ludcmp, ierr)
 if(ierr /= 0) call error_mesg('vegn_crop_init in vegn_crop_mod','Error in subroutine ludcmp. Probably because matrix is singular.', FATAL)
!------------------------------------------------------------------------------------------------------------------------------------------------------
! compute day length at 5 day intervals
 call get_orbital_parameters(ecc, obliq, per)
 call get_ref_date_of_ae(day_ae, month_ae, year_ae, second_ae, minute_ae, hour_ae)
 autumnal_eq_ref = set_date(year_ae, month_ae, day_ae, hour_ae, minute_ae, second_ae)
 period_time_type = length_of_year()
 call orbit

 allocate(day_length(0:num_test_days+1,lnd%ls:lnd%le))

 ce = first_elmt(land_tile_map, ls=lnd%ls)
 do while(loop_over_tiles(ce,tile,L))
   do k=1,num_test_days
     call compute_day_length(year_ae, 5*k, lnd%ug_lat(L), day_length(k,L))
   enddo
!  extend day_length array one step beyond either end of the year to facilitate interpolation between 5 day intervals
   day_length(0,L) = day_length(num_test_days,L)
   day_length(num_test_days+1,L) = day_length(1,L)
 enddo
!------------------------------------------------------------------------------------------------------------------------------------------------------
! compute coefficients used by function interp_between_mid_mths
 t_mid_month(0) = -.5*days_in_month(12)
 t_mid_month(1) =  .5*days_in_month(1)
 do m=2,12
   t_mid_month(m) = t_mid_month(m-1) + .5*(days_in_month(m-1)+days_in_month(m))
 enddo
 t_mid_month(13) = t_mid_month(12) + .5*(days_in_month(12)+days_in_month(1))
!------------------------------------------------------------------------------------------------------------------------------------------------------
 call crop_diag_init(id_ug)
 allocate(potential_crop(lnd%ls:lnd%le,num_crop_types))
 ce = first_elmt(land_tile_map, ls=lnd%ls)
 do while(loop_over_tiles(ce,tile,L,k))
   if(.not.associated(tile%vegn)) cycle
   potential_crop(L,:) = tile%vegn%Crop%potential_crop
 enddo
 if(id_potential_crop > 0) used = send_data(id_potential_crop, real(potential_crop), lnd%time)
 deallocate(potential_crop)
 crop_mod_initialized = .TRUE.
 end subroutine vegn_crop_init
!======================================================================================================================================================
 subroutine compute_day_length(year_ae, doy, lat, daylen)
 integer, intent(in) :: year_ae, doy
 real, intent(in) :: lat
 real, intent(out) :: daylen

 type(time_type) :: Jan01_year_of_ae
 real :: otime, ang, sindec, cosdec, tandec, cos_half_day

 Jan01_year_of_ae = set_date(year_ae, 1, 1, 0, 0, 0)
 otime = orbital_time(Jan01_year_of_ae + set_time(0,doy) )
 ang = angle(otime)
 sindec = -sin(PI*obliq/180.)*sin(ang)
 cosdec = sqrt(1.0 - sindec**2)
 tandec = sindec/cosdec
 if(lat == -0.5*PI) then
   if(sindec > 0.0) then
     daylen = 0.0
   else
     daylen = 1.0
   endif
 else if(lat == 0.5*PI) then
   if(sindec > 0.0) then
     daylen = 1.0
   else
     daylen = 0.0
   endif
 else
   cos_half_day = -tan(lat)*tandec
   if(cos_half_day <= -1.0) then
     daylen = 1.0
   else if(cos_half_day >= 1.0) then
     daylen = 0.0
   else
     daylen = acos(cos_half_day)/PI
   endif
 endif

 end subroutine compute_day_length
!======================================================================================================================================================
 subroutine orbit
 integer :: n
 real :: d1, d2, d3, d4, d5, dt, norm

 ! Solves for orbital angle as a function of time via Runge-Kutta.
 ! The equation solved is the conservation of angular momentum.
 orb_angle(0) = 0.0
 dt = 2*PI/float(num_angles)
 norm = sqrt(1.0 - ecc**2)
 dt = dt*norm
 do n = 1,num_angles
   d1 = dt*r_inv_squared(orb_angle(n-1))
   d2 = dt*r_inv_squared(orb_angle(n-1)+0.5*d1)
   d3 = dt*r_inv_squared(orb_angle(n-1)+0.5*d2)
   d4 = dt*r_inv_squared(orb_angle(n-1)+d3)
   d5 = d1/6.0 + d2/3.0 + d3/3.0 + d4/6.0
   orb_angle(n) = orb_angle(n-1) + d5
 end do

 end subroutine orbit
!======================================================================================================================================================
 function r_inv_squared(ang)
 real, intent(in) :: ang
 real :: r, r_inv_squared, rad_per

 rad_per = PI*per/180.
 r = (1. - ecc**2)/(1.0 + ecc*cos(ang - rad_per))
 r_inv_squared = 1.0/r**2
 end function r_inv_squared
!======================================================================================================================================================
 function orbital_time(time) result(otime)
 type(time_type), intent(in) :: time
 real :: otime

 otime = real((time - autumnal_eq_ref)//period_time_type) ! What is the purpose of "real". Looks like it's not needed.
 otime = 2*PI*(otime - floor(otime))
 if(time < autumnal_eq_ref) otime = 2*PI - otime ! This is necessary because time_type is always positive.
                                                 ! Therefore, otime is positive when time is before autumnal_eq_ref.
 end function orbital_time
!======================================================================================================================================================
 function angle(otime) result(ang)
 real, intent(in) :: otime
 real :: ang, index_as_real, x
 integer :: index_as_int

 index_as_real = otime*float(num_angles)/(2*PI)
 index_as_int = floor(index_as_real)
 index_as_int = modulo(index_as_int,num_angles)
 x = index_as_real - index_as_int
 ang = (1.0 -x)*orb_angle(index_as_int) + x*orb_angle(index_as_int+1)
 ang = modulo(ang, 2*PI)
 end function angle
!======================================================================================================================================================
 subroutine crop_diag_init(id_ug)
 integer, intent(in) :: id_ug
 integer :: id_month, mth, id_crop_num, ical, id_season, icrop, iseason, iwater, iph
 character(len=256) :: diag_fieldname

 id_month    = diag_axis_init('month', (/(float(mth),mth=1,12)/),'none','Z','month of year')

 diag_fieldname = trim(crop_name(1))
 do icrop=2,num_crop_types
   diag_fieldname = trim(diag_fieldname)//', '//trim(crop_name(icrop))
 enddo
 id_crop_num = diag_axis_init('crop_num',(/(float(icrop),icrop=1,num_crop_types)/),'none','Z',trim(diag_fieldname))

 diag_fieldname = trim(season_name(1))
 do iseason=2,num_crop_seasons
   diag_fieldname = trim(diag_fieldname)//', '//trim(season_name(iseason))
 enddo
 id_season = diag_axis_init('crop_seasons',(/(float(iseason),iseason=1,num_crop_seasons)/),'none','Z',trim(diag_fieldname))

 call set_default_diag_filter('soil')
 id_T_ave = register_tiled_diag_field(module_name,'tc_av_climate',    (/id_ug,id_month/),lnd%time,'climatological monthly mean temperature','deg K', missing_value=-1.0)
 id_P_ave = register_tiled_diag_field(module_name,'precip_av_climate',(/id_ug,id_month/),lnd%time,'climatological monthly mean precip rate','Kg/s*m^2',missing_value=-1.0)
 do iwater=1,num_crop_water_sources
   do iseason=1,num_crop_seasons
     do iph=1,2
       if(iph == 1) then
         diag_fieldname = trim(water_source_name(iwater))//'_'//trim(season_name(iseason))//'_planting_date'
       else
         diag_fieldname = trim(water_source_name(iwater))//'_'//trim(season_name(iseason))//'_harvest_date'
       endif
       id_crop_calendars(iph,iseason,iwater) = register_tiled_diag_field(module_name,trim(diag_fieldname),(/id_ug,id_crop_num/),lnd%time,trim(diag_fieldname),missing_value= 0.0)
     enddo
   enddo
 enddo
 do iseason=1,num_crop_seasons
   diag_fieldname = trim(season_name(iseason))//'_planting_and_harvest_dates'
   id_chosen_calendars(iseason) = register_tiled_diag_field(module_name,trim(diag_fieldname),(/id_ug,id_season/),lnd%time,trim(diag_fieldname),missing_value= 0.0)
 enddo

 id_potential_crop  = register_static_field(module_name,'potential_crop',(/id_ug,id_crop_num/), 'crops sorted by area in the MIRCA2000 data set',missing_value=0.0)
 call set_default_diag_filter('crop')
 id_status = register_tiled_diag_field(module_name,'status',(/id_ug/),lnd%time,'IDLE = 0,ACTIVE_ON_CROP_SCHEDULE = 1,ACTIVE_ON_LM3_SCHEDULE = 2',missing_value= -1.0)
 id_chosen_crop = register_tiled_diag_field(module_name,'chosen_crop', (/id_ug,id_season/), lnd%time, 'number of chosen crop', missing_value=0.0)
 end subroutine crop_diag_init
 !======================================================================================================================================================
 subroutine save_crop_restart(tile_dim_length, timestamp)
 integer, intent(in) :: tile_dim_length
 character(*), intent(in) :: timestamp
 type(land_restart_type) :: restart
 character(len=256) :: filename
 integer :: idate, iwater, icrop, iseason, mth, iperiod, iph

 if(.not. crop_mod_initialized) call error_mesg('save_crop_restart','vegn_crop_init has not been called', FATAL)
 filename = 'RESTART/'//trim(timestamp)//trim(restart_file_name)
 call error_mesg('save_crop_restart', 'writing NetCDF restart "'//trim(filename)//'"', NOTE)
 call init_land_restart(restart, filename, vegn_tile_exists, tile_dim_length)
 call add_restart_axis(restart,'month', (/(float(mth),mth=1,12)/),.false.,longname='calendar month')
 call add_restart_axis(restart,'crop_cal_date',(/(float(idate),idate=1,num_crop_cal)/),.false.,longname='crop calendar dates as day of year')
 call add_restart_axis(restart,'crop_cal_periods',(/(float(iperiod),iperiod=1,num_crop_periods)/),.false.,longname='opt beg end')
 call add_restart_axis(restart,'plant_harvest',(/(float(iph),iph=1,2)/),.false.,longname='plant harvest')
 call add_restart_axis(restart,'crop_water',(/(float(iwater),iwater=1,num_crop_water_sources)/),.false.,longname='irrigated, rainfed')
 call add_restart_axis(restart,'crop_seasons',(/(float(iseason),iseason=1,num_crop_seasons)/),.false.,longname='crop season')
 call add_restart_axis(restart,'crop_types',(/(float(icrop),icrop=1,num_crop_types)/),.false.,longname='Maize, Soybean, Rice, Spring Wheat, Winter Wheat')
 call add_tile_data(restart,'tc_av_climate', 'month', vegn_tc_av_climate_ptr, 'climatological monthly average canopy air temperature','degK')
 call add_tile_data(restart,'precip_av_climate', 'month', vegn_precip_av_climate_ptr,'climatological monthly average precipitation rate','Kg/s*m^2')
 call add_tile_data(restart,'T_mid_mth', 'month', vegn_T_mid_mth_ptr, 'climatological average mid-month canopy air temperature','degK')
 call add_tile_data(restart,'P_mid_mth', 'month', vegn_P_mid_mth_ptr, 'climatological average mid-month precipitation rate','Kg/s*m^2')
 call add_int_tile_data(restart, trim(restart_fieldname(1,1,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_1_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,1,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_1_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,2,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_2_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,2,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_2_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,3,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_3_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,3,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_3_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,4,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_4_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,4,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_4_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,5,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_5_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,5,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_5_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,1,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_1_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,1,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_1_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,2,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_2_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,2,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_2_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,3,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_3_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,3,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_3_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,4,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_4_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,4,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_4_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,5,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_5_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,5,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_5_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')

 call add_int_tile_data(restart, 'chosen_calendars', 'plant_harvest', 'crop_seasons', vegn_chosen_calendars_ptr, 'popt, hopt, 1st and 2nd seasons','day of year')
 call add_int_tile_data(restart, 'potential_crop', 'crop_types',vegn_potential_crop_ptr, 'crop number of potential crop','crop number')
 call add_int_tile_data(restart, 'chosen_crop', 'crop_seasons', vegn_chosen_crop_ptr, 'crop number of chosen crop',      'crop number')
 call add_int_tile_data(restart, 'status',  vegn_status_ptr, 'IDLE = 0, ACTIVE_ON_CROP_SCHEDULE = 1, ACTIVE_ON_LM3_SCHEDULE = 2','dimensionless')
 call save_land_restart(restart)
 call free_land_restart(restart)
 end subroutine save_crop_restart
!======================================================================================================================================================
 subroutine vegn_crop_end()
  crop_mod_initialized = .FALSE.
 end subroutine vegn_crop_end
!======================================================================================================================================================
 logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
 end function vegn_tile_exists
!======================================================================================================================================================
 function modulo_no_zero(nn,cycle_len) result(nn_within_cycle)
 integer, intent(in) :: nn,cycle_len
 integer :: nn_within_cycle

 nn_within_cycle = modulo(nn,cycle_len)
 if(nn_within_cycle == 0) nn_within_cycle = cycle_len
 end function modulo_no_zero
!======================================================================================================================================================
 integer function get_crop_index_from_name(cropname) result(cropindex)
 character(len=*), intent(in) :: cropname
 integer :: icrop

 do icrop=0,num_crop_types
   if(trim(crop_name(icrop)) == trim(cropname)) then
     cropindex = icrop
     exit
   endif
 enddo
 end function get_crop_index_from_name
!======================================================================================================================================================
subroutine vegn_T_mid_mth_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%T_mid_mth(n)
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_P_mid_mth_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%P_mid_mth(n)
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_tc_av_climate_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%tc_av_climate(n)
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_precip_av_climate_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%precip_av_climate(n)
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_status_ptr(t,p)
 type(land_tile_type),pointer::t
 integer,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%status
 endif
 end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_1_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,1,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_1_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,1,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_2_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,2,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_2_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,2,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_3_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,3,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_3_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,3,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_4_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,4,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_4_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,4,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_5_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,5,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_5_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,5,1)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_1_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,1,2)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_1_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,1,2)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_2_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,2,2)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_2_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,2,2)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_3_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,3,2)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_3_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,3,2)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_4_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,4,2)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_4_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,4,2)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_1_5_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,5,2)
  endif
  end subroutine
!======================================================================================================================================================
 subroutine crop_calendar_2_5_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,5,2)
  endif
  end subroutine
!======================================================================================================================================================
subroutine vegn_potential_crop_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 integer,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%potential_crop(n)
 endif
end subroutine vegn_potential_crop_ptr
!======================================================================================================================================================
subroutine vegn_chosen_calendars_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%chosen_calendars(m,n)
  endif
end subroutine vegn_chosen_calendars_ptr
!======================================================================================================================================================
subroutine vegn_chosen_crop_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 integer,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%chosen_crop(n)
 endif
end subroutine vegn_chosen_crop_ptr
!======================================================================================================================================================
 end module vegn_crop_mod
