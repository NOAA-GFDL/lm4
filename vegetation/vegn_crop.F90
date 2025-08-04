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
! precipitation is adaquate. In all cases, planting and harvests dates are set to the value of
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
 use mpp_mod, only: input_nml_file, mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
 use fms_mod, only: error_mesg, NOTE, FATAL, check_nml_error, stdlog, CLOCK_FLAG_DEFAULT
 use time_manager_mod, only: time_type, set_date, get_date, operator(-), set_time, operator(+), length_of_year, operator(//), operator(<)
 use constants_mod, only: TFREEZE, SECONDS_PER_DAY, PI
 use land_tile_mod, only: land_tile_type, land_tile_enum_type, first_elmt, loop_over_tiles, land_tile_map
 use vegn_tile_mod, only: vegn_tile_type
 use vegn_data_mod, only: NO_DATE, NO_CROP, MAIN_SEASON, SECOND_SEASON, season_name, period_name, &
                          IRRIGATED_MAIZE, IRRIGATED_SOYBEAN, IRRIGATED_RICE, IRRIGATED_SPRING_WHEAT, IRRIGATED_WINTER_WHEAT, &
                          RAINFED_MAIZE,   RAINFED_SOYBEAN,   RAINFED_RICE,   RAINFED_SPRING_WHEAT,   RAINFED_WINTER_WHEAT, &
                          crop_name, num_crop_types, num_crop_cal, num_crop_seasons, num_crop_periods, &
                          landuse_name, water_source_name, landuse_longname, LU_IRRIG, LU_RAINF
 use land_data_mod, only: lnd
 use land_tile_io_mod, only: land_restart_type, init_land_restart, open_land_restart, save_land_restart, &
                             free_land_restart, add_restart_axis, add_tile_data, get_tile_data, field_exists, add_int_tile_data, get_int_tile_data
 use astronomy_mod, only: get_orbital_parameters, get_ref_date_of_ae
 use diag_manager_mod,  only: diag_axis_init, send_data, register_static_field, diag_field_add_attribute
 use land_tile_diag_mod,only: register_tiled_diag_field, send_tile_data, diag_buff_type, set_default_diag_filter
 use land_numerics_mod, only: ludcmp, lubksb
 use land_io_mod, only: init_cover_field, read_field
 use fms2_io_mod, only: close_file, FmsNetcdfFile_t, open_file
 use land_debug_mod, only: is_watch_cell     ! watchpoint_code
 use vegn_debug_crop_mod, only: debug_crop_2 ! watchpoint_code

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
 real, parameter :: unsuitable = 1000. ! suitablity index is set to unsuitable whenever it exceeds SI_crit. (Is this necessary?)
 real, parameter :: aPTT_interval = 200. ! Wheat suitablity is tested at intervals of aPTT_interval units of photo-thermal time.
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
 character(len=64) :: restart_fieldname(num_crop_seasons,num_crop_types)

 real :: t_mid_month(0:13)
! Note that t_mid_month is dimensioned (0:13) where t_mid_month(0) is negative because it is the middle of Dec
! of the previous year and t_mid_month(13) is > 365. because it is the middle of Jan of the following year.

 integer :: id_crop_calendars(2,num_crop_periods,num_crop_seasons,num_crop_types), id_chosen_calendars(num_crop_seasons)
 integer :: id_T_ave, id_P_ave, id_potential_crop, id_chosen_crop, cropclock

 real :: weight_climate=.10
 real :: max_planting_SI_SW = 9.75
 real :: Tbase_Wheat = 5.0 + TFREEZE
 integer :: length_of_vernalization_period = 40
 real :: max_T_for_vernalization = 8.0 + TFREEZE
 real :: min_T_GP_SW = 5.0 + TFREEZE
 real :: absolute_min_T_for_Wheat = -8.0 + TFREEZE
 real :: SI_crit(num_crop_types)
 integer :: GP(num_crop_types)
 real, dimension(0:num_m,num_crop_types) :: central_T, central_P, central_D, variance_T, variance_P, variance_D
 logical :: match_crop_water_source_with_landuse = .FALSE. ! If TRUE, the crop
        !  selection only allows crops with matching water source on the land tiles marked
        !  as irrigated or rain-fed by the land use data set: that is, rain-fed crops are
        !  not allowed on irrigated land use tiles, and irrigated crops are not allowed on
        !  rain-fed tiles.
        !
        ! By default (FALSE), any crop (rain-fed or irrigated) can be selected for any
        ! cropland land use type (irrigated or rain-fed). This is appropriate for the case
        ! when irrigation transitions are not turned on in the land use module: all crop
        ! tiles are tagged as "rain-fed" in this case, so excluding certain crops from
        ! selection would likely result in underestimation of the crop coverage compared
        ! to the real world.

 namelist / vegn_crop_nml / weight_climate, &
            max_planting_SI_SW, Tbase_Wheat, length_of_vernalization_period, &
            max_T_for_vernalization, min_T_GP_SW, absolute_min_T_for_Wheat, &
            match_crop_water_source_with_landuse
 contains
!============================================================================
 subroutine compute_crop_calendars(vegn, diag, L)
 type(vegn_tile_type), intent(inout) :: vegn
 type(diag_buff_type), intent(inout) :: diag
 integer, intent(in) :: L ! index of grid cell which contains this tile
 integer :: m, second, minute, hour, day0, day1, month0, month1, year0, year1, nn
 integer :: iph, iseason, mth, iwater, ipref, pot_crop, iperiod, day_beg, day_opt, day_end
 logical :: new_month
 real, dimension(12) :: rhs
 character(len=512) :: outname, text
 character(len=24) :: crp_name

 integer, dimension(2) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end
 real :: aPTTtH_range(2)
 character(len=2) :: Wheat_type

 if(.not.crop_mod_initialized) call error_mesg('compute_crop_calendars','vegn_crop_init has not been called', FATAL)
 call mpp_clock_begin(cropclock)
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
    crop_loop_1: do ipref=1,num_crop_types
      pot_crop = vegn%Crop%potential_crop(ipref)
      if(is_watch_cell()) then             ! watchpoint_code
        text = ' watchpoint subroutine compute_crop_calendars 0: ipref =   , pot_crop =    = '//trim(crop_name(pot_crop)) ! watchpoint_code
        write(text(57:59),'(i3)') ipref    ! watchpoint_code
        write(text(72:74),'(i3)') pot_crop ! watchpoint_code
        call debug_crop_2(vegn, text)      ! watchpoint_code
      endif                                ! watchpoint_code
      if(pot_crop == NO_CROP) exit crop_loop_1
      if(pot_crop == IRRIGATED_SPRING_WHEAT .or. pot_crop == RAINFED_SPRING_WHEAT) then
        aPTTtH_range = aPTTtH_range_SW
        Wheat_type = 'SW'
      else if(pot_crop == IRRIGATED_WINTER_WHEAT .or. pot_crop == RAINFED_WINTER_WHEAT) then
        aPTTtH_range = aPTTtH_range_WW
        Wheat_type = 'WW'
      endif
      if(pot_crop == IRRIGATED_SPRING_WHEAT .or. pot_crop == IRRIGATED_WINTER_WHEAT .or. pot_crop == RAINFED_SPRING_WHEAT .or. pot_crop == RAINFED_WINTER_WHEAT) then
        call CCA_Wheat(L, vegn, Wheat_type, central_T(:,pot_crop), variance_T(:,pot_crop), central_P(:,pot_crop), & ! intent(in)
                       variance_P(:,pot_crop), central_D(:,pot_crop), variance_D(:,pot_crop), & ! intent(in)
                       SI_crit(pot_crop), max_planting_SI_SW, Tbase_Wheat, aPTTtH_range,      & ! intent(in)
                       length_of_vernalization_period, max_T_for_vernalization, min_T_GP_SW,  & ! intent(in)
                       pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)
        if(pot_crop == IRRIGATED_SPRING_WHEAT .or. pot_crop == IRRIGATED_WINTER_WHEAT) then
          vegn%Crop%crop_calendars(:,1,MAIN_SEASON,ipref) = (/pday(1), hday(1)/)
          vegn%Crop%crop_calendars(:,2,MAIN_SEASON,ipref) = (/pday_beg(1), hday_beg(1)/)
          vegn%Crop%crop_calendars(:,3,MAIN_SEASON,ipref) = (/pday_end(1), hday_end(1)/)
        else if(pot_crop == RAINFED_SPRING_WHEAT .or. pot_crop == RAINFED_WINTER_WHEAT) then
          vegn%Crop%crop_calendars(:,1,MAIN_SEASON,ipref) = (/pday(2), hday(2)/)
          vegn%Crop%crop_calendars(:,2,MAIN_SEASON,ipref) = (/pday_beg(2), hday_beg(2)/)
          vegn%Crop%crop_calendars(:,3,MAIN_SEASON,ipref) = (/pday_end(2), hday_end(2)/)
        endif
        ! Two seasons of Wheat on the same land is rare. Assume it doesn't exist.
        vegn%Crop%crop_calendars(:,:,SECOND_SEASON,ipref) = NO_DATE
      else
        if(pot_crop == IRRIGATED_MAIZE .or. pot_crop == IRRIGATED_SOYBEAN .or. pot_crop == IRRIGATED_RICE) then
          iwater = 1
        else if(pot_crop == RAINFED_MAIZE .or. pot_crop == RAINFED_SOYBEAN .or. pot_crop == RAINFED_RICE) then
          iwater = 2
        endif
        call CCA_Maize_Soybean_Rice(L, vegn, trim(water_source_name(iwater)), GP(pot_crop), central_T(:,pot_crop), variance_T(:,pot_crop), & ! intent(in)
                      central_P(:,pot_crop), variance_P(:,pot_crop), central_D(:,pot_crop), variance_D(:,pot_crop), SI_crit(pot_crop), & ! intent(in)
                      pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)
        vegn%Crop%crop_calendars(:,1,MAIN_SEASON,  ipref) = (/ pday(1), hday(1)/)
        vegn%Crop%crop_calendars(:,1,SECOND_SEASON,ipref) = (/ pday(2), hday(2)/)
        vegn%Crop%crop_calendars(:,2,MAIN_SEASON,  ipref) = (/ pday_beg(1), hday_beg(1)/)
        vegn%Crop%crop_calendars(:,2,SECOND_SEASON,ipref) = (/ pday_beg(2), hday_beg(2)/)
        vegn%Crop%crop_calendars(:,3,MAIN_SEASON,  ipref) = (/ pday_end(1), hday_end(1)/)
        vegn%Crop%crop_calendars(:,3,SECOND_SEASON,ipref) = (/ pday_end(2), hday_end(2)/)
      endif
      if(is_watch_cell()) then                                                                                                 ! watchpoint_code
        if(vegn%Crop%crop_calendars(1,1,MAIN_SEASON,ipref) == NO_DATE) then                                                    ! watchpoint_code
          text = ' watchpoint subroutine compute_crop_calendars 1:'// &                                                        ! watchpoint_code
                 ' The CCA has determined that conditions are unsuitable for cultivation of '//trim(crop_name(pot_crop))       ! watchpoint_code
          call debug_crop_2(vegn, text)                                                                                        ! watchpoint_code
        else                                                                                                                   ! watchpoint_code
          text = ' watchpoint subroutine compute_crop_calendars 2:'// &                                                        ! watchpoint_code
                 ' planting date range for 1st season '//trim(crop_name(pot_crop))//' ='                                       ! watchpoint_code
          nn = len_trim(text)                                                                                                  ! watchpoint_code
          write(text(nn+1:nn+8),'(2i4)') vegn%Crop%crop_calendars(1,2:3,MAIN_SEASON,ipref)                                     ! watchpoint_code
          call debug_crop_2(vegn, text)                                                                                        ! watchpoint_code
          text = ' watchpoint subroutine compute_crop_calendars 3:'// &                                                        ! watchpoint_code
                 ' harvest  date range for 1st season '//trim(crop_name(pot_crop))//' ='                                       ! watchpoint_code
          nn = len_trim(text)                                                                                                  ! watchpoint_code
          write(text(nn+1:nn+8),'(2i4)') vegn%Crop%crop_calendars(2,2:3,MAIN_SEASON,ipref)                                     ! watchpoint_code
          call debug_crop_2(vegn, text)                                                                                        ! watchpoint_code
          if(vegn%Crop%crop_calendars(1,1,SECOND_SEASON,ipref) == NO_DATE) then                                                ! watchpoint_code
            text = ' watchpoint subroutine compute_crop_calendars 4:'// &                                                      ! watchpoint_code
                   ' The CCA has determined that conditions are unsuitable for a second season of '//trim(crop_name(pot_crop)) ! watchpoint_code
            call debug_crop_2(vegn, text)                                                                                      ! watchpoint_code
          else                                                                                                                 ! watchpoint_code
            text = ' watchpoint subroutine compute_crop_calendars 5:'// &                                                      ! watchpoint_code
                   ' planting date range for 2nd season '//trim(crop_name(pot_crop))//' ='                                     ! watchpoint_code
            nn = len_trim(text)                                                                                                ! watchpoint_code
            write(text(nn+1:nn+8),'(2i4)') vegn%Crop%crop_calendars(1,2:3,SECOND_SEASON,ipref)                                 ! watchpoint_code
            call debug_crop_2(vegn, text)                                                                                      ! watchpoint_code
            text = ' watchpoint subroutine compute_crop_calendars 6:'// &                                                      ! watchpoint_code
                   ' harvest  date range for 2nd season '//trim(crop_name(pot_crop))//' ='                                     ! watchpoint_code
            nn = len_trim(text)                                                                                                ! watchpoint_code
            write(text(nn+1:nn+8),'(2i4)') vegn%Crop%crop_calendars(2,2:3,SECOND_SEASON,ipref)                                 ! watchpoint_code
            call debug_crop_2(vegn, text)                                                                                      ! watchpoint_code
          endif                                                                                                                ! watchpoint_code
        endif                                                                                                                  ! watchpoint_code
      endif                                                                                                                    ! watchpoint_code
    enddo crop_loop_1

    call crop_selection(vegn)

    if(is_watch_cell()) then                                                                         ! watchpoint_code
      text = ' watchpoint subroutine compute_crop_calendars 7: chosen_crops = '// &                  ! watchpoint_code
      trim(crop_name(vegn%Crop%chosen_crop(1)))//' '//trim(crop_name(vegn%Crop%chosen_crop(2)))      ! watchpoint_code
      call debug_crop_2(vegn, text)                                                                  ! watchpoint_code
      text = ' watchpoint subroutine compute_crop_calendars 8: planting day = '                      ! watchpoint_code
      nn = len_trim(text)                                                                            ! watchpoint_code
      write(text(nn+1:nn+8),'(2i4)') vegn%Crop%chosen_calendars(1,1),vegn%Crop%chosen_calendars(1,2) ! watchpoint_code
      call debug_crop_2(vegn, text)                                                                  ! watchpoint_code
      text = ' watchpoint subroutine compute_crop_calendars 9: harvest  day = '                      ! watchpoint_code
      nn = len_trim(text)                                                                            ! watchpoint_code
      write(text(nn+1:nn+8),'(2i4)') vegn%Crop%chosen_calendars(2,1),vegn%Crop%chosen_calendars(2,2) ! watchpoint_code
      call debug_crop_2(vegn, text)                                                                  ! watchpoint_code
    endif                                                                                            ! watchpoint_code
 endif ! if(new_month)

 call send_tile_data(id_T_ave, vegn%Crop%tc_av_climate,diag)
 call send_tile_data(id_P_ave, vegn%Crop%precip_av_climate,diag)
 call send_tile_data(id_chosen_calendars(1), real(vegn%Crop%chosen_calendars(:,1)),diag)
 call send_tile_data(id_chosen_calendars(2), real(vegn%Crop%chosen_calendars(:,2)),diag)
 call send_tile_data(id_chosen_crop, real(vegn%Crop%chosen_crop),diag)
 do ipref=1,num_crop_types
   do iseason=1,num_crop_seasons
     do iperiod=1,num_crop_periods
       do iph=1,2
         call send_tile_data(id_crop_calendars(iph,iperiod,iseason,ipref), real(vegn%Crop%crop_calendars(iph,iperiod,iseason,ipref)),diag)
       enddo
     enddo
   enddo
 enddo
 call mpp_clock_end(cropclock)

 end subroutine compute_crop_calendars
!======================================================================================================================================================
 subroutine crop_selection(vegn)
 integer :: calendars(2, num_crop_periods, num_crop_seasons, num_crop_types)
 integer :: potential_crop(num_crop_types)
 integer :: potential_chosen_crops(num_crop_seasons)
 integer :: potential_chosen_calendars(2,num_crop_seasons)
 type(vegn_tile_type), intent(inout) :: vegn
 integer :: iseason

 calendars      = vegn%Crop%crop_calendars
 potential_crop = vegn%Crop%potential_crop
 call crop_selection_sub(calendars, potential_crop, vegn%landuse, potential_chosen_crops, potential_chosen_calendars)

 ! Change the chosen crop and it's calendar only if it is not active
 do iseason=1,num_crop_seasons
   if(.not.vegn%Crop%chosen_crop_is_active(iseason)) then
     vegn%Crop%chosen_calendars(:,iseason) = potential_chosen_calendars(:,iseason)
     vegn%Crop%chosen_crop(iseason)        = potential_chosen_crops(iseason)
   endif
 enddo

 end subroutine crop_selection
!======================================================================================================================================================
 subroutine crop_selection_sub(calendars, potential_crop, landuse, chosen_crops, chosen_calendars)
 integer, intent(in)  :: calendars(2, num_crop_periods, num_crop_seasons, num_crop_types)
 integer, intent(in)  :: potential_crop(num_crop_types)
 integer, intent(in)  :: landuse
 integer, intent(out) :: chosen_crops(num_crop_seasons)
 integer, intent(out) :: chosen_calendars(2,num_crop_seasons)
 integer :: ipref, ipref_beg, ipref_end
 integer, dimension(2,num_crop_seasons) :: dble_cropping_calendar
 character(len=256) :: text
 character(len=16) :: cn1, cn2

 integer :: vcal_1, vcal_2, vcal_3, vcal_4, vcal_5, vcal_6, vcal_7, vcal_8, vcal_9, vcal_10
! vcal_1 = The highest crop preference among crops with valid crop calendars. Valid calendar number 1.
! vcal_2 = The next highest. Valid calendar number 2.
! vcal_3 = The next. And so on for vcal_4 through vcal_10
! Example
! potential_crop(1) = RAINFED_RICE           calendars(1,1,1,1) = NO_DATE
! potential_crop(2) = IRRIGATED_SPRING_WHEAT calendars(1,1,1,2) = (a valid day of year)
! potential_crop(3) = RAINFED_MAIZE          calendars(1,1,1,3) = NO_DATE
! potential_crop(4) = RAINFED_SOYBEAN        calendars(1,1,1,4) = (a valid day of year)
! potential_crop(icrop,icrop=5,num_crop_types) = NO_CROP

! In the example above vcal_1 = 2, vcal_2 = 4, vcal_3 to vcal_10 = 0
! potential_crop(vcal_1) = IRRIGATED_SPRING_WHEAT, potential_crop(vcal_2) = RAINFED_SOYBEAN

! Find the two crops of highest preference which have valid optimal planting dates for the 1st season.
! Valid optimal planting dates for the 1st season exist when the CCA has determined that conditions are suitable.
 vcal_1 = 0; vcal_2 = 0; vcal_3 = 0; vcal_4 = 0; vcal_5 = 0
 vcal_6 = 0; vcal_7 = 0; vcal_8 = 0; vcal_9 = 0; vcal_10 = 0
 if(match_crop_water_source_with_landuse) then
   ! In this case, irrigated and rain-fed crops are distinct land use types:
   ! the result of land use transitions that distinguish them.
   select case (landuse)
   case (LU_IRRIG) ! irrigated crop
     ipref_beg = 1
     ipref_end = 5
   case (LU_RAINF) ! rain-fed crop
     ipref_beg = 6
     ipref_end = 10
   case default
     call error_mesg('crop_selection','landuse argument must be LU_IRRIG or LU_RAINF', FATAL)
   end select
 else
   ! Irrigated or rain-fed crops could be selected for any crop tile, regardless of water
   ! source. NOTE: if do_irrigation is FALSE in transition, all crop tiles have LU_RAINF
   ! land use type.
   ipref_beg = 1
   ipref_end = 10
 endif
 ipref_loop_1: do ipref=ipref_beg,ipref_end
   if(calendars(1,1,1,ipref) == NO_DATE) cycle ipref_loop_1
   if(vcal_1 == 0) then
     vcal_1 = ipref
   else if(vcal_2 == 0) then
     vcal_2 = ipref
   else if(vcal_3 == 0) then
     vcal_3 = ipref
   else if(vcal_4 == 0) then
     vcal_4 = ipref
   else if(vcal_5 == 0) then
     vcal_5 = ipref
   else if(vcal_6 == 0) then
     vcal_6 = ipref
   else if(vcal_7 == 0) then
     vcal_7 = ipref
   else if(vcal_8 == 0) then
     vcal_8 = ipref
   else if(vcal_9 == 0) then
     vcal_9 = ipref
   else if(vcal_10 == 0) then
     vcal_10 = ipref
   endif
 enddo ipref_loop_1
 if(vcal_1 == 0) then
   ! There are no crops for which conditions are suitable.
   chosen_crops(:) = NO_CROP
   chosen_calendars(:,1) = (/NO_DATE,NO_DATE/)
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

 if(vcal_2 == 0) then
   ! There is only one crop for which conditions are suitable.
   if(calendars(1,1,2,vcal_1) == NO_DATE) then
     ! In addition, this crop has only one season.
     ! There is no need to test for the possibility of double-cropping.
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
     chosen_crops(2) = NO_CROP
     chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
     return
   endif
 endif

! Execution reaches this point only if there is the possibility of double-cropping.
! First test if double-cropping of the same crop is possible.
 if(calendars(1,1,2,vcal_1) /= NO_CROP) then
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_1), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_1)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if double-cropping of the same crop is not possible.

 if(vcal_2 == 0) then
   ! There is no crop vcal_2. Double-cropping is not possible.
   chosen_crops(1) = potential_crop(vcal_1)
   chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
   chosen_crops(2) = NO_CROP
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 1st season of crop number potential_crop(vcal_2).
 if(calendars(1,1,1,vcal_2) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_2))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,1,vcal_2), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_2)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 2nd season of crop number potential_crop(vcal_2).
 if(calendars(1,1,2,vcal_2) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_2))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_2), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_2)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if double-cropping of crop vcal_1 with vcal_2 is not possible. Try double-cropping with crop vcal_3.

 if(vcal_3 == 0) then
   ! There is no crop vcal_3. Double-cropping is not possible.
   chosen_crops(1) = potential_crop(vcal_1)
   chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
   chosen_crops(2) = NO_CROP
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 1st season of crop number potential_crop(vcal_3).
 if(calendars(1,1,1,vcal_3) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_3))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,1,vcal_3), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_3)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 2nd season of crop number potential_crop(vcal_3).
 if(calendars(1,1,2,vcal_3) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_3))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_3), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_3)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if double-cropping of crop vcal_1 with vcal_3 is not possible. Try double-cropping with crop vcal_4.

 if(vcal_4 == 0) then
   ! There is no crop vcal_4. Double-cropping is not possible.
   chosen_crops(1) = potential_crop(vcal_1)
   chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
   chosen_crops(2) = NO_CROP
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 1st season of crop number potential_crop(vcal_4).
 if(calendars(1,1,1,vcal_4) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_4))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,1,vcal_4), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_4)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 2nd season of crop number potential_crop(vcal_4).
 if(calendars(1,1,2,vcal_4) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_4))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_4), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_4)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if double-cropping of crop vcal_1 with vcal_4 is not possible. Try double-cropping with crop vcal_5.

 if(vcal_5 == 0) then
   ! There is no crop vcal_5. Double-cropping is not possible.
   chosen_crops(1) = potential_crop(vcal_1)
   chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
   chosen_crops(2) = NO_CROP
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 1st season of crop number potential_crop(vcal_5).
 if(calendars(1,1,1,vcal_5) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_5))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,1,vcal_5), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_5)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 2nd season of crop number potential_crop(vcal_5).
 if(calendars(1,1,2,vcal_5) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_5))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_5), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_5)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if double-cropping of crop vcal_1 with vcal_5 is not possible. Try double-cropping with crop vcal_6.

 if(vcal_6 == 0) then
   ! There is no crop vcal_6. Double-cropping is not possible.
   chosen_crops(1) = potential_crop(vcal_1)
   chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
   chosen_crops(2) = NO_CROP
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 1st season of crop number potential_crop(vcal_6).
 if(calendars(1,1,1,vcal_6) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_6))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,1,vcal_6), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_6)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 2nd season of crop number potential_crop(vcal_6).
 if(calendars(1,1,2,vcal_6) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_6))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_6), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_6)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if double-cropping of crop vcal_1 with vcal_6 is not possible. Try double-cropping with crop vcal_7.

 if(vcal_7 == 0) then
   ! There is no crop vcal_7. Double-cropping is not possible.
   chosen_crops(1) = potential_crop(vcal_1)
   chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
   chosen_crops(2) = NO_CROP
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 1st season of crop number potential_crop(vcal_7).
 if(calendars(1,1,1,vcal_7) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_7))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,1,vcal_7), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_7)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 2nd season of crop number potential_crop(vcal_7).
 if(calendars(1,1,2,vcal_7) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_7))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_7), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_7)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if double-cropping of crop vcal_1 with vcal_7 is not possible. Try double-cropping with crop vcal_8.

 if(vcal_8 == 0) then
   ! There is no crop vcal_8. Double-cropping is not possible.
   chosen_crops(1) = potential_crop(vcal_1)
   chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
   chosen_crops(2) = NO_CROP
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 1st season of crop number potential_crop(vcal_8).
 if(calendars(1,1,1,vcal_8) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_8))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,1,vcal_8), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_8)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 2nd season of crop number potential_crop(vcal_8).
 if(calendars(1,1,2,vcal_8) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_8))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_8), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_8)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if double-cropping of crop vcal_1 with vcal_8 is not possible. Try double-cropping with crop vcal_9.

 if(vcal_9 == 0) then
   ! There is no crop vcal_9. Double-cropping is not possible.
   chosen_crops(1) = potential_crop(vcal_1)
   chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
   chosen_crops(2) = NO_CROP
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 1st season of crop number potential_crop(vcal_9).
 if(calendars(1,1,1,vcal_9) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_9))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,1,vcal_9), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_9)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 2nd season of crop number potential_crop(vcal_9).
 if(calendars(1,1,2,vcal_9) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_9))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_9), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_9)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if double-cropping of crop vcal_1 with vcal_9 is not possible. Try double-cropping with crop vcal_10.

 if(vcal_10 == 0) then
   ! There is no crop vcal_10. Double-cropping is not possible.
   chosen_crops(1) = potential_crop(vcal_1)
   chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
   chosen_crops(2) = NO_CROP
   chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
   return
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 1st season of crop number potential_crop(vcal_10).
 if(calendars(1,1,1,vcal_10) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_10))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,1,vcal_10), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_10)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif

! Test if double-cropping is possible by combining the 1st season of crop number
! potential_crop(vcal_1) with the 2nd season of crop number potential_crop(vcal_10).
 if(calendars(1,1,2,vcal_10) /= NO_CROP) then
   cn1 = crop_name(potential_crop(vcal_1))
   cn2 = crop_name(potential_crop(vcal_10))
   call search_for_non_overlapping_seasons(calendars(:,:,1,vcal_1), calendars(:,:,2,vcal_10), dble_cropping_calendar)
   if(dble_cropping_calendar(1,2) /= NO_DATE) then
     chosen_crops(1) = potential_crop(vcal_1)
     chosen_crops(2) = potential_crop(vcal_10)
     chosen_calendars = dble_cropping_calendar
     return
   endif
 endif
!--------------------------------------------------------------------------------------------------------------------------------------
! Execution reaches this point only if all attempts to double-crop have failed.
! The 1st season of crop number potential_crop(vcal_1) is the only possible cropping season.
 chosen_crops(1) = potential_crop(vcal_1)
 chosen_calendars(:,1) = (/calendars(1,1,1,vcal_1),calendars(2,1,1,vcal_1)/)
 chosen_crops(2) = NO_CROP
 chosen_calendars(:,2) = (/NO_DATE,NO_DATE/)
 return
 end subroutine crop_selection_sub
!==================================================================================================================================================
 subroutine search_for_non_overlapping_seasons(p_and_h_date_ranges_1, p_and_h_date_ranges_2, dble_cropping_calendar)
 integer, dimension(2,num_crop_periods), intent(in) :: p_and_h_date_ranges_1, p_and_h_date_ranges_2
 integer, dimension(2,num_crop_seasons), intent(out) :: dble_cropping_calendar
 logical :: equal_GP_1, equal_GP_2
 integer :: GP_1, GP_2
 integer :: dopt1, dbeg1, dend1, pday1, dopt2, dbeg2, dend2, pday2, hday1, hday2, iperiod
 integer :: K1, K2, DE1, DE2, DO1, DO2
 logical :: B1, B2, E1, E2, I1, I2
 character(len=256) :: text

 ! Test for a pair of planting dates, one from p_and_h_date_ranges_1 and one from p_and_h_date_ranges_2, that do not result in overlapping growing periods.
 ! There may be no such pairs, in which case double-cropping is not possible.
 ! There may be multiple pairs, in which case the pair chosen is the first one found that results in no overlap.

 ! Function flip_flop needs to remember the state it was in after the previous execution.
 ! Several variables exist for this purpose. They are:
 ! K1, DE1, DO1, B1, E1, I1 and K2, DE2, DO2, B2, E2, I2
 ! They are re-assigned every execution of flip_flop and passed back in on the subsequent execution.

 dble_cropping_calendar = NO_DATE ! Will be over-written if a non-overlapping pair of dates is found

 equal_GP_1 = equal_GP(p_and_h_date_ranges_1)
 if(equal_GP_1) then
   GP_1 = p_and_h_date_ranges_1(2,1) - p_and_h_date_ranges_1(1,1)
   if(GP_1 < 0) GP_1 = GP_1 + 365
 endif

 equal_GP_2 = equal_GP(p_and_h_date_ranges_2)
 if(equal_GP_2) then
   GP_2 = p_and_h_date_ranges_2(2,1) - p_and_h_date_ranges_2(1,1)
   if(GP_2 < 0) GP_2 = GP_2 + 365
 endif

 if(.not.equal_GP_1 .and. .not.equal_GP_2) then
   dble_cropping_calendar(1,1) = p_and_h_date_ranges_1(1,1)
   dble_cropping_calendar(2,1) = p_and_h_date_ranges_1(2,1)
   dble_cropping_calendar(1,2) = NO_DATE
   dble_cropping_calendar(2,2) = NO_DATE
   return
 endif

 if(equal_GP_1 .and. .not.equal_GP_2) then
   dopt1 = p_and_h_date_ranges_1(1,1)
   dbeg1 = p_and_h_date_ranges_1(1,2)
   dend1 = p_and_h_date_ranges_1(1,3)
   period_loop_1: do iperiod=1,num_crop_periods
     pday2 = p_and_h_date_ranges_2(1,iperiod)
     hday2 = p_and_h_date_ranges_2(2,iperiod)
     I1 = .false. ! Tells function flip_flop that it needs to initialize. Flipped to .true. after the first execution of flip_flop.
     loop1: do
       pday1 = flip_flop(NO_DATE, dopt1, dbeg1, dend1, K1, DE1, DO1, B1, E1, I1)
       if(B1 .and. E1) exit loop1
       if(pday1 == 0) cycle loop1
       hday1 = pday1 + GP_1
       hday1 = modulo_no_zero(hday1,365)
       if(no_overlap(pday1, hday1, pday2, hday2)) then
         !----------------------------------
          dble_cropping_calendar(1,1) = pday1
          hday1 = pday1 + GP_1
          hday1 = modulo_no_zero(hday1,365)
          dble_cropping_calendar(2,1) = hday1
         !----------------------------------
          dble_cropping_calendar(1,2) = pday2
          dble_cropping_calendar(2,2) = hday2
         !----------------------------------
         return
       else
       endif
     enddo loop1
   enddo period_loop_1
 endif

 if(equal_GP_2 .and. .not.equal_GP_1) then
   dopt2 = p_and_h_date_ranges_2(1,1)
   dbeg2 = p_and_h_date_ranges_2(1,2)
   dend2 = p_and_h_date_ranges_2(1,3)
   period_loop_2: do iperiod=1,num_crop_periods
     pday1 = p_and_h_date_ranges_1(1,iperiod)
     hday1 = p_and_h_date_ranges_1(2,iperiod)
     I2 = .false. ! Tells function flip_flop that it needs to initialize. Flipped to .true. after the first execution of flip_flop.
     loop2: do
       pday2 = flip_flop(NO_DATE, dopt2, dbeg2, dend2, K2, DE2, DO2, B2, E2, I2)
       if(B2 .and. E2) exit loop2
       if(pday2 == 0) cycle loop2
       hday2 = pday2 + GP_2
       hday2 = modulo_no_zero(hday2,365)
       if(no_overlap(pday1, hday1, pday2, hday2)) then
         !----------------------------------
          dble_cropping_calendar(1,2) = pday2
          hday2 = pday2 + GP_2
          hday2 = modulo_no_zero(hday2,365)
          dble_cropping_calendar(2,2) = hday2
         !----------------------------------
          dble_cropping_calendar(1,1) = pday1
          dble_cropping_calendar(2,1) = hday1
         !----------------------------------
         return
       endif
     enddo loop2
   enddo period_loop_2
 endif

 if(equal_GP_1 .and. equal_GP_2) then
   dopt1 = p_and_h_date_ranges_1(1,1)
   dbeg1 = p_and_h_date_ranges_1(1,2)
   dend1 = p_and_h_date_ranges_1(1,3)
   I1 = .false. ! Tells function flip_flop that it needs to initialize each iteration of loop3. Flipped to .true. after the first execution of flip_flop.
   loop3: do
     pday1 = flip_flop(NO_DATE, dopt1, dbeg1, dend1, K1, DE1, DO1, B1, E1, I1)
     hday1 = pday1 + GP_1
     hday1 = modulo_no_zero(hday1,365)
     if(B1 .and. E1) exit loop3
     if(pday1 == 0) cycle loop3
     I2 = .false. ! Tells function flip_flop within loop4 that it needs to initialize. It must initialize each interation of loop3.
     loop4: do
       dopt2 = p_and_h_date_ranges_2(1,1)
       dbeg2 = p_and_h_date_ranges_2(1,2)
       dend2 = p_and_h_date_ranges_2(1,3)
       pday2 = flip_flop(NO_DATE, dopt2, dbeg2, dend2, K2, DE2, DO2, B2, E2, I2)
       if(B2 .and. E2) exit loop4
       if(pday2 == 0) cycle loop4
       hday2 = pday2 + GP_2
       hday2 = modulo_no_zero(hday2,365)
       if(no_overlap(pday1, hday1, pday2, hday2)) then
        !----------------------------------
         dble_cropping_calendar(1,1) = pday1
         dble_cropping_calendar(2,1) = hday1
        !----------------------------------
         dble_cropping_calendar(1,2) = pday2
         hday2 = pday2 + GP_2
         hday2 = modulo_no_zero(hday2,365)
         dble_cropping_calendar(2,2) = hday2
        !----------------------------------
         return
       endif
     enddo loop4
   enddo loop3
 endif

! Execution reaches this point only if all days tested result in overlap of growing seasons.
! Assign p_and_h_date_ranges_1 to the 1st season and NO_DATE to the 2nd.
 dble_cropping_calendar(:,1) = p_and_h_date_ranges_1(:,1)
 dble_cropping_calendar(:,2) = NO_DATE
 end subroutine search_for_non_overlapping_seasons
!==================================================================================================================================================
 logical function equal_GP(p_and_h_date_ranges)
 integer, dimension(2,num_crop_periods), intent(in) :: p_and_h_date_ranges
 integer, dimension(num_crop_periods) :: GP_tmp
 integer :: iperiod

 equal_GP = .true.
 GP_tmp(1) = p_and_h_date_ranges(2,1) - p_and_h_date_ranges(1,1)
 if(GP_tmp(1) < 0) GP_tmp(1) = GP_tmp(1) + 365
 do_loop: do iperiod=2,num_crop_periods
   GP_tmp(iperiod) = p_and_h_date_ranges(2,iperiod) - p_and_h_date_ranges(1,iperiod)
   if(GP_tmp(iperiod) < 0) GP_tmp(iperiod) = GP_tmp(iperiod) + 365
   if(GP_tmp(iperiod) /= GP_tmp(1)) then
     equal_GP = .false.
     exit do_loop
   endif
 enddo do_loop
 end function equal_GP
!==================================================================================================================================================
 integer function flip_flop(NO_DATE, dopt, dbeg, dend, K, dend_local, dopt_local, before_beg, after_end, initialized) result(day)
 integer, intent(in) :: NO_DATE, dopt, dbeg, dend
 integer, intent(inout) :: K, dend_local, dopt_local
 logical, intent(inout) :: before_beg, after_end, initialized

 if(.not.initialized) then
   if(dend < dbeg) then
     dend_local = dend + 365
   else
     dend_local = dend
   endif
   if(dopt < dbeg) then
     dopt_local = dopt + 365
   else
     dopt_local = dopt
   endif
   K = 0
   before_beg = .false.
   after_end  = .false.
   initialized = .true.
 else
   if(K == 0) then
     K = -1
   else if(K > 0) then
     K = -(K+1)
   else if(K < 0) then
     K = abs(K)
   endif
 endif

 if(before_beg .and. after_end) then
   day = NO_DATE
   K = 0
   dend_local = 0
   return
 endif

 day = dopt_local + 5*K
 if(day < dbeg) then
   day = NO_DATE
   before_beg = .true.
   return
 endif

 if(day > dend_local) then
   day = NO_DATE
   after_end  = .true.
   return
 endif
 day = modulo_no_zero(day,365)
 end function flip_flop
!======================================================================================================================================================
 logical function no_overlap(pday1, hday1, pday2, hday2)
 integer, intent(in) :: pday1, hday1, pday2, hday2
 integer :: tbeg1_local, tend1_local, pday2_local, tend2_local, beg_doy

 if(any((/pday1,hday1,pday2,hday2/) == (/NO_DATE,NO_DATE,NO_DATE,NO_DATE/))) then
   no_overlap = .false.
   return
 endif
 if(hday1 < pday1) then
   tend1_local = hday1 + 365
 else
   tend1_local = hday1
 endif
 if(hday2 < pday2) then
   tend2_local = hday2 + 365
 else
   tend2_local = hday2
 endif
 beg_doy = max(tend1_local,tend2_local) - 365
 if(pday1 <= beg_doy) then
   tbeg1_local = pday1 + 365
   tend1_local = tend1_local + 365
 else
   tbeg1_local = pday1
   tend1_local = tend1_local
 endif
 if(pday2 <= beg_doy) then
   pday2_local = pday2 + 365
   tend2_local = tend2_local + 365
 else
   pday2_local = pday2
   tend2_local = tend2_local
 endif
 no_overlap = .not.(tend1_local >= pday2_local .and. tend2_local >= tbeg1_local)
 end function no_overlap
!======================================================================================================================================================
 subroutine CCA_Maize_Soybean_Rice(L, vegn, water, GP_in, central_T, variance_T, central_P, variance_P, &
                     central_D, variance_D, SI_crit, pday, pday_beg, pday_end, hday, hday_beg, hday_end)
 integer, intent(in) :: L
 type(vegn_tile_type), intent(in) :: vegn
 character(len=*), intent(in) :: water
 integer, intent(in) :: GP_in
 real, intent(in) :: central_T(0:), variance_T(0:), central_P(0:), variance_P(0:), central_D(0:), variance_D(0:), SI_crit
 integer, dimension(num_crop_seasons), intent(out) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end
 integer :: k, km, kp, daybeg, mths_after, doy, k_at_SI_min, ktest, day, iseason, max_range_length
 real :: Temp, Prec, annual_SI_min, dlen
 real :: TSI, DSI, PSI
 real, dimension(num_test_days) :: SI
 character(len=24) :: crp_name
!---------------------------------------------------------------------
 if(trim(water) /= 'irrigated' .and. trim(water) /= 'rainfed') then
   call error_mesg('CCA_Maize_Soybean_Rice ERROR: '//trim(water), 'is not a valid value of water', FATAL)
 endif
 k_loop_1: do k=1,num_test_days ! Compute the suitability index at 5 day intervals, starting with Jan 5
   TSI = 0.0
   PSI = 0.0
   DSI = 0.0
   daybeg = 5*k
   mths_loop: do mths_after=0,num_m
     day = daybeg + 30*mths_after
     doy = modulo_no_zero(day,365)
     km = doy/5
     kp = km+1
     Temp = interp_between_mid_mths(doy, vegn%Crop%T_mid_mth)
     TSI = TSI + (Temp - central_T(mths_after))**2/variance_T(mths_after)
     if(mths_after < 4) then
       ! Month 4 is not tested for precip or day length
       Prec = interp_between_mid_mths(doy, vegn%Crop%P_mid_mth)
       if(trim(water) == 'irrigated') Prec = max(Prec,central_P(mths_after))
       PSI = PSI + (Prec - central_P(mths_after))**2/variance_P(mths_after)
       dlen = .2*((doy-5*km)*day_length(kp,L) + (5*kp-doy)*day_length(km,L))
       DSI = DSI + (dlen - central_D(mths_after))**2/variance_D(mths_after)
     endif
   enddo mths_loop
   SI(k) = TSI + DSI + PSI
   if(SI(k) > SI_crit) cycle k_loop_1
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
   max_range_length = 179 - GP_in
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
     hday_beg(iseason) = modulo_no_zero(pday_beg(iseason) + GP_in, 365)
     hday(iseason)     = modulo_no_zero(pday(iseason)     + GP_in, 365)
     hday_end(iseason) = modulo_no_zero(pday_end(iseason) + GP_in, 365)
   else
     hday_beg(iseason) = NO_DATE
     hday(iseason)     = NO_DATE
     hday_end(iseason) = NO_DATE
   endif
 enddo

 end subroutine CCA_Maize_Soybean_Rice
!======================================================================================================================================================
 subroutine CCA_Wheat(L, vegn, Wtype, central_T, variance_T, central_P, variance_P, central_D, variance_D, & ! intent(in)
                      SI_crit, max_planting_SI_SW, Tbase, aPTTtH_range, length_of_vernalization_period, & ! intent(in)
                      max_T_for_vernalization, min_T_GP_SW, & ! intent(in)
                      pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)

! 1. Compute dates of accumulated photo-thermal time at intervals of 200 units, from zero to 800, for each
!    candidate Optimal Planting Date (OPD) starting with Jan 5 and at five day intervals throughout the year.

! 2. If the accumulated photo-thermal time does not exceed 800 units starting from any date then the climate is deemed unsuitable for wheat.

! 3. Compute suitability index using climatic conditions at intervals of 200 units of accumulated photo-thermal time
!    The suitability index is specific to the water source and variety:
!    irrigated winter wheat, rainfed winter wheat, irrigated spring wheat, rainfed spring wheat

! 4. If the suitability index for the specific type of wheat being tested exceeds the critical
!    value at all tested dates throughout the year then the climate is deemed unsuitable.

! 5. Reduce the candidate OPDs to those for which the suitability index is below the critical value.
!    The corresponding harvest dates are the dates when the accumulated photo-thermal time reaches 837 units or the maximum, starting from the candidate OPD.

! 6. Reduce the candidate OPDs to those for which the temperature never drops below absolute_min_T_for_Wheat before the corresponding harvest date.

! 7. Reduce the candidate OPDs to those which are warmer than min_T_GP_SW.

! 8. For winter wheat: Reduce the candidate OPDs to those for which temperature drops below
!    max_T_for_vernalization for at least 40 days between the planting and harvest dates.
!    For spring wheat: Reduce the candidate OPDs to those for which temperature remains above min_T_GP_SW between the planting and harvest dates.

! 9. For winter wheat: The predicted OPD is the date of minimum suitability index among the remaining candidate OPDs.
!    The predicted harvest date is the date when the accumulated photo-thermal time reaches 837 units or the maximum, starting from the candidate OPD.
!    For spring wheat: The predicted OPD is the date of minimum suitability index among the remaining candidate OPDs if the minimum is between 3.25 and 9.0
!    or, if the minimum is below 3.25, the date prior to the the date of the minimum when it reaches 3.25

!10. The range of suitable dates includes all contiguous dates having a suitability index below critial before and after the OPD.
!======================================================================================================================================================
 integer, intent(in) :: L ! index of grid cell which contains this tile
 type(vegn_tile_type), intent(inout) :: vegn
 character(len=2), intent(in) :: Wtype ! 'SW' or 'WW'
 real, intent(in) :: central_T(0:), variance_T(0:), central_P(0:), variance_P(0:), central_D(0:), variance_D(0:) ! At intervals of 200 aPTT units after planting. (0) is planting day.
 real, intent(in) :: SI_crit, max_planting_SI_SW, Tbase, aPTTtH_range(2)
 integer, intent(in) :: length_of_vernalization_period
 real, intent(in) :: max_T_for_vernalization, min_T_GP_SW
 integer, dimension(num_water), intent(out) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end

 integer :: k, k2, k2m, k2p, km, kp, daybeg, crossing_point, kautumn, k_of_ann_SI_min, kk, kkp, k_of_ann_SI_max, doy, iwater
 real :: Temp, Prec, TSI_test, DSI_test, dlen
 real :: PSI_test(2), SI_test(2) ! first element for irrigated, second for rainfed
 real :: annual_SI_max, annual_SI_min, min_T_planting
 real :: SI(num_test_days,2) ! Suitability Index. Computed at 5 day intervals from Jan 5 to Dec 31.
 integer :: crossing_day_400(num_test_days) ! Date at which accumulated photo-thermal time (aPTT) since planting reaches 400 units
 integer :: crossing_days(0:num_m)
 integer :: hday_list(num_test_days) ! Remember the values for each test date then choose the one that corresponds to the annual minimum suitability index.
 real :: aPTTtH_list(num_test_days)  ! Remember the values for each test date then choose the one that corresponds to the annual minimum suitability index.
 logical :: passes_other_criteria
 character(len=256) :: text ! watchpoint_code

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
     cycle k_loop_1 ! cycle k loop if aPTT never reaches aPTTtH_range(1) or if the temperature drops
                    ! below absolute_min_T_for_Wheat before aPTTtH_range(1) units of aPTT is reached.
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
       SI(k,:) = unsuitable ! If irrigated Wheat exceeds critical, then so does rainfed.
       if(is_watch_cell()) then                                                          ! watchpoint_code
         text = ' Conditions are unsuitable for either irrigated or rainfed '//Wtype// & ! watchpoint_code
                ' if planted on day     because the suitability index exceeds critical'  ! watchpoint_code
         write(text(81:83),'(i3)') 5*k                                                   ! watchpoint_code
         call debug_crop_2(vegn, text)                                                   ! watchpoint_code
       endif                                                                             ! watchpoint_code
       cycle k_loop_1
     endif
   enddo ! do crossing_point=0,num_m
   SI(k,1) = SI_test(1) ! Conditions are suitable for planting irrigated Wheat on day of the year 5*k, provided it passes the tests in k_loop_2 and k_loop_3
   if(SI_test(2) > SI_crit) then
     SI(k,2) = unsuitable
     if(is_watch_cell()) then                                                         ! watchpoint_code
       text = ' Conditions are unsuitable for rainfed '//Wtype// &                    ! watchpoint_code
              ' if planted on day     because the suitability index exceeds critical' ! watchpoint_code
       write(text(61:63),'(i3)') 5*k                                                  ! watchpoint_code
       call debug_crop_2(vegn, text)                                                  ! watchpoint_code
     endif                                                                            ! watchpoint_code
   else
     SI(k,2) = SI_test(2) ! Conditions are suitable for planting rainfed Wheat on day of the year 5*k, provided it passes the tests in k_loop_2 and k_loop_3
   endif
 enddo k_loop_1

 water_loop: do iwater=1,num_water
   passes_other_criteria = .true.
   k_loop_2: do k=1,num_test_days ! Step 8: check that vernalization is possible for winter wheat and that temperature remains above min_T_GP_SW for spring wheat.
     if(SI(k,iwater)==unsuitable) cycle k_loop_2 ! Step 4 If the suitability index for the specific type of wheat being tested exceeds the critical value at
                                                 ! all tested dates throughout the year then the tests within k_loop_2 are not necessary and will be skipped.
     if(Wtype == 'SW') then
       ! If the temperature drops below min_T_GP_SW during the growing period then flag it as unsuitable for planting.
       if(T_too_cold_during_GP(5*k, hday_list(k), vegn%Crop%T_mid_mth, min_T_GP_SW)) then
         SI(k,iwater) = unsuitable ! Step 8
         passes_other_criteria = .false.
         text = ' Conditions are unsuitable for either irrigated or rainfed '//Wtype// & ! watchpoint_code
              ' if planted on day     because the climatological mean temperature'// &   ! watchpoint_code
              ' drops below min_T_GP_SW during what would othwise be a suitable growing period'  ! watchpoint_code
         write(text(81:83),'(i3)') 5*k                                                   ! watchpoint_code
         call debug_crop_2(vegn, text)                                                   ! watchpoint_code
       endif
     endif
     if(Wtype == 'WW') then
       if(.not.vernalization_is_possible(5*k, crossing_day_400(k), vegn%Crop%T_mid_mth, length_of_vernalization_period, max_T_for_vernalization)) then
         SI(k,iwater) = unsuitable ! Step 8
         passes_other_criteria = .false.
         text = ' Conditions are unsuitable for either irrigated or rainfed '//Wtype// &                                           ! watchpoint_code
                ' if planted on day     because the evolution of climatological mean temperature does not allow for vernalization' ! watchpoint_code
         write(text(81:83),'(i3)') 5*k                                                                                             ! watchpoint_code
         call debug_crop_2(vegn, text)                                                                                             ! watchpoint_code
       endif
     endif
   enddo k_loop_2

   min_T_planting = min_T_GP_SW
   k_loop_3: do k=1,num_test_days ! Do not plant when the temperature is below min_T_planting
     if(SI(k,iwater) == unsuitable) cycle k_loop_3
     Temp = interp_between_mid_mths(5*k, vegn%Crop%T_mid_mth)
     if(Temp < min_T_planting) then
       SI(k,iwater) = unsuitable ! Step 7
       passes_other_criteria = .false.
       text = ' Conditions are unsuitable for either irrigated or rainfed '//Wtype// &                                                  ! watchpoint_code
              ' because the climatological mean temperature is below min_T_GP_SW on day    , which would otherwise be a suitable planting date' ! watchpoint_code
       write(text(127:129),'(i3)') 5*k                                                                                                  ! watchpoint_code
       call debug_crop_2(vegn, text)                                                                                                    ! watchpoint_code
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
 logical function vernalization_is_possible(pday, day_aPPT_400, T_mid_mth, length_of_vernalization_period, max_T_for_vernalization)
 integer, intent(in) :: pday, day_aPPT_400
 real, intent(in) :: T_mid_mth(12)
 integer, intent(in) :: length_of_vernalization_period
 real, intent(in) :: max_T_for_vernalization
 integer :: ndays, day, num_vern_days, doy
 real :: Temp

 vernalization_is_possible = .FALSE.
 ndays = day_aPPT_400 - pday
 if(ndays <= 0) ndays = ndays + 365
 num_vern_days = 0
 day_loop: do day=pday,pday+ndays
   doy = modulo_no_zero(day,365)
   Temp = interp_between_mid_mths(doy, T_mid_mth)
   if(Temp > TFREEZE .and. Temp < max_T_for_vernalization) then
     num_vern_days = num_vern_days + 1
     if(num_vern_days > length_of_vernalization_period) then
       vernalization_is_possible = .true.
       exit day_loop
     endif
   endif
 enddo day_loop
 end function vernalization_is_possible
!======================================================================================================================================================
 logical function T_too_cold_during_GP(pday, hday, T_mid_mth, min_T_GP_SW)
 integer, intent(in) :: pday, hday
 real, intent(in) :: T_mid_mth(12), min_T_GP_SW
 integer :: GP, day, doy
 real :: Temp

 T_too_cold_during_GP = .FALSE.
 GP = hday - pday
 if(GP <= 0) GP = GP + 365
 day_loop: do day=pday,pday+GP
   doy = modulo_no_zero(day,365)
   Temp = interp_between_mid_mths(doy, T_mid_mth)
   if(Temp < min_T_GP_SW) then
     T_too_cold_during_GP = .TRUE.
     exit day_loop
   endif
 enddo day_loop
 end function T_too_cold_during_GP
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
 ! crossing_days(m) = NO_DATE if m*aPTT_interval is never reached or if temperature drops
 ! below absolute_min_T_for_Wheat before aPTT reaches a value of num_m*aPTT_interval.
 ! The date being tested is not suitable for planting either Spring or Winter Wheat if any of crossing_days(:) returned is NO_DATE

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
  day_loop: do dd=daybeg,daybeg+363
   doy_today    = modulo_no_zero(dd,365)
   doy_tomorrow = modulo_no_zero(dd+1,365)
   if(aPTT(doy_today) > aPTTtH_range(2)) exit day_loop
   T_today = interp_between_mid_mths(doy_today, T_mid_mth)
   if(T_today < absolute_min_T_for_Wheat) exit day_loop
   if(T_today < Tbase .and. aPTT(doy_today) > aPTTtH_range(1)) exit day_loop
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
 integer :: logunit, io, ierr

  read(input_nml_file, nml=vegn_crop_nml, iostat=io)
  ierr = check_nml_error(io, 'vegn_crop_nml')
  logunit = stdlog()
  write(logunit, nml=vegn_crop_nml)
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
 integer :: year0, month0, day0, year1, month1, day1, hour, minute, second
 type(land_tile_enum_type) :: ce
 type(land_tile_type), pointer :: tile
 integer :: day_ae, month_ae, year_ae, hour_ae, minute_ae, second_ae, ierr, icrop, iwater, iseason
 real, allocatable, dimension(:,:) :: MIRCA_crop_frac, crop_frac_tmp, Tclim, Pclim
 integer, allocatable, dimension(:,:) :: potential_crop
 character(len=256) :: infile
 integer, dimension(num_crop_types) :: input_cover_types

 character(len=24) :: crp_name
 integer :: ipref
 character(len=256) :: text

 cropclock = mpp_clock_id('compute_crop_calendars', CLOCK_FLAG_DEFAULT, CLOCK_ROUTINE)

 call get_date(lnd%time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
 call get_date(lnd%time, year0,month0,day0,hour,minute,second)

 call read_crop_namelist

 do icrop=1,num_crop_types
 do iseason=1,num_crop_seasons
   restart_fieldname(iseason,icrop) = trim(crop_name(icrop))//'_'//trim(season_name(iseason))
 enddo
 enddo
 infile = ''
! Read the restart data
 infile = 'INPUT/'//trim(restart_file_name)
 call open_land_restart(restart,trim(infile),restart_exists)
 if(restart_exists) then
   call error_mesg('vegn_crop_init', 'reading NetCDF restart', NOTE)
   call get_tile_data(restart, 'tc_av_climate', 'month', vegn_tc_av_climate_ptr)                      !             +-- season number
   call get_tile_data(restart, 'precip_av_climate', 'month', vegn_precip_av_climate_ptr)              !             | +-- crop number
   call get_tile_data(restart, 'T_mid_mth', 'month', vegn_T_mid_mth_ptr)                              !             | |
   call get_tile_data(restart, 'P_mid_mth', 'month', vegn_P_mid_mth_ptr)                              !             v v
   call get_int_tile_data(restart, trim(restart_fieldname(1,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_1_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_2_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,3)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_3_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,3)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_3_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,4)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_4_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,4)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_4_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,5)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_5_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,5)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_5_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,6)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_6_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,6)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_6_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,7)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_7_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,7)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_7_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,8)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_8_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,8)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_8_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,9)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_9_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,9)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_9_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(1,10)),'plant_harvest', 'crop_cal_periods', crop_calendar_1_10_ptr)
   call get_int_tile_data(restart, trim(restart_fieldname(2,10)),'plant_harvest', 'crop_cal_periods', crop_calendar_2_10_ptr)
   call get_int_tile_data(restart, 'chosen_calendars', 'plant_harvest', 'crop_seasons', vegn_chosen_calendars_ptr)
   call get_int_tile_data(restart, 'potential_crop', 'crop_types',vegn_potential_crop_ptr)
   call get_int_tile_data(restart, 'chosen_crop', 'crop_seasons', vegn_chosen_crop_ptr)
   call get_int_tile_data(restart, 'chosen_crop_is_active', 'crop_seasons', vegn_chosen_crop_is_active_ptr)
   call get_int_tile_data(restart, 'grass_is_active', vegn_grass_is_active_ptr)
   ce = first_elmt(land_tile_map, ls=lnd%ls)
   do while(loop_over_tiles(ce,tile,L,k))
     if(.not.associated(tile%vegn)) cycle
     do iseason=1,num_crop_seasons
       if(tile%vegn%Crop%chosen_crop_is_active_int(iseason) == 0) then
         tile%vegn%Crop%chosen_crop_is_active(iseason) = .false.
       else if(tile%vegn%Crop%chosen_crop_is_active_int(iseason) == -1) then
         tile%vegn%Crop%chosen_crop_is_active(iseason) = .true.
       else
         call error_mesg('vegn_crop_init', 'chosen_crop_is_active is neither .true. or .false. This should never happen. Contact developer.', FATAL)
       endif
     enddo
     if(tile%vegn%Crop%grass_is_active_int == 0) then
       tile%vegn%Crop%grass_is_active = .false.
     else if(tile%vegn%Crop%grass_is_active_int == -1) then
       tile%vegn%Crop%grass_is_active = .true.
     else
       call error_mesg('vegn_crop_init', 'grass_is_active is neither .true. or .false. This should never happen. Contact developer.', FATAL)
     endif
   enddo

 else ! if(restart_exists) then
   call error_mesg('vegn_crop_init', 'cold starting vegn_crop_mod', NOTE)
   allocate(Tclim(lnd%ls:lnd%le,12), Pclim(lnd%ls:lnd%le,12))
   init_clim_exists = open_file(fileobj, "INPUT/initial_climatology.nc", "read")
   if (.not. init_clim_exists) then
     call error_mesg('vegn_crop_init', 'INPUT/initial_climatology.nc does not exist', FATAL)
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
     tile%vegn%Crop%chosen_crop_is_active = .FALSE.
     tile%vegn%Crop%chosen_crop_is_active_int = 0
   enddo
   deallocate(Tclim, Pclim)

! Read the MIRCA crop fractions. The crop with the largest fraction becomes the potential_crop.
! These crop fractions are the sum of irrigated and rainfed fractions from the MIRCA2000 data.
   allocate(MIRCA_crop_frac(lnd%ls:lnd%le,num_crop_types), crop_frac_tmp(lnd%ls:lnd%le,num_crop_types))
   allocate(potential_crop(lnd%ls:lnd%le,num_crop_types))
   dummyi = 0

   input_cover_types = (/ 1,-1,-1,-1,-1,-1,-1,-1,-1,-1/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,1) = crop_frac_tmp(:,1) ! IRRIGATED_MAIZE

   input_cover_types = (/-1, 2,-1,-1,-1,-1,-1,-1,-1,-1/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,2) = crop_frac_tmp(:,2) ! IRRIGATED_SOYBEAN

   input_cover_types = (/-1,-1, 3,-1,-1,-1,-1,-1,-1,-1/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,3) = crop_frac_tmp(:,3) ! IRRIGATED_RICE

   input_cover_types = (/-1,-1,-1, 4,-1,-1,-1,-1,-1,-1/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,4) = crop_frac_tmp(:,4) ! IRRIGATED_SPRING_WHEAT

   input_cover_types = (/-1,-1,-1,-1, 5,-1,-1,-1,-1,-1/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,5) = crop_frac_tmp(:,5) ! IRRIGATED_WINTER_WHEAT

   input_cover_types = (/-1,-1,-1,-1,-1, 6,-1,-1,-1,-1/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,6) = crop_frac_tmp(:,6) ! RAINFED_MAIZE

   input_cover_types = (/-1,-1,-1,-1,-1,-1, 7,-1,-1,-1/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,7) = crop_frac_tmp(:,7) ! RAINFED_SOYBEAN

   input_cover_types = (/-1,-1,-1,-1,-1,-1,-1, 8,-1,-1/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,8) = crop_frac_tmp(:,8) ! RAINFED_RICE

   input_cover_types = (/-1,-1,-1,-1,-1,-1,-1,-1, 9,-1/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,9) = crop_frac_tmp(:,9) ! RAINFED_SPRING_WHEAT

   input_cover_types = (/-1,-1,-1,-1,-1,-1,-1,-1,-1, 10/)
   call init_cover_field('single-tile', 'INPUT/crop_frac_by_water_source.nc', 'cover', 'crop_frac', lnd%sg_lonb, lnd%sg_latb, dummyi, input_cover_types, crop_frac_tmp)
   MIRCA_crop_frac(:,10) = crop_frac_tmp(:,10) ! RAINFED_WINTER_WHEAT

   do L=lnd%ls,lnd%le
     call compute_potential_crop(MIRCA_crop_frac(L,:), potential_crop(L,:))
   enddo

   ce = first_elmt(land_tile_map, ls=lnd%ls)
   do while(loop_over_tiles(ce,tile,L,k))
     if(.not.associated(tile%vegn)) cycle
     tile%vegn%Crop%potential_crop = potential_crop(L,:)
     call crop_selection(tile%vegn)
   enddo
   deallocate(MIRCA_crop_frac, crop_frac_tmp, potential_crop)
 endif ! if(restart_exists)

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
   if(icrop == IRRIGATED_MAIZE .or. icrop == RAINFED_MAIZE) then
     SI_crit(icrop) = SI_crit_Maize
     GP(icrop) = GP_Maize
     central_T(:,icrop)  = central_T_Maize
     central_P(:,icrop)  = central_P_Maize
     central_D(:,icrop)  = central_D_Maize
     variance_T(:,icrop) = variance_T_Maize
     variance_P(:,icrop) = variance_P_Maize
     variance_D(:,icrop) = variance_D_Maize
   else if(icrop == IRRIGATED_SOYBEAN .or. icrop == RAINFED_SOYBEAN) then
     SI_crit(icrop) = SI_crit_Soy
     GP(icrop) = GP_Soy
     central_T(:,icrop)  = central_T_Soy
     central_P(:,icrop)  = central_P_Soy
     central_D(:,icrop)  = central_D_Soy
     variance_T(:,icrop) = variance_T_Soy
     variance_P(:,icrop) = variance_P_Soy
     variance_D(:,icrop) = variance_D_Soy
   else if(icrop == IRRIGATED_RICE .or. icrop == RAINFED_RICE) then
     SI_crit(icrop) = SI_crit_Rice
     GP(icrop) = GP_Rice
     central_T(:,icrop)  = central_T_Rice
     central_P(:,icrop)  = central_P_Rice
     central_D(:,icrop)  = central_D_Rice
     variance_T(:,icrop) = variance_T_Rice
     variance_P(:,icrop) = variance_P_Rice
     variance_D(:,icrop) = variance_D_Rice
   else if(icrop == IRRIGATED_SPRING_WHEAT .or. icrop == RAINFED_SPRING_WHEAT) then
     GP(icrop) = 0
     SI_crit(icrop) = SI_crit_SW
     central_T(:,icrop)  = central_T_SW
     central_P(:,icrop)  = central_P_SW
     central_D(:,icrop)  = central_D_SW
     variance_T(:,icrop) = variance_T_SW
     variance_P(:,icrop) = variance_P_SW
     variance_D(:,icrop) = variance_D_SW
   else if(icrop == IRRIGATED_WINTER_WHEAT .or. icrop == RAINFED_WINTER_WHEAT) then
     GP(icrop) = 0
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
 potential_crop = NO_CROP
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
 subroutine compute_potential_crop(MIRCA_crop_frac, potential_crop)
 real,    intent(in)  :: MIRCA_crop_frac(num_crop_types)
 integer, intent(out) :: potential_crop(num_crop_types)
 real :: frac_tmp(num_crop_types), max_frac
 integer :: m, k, k_of_max_frac

 frac_tmp = MIRCA_crop_frac
 do m=1,num_crop_types
   max_frac = 0.0
   k_of_max_frac = 0
   do k=1,num_crop_types
     if(frac_tmp(k) > max_frac) then
       max_frac = frac_tmp(k)
       k_of_max_frac = k
     endif
   enddo
   if(k_of_max_frac == 0) then
     potential_crop(m) = NO_CROP
   else
     potential_crop(m) = k_of_max_frac
     frac_tmp(k_of_max_frac) = 0.0
   endif
 enddo

 end subroutine compute_potential_crop
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
 integer :: id_month, mth, id_crop_num, ical, id_season, icrop, iseason, iph, ipref, iperiod, nn, id_plant_harvest
 character(len=256) :: diag_fieldname

 id_month = diag_axis_init('month', (/(float(mth),mth=1,12)/),'none','Z','month of year')
 id_plant_harvest = diag_axis_init('plant_harvest', (/(float(iph),iph=1,2)/), 'none','Z','plant harvest')

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

 do ipref=1,num_crop_types
   do iseason=1,num_crop_seasons
     do iperiod=1,num_crop_periods
       do iph=1,2
         if(iph == 1) then
           diag_fieldname = trim(season_name(iseason))//'_'//trim(period_name(iperiod))//'_planting_date_for_crop_'
           nn = len_trim(diag_fieldname)
           write(diag_fieldname(nn+1:nn+2),'(i2.2)') ipref
         else
           diag_fieldname = trim(season_name(iseason))//'_'//trim(period_name(iperiod))//'_harvest_date_for_crop_'
           nn = len_trim(diag_fieldname)
           write(diag_fieldname(nn+1:nn+2),'(i2.2)') ipref
         endif
         id_crop_calendars(iph,iperiod,iseason,ipref) = register_tiled_diag_field(module_name,trim(diag_fieldname),(/id_ug/),lnd%time,trim(diag_fieldname),missing_value=0.0)
       enddo
     enddo
   enddo
 enddo

 do iseason=1,num_crop_seasons
   diag_fieldname = trim(season_name(iseason))//'_planting_and_harvest_dates'
   id_chosen_calendars(iseason) = register_tiled_diag_field(module_name,trim(diag_fieldname),(/id_ug,id_plant_harvest/),lnd%time,trim(diag_fieldname),missing_value=0.0)
 enddo

 id_potential_crop  = register_static_field(module_name,'potential_crop',(/id_ug,id_crop_num/), 'crops sorted by area in the MIRCA2000 data set',missing_value=0.0)
 call set_default_diag_filter('crop')
 id_chosen_crop = register_tiled_diag_field(module_name,'chosen_crop', (/id_ug,id_season/), lnd%time, 'number of chosen crop', missing_value=0.0)
 end subroutine crop_diag_init
 !======================================================================================================================================================
 subroutine save_crop_restart(tile_dim_length, timestamp)
 integer, intent(in) :: tile_dim_length
 character(*), intent(in) :: timestamp
 type(land_restart_type) :: restart
 character(len=256) :: filename
 integer :: idate, icrop, iseason, mth, iperiod, iph, L, k
 type(land_tile_enum_type) :: ce
 type(land_tile_type), pointer :: tile

 if(.not. crop_mod_initialized) call error_mesg('save_crop_restart','vegn_crop_init has not been called', FATAL)
 filename = 'RESTART/'//trim(timestamp)//trim(restart_file_name)
 call error_mesg('save_crop_restart', 'writing NetCDF restart "'//trim(filename)//'"', NOTE)
 call init_land_restart(restart, filename, vegn_tile_exists, tile_dim_length)
 call add_restart_axis(restart,'month', (/(float(mth),mth=1,12)/),.false.,longname='calendar month')
 call add_restart_axis(restart,'crop_cal_date',(/(float(idate),idate=1,num_crop_cal)/),.false.,longname='crop calendar dates as day of year')
 call add_restart_axis(restart,'crop_cal_periods',(/(float(iperiod),iperiod=1,num_crop_periods)/),.false.,longname='opt beg end')
 call add_restart_axis(restart,'plant_harvest',(/(float(iph),iph=1,2)/),.false.,longname='plant harvest')
 call add_restart_axis(restart,'crop_seasons',(/(float(iseason),iseason=1,num_crop_seasons)/),.false.,longname='crop season')
 call add_restart_axis(restart,'crop_types',(/(float(icrop),icrop=1,num_crop_types)/),.false.,longname='Maize, Soybean, Rice, Spring Wheat, Winter Wheat')
 call add_tile_data(restart,'tc_av_climate', 'month', vegn_tc_av_climate_ptr, 'climatological monthly average canopy air temperature','degK')
 call add_tile_data(restart,'precip_av_climate', 'month', vegn_precip_av_climate_ptr,'climatological monthly average precipitation rate','Kg/s*m^2')
 call add_tile_data(restart,'T_mid_mth', 'month', vegn_T_mid_mth_ptr, 'climatological average mid-month canopy air temperature','degK')
 call add_tile_data(restart,'P_mid_mth', 'month', vegn_P_mid_mth_ptr, 'climatological average mid-month precipitation rate','Kg/s*m^2')
 call add_int_tile_data(restart, trim(restart_fieldname(1,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,1)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_1_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,2)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_2_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,3)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_3_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,3)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_3_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,4)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_4_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,4)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_4_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,5)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_5_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,5)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_5_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,6)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_6_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,6)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_6_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,7)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_7_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,7)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_7_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,8)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_8_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,8)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_8_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,9)), 'plant_harvest', 'crop_cal_periods', crop_calendar_1_9_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,9)), 'plant_harvest', 'crop_cal_periods', crop_calendar_2_9_ptr, 'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(1,10)),'plant_harvest', 'crop_cal_periods', crop_calendar_1_10_ptr,'day of year: popt,hopt,pbeg,hbeg,pend,hend')
 call add_int_tile_data(restart, trim(restart_fieldname(2,10)),'plant_harvest', 'crop_cal_periods', crop_calendar_2_10_ptr,'day of year: popt,hopt,pbeg,hbeg,pend,hend')

 call add_int_tile_data(restart, 'chosen_calendars', 'plant_harvest', 'crop_seasons', vegn_chosen_calendars_ptr, 'popt, hopt, 1st and 2nd seasons','day of year')
 call add_int_tile_data(restart, 'potential_crop', 'crop_types',vegn_potential_crop_ptr, 'crop number of potential crop','crop number')
 call add_int_tile_data(restart, 'chosen_crop', 'crop_seasons', vegn_chosen_crop_ptr, 'crop number of chosen crop',      'crop number')
 ce = first_elmt(land_tile_map, ls=lnd%ls)
 do while(loop_over_tiles(ce,tile,L,k))
   if(.not.associated(tile%vegn)) cycle
   do iseason=1,num_crop_seasons
     if(tile%vegn%Crop%chosen_crop_is_active(iseason)) then
       tile%vegn%Crop%chosen_crop_is_active_int(iseason) = -1
     else
       tile%vegn%Crop%chosen_crop_is_active_int(iseason) = 0
     endif
   enddo
 enddo
 call add_int_tile_data(restart, 'chosen_crop_is_active', 'crop_seasons', vegn_chosen_crop_is_active_ptr, 'true when actively growing')
 call add_int_tile_data(restart, 'grass_is_active', vegn_grass_is_active_ptr, 'true when actively growing') ! Convert to integer
 call save_land_restart(restart)
 call free_land_restart(restart)
 end subroutine save_crop_restart
!=============================================================
 subroutine vegn_crop_end()
  crop_mod_initialized = .FALSE.
 end subroutine vegn_crop_end
!=============================================================
 logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
 end function vegn_tile_exists
!=============================================================
 function modulo_no_zero(nn,cycle_len) result(nn_within_cycle)
 integer, intent(in) :: nn,cycle_len
 integer :: nn_within_cycle

 nn_within_cycle = modulo(nn,cycle_len)
 if(nn_within_cycle == 0) nn_within_cycle = cycle_len
 end function modulo_no_zero
!=============================================================
subroutine vegn_T_mid_mth_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%T_mid_mth(n)
 endif
 end subroutine vegn_T_mid_mth_ptr
!=============================================================
subroutine vegn_P_mid_mth_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%P_mid_mth(n)
 endif
 end subroutine vegn_P_mid_mth_ptr
!=============================================================
subroutine vegn_tc_av_climate_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%tc_av_climate(n)
 endif
 end subroutine vegn_tc_av_climate_ptr
!=============================================================
subroutine vegn_precip_av_climate_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%precip_av_climate(n)
 endif
 end subroutine vegn_precip_av_climate_ptr
!=============================================================
!                         +-- season number
!                         | +-- crop number
!                         | |
!                         v v
 subroutine crop_calendar_1_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n             !              +-- plant, harvest
  integer,pointer::p                   !              | +-- period number (opt, beg, end)
  p=>NULL()                            !              | |
  if(associated(t))then                !              v v
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,1)
  endif
  end subroutine crop_calendar_1_1_ptr
!=============================================================
 subroutine crop_calendar_2_1_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,1)
  endif
  end subroutine crop_calendar_2_1_ptr
!=============================================================
 subroutine crop_calendar_1_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,2)
  endif
  end subroutine crop_calendar_1_2_ptr
!=============================================================
 subroutine crop_calendar_2_2_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,2)
  endif
  end subroutine crop_calendar_2_2_ptr
!=============================================================
 subroutine crop_calendar_1_3_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,3)
  endif
  end subroutine crop_calendar_1_3_ptr
!=============================================================
 subroutine crop_calendar_2_3_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,3)
  endif
  end subroutine crop_calendar_2_3_ptr
!=============================================================
 subroutine crop_calendar_1_4_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,4)
  endif
  end subroutine crop_calendar_1_4_ptr
!=============================================================
 subroutine crop_calendar_2_4_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,4)
  endif
  end subroutine crop_calendar_2_4_ptr
!=============================================================
 subroutine crop_calendar_1_5_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,5)
  endif
  end subroutine crop_calendar_1_5_ptr
!=============================================================
 subroutine crop_calendar_2_5_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,5)
  endif
  end subroutine crop_calendar_2_5_ptr
!=============================================================
 subroutine crop_calendar_1_6_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,6)
  endif
  end subroutine crop_calendar_1_6_ptr
!=============================================================
 subroutine crop_calendar_2_6_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,6)
  endif
  end subroutine crop_calendar_2_6_ptr
!=============================================================
 subroutine crop_calendar_1_7_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,7)
  endif
  end subroutine crop_calendar_1_7_ptr
!=============================================================
 subroutine crop_calendar_2_7_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,7)
  endif
  end subroutine crop_calendar_2_7_ptr
!=============================================================
 subroutine crop_calendar_1_8_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,8)
  endif
  end subroutine crop_calendar_1_8_ptr
!=============================================================
 subroutine crop_calendar_2_8_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,8)
  endif
  end subroutine crop_calendar_2_8_ptr
!=============================================================
 subroutine crop_calendar_1_9_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,9)
  endif
  end subroutine crop_calendar_1_9_ptr
!=============================================================
 subroutine crop_calendar_2_9_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,9)
  endif
  end subroutine crop_calendar_2_9_ptr
!=============================================================
 subroutine crop_calendar_1_10_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,1,10)
  endif
  end subroutine crop_calendar_1_10_ptr
!=============================================================
 subroutine crop_calendar_2_10_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%crop_calendars(m,n,2,10)
  endif
  end subroutine crop_calendar_2_10_ptr
!=============================================================
subroutine vegn_potential_crop_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 integer,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%potential_crop(n)
 endif
end subroutine vegn_potential_crop_ptr
!=============================================================
subroutine vegn_chosen_calendars_ptr(t,m,n,p)
  type(land_tile_type),pointer::t
  integer,intent(in):: m,n
  integer,pointer::p
  p=>NULL()
  if(associated(t))then
  if(associated(t%vegn))p=>t%vegn%Crop%chosen_calendars(m,n)
  endif
end subroutine vegn_chosen_calendars_ptr
!=============================================================
subroutine vegn_chosen_crop_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 integer,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%chosen_crop(n)
 endif
end subroutine vegn_chosen_crop_ptr
!=============================================================
subroutine vegn_chosen_crop_is_active_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 integer,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%chosen_crop_is_active_int(n)
 endif
end subroutine vegn_chosen_crop_is_active_ptr
!=============================================================
subroutine vegn_grass_is_active_ptr(t,p)
 type(land_tile_type),pointer::t
 integer,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%grass_is_active_int
 endif
end subroutine vegn_grass_is_active_ptr
!=============================================================
 end module vegn_crop_mod
