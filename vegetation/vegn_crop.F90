 module vegn_crop_mod

! This module contains routines that determine the most favorable dates for planting, and a range of dates
! when conditions are favorable, the most likely dates of harvest, and a range of likely harvest dates,
! for the four agricultural crops that cover the most land surface area. These are maize, soybean, wheat and rice.
! These dates, taken together, are referred to as a crop calendar.
! The crop calendars are determined separately for irrigated and rainfed crops, although
! they do not have to differ, and often do not in regions where precipitation is adequate.

! A single routine, subroutine CCA_Maize_Soy, handles both maize and soybean because the algorithm
! is identical for both these crops except for the numerical values of certain parameters.

! Subroutine CCA_Wheat executes the algorithm for wheat, but is called separately for spring and winter wheat
! because the algorithm differs for each type of wheat. The type of wheat is specified by an input argument to the routine.

! Subroutine CCA_Rice executes the algorithm for rice. Because rice is often double-cropped, CCA_Rice differs from
! the other routines in that it returns dates for main and second season crops instead of irrigated and rainfed.
! Therefore, CCA_Rice must be called twice, once each for irrigated and rainfed rice.
! The water source (irrigated or rainfed) is specified by an input argument to the routine.

! The crop calendars are computed for all tile types, whether LU_CROP or not.

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

 use mpp_mod, only: get_unit, mpp_pe
 use fms_mod, only: error_mesg, NOTE, WARNING, FATAL, file_exist, close_file, check_nml_error, stdlog
 use time_manager_mod, only: time_type, set_date, get_date, operator(-), set_time, operator(+), length_of_year, operator(//), operator(<), print_date, get_time
 use field_manager_mod, only: fm_field_name_len
 use constants_mod, only: TFREEZE, SECONDS_PER_DAY, PI
 use land_tile_mod, only: land_tile_type, land_tile_enum_type, first_elmt, loop_over_tiles, land_tile_map
 use vegn_tile_mod, only: vegn_tile_type, ITRUE
 use vegn_data_mod, only: LEAF_ON, LEAF_OFF, LU_CROP, MAIZE, SOYBEAN, RICE, SPRING_WHEAT, WINTER_WHEAT, NO_CROP
 use land_data_mod, only: lnd
 use land_tile_io_mod, only: land_restart_type, init_land_restart, open_land_restart, save_land_restart, &
                              free_land_restart, add_restart_axis, add_tile_data, get_tile_data, field_exists, add_int_tile_data, get_int_tile_data
 use astronomy_mod, only: get_orbital_parameters, get_ref_date_of_ae
 use diag_manager_mod, only: diag_axis_init, send_data, register_static_field, diag_field_add_attribute
 use land_tile_diag_mod,only: register_tiled_diag_field, send_tile_data, diag_buff_type, set_default_diag_filter
 use land_numerics_mod, only: ludcmp, lubksb
 use land_io_mod, only: init_cover_field
 use crop_debug_mod, only: debug_crop

 implicit none
 private

 public :: crop_init, crop_end, save_crop_restart, crop_calendar
 character(len=4), private, parameter :: module_name = 'crop'
 integer, parameter :: num_crop_types=5, num_watch=6, num_angles=3600
 integer, parameter :: days_in_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
 character(len=3), parameter :: month_name(12) = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)
 integer, parameter :: num_m = 4
 real, parameter :: aPTTtH_range_SW(2) = (/800.,967./)
 real, parameter :: aPTTtH_range_WW(2) = (/800.,837./)
 real, parameter :: unsuitable = 1000. ! suitablity index is set to unsuitable whenever it exeeds SI_crit. (Is this necessary?)
 real, parameter :: aPTT_interval = 200. ! Wheat suitablity is tested at intervals of aPTT_interval units of photo-thermal time.
 integer, parameter :: GP_Maize = 149, GP_Soy = 142, GP_Rice = 137 ! Time from planting to harvest
 character(len=10), parameter :: cwater(2) = (/' irrigated',' rainfed  '/)
 character(len=5), parameter :: cseason(2) = (/' main'," 2'nd"/)
 character(len=16), parameter :: restart_file_name = 'crop.nc'

!-----------------------------------------------------------
! used by routines that compute_day_length
 real :: ecc, obliq, per
 real :: orb_angle(0:num_angles)
 type(time_type) :: period_time_type, autumnal_eq_ref
!-----------------------------------------------------------

 real, parameter, dimension(0:num_m) :: central_T_Maize_orig_units  = (/18.96, 21.80, 23.41, 23.39, 21.12/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_Maize            = (/30.87, 13.87, 8.87, 7.92, 14.10/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_Maize_orig_units  = (/ 3.41, 4.26, 4.46, 4.34, 3.82/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_Maize_orig_units = (/ 3.37, 4.48, 4.49, 5.60, 5.62/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_Maize             = (/.5618, .5842, .5819, .5563, .5188/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_Maize            = (/.001449, .002331, .002202, .001230, .000514/) ! fraction^2

 real, parameter, dimension(0:num_m) :: central_T_Soy_orig_units  = (/21.19, 23.58, 24.35, 23.06, 19.63/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_Soy            = (/21.39, 7.78, 4.41, 7.08, 18.95/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_Soy_orig_units  = (/ 4.45, 4.87, 4.65, 4.02, 3.24/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_Soy_orig_units = (/ 3.75, 5.76, 4.56, 3.06, 3.32/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_Soy             = (/.5837, .5917, .5724, .5356, .4926/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_Soy            = (/.001053, .001243, .000980, .000597, .000441/) ! fraction^2

 real, parameter, dimension(0:num_m) :: central_T_SW_orig_units  = (/12.86, 15.96, 18.35, 20.31, 21.11/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_SW            = (/24.11, 8.33, 5.54, 4.18, 10.62/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_SW_orig_units  = (/ 1.59, 1.94, 1.89, 1.70, 1.49/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_SW_orig_units = (/ 0.58, 0.75, 0.89, 1.02, 1.14/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_SW             = (/.4755, .5291, .5519, .5514, .5331/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_SW            = (/.014516, .017204, .010947, .005253, .002207/) ! fraction^2

 real, parameter, dimension(0:num_m) :: central_T_WW_orig_units  = (/14.47, 12.38, 17.31, 20.49, 22.43/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_WW            = (/17.93, 5.30, 2.57, 2.91, 4.52/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_WW_orig_units  = (/ 1.80, 1.90, 2.04, 2.12, 2.16/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_WW_orig_units = (/ 0.71, 0.48, 0.74, 1.00, 1.43/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_WW             = (/.3992, .4994, .5711, .5860, .5834/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_WW            = (/.002711, .009490, .004064, .002138, .001261/) ! fraction^2

 real, parameter, dimension(0:num_m) :: central_T_Rice_orig_units  = (/23.77, 24.85, 25.90, 26.27, 25.07/) ! deg C
 real, parameter, dimension(0:num_m) :: variance_T_Rice            = (/32.11, 16.55, 9.44, 7.75, 13.79/) ! deg^2
 real, parameter, dimension(0:num_m) :: central_P_Rice_orig_units  = (/ 4.97, 5.43, 5.75, 5.22, 3.92/) ! mm/day
 real, parameter, dimension(0:num_m) :: variance_P_Rice_orig_units = (/14.42, 17.87, 15.77, 9.58, 8.50/) ! (mm/day)^2
 real, parameter, dimension(0:num_m) :: central_D_Rice             = (/.5306, .5364, .5325, .5198, .5022/) ! fraction of 24 hour day
 real, parameter, dimension(0:num_m) :: variance_D_Rice            = (/.001748, .002061, .001851, .001230, .000910/) ! fraction^2

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
 character(len=6) :: date_from_doy(0:365), date1, date2

 real :: t_mid_month(0:13)
! Note that t_mid_month is dimensioned (0:13) where t_mid_month(0) is negative because it's the middle of Dec
! of the previous year and t_mid_month(13) is > 365. because it's the middle of Jan of the following year.

 character(len=256) :: text
 integer :: id_T_ave, id_P_ave
 integer :: id_crop_cal_Maize, id_crop_cal_Soy, id_crop_cal_SW, id_crop_cal_WW, id_crop_cal_Rice_1, id_crop_cal_Rice_2
 integer :: id_MIRCA_crop_frac, id_current_crop

 real :: weight_climate=.10
 real :: Twt_Maize=1./3., Pwt_Maize=1./3.
 real :: Twt_Soy=1./3., Pwt_Soy=1./3.
 real :: Twt_SW=1./3., Pwt_SW=1./3.
 real :: Twt_WW=1./3., Pwt_WW=1./3.
 real :: Twt_Rice=1./3., Pwt_Rice=1./3.
 real :: SI_crit_Maize = 7.0
 real :: SI_crit_Soy = 9.0
 real :: SI_crit_SW = 7.0
 real :: SI_crit_WW = 13.0
 real :: SI_crit_Rice = 8.2
 real :: max_planting_SI_SW = 3.25
 real :: Tbase_Wheat = 5.0 + TFREEZE
 integer :: length_of_vernalization_period = 40
 real :: max_T_for_vernalization = 7.0 + TFREEZE
 real :: min_planting_T_Wheat = 5.0 + TFREEZE
 real :: absolute_min_T_for_Wheat = -7.0 + TFREEZE
 real :: lat_watch = 100.0, lon_watch = 400.0
 logical :: all_potential_growing_areas=.FALSE.
 character(len=9) :: crop_calendar_option = 'rainfed  ' ! valid options are 'rainfed' and 'irrigated'
 logical :: cold_start_idle = .FALSE.
 ! all_potential_growing_areas affects only diagnostic output.
 ! If .TRUE., then the diagnostic output will include crop calendars for all tiles,
 ! whether or not they are crop tiles and whether or not the MIRCA data shows any crop area.
 ! Regardless of the value of all_potential_growing_areas, only the calendar corresponding to the crop
 ! with the largest MIRCA area will be assigned to vegn%Crop%plant_beg, plant_opt, etc.

 namelist / vegn_crop_nml / weight_climate, &
            Twt_Maize, Twt_Soy, Twt_SW, Twt_WW, Twt_Rice, &
            Pwt_Maize, Pwt_Soy, Pwt_SW, Pwt_WW, Pwt_Rice, &
            SI_crit_Maize, SI_crit_Soy, SI_crit_SW, SI_crit_WW, SI_crit_Rice , &
            max_planting_SI_SW, Tbase_Wheat, length_of_vernalization_period, &
            max_T_for_vernalization, min_planting_T_Wheat, absolute_min_T_for_Wheat, &
            lat_watch, lon_watch, all_potential_growing_areas, crop_calendar_option, cold_start_idle
 contains
!============================================================================
 subroutine crop_calendar(vegn,diag,L)
 type(vegn_tile_type), intent(inout) :: vegn
 type(diag_buff_type), intent(inout) :: diag
 integer, intent(in) :: L ! index of grid cell which contains this tile
 integer :: n1, n2, second, minute, hour, day0, day1, month0, month1, year0, year1, watch_unit, iseason, mth, iwater
 logical :: new_month, compute_calendars
 real, dimension(12) :: rhs
 character(len=33) :: outname

 ! The routines for Maize, Soybean and Wheat return dates for irrigated crop
 ! as the first element of pday etc. and rainfed as the second element.
 ! CCA_Rice returns dates for the main and 2'nd season crops as the first element and second elements.
 ! CCA_Rice is called twice. Once each for irrigated and rainfed Rice.
 integer, dimension(2) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end

 if(.not.crop_mod_initialized) call error_mesg('crop_calendar','crop_init has not been called', FATAL)
 call get_date(lnd%time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
 call get_date(lnd%time, year0,month0,day0,hour,minute,second)
 new_month = month0 /= month1
 watch_unit = 0
 vegn%Crop%watchpoint = .false.
 if(vegn%landuse == LU_CROP) then
   if(abs(lon_watch - 180*lnd%ug_lon(L)/PI) < 0.7 .and. abs(lat_watch - 180*lnd%ug_lat(L)/PI) < 0.7) then
     outname = 'crop_calendar.      E.      N.out'
     write(outname(15:20),'(f6.2)') 180*lnd%ug_lon(L)/PI
     write(outname(23:28),'(f6.2)') 180*lnd%ug_lat(L)/PI
     do n1=1,4
       n2 = scan(outname,' ')
       if(n2 == 0) exit
       outname(n2:n2) = '0'
     enddo
     watch_unit = get_unit()
     open(unit=watch_unit, file=outname, action='write', form='formatted', position='append')
     text = ' L=       longitude=          latitude=        '
     write(text( 4: 8),'(i5)') L
     write(text(21:28),'(f8.3)') 180*lnd%ug_lon(L)/PI
     write(text(40:47),'(f8.3)') 180*lnd%ug_lat(L)/PI
     write(watch_unit,'(a)') ' Watchpoint: HelloZ '//trim(text)
     call error_mesg('Watchpoint: HelloZ',trim(text),NOTE)
     vegn%Crop%watchpoint = .true.
   endif
 endif
 if(new_month) then
    vegn%Crop%tc_av_climate(month1) = weight_climate*vegn%tc_av + (1-weight_climate)*vegn%Crop%tc_av_climate(month1)
    vegn%Crop%precip_av_climate(month1) = weight_climate*vegn%precip_av + (1-weight_climate)*vegn%Crop%precip_av_climate(month1)
    if(vegn%landuse == LU_CROP .and. watch_unit > 0) then
      write(watch_unit,'(2(a,3i4))') ' year1,month1,day1=',year1,month1,day1,' year0,month0,day0=',year0,month0,day0
      write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
      write(watch_unit,'(a,f7.3,a)') ' '//month_name(month1)//' Ave precipitation =',3401.6*days_in_month(month1)*vegn%Crop%precip_av_climate(month1),' inches'
      write(watch_unit,'(a,f7.3,a)') ' '//month_name(month1)//' Ave temperature   =',(vegn%Crop%tc_av_climate(month1)-273.15),'◦C'
    endif
 endif
 if(all_potential_growing_areas) then
   compute_calendars = new_month
 else
   compute_calendars = new_month .and. vegn%landuse == LU_CROP
 endif
 if(compute_calendars) then ! compute new climatological temperature and precipitation rates and new crop calendars once per month
    vegn%Crop%T_mid_mth = 4*vegn%Crop%tc_av_climate ! 4*tc_av_climate is the rhs. lubksb overwrites it with the solution.
    call lubksb(X_ludcmp, indx_ludcmp, vegn%Crop%T_mid_mth)
    vegn%Crop%P_mid_mth = 4*vegn%Crop%precip_av_climate ! 4*tcprecip_av_climate is the rhs. lubksb overwrites it with the solution.
    call lubksb(X_ludcmp, indx_ludcmp, vegn%Crop%P_mid_mth)
    if(watch_unit > 0) then
      write(watch_unit,'(a)') ''
      write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
      do mth=1,12
        write(watch_unit,'(a,2(f6.2,a))') ' '//month_name(mth)//' T_mid_mth =',vegn%Crop%T_mid_mth(mth)-273.15,'C  tc_av_climate=',vegn%Crop%tc_av_climate(mth)-273.15,'C'
      enddo
      do mth=1,12
        write(watch_unit,'(a,2(f6.2,a))') ' '//month_name(mth)//' P_mid_mth =',3401.6*days_in_month(mth)*vegn%Crop%P_mid_mth(mth),' inches  vegn%Crop%precip_av_climate=', &
                                             3401.6*days_in_month(mth)*vegn%Crop%precip_av_climate(mth),' inches'
      enddo
      if(vegn%Crop%current_crop == NO_CROP)      write(watch_unit,'(a)') ' current_crop = NO_CROP'
      if(vegn%Crop%current_crop == MAIZE)        write(watch_unit,'(a)') ' current_crop = MAIZE'
      if(vegn%Crop%current_crop == SOYBEAN)      write(watch_unit,'(a)') ' current_crop = SOYBEAN'
      if(vegn%Crop%current_crop == RICE)         write(watch_unit,'(a)') ' current_crop = RICE'
      if(vegn%Crop%current_crop == SPRING_WHEAT) write(watch_unit,'(a)') ' current_crop = SPRING_WHEAT'
      if(vegn%Crop%current_crop == WINTER_WHEAT) write(watch_unit,'(a)') ' current_crop = WINTER_WHEAT'
    endif
    if(all_potential_growing_areas .or. vegn%Crop%current_crop == MAIZE) then
      call CCA_Maize_Soy('Maize', watch_unit, vegn, L, central_T_Maize, variance_T_Maize, central_P_Maize, variance_P_Maize, &
                         central_D_Maize, variance_D_Maize, Twt_Maize, Pwt_Maize, SI_crit_Maize, pday, pday_beg, pday_end, hday, hday_beg, hday_end)
      if(vegn%Crop%idle == ITRUE) then
        ! The crop calendar should not be changed while the crop is actively growing.
        ! It should be changed only in the off season. It is the harvest date that is at issue here.
        ! Once a crop is planted one is committed to the harvest date that was determined at the time of planting.
        vegn%Crop%crop_cal_Maize(1: 6) = (/pday_beg(1), pday(1), pday_end(1), hday_beg(1), hday(1), hday_end(1)/)
        vegn%Crop%crop_cal_Maize(7:12) = (/pday_beg(2), pday(2), pday_end(2), hday_beg(2), hday(2), hday_end(2)/)
        call debug_crop(vegn,'HelloZ '//trim(crop_calendar_option)//' Maize cropland is idle')   ! debug
      else
        call debug_crop(vegn,'HelloZ '//trim(crop_calendar_option)//' Maize cropland is active') ! debug
      endif ! if(vegn%Crop%idle == ITRUE) then
      if(watch_unit > 0) then
        write(watch_unit,'(a)') ''
        write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
        write(watch_unit,'(a)') " 1'st and 3'rd columns are for irrigated Maize, 2'nd and 4'th for rainfed"
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(1)),'Hello01') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(7)),'Hello02')
        write(watch_unit,'(a,2i4,a)') ' begin planting period =',nint(vegn%Crop%crop_cal_Maize(1)),nint(vegn%Crop%crop_cal_Maize(7)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(2)),'Hello03') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(8)),'Hello04')
        write(watch_unit,'(a,2i4,a)') ' optimal planting date =',nint(vegn%Crop%crop_cal_Maize(2)),nint(vegn%Crop%crop_cal_Maize(8)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(3)),'Hello05') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(9)),'Hello06')
        write(watch_unit,'(a,2i4,a)') ' end   planting period =',nint(vegn%Crop%crop_cal_Maize(3)),nint(vegn%Crop%crop_cal_Maize(9)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(4)),'Hello07') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(10)),'Hello08')
        write(watch_unit,'(a,2i4,a)') ' begin harvest  period =',nint(vegn%Crop%crop_cal_Maize(4)),nint(vegn%Crop%crop_cal_Maize(10)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(5)),'Hello09') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(11)),'Hello10')
        write(watch_unit,'(a,2i4,a)') ' optimal harvest  date =',nint(vegn%Crop%crop_cal_Maize(5)),nint(vegn%Crop%crop_cal_Maize(11)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(6)),'Hello11') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Maize(12)),'Hello12')
        write(watch_unit,'(a,2i4,a)') ' end   harvest  period =',nint(vegn%Crop%crop_cal_Maize(6)),nint(vegn%Crop%crop_cal_Maize(12)),' '//date1//' '//date2
      endif
    endif ! if(all_potential_growing_areas .or. vegn%Crop%current_crop == MAIZE)
    if(all_potential_growing_areas .or. vegn%Crop%current_crop == SOYBEAN) then
      call CCA_Maize_Soy('Soybean', watch_unit, vegn, L, central_T_Soy, variance_T_Soy, central_P_Soy, variance_P_Soy, central_D_Soy, variance_D_Soy, &
                         Twt_Soy, Pwt_Soy, SI_crit_Soy, pday, pday_beg, pday_end, hday, hday_beg, hday_end)
      if(vegn%Crop%idle == ITRUE) then
        vegn%Crop%crop_cal_Soy(1: 6) = (/pday_beg(1), pday(1), pday_end(1), hday_beg(1), hday(1), hday_end(1)/)
        vegn%Crop%crop_cal_Soy(7:12) = (/pday_beg(2), pday(2), pday_end(2), hday_beg(2), hday(2), hday_end(2)/)
        call debug_crop(vegn,'HelloZ '//trim(crop_calendar_option)//' Soybean cropland is idle')   ! debug
      else
        call debug_crop(vegn,'HelloZ '//trim(crop_calendar_option)//' Soybean cropland is active') ! debug
      endif
      if(watch_unit > 0) then
        write(watch_unit,'(a)') ''
        write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
        write(watch_unit,'(a)') " 1'st and 3'rd columns are for irrigated Soybean, 2'nd and 4'th for rainfed"
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(1)),'Hello13') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(7)),'Hello14')
        write(watch_unit,'(a,2i4,a)') ' begin planting period =',nint(vegn%Crop%crop_cal_Soy(1)),nint(vegn%Crop%crop_cal_Soy(7)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(2)),'Hello15') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(8)),'Hello16')
        write(watch_unit,'(a,2i4,a)') ' optimal planting date =',nint(vegn%Crop%crop_cal_Soy(2)),nint(vegn%Crop%crop_cal_Soy(8)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(3)),'Hello17') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(9)),'Hello18')
        write(watch_unit,'(a,2i4,a)') ' end   planting period =',nint(vegn%Crop%crop_cal_Soy(3)),nint(vegn%Crop%crop_cal_Soy(9)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(4)),'Hello19') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(10)),'Hello20')
        write(watch_unit,'(a,2i4,a)') ' begin harvest  period =',nint(vegn%Crop%crop_cal_Soy(4)),nint(vegn%Crop%crop_cal_Soy(10)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(5)),'Hello21') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(11)),'Hello22')
        write(watch_unit,'(a,2i4,a)') ' optimal harvest  date =',nint(vegn%Crop%crop_cal_Soy(5)),nint(vegn%Crop%crop_cal_Soy(11)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(6)),'Hello23') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Soy(12)),'Hello24')
        write(watch_unit,'(a,2i4,a)') ' end   harvest  period =',nint(vegn%Crop%crop_cal_Soy(6)),nint(vegn%Crop%crop_cal_Soy(12)),' '//date1//' '//date2
      endif
    endif ! if(all_potential_growing_areas .or. vegn%Crop%current_crop == SOYBEAN)
    if(all_potential_growing_areas .or. vegn%Crop%current_crop == SPRING_WHEAT) then
      call CCA_Wheat(watch_unit, vegn, L, 'SW', central_T_SW, variance_T_SW, central_P_SW, variance_P_SW, central_D_SW, variance_D_SW, & ! intent(in)
                     Twt_SW, Pwt_SW, SI_crit_SW, max_planting_SI_SW, Tbase_Wheat, aPTTtH_range_SW, & ! intent(in)
                     length_of_vernalization_period, max_T_for_vernalization, min_planting_T_Wheat, & ! intent(in)
                     pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)
      if(vegn%Crop%idle == ITRUE) then
        vegn%Crop%crop_cal_SW(1: 6) = (/pday_beg(1), pday(1), pday_end(1), hday_beg(1), hday(1), hday_end(1)/)
        vegn%Crop%crop_cal_SW(7:12) = (/pday_beg(2), pday(2), pday_end(2), hday_beg(2), hday(2), hday_end(2)/)
        call debug_crop(vegn,'HelloZ '//trim(crop_calendar_option)//' SW cropland is idle')   ! debug
      else
        call debug_crop(vegn,'HelloZ '//trim(crop_calendar_option)//' SW cropland is active') ! debug
      endif
      if(watch_unit > 0) then
        write(watch_unit,'(a)') ''
        write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
        write(watch_unit,'(a)') " 1'st and 3'rd columns are for irrigated Spring Wheat, 2'nd and 4'th for rainfed"
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(1)),'Hello25') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(7)),'Hello26')
        write(watch_unit,'(a,2i4,a)') ' begin planting period =',nint(vegn%Crop%crop_cal_SW(1)),nint(vegn%Crop%crop_cal_SW(7)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(2)),'Hello27') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(8)),'Hello28')
        write(watch_unit,'(a,2i4,a)') ' optimal planting date =',nint(vegn%Crop%crop_cal_SW(2)),nint(vegn%Crop%crop_cal_SW(8)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(3)),'Hello29') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(9)),'Hello30')
        write(watch_unit,'(a,2i4,a)') ' end   planting period =',nint(vegn%Crop%crop_cal_SW(3)),nint(vegn%Crop%crop_cal_SW(9)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(4)),'Hello31') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(10)),'Hello32')
        write(watch_unit,'(a,2i4,a)') ' begin harvest  period =',nint(vegn%Crop%crop_cal_SW(4)),nint(vegn%Crop%crop_cal_SW(10)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(5)),'Hello33') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(11)),'Hello34')
        write(watch_unit,'(a,2i4,a)') ' optimal harvest  date =',nint(vegn%Crop%crop_cal_SW(5)),nint(vegn%Crop%crop_cal_SW(11)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(6)),'Hello35') ; date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_SW(12)),'Hello36')
        write(watch_unit,'(a,2i4,a)') ' end   harvest  period =',nint(vegn%Crop%crop_cal_SW(6)),nint(vegn%Crop%crop_cal_SW(12)),' '//date1//' '//date2
      endif
    endif ! if(all_potential_growing_areas .or. vegn%Crop%current_crop == SPRING_WHEAT)
    if(all_potential_growing_areas .or. vegn%Crop%current_crop == WINTER_WHEAT) then
      call CCA_Wheat(watch_unit, vegn, L, 'WW', central_T_WW, variance_T_WW, central_P_WW, variance_P_WW, central_D_WW, variance_D_WW, & ! intent(in)
                     Twt_WW, Pwt_WW, SI_crit_WW, max_planting_SI_SW, Tbase_Wheat, aPTTtH_range_WW, & ! intent(in)
                     length_of_vernalization_period, max_T_for_vernalization, min_planting_T_Wheat, & ! intent(in)
                     pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)
      if(vegn%Crop%idle == ITRUE) then
        vegn%Crop%crop_cal_WW(1: 6) = (/pday_beg(1), pday(1), pday_end(1), hday_beg(1), hday(1), hday_end(1)/)
        vegn%Crop%crop_cal_WW(7:12) = (/pday_beg(2), pday(2), pday_end(2), hday_beg(2), hday(2), hday_end(2)/)
        call debug_crop(vegn,'HelloZ '//trim(crop_calendar_option)//' WW cropland is idle')   ! debug
      else
        call debug_crop(vegn,'HelloZ '//trim(crop_calendar_option)//' WW cropland is active') ! debug
      endif
      if(watch_unit > 0) then
        write(watch_unit,'(a)') ''
        write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
        write(watch_unit,'(a)') " 1'st and 3'rd columns are for irrigated Winter Wheat, 2'nd and 4'th for rainfed"
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(1)),'Hello37'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(7)),'Hello38')
        write(watch_unit,'(a,2i4,a)') ' begin planting period =',nint(vegn%Crop%crop_cal_WW(1)),nint(vegn%Crop%crop_cal_WW(7)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(2)),'Hello39'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(8)),'Hello40')
        write(watch_unit,'(a,2i4,a)') ' optimal planting date =',nint(vegn%Crop%crop_cal_WW(2)),nint(vegn%Crop%crop_cal_WW(8)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(3)),'Hello41'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(9)),'Hello42')
        write(watch_unit,'(a,2i4,a)') ' end   planting period =',nint(vegn%Crop%crop_cal_WW(3)),nint(vegn%Crop%crop_cal_WW(9)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(4)),'Hello43'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(10)),'Hello44')
        write(watch_unit,'(a,2i4,a)') ' begin harvest  period =',nint(vegn%Crop%crop_cal_WW(4)),nint(vegn%Crop%crop_cal_WW(10)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(5)),'Hello45'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(11)),'Hello46')
        write(watch_unit,'(a,2i4,a)') ' optimal harvest  date =',nint(vegn%Crop%crop_cal_WW(5)),nint(vegn%Crop%crop_cal_WW(11)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(6)),'Hello47'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_WW(12)),'Hello48')
        write(watch_unit,'(a,2i4,a)') ' end   harvest  period =',nint(vegn%Crop%crop_cal_WW(6)),nint(vegn%Crop%crop_cal_WW(12)),' '//date1//' '//date2
      endif
    endif ! if(all_potential_growing_areas .or. vegn%Crop%current_crop == WINTER_WHEAT)
    if(all_potential_growing_areas .or. vegn%Crop%current_crop == RICE) then
      water_loop: do iwater=1,2
        call CCA_Rice(watch_unit, vegn, L, trim(cwater(iwater)), Twt_Rice, Pwt_Rice, central_T_Rice, variance_T_Rice, & ! intent(in)
                      central_P_Rice, variance_P_Rice, central_D_Rice, variance_D_Rice, SI_crit_Rice, & ! intent(in)
                      pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)
        if(vegn%Crop%idle == ITRUE) then
          if(iwater==1) then
            vegn%Crop%crop_cal_Rice_1(1:6) = (/pday_beg(1), pday(1), pday_end(1), hday_beg(1), hday(1), hday_end(1)/)
            vegn%Crop%crop_cal_Rice_2(1:6) = (/pday_beg(2), pday(2), pday_end(2), hday_beg(2), hday(2), hday_end(2)/)
          endif
          if(iwater==2) then
            vegn%Crop%crop_cal_Rice_1(7:12) = (/pday_beg(1), pday(1), pday_end(1), hday_beg(1), hday(1), hday_end(1)/)
            vegn%Crop%crop_cal_Rice_2(7:12) = (/pday_beg(2), pday(2), pday_end(2), hday_beg(2), hday(2), hday_end(2)/)
          endif
          call debug_crop(vegn,'HelloZ '//trim(cwater(iwater))//' Rice cropland is idle')   ! debug
        else
          call debug_crop(vegn,'HelloZ '//trim(cwater(iwater))//' Rice cropland is active') ! debug
        endif
      enddo water_loop

      if(watch_unit > 0) then
        write(watch_unit,'(a)') ''
        write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
        write(watch_unit,'(a)') ' main season rice'
        write(watch_unit,'(a)') " 1'st and 3'rd columns are for irrigated Rice, 2'nd and 4'th for rainfed"
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(1)),'Hello49'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(7)),'Hello50')
        write(watch_unit,'(a,2i4,a)') ' begin planting period =',nint(vegn%Crop%crop_cal_Rice_1(1)),nint(vegn%Crop%crop_cal_Rice_1(7)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(2)),'Hello51'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(8)),'Hello52')
        write(watch_unit,'(a,2i4,a)') ' optimal planting date =',nint(vegn%Crop%crop_cal_Rice_1(2)),nint(vegn%Crop%crop_cal_Rice_1(8)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(3)),'Hello53'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(9)),'Hello54')
        write(watch_unit,'(a,2i4,a)') ' end   planting period =',nint(vegn%Crop%crop_cal_Rice_1(3)),nint(vegn%Crop%crop_cal_Rice_1(9)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(4)),'Hello55'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(10)),'Hello56')
        write(watch_unit,'(a,2i4,a)') ' begin harvest  period =',nint(vegn%Crop%crop_cal_Rice_1(4)),nint(vegn%Crop%crop_cal_Rice_1(10)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(5)),'Hello57'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(11)),'Hello58')
        write(watch_unit,'(a,2i4,a)') ' optimal harvest  date =',nint(vegn%Crop%crop_cal_Rice_1(5)),nint(vegn%Crop%crop_cal_Rice_1(11)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(6)),'Hello59'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_1(12)),'Hello60')
        write(watch_unit,'(a,2i4,a)') ' end   harvest  period =',nint(vegn%Crop%crop_cal_Rice_1(6)),nint(vegn%Crop%crop_cal_Rice_1(12)),' '//date1//' '//date2
        write(watch_unit,'(a)') ''
        write(watch_unit,'(a)') ' second season rice'
        write(watch_unit,'(a)') " 1'st and 3'rd columns are for irrigated Rice, 2'nd and 4'th for rainfed"
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(1)),'Hello61'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(7)),'Hello62')
        write(watch_unit,'(a,2i4,a)') ' begin planting period =',nint(vegn%Crop%crop_cal_Rice_2(1)),nint(vegn%Crop%crop_cal_Rice_2(7)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(2)),'Hello63'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(8)),'Hello64')
        write(watch_unit,'(a,2i4,a)') ' optimal planting date =',nint(vegn%Crop%crop_cal_Rice_2(2)),nint(vegn%Crop%crop_cal_Rice_2(8)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(3)),'Hello65'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(9)),'Hello66')
        write(watch_unit,'(a,2i4,a)') ' end   planting period =',nint(vegn%Crop%crop_cal_Rice_2(3)),nint(vegn%Crop%crop_cal_Rice_2(9)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(4)),'Hello67'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(10)),'Hello68')
        write(watch_unit,'(a,2i4,a)') ' begin harvest  period =',nint(vegn%Crop%crop_cal_Rice_2(4)),nint(vegn%Crop%crop_cal_Rice_2(10)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(5)),'Hello69'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(11)),'Hello70')
        write(watch_unit,'(a,2i4,a)') ' optimal harvest  date =',nint(vegn%Crop%crop_cal_Rice_2(5)),nint(vegn%Crop%crop_cal_Rice_2(11)),' '//date1//' '//date2
        date1 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(6)),'Hello71'); date2 = get_date_from_doy(nint(vegn%Crop%crop_cal_Rice_2(12)),'Hello72')
        write(watch_unit,'(a,2i4,a)') ' end   harvest  period =',nint(vegn%Crop%crop_cal_Rice_2(6)),nint(vegn%Crop%crop_cal_Rice_2(12)),' '//date1//' '//date2
      endif
    endif ! if(all_potential_growing_areas .or. vegn%Crop%current_crop == RICE)
 endif ! if(compute_calendars)
 call set_crop_calendar(vegn) ! Select the crop calendar of chosen water source (specified via namelist switch crop_calendar_option)
 call send_tile_data(id_T_ave, vegn%Crop%tc_av_climate, diag)
 call send_tile_data(id_P_ave, vegn%Crop%precip_av_climate,diag)
 call send_tile_data(id_current_crop,real(vegn%Crop%current_crop),  diag)
 call send_tile_data(id_crop_cal_Maize,  vegn%Crop%crop_cal_Maize,  diag)
 call send_tile_data(id_crop_cal_Soy,    vegn%Crop%crop_cal_Soy,    diag)
 call send_tile_data(id_crop_cal_SW,     vegn%Crop%crop_cal_SW,     diag)
 call send_tile_data(id_crop_cal_WW,     vegn%Crop%crop_cal_WW,     diag)
 call send_tile_data(id_crop_cal_Rice_1, vegn%Crop%crop_cal_Rice_1, diag)
 call send_tile_data(id_crop_cal_Rice_2, vegn%Crop%crop_cal_Rice_2, diag)
 if(watch_unit > 0) close(watch_unit)

 end subroutine crop_calendar
!======================================================================================================================================================
 subroutine set_crop_calendar(vegn)
 type(vegn_tile_type), intent(inout) :: vegn
 integer :: offset

 if(trim(crop_calendar_option) == 'rainfed') then
   offset = 6
 else if(trim(crop_calendar_option) == 'irrigated') then
   offset = 0
 else
   call error_mesg('set_crop_calendar',trim(crop_calendar_option)//' is an invalid value of crop_calendar_option', FATAL)
 endif
 if(vegn%Crop%current_crop == MAIZE) then
   vegn%Crop%plant_beg = vegn%Crop%crop_cal_Maize(1+offset)
   vegn%Crop%plant_opt = vegn%Crop%crop_cal_Maize(2+offset)
   vegn%Crop%plant_end = vegn%Crop%crop_cal_Maize(3+offset)
   vegn%Crop%harvest_beg = vegn%Crop%crop_cal_Maize(4+offset)
   vegn%Crop%harvest_opt = vegn%Crop%crop_cal_Maize(5+offset)
   vegn%Crop%harvest_end = vegn%Crop%crop_cal_Maize(6+offset)
 else if(vegn%Crop%current_crop == SOYBEAN) then
   vegn%Crop%plant_beg = vegn%Crop%crop_cal_Soy(1+offset)
   vegn%Crop%plant_opt = vegn%Crop%crop_cal_Soy(2+offset)
   vegn%Crop%plant_end = vegn%Crop%crop_cal_Soy(3+offset)
   vegn%Crop%harvest_beg = vegn%Crop%crop_cal_Soy(4+offset)
   vegn%Crop%harvest_opt = vegn%Crop%crop_cal_Soy(5+offset)
   vegn%Crop%harvest_end = vegn%Crop%crop_cal_Soy(6+offset)
 else if(vegn%Crop%current_crop == SPRING_WHEAT) then
   vegn%Crop%plant_beg = vegn%Crop%crop_cal_SW(1+offset)
   vegn%Crop%plant_opt = vegn%Crop%crop_cal_SW(2+offset)
   vegn%Crop%plant_end = vegn%Crop%crop_cal_SW(3+offset)
   vegn%Crop%harvest_beg = vegn%Crop%crop_cal_SW(4+offset)
   vegn%Crop%harvest_opt = vegn%Crop%crop_cal_SW(5+offset)
   vegn%Crop%harvest_end = vegn%Crop%crop_cal_SW(6+offset)
 else if(vegn%Crop%current_crop == WINTER_WHEAT) then
   vegn%Crop%plant_beg = vegn%Crop%crop_cal_WW(1+offset)
   vegn%Crop%plant_opt = vegn%Crop%crop_cal_WW(2+offset)
   vegn%Crop%plant_end = vegn%Crop%crop_cal_WW(3+offset)
   vegn%Crop%harvest_beg = vegn%Crop%crop_cal_WW(4+offset)
   vegn%Crop%harvest_opt = vegn%Crop%crop_cal_WW(5+offset)
   vegn%Crop%harvest_end = vegn%Crop%crop_cal_WW(6+offset)
 else if(vegn%Crop%current_crop == RICE) then
   vegn%Crop%plant_beg = vegn%Crop%crop_cal_Rice_1(1+offset)
   vegn%Crop%plant_opt = vegn%Crop%crop_cal_Rice_1(2+offset)
   vegn%Crop%plant_end = vegn%Crop%crop_cal_Rice_1(3+offset)
   vegn%Crop%harvest_beg = vegn%Crop%crop_cal_Rice_1(4+offset)
   vegn%Crop%harvest_opt = vegn%Crop%crop_cal_Rice_1(5+offset)
   vegn%Crop%harvest_end = vegn%Crop%crop_cal_Rice_1(6+offset)
 else if(vegn%Crop%current_crop == NO_CROP) then
   vegn%Crop%plant_beg = 0.0
   vegn%Crop%plant_opt = 0.0
   vegn%Crop%plant_end = 0.0
   vegn%Crop%harvest_beg = 0.0
   vegn%Crop%harvest_opt = 0.0
   vegn%Crop%harvest_end = 0.0
 endif
 end subroutine set_crop_calendar
!======================================================================================================================================================
 subroutine CCA_Maize_Soy(crop_name, watch_unit, vegn, L, central_T, variance_T, central_P, variance_P, central_D, variance_D, Twt, Pwt, SI_crit, &
                          pday, pday_beg, pday_end, hday, hday_beg, hday_end)
 character(len=*), intent(in) :: crop_name
 integer, intent(in) :: watch_unit
 type(vegn_tile_type), intent(in) :: vegn
 integer, intent(in) :: L ! index of grid cell which contains this tile
 real, intent(in), dimension(0:) :: central_T, variance_T, central_P, variance_P, central_D, variance_D ! dimension is for months after planting. (0) is planting day.
 real, intent(in) :: Twt, Pwt, SI_crit
 integer, dimension(2), intent(out) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end

 integer :: k, km, kp, k2, k2m, k2p, daybeg, mths_after, doy, day, water_source, GP
 real :: Temp, Prec, TSI_test, DSI_test, dlen, Dwt, Ttmp, Dtmp, Ptmp(2)
 real :: PSI_test(2), SI_test(2) ! first element for irrigated, second for rainfed
 real :: annual_SI_max, annual_SI_min
 integer :: k_of_ann_SI_min
 real :: SI(num_test_days,2) ! Suitability Index. Computed at 5 day intervals from Jan 5 to Dec 31. SI(:,1) for irrigated, SI(:,2) for rainfed

 Dwt = 1.0 - Twt - Pwt
 if(vegn%landuse == LU_CROP .and. watch_unit > 0) then
   write(watch_unit,'(a)') ''
   write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
   write(watch_unit,'(a)') ' '//trim(crop_name)
   write(watch_unit,'(a)') ' doy  date   TSI  PSI_i PSI_r  DSI  SI_i  SI_r'
 endif
 ! compute SI at 5 day intervals from Jan 5 to Dec 31. SI_test is used for this.
 ! SI_test is not loaded into SI unless it passes the suitability test.
 k_loop_1: do k=1,num_test_days
   SI_test = 0.0
   daybeg = 5*k
   Ttmp = 0.0
   Dtmp = 0.0
   Ptmp = 0.0
   m_loop: do mths_after=0,num_m
     day = daybeg + 30*mths_after
     doy = modulo_no_zero(day,365)
     km = doy/5
     kp = km+1
     Temp = interp_between_mid_mths(doy, vegn%Crop%T_mid_mth)
     TSI_test = (Temp - central_T(mths_after))**2/variance_T(mths_after)
     Ttmp = Ttmp + Twt*TSI_test
     if(mths_after == num_m) then
       ! Month 4 is not tested for precip or day_length
       SI_test(1) = SI_test(1) + Twt*TSI_test
       SI_test(2) = SI_test(2) + Twt*TSI_test
     else
       Prec = interp_between_mid_mths(doy, vegn%Crop%P_mid_mth)
       PSI_test(1) = (max(Prec,central_P(mths_after)) - central_P(mths_after))**2/variance_P(mths_after)! Irrigation is equivalent to a minimum recipitation rate of central_P
       PSI_test(2) = (Prec - central_P(mths_after))**2/variance_P(mths_after)
       Ptmp = Ptmp + Pwt*PSI_test
       dlen = .2*((doy-5*km)*day_length(kp,L) + (5*kp-doy)*day_length(km,L))
       DSI_test = (dlen - central_D(mths_after))**2/variance_D(mths_after)
       Dtmp = Dtmp + Dwt*DSI_test
       SI_test(1) = SI_test(1) + Twt*TSI_test + Dwt*DSI_test + Pwt*PSI_test(1)
       SI_test(2) = SI_test(2) + Twt*TSI_test + Dwt*DSI_test + Pwt*PSI_test(2)
     endif
     if(SI_test(1) > SI_crit) then
       if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(i4,a,6f6.2,a)') 5*k,' '//get_date_from_doy(5*k,'Hello73'),Ttmp,Ptmp,Dtmp,SI_test(:),' (unsuitable)'
       exit m_loop ! If SI_test(1) exceeds critical, then so does SI_test(2)
     endif
   enddo m_loop
   SI(k,1) = SI_test(1)
   SI(k,2) = SI_test(2)
   if(SI(k,1) < SI_crit) then
     if(SI(k,2) > SI_crit) then
       if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(i4,a,6f6.2,a)') 5*k,' '//get_date_from_doy(5*k,'Hello74'),Ttmp,Ptmp,Dtmp,SI_test(:),' (only rainfed suitability index exceeds critical)'
     else
       if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(i4,a,6f6.2)') 5*k,' '//get_date_from_doy(5*k,'Hello75'),Ttmp,Ptmp,Dtmp,SI_test(:)
     endif
   endif
 enddo k_loop_1

 if(trim(crop_name) == 'Maize') then
   GP = GP_Maize
 else if(trim(crop_name) == 'Soybean') then
   GP = GP_Soy
 else
   call error_mesg('CCA_Maize_Soy in vegn_crop_mod',trim(crop_name)//' is not a valid value crop name', FATAL)
 endif 
 water_loop: do water_source=1,2
   annual_SI_min = HUGE(1.0)
   annual_SI_max = 0.0
   k_of_ann_SI_min = 0
   k_loop_2: do k=1,num_test_days ! Find the optimal time to plant this crop.
     if(SI(k,water_source) < annual_SI_min) then
       annual_SI_min = SI(k,water_source)
       k_of_ann_SI_min = k
     endif
     if(SI(k,water_source) > annual_SI_max) then
       annual_SI_max = SI(k,water_source)
     endif
   enddo k_loop_2
   if(annual_SI_min < SI_crit) then
     pday(water_source) = 5*k_of_ann_SI_min
   else
     pday(water_source) = 0
   endif
   if(annual_SI_min > SI_crit) then ! Find the range of dates over which conditions are suitable for planting
     ! SI is above SI_crit all year
     pday_beg(water_source) = 0
     pday_end(water_source) = 0
   else if(annual_SI_max < SI_crit) then
     ! SI is below SI_crit all year
     pday_beg(water_source) = 5
     pday_end(water_source) = 365
   else
     k_loop_3: do k=k_of_ann_SI_min-1,k_of_ann_SI_min-num_test_days+1,-1 ! Go back in time to find the first date where SI < SI_crit
       k2 = modulo_no_zero(k,num_test_days)
       if(SI(k2,water_source) > SI_crit) then
         k2p = modulo_no_zero(k2+1,num_test_days) ! move forward one step to get first date where SI < SI_crit
         pday_beg(water_source) = 5*k2p
         exit k_loop_3
       endif
     enddo k_loop_3
     k_loop_4: do k=k_of_ann_SI_min+1,k_of_ann_SI_min+num_test_days-1 ! Go forward in time to find the last date where SI < SI_crit
       k2 = modulo_no_zero(k,num_test_days)
       if(SI(k2,water_source) > SI_crit) then
         k2m = modulo_no_zero(k2-1,num_test_days) ! back up one step to get last date where SI < SI_crit
         pday_end(water_source) = 5*k2m
         exit k_loop_4
       endif
     enddo k_loop_4
   endif
   if(pday(water_source) > 0) then
     hday_beg(water_source) = modulo_no_zero(pday_beg(water_source) + GP, 365)
     hday(water_source)     = modulo_no_zero(pday(water_source)     + GP, 365)
     hday_end(water_source) = modulo_no_zero(pday_end(water_source) + GP, 365)
   else
     hday_beg(water_source) = 0
     hday(water_source)     = 0
     hday_end(water_source) = 0
   endif
 enddo water_loop
 end subroutine CCA_Maize_Soy
!======================================================================================================================================================
 subroutine CCA_Wheat(watch_unit, vegn, L, Wtype, central_T, variance_T, central_P, variance_P, central_D, variance_D, & ! intent(in)
                      Twt, Pwt, SI_crit, max_planting_SI_SW, Tbase, aPTTtH_range, & ! intent(in)
                      length_of_vernalization_period, max_T_for_vernalization, min_planting_T, & ! intent(in)
                      pday, pday_beg, pday_end, hday, hday_beg, hday_end) ! intent(out)
 integer, intent(in) :: watch_unit

! 1. Compute dates of accumulated photo-thermal time at intervals of 200 units, from zero to 800, for each
! candidate Optimal Planting Date (OPD) starting with Jan 5 and at five day intervals throughout the year.

! 2. If the accumulated photo-thermal time does not exceed 800 units starting from any date then the climate is deemed unsuitable for wheat.

! 3. Compute suitability index using climatic conditions at intervals of 200 units of accumulated photo-thermal time
! The suitability index is specific to the water source and variety: irrigated winter wheat, rainfed winter wheat, irrigated spring wheat, rainfed spring wheat

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
!===============================================================================================================================================================================
 type(vegn_tile_type), intent(inout) :: vegn
 integer :: L ! index of grid cell which contains this tile
 character(len=2) :: Wtype ! 'SW' or 'WW'
 real, intent(in) :: central_T(0:), variance_T(0:), central_P(0:), variance_P(0:), central_D(0:), variance_D(0:) ! At intervals of 200 aPTT units after planting. (0) is planting day.
 real, intent(in) :: Twt, Pwt, SI_crit, max_planting_SI_SW, Tbase, aPTTtH_range(2)
 integer, intent(in) :: length_of_vernalization_period
 real, intent(in) :: max_T_for_vernalization, min_planting_T ! default values are 7°C and 6°C
 integer, dimension(2), intent(out) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end

 character(len=12) :: cwater(2)
 integer :: k, k2, k2m, k2p, km, kp, daybeg, crossing_point, kautumn, k_of_ann_SI_min, kk, kkp, k_of_ann_SI_max, doy, pdayy, hdayy, water_source
 real :: Temp, Prec, TSI_test, DSI_test, aPTTtH, dlen, Dwt, Ttmp, Dtmp, Ptmp(2)
 real :: PSI_test(2), SI_test(2) ! first element for irrigated, second for rainfed
 real :: annual_SI_max, annual_SI_min
 real :: SI(num_test_days,2) ! Suitability Index. Computed at 5 day intervals from Jan 5 to Dec 31.
 integer :: crossing_day_400(num_test_days) ! Date at which accumulated photo-thermal time (aPTT) since planting reaches 400 units
 integer :: crossing_days(0:num_m)
 integer :: hday_list(num_test_days) ! Remember the values for each test date then choose the one that correspond to the annual minimum suitability index.
 real :: aPTTtH_list(num_test_days) ! Remember the values for each test date then choose the one that correspond to the annual minimum suitability index.
 logical :: passes_other_criteria

 if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(a)') ''
 cwater(1) = 'irrigated '//Wtype
 cwater(2) = 'rainfed '//Wtype
 Dwt = 1.0 - Twt - Pwt
 pday = 0
 hday = 0
 aPTTtH = 0.0
 hday_list = 0
 aPTTtH_list = 0.0
 ! compute SI at 5 day intervals from Jan 5 to Dec 31. SI_test is used for this.
 ! SI_test is not loaded into SI unless it passes the suitability test.
 if(vegn%landuse == LU_CROP .and. watch_unit > 0) then
   write(watch_unit,'(a)') ''
   write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
   write(watch_unit,'(a)') ' '//Wtype
   write(watch_unit,'(a)') ' doy  date   TSI  PSI_i PSI_r  DSI  SI_i  SI_r'
 endif
 k_loop_1: do k=1,num_test_days
   daybeg = 5*k
! Steps 1, 2, 6 and the harvest day of Step 5 are all handled within subroutine days_of_aPTT_crossings
   call days_of_aPTT_crossings(L, daybeg, num_m, Tbase, aPTT_interval, aPTTtH_range, vegn%Crop%T_mid_mth, & ! intent(in)
                               hday_list(k), aPTTtH_list(k), crossing_days) ! intent(out)
   if(any(crossing_days(:) == (/0,0,0,0,0/))) then
     SI(k,:) = unsuitable ! Steps 2, 6 and planting day of Step 5
     if(vegn%landuse == LU_CROP .and. watch_unit > 0) then
       write(watch_unit,'(i4,a)') 5*k,' '//get_date_from_doy(5*k,'Hello76')//' (unsuitable: aPTT never reaches 800 or if the temperature drops below -7°C before 800 units of aPTT is reached)'
     endif
     cycle k_loop_1 ! cycle k loop if aPTT never reaches 800 or if the temperature drops below -7°C before 800 units of aPTT is reached.
   endif
   crossing_day_400(k) = crossing_days(2)
   Ttmp = 0.0
   Dtmp = 0.0
   Ptmp = 0.0
   SI_test = 0.0
   do crossing_point=0,num_m ! Step 3
     doy = crossing_days(crossing_point)
     Temp = interp_between_mid_mths(doy, vegn%Crop%T_mid_mth)
     TSI_test = (Temp - central_T(crossing_point))**2/variance_T(crossing_point)
     Ttmp = Ttmp + Twt*TSI_test
     Prec = interp_between_mid_mths(doy, vegn%Crop%P_mid_mth)
     if(vegn%landuse == LU_CROP .and. watch_unit > 0) then
       write(watch_unit,'(i4,a,2(f6.2,a))') 200*crossing_point,' units of aPPT starting on '//get_date_from_doy(5*k,'Hello76.1')//' is reached on '// &
       get_date_from_doy(doy,'Hello76.2')//'  temperature on that date is',Temp-273.15,'C  precip rate on that date is',1.034e5*Prec,' inches/month'
     endif
     PSI_test(1) = (max(Prec,central_P(crossing_point)) - central_P(crossing_point))**2/variance_P(crossing_point) ! Irrigation is equivalent to a minimum recipitation rate of central_P
     PSI_test(2) = (Prec - central_P(crossing_point))**2/variance_P(crossing_point)
     Ptmp = Ptmp + Pwt*PSI_test
     km = doy/5
     kp = km+1
     dlen = .2*((doy-5*km)*day_length(kp,L) + (5*kp-doy)*day_length(km,L))
     DSI_test = (dlen - central_D(crossing_point))**2/variance_D(crossing_point)
     Dtmp = Dtmp + Dwt*DSI_test
     SI_test(1) = SI_test(1) + Twt*TSI_test + Dwt*DSI_test + Pwt*PSI_test(1)
     SI_test(2) = SI_test(2) + Twt*TSI_test + Dwt*DSI_test + Pwt*PSI_test(2)
     if(SI_test(1) > SI_crit) then
       if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(i4,a,6f6.2,a)') 5*k,' '//get_date_from_doy(5*k,'Hello77'),Ttmp,Ptmp,Dtmp,SI_test(:),' (unsuitable)'
       SI(k,:) = unsuitable ! If SI_test(1) exceeds critical, then so does SI_test(2)
       cycle k_loop_1
     endif
   enddo ! do crossing_point=0,num_m
   if(SI_test(2) < SI_crit) then
     if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(i4,a,6f6.2)') 5*k,' '//get_date_from_doy(5*k,'Hello78'),Ttmp,Ptmp,Dtmp,SI_test(:)
   else
     if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(i4,a,6f6.2,a)') 5*k,' '//get_date_from_doy(5*k,'Hello79'),Ttmp,Ptmp,Dtmp,SI_test(:),' (only rainfed suitability index exceeds critical)'
   endif
   SI(k,1) = SI_test(1) ! Conditions are suitable for planting irrigated Wheat on day of the year 5*k, provided it passes the tests in k_loop_2 and k_loop_3
   if(SI_test(2) > SI_crit) then
     SI(k,2) = unsuitable
   else
     SI(k,2) = SI_test(2) ! Conditions are suitable for planting rainfed Wheat on day of the year 5*k, provided it passes the tests in k_loop_2 and k_loop_3
   endif
 enddo k_loop_1

 water_loop: do water_source=1,2
   if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(a)')
   if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(a)') ' Algorithm details for '//trim(cwater(water_source))
   passes_other_criteria = .true.
   k_loop_2: do k=1,num_test_days ! Step 8: check that vernalization is possible for winter wheat and that temperature remains above 5C for spring wheat.
     if(SI(k,water_source) == unsuitable) cycle k_loop_2 ! Step 4 If the suitability index for the specific type of wheat being tested exceeds the critical value at
                                            ! all tested dates throughout the year then the tests within k_loop_2 are not necessary and will be skipped.
     if(Wtype == 'SW') then
       ! If the temperature drops below 5°C during the growing period then flag it as unsuitable for planting.
       if(T_goes_below_5C_during_GP(5*k, hday_list(k), vegn%Crop%T_mid_mth)) then
         SI(k,water_source) = unsuitable ! Step 8
         if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(a)') get_date_from_doy(5*k,'Hello80')//' is an unsuitable planting date for '//trim(cwater(water_source))//' because the temperature drops below 5C before the harvest date of '//get_date_from_doy(hday_list(k),'Hello81')
         passes_other_criteria = .false.
       endif
     endif
     if(Wtype == 'WW') then
       if(.not.vernalization_is_possible(5*k, crossing_day_400(k), vegn%Crop%T_mid_mth, length_of_vernalization_period, max_T_for_vernalization)) then
         SI(k,water_source) = unsuitable ! Step 8
         if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(a)') get_date_from_doy(5*k,'Hello82')//' is an unsuitable planting date for '//trim(cwater(water_source))//' because vernalization is not posible'
         passes_other_criteria = .false.
       endif
     endif
   enddo k_loop_2

   if(vegn%landuse == LU_CROP .and. watch_unit > 0) write(watch_unit,'(a)')
   k_loop_3: do k=1,num_test_days ! Do not plant when the temperature is below min_planting_T (default value is 5°C)
     if(SI(k,water_source) == unsuitable) cycle k_loop_3
     Temp = interp_between_mid_mths(5*k, vegn%Crop%T_mid_mth)
     if(Temp < min_planting_T) then
       SI(k,water_source) = unsuitable ! Step 7
       if(vegn%landuse == LU_CROP .and. watch_unit > 0) then
         write(watch_unit,'(a)') get_date_from_doy(5*k,'Hello83')//' is an unsuitable planting date for '//trim(cwater(water_source))//' because the temperature is below 5C on this date'
       endif
       passes_other_criteria = .false.
     endif
   enddo k_loop_3

   if(vegn%landuse == LU_CROP .and. watch_unit > 0) then
     write(watch_unit,'(a)')
     if(passes_other_criteria) then
        write(watch_unit,'(a)') ' No other planting criteria other than suitability index are violated for '//trim(cwater(water_source))
     else
       write(watch_unit,'(a)') ' After elimination of suitable planting dates above, the remaining suitability indices for '//trim(cwater(water_source))//' are:'
       do k=1,num_test_days
         if(SI(k,water_source) > SI_crit) then
           write(watch_unit,'(i4,a)') 5*k,' '//get_date_from_doy(5*k,'Hello84')//' unsuitable for '//trim(cwater(water_source))
         else
           write(watch_unit,'(i4,a,f6.2)') 5*k,' '//get_date_from_doy(5*k,'Hello85'),SI(k,water_source)
         endif
       enddo
     endif
   endif

   annual_SI_min = unsuitable
   pday(water_source) = 0
   hday(water_source) = 0
   aPTTtH = 0.0
   k_of_ann_SI_min = index_of_annual_SI_min(SI(:,water_source)) ! returns -1 if unsuitable all year
   k_of_ann_SI_max = index_of_annual_SI_max(SI(:,water_source)) ! never returns -1, returns the index of an unsuitable date if there are any.
   if(k_of_ann_SI_min > 0) then ! compute predicted planting and harvest days. Both remain zero if a planting date cannot be found.
     if(SI(k_of_ann_SI_min,water_source) /= unsuitable) then
       if(Wtype == 'SW') then
         annual_SI_min = SI(k_of_ann_SI_min,water_source)
         if(annual_SI_min > max_planting_SI_SW) then
           pday(water_source) = 5*k_of_ann_SI_min ! Step 9 If the annual minimum SI is above max_planting_SI_SW, then the predicted optimal planting date is the date of the minimum.
           hday(water_source) = hday_list(k_of_ann_SI_min) ! Step 9
           aPTTtH = aPTTtH_list(k_of_ann_SI_min)
         else
           k_loop_4: do k=k_of_ann_SI_min-1,k_of_ann_SI_min-36,-1
             ! Go back in time until SI reaches max_planting_SI_SW or until the end of the suitable period is reached.
             kk = modulo_no_zero(k, num_test_days)
             if(SI(kk,water_source) > max_planting_SI_SW .or. SI(kk,water_source) == unsuitable) then
               kkp = modulo_no_zero(k+1,num_test_days)
               pday(water_source) = 5*kkp ! Step 9 If the annual minimum SI is below 3.25 then the predicted optimal
                                     ! planting date is the date before the minimum when SI first drops below 3.25
               hday(water_source) = hday_list(kkp)
               aPTTtH = aPTTtH_list(kkp)
               exit k_loop_4
             endif
           enddo k_loop_4
         endif
       endif
       if(Wtype == 'WW') then
         annual_SI_min = SI(k_of_ann_SI_min,water_source)
         pday(water_source) = 5*k_of_ann_SI_min
         hday(water_source) = hday_list(k_of_ann_SI_min)
         aPTTtH = aPTTtH_list(k_of_ann_SI_min)
       endif
     endif
   endif

   annual_SI_max = SI(k_of_ann_SI_max,water_source)
   if(annual_SI_min > SI_crit) then ! Step 10 Find the range of dates over which conditions are suitable for planting
     ! SI is above SI_crit all year
     pday_beg(water_source) = 0
     pday_end(water_source) = 0
     hday_beg(water_source) = 0
     hday_end(water_source) = 0
   else if(annual_SI_max < SI_crit) then
     ! SI is below SI_crit all year
     pday_beg(water_source) = 5
     pday_end(water_source) = 365
     hday_beg(water_source) = 5
     hday_end(water_source) = 365
   else
     k_loop_5: do k=k_of_ann_SI_min-1,k_of_ann_SI_min-num_test_days+1,-1 ! Go back in time to find the first date where SI < SI_crit
       k2 = modulo_no_zero(k,num_test_days)
       if(SI(k2,water_source) > SI_crit) then
         k2p = modulo_no_zero(k2+1,num_test_days)
         pday_beg(water_source) = 5*k2p
         hday_beg(water_source) = hday_list(k2p)
         exit k_loop_5
       endif
     enddo k_loop_5
     k_loop_6: do k=k_of_ann_SI_min+1,k_of_ann_SI_min+num_test_days-1 ! Go forward in time to find the last date where SI < SI_crit
       k2 = modulo_no_zero(k,num_test_days)
       if(SI(k2,water_source) > SI_crit) then
         k2m = modulo_no_zero(k2-1,num_test_days)
         pday_end(water_source) = 5*k2m
         hday_end(water_source) = hday_list(k2m)
         exit k_loop_6
       endif
     enddo k_loop_6
   endif
 enddo water_loop

 end subroutine CCA_Wheat
!======================================================================================================================================================
 function vernalization_is_possible(pday, hday, T_mid_mth, length_of_vernalization_period, max_T_for_vernalization) result(It_is)
 integer, intent(in) :: pday, hday
 real, intent(in) :: T_mid_mth(12)
 integer, intent(in) :: length_of_vernalization_period
 real, intent(in) :: max_T_for_vernalization
 integer :: day, num_cold_days, hdayy
 logical :: It_is
 real :: Temp

 It_is = .false.
 if(hday < pday) then
   hdayy = hday + 365
 else
   hdayy = hday
 endif
 num_cold_days = 0
 day_loop: do day=pday,hdayy
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
 integer, intent(in) :: L ! index of grid cell which contains this tile
 integer, intent(in) :: daybeg, num_m
 real, intent(in) :: Tbase, aPTT_interval, aPTTtH_range(2)
 real, intent(in) :: T_mid_mth(12)

 ! crossing_days(m) = day of year when aPTT reaches m*aPTT_interval
 ! crossing_days(m) = zero if m*aPTT_interval is never reached or if temperature drops below -7C before aPTT reaches a value of num_m*aPTT_interval
 ! The date being tested is not suitable for planting either Spring or Winter Wheat if any of crossing_days(:) returned is zero

 integer, intent(out) :: harvestday
 real, intent(out) :: harvest_aPTT
 integer, intent(out) :: crossing_days(0:)

 integer :: day, m, doy, km, kp
 real :: TmTbase, aPTT_today, aPTT_yesterday, aPTT_target, T_today, dlen

 aPTT_yesterday = 0.0
 aPTT_target = 0.0
 harvest_aPTT = 0.0
 harvestday = 0
 crossing_days = 0
 m = -1
 day_loop: do day=daybeg,daybeg+364
   if(aPTT_yesterday > aPTTtH_range(2)) exit day_loop
   doy = modulo_no_zero(day,365)
   T_today = interp_between_mid_mths(doy, T_mid_mth)
   if(T_today < absolute_min_T_for_Wheat) exit day_loop ! Step 6 If the temperature never drops below -7°C before the harvest date then
                                    ! one or more elements of crossing_days will be returned with the full value.
   TmTbase = T_today - Tbase
   km = doy/5
   kp = km+1
   dlen = .2*((doy-5*km)*day_length(kp,L) + (5*kp-doy)*day_length(km,L))
   aPTT_today = aPTT_yesterday + dlen*max(0.0,TmTbase)
   if(aPTT_yesterday <= aPTT_target .and. aPTT_today > aPTT_target .and. m < num_m) then
     m = m + 1 ! Step 2 If the accumulated photo-thermal time does not exceed 800 units then m will not
               ! reach num_m and one or more elements of crossing_days will be returned with the full value.
     crossing_days(m) = doy ! Step 1
     aPTT_target = aPTT_target + aPTT_interval
   endif
   if(aPTT_today > aPTTtH_range(1) .and. aPTT_today < aPTTtH_range(2)) then
     harvestday = doy ! Step 5
     harvest_aPTT = aPTT_today
   endif
   aPTT_yesterday = aPTT_today
 enddo day_loop

 end subroutine days_of_aPTT_crossings
!======================================================================================================================================================
 subroutine CCA_Rice(watch_unit, vegn, L, water, Twt, Pwt, central_T, variance_T, central_P, variance_P, &
                                 central_D, variance_D, SI_crit, pday, pday_beg, pday_end, hday, hday_beg, hday_end)
 integer, intent(in) :: watch_unit
 type(vegn_tile_type), intent(in) :: vegn
 integer, intent(in) :: L ! index of grid cell which contains this tile
 character(len=*), intent(in) :: water ! Valid options are 'irrigated' 'rainfed'
 real, intent(in) :: Twt, Pwt
 real, intent(in) :: central_T(0:), variance_T(0:), central_P(0:), variance_P(0:), central_D(0:), variance_D(0:), SI_crit
 integer, dimension(2), intent(out) :: pday, pday_beg, pday_end, hday, hday_beg, hday_end
 character(len=256) :: mesg
 integer :: k, km, kp, daybeg, mths_after, doy, k_at_SI_min, ktest, iret, day, iseason
 real :: Temp, Prec, annual_SI_min, annual_SI_max, Dwt, dlen
 real, dimension(num_test_days) :: TSI, DSI, PSI, SI
!---------------------------------------------------------------------
 if(trim(water) == ' irrigated' .or. trim(water) == ' rainfed') then
   if(vegn%landuse == LU_CROP .and. watch_unit > 0) then
     write(watch_unit,'(a)') ''
     write(watch_unit,'(2(a,f7.3))') ' lon=',180*lnd%ug_lon(L)/PI,' lat=',180*lnd%ug_lat(L)/PI
     write(watch_unit,'(a)') ' '//trim(water)//' Rice'
     write(watch_unit,'(a)') ' doy  date   TSI   PSI   DSI   SI'
   endif
 else
   call error_mesg('CCA_Rice in vegn_crop_mod',trim(water)//' is not a valid value of water', FATAL)
 endif
 Dwt = 1.0 - Twt - Pwt
 k_loop: do k=1,num_test_days ! Compute the suitability index at 5 day intervals, starting with Jan 1
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
     TSI(k) = TSI(k) + Twt*(Temp - central_T(mths_after))**2/variance_T(mths_after)
     if(mths_after < 4) then
       ! Month 4 is not tested for precip or daylength
       Prec = interp_between_mid_mths(doy, vegn%Crop%P_mid_mth)
       if(trim(water) == ' irrigated') Prec = max(Prec,central_P(mths_after))
       PSI(k) = PSI(k) + Pwt*(Prec - central_P(mths_after))**2/variance_P(mths_after)
       km = doy/5
       kp = km+1
       dlen = .2*((doy-5*km)*day_length(kp,L) + (5*kp-doy)*day_length(km,L))
       DSI(k) = DSI(k) + Dwt*(dlen - central_D(mths_after))**2/variance_D(mths_after)
     endif
   enddo mths_loop
   SI(k) = TSI(k) + DSI(k) + PSI(k)
   if(vegn%landuse == LU_CROP .and. watch_unit > 0) then
     if(SI(k) > SI_crit) then
       write(watch_unit,'(i4,a,4f6.2,a)') 5*k,' '//get_date_from_doy(5*k,'Hello86'),TSI(k),PSI(k),DSI(k),SI(k),' (unsuitable)'
     else
       write(watch_unit,'(i4,a,4f6.2)') 5*k,' '//get_date_from_doy(5*k,'Hello87'),TSI(k),PSI(k),DSI(k),SI(k)
     endif
   endif
 enddo k_loop

 annual_SI_min = HUGE(1.0)
 do k=1,num_test_days ! Find the k index of the annual minimum of SI
   if(SI(k) < annual_SI_min) then
     annual_SI_min = SI(k)
     k_at_SI_min = k
   endif
 enddo

 if(annual_SI_min < SI_crit) then
   pday(1) = 5*k_at_SI_min
 else
   pday(:) = 0
   pday_beg(:) = 0
   pday_end(:) = 0
   return
 endif

 ktest = k_at_SI_min + 36
 if(ktest > num_test_days) ktest = ktest - num_test_days
 if(SI(ktest) < SI_crit) then
   pday(2) = 5*ktest
   do iseason=1,2
     call find_planting_date_range_Rice(pday(iseason), SI_crit, SI, 40, pday_beg(iseason), pday_end(iseason))
   enddo
 else
   call find_planting_date_range_Rice(pday(1), SI_crit, SI, 360, pday_beg(1), pday_end(1))
   pday(2) = 0
   pday_beg(2) = 0
   pday_end(2) = 0
 endif

 do iseason=1,2
   if(pday(iseason) > 0) then
     hday_beg(iseason) = modulo_no_zero(pday_beg(iseason) + GP_Rice, 365)
     hday(iseason)     = modulo_no_zero(pday(iseason)     + GP_Rice, 365)
     hday_end(iseason) = modulo_no_zero(pday_end(iseason) + GP_Rice, 365)
   else
     hday_beg(iseason) = 0
     hday(iseason)     = 0
     hday_end(iseason) = 0
   endif
 enddo

 end subroutine CCA_Rice
!======================================================================================================================================================
 subroutine find_planting_date_range_Rice(pday, SI_crit, SI, max_range_length, pday_beg, pday_end)
 integer, intent(in) :: pday
 real, intent(in) :: SI_crit, SI(num_test_days)
 integer, intent(in) :: max_range_length
 integer, intent(out) :: pday_beg, pday_end
 character(len=256) :: mesg
 integer :: k, k2, k2m, k2p, k_at_pday, krange

 k_at_pday = pday/5
 if(k_at_pday < 1 .or. k_at_pday > num_test_days) then
   mesg = 'ERROR1 in subroutine find_planting_date_range_Rice: invalid value of pday. pday='
   write(mesg(78:83),'(i6)') pday
   call error_mesg('find_planting_date_range_Rice in vegn_crop_mod',trim(mesg), FATAL)
 endif
 if(SI(k_at_pday) > SI_crit) then
   mesg = 'ERROR2 in subroutine find_planting_date_range_Rice: suitability index at planting day exceeds SI_crit'
   call error_mesg('find_planting_date_range_Rice in vegn_crop_mod',trim(mesg), FATAL)
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
 end subroutine find_planting_date_range_Rice
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
! Note that t_mid_month is dimensioned (0:13) where t_mid_month(0) is negative because it's the middle of Dec
! of the previous year and t_mid_month(13) is > 365. because it's the middle of Jan of the following year.
 w0 = ( day_of_year_local - t_mid_month(mon)) / (t_mid_month(mon+1) - t_mid_month(mon))
 w1 = (t_mid_month(mon+1) - day_of_year_local) / (t_mid_month(mon+1) - t_mid_month(mon))
 interp_between_mid_mths = w1*tmp(mon) + w0*tmp(mon+1)
 end function interp_between_mid_mths
!======================================================================================================================================================
 subroutine read_crop_namelist
 integer :: outunit, io, ierr

#ifdef INTERNAL_FILE_NML
    read(input_nml_file, nml=vegn_crop_nml, iostat=io)
    ierr = check_nml_error(io, 'vegn_crop_nml')
#else
  if (file_exist('input.nml')) then
     outunit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (outunit, nml=vegn_crop_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'vegn_crop_nml')
     enddo
10   continue
     call close_file(outunit)
  endif
#endif
  outunit = stdlog()
  write(outunit, nml=vegn_crop_nml)
 end subroutine read_crop_namelist
!======================================================================================================================================================
 subroutine crop_init(id_ug)
 integer, intent(in) :: id_ug
 type(land_restart_type) :: restart
 logical :: restart_exists, used
 real :: Dmm, Dm0, Dmp
 real, dimension(12,12) :: X
 real :: cosz, fracday1, fracday2, rrsun, max_frac
 integer :: dummyi, L, m, k, dom, doy ! L = index of grid cell which contains this tile
 type(land_tile_enum_type) :: ce
 type(land_tile_type), pointer :: tile
 character(len=3) :: month_name(12) = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)
 character(len=64) :: outfile
 integer :: ierr, outunit
 integer :: day_ae, month_ae, year_ae, hour_ae, minute_ae, second_ae
 real, allocatable, dimension(:,:) :: MIRCA_crop_frac, crop_frac_tmp
 integer, allocatable, dimension(:) :: current_crop

 integer :: i,j
 character(len=256) :: text

 call read_crop_namelist
!------------------------------------------------------------------------------------------------------------------------------------------------------
! Read the restart data
 text = 'INPUT/'//trim(restart_file_name)
 call open_land_restart(restart,trim(text),restart_exists)
 if(restart_exists) then
   call error_mesg('crop_init', 'reading NetCDF restart', NOTE)
   call get_tile_data(restart, 'T_mid_mth', 'month', vegn_T_mid_mth_ptr)
   call get_tile_data(restart, 'P_mid_mth', 'month', vegn_P_mid_mth_ptr)
   call get_tile_data(restart, 'tc_av_climate', 'month', vegn_tc_av_climate_ptr)
   call get_tile_data(restart, 'precip_av_climate', 'month', vegn_precip_av_climate_ptr)
   call get_int_tile_data(restart, 'current_crop', vegn_current_crop_ptr)
   call get_tile_data(restart, 'crop_cal_Maize', 'crop_cal', vegn_crop_cal_Maize_ptr)
   call get_tile_data(restart, 'crop_cal_Soy', 'crop_cal', vegn_crop_cal_Soy_ptr)
   call get_tile_data(restart, 'crop_cal_SW', 'crop_cal', vegn_crop_cal_SW_ptr)
   call get_tile_data(restart, 'crop_cal_WW', 'crop_cal', vegn_crop_cal_WW_ptr)
   call get_tile_data(restart, 'crop_cal_Rice_1', 'crop_cal', vegn_crop_cal_Rice_1_ptr)
   call get_tile_data(restart, 'crop_cal_Rice_2', 'crop_cal', vegn_crop_cal_Rice_2_ptr)
   if(cold_start_idle) then
     call error_mesg('crop_init', 'cold starting Crop%idle', NOTE)
    ce = first_elmt(land_tile_map, ls=lnd%ls)
    do while(loop_over_tiles(ce,tile,L,k))
      if(.not.associated(tile%vegn)) cycle
      tile%vegn%Crop%idle = ITRUE
    enddo
   else
     call get_int_tile_data(restart, 'idle', vegn_idle_ptr)
   endif
 else
   call error_mesg('crop_init', 'cold starting vegn_crop_mod', NOTE)
   ce = first_elmt(land_tile_map, ls=lnd%ls)
   do while(loop_over_tiles(ce,tile,L,k))
     if(.not.associated(tile%vegn)) cycle
     tile%vegn%Crop%T_mid_mth = 288.
     tile%vegn%Crop%P_mid_mth = 3.0e-5
     tile%vegn%Crop%tc_av_climate = 288.
     tile%vegn%Crop%precip_av_climate = 3.0e-5
     tile%vegn%Crop%current_crop = NO_CROP
     tile%vegn%Crop%crop_cal_Maize = 0.0
     tile%vegn%Crop%crop_cal_Soy = 0.0
     tile%vegn%Crop%crop_cal_SW = 0.0
     tile%vegn%Crop%crop_cal_WW = 0.0
     tile%vegn%Crop%crop_cal_Rice_1 = 0.0
     tile%vegn%Crop%crop_cal_Rice_2 = 0.0
     tile%vegn%Crop%idle = ITRUE
   enddo
 endif

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
 if(ierr /= 0) call error_mesg('crop_init in vegn_crop_mod','Error in subroutine ludcmp. Probably because matrix is singular.', FATAL)
!------------------------------------------------------------------------------------------------------------------------------------------------------
! compute day length at 5 day intervals
 call get_orbital_parameters(ecc, obliq, per)
 call get_ref_date_of_ae(day_ae, month_ae, year_ae, second_ae, minute_ae, hour_ae)
 autumnal_eq_ref = set_date(year_ae, month_ae, day_ae, hour_ae, minute_ae, second_ae)
 period_time_type = length_of_year()
 call orbit

! extend day_length array one step beyond either end of the year to facilitate interpolation between 5 day intervals
 allocate(day_length(0:num_test_days+1,lnd%ls:lnd%le))

 ce = first_elmt(land_tile_map, ls=lnd%ls)
 do while(loop_over_tiles(ce,tile,L))
   do k=1,num_test_days
     call compute_day_length(year_ae, 5*k, lnd%ug_lat(L), day_length(k,L))
   enddo
   day_length(0,L) = day_length(num_test_days,L)
   day_length(num_test_days+1,L) = day_length(1,L)
 enddo
!------------------------------------------------------------------------------------------------------------------------------------------------------
! compute coefficients used by function interp_between_mid_mths
 t_mid_month(0) = -.5*days_in_month(12)
 t_mid_month(1) = .5*days_in_month(1)
 do m=2,12
   t_mid_month(m) = t_mid_month(m-1) + .5*(days_in_month(m-1)+days_in_month(m))
 enddo
 t_mid_month(13) = t_mid_month(12) + .5*(days_in_month(12)+days_in_month(1))
!------------------------------------------------------------------------------------------------------------------------------------------------------
! Convert units of means and variances
 central_T_Maize = central_T_Maize_orig_units + TFREEZE
 central_P_Maize = central_P_Maize_orig_units/SECONDS_PER_DAY
 variance_P_Maize = variance_P_Maize_orig_units/SECONDS_PER_DAY**2
 central_T_Soy = central_T_Soy_orig_units + TFREEZE
 central_P_Soy = central_P_Soy_orig_units/SECONDS_PER_DAY
 variance_P_Soy = variance_P_Soy_orig_units/SECONDS_PER_DAY**2
 central_T_SW = central_T_SW_orig_units + TFREEZE
 central_P_SW = central_P_SW_orig_units/SECONDS_PER_DAY
 variance_P_SW = variance_P_SW_orig_units/SECONDS_PER_DAY**2
 central_T_WW = central_T_WW_orig_units + TFREEZE
 central_P_WW = central_P_WW_orig_units/SECONDS_PER_DAY
 variance_P_WW = variance_P_WW_orig_units/SECONDS_PER_DAY**2
 central_T_Rice = central_T_Rice_orig_units + TFREEZE
 central_P_Rice = central_P_Rice_orig_units/SECONDS_PER_DAY
 variance_P_Rice = variance_P_Rice_orig_units/SECONDS_PER_DAY**2
!------------------------------------------------------------------------------------------------------------------------------------------------------
! date_from_doy is returned by function get_date_from_doy.
! function get_date_from_doy is used only to facilitate watchpoint output.
 doy = 0
 date_from_doy(0) = '------'
 do m=1,12
   do dom=1,days_in_month(m)
     doy = doy + 1
     date_from_doy(doy)(1:3) = month_name(m)
     write(date_from_doy(doy)(4:6),'(i3)') dom
   enddo
 enddo
!------------------------------------------------------------------------------------------------------------------------------------------------------
! Read the MIRCA crop fractions. The crop with the largest fraction becomes the current_crop.
! These crop fractions are the sum of irrigated and rainfed fractions from the MIRCA2000 data.
 allocate(MIRCA_crop_frac(lnd%ls:lnd%le,num_crop_types), crop_frac_tmp(lnd%ls:lnd%le,num_crop_types))
 allocate(current_crop(lnd%ls:lnd%le))
!outfile = 'ls_le_pe.    .out' ! debug
!write(outfile(10:13),'(i4.4)') mpp_pe() ! debug
!outunit = get_unit() ! debug
!open(unit=outunit, file=trim(outfile), action='write', form='formatted') ! debug
!write(outunit,'(2(a,i6))') 'lnd%ls=',lnd%ls,' lnd%le=',lnd%le ! debug
!close(outunit) ! debug
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

 do L=lnd%ls,lnd%le
   max_frac = 0.0
   do k=1,5
     if(MIRCA_crop_frac(L,k) > max_frac) then
       max_frac = MIRCA_crop_frac(L,k)
       current_crop(L) = k
     endif
   enddo
   if(max_frac == 0.0) current_crop(L) = 0
 enddo

 ce = first_elmt(land_tile_map, ls=lnd%ls)
 do while(loop_over_tiles(ce,tile,L,k))
   if(.not.associated(tile%vegn)) cycle
   tile%vegn%Crop%current_crop = current_crop(L)
   call set_crop_calendar(tile%vegn)
 enddo

 call crop_diag_init(id_ug)
 if(id_MIRCA_crop_frac>0) used = send_data(id_MIRCA_crop_frac, MIRCA_crop_frac, lnd%time)
 deallocate(current_crop, crop_frac_tmp, MIRCA_crop_frac)

!outfile = 'lat_lon_pe.    .out' ! debug
!write(outfile(12:15),'(i4.4)') mpp_pe() ! debug
!outunit = get_unit() ! debug
!open(unit=outunit, file=trim(outfile), action='write', form='formatted') ! debug
!do j=lnd%js,lnd%je ! debug
!do i=lnd%is,lnd%ie ! debug
!  text = 'face= , i=  , j=  , lon=        , lat=       ' ! debug
!  write(text( 6: 6),'(i1)') lnd%sg_face ! debug
!  write(text(11:12),'(i2)') i ! debug
!  write(text(17:18),'(i2)') j ! debug
!  write(text(25:32),'(f8.3)') lnd%sg_lon(i,j) ! debug
!  write(text(39:45),'(f7.3)') lnd%sg_lat(i,j) ! debug
!  write(outunit,'(a)') trim(text) ! debug
!enddo ! debug
!enddo ! debug
!close(outunit) ! debug

 crop_mod_initialized = .TRUE.
 end subroutine crop_init
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

 ! Solve via Runge-Kutta
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
 function get_date_from_doy(doy,msgid) result(date_string)
 integer, intent(in) :: doy
 character(len=*), intent(in) :: msgid
 character(len=6) :: date_string
 character(len=128) :: error_message

 if(doy < 0 .or. doy > 365) then
   error_message = 'Day of year is out of bounds. doy=                 '//trim(msgid)
   write(error_message(35:50),'(i12)') doy
   call error_mesg('get_date_from_doy',trim(error_message), FATAL)
 endif
 date_string = date_from_doy(doy)
 end function get_date_from_doy
!======================================================================================================================================================
 subroutine crop_diag_init(id_ug)
 integer, intent(in) :: id_ug
 integer :: id_month, mth, id_crop_cal, id_crop_num, ical, iseason, icrop
 integer :: year,month,day,hour,minute,second

 id_month = diag_axis_init('month', (/(float(mth),mth=1,12)/),'none','Z','month of year')
 id_crop_cal = diag_axis_init('crop_cal',(/(float(ical),ical=1,12)/),'none','Z','plant beg, plant optimal, plant end, harvest beg, harvest optimal, harvest end')
 id_crop_num = diag_axis_init('crop_num',(/(float(icrop),icrop=1,num_crop_types)/),'none','Z','Maize, Soybean, Rice, Spring Wheat, Winter Wheat')
 call set_default_diag_filter('crop')
 id_T_ave = register_tiled_diag_field(module_name,'tc_av_climate', (/id_ug,id_month/),lnd%time,'climatological monthly mean temperature','deg K', missing_value=-1.0)
 id_P_ave = register_tiled_diag_field(module_name,'precip_av_climate',(/id_ug,id_month/),lnd%time,'climatological monthly mean precipitation','Kg/s*m^2', missing_value=-1.0)
 id_current_crop    = register_tiled_diag_field(module_name,'current_crop', (/id_ug/), lnd%time, 'crop of maximum area', missing_value= 0.0)
 id_crop_cal_Maize  = register_tiled_diag_field(module_name,'crop_cal_Maize', (/id_ug,id_crop_cal/),lnd%time,'maize crop calendar', missing_value= 0.0)
 id_crop_cal_Soy    = register_tiled_diag_field(module_name,'crop_cal_Soy', (/id_ug,id_crop_cal/),lnd%time,'soybean crop calendar', missing_value= 0.0)
 id_crop_cal_SW     = register_tiled_diag_field(module_name,'crop_cal_SW', (/id_ug,id_crop_cal/),lnd%time,'spring wheat crop calendar', missing_value= 0.0)
 id_crop_cal_WW     = register_tiled_diag_field(module_name,'crop_cal_WW', (/id_ug,id_crop_cal/),lnd%time,'winter wheat crop calendar', missing_value= 0.0)
 id_crop_cal_Rice_1 = register_tiled_diag_field(module_name,'crop_cal_Rice_1',(/id_ug,id_crop_cal/),lnd%time,'rice crop calendar. Main crop.', missing_value= 0.0)
 id_crop_cal_Rice_2 = register_tiled_diag_field(module_name,'crop_cal_Rice_2',(/id_ug,id_crop_cal/),lnd%time,'rice crop calendar. Second crop.',missing_value= 0.0)
 id_MIRCA_crop_frac = register_static_field(module_name,'MIRCA_crop_frac',(/id_ug,id_crop_num/),'crop fraction','unitless')
 call diag_field_add_attribute(id_MIRCA_crop_frac,'ocean_fillvalue',0.0)

 end subroutine crop_diag_init
!======================================================================================================================================================
 subroutine save_crop_restart(tile_dim_length, timestamp)
 integer, intent(in) :: tile_dim_length
 character(*), intent(in) :: timestamp
 type(land_restart_type) :: restart
 character(len=256) :: filename
 integer :: mn

 if(.not. crop_mod_initialized) call error_mesg('save_crop_restart','crop_init has not been called', FATAL)
 filename = 'RESTART/'//trim(timestamp)//trim(restart_file_name)
 call error_mesg('save_crop_restart', 'writing NetCDF restart "'//trim(filename)//'"', NOTE)
 call init_land_restart(restart, filename, vegn_tile_exists, tile_dim_length)
 call add_restart_axis(restart,'month', (/(float(mn),mn=1,12)/),.false.,longname='calendar month')
 call add_restart_axis(restart,'crop_cal',(/(float(mn),mn=1,12)/),.false.,longname='crop calendar')
 call add_tile_data(restart,'T_mid_mth', 'month', vegn_T_mid_mth_ptr, 'climatological average mid-month canopy air temperature','degK')
 call add_tile_data(restart,'P_mid_mth', 'month', vegn_P_mid_mth_ptr, 'climatological average mid-month precipitation rate','mm/sec')
 call add_tile_data(restart,'tc_av_climate', 'month', vegn_tc_av_climate_ptr, 'climatological monthly average canopy air temperature','degK')
 call add_tile_data(restart,'precip_av_climate', 'month', vegn_precip_av_climate_ptr,'climatological monthly average precipitation rate','mm/sec')
 call add_int_tile_data(restart,'current_crop', vegn_current_crop_ptr, '1=Maize 2=Soybean 3=Rice 4=Spring Wheat 5=Winter Wheat','dimensionless')
 call add_tile_data(restart,'crop_cal_Maize', 'crop_cal', vegn_crop_cal_Maize_ptr, 'rainfed maize crop calendar', 'day_of_year')
 call add_tile_data(restart,'crop_cal_Soy', 'crop_cal', vegn_crop_cal_Soy_ptr, 'rainfed soybean crop calendar', 'day_of_year')
 call add_tile_data(restart,'crop_cal_SW', 'crop_cal', vegn_crop_cal_SW_ptr, 'rainfed spring wheat crop calendar', 'day_of_year')
 call add_tile_data(restart,'crop_cal_WW', 'crop_cal', vegn_crop_cal_WW_ptr, 'rainfed winter wheat crop calendar', 'day_of_year')
 call add_tile_data(restart,'crop_cal_Rice_1', 'crop_cal', vegn_crop_cal_Rice_1_ptr, 'rainfed Rice crop calendar. Main crop.', 'day_of_year')
 call add_tile_data(restart,'crop_cal_Rice_2', 'crop_cal', vegn_crop_cal_Rice_2_ptr, 'rainfed Rice crop calendar. Second crop.', 'day_of_year')
 call add_int_tile_data(restart,'idle', vegn_idle_ptr, '-1 = .true.  0 = .false.','dimensionless')
 call save_land_restart(restart)
 call free_land_restart(restart)
 end subroutine save_crop_restart
!======================================================================================================================================================
 subroutine crop_end()
  crop_mod_initialized = .FALSE.
 end subroutine crop_end
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
subroutine vegn_current_crop_ptr(t,p)
 type(land_tile_type),pointer::t
 integer,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%current_crop
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_idle_ptr(t,p)
 type(land_tile_type),pointer::t
 integer,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%idle
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_crop_cal_Maize_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%crop_cal_Maize(n)
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_crop_cal_Soy_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%crop_cal_Soy(n)
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_crop_cal_SW_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%crop_cal_SW(n)
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_crop_cal_WW_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%crop_cal_WW(n)
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_crop_cal_Rice_1_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%crop_cal_Rice_1(n)
 endif
 end subroutine
!======================================================================================================================================================
subroutine vegn_crop_cal_Rice_2_ptr(t,n,p)
 type(land_tile_type),pointer::t
 integer,intent(in)::n
 real,pointer::p
 p=>NULL()
 if(associated(t))then
 if(associated(t%vegn))p=>t%vegn%Crop%crop_cal_Rice_2(n)
 endif
 end subroutine
!======================================================================================================================================================
 end module vegn_crop_mod
