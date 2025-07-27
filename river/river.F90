module river_mod

#include "../shared/debug.inc"
!-----------------------------------------------------------------------
!                   GNU General Public License
!
! This program is free software; you can redistribute it and/or modify it and
! are expected to follow the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2 of
! the License, or (at your option) any later version.
!
! MOM is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.
! or see:   http://www.gnu.org/licenses/gpl.html
!-----------------------------------------------------------------------
! <CONTACT EMAIL="Kirsten.Findell@@noaa.gov"> Kirsten Findell </CONTACT>
! <CONTACT EMAIL="Zhi.Liang@@noaa.gov"> Zhi Liang </CONTACT>
! <NAMELIST NAME="river_nml">
! <DATA NAME="layout" TYPE="integer, dimension(2)">
!  Processor domain layout for river model. If layout(1)*layout(2) is not equal
!  to mpp_npes, the river model layout will be assigned the layout of land model
!  passed through river_init.
!  </DATA>
! <DATA NAME="do_rivers" TYPE="logical">
!   set true to run river model ( default is true). If FALSE, rivers are
!   essentially turned off to save computing time
!  </DATA>
! <DATA NAME="dt_slow" TYPE="real">
!   slow time step for river model. dt_slow must be integer multiplier of dt_fast passed
!   from land model ( land model time step).
!  </DATA>
! <DATA NAME="diag_freq" TYPE="integer">
!   Number of slow time steps between sending out diagnositics data(default is 1). Please
!   note that diagnostic output frequency ( specified in diag_table ) must be divided by
!   diag_freq*dt_slow.
!  </DATA>
! </NAMELIST>

use mpp_mod,             only : CLOCK_SUBCOMPONENT, CLOCK_ROUTINE, &
   mpp_error, FATAL, WARNING, NOTE, stdout, stdlog, input_nml_file, &
   mpp_pe, mpp_chksum, mpp_max, &
   mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_domains_mod,     only : domain2d, mpp_get_compute_domain, mpp_get_global_domain, &
   mpp_get_data_domain, mpp_update_domains, mpp_get_ntile_count, mpp_get_tile_id, &
   domainUG, mpp_get_UG_compute_domain, mpp_pass_ug_to_sg, mpp_pass_sg_to_ug
use fms_mod,             only : check_nml_error, string, CLOCK_FLAG_DEFAULT, error_mesg
use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, register_axis, &
   register_restart_field, variable_exists, register_field, &
   read_restart, close_file, write_data, &
   get_global_io_domain_indices, FmsNetcdfFile_t, &
   get_variable_size, read_data, get_variable_num_dimensions, unlimited, &
   get_instance_filename
use diag_manager_mod,    only : diag_axis_init, register_diag_field, &
   register_static_field, send_data, diag_field_add_attribute
use time_manager_mod,    only : time_type, increment_time, get_time
use data_override_mod,   only : data_override
use tracer_manager_mod, only : NO_TRACER

use river_type_mod,      only : river_type, Leo_Mad_trios, NO_RIVER_FLAG
use river_tracers_mod,   only : num_phys, num_species, trdata, river_tracer_index
use river_physics_mod,   only : river_physics_step, river_physics_init, &
   river_impedes_lake, river_impedes_large_lake
use constants_mod,       only : PI, RADIAN, tfreeze, DENS_H2O, hlf
use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT
use land_io_mod,         only : register_variable_string_attribute, read_field
use land_tile_mod,       only : land_tile_map, land_tile_type, land_tile_enum_type, &
   first_elmt, loop_over_tiles, nitems, elmt_at_index
use land_data_mod,       only : land_data_type, log_version, lnd
use land_debug_mod, only : is_watch_point, is_watch_cell, get_current_point, &
   set_current_point, check_var_range, check_conservation, land_error_message, &
   carbon_cons_tol, nitrogen_cons_tol
use lake_tile_mod,       only : num_l
use soil_tile_mod,      only : soil_tile_type, num_soil=>num_l, dz_soil=>dz
use lake_mod,           only : use_reservoir
use land_numerics_mod,  only : rank_descending

implicit none
private

!--- version information ---------------------------------------------
character(len=*), parameter :: module_name = 'river_mod'
#include "../shared/version_variable.inc"

!--- public interface ------------------------------------------------
public :: river_init, river_end, river_type, update_river, river_stock_pe
public :: save_river_restart
public :: get_river_water

!--- namelist interface ----------------------------------------------
logical            :: do_rivers       = .TRUE.  ! if FALSE, rivers are essentially turned off to save computing time
real               :: dt_slow
integer            :: diag_freq       = 1       ! Number of slow time steps between sending out diagnostics data.
logical            :: debug_river     = .FALSE.
real               :: Somin           = 0.00005 ! There are 7 points with So = -9.999 but basinid > 0....
real               :: outflowmean_min = 1.      ! temporary fix, should not allow zero in input file
logical            :: land_area_called_cellarea = .false.
logical            :: all_big_outlet_ctn0 = .false.

real, dimension(3) :: ave_DHG_exp = (/0.49,0.33,0.18/)  ! (/B, F, M for avg of many rivers, 15Nov05/)
real, dimension(3) :: ave_AAS_exp = (/0.19,0.39,0.42/)  ! (/b, f, m for avg of many rivers, 15Nov05/)
real, dimension(3) :: ave_DHG_coef = (/4.62,0.26,0.82/) ! (/A, C, K for avg of many rivers, 15Nov05/)
real               :: sinuosity = 1.3
real               :: channel_tau = 86400*365.25*10     ! channel geometry reflects average flow over O(10 y)
logical :: lake_area_bug = .FALSE. ! if set to true, reverts to buggy (quebec)
    ! behavior, where by mistake cell area was used instead of land area to
    ! compute the area of lakes.
logical :: stop_on_mask_mismatch = .TRUE. ! If set to false, then the data mismatches
    ! (mismatch of land and river masks, and discharges in points where there is no
    ! ocean) are reported, but do not cause the abort of the program.

! ZMS
logical :: tracers_from_runoff = .false. ! if true, use runoff_c(:,:,num_phys+1:num_species)
        ! rather than source concentration and flux files
logical :: do_groundwater_abstraction = .false.
logical :: do_deep_gw_abst = .false. ! If true, water is borrowed from imaginary
        ! "deep aquifers" of infinite capacity, violating water conservation
        ! in the system.

namelist /river_nml/ dt_slow, diag_freq, debug_river,                      &
                     Somin, outflowmean_min, ave_DHG_exp, ave_AAS_exp,     &
                     ave_DHG_coef, do_rivers, sinuosity, channel_tau,      &
                     land_area_called_cellarea, all_big_outlet_ctn0,       &
                     lake_area_bug, stop_on_mask_mismatch,                 &
                     tracers_from_runoff, &
                     do_groundwater_abstraction, do_deep_gw_abst

character(len=128) :: river_src_file   = 'INPUT/river_data.nc'
character(len=128) :: river_Omean_file = 'INPUT/river_Omean.nc'
character(len=128) :: river_threshold_file = 'INPUT/threshold.nc'
character(len=128) :: env_flow_file = 'INPUT/env_flow.nc'

!---------------------------------------------------------------------
logical :: module_is_initialized = .FALSE.
integer :: isc, iec, jsc, jec                         ! compute domain decomposition
integer :: isd, ied, jsd, jed                         ! data domain decomposition
integer :: lsc, lec                                   ! unstructured domain decomposition
integer :: nlon, nlat                                 ! size of computational river grid
integer :: num_lake_lev
integer :: id_outflowmean, id_lake_depth_sill
integer :: id_dx, id_basin, id_So, id_depth, id_width, id_vel
integer :: id_lake_abst, id_lake_habst
integer :: id_rsv_outflow
integer :: id_irr_full, id_irr_met, id_irr_unmet
integer :: id_gw_s_abst, id_gw_d_abst, id_gw_s_habst, id_gw_d_habst
integer :: id_LWSr, id_FWSr, id_HSr, id_meltr
integer :: id_travel, id_elev, id_tocell
integer :: maxtravel
real    :: missing = -1.e8

real,    parameter :: CONST_OMEAN = 80000
real,    parameter :: epsln = 1.e-6
real,    parameter :: sec_in_day = 86400.

real     :: discharge_tol=0.0, clw=0.0, csw=0.0  ! will get these values from land model
integer  :: i_river_ice, i_river_heat, i_river_DOC
logical, allocatable, dimension(:,:) :: missing_rivers
real,  allocatable, dimension(:,:)   :: discharge2ocean_next   ! store discharge value
real,  allocatable, dimension(:,:,:) :: discharge2ocean_next_c ! store discharge value
! IDs of diag fields normalized per land area
integer, allocatable, dimension(:)   :: id_infloc,  id_storage, id_stordis, id_inflow, &
      id_run_stor, id_outflow, id_removal, id_dis, id_lake_outflow, id_abstflow
! IDs of diag fields normalized per cell area
integer, allocatable, dimension(:)   :: id_infloc_c,  id_storage_c, id_stordis_c, id_inflow_c, &
      id_run_stor_c, id_outflow_c, id_removal_c, id_dis_c, id_lake_outflow_c, id_abstflow_c
integer :: id_dis_liq,  id_dis_ice,  id_dis_heat, id_dis_sink, id_dis_DOC, id_no_riv
integer, public :: num_fast_calls !public for soil_mod, this is not good, but in original code of lm4p1, soil_mod calls river_mod
integer :: slow_step = 0          ! record number of slow time step run.
type(domain2d), pointer :: domain    => NULL()
type(domainUG), pointer :: UG_domain => NULL()
type(river_type), save :: River

!--- clock id variable
integer :: slowclock, bndslowclock, physicsclock, diagclock, riverclock

character(len=8),parameter :: river_res_xdim = "xaxis_1"
character(len=8),parameter :: river_res_ydim = "yaxis_1"
character(len=8),parameter :: river_res_zdim = "zaxis_1"

contains ! ===--------------------------------------------------------

!#####################################################################
  subroutine river_init( land_lon, land_lat, time, dt_fast, land_domain, land_UG_domain, &
                         land_frac, discharge_tol_in, clw_in, csw_in )
    real,            intent(in) :: land_lon(:,:)     ! geographical longitude of cell center
    real,            intent(in) :: land_lat(:,:)     ! geographical latitude of cell center
    type(time_type), intent(in) :: time              ! current time
    type(time_type), intent(in) :: dt_fast           ! fast time step
    type(domain2d),  intent(in), target :: land_domain       ! land domain
    type(domainUG), intent(in), target :: land_UG_domain
    real,            intent(in) :: land_frac(:,:)       ! land area fraction from land model on UG_domain
    real,            intent(in) :: discharge_tol_in, clw_in, csw_in

    integer              :: unit, io_status, ierr, id_restart
    integer              :: sec, day, i, j, i_species
    integer              :: nxc, nyc
    character(len=*), parameter :: filename = "INPUT/river.nc"
    integer              :: id_lon, id_lat, id_lonb, id_latb
    type(Leo_Mad_trios)   :: DHG_exp            ! downstream equation exponents
    type(Leo_Mad_trios)   :: DHG_coef           ! downstream equation coefficients
    type(Leo_Mad_trios)   :: AAS_exp            ! at-a-station equation exponents

    type(FmsNetcdfDomainFile_t) :: river_restart
    logical :: exists
    type(FmsNetcdfFile_t) :: river_input

    riverclock = mpp_clock_id('update_river'           , CLOCK_FLAG_DEFAULT, CLOCK_SUBCOMPONENT)
    slowclock = mpp_clock_id('update_river_slow'       , CLOCK_FLAG_DEFAULT, CLOCK_ROUTINE)
    bndslowclock = mpp_clock_id('update_river_bnd_slow', CLOCK_FLAG_DEFAULT, CLOCK_ROUTINE)
    physicsclock = mpp_clock_id('river phys'           , CLOCK_FLAG_DEFAULT, CLOCK_ROUTINE)
    diagclock    = mpp_clock_id('river diag'           , CLOCK_FLAG_DEFAULT, CLOCK_ROUTINE)

!--- read namelist -------------------------------------------------
    read (input_nml_file, nml=river_nml, iostat=io_status)
    ierr = check_nml_error(io_status, 'river_nml')

!--- write version and namelist info to logfile --------------------
    call log_version(version, module_name, __FILE__)
    unit=stdlog()
    write(unit, river_nml)

    if(.not.do_rivers) return ! do nothing further if the rivers are turned off

!--- check name list variables
    if(diag_freq .le. 0) call mpp_error(FATAL,'river_mod: diag_freq should be a positive integer')

! set up time-related values
    River%time = time
    call get_time(dt_fast, sec, day)
    River%dt_fast = day*sec_in_day+sec

    River%dt_slow = dt_slow
    River%channel_tau = channel_tau

    num_fast_calls = River%dt_slow/River%dt_fast

    call mpp_error(NOTE,'river_mod: tracer numbers: num_phys='//string(num_phys)//' num_species='//string(num_species))

    if(River%dt_slow .lt. River%dt_fast) call mpp_error(FATAL, &
         'river_mod: river slow time step dt_slow should be no less than land model fast time step dt_fast')

    if ( mod(River%dt_slow,River%dt_fast) .ne. 0  ) call mpp_error(FATAL, &
         'river_mod: river slow time step dt_slow should be multiple of land model fast time step dt_fast')

    discharge_tol = discharge_tol_in
    clw = clw_in
    csw = csw_in

!--- get the domain decomposition, river and land will be on the same grid and have the same domain decomposition.
    domain => land_domain
    UG_domain => land_UG_domain
    call mpp_get_global_domain (domain, xsize=River%nlon, ysize=River%nlat)
    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain   (domain, isd, ied, jsd, jed)
    call mpp_get_UG_compute_domain(UG_domain, lsc, lec)

!---- make sure the halo size is 1
    if( ied-iec .NE. 1 .OR. isc-isd .NE. 1 .OR. jed-jec .NE. 1 .OR. jsc-jsd .NE. 1 ) &
      call mpp_error(FATAL, "river_mod: halo size in four direction should all be 1")

    nxc = iec - isc + 1; nyc = jec - jsc + 1
    !--- make sure land_lon, land_lat is on the compute domain
    if(size(land_lon,1) .NE. nxc .OR. size(land_lon,2) .NE. nyc ) call mpp_error(FATAL, &
        "river_mod: land_lon should be on the compute domain")
    if(size(land_lat,1) .NE. nxc .OR. size(land_lat,2) .NE. nyc ) call mpp_error(FATAL, &
        "river_mod: land_lat should be on the compute domain")

    allocate(missing_rivers        (isc:iec,jsc:jec            ))
    allocate(discharge2ocean_next  (isc:iec,jsc:jec            ))
    allocate(discharge2ocean_next_c(isc:iec,jsc:jec,num_species))
    allocate(id_infloc (0:num_species), id_storage(0:num_species))
    allocate(id_inflow (0:num_species), id_outflow(0:num_species))
    allocate(id_dis    (0:num_species), id_lake_outflow (0:num_species))
    allocate(id_removal(0:num_species), id_stordis(0:num_species), id_run_stor(0:num_species))
    allocate(id_abstflow(0:num_species))
    discharge2ocean_next = 0
    discharge2ocean_next_c = 0
    ! IDs of diag fields normalized per cell area
    allocate(id_infloc_c (0:num_species), id_storage_c(0:num_species))
    allocate(id_inflow_c (0:num_species), id_outflow_c(0:num_species))
    allocate(id_dis_c    (0:num_species), id_lake_outflow_c (0:num_species))
    allocate(id_removal_c(0:num_species), id_stordis_c(0:num_species), id_run_stor_c(0:num_species))
    allocate(id_abstflow_c(0:num_species))

!--- read the data from the file river_src_file -- has all static river network data
    call get_river_data(land_lon, land_lat, land_frac)

    missing_rivers = (lnd%sg_landfrac .gt. 0 .and. .not. River%mask)

    i_river_ice  = river_tracer_index('ice')
    i_river_heat = river_tracer_index('het')
    i_river_DOC  = river_tracer_index('doc')
    if (i_river_ice  == NO_TRACER) call mpp_error(FATAL, 'river_mod: required river tracer for ice not found')
    if (i_river_heat == NO_TRACER) call mpp_error(FATAL, 'river_mod: required river tracer for heat not found')

!--- register diag field
    if(mpp_get_ntile_count(domain)==1) then
       ! grid has just one tile, so we assume that the grid is regular lat-lon
       ! define longitude axes and its edges
       id_lonb = diag_axis_init ( &
            'lonb', lnd%coord_glonb, 'degrees_E', 'X', 'longitude edges', &
            set_name='river', domain2=domain )
       id_lon  = diag_axis_init (                                                &
            'lon',  lnd%coord_glon, 'degrees_E', 'X',  &
            'longitude', set_name='river',  edges=id_lonb, domain2=domain )

       ! define latitude axes and its edges
       id_latb = diag_axis_init ( &
            'latb', lnd%coord_glatb, 'degrees_N', 'Y', 'latitude edges',  &
            set_name='river',  domain2=domain   )
       id_lat = diag_axis_init (                                                &
            'lat',  lnd%coord_glat, 'degrees_N', 'Y', &
            'latitude', set_name='river', edges=id_latb, domain2=domain   )
    else
       id_lon = diag_axis_init ( 'grid_xt', (/(real(i),i=1,River%nlon)/), 'degrees_E', 'X', &
            'T-cell longitude', set_name='river',  domain2=domain, aux='geolon_t' )
       id_lat = diag_axis_init ( 'grid_yt', (/(real(i),i=1,River%nlat)/), 'degrees_N', 'Y', &
            'T-cell latitude', set_name='river',  domain2=domain, aux='geolat_t' )
    endif

    call river_diag_init (id_lon, id_lat)

!--- read restart file
    exists = open_file(river_restart, filename, "read", domain, &
                       is_restart=.true.)
    if (exists) then
        call mpp_error(NOTE, 'river_init : Read restart files '//trim(filename))
        call register_axis(river_restart, river_res_xdim, "x")
        call register_axis(river_restart, river_res_ydim, "y")
        call register_restart_field(river_restart, "storage", river%storage)
        call register_restart_field(river_restart, "discharge2ocean", &
                                    discharge2ocean_next)
        if (variable_exists(river_restart, "discharge2ocean_c")) then
            call register_restart_field(river_restart, "discharge2ocean_c", discharge2ocean_next_c)
            call register_restart_field(river_restart, "storage_c", river%storage_c)
        else
            do i_species = 1, num_species
                if (variable_exists(river_restart, "disch2ocn_"//trdata(i_species)%name)) then
                    call register_restart_field(river_restart, "disch2ocn_"//trdata(i_species)%name, &
                                                discharge2ocean_next_c(:,:,i_species))
                else
                    call mpp_error(NOTE, "river_init: disch2ocn_"//trim(trdata(i_species)%name)//" does not exist in "//trim(filename))
                endif
                if (variable_exists(river_restart, "storage_"//trdata(i_species)%name)) then
                    call register_restart_field(river_restart, "storage_"//trdata(i_species)%name, &
                                                river%storage_c(:,:,i_species))
                else
                    call mpp_error(NOTE, "river_init: storage_"//trim(trdata(i_species)%name)//" does not exist in "//trim(filename))
                endif
              enddo
           endif
        call register_restart_field(river_restart, "Omean", river%outflowmean)
        if (variable_exists(river_restart, "depth")) then
            call register_restart_field(river_restart, "depth", river%depth)
        endif
        call read_restart(river_restart)
        call close_file(river_restart)
    else
        call mpp_error(NOTE, 'river_init : cold start, set data to 0')
        River%storage    = 0.0
        River%storage_c  = 0.0
        discharge2ocean_next   = 0.0
        discharge2ocean_next_c = 0.0
        exists = open_file(river_input, river_Omean_file, "read")
        if (exists) then
           call read_data(river_restart, 'Omean', River%outflowmean)
           call close_file(river_input)
        else
           River%outflowmean = CONST_OMEAN
        end if
    endif
    River%stordis_c = River%dt_slow * discharge2ocean_next_c/DENS_H2O
    River%stordis   = River%dt_slow *(discharge2ocean_next + &
                                      discharge2ocean_next_c(:,:,1))/DENS_H2O
    where(River%outflowmean .le. outflowmean_min) River%outflowmean=outflowmean_min

    maxtravel = maxval(River%travel)
    call mpp_max(maxtravel)

    call river_physics_init(River, domain, id_lon, id_lat)
    call get_Leo_Mad_params(DHG_exp, DHG_coef, AAS_exp)
    River%o_exp  = 1./ (AAS_exp%on_w + AAS_exp%on_d)
    do j = jsc, jec
       do i = isc, iec
          if ( River%reach_length(i,j) > 0.0) then
              River%o_coef(i,j) = River%outflowmean(i,j) / &
                   ((sinuosity*River%reach_length(i,j))*DHG_coef%on_w*DHG_coef%on_d &
                   *(River%outflowmean(i,j)**(DHG_exp%on_w+DHG_exp%on_d)))**River%o_exp
          endif
       enddo
    enddo
    River%d_exp  = AAS_exp%on_d
    River%d_coef = DHG_coef%on_d                        &
         *(River%outflowmean**(DHG_exp%on_d-AAS_exp%on_d))
    River%w_exp  = AAS_exp%on_w
    River%w_coef = DHG_coef%on_w                        &
         *(River%outflowmean**(DHG_exp%on_w-AAS_exp%on_w))

    num_lake_lev = num_l
    module_is_initialized = .TRUE.

  end subroutine river_init

!#####################################################################
  subroutine update_river ( runoff, runoff_c, land2cplr)
    real, dimension(:,:),   intent(in)  :: runoff
    real, dimension(:,:,:), intent(in)  :: runoff_c
    type(land_data_type), intent(inout) :: land2cplr

    ! --- local vars
    real, dimension(size(runoff,1),size(runoff,2)) :: &
        heat_frac_liq,    & ! fraction of runoff heat in liquid
        discharge_l,      & ! discharge of liquid water to ocean
        discharge_sink      ! container to collect small/negative values for later accounting
    real, dimension(size(runoff,1),size(runoff,2),num_species) ::  &
        discharge_c    ! runoff of tracers accumulated over tiles in cell (including ice and heat)

    integer, save :: n = 0  ! fast time step with each slow time step
    integer       :: i_species
    logical       :: used

    call mpp_clock_begin(riverclock)
    if (.not.do_rivers) then
        call mpp_clock_end(riverclock)  !needed for early exit when do_rivers=.false.
        return
    endif

    discharge_l   = discharge2ocean_next
    discharge_c(:,:,1:num_species) = discharge2ocean_next_c(:,:,1:num_species)
! deplete the discharge storage pools
    River%stordis_c = River%stordis_c &
         - River%dt_fast * discharge_c/DENS_H2O
    River%stordis   = River%stordis   &
         - River%dt_fast *(discharge_l + &
                           discharge_c(:,:,1))/DENS_H2O

!  increment time
    River%Time = increment_time(River%Time, River%dt_fast, 0)
    n = n + 1
!--- accumulate runoff ---------------------
    River%run_stor   = River%run_stor   + runoff
    River%run_stor_c = River%run_stor_c + runoff_c

    if(n == num_fast_calls) then
        call mpp_clock_begin(slowclock)
        call update_river_slow(River%run_stor/real(num_fast_calls), &
             River%run_stor_c(:,:,:)/real(num_fast_calls))
        call mpp_clock_end(slowclock)
        call mpp_clock_begin(bndslowclock)
        call update_river_bnd_slow
        call mpp_clock_end(bndslowclock)
        n = 0
        River%run_stor = 0
        River%run_stor_c = 0
    endif

    discharge_l = discharge_l/lnd%sg_cellarea
    do i_species = 1, num_species
       discharge_c(:,:,i_species) =  discharge_c(:,:,i_species)/lnd%sg_cellarea
    enddo

    ! pass through to ocean the runoff that was not seen by river module because of land_frac diffs.
    ! need to multiply by gfrac to spread over whole cell
    where (missing_rivers) discharge_l = (runoff-runoff_c(:,:,i_river_ice))*lnd%sg_landfrac
    do i_species = 1, num_species
       where (missing_rivers) &
          discharge_c(:,:,i_species) = runoff_c(:,:,i_species)*lnd%sg_landfrac
    enddo

    ! do not send negatives or insignificant values to ocean. put them in the sink instead.
    ! this code does not seem necessary, and default discharge_tol value should be used.
    discharge_sink = 0.0
    where (discharge_l.le.discharge_tol)
       discharge_sink = discharge_sink + discharge_l
       discharge_l    = 0.0
    end where
    where (discharge_c(:,:,i_river_ice).le.discharge_tol)
       discharge_sink     = discharge_sink + discharge_c(:,:,i_river_ice)
       discharge_c(:,:,i_river_ice) = 0.0
    end where

    ! find phase partitioning ratio for discharge sensible heat flux
    where (discharge_l.gt.0. .or. discharge_c(:,:,i_river_ice).gt.0.)
       heat_frac_liq = clw*discharge_l / (clw*discharge_l+csw*discharge_c(:,:,i_river_ice))
    elsewhere
       heat_frac_liq = 1.0
    end where

    ! scale up fluxes sent to ocean to compensate for non-ocean fraction of discharge cell.
    ! split heat into liquid and solid streams
    where (lnd%sg_landfrac.lt.1.)
       land2cplr%discharge           = discharge_l        / (1-lnd%sg_landfrac)
       land2cplr%discharge_snow      = discharge_c(:,:,i_river_ice) / (1-lnd%sg_landfrac)
       land2cplr%discharge_heat      = heat_frac_liq*discharge_c(:,:,i_river_heat) / (1-lnd%sg_landfrac)
       land2cplr%discharge_snow_heat =               discharge_c(:,:,i_river_heat) / (1-lnd%sg_landfrac) &
                                      - land2cplr%discharge_heat
    end where

#ifdef ZMSDEBUG
    do j=lnd%js,lnd%je
    do i=lnd%is,lnd%ie
       call set_current_point(i, j, 1)
       call check_var_range(land2cplr%discharge(i,j), -10., 10., 'gcell DISCHARGE CHECK', &
                            'Liquid Discharge (mm/s)', WARNING)
       call check_var_range(land2cplr%discharge_snow(i,j), -10., 10., 'gcell DISCHARGE CHECK', &
                            'Snow Discharge (mm/s)', WARNING)
       call check_var_range(land2cplr%discharge_heat(i,j), -1000., 1000., 'gcell DISCHARGE CHECK', &
                            'Discharge Heat (W/m^2/s)', WARNING)
       call check_var_range(land2cplr%discharge_snow_heat(i,j), -1000., 1000., 'gcell DISCHARGE CHECK', &
                            'Discharge Snow Heat (W/m^2/s)', WARNING)
    end do
    end do
#endif

    if (id_dis_liq > 0)  used = send_data (id_dis_liq,  discharge_l, lnd%time)
    if (id_dis_ice > 0)  used = send_data (id_dis_ice,  discharge_c(:,:,i_river_ice), lnd%time)
    if (id_dis_heat > 0) used = send_data (id_dis_heat, discharge_c(:,:,i_river_heat), lnd%time)
    if (id_dis_sink > 0) used = send_data (id_dis_sink, discharge_sink, lnd%time)
    if (id_dis_DOC > 0.and.i_river_DOC/=NO_TRACER) &
                         used = send_data (id_dis_DOC,  discharge_c(:,:,i_river_DOC), lnd%time)

    call mpp_clock_end(riverclock)

  end subroutine update_river

!#####################################################################
  subroutine update_river_slow(runoff, runoff_c)
    real, dimension(:,:),   intent(in)  :: runoff
    real, dimension(:,:,:), intent(in)  :: runoff_c

    real, dimension(isd:ied,jsd:jed) :: &
                             lake_sfc_A, lake_sfc_bot, lake_conn, &
                             Afrac_rsv, Vfrac_rsv
    real, dimension(isd:ied,jsd:jed,num_lake_lev) :: &
                             lake_wl, lake_ws, lake_dz, lake_dhcap
    real, dimension(isc:iec,jsc:jec) :: &
                             lake_depth_sill, lake_width_sill, lake_backwater, &
                             lake_backwater_1, &
                             lake_whole_area, &
                             rivr_LMASS,       & ! mass of liquid water in rivers in cell
                             rivr_FMASS,       & ! mass of ice in rivers in cell
                             rivr_MELT,        & ! net mass melt in rivers in cell
                             rivr_HEAT,        & ! sensible heat content of rivers in cell
                             irr_demand,       & ! left irrigation demand
                             rsv_depth,        &
                             rsv_outflow
    real, dimension(isc:iec,jsc:jec,num_lake_lev) :: &
                             lake_T

    real, dimension(lnd%ls:lnd%le) :: &
                             lake_sfc_A_ug, lake_sfc_bot_ug, lake_conn_ug
    real, dimension(lnd%ls:lnd%le,num_lake_lev) :: &
                             lake_wl_ug, lake_ws_ug, lake_dz_ug, lake_dhcap_ug
    real, dimension(lnd%ls:lnd%le) :: &
                             lake_depth_sill_ug, lake_width_sill_ug, lake_backwater_ug, &
                             lake_backwater_1_ug, &
                             lake_whole_area_ug, &
                             irr_demand_ug, lake_abst_ug, lake_habst_ug, river_abst_ug, &
                             rsv_depth_ug, Afrac_rsv_ug, Vfrac_rsv_ug
    real, dimension(lnd%ls:lnd%le,num_lake_lev) :: &
                             lake_T_ug
    real, dimension(lnd%ls:lnd%le,num_species) :: river_abstflow_c_ug

    integer                             :: travelnow, lev, l
    type(Leo_Mad_trios)   :: DHG_exp
    type(Leo_Mad_trios)   :: DHG_coef
    type(Leo_Mad_trios)   :: AAS_exp
    integer i,j,k, i_next, j_next, i_species
    type(land_tile_enum_type)     :: ce    ! land tile enumerator
    type(land_tile_type), pointer :: tile  ! pointer to current tile
    type(soil_tile_type), pointer :: soil
    logical :: used
    integer :: ntiles, nlow
    real,    allocatable :: priority(:) ! priority of the gw withdrawal for each tile
    integer, allocatable :: idx(:)      ! array of tile indices in the descending priority order
    real :: abst_thres = 1.e-15 !m3
    real,    allocatable :: hlsp_irr_demand_gw(:)
    real :: irr_demand_gw
    integer :: nk_g, hidxk

    ! variables for data override
    real, dimension(isc:iec,jsc:jec) :: src_conc, src_flux
    logical :: src_flux_overridden, src_conc_overridden
    real, dimension(lnd%ls:lnd%le) :: tot_demand_full !kg
    real :: tile_demand_full !kg
    real :: demand_left_tile !kg/m2
    real, dimension(lnd%ls:lnd%le) :: demand_full_ug, demand_met_ug, demand_unmet_ug !m3
    real, dimension(lnd%ls:lnd%le) :: gw_s_abst_ug, gw_d_abst_ug !m3
    real, dimension(lnd%ls:lnd%le) :: gw_s_habst_ug, gw_d_habst_ug !J
    real, dimension(isc:iec,jsc:jec) :: demand_full, demand_met, demand_unmet !m3
    real, dimension(isc:iec,jsc:jec) :: gw_s_abst, gw_d_abst !m3
    real, dimension(isc:iec,jsc:jec) :: gw_s_habst, gw_d_habst !J
    real :: tot_abst, tot_habst, shallow_abst, shallow_habst, deep_abst, deep_habst, frac

    slow_step = slow_step + 1

    River%infloc   = River%land_area*runoff  /DENS_H2O !m2 * kg/(m2 s) / kg/m3 = m3/s
    River%infloc_c = 0
    do i_species = 1, num_species
       River%infloc_c(:,:,i_species) = River%land_area*runoff_c(:,:,i_species)/DENS_H2O
       src_conc = 0.0
       src_flux = 0.0
       call data_override('LND','river_src_flux_'//trdata(i_species)%name, src_flux, River%time, override=src_flux_overridden)
       call data_override('LND','river_src_conc_'//trdata(i_species)%name, src_conc, River%time, override=src_conc_overridden)
       if (src_conc_overridden.OR.src_flux_overridden) then
          where (River%land_area.gt.0)  &
               River%infloc_c(:,:,i_species) = River%infloc*src_conc + src_flux
       endif
    enddo

    River%inflow   = 0
    River%inflow_c = 0
    River%lake_outflow   = 0
    River%lake_outflow_c = 0
    River%disw2o = 0.
    River%disc2o = 0.
    River%melt   = 0.
    lake_sfc_A  = 0
    lake_sfc_bot= 0
    lake_T  = 0
    lake_wl = 0
    lake_ws = 0
    lake_dz = 0
    lake_depth_sill  = 0
    lake_width_sill  = 0
    lake_whole_area  = 0
    lake_conn   = 0
    lake_backwater = 0
    lake_backwater_1 = 0
    irr_demand = 0
    rsv_depth = 0
    Afrac_rsv = 0
    Vfrac_rsv = 0
    rsv_outflow = 0
    lake_sfc_A_ug  = 0
    lake_sfc_bot_ug= 0
    lake_T_ug  = 0
    lake_wl_ug = 0
    lake_ws_ug = 0
    lake_dz_ug = 0
    lake_dhcap_ug = 0
    lake_depth_sill_ug  = 0
    lake_width_sill_ug  = 0
    lake_whole_area_ug  = 0
    lake_conn_ug   = 0
    lake_backwater_ug = 0
    lake_backwater_1_ug = 0
    irr_demand_ug = 0
    lake_abst_ug = 0
    lake_habst_ug = 0
    river_abst_ug = 0
    river_abstflow_c_ug = 0
    rsv_depth_ug = 0
    Afrac_rsv_ug = 0
    Vfrac_rsv_ug = 0

    ce = first_elmt(land_tile_map, ls=lnd%ls)
    do while(loop_over_tiles(ce, tile, l,k))
       if (.not.associated(tile%lake)) cycle
       if (lake_area_bug) then
          lake_sfc_A_ug (l) = tile%frac * lnd%ug_cellarea(l)
       else
          lake_sfc_A_ug (l) = tile%frac * lnd%ug_area(l)
       endif
       do lev = 1, num_lake_lev
         lake_T_ug (l,lev)   = tile%lake%T(lev)
         lake_wl_ug(l,lev)   = tile%lake%wl(lev)
         lake_ws_ug(l,lev)   = tile%lake%ws(lev)
         lake_dz_ug(l,lev)   = tile%lake%dz(lev)
         lake_dhcap_ug(l,lev)= tile%lake%heat_capacity_dry(lev)
       enddo
       if(use_reservoir)then
         rsv_depth_ug(l)       = tile%lake%rsv_depth
         Afrac_rsv_ug(l)       = tile%lake%Afrac_rsv
         Vfrac_rsv_ug(l)       = tile%lake%Vfrac_rsv
         !this is still an approximation, because we didn't consider reservoir area in other gridcells with the same lake
         !if((.not.do_lake_change))then
         !  lake_whole_area_ug(l) = max(0., tile%lake%pars%whole_area-Afrac_rsv_ug(l)*tile%frac*lnd%ug_area(l))
         !else
         !  lake_whole_area_ug(l) = tile%lake%pars%whole_area
         !endif
         if(Afrac_rsv_ug(l)<1.)then
           lake_sfc_bot_ug(l) = (1.-Vfrac_rsv_ug(l))*lake_sfc_A_ug(l)*(sum(tile%lake%wl(:)+tile%lake%ws(:))-tile%lake%wl(1)-tile%lake%ws(1))/DENS_H2O & !m2 * kg/m2 / (kg/m3) = m3
                               /((1.-Afrac_rsv_ug(l))*lake_sfc_A_ug(l)) !m2
         else
           lake_sfc_bot_ug(l) = 0.
         endif
       else
         rsv_depth_ug(l)       = 0.
         Afrac_rsv_ug(l)       = 0.
         Vfrac_rsv_ug(l)       = 0.
         !lake_whole_area_ug(l) = tile%lake%pars%whole_area !+ Afrac_rsv_ug(l)*tile%frac*lnd%ug_area(l)
         lake_sfc_bot_ug(l)    = (sum(tile%lake%wl(:)+tile%lake%ws(:)) &
                                 -tile%lake%wl(1)-tile%lake%ws(1) ) &
                                      / DENS_H2O
       endif
       lake_whole_area_ug(l)  = tile%lake%pars%whole_area
       lake_depth_sill_ug(l)  = tile%lake%pars%depth_sill
       lake_width_sill_ug(l)  = tile%lake%pars%width_sill
       lake_conn_ug (l)       = tile%lake%pars%connected_to_next
       lake_backwater_ug(l)   = tile%lake%pars%backwater
       lake_backwater_1_ug(l) = tile%lake%pars%backwater_1
    enddo
    ! get irrigation demand
    tot_demand_full(:) = 0. !kg
    do l = lnd%ls,lnd%le
      ce = first_elmt(land_tile_map(l))
      do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%soil)) cycle
        tot_demand_full(l) = tot_demand_full(l) + tile%frac*lnd%ug_area(l) * tile%soil%irr_demand_ac !m2 * kg/m2 = kg
      enddo
    enddo

!z1l: The following might be changed for performance issue. This might be a temporary solution.
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_sfc_A_ug, lake_sfc_A)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_T_ug, lake_T)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_wl_ug, lake_wl)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_ws_ug, lake_ws)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_dz_ug, lake_dz)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_dhcap_ug, lake_dhcap)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_sfc_bot_ug, lake_sfc_bot)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_depth_sill_ug, lake_depth_sill)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_width_sill_ug, lake_width_sill)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_whole_area_ug, lake_whole_area)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_conn_ug, lake_conn)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_backwater_ug, lake_backwater)
    call mpp_pass_UG_to_SG(lnd%ug_domain, lake_backwater_1_ug, lake_backwater_1)
    call mpp_pass_UG_to_SG(lnd%ug_domain, tot_demand_full, irr_demand) !kg
    irr_demand = irr_demand/DENS_H2O !kg / kg/m3 = m3
    call mpp_pass_UG_to_SG(lnd%ug_domain, rsv_depth_ug, rsv_depth)
    call mpp_pass_UG_to_SG(lnd%ug_domain, Afrac_rsv_ug, Afrac_rsv)
    call mpp_pass_UG_to_SG(lnd%ug_domain, Vfrac_rsv_ug, Vfrac_rsv)
    call mpp_update_domains (lake_sfc_A,  domain)
    call mpp_update_domains (lake_sfc_bot,domain)
    call mpp_update_domains (lake_wl, domain)
    call mpp_update_domains (lake_ws, domain)
    call mpp_update_domains (lake_dz, domain)
    call mpp_update_domains (lake_dhcap,  domain)
    call mpp_update_domains (lake_conn,   domain)
    call mpp_update_domains (Afrac_rsv,   domain)
    call mpp_update_domains (Vfrac_rsv,   domain)
    do i=isc,iec
       do j=jsc,jec
          if (River%i_tocell(i,j)/=NO_RIVER_FLAG) then
             i_next = River%i_tocell(i,j)
             j_next = River%j_tocell(i,j)
          else
             ! to avoid indices out of bounds in the lake_sfc_A check
             i_next = i; j_next=j
          endif

          if (lake_backwater(i,j).gt.0.5 .and. lake_sfc_A(i,j).gt.0. .and. &
                                            lake_sfc_A(i_next,j_next).gt.0. ) then
             ! because of river backwater, lake in this cell relaxes toward level of
             ! lake in next cell downstream. (river depth is still simple function
             ! of discharge though.)
             if(use_reservoir.and.Afrac_rsv(i_next,j_next)>0.)then
               lake_depth_sill(i,j) = lake_sfc_bot(i_next,j_next) &
                +(lake_wl(i_next,j_next,1)+lake_ws(i_next,j_next,1))/DENS_H2O*(Vfrac_rsv(i_next,j_next)/Afrac_rsv(i_next,j_next))
             else
               lake_depth_sill(i,j) = lake_sfc_bot(i_next,j_next) &
                +(lake_wl(i_next,j_next,1)+lake_ws(i_next,j_next,1))/DENS_H2O
             endif
          elseif (lake_backwater_1(i,j).gt.0.5) then
             ! to determine depth of backwater, lake at coastal cell has base level
             ! set to river depth in same cell
             lake_depth_sill(i,j) = lake_depth_sill(i,j) + River%depth(i,j)
          elseif (lake_conn(i,j).gt.0.5 ) then
             ! for all but furthest downstream cell of a multi-cell lake,
             ! relax toward level in next cell (same lake) downstream
             if (lake_conn(i_next,j_next).gt.0.5 .or. all_big_outlet_ctn0) then
                  lake_depth_sill(i,j) = lake_sfc_bot(i_next,j_next) &
                   +(lake_wl(i_next,j_next,1)+lake_ws(i_next,j_next,1))/DENS_H2O
             endif
          elseif (river_impedes_lake) then
              if (lake_width_sill(i,j).lt.0..or.river_impedes_large_lake) then
                  ! lake level in cell relaxes toward river level in cell
                  lake_depth_sill(i,j) = lake_depth_sill(i,j) + River%depth(i,j)
              endif
          endif
       enddo
    enddo

! leftovers from horizontal mixing option, now gone
!call mpp_update_domains (lake_T,  domain)
!call mpp_update_domains (lake_depth_sill, domain)
!call mpp_update_domains (lake_tau, domain)

    travelnow = maxtravel
    do travelnow = maxtravel, 0, -1
       call mpp_clock_begin(physicsclock)
!***************************************************************
       call river_physics_step (River, travelnow, &
         lake_sfc_A, lake_sfc_bot, lake_depth_sill, &
         lake_width_sill, lake_whole_area,         &
         lake_T, lake_wl, lake_ws, lake_dz, lake_dhcap, irr_demand, &
         rsv_depth, Afrac_rsv, Vfrac_rsv, rsv_outflow )
!***************************************************************
       call mpp_clock_end(physicsclock)
    enddo

    call mpp_pass_SG_to_UG(lnd%ug_domain, lake_T, lake_T_ug)
    call mpp_pass_SG_to_UG(lnd%ug_domain, lake_wl, lake_wl_ug)
    call mpp_pass_SG_to_UG(lnd%ug_domain, lake_ws, lake_ws_ug)
    call mpp_pass_SG_to_UG(lnd%ug_domain, lake_dz, lake_dz_ug)
    call mpp_pass_SG_to_UG(lnd%ug_domain, lake_dhcap, lake_dhcap_ug)
    call mpp_pass_SG_to_UG(lnd%ug_domain, irr_demand, irr_demand_ug) !m3
    call mpp_pass_SG_to_UG(lnd%ug_domain, River%lake_abst, lake_abst_ug) !m3
    call mpp_pass_SG_to_UG(lnd%ug_domain, River%lake_habst, lake_habst_ug) !J
    call mpp_pass_SG_to_UG(lnd%ug_domain, River%abst, river_abst_ug) !m3
    call mpp_pass_SG_to_UG(lnd%ug_domain, River%abstflow_c, river_abstflow_c_ug)  ! m3/s, J m3/kg / s
    call mpp_pass_SG_to_UG(lnd%ug_domain, Vfrac_rsv, Vfrac_rsv_ug)

    ce = first_elmt(land_tile_map, ls=lnd%ls)
    do while(loop_over_tiles(ce, tile, l,k))
       if (.not.associated(tile%lake)) cycle
       do lev = 1, num_lake_lev
         tile%lake%T(lev)  = lake_T_ug (l,lev)
         tile%lake%wl(lev) = lake_wl_ug(l,lev)
         tile%lake%ws(lev) = lake_ws_ug(l,lev)
         tile%lake%dz(lev) = lake_dz_ug(l,lev)
       enddo
       if(use_reservoir) tile%lake%Vfrac_rsv = Vfrac_rsv_ug(l)
    enddo

    ! account for groundwater abstraction and calculate irrigation rate for next dt_slow
    demand_full_ug(:) = 0.  !m3
    demand_met_ug(:) = 0.   !m3
    demand_unmet_ug(:) = 0. !m3
    gw_s_abst_ug(:) = 0. ; gw_d_abst_ug(:) = 0. !m3
    gw_s_habst_ug(:) = 0. ; gw_d_habst_ug(:) = 0. !J

    ! Loop over the entire land domain by grid cell
    do l=lnd%ls, lnd%le
        ce = first_elmt(land_tile_map(l))
        ! Loop over tiles within the current grid cell
        do while(loop_over_tiles(ce,tile,k=k))
            call set_current_point(l,k)
            if (.not.associated(tile%soil)) cycle
            soil => tile%soil
            tile_demand_full = tile%frac*lnd%ug_area(l) * soil%irr_demand_ac !m2 * kg/m2 = kg
            frac = 0.
            if(tot_demand_full(l)>0.) frac = tile_demand_full/tot_demand_full(l)
            demand_left_tile = frac*irr_demand_ug(l) * DENS_H2O/(tile%frac*lnd%ug_area(l)) !m3 * kg/m3 / m2 = kg/m2
            if(is_watch_point()) then
                write(*,*) '########### irrigation checkpoint 2 ###########'
                __DEBUG1__(tot_demand_full(l))
                __DEBUG1__(tile_demand_full)
                __DEBUG1__(irr_demand_ug(l))
                __DEBUG1__(frac)
                __DEBUG1__(demand_left_tile)
            end if
            call groundwater_abstraction(soil, demand_left_tile, shallow_abst, shallow_habst, deep_abst, deep_habst)
            soil%abst_s = shallow_abst/River%dt_slow
            soil%habst_s = shallow_habst/River%dt_slow
            soil%abst_d = deep_abst/River%dt_slow
            soil%habst_d = deep_habst/River%dt_slow
            ! call check_var_range(frac*river_abst_ug(l), -0.1, 0.0, 'River abstraction check (pos frac)', 'frac*river_abst_ug(l)', WARNING)
            tot_abst = frac*lake_abst_ug(l) * DENS_H2O/(tile%frac*lnd%ug_area(l)) & !m3 * kg/m3 / m2 = kg/m2
                    +frac*river_abst_ug(l) * DENS_H2O/(tile%frac*lnd%ug_area(l)) & !kg/m2
                    +shallow_abst &  !kg/m2
                    +deep_abst !kg/m2
            if(is_watch_point()) then
                write(*,*) '########### irrigation checkpoint 2.1 (abstraction) ###########'
                __DEBUG1__(lake_abst_ug(l))
                __DEBUG1__(river_abst_ug(l))
                __DEBUG1__(shallow_abst)
                __DEBUG1__(deep_abst)
            end if
            soil%irr_rate = tot_abst/River%dt_slow !kg/(m2 s)
            tot_habst = frac*lake_habst_ug(l)/(tile%frac*lnd%ug_area(l)) & !J/m2
                        +frac*(river_abstflow_c_ug(l,2)*DENS_H2O*River%dt_slow)/(tile%frac*lnd%ug_area(l)) & ! (J m3/kg / s) * kg/m3 * s / m2 = J/m2
                        +shallow_habst & !J/m2
                        +deep_habst !J/m2
            soil%hirr_rate = tot_habst/River%dt_slow !W/m2
            if(is_watch_point()) then
                write(*,*) '########### irrigation checkpoint 3 ###########'
                __DEBUG1__(soil%abst_s)
                __DEBUG1__(soil%abst_d)
                __DEBUG1__(tot_abst)
                __DEBUG1__(soil%irr_rate)
            end if
            gw_s_abst_ug(l) = gw_s_abst_ug(l) + shallow_abst * (tile%frac*lnd%ug_area(l))/DENS_H2O !kg/m2 * m2 / kg/m3 = m3
            gw_d_abst_ug(l) = gw_d_abst_ug(l) + deep_abst * (tile%frac*lnd%ug_area(l))/DENS_H2O !kg/m2 * m2 / kg/m3 = m3
            gw_s_habst_ug(l) = gw_s_habst_ug(l) + shallow_habst * (tile%frac*lnd%ug_area(l)) !J/m2 * m2 = J
            gw_d_habst_ug(l) = gw_d_habst_ug(l) + deep_habst * (tile%frac*lnd%ug_area(l)) !J/m2 * m2 = J
            demand_full_ug(l) =  demand_full_ug(l) + soil%irr_demand_ac * (tile%frac*lnd%ug_area(l))/DENS_H2O !kg/m2 * m2 / kg/m3 = m3
            demand_met_ug(l) = demand_met_ug(l) + soil%irr_rate*River%dt_slow * (tile%frac*lnd%ug_area(l))/DENS_H2O !kg/(m2 s) * s * m2 / kg/m3 = m3
            demand_unmet_ug(l) = demand_unmet_ug(l) &
                                +(demand_left_tile-shallow_abst-deep_abst) * (tile%frac*lnd%ug_area(l))/DENS_H2O !kg/m2 * m2 / kg/m3 = m3
            soil%hlsp%irrrate_soil = tot_abst/River%dt_slow !kg/(m2 s)
            soil%hlsp%hirrrate_soil = tot_habst/River%dt_slow !W/m2
            soil%hlsp%absts_soil = shallow_abst/River%dt_slow
            soil%hlsp%habsts_soil = shallow_habst/River%dt_slow
            soil%hlsp%abstd_soil = deep_abst/River%dt_slow
            soil%hlsp%habstd_soil = deep_habst/River%dt_slow
            if(is_watch_point()) then
                write(*,*) '########### irrigation checkpoint 4 ###########'
                __DEBUG1__(demand_full_ug(l))
                __DEBUG1__(demand_met_ug(l))
                __DEBUG1__(demand_unmet_ug(l))
                __DEBUG1__(soil%hlsp%irrrate_soil)
            end if
        enddo
    enddo

    demand_full(:,:) = 0. ; demand_met(:,:) = 0. ; demand_unmet(:,:) = 0. !m3
    gw_s_abst(:,:) = 0. ; gw_d_abst(:,:) = 0. !m3
    gw_s_habst(:,:) = 0. ; gw_d_habst(:,:) = 0. !J
    call mpp_pass_UG_to_SG(lnd%ug_domain, demand_full_ug, demand_full)  !m3
    call mpp_pass_UG_to_SG(lnd%ug_domain, demand_met_ug, demand_met)  !m3
    call mpp_pass_UG_to_SG(lnd%ug_domain, demand_unmet_ug, demand_unmet)  !m3
    call mpp_pass_UG_to_SG(lnd%ug_domain, gw_s_abst_ug, gw_s_abst) !m3
    call mpp_pass_UG_to_SG(lnd%ug_domain, gw_d_abst_ug, gw_d_abst) !m3
    call mpp_pass_UG_to_SG(lnd%ug_domain, gw_s_habst_ug, gw_s_habst)  !J
    call mpp_pass_UG_to_SG(lnd%ug_domain, gw_d_habst_ug, gw_d_habst)  !J


    River%outflowmean = River%outflowmean + &
       (River%outflow-River%outflowmean)*River%dt_slow/River%channel_tau
    where(River%outflowmean .le. outflowmean_min) River%outflowmean=outflowmean_min
    call get_Leo_Mad_params(DHG_exp, DHG_coef, AAS_exp)
    do j = jsc, jec
       do i = isc, iec
          if ( River%reach_length(i,j) > 0.0) then
              River%o_coef(i,j) = River%outflowmean(i,j) / &
                   ((sinuosity*River%reach_length(i,j))*DHG_coef%on_w*DHG_coef%on_d &
                   *(River%outflowmean(i,j)**(DHG_exp%on_w+DHG_exp%on_d)))**River%o_exp
          endif
       enddo
    enddo
    River%d_coef = DHG_coef%on_d                        &
         *(River%outflowmean**(DHG_exp%on_d-AAS_exp%on_d))
    River%w_coef = DHG_coef%on_w                        &
         *(River%outflowmean**(DHG_exp%on_w-AAS_exp%on_w))

    River%stordis = River%dt_slow*River%disw2o
    do i_species = 1, num_species
       River%stordis_c(:,:,i_species) = River%dt_slow*River%disc2o(:,:,i_species)
    enddo

    rivr_FMASS = DENS_H2O * (River%storage_c(:,:,1) + River%stordis_c(:,:,1)) !kg
    rivr_LMASS = DENS_H2O * (River%storage + River%stordis) - rivr_FMASS !kg
    rivr_MELT  = DENS_H2O *  River%melt / River%dt_slow ! J m3/kg * kg/m3 / s = J/s = W
    rivr_HEAT  = DENS_H2O * (River%storage_c(:,:,2) + River%stordis_c(:,:,2)) & !kg/m3 * J m3/kg = J
                      - hlf*rivr_FMASS

    call mpp_clock_begin(diagclock)
  ! convert area-integrated river outputs to unit-area quantities,
  ! using land area for stores, cell area for fluxes to ocean
    if (id_LWSr > 0) then
       where (lnd%sg_area > 0) &
            rivr_LMASS = rivr_LMASS / lnd%sg_area !kg/m2
       used = send_data (id_LWSr, rivr_LMASS, River%time, mask=lnd%sg_area>0)
    endif
    if (id_FWSr > 0) then
       where (lnd%sg_area > 0) &
            rivr_FMASS = rivr_FMASS / lnd%sg_area !kg/m2
       used = send_data (id_FWSr, rivr_FMASS, River%time, mask=lnd%sg_area>0)
    endif
    if (id_HSr > 0) then
       where (lnd%sg_area > 0) &
            rivr_HEAT = rivr_HEAT / lnd%sg_area !J/m2
       used = send_data (id_HSr, rivr_HEAT, River%time, mask=lnd%sg_area>0)
    endif
    if (id_meltr > 0) then
       where (lnd%sg_area > 0) &
            rivr_MELT = rivr_MELT / lnd%sg_area !W/m2
       used = send_data (id_meltr, rivr_MELT, River%time, mask=lnd%sg_area>0)
    end if

    if (id_gw_s_abst > 0) then
       gw_s_abst = gw_s_abst*DENS_H2O / (River%land_area*River%dt_slow)  ! m3 * kg/m3 / (m2 s) = kg/(m2 s)
       used = send_data (id_gw_s_abst, gw_s_abst, River%time, mask=River%mask)
    end if

    if (id_gw_d_abst > 0) then
       gw_d_abst = gw_d_abst*DENS_H2O / (River%land_area*River%dt_slow)  ! m3 * kg/m3 / (m2 s) = kg/(m2 s)
       used = send_data (id_gw_d_abst, gw_d_abst, River%time, mask=River%mask)
    end if

    if (id_gw_s_habst > 0) then
       gw_s_habst = gw_s_habst / (River%land_area*River%dt_slow)  ! J / (m2 s) = W/m2
       used = send_data (id_gw_s_habst, gw_s_habst, River%time, mask=River%mask)
    end if

    if (id_gw_d_habst > 0) then
       gw_d_habst = gw_d_habst / (River%land_area*River%dt_slow)  ! J / (m2 s) = W/m2
       used = send_data (id_gw_d_habst, gw_d_habst, River%time, mask=River%mask)
    end if

    if (id_irr_full > 0) then
       demand_full = demand_full*DENS_H2O / (River%land_area*River%dt_slow)  ! m3 * kg/m3 / (m2 s) = kg/(m2 s)
       used = send_data (id_irr_full, demand_full, River%time, mask=River%mask)
    end if

    if (id_irr_met > 0) then
       demand_met = demand_met*DENS_H2O / (River%land_area*River%dt_slow)  ! m3 * kg/m3 / (m2 s) = kg/(m2 s)
       used = send_data (id_irr_met, demand_met, River%time, mask=River%mask)
    end if

    if (id_irr_unmet > 0) then
       demand_unmet = demand_unmet*DENS_H2O / (River%land_area*River%dt_slow)  ! m3 * kg/m3 / (m2 s) = kg/(m2 s)
       used = send_data (id_irr_unmet, demand_unmet, River%time, mask=River%mask)
    end if

    if (id_rsv_outflow > 0) then
       rsv_outflow = rsv_outflow / (River%land_area*River%dt_slow)  ! kg / (m2 s) = kg/(m2 s)
       used = send_data (id_rsv_outflow, rsv_outflow, River%time, mask=River%mask)
    end if


    if(mod(slow_step, diag_freq) == 0)  call river_diag(lake_depth_sill)
    call mpp_clock_end(diagclock)


  end subroutine update_river_slow

!--------------------------------------------------------
subroutine groundwater_abstraction(soil,irr_demand, abst_s, habst_s, abst_d, habst_d)

  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: irr_demand !kg/m2
  real, intent(out) :: abst_s, abst_d !kg/m2
  real, intent(out) :: habst_s, habst_d !J/m2

  integer :: lev
  real :: avail, abst_lev
  real :: abst_thres = 1.e-20 !kg/m2

  abst_s = 0.; abst_d = 0. !kg/m2
  habst_s = 0.; habst_d = 0. !J/m2
  if(.not.do_groundwater_abstraction) return
  if(irr_demand<=abst_thres) return
  do lev = 1, num_soil
    avail = soil%wl(lev)-soil%w_fc(lev)*(DENS_H2O*dz_soil(lev)) !1 * kg/m3 * m = kg/m2
    if (is_watch_point()) then
    write(*,*) '########### gw_abstraction 1 ###########'
     __DEBUG1__(irr_demand)
     __DEBUG1__(lev)
     __DEBUG1__(soil%wl(lev))
     __DEBUG1__(soil%w_fc(lev))
     __DEBUG1__(dz_soil(lev))
     __DEBUG1__(avail)
    endif
    if(avail <= 0.) cycle
    abst_lev = max(0., min(irr_demand-abst_s, avail)) !kg/m2
    soil%wl(lev) = soil%wl(lev) - abst_lev !kg/m2
    abst_s = abst_s + abst_lev !kg/m2
    habst_s = habst_s + clw*(soil%T(lev)-tfreeze)*abst_lev  !J/m2
    if (is_watch_point()) then
    write(*,*) '########### gw_abstraction 2 ###########'
     __DEBUG1__(abst_lev)
     __DEBUG1__(soil%wl(lev))
     __DEBUG1__(abst_s)
    endif
    if((irr_demand-abst_s)<=abst_thres) exit
  enddo
  if((irr_demand-abst_s)>abst_thres .and. do_deep_gw_abst)then
    abst_d = max(0., irr_demand-abst_s) !kg/m2
    habst_d = clw*(soil%T(num_soil)-tfreeze)*abst_d !J/m2
    if (is_watch_point()) then
      write(*,*) '########### gw_abstraction deep ###########'
       __DEBUG1__(abst_d)
       __DEBUG1__(irr_demand)
       __DEBUG1__(abst_s)
      endif
  endif

end subroutine groundwater_abstraction

!#####################################################################

  subroutine update_river_bnd_slow
    integer :: i_species
! note that land_area is not the total area of the cell, but just the land area
! within the cell, so it cannot be used to normalize fluxes to all-ocean cells.
! we need a true cell area for normalization, so river will
! just return mass flux per unit time and let land_model divide by area
    discharge2ocean_next = DENS_H2O*(River%disw2o - River%disc2o(:,:,1))

    do i_species = 1, num_species
       discharge2ocean_next_c(:,:,i_species) = DENS_H2O*River%disc2o(:,:,i_species)
    enddo

  end subroutine update_river_bnd_slow

!#####################################################################

  subroutine river_end
    integer :: outunit ! unit number for stdout

    if(.not.do_rivers) return ! do nothing further if rivers are turned off

!--- write out checksum
    outunit=stdout()
    write(outunit,*)"Chksum for storage ==> ", mpp_chksum(River%storage(isc:iec,jsc:jec))
    write(outunit,*)"Chksum for storage_c ==> ", mpp_chksum(River%storage_c(isc:iec,jsc:jec,:))
    write(outunit,*)"Chksum for discharge2ocean_next ==> ", mpp_chksum(discharge2ocean_next(isc:iec,jsc:jec))
    write(outunit,*)"Chksum for discharge2ocean_next_c ==> ", mpp_chksum(discharge2ocean_next_c(isc:iec,jsc:jec,:))

!--- release memory
    deallocate(discharge2ocean_next, discharge2ocean_next_c,&
         River%run_stor, River%run_stor_c)

    deallocate( River%lon, River%lat)
    deallocate(River%land_area ,     River%basinid        )
    deallocate(River%landfrac )
    deallocate(River%tocell )
    deallocate(River%travel )
    deallocate(River%outflow  )
    deallocate(River%inflow  )
    deallocate(River%lake_outflow)
    deallocate(River%lake_outflow_c)
    deallocate(River%storage        )
    deallocate(River%stordis        )
    deallocate(River%melt           )
    deallocate(River%disw2o        )
    deallocate(River%disc2o        )
    deallocate(River%infloc   )
    deallocate(River%reach_length    )
    deallocate(River%mask        )
    deallocate(River%So        )
    deallocate(River%depth     )
    deallocate(River%width     )
    deallocate(River%vel       )
    deallocate(River%infloc_c ,     River%storage_c ,     River%stordis_c    )
    deallocate(River%inflow_c, River%outflow_c )
    deallocate(River%removal_c )
    deallocate(River%d_coef,River%o_coef,River%w_coef)
    deallocate(River%threshold)
    deallocate(River%env_flow)
    deallocate(River%abst)
    deallocate(River%abstflow_c)
    deallocate(River%lake_abst)
    deallocate(River%lake_habst)

    module_is_initialized = .FALSE.

  end subroutine river_end


!#####################################################################
  !--- write to restart file
  subroutine save_river_restart(timestamp)
    character(*), intent(in) :: timestamp

    type(FmsNetcdfDomainFile_t) :: river_restart
    logical :: s
    integer :: tr
    integer, dimension(:), allocatable :: buffer
    integer :: starting, ending, i
    real, dimension(:,:,:), allocatable :: buffer3d

    if (.not. do_rivers) return ! do nothing further if rivers are turned off
    s = open_file(river_restart, 'RESTART/'//trim(timestamp)//"river.nc", &
                  "overwrite", domain, is_restart=.true.)

    call register_axis(river_restart, river_res_xdim, "x")
    call register_field(river_restart, river_res_xdim, "double", (/river_res_xdim/))
    call register_variable_string_attribute(river_restart, river_res_xdim, "long_name", river_res_xdim)
    call register_variable_string_attribute(river_restart, river_res_xdim, "units", "none")
    call register_variable_string_attribute(river_restart, river_res_xdim, "cartesian_axis", "X")
    call get_global_io_domain_indices(river_restart, river_res_xdim, starting, ending)
    allocate(buffer(ending-starting+1))
    do i = starting, ending
       buffer(i-starting+1) = i
    end do
    call write_data(river_restart, river_res_xdim, buffer)
    deallocate(buffer)

    call register_axis(river_restart, river_res_ydim, "y")
    call register_field(river_restart, river_res_ydim, "double", (/river_res_ydim/))
    call register_variable_string_attribute(river_restart, river_res_ydim, "long_name", river_res_ydim)
    call register_variable_string_attribute(river_restart, river_res_ydim, "units", "none")
    call register_variable_string_attribute(river_restart, river_res_ydim, "cartesian_axis", "Y")
    call get_global_io_domain_indices(river_restart, river_res_ydim, starting, ending)
    allocate(buffer(ending-starting+1))
    do i = starting, ending
       buffer(i-starting+1) = i
    end do
    call write_data(river_restart, river_res_ydim, buffer)
    deallocate(buffer)

    call register_axis(river_restart, river_res_zdim, 1)
    call register_field(river_restart, river_res_zdim, "double", (/river_res_zdim/))
    call register_variable_string_attribute(river_restart, river_res_zdim, "long_name", river_res_zdim)
    call register_variable_string_attribute(river_restart, river_res_zdim, "units", "none")
    call register_variable_string_attribute(river_restart, river_res_zdim, "cartesian_axis", "Z")
    call write_data(river_restart, river_res_zdim, 1)

    call register_axis(river_restart, "Time", unlimited)
    call register_field(river_restart, "Time", "double", (/"Time"/))
    call register_variable_string_attribute(river_restart, "Time", "long_name", "Time")
    call register_variable_string_attribute(river_restart, "Time", "units", "time level")
    call register_variable_string_attribute(river_restart, "Time", "cartesian_axis", "T")
    call write_data(river_restart, "Time", 1)

    allocate(buffer3d(size(river%storage,1), size(river%storage,2), 1))
    buffer3d(:,:,1) = river%storage
    call register_restart_field(river_restart, "storage", buffer3d, (/river_res_xdim, river_res_ydim, river_res_zdim, "Time"/))
    call register_variable_string_attribute(river_restart, "storage", "long_name", "storage")
    call register_variable_string_attribute(river_restart, "storage", "units", "none")
    call write_data(river_restart, "storage", buffer3d)
    deallocate(buffer3d)

    allocate(buffer3d(size(discharge2ocean_next,1), size(discharge2ocean_next,2), 1))
    buffer3d(:,:,1) = discharge2ocean_next
    call register_restart_field(river_restart, "discharge2ocean", buffer3d, (/river_res_xdim, river_res_ydim, river_res_zdim, "Time"/))
    call register_variable_string_attribute(river_restart, "discharge2ocean", "long_name", "discharge2ocean")
    call register_variable_string_attribute(river_restart, "discharge2ocean", "units", "none")
    call write_data(river_restart, "discharge2ocean", buffer3d)
    deallocate(buffer3d)

    do tr = 1, num_species
       allocate(buffer3d(size(river%storage_c(:,:,tr),1), size(river%storage_c(:,:,tr),2), 1))
       buffer3d(:,:,1) = river%storage_c(:,:,tr)
       call register_restart_field(river_restart, "storage_"//trdata(tr)%name, buffer3d, (/river_res_xdim, river_res_ydim, river_res_zdim, "Time"/))
       call register_variable_string_attribute(river_restart, "storage_"//trdata(tr)%name, "long_name", "storage_"//trdata(tr)%name)
       call register_variable_string_attribute(river_restart, "storage_"//trdata(tr)%name, "units", "none")
       call write_data(river_restart, "storage_"//trdata(tr)%name, buffer3d)
       deallocate(buffer3d)

       allocate(buffer3d(size(discharge2ocean_next_c(:,:,tr),1), size(discharge2ocean_next_c(:,:,tr),2), 1))
       buffer3d(:,:,1) = discharge2ocean_next_c(:,:,tr)
       call register_restart_field(river_restart, "disch2ocn_"//trdata(tr)%name, buffer3d, (/river_res_xdim, river_res_ydim, river_res_zdim, "Time"/))
       call register_variable_string_attribute(river_restart, "disch2ocn_"//trdata(tr)%name, "long_name", "disch2ocn_"//trdata(tr)%name)
       call register_variable_string_attribute(river_restart, "disch2ocn_"//trdata(tr)%name, "units", "none")
       call write_data(river_restart, "disch2ocn_"//trdata(tr)%name, buffer3d)
       deallocate(buffer3d)
    enddo

    allocate(buffer3d(size(river%outflowmean,1), size(river%outflowmean,2), 1))
    buffer3d(:,:,1) = river%outflowmean
    call register_restart_field(river_restart, "Omean", buffer3d, (/river_res_xdim, river_res_ydim, river_res_zdim, "Time"/))
    call register_variable_string_attribute(river_restart, "Omean", "long_name", "Omean")
    call register_variable_string_attribute(river_restart, "Omean", "units", "none")
    call write_data(river_restart, "Omean", buffer3d)
    deallocate(buffer3d)

    allocate(buffer3d(size(river%depth,1), size(river%depth,2), 1))
    buffer3d(:,:,1) = river%depth
    call register_restart_field(river_restart, "depth", buffer3d, (/river_res_xdim, river_res_ydim, river_res_zdim, "Time"/))
    call register_variable_string_attribute(river_restart, "depth", "long_name", "depth")
    call register_variable_string_attribute(river_restart, "depth", "units", "none")
    call write_data(river_restart, "depth", buffer3d)
    deallocate(buffer3d)

    call close_file(river_restart)
  end subroutine save_river_restart

!#####################################################################
  subroutine get_river_data(land_lon, land_lat, land_frac)
    real,            intent(in) :: land_lon(isc:,jsc:)  ! geographical longitude of cell center
    real,            intent(in) :: land_lat(isc:,jsc:)  ! geographical latitude of cell center
    real,            intent(in) :: land_frac(isc:,jsc:) ! land area fraction of land grid.

    integer                           :: ni, nj, i, j, ntiles
    real, dimension(:,:), allocatable :: xt, yt, frac, glon, glat, lake_frac
    integer :: nerrors ! number of errors detected during initialization
    type(FmsNetcdfFile_t) :: fileobj
    logical :: exists
    integer, dimension(:), allocatable :: siz
    integer :: isize, jsize
    integer :: ndims, L
    integer, dimension(1) :: tile_id
    real, dimension(lnd%ls:lnd%le) :: threshold, env_flow

    ntiles = mpp_get_ntile_count(domain)
    tile_id = mpp_get_tile_id(domain)

    if (ntiles>1) then
        L = len(trim(river_src_file))
        write(river_src_file, '(a,a,i1,a)') trim(river_src_file(1:L-2)), 'tile', tile_id(1), '.nc'
    endif

    exists = open_file(fileobj, river_src_file, "read")
    if (.not. exists) then
      call error_mesg("get_river_data", &
                      "file "//trim(river_src_file)//" does not exist.", &
                      FATAL)
    endif
    ndims = get_variable_num_dimensions(fileobj, "basin")
    allocate(siz(ndims))
    call get_variable_size(fileobj, "basin", siz)
    ni = siz(1)
    nj = siz(2)
    deallocate(siz)
    if(ni .NE. River%nlon .OR. nj .NE. River%nlat) call mpp_error(FATAL, &
       "river_mod: size mismatch between river grid and land grid")

    allocate(glon(ni,nj), glat(ni, nj))
    allocate(xt(isc:iec, jsc:jec), yt(isc:iec, jsc:jec), frac(isc:iec, jsc:jec) )
    allocate(lake_frac(isc:iec, jsc:jec))

    if (ntiles == 1) then
        call read_data(fileobj, "x", glon)
        call read_data(fileobj, "y", glat)
      endif
     isize = iec - isc + 1
    jsize = jec - jsc + 1
    call read_data(fileobj, "x", xt, corner=(/isc, jsc/), edge_lengths=(/isize, jsize/))
    call read_data(fileobj, "y", yt, corner=(/isc, jsc/), edge_lengths=(/isize, jsize/))
    call read_data(fileobj, "land_frac", frac, corner=(/isc, jsc/), &
                   edge_lengths=(/isize, jsize/))
    call read_data(fileobj, "lake_frac", lake_frac, &
                   corner=(/isc, jsc/), edge_lengths=(/isize, jsize/))

    !--- the following will be changed when the river data sets is finalized.
    xt = land_lon
    yt = land_lat
!--- transform to radians, since land model grid use radians and compare with land grid.

    allocate(River%lon_1d    (1:ni            ) )
    allocate(River%lat_1d    (1:nj            ) )
    allocate(River%lon       (isc:iec, jsc:jec) )
    allocate(River%lat       (isc:iec, jsc:jec) )
    allocate(River%land_area  (isc:iec, jsc:jec) )
    allocate(River%basinid   (isc:iec, jsc:jec) )
    allocate(River%landfrac  (isc:iec, jsc:jec) )
    allocate(River%mask      (isc:iec, jsc:jec) )
    allocate(River%tocell    (isc:iec, jsc:jec) )
    allocate(River%i_tocell  (isc:iec, jsc:jec) )
    allocate(River%j_tocell  (isc:iec, jsc:jec) )
    allocate(River%travel    (isd:ied, jsd:jed) )
    allocate(River%inflow    (isc:iec, jsc:jec) )
    allocate(River%outflow   (isc:iec, jsc:jec) )
    allocate(River%lake_outflow(isc:iec, jsc:jec) )
    allocate(River%storage   (isc:iec, jsc:jec) )
    allocate(River%stordis   (isc:iec, jsc:jec) )
    allocate(River%run_stor  (isc:iec, jsc:jec) )
    allocate(River%melt      (isc:iec, jsc:jec) )
    allocate(River%disw2o    (isc:iec, jsc:jec) )
    allocate(River%infloc    (isc:iec, jsc:jec))
    allocate(River%reach_length(isc:iec, jsc:jec) )
    allocate(River%So        (isc:iec, jsc:jec) )
    allocate(River%depth     (isc:iec, jsc:jec) )
    allocate(River%width    (isc:iec, jsc:jec) )
    allocate(River%vel      (isc:iec, jsc:jec) )
    allocate(River%infloc_c  (isc:iec, jsc:jec, num_species) )
    allocate(River%storage_c (isc:iec, jsc:jec, num_species) )
    allocate(River%stordis_c (isc:iec, jsc:jec, num_species) )
    allocate(River%run_stor_c (isc:iec, jsc:jec, num_species) )
    allocate(River%outflow_c (isc:iec, jsc:jec, num_species) )
    allocate(River%lake_outflow_c (isc:iec, jsc:jec, num_species) )
    allocate(River%removal_c (isc:iec, jsc:jec, num_species) )
    allocate(River%inflow_c  (isc:iec, jsc:jec, num_species) )
    allocate(River%disc2o    (isc:iec, jsc:jec, num_species))
    allocate(River%d_coef    (isc:iec, jsc:jec) )
    allocate(River%o_coef    (isc:iec, jsc:jec) )
    allocate(River%w_coef    (isc:iec, jsc:jec) )
    allocate(River%outflowmean(isc:iec, jsc:jec) )
    allocate(River%threshold  (isc:iec, jsc:jec) )
    allocate(River%env_flow  (isc:iec, jsc:jec) )
    allocate(River%abst  (isc:iec, jsc:jec) )
    allocate(River%abstflow_c (isc:iec, jsc:jec, num_species) )
    allocate(River%lake_abst (isc:iec, jsc:jec) )
    allocate(River%lake_habst (isc:iec, jsc:jec) )

    if(ntiles == 1) then   ! lat-lon grid, use actual grid location
       River%lon_1d(:)      = glon(:,1)
       River%lat_1d(:)      = glat(1,:)
    else                   ! cubic grid, use index.
       River%lon_1d(:)      = (/ (i, i=1,River%nlon) /)
       River%lat_1d(:)      = (/ (i, i=1,River%nlat) /)
    end if
    deallocate(glon, glat)

    River%lon(:,:)       = land_lon(:,:)
    River%lat(:,:)       = land_lat(:,:)
!!$    River%landfrac(:,:)  = land_frac(:,:)
    River%landfrac(:,:)  = frac(:,:)
    River%infloc    = 0.0
    River%infloc_c  = 0.0
    River%storage   = 0.0
    River%storage_c = 0.0
    River%stordis   = 0.0
    River%run_stor  = 0.0
    River%stordis_c = 0.0
    River%run_stor_c= 0.0
    River%removal_c = 0.0
    River%depth     = 0.
    River%width     = 0.
    River%vel       = 0.
    River%outflow   = 0.
    River%outflow_c = 0.
    River%inflow    = 0.
    River%inflow_c  = 0.
    River%threshold = 0.
    River%env_flow  = 0.
    River%abst      = 0.
    River%abstflow_c= 0.
    River%lake_abst = 0.
    River%lake_habst = 0.

!--- read the data from the source file
    call read_data(fileobj, "tocell", River%tocell, corner=(/isc, jsc/), &
                   edge_lengths=(/isize, jsize/))

    where (River%tocell(:,:).eq.  4) River%tocell(:,:)=3
    where (River%tocell(:,:).eq.  8) River%tocell(:,:)=4
    where (River%tocell(:,:).eq. 16) River%tocell(:,:)=5
    where (River%tocell(:,:).eq. 32) River%tocell(:,:)=6
    where (River%tocell(:,:).eq. 64) River%tocell(:,:)=7
    where (River%tocell(:,:).eq.128) River%tocell(:,:)=8

    nerrors = 0
    do j = jsc, jec
    do i = isc, iec
!!$          if(abs(xt(i,j) - land_lon(i,j)) > epsln) call mpp_error(FATAL, &
!!$             "get_river_data: longitude mismatch between river grid and land grid")
!!$          if(abs(yt(i,j) - land_lat(i,j)) > epsln) call mpp_error(FATAL, &
!!$             "get_river_data: latitude mismatch between river grid and land grid")
!!$          if(abs(frac(i,j) - land_frac(i,j)) > epsln) call mpp_error(FATAL, &
!!$             "get_river_data: area fraction mismatch between river grid and land grid")

       ! check that river and land masks match
       if ((frac(i,j)>0).neqv.(land_frac(i,j)>0)) then
          call mpp_error(WARNING,'get_river_data: land and river masks do not match at '//&
               trim(coordinates(i,j)))
          nerrors = nerrors+1
       endif

       ! check that the rivers do not discarge in the middle of the continents
       if ((River%tocell(i,j)==0).and.(land_frac(i,j)>1.0-epsln)) then
          print*,i,j,River%tocell(i,j),land_frac(i,j)
          call mpp_error(WARNING, &
               'get_river_data: river discharges into a land point '&
               //trim(coordinates(i,j))//' where there is no ocean')
          nerrors = nerrors+1
       endif
    end do
    end do

    if (nerrors>0.and.stop_on_mask_mismatch) call mpp_error(FATAL,&
        'get_river_data: river/land mask-related mismatch detected during river data initialization')

    call read_data(fileobj, "basin", River%basinid, corner=(/isc, jsc/), &
                   edge_lengths=(/isize, jsize/))
    where (River%basinid >0)
       River%mask = .true.
    elsewhere
       River%mask = .false.
    endwhere

    River%travel = 0
    call read_data(fileobj, "travel", River%travel(isc:iec,jsc:jec), &
                   corner=(/isc, jsc/), edge_lengths=(/isize, jsize/))
    call mpp_update_domains(River%travel, domain)
    call read_data(fileobj, "celllength", River%reach_length, &
                   corner=(/isc, jsc/), edge_lengths=(/isize, jsize/))
    River%reach_length = River%reach_length * River%landfrac * (1.-lake_frac)
    if (land_area_called_cellarea) then
        call read_data(fileobj, "cellarea", River%land_area, &
                       corner=(/isc, jsc/), edge_lengths=(/isize, jsize/))
      else
        call read_data(fileobj, "land_area", River%land_area, &
                       corner=(/isc, jsc/), edge_lengths=(/isize, jsize/))
      endif
!    call read_data(fileobj, "So", River%So)
    River%So = 0.0
    where (River%So .LT. 0.0) River%So = Somin
    call close_file(fileobj)

    exists = open_file(fileobj, river_threshold_file, "read")
    if(exists)then
       call read_field(fileobj, 'Threshold', threshold,fill=-1e8) !kg/m2
    !    write(*,*) 'Threshold:', threshold
    !    write(*,*) 'lnd%ug_cellarea:', lnd%ug_cellarea
    !    write(*,*) 'DENS_H2O:', DENS_H2O
       threshold = threshold*lnd%ug_cellarea/DENS_H2O !kg/m2 * m2 / kg/m3 = m3
       where (threshold<0.0) threshold = 0.0
       call mpp_pass_UG_to_SG(lnd%ug_domain,threshold,River%threshold)
     else
       River%threshold = 0.0
    endif
    call close_file(fileobj)

    exists = open_file(fileobj, env_flow_file, "read")
    if(exists)then
       call read_field(fileobj, 'Env_flow', env_flow,fill=-1e8) !kg/(m2 s)
       env_flow = env_flow*lnd%ug_cellarea/DENS_H2O !kg/(m2 s) * m2 / kg/m3 = m3/s
       where (env_flow<0.0) env_flow = 0.0
       call mpp_pass_UG_to_SG(lnd%ug_domain,env_flow,River%env_flow)
    else
       River%env_flow = 0.0
    endif
    call close_file(fileobj)

    deallocate(lake_frac)

  end subroutine get_river_data

!#####################################################################

  subroutine river_diag_init(id_lon, id_lat)
    integer, intent(in) :: id_lon  ! ID of land longitude (X) diag axis
    integer, intent(in) :: id_lat  ! ID of land latitude (Y) diag axis

    character(len=11)                :: mod_name = 'river'
    real, dimension(isc:iec,jsc:jec) :: tmp, no_riv
    logical                          :: sent
    integer                          :: i
    integer :: id_geolon_t, id_geolat_t, id_area_land, id_cellarea

! static fields
    id_geolon_t = register_static_field ( mod_name, 'geolon_t', (/id_lon,id_lat/), &
         'longitude of grid cell centers', 'degrees_E', missing_value = -1.0e+20 )
    id_geolat_t = register_static_field ( mod_name, 'geolat_t', (/id_lon,id_lat/), &
         'latitude of grid cell centers', 'degrees_N', missing_value = -1.0e+20 )
    id_cellarea = register_static_field ( mod_name, 'cell_area', (/id_lon, id_lat/), &
         'River Model Grid-Cell Area', 'm2', standard_name='cell_area', missing_value = -1.0e+20 )
    call diag_field_add_attribute(id_cellarea,'cell_methods','area: sum')

    id_area_land = register_static_field ( mod_name, 'area_land', (/id_lon, id_lat/), &
         'land area', 'm2', missing_value = -1.0e+20 )
    call diag_field_add_attribute(id_area_land,'cell_methods','area: sum')

    ! regular diagnostic fields normalized per land area; values outside of land are zeroed
    ! out. This is the traditional way of saving the river diagnostics.
    ! NOTE that for some fields it is problematic, for example when river discharges into
    ! a grid cell where there is no land, and the value is zeroed-out in the output.
    do i = 0, num_species
      id_inflow(i) = register_diag_field ( mod_name, 'rv_i_'//trim(trdata(i)%name),      &
           (/id_lon, id_lat/), River%Time, 'river inflow, '//trim(trdata(i)%longname),   &
           trdata(i)%flux_units, missing_value=missing )
      id_outflow(i) = register_diag_field ( mod_name, 'rv_o_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river outflow, '//trim(trdata(i)%longname),  &
           trdata(i)%flux_units, missing_value=missing )
      id_abstflow(i) = register_diag_field ( mod_name, 'rv_a_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river abstraction flow, '//trim(trdata(i)%longname),  &
           trdata(i)%flux_units, missing_value=missing )
      id_dis(i)     = register_diag_field ( mod_name, 'rv_d_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'ocean_discharge, '//trim(trdata(i)%longname),&
           trdata(i)%flux_units, missing_value=missing, area=id_area_land )
      call diag_field_add_attribute(id_dis(i),'cell_methods', 'area: mean')
      id_lake_outflow(i) = register_diag_field ( mod_name, 'rv_l_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'lake outflow, '//trim(trdata(i)%longname), &
           trdata(i)%flux_units, missing_value=missing )
      id_infloc(i) = register_diag_field ( mod_name, 'rv_r_'//trim(trdata(i)%name),      &
           (/id_lon, id_lat/), River%Time, 'local runoff, '//trim(trdata(i)%longname),   &
           trdata(i)%flux_units, missing_value=missing )
      id_removal(i) = register_diag_field ( mod_name, 'rv_m_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river removal, '//trim(trdata(i)%longname),  &
           trdata(i)%flux_units, missing_value=missing )
      id_storage(i) = register_diag_field ( mod_name, 'rv_s_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river storage, '//trim(trdata(i)%longname),  &
           trdata(i)%store_units, missing_value=missing, area=id_area_land )
      call diag_field_add_attribute(id_storage(i),'cell_methods', 'area: mean')
      id_stordis(i) = register_diag_field ( mod_name, 'rv_n_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river discharge lag (numerical) storage, '//trim(trdata(i)%longname), &
           trdata(i)%store_units, missing_value=missing )
      id_run_stor(i) = register_diag_field ( mod_name, 'rv_u_'//trim(trdata(i)%name),    &
           (/id_lon, id_lat/), River%Time, 'river runoff lag (numerical) storage, '//trim(trdata(i)%longname), &
           trdata(i)%store_units, missing_value=missing )
    enddo

    ! register fields normalized per cell area
    do i = 0, num_species
      id_inflow_c(i) = register_diag_field ( mod_name, 'rv_i_c_'//trim(trdata(i)%name),    &
           (/id_lon, id_lat/), River%Time, 'river inflow, '//trim(trdata(i)%longname)//', per unit cell area',   &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea)
      call diag_field_add_attribute(id_inflow_c(i),'cell_methods', 'area: mean')

      id_outflow_c(i) = register_diag_field ( mod_name, 'rv_o_c_'//trim(trdata(i)%name),   &
           (/id_lon, id_lat/), River%Time, 'river outflow, '//trim(trdata(i)%longname)//', per unit cell area',  &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_outflow_c(i),'cell_methods', 'area: mean')

      id_abstflow_c(i) = register_diag_field ( mod_name, 'rv_a_c_'//trim(trdata(i)%name),   &
           (/id_lon, id_lat/), River%Time, 'river abstraction flow, '//trim(trdata(i)%longname)//', per unit cell area',  &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_abstflow_c(i),'cell_methods', 'area: mean')

      id_dis_c(i)     = register_diag_field ( mod_name, 'rv_d_c_'//trim(trdata(i)%name),   &
           (/id_lon, id_lat/), River%Time, 'ocean_discharge, '//trim(trdata(i)%longname)//', per unit cell area',&
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_dis_c(i),'cell_methods', 'area: mean')

      id_lake_outflow_c(i) = register_diag_field ( mod_name, 'rv_l_c_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'lake outflow, '//trim(trdata(i)%longname)//', per unit cell area', &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_lake_outflow_c(i),'cell_methods', 'area: mean')

      id_infloc_c(i) = register_diag_field ( mod_name, 'rv_r_c_'//trim(trdata(i)%name),    &
           (/id_lon, id_lat/), River%Time, 'local runoff, '//trim(trdata(i)%longname)//', per unit cell area',   &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_infloc_c(i),'cell_methods', 'area: mean')

      id_removal_c(i) = register_diag_field ( mod_name, 'rv_m_c_'//trim(trdata(i)%name),   &
           (/id_lon, id_lat/), River%Time, 'river removal, '//trim(trdata(i)%longname)//', per unit cell area',  &
           trdata(i)%flux_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_removal_c(i),'cell_methods', 'area: mean')

      id_storage_c(i) = register_diag_field ( mod_name, 'rv_s_c_'//trim(trdata(i)%name),   &
           (/id_lon, id_lat/), River%Time, 'river storage, '//trim(trdata(i)%longname)//', per unit cell area',  &
           trdata(i)%store_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_storage_c(i),'cell_methods', 'area: mean')

      id_stordis_c(i) = register_diag_field ( mod_name, 'rv_n_c_'//trim(trdata(i)%name),     &
           (/id_lon, id_lat/), River%Time, 'river discharge lag (numerical) storage, '//trim(trdata(i)%longname)//', per unit cell area', &
           trdata(i)%store_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_stordis_c(i),'cell_methods', 'area: mean')

      id_run_stor_c(i) = register_diag_field ( mod_name, 'rv_u_c_'//trim(trdata(i)%name),    &
           (/id_lon, id_lat/), River%Time, 'river runoff lag (numerical) storage, '//trim(trdata(i)%longname)//', per unit cell area', &
           trdata(i)%store_units, missing_value=missing, area = id_cellarea )
      call diag_field_add_attribute(id_run_stor_c(i),'cell_methods', 'area: mean')
    enddo

    id_lake_depth_sill= register_diag_field ( mod_name, 'rv_dsill', (/id_lon, id_lat/), &
         River%Time, 'effective lake sill depth', 'm', missing_value=missing )
    id_outflowmean   = register_diag_field ( mod_name, 'rv_Qavg', (/id_lon, id_lat/), &
         River%Time, 'long-time average vol. flow', 'm3/s', missing_value=missing )
    id_depth     = register_diag_field ( mod_name, 'rv_depth', (/id_lon, id_lat/), &
         River%Time, 'river flow depth', 'm', missing_value=missing )
    id_width     = register_diag_field ( mod_name, 'rv_width', (/id_lon, id_lat/), &
         River%Time, 'river flow width', 'm', missing_value=missing )
    id_vel       = register_diag_field ( mod_name, 'rv_veloc', (/id_lon, id_lat/), &
         River%Time, 'river flow velocity', 'm/s', missing_value=missing )


    id_lake_abst   = register_diag_field ( mod_name, 'lake_abst', (/id_lon, id_lat/), &
         River%Time, 'lake abstraction rate over land', 'kg/(m2 s)', missing_value=missing )
    id_lake_habst   = register_diag_field ( mod_name, 'lake_habst', (/id_lon, id_lat/), &
         River%Time, 'heat associated with lake abstraction over land', 'W/m2', missing_value=missing )
    id_gw_s_abst   = register_diag_field ( mod_name, 'gw_s_abst', (/id_lon, id_lat/), &
         River%Time, 'shallow groundwater abstraction rate over land', 'kg/(m2 s)', missing_value=missing )
    id_gw_d_abst   = register_diag_field ( mod_name, 'gw_d_abst', (/id_lon, id_lat/), &
         River%Time, 'deep groundwater abstraction rate over land', 'kg/(m2 s)', missing_value=missing )
    id_gw_s_habst   = register_diag_field ( mod_name, 'gw_s_habst', (/id_lon, id_lat/), &
         River%Time, 'heat associated with shallow groundwater abstraction over land', 'W/m2', missing_value=missing )
    id_gw_d_habst   = register_diag_field ( mod_name, 'gw_d_habst', (/id_lon, id_lat/), &
         River%Time, 'heat associated with deep groundwater abstraction over land', 'W/m2', missing_value=missing )
    id_irr_full   = register_diag_field ( mod_name, 'irr_full', (/id_lon, id_lat/), &
         River%Time, 'needed irrigation rate over land', 'kg/(m2 s)', missing_value=missing )
    id_irr_met   = register_diag_field ( mod_name, 'irr_met', (/id_lon, id_lat/), &
         River%Time, 'met irrigation rate over land', 'kg/(m2 s)', missing_value=missing )
    id_irr_unmet   = register_diag_field ( mod_name, 'irr_unmet', (/id_lon, id_lat/), &
         River%Time, 'unmet irrigation rate over land', 'kg/(m2 s)', missing_value=missing )

    id_rsv_outflow   = register_diag_field ( mod_name, 'rsv_outflow', (/id_lon, id_lat/), &
         River%Time, 'reservoir outflow', 'kg/(m2 s)', missing_value=missing )


    id_LWSr   = register_diag_field ( mod_name, 'LWSr', (/id_lon, id_lat/), &
         River%Time, 'river liquid mass storage', 'kg/m2', missing_value=-1.0e+20, area=id_area_land )
    call diag_field_add_attribute(id_LWSr,'cell_methods', 'area: mean')

    id_FWSr   = register_diag_field ( mod_name, 'FWSr', (/id_lon, id_lat/), &
         River%Time, 'river ice mass storage', 'kg/m2', missing_value=-1.0e+20, area=id_area_land )
    call diag_field_add_attribute(id_FWSr,'cell_methods', 'area: mean')

    id_HSr   = register_diag_field ( mod_name, 'HSr', (/id_lon, id_lat/), &
         River%Time, 'river heat storage', 'J/m2', missing_value=-1.0e+20, area=id_area_land )
    call diag_field_add_attribute(id_HSr,'cell_methods', 'area: mean')

    id_meltr   = register_diag_field ( mod_name, 'meltr', (/id_lon, id_lat/), &
         River%Time, 'melt in river system', 'kg/m2/s', missing_value=-1.0e+20, area=id_area_land )
    call diag_field_add_attribute(id_meltr,'cell_methods', 'area: mean')

    id_dis_liq   = register_diag_field ( mod_name, 'dis_liq', (/id_lon, id_lat/), &
         River%time, 'liquid discharge to ocean', 'kg/(m2 s)', missing_value=-1.0e+20, area=id_cellarea )
    call diag_field_add_attribute(id_dis_liq,'cell_methods', 'area: mean')
    id_dis_ice   = register_diag_field ( mod_name, 'dis_ice', (/id_lon, id_lat/), &
         River%time, 'ice discharge to ocean', 'kg/(m2 s)', missing_value=-1.0e+20, area=id_cellarea )
    call diag_field_add_attribute(id_dis_ice,'cell_methods', 'area: mean')
    id_dis_heat   = register_diag_field ( mod_name, 'dis_heat', (/id_lon, id_lat/), &
         River%time, 'heat of mass discharge to ocean', 'W/m2', missing_value=-1.0e+20, area=id_cellarea )
    call diag_field_add_attribute(id_dis_heat,'cell_methods', 'area: mean')
    id_dis_sink   = register_diag_field ( mod_name, 'dis_sink', (/id_lon, id_lat/), &
         River%time, 'burial rate of small/negative discharge', 'kg/(m2 s)', missing_value=-1.0e+20, area=id_cellarea )
    call diag_field_add_attribute(id_dis_sink,'cell_methods', 'area: mean')
    id_dis_DOC    = register_diag_field ( mod_name, 'dis_DOC', (/id_lon, id_lat/), &
         River%time, 'DOC discharge to ocean', 'kgC/m^2/s', missing_value=-1.0e+20 )
    call diag_field_add_attribute(id_dis_sink,'cell_methods', 'area: mean')

    ! static fields

    id_dx = register_static_field ( mod_name, 'rv_length', (/id_lon, id_lat/), &
           'river reach length', 'm', missing_value=missing )
    id_basin = register_static_field ( mod_name, 'rv_basin', (/id_lon, id_lat/), &
           'river basin id', 'none', missing_value=missing )
    id_So = register_static_field ( mod_name, 'So', (/id_lon, id_lat/), &
           'Slope', 'none', missing_value=missing )
    id_travel = register_static_field ( mod_name, 'rv_trav', (/id_lon, id_lat/), &
           'cells left to travel before reaching ocean', 'none', missing_value=missing )
    id_tocell = register_static_field ( mod_name, 'rv_dir', (/id_lon, id_lat/), &
           'outflow direction code', 'none', missing_value=missing )
    id_no_riv = register_static_field ( mod_name, 'no_riv', (/id_lon, id_lat/), &
         'indicator of land without rivers','unitless', missing_value=-1.0 )

    if(id_geolon_t>0) sent=send_data(id_geolon_t, River%lon*180.0/PI, River%time )
    if(id_geolat_t>0) sent=send_data(id_geolat_t, River%lat*180.0/PI, River%time )
    if(id_area_land>0) sent=send_data(id_area_land, lnd%sg_area, River%time )
    if(id_cellarea>0) sent=send_data(id_cellarea, lnd%sg_cellarea, River%time )

    if (id_dx>0) sent=send_data(id_dx, River%reach_length, River%Time, mask=River%mask )
    if (id_basin>0) then
        tmp = River%basinid(isc:iec,jsc:jec)
        sent=send_data(id_basin, tmp, River%Time, mask=River%mask )
      end if
    if (id_So>0) sent=send_data(id_So, River%So, River%Time, mask=River%mask )
    if (id_travel>0) then
        tmp = River%travel(isc:iec,jsc:jec)
        sent=send_data(id_travel, tmp, River%Time, mask=River%mask )
      end if
    if (id_tocell>0) then
        tmp = River%tocell(isc:iec,jsc:jec)
        sent=send_data(id_tocell, tmp, River%Time, mask=River%mask )
      end if

    no_riv = 0.
    where ( lnd%sg_landfrac .gt. 0 .and. .not. River%mask ) no_riv = 1.
    if ( id_no_riv > 0 ) sent = send_data( id_no_riv, no_riv, River%time )


  end subroutine river_diag_init

!#####################################################################

  subroutine river_diag(lake_depth_sill)
    real, dimension(isc:iec,jsc:jec), intent(in) :: lake_depth_sill
    logical :: used   ! logical for send_data
    real diag_factor  (isc:iec,jsc:jec)
    real diag_factor_2(isc:iec,jsc:jec)
    integer :: tr ! iterator over river tracers

    diag_factor   = DENS_H2O/lnd%sg_cellarea(:,:) !kg/m3 / m2
    diag_factor_2 = 1.0/(lnd%sg_cellarea(:,:)*River%dt_slow) ! 1/(m2 s)

    if (id_inflow_c(0) > 0) used = send_data (id_inflow_c(0), &
            diag_factor*River%inflow(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_outflow_c(0) > 0) used = send_data (id_outflow_c(0), &
            diag_factor*River%outflow(isc:iec,jsc:jec), River%Time, mask=River%mask ) !kg/m3 / m2 * m3/s = kg/(m2 s)
    if (id_storage_c(0) > 0) used = send_data (id_storage_c(0), &
            diag_factor*River%storage(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_stordis_c(0) > 0) used = send_data (id_stordis_c(0), &
            diag_factor*River%stordis(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_run_stor_c(0) > 0) used = send_data (id_run_stor_c(0), &
            River%dt_fast*River%run_stor(isc:iec,jsc:jec)/lnd%sg_cellarea, River%Time, mask=River%mask )
    if (id_infloc_c(0) > 0) used = send_data (id_infloc_c(0), &
            diag_factor*River%infloc(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_dis_c(0) > 0)    used = send_data (id_dis_c(0), &
            diag_factor*River%disw2o(isc:iec,jsc:jec), River%Time)
    if (id_lake_outflow_c(0) > 0) used = send_data (id_lake_outflow_c(0), &
            diag_factor_2*River%lake_outflow(isc:iec,jsc:jec), River%Time, mask=River%mask ) !kg/(m2 s), River%lake_outflow: kg
    if (id_abstflow_c(0) > 0) used = send_data (id_abstflow_c(0), &
            diag_factor_2*River%abst(isc:iec,jsc:jec)*DENS_H2O, River%Time, mask=River%mask ) !m3 * kg/m3  /(m2 s) = kg/(m2 s)

    do tr = 1, num_species
       if (id_outflow_c(tr) > 0) used = send_data (id_outflow_c(tr), &
         diag_factor*River%outflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask ) !m3/s * kg/m3 / m2 = kg/(m2 s), J m3/kg / s  * kg/m3 / m2 = W/m2
       if (id_abstflow_c(tr) > 0) used = send_data (id_abstflow_c(tr), &
         diag_factor*River%abstflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )  !abstflow_c units are same as River%outflow_c
       if (id_lake_outflow_c(tr) > 0) used = send_data (id_lake_outflow_c(tr), &
         diag_factor_2*River%lake_outflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_inflow_c(tr) > 0) used = send_data (id_inflow_c(tr), &
         diag_factor*River%inflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_storage_c(tr) > 0) used = send_data (id_storage_c(tr), &
         diag_factor*River%storage_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_stordis_c(tr) > 0) used = send_data (id_stordis_c(tr), &
         diag_factor*River%stordis_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_run_stor_c(tr) > 0) used = send_data (id_run_stor_c(tr), &
         River%dt_fast*River%run_stor_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_infloc_c(tr) > 0) used = send_data (id_infloc_c(tr), &
         diag_factor*River%infloc_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_removal_c(tr) > 0) used = send_data (id_removal_c(tr), &
         diag_factor*River%removal_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_dis_c(tr) > 0)    used = send_data (id_dis_c(tr), &
         diag_factor*River%disc2o(isc:iec,jsc:jec,tr), River%Time)
    enddo

    ! recalculate area normalization factors and send data normalized per land area.
    diag_factor = 0.
    diag_factor_2 = 0.
    where (River%land_area(isc:iec,jsc:jec).gt.0.) &
                     diag_factor=DENS_H2O/River%land_area(isc:iec,jsc:jec)
    where (River%land_area(isc:iec,jsc:jec).gt.0.) &
                     diag_factor_2=1./(River%land_area(isc:iec,jsc:jec)*River%dt_slow)

    if (id_inflow(0) > 0) used = send_data (id_inflow(0), &
            diag_factor*River%inflow(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_outflow(0) > 0) used = send_data (id_outflow(0), &
            diag_factor*River%outflow(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_storage(0) > 0) used = send_data (id_storage(0), &
            diag_factor*River%storage(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_stordis(0) > 0) used = send_data (id_stordis(0), &
            diag_factor*River%stordis(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_run_stor(0) > 0) used = send_data (id_run_stor(0), &
            River%dt_fast*River%run_stor(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_infloc(0) > 0) used = send_data (id_infloc(0), &
            diag_factor*River%infloc(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_dis(0) > 0)    used = send_data (id_dis(0), &
            diag_factor*River%disw2o(isc:iec,jsc:jec), River%Time)
    if (id_lake_outflow(0) > 0) used = send_data (id_lake_outflow(0), &
            diag_factor_2*River%lake_outflow(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_abstflow(0) > 0) used = send_data (id_abstflow(0), &
            diag_factor_2*River%abst(isc:iec,jsc:jec)*DENS_H2O, River%Time, mask=River%mask ) !m3 * kg/m3  /(m2 s) = kg/(m2 s)

    do tr = 1, num_species
       if (id_outflow(tr) > 0) used = send_data (id_outflow(tr), &
         diag_factor*River%outflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_abstflow(tr) > 0) used = send_data (id_abstflow(tr), &
         diag_factor*River%abstflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )  !abstflow_c units are same as River%outflow_c
       if (id_lake_outflow(tr) > 0) used = send_data (id_lake_outflow(tr), &
         diag_factor_2*River%lake_outflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_inflow(tr) > 0) used = send_data (id_inflow(tr), &
         diag_factor*River%inflow_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_storage(tr) > 0) used = send_data (id_storage(tr), &
         diag_factor*River%storage_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_stordis(tr) > 0) used = send_data (id_stordis(tr), &
         diag_factor*River%stordis_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_run_stor(tr) > 0) used = send_data (id_run_stor(tr), &
         River%dt_fast*River%run_stor_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_infloc(tr) > 0) used = send_data (id_infloc(tr), &
         diag_factor*River%infloc_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_removal(tr) > 0) used = send_data (id_removal(tr), &
         diag_factor*River%removal_c(isc:iec,jsc:jec,tr), River%Time, mask=River%mask )
       if (id_dis(tr) > 0)    used = send_data (id_dis(tr), &
         diag_factor*River%disc2o(isc:iec,jsc:jec,tr), River%Time)
    enddo

    if (id_lake_depth_sill > 0) used = send_data (id_lake_depth_sill, &
            lake_depth_sill, River%Time, mask=River%mask )
    if (id_outflowmean > 0) used = send_data (id_outflowmean, &
            River%outflowmean(isc:iec,jsc:jec), River%Time, mask=River%mask ) !m3/s
    if (id_width > 0) used = send_data (id_width, &
            River%width(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_depth > 0) used = send_data (id_depth, &
            River%depth(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_vel > 0) used = send_data (id_vel, &
            River%vel(isc:iec,jsc:jec), River%Time, mask=River%mask )
    if (id_lake_abst > 0) used = send_data (id_lake_abst, &
            diag_factor_2*River%lake_abst(isc:iec,jsc:jec)*DENS_H2O, River%Time, mask=River%mask )   !kg/(m2 s), River%lake_abst: m3
    if (id_lake_habst > 0) used = send_data (id_lake_habst, &
            diag_factor_2*River%lake_habst(isc:iec,jsc:jec), River%Time, mask=River%mask )   ! J / (m2 s) = W/m2 River%lake_habst: J

  end subroutine river_diag

!#####################################################################

  subroutine get_Leo_Mad_params(DHG_exp, DHG_coef, AAS_exp)

    type(Leo_Mad_trios), intent(inout) :: DHG_exp  ! Exponents for downstream equations
    type(Leo_Mad_trios), intent(inout) :: DHG_coef ! Coefficients for downstream equations
    type(Leo_Mad_trios), intent(inout) :: AAS_exp  ! Exponents for at-a-station equations

!!! Exponents for the downstream hydraulic geometry equations
    DHG_exp%on_w = ave_DHG_exp(1)
    DHG_exp%on_d = ave_DHG_exp(2)
    DHG_exp%on_V = ave_DHG_exp(3)

!!! Coefficients for the downstream hydraulic geometry equations
    DHG_coef%on_w = ave_DHG_coef(1)
    DHG_coef%on_d = ave_DHG_coef(2)
    DHG_coef%on_V = ave_DHG_coef(3)

!!! Exponents for the at-a-station hydraulic geometry equations
    AAS_exp%on_w = ave_AAS_exp(1)
    AAS_exp%on_d = ave_AAS_exp(2)
    AAS_exp%on_V = ave_AAS_exp(3)

  end subroutine get_Leo_Mad_params

!#####################################################################

subroutine river_stock_pe(index, value)
integer, intent(in)  :: index
real   , intent(out) :: value ! Domain water (Kg) or heat (Joules)

value = 0.0
if (.not.do_rivers) return

select case(index)
case(ISTOCK_WATER)
  value = DENS_H2O*(sum(River%storage)+sum(River%stordis)) &
        + sum(River%run_stor*River%land_area)*River%dt_fast

case(ISTOCK_HEAT)
! heat stock not yet implemented
  value = 0
case default
! Lnd_stock_pe issues a FATAL error message if index is invalid
end select

end subroutine river_stock_pe

!#####################################################################
! returns total amount of water (liquid and frozen) in rivers, kg/m2 of land
subroutine get_river_water(water)
  real, intent(out) :: water(lnd%is:lnd%ie,lnd%js:lnd%je)

  if (do_rivers) then
     water(:,:) = River%run_stor*River%dt_fast
     where (lnd%sg_area > 0.0) &
         water(:,:) = water(:,:) + DENS_H2O*(River%storage+River%stordis)/lnd%sg_area
  else
     water(:,:) = 0.0
  endif
end subroutine get_river_water

!#####################################################################
! returns string indicating the coordinates of the point i,j
function coordinates(i,j) result(s); character(128) :: s
   integer, intent(in) :: i,j
   s ='('//trim(string(i))//','//trim(string(j))//')'
   if (lnd%nfaces>1) s=trim(s)//' on cubic sphere face '//string(lnd%sg_face)
end function coordinates

end module river_mod
