! this is a copy of the testing program that was originally embedded at the
! end of river.F90. This has not been used for long time, and is unlikely
! to be functional.
program river_solo
  use mpp_mod,                  only : mpp_error, mpp_pe, mpp_root_pe, mpp_npes, FATAL
  use mpp_mod,                  only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_domains_mod,          only : mpp_define_layout, mpp_define_domains
  use mpp_domains_mod,          only : mpp_get_compute_domain, domain2d, CYCLIC_GLOBAL_DOMAIN
  use mpp_domains_mod,          only : mpp_get_current_ntile, mpp_get_tile_id
  use fms_mod,                  only : fms_init, fms_end, stdlog
  use fms_mod,                  only : check_nml_error, stdout, input_nml_file
  use time_manager_mod,         only : time_type, increment_time, set_date, increment_date, set_time
  use time_manager_mod,         only : set_calendar_type, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only : operator(/), operator(-), operator( + ), month_name, get_date
  use diag_manager_mod,         only : diag_manager_init, diag_manager_end
  use river_mod,                only : river_init, river_end, update_river
  use constants_mod,            only : constants_init, PI, radius
  use grid_mod,                 only : get_grid_size, get_grid_cell_vertices
  use grid_mod,                 only : get_grid_cell_centers, get_grid_cell_area, get_grid_comp_area
  use grid_mod,                 only : define_cube_mosaic, get_grid_ntiles
  use river_mod,                only : num_species
  use fms2_io_mod, only: close_file, file_exists, FmsNetcdfFile_t, open_file

  implicit none

  real, parameter       :: CONST_RUNOFF = 200.0

!--- namelist -----------------------------------------

  integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=16)     :: calendar = 'julian'
  integer               :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0
  integer               :: dt_fast     = 0
  integer               ::layout(2) = (/1,0/)
  namelist /river_solo_nml/ current_date, dt_fast, years, months, days, &
       hours, minutes, seconds, calendar, layout

!--------------------------------------------------------------------
  type(time_type)      :: Time, Time_start, Time_end, Run_len, Time_step_fast
  integer              :: nf, num_fast_step, unit
  integer              :: yr,mon,day,hr,min,sec, calendar_type=-1
  integer              :: outunit
  real, allocatable    :: runoff(:,:), discharge(:,:)
  integer              :: initClock, mainClock, termClock, updateClock
!Balaji
  real, allocatable :: runoff_c(:,:,:)
  real, allocatable :: discharge2ocean(:,:)
  real, allocatable :: discharge2ocean_c(:,:,:)
  real, allocatable :: liq(:,:), sol(:,:), mel(:,:), hea(:,:)

  call fms_init

  initClock = mpp_clock_id( 'Initialization' )
  mainClock = mpp_clock_id( 'Main loop' )
  termClock = mpp_clock_id( 'Termination' )
  updateClock = mpp_clock_id( 'update river')

  call mpp_clock_begin(initClock)
  call river_solo_init
  call mpp_clock_end (initClock) !end initialization

  call mpp_clock_begin(mainClock) !begin main loop
  outunit=stdout()
  do nf = 1, num_fast_step
     write(outunit,*)' at river fast time step ', nf
     call mpp_clock_begin(updateClock)
!Balaji
     call update_river( runoff, runoff_c, discharge2ocean, discharge2ocean_c, &
     liq, sol, mel, hea )
!     call update_river ( runoff, discharge )
     call mpp_clock_end (updateClock)
     Time = Time + Time_step_fast
  enddo
  call mpp_clock_end(mainClock)

  call mpp_clock_begin(termClock)
  call river_end
  call diag_manager_end(Time)
  call get_date(Time,yr,mon,day,hr,min,sec)

  if (mpp_pe() == mpp_root_pe()) then
      open(newunit=unit, file="RESTART/river_solo.res", action="write")
      write(unit,*) yr, mon, day, hr, min, sec
      write(unit,*) calendar_type
      close(unit)
  endif

  call fms_end


contains

!#######################################################################
  subroutine river_solo_init

    integer                     :: unit, ierr, io
    integer                     :: ni, nj, npes, isc, iec, jsc, jec, ntiles, tile
    type(domain2d)              :: Domain
    integer, allocatable        :: tile_ids(:)             ! mosaic tile IDs for the current PE
    real,   allocatable         :: lon(:,:), lat(:,:)
    real,   allocatable         :: area_lnd(:,:), area_lnd_cell(:,:), gfrac(:,:)
    integer                     :: date(6)
    character(len=9)            :: month
    type(FmsNetcdfFile_t ) :: fileobj
    logical :: exists

    call constants_init

    read (input_nml_file, nml=river_solo_nml, iostat=io)
    ierr = check_nml_error(io, 'river_solo_nml')

    unit=stdlog()
    write(unit, nml= river_solo_nml)

! set the calendar
    if (calendar(1:6) == 'julian') then
        calendar_type = julian
    else if (calendar(1:6) == 'NOLEAP') then
        calendar_type = NOLEAP
    else if (calendar(1:10) == 'thirty_day') then
        calendar_type = THIRTY_DAY_MONTHS
    else if (calendar(1:11) == 'no_calendar') then
        calendar_type = NO_CALENDAR
    else if (calendar(1:1) /= ' ') then
        call mpp_error (FATAL,'==>Error from ocean_solo_mod: invalid namelist value for calendar')
    else
        call mpp_error (FATAL,'==>Error from ocean_solo_mod: no namelist value for calendar')
    endif

! get river_solo restart
    if (file_exists("INPUT/river_solo.res")) then
        open(newunit=unit, file="INPUT/river_solo.res", action="read")
        read(unit,*) date
        read(unit,*) calendar_type
        close(unit)
    endif

    call set_calendar_type (calendar_type)

    call diag_manager_init

    if (sum(current_date) <= 0) then
        call mpp_error(FATAL,'==>Error from river_solo_mod: no namelist value for current date')
    else
        Time_start  = set_date(current_date(1),current_date(2), current_date(3), &
             current_date(4),current_date(5),current_date(6))
    endif

    if (file_exists('INPUT/river_solo.res')) then
        Time_start =  set_date(date(1),date(2),date(3),date(4),date(5),date(6))
    else
        Time_start = Time_start
        date = current_date
    endif

    Time           = Time_start
    Time_end       = increment_date(Time_start, years, months, days, hours, minutes, seconds)
    Run_len        = Time_end - Time_start
    Time_step_fast = set_time(dt_fast, 0)
    num_fast_step  = Run_len/Time_step_fast

      if (mpp_pe() .eq. mpp_root_pe()) then
      open(newunit=unit, file="time_stamp.out", action="write")
      month = month_name(current_date(2))
      write (unit,'(6i4,2x,a3)') date, month(1:3)
      call get_date (Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
      month = month_name(date(2))
      write (unit,'(6i4,2x,a3)') date, month(1:3)
      close(unit)
    endif

!--- get the land grid and set up domain decomposition
    call get_grid_size('LND', 1, ni, nj)

    npes = mpp_npes()
!--- define domain ------------------------------------------------
    call get_grid_ntiles('LND',ntiles)
    if(layout(1)*layout(2)*ntiles .NE. npes) call mpp_define_layout((/1,ni,1,nj/),npes/ntiles,layout)
    if (ntiles==1) then
       call mpp_define_domains ((/1,ni, 1, nj/), layout, domain, &
            xflags = CYCLIC_GLOBAL_DOMAIN, whalo=1, ehalo=1, shalo=1, nhalo=1, name = 'LAND MODEL')
    else
       call define_cube_mosaic ('LND', domain, layout, halo=1 )
    endif
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)

    allocate(tile_ids(mpp_get_current_ntile(Domain)))
    tile_ids = mpp_get_tile_id(domain)
    tile = tile_ids(1)
    deallocate(tile_ids)

!--- get grid information
    allocate( lon(isc:iec,jsc:jec), lat(isc:iec,jsc:jec) )
    allocate( area_lnd(ni,nj), area_lnd_cell(ni,nj), gfrac(ni,nj) )
    call get_grid_cell_centers ('LND', tile, lon, lat, domain)
!!$    call get_grid_cell_area    ('LND',tile, area_lnd_cell)
!!$    call get_grid_comp_area    ('LND',tile, area_lnd)

    lon = lon * PI/180.
    lat = lat * PI/180.
!!$    gfrac = area_lnd/area_lnd_cell
    gfrac = 1
    npes = mpp_npes()

    call river_init( lon, lat, Time_start, Time_step_fast, Domain, gfrac(isc:iec,jsc:jec)  )

    allocate(runoff_c(isc:iec,jsc:jec,num_species) )
    allocate(runoff(isc:iec,jsc:jec), discharge(isc:iec,jsc:jec) )
   exists = open_file(fileobj, "INPUT/runoff.nc", "read")
    if (exists) then
       call read_data(fileobj, "runoff", runoff)
       call read_data(fileobj, "runoff_c", runoff_c)
       call close_file(fileobj)
    else
       runoff = CONST_RUNOFF
       runoff_c = CONST_RUNOFF
    end if
    allocate( discharge2ocean(isc:iec,jsc:jec) )
    allocate( discharge2ocean_c(isc:iec,jsc:jec,num_species) )

  end subroutine river_solo_init

!#####################################################################

end program river_solo
