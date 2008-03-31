module river_mod
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
  ! <CONTACT EMAIL="klf@gfdl.noaa.gov"> Kirsten Findell </CONTACT> 
  ! <CONTACT EMAIL="z1l@gfdl.noaa.gov"> Zhi Liang </CONTACT> 

  use mpp_mod,             only : mpp_error, mpp_chksum, FATAL, WARNING, stdlog, mpp_npes, mpp_sum
  use mpp_mod,             only : mpp_pe, mpp_root_pe, stdout, mpp_get_current_pelist
  use mpp_mod,             only : mpp_sync_self, mpp_broadcast, mpp_sync
  use mpp_mod,             only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, MPP_CLOCK_DETAILED
  use mpp_domains_mod,     only : domain2d, mpp_get_compute_domain, mpp_define_layout, mpp_define_domains
  use mpp_io_mod,          only : mpp_open, mpp_read, axistype, fieldtype, mpp_close, mpp_get_axes
  use mpp_io_mod,          only : mpp_get_info, mpp_get_atts, mpp_get_axis_data, mpp_get_fields
  use mpp_io_mod,          only : MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE
  use fms_mod,             only : write_version_number, open_namelist_file, check_nml_error
  use fms_mod,             only : close_file, file_exist, field_size, read_data, write_data, lowercase
  use fms_mod,             only : field_exist
  use axis_utils_mod,      only : get_axis_cart, nearest_index
  use diag_manager_mod,    only : diag_axis_init, register_diag_field, register_static_field, send_data
  use time_manager_mod,    only : time_type, increment_time, get_time
  use horiz_interp_mod,    only : horiz_interp_type, horiz_interp_new, horiz_interp, horiz_interp_del
  use river_type_mod,      only : river_type, Leo_Mad_trios
  use river_physics_mod,   only : river_physics_step, river_physics_init
  use constants_mod,       only : PI, RADIAN, tfreeze, DENS_H2O
  use time_interp_external_mod, only : init_external_field, time_interp_external, time_interp_external_init

  implicit none
  private

  !--- version information ---------------------------------------------
  character(len=128) :: version = '$Id: river.F90,v 15.0.2.7 2007/12/19 23:25:41 slm Exp $'
  character(len=128) :: tagname = '$Name: omsk_2008_03 $'

  !--- public interface ------------------------------------------------
  public :: river_init, river_end, river_type, update_river

  !--- namelist interface ----------------------------------------------
  logical            :: do_rivers   = .TRUE. ! if FALSE, rivers are essentially turned off to save computing time
  real               :: dt_slow
  integer            :: diag_freq   = 1                 ! Number of slow time steps between sending out diagnositics data.
  character(len=32)  :: runoff_type = "from_land_model" ! with value "from_land_model", "read_on_land" or "read_on_river"
  logical            :: debug_river      = .FALSE.
  logical            :: do_age           = .false.
  real :: Somin = 0.00005 ! There are 7 points with So = -9.999 but basinid > 0....
  real :: outflowmean_min = 1. ! temporary fix, should not allow zero in input file
  integer                         :: num_c, num_species
  character(len=6), dimension(10) :: rt_c_name
  character(len=128),dimension(10) :: rt_source_conc_file, rt_source_flux_file
  character(len=128),dimension(10) :: rt_source_conc_name, rt_source_flux_name
  real,             dimension(10) :: rt_t_ref, rt_vf_ref, rt_q10, rt_kinv
  character(len=6),  allocatable, dimension(:) :: c_name
  character(len=128),allocatable, dimension(:) :: source_conc_file, source_flux_file
  character(len=128),allocatable, dimension(:) :: source_conc_name, source_flux_name
  real, dimension(3) :: ave_DHG_exp = (/0.49,0.33,0.18/)  ! (/B, F, M for avg of many rivers, 15Nov05/)
  real, dimension(3) :: ave_AAS_exp = (/0.19,0.39,0.42/)  ! (/b, f, m for avg of many rivers, 15Nov05/)
  real, dimension(3) :: ave_DHG_coef = (/4.62,0.26,0.82/) ! (/A, C, K for avg of many rivers, 15Nov05/)
  real               :: sinuosity = 1.3
  namelist /river_nml/ dt_slow, diag_freq, runoff_type, debug_river, do_age, &
                         Somin, outflowmean_min, &
                         num_c, rt_c_name, rt_t_ref, rt_vf_ref, rt_q10, &
                                     rt_kinv, &
                        rt_source_conc_file, rt_source_flux_file, &
                        rt_source_conc_name, rt_source_flux_name, &
                        ave_DHG_exp, ave_AAS_exp, ave_DHG_coef, &
                        do_rivers, &
                        sinuosity

  character(len=128) :: runoff_file      = "INPUT/river_runoff.nc"
  character(len=128) :: river_src_file   = 'INPUT/river_data.nc'
  character(len=128) :: river_Omean_file = 'INPUT/river_Omean.nc'
!  character(len=128) :: river_Nload_file = 'INPUT/river_Nload.nc'
  !---------------------------------------------------------------------
  integer :: is_runoff, ie_runoff, js_runoff, je_runoff ! index in runoff data grid.
  integer :: is_model,  ie_model,  js_model,  je_model  ! index in river model grid.
  integer :: id_runoff_forcing                          ! id related to time_interp_external.
  logical :: module_is_initialized = .FALSE.
  integer :: isc, iec, jsc, jec                         ! domain decomposition of river grid
  integer :: isc_lnd, iec_lnd, jsc_lnd, jec_lnd         ! domain decomposition of land grid
  integer :: nlon_lnd, nlat_lnd                         ! size of land grid
  integer :: nlon, nlat                                 ! size of computational river grid 
  integer :: id_inflow, id_outflow, id_infloc
  integer :: id_storage, id_dx, id_basin, id_So, id_depth !, id_width, id_vel
  integer :: id_disw2o, id_disw2l
  integer :: i_species
  integer :: id_travel, id_elev, id_tocell
  integer :: maxtravel
  real    :: missing = -1.e8

  integer            :: runoff_source                   ! indicate where to obtain runoff data.
  integer, parameter :: FROM_LAND_MODEL = 1
  integer, parameter :: READ_ON_LAND    = 2
  integer, parameter :: READ_ON_RIVER   = 3
  integer, parameter :: num_phys = 2
  real,    parameter :: time_interp_missing = -1e99
  real,    parameter :: epsln = 1.e-10
  real,    parameter :: sec_in_day = 86400.
  real,  allocatable :: gland_mask(:,:)        ! land mask on global domain
  real,  allocatable :: discharge2ocean_next(:,:)         ! store discharge value
  real,  allocatable :: discharge2ocean_next_c(:,:,:)     ! store discharge value
  real,  allocatable :: discharge2land_next(:,:)         ! store discharge value
  real,  allocatable :: discharge2land_next_c(:,:,:)     ! store discharge value
  real,  allocatable :: runoff_total(:,:)      ! store runoff accumulation
  real,  allocatable :: runoff_c_total(:,:,:)  ! store runoff_c accumulation
  real,  allocatable :: land_area(:,:)         ! land area on computation domain
  real,  allocatable :: wrk1(:,:)              ! on river global domain
  real,  allocatable :: wrk1_c(:,:,:)
  real,  allocatable :: wrk2(:,:)              ! on land global domain
  real,  allocatable :: wrk3(:,:)              ! on land computation domain
  real,  allocatable :: wrk4(:,:)              ! hold runoff forcing data
  integer,  allocatable :: id_infloc_c(:), id_storage_c(:), &
                           id_inflow_c(:), id_outflow_c(:), id_removal_c(:), id_disc2o(:), id_disc2l(:)
  character(len=8), allocatable:: cname(:), ifname(:), ofname(:), stname(:), rfname(:), rmname(:)
  character(len=8), allocatable:: doname(:), dlname(:)
  character(len=5), allocatable:: cunits(:), ifunits(:), ofunits(:), stunits(:), rfunits(:), rmunits(:)
  character(len=5), allocatable:: dounits(:), dlunits(:)
  integer            :: num_fast_calls 
  integer            :: slow_step = 0          ! record number of slow time step run.
  real               :: D2R

  type(horiz_interp_type), save :: Interp_land_to_river, Interp_river_to_land
  type(domain2d),          save :: domain_lnd
  type(domain2d),          save :: domain_river            
  type(river_type) ,       save :: River

  !--- these variables are for communication purpose
  integer              :: pe, root_pe, npes

  integer, allocatable :: petravel(:,:)        ! travel on each pe.
  integer, allocatable :: ncells(:)            ! number of points with each travel value

  !--- clock id variable 
  integer :: slowclock, bndslowclock, physicsclock, diagclock

contains


  !#####################################################################
  subroutine river_init( gblon, gblat, time, dt_fast, domain, land_area_in, land_mask )
    real,            intent(in) :: gblon(:)          ! lon boundaries of the grid cells
    real,            intent(in) :: gblat(:)          ! lat boundaries of the grid cells
    type(time_type), intent(in) :: time              ! current time
    type(time_type), intent(in) :: dt_fast           ! fast time step
    type(domain2d),  intent(in) :: domain            ! our domain
    real,            intent(in) :: land_area_in(:,:) 
    real,            intent(in) :: land_mask(:,:)

    integer            :: layout(2) = (/1,0/)    
    integer            :: unit, io_status, ierr, siz(4)
    integer            :: sec, day, i, j
    character(len=128) :: filename

    type(Leo_Mad_trios)   :: DHG_exp            ! downstream equation exponents
    type(Leo_Mad_trios)   :: DHG_coef           ! downstream equation coefficients
    type(Leo_Mad_trios)   :: AAS_exp            ! at-a-station equation exponents 

    D2R = PI/180.
    slowclock = mpp_clock_id('update_river_slow')
    bndslowclock = mpp_clock_id('update_river_bnd_slow')    
    physicsclock = mpp_clock_id('river phys')
    diagclock    = mpp_clock_id('river diag')
    !--- read namelist -------------------------------------------------
    unit = open_namelist_file()
    ierr = 1;
    do while (ierr /= 0)
       read  (unit, nml=river_nml, iostat=io_status, end=10)
    ierr = check_nml_error(io_status,'river_nml')
    enddo
10  call close_file (unit)

    !--- write version and namelist info to logfile --------------------
    call write_version_number(version,tagname)
    write (stdlog(), river_nml)  

    if(.not.do_rivers) return ! do nothing further if the rivers are turned off

    !--- check name list variables 

    if(diag_freq .le. 0) call mpp_error(FATAL,'river_mod: diag_freq should be a positive integer')

    pe      = mpp_pe()
    root_pe = mpp_root_pe()
    npes    = mpp_npes()

    ! set up time-related values
    River % time = time
    call get_time(dt_fast, sec, day)
    River%dt_fast = day*sec_in_day+sec

    River%dt_slow = dt_slow

    num_fast_calls = River%dt_slow/River%dt_fast
    num_species = num_phys + num_c
    if (do_age) num_species = num_species + 1
    River%num_species = num_species
    River%num_c = num_c
    River%do_age = do_age
    River%num_phys = num_phys

    if(River%dt_slow .lt. River%dt_fast) call mpp_error(FATAL, &
           'river_mod: river slow time step dt_slow should be no less than land model fast time step dt_fast')

    if ( mod(River%dt_slow,River%dt_fast) .ne. 0  ) call mpp_error(FATAL, &
           'river_mod: river slow time step dt_slow should be multiple of land model fast time step dt_fast')

    select case(trim(runoff_type))
    case ('from_land_model' )
       runoff_source = FROM_LAND_MODEL
    case ('read_on_land')
       runoff_source = READ_ON_LAND
    case ('read_on_river')
       runoff_source = READ_ON_RIVER
    case default
       call mpp_error(FATAL, 'river_mod: runoff_type = '//trim(runoff_type)//' is not a valid namelist option')
    end select


    !--- get the domain decompsition of land model
    nlon_lnd = size(gblon(:)) -1
    nlat_lnd = size(gblat(:)) -1
    call mpp_get_compute_domain(domain, isc_lnd, iec_lnd, jsc_lnd, jec_lnd)
    allocate(discharge2ocean_next(isc_lnd:iec_lnd,jsc_lnd:jec_lnd))
    allocate(discharge2ocean_next_c(isc_lnd:iec_lnd,jsc_lnd:jec_lnd,num_species))
    allocate(discharge2land_next(isc_lnd:iec_lnd,jsc_lnd:jec_lnd))
    allocate(discharge2land_next_c(isc_lnd:iec_lnd,jsc_lnd:jec_lnd,num_species))
    allocate(runoff_total(isc_lnd:iec_lnd,jsc_lnd:jec_lnd))
    allocate(runoff_c_total(isc_lnd:iec_lnd,jsc_lnd:jec_lnd,num_species))
    allocate(land_area(isc_lnd:iec_lnd,jsc_lnd:jec_lnd))
    allocate(gland_mask(nlon_lnd, nlat_lnd))
    allocate(id_infloc_c(num_species),id_storage_c(num_species))
    allocate(id_inflow_c(num_species),id_outflow_c(num_species),id_removal_c(num_species))
    allocate(id_disc2o(num_species), id_disc2l(num_species))
    allocate(cname(num_species),ifname(num_species),ofname(num_species), &
             stname(num_species),rfname(num_species),rmname(num_species))
    allocate(doname(num_species),dlname(num_species))
    allocate(cunits(num_species),ifunits(num_species),ofunits(num_species),&
             stunits(num_species),rfunits(num_species),rmunits(num_species))
    allocate(dounits(num_species), dlunits(num_species))
    allocate(source_conc_file(num_species-num_c+1:num_species))
    allocate(source_flux_file(num_species-num_c+1:num_species))
    allocate(source_conc_name(num_species-num_c+1:num_species))
    allocate(source_flux_name(num_species-num_c+1:num_species))
    allocate(          c_name(num_species-num_c+1:num_species))

    land_area  = land_area_in
    domain_lnd = domain
    runoff_total = 0.0
    runoff_c_total = 0.0
    gland_mask = land_mask
    !--- check the compatibility of land model grid with ocean grid used in river_regrid.
    call compare_ocean_grid()

    !--- read the data from the file river_src_file -- has all static river network data
    call get_river_data(gblon, gblat)

    c_name       = rt_c_name(1:num_c)
    River%t_ref  = rt_t_ref (1:num_c)
    River%vf_ref = rt_vf_ref(1:num_c)
    River%q10    = rt_q10   (1:num_c)
    River%kinv   = rt_kinv  (1:num_c)
    source_conc_file(num_species-num_c+1:num_species) = rt_source_conc_file(1:num_c)
    source_flux_file(num_species-num_c+1:num_species) = rt_source_flux_file(1:num_c)
    source_conc_name(num_species-num_c+1:num_species) = rt_source_conc_name(1:num_c)
    source_flux_name(num_species-num_c+1:num_species) = rt_source_flux_name(1:num_c)

    cname (1)='frazil'; cunits (1)='  -  '
    cname (2)='t_rivr'; cunits (2)='  K  '
    if (do_age) cname (3)='age_rv'
    if (do_age) cunits (3)='days '
    do i_species = num_species-num_c+1, num_species
      cname(i_species)=trim(c_name(i_species))
      enddo
    cunits(num_species-num_c+1:num_species)='kg/m3'
    ifunits(1)  ='m3/s '
    ifunits(2)  =' J/s '
    if (do_age) ifunits(3)  =' m3  '
    ifunits(num_species-num_c+1:num_species)='kg/s'
    ofunits(1)  ='m3/s '
    ofunits(2)  =' J/s '
    if (do_age) ofunits(3)  =' m3  '
    ofunits(num_species-num_c+1:num_species)='kg/s'
    stunits(1)  ='m3   '
    stunits(2)  =' J   '
    if (do_age) stunits(3)  ='s-m3 '
    stunits(num_species-num_c+1:num_species)='kg'
    rfunits(1)  ='m/s  '
    rfunits(2)  =' J/s '
    if (do_age) rfunits(3)  ='  m3 '
    rfunits(num_species-num_c+1:num_species)='kg/s'
    rmunits(1)  ='m3/s '
    rmunits(2)  =' J/s '
    if (do_age) rmunits(3)  ='  m3 '
    rmunits(num_species-num_c+1:num_species)='kg/s'
    dounits(1)  ='m3/s '
    dounits(2)  =' J/s '
    if (do_age) dounits(3)  =' m3  '
    dounits(num_species-num_c+1:num_species)='kg/s'
    dlunits(1)  ='m3/s '
    dlunits(2)  =' J/s '
    if (do_age) dlunits(3)  =' m3  '
    dlunits(num_species-num_c+1:num_species)='kg/s'

    do i_species = 1, num_species
      ifname(i_species)='i_'//trim(cname(i_species))
      ofname(i_species)='o_'//trim(cname(i_species))
      stname(i_species)='s_'//trim(cname(i_species))
      rfname(i_species)='r_'//trim(cname(i_species))
      rmname(i_species)='l_'//trim(cname(i_species))
      dlname(i_species)='dl'//trim(cname(i_species))
      doname(i_species)='do'//trim(cname(i_species))
      enddo

    !--- if runoff_type is not "from_land_model", need to read the axis data of the
    !--- runoff file and compare with the river grid.
    if(runoff_source .NE. FROM_LAND_MODEL) call runoff_setup(gblon, gblat)

    !--- river grid domain decomposition
    call mpp_define_layout((/1,nlon,1,nlat/),npes,layout)
    call mpp_define_domains((/1,nlon,1,nlat/),layout, domain_river)
    call mpp_get_compute_domain(domain_river,isc,iec,jsc,jec)

    !--- working array memory allocation
    allocate(wrk1(nlon,nlat), &
           wrk1_c(nlon,nlat,num_species), &
           wrk2(nlon_lnd, nlat_lnd), wrk3(isc_lnd:iec_lnd,jsc_lnd:jec_lnd))
    wrk1 = 0.0
    wrk1_c = 0.
    wrk2 = 0.0
    wrk3 = 0.0

    !--- set up the interp relation between land and river
!    call horiz_interp_new(Interp_land_to_river, gblon, gblat,                &
!                          River%lonb(isc:iec+1), River%latb(jsc:jec+1) )
    call horiz_interp_new(Interp_land_to_river, gblon, gblat,                &
                          River%lonb, River%latb )
    call horiz_interp_new(Interp_river_to_land, River%lonb, River%latb,      &
                          gblon(isc_lnd:iec_lnd+1), gblat(jsc_lnd:jec_lnd+1) )

    call sort_basins

    !--- register diag field
    call river_diag_init

    !--- read outflow mean data 
    allocate(River%outflowmean(nlon,nlat))
    call read_data(river_Omean_file, 'Omean', River%outflowmean, no_domain = .true.)
    where(River%outflowmean .le. outflowmean_min) River%outflowmean=outflowmean_min

    allocate(River%source_conc(nlon,nlat,num_species-num_c+1:num_species))
    allocate(River%source_flux(nlon,nlat,num_species-num_c+1:num_species))
    do i_species = num_species-num_c+1, num_species
      if (trim(source_conc_file(i_species)).eq.'') then
          River%source_conc(:,:,i_species)=0
          if (trim(source_conc_name(i_species)).eq.'one') River%source_conc(:,:,i_species)=1
        else if (trim(source_conc_name(i_species)).ne.'') then
          call read_data(trim(source_conc_file(i_species)), trim(source_conc_name(i_species)), &
                River%source_conc(:,:,i_species), no_domain=.true.)
        else
          River%source_conc(:,:,i_species) = 0
        endif
      if (trim(source_flux_file(i_species)).eq.'') then
          River%source_flux(:,:,i_species)=0
          if (trim(source_flux_name(i_species)).eq.'one') River%source_flux(:,:,i_species)=1
        else if (trim(source_flux_name(i_species)).ne.'') then
          call read_data(trim(source_flux_file(i_species)), &
               trim(source_flux_name(i_species)), &
                River%source_flux(:,:,i_species), no_domain=.true.)
        else
          River%source_flux(:,:,i_species) = 0
        endif
      enddo

  River%source_conc = max(River%source_conc, 0.)
  River%source_flux = max(River%source_flux, 0.)
    River%depth     = 0.
!    River%width     = 0.
!    River%vel       = 0.
    River%outflow   = 0.
    River%outflow_c = 0.
    River%inflow    = 0.
    River%inflow_c  = 0.

    !--- read restart file 
    filename = 'INPUT/river.res.nc'
    if(file_exist(trim(filename)) ) then
       call field_size(filename, 'storage', siz)
       if(siz(1) .ne. nlon .or. siz(2) .ne. nlat) then
          call mpp_error(FATAL,'river_mod: size mismatch between file INPUT/river.res.nc and '//trim(river_src_file) )
       endif

       call read_data(filename,'storage', River%storage, no_domain = .true.)
       where( .not. River%pemask)
          River%storage = 0.0
       end where

       call read_data(filename,'storage_c', River%storage_c, no_domain = .true.)
       do i_species = 1, num_species
       where( .not. River%pemask)
          River%storage_c(:,:,i_species) = 0.0
       end where
       enddo

       call read_data(filename,'discharge2ocean',  discharge2ocean_next,   domain_lnd)
       call read_data(filename,'discharge2ocean_c',discharge2ocean_next_c, domain_lnd)
       call read_data(filename,'discharge2land',   discharge2land_next,    domain_lnd)
       call read_data(filename,'discharge2land_c', discharge2land_next_c,  domain_lnd)

       if( pe == root_pe )write(stdout(),*) 'Read restart files INPUT/river.res.nc'
    else
       River%storage    = 0.0
       River%storage_c  = 0.0
       discharge2ocean_next   = 0.0
       discharge2ocean_next_c = 0.0
       discharge2land_next    = 0.0
       discharge2land_next_c  = 0.0
       if( pe == root_pe )write(stdout(),*) 'cold restart, set data to 0 '
    endif

    allocate(petravel(nlon, nlat), ncells(0:maxtravel) )
    petravel = -1
    do j = 1, nlat
       do i = 1, nlon
          if(River%pemask(i,j)) petravel(i,j) = River%travel(i,j)
       enddo
    enddo

    do i = 0, maxtravel
       ncells(i) = count(petravel == i)
    enddo

    call river_physics_init()
    call get_Leo_Mad_params(DHG_exp, DHG_coef, AAS_exp)
    River%o_exp  = 1./ (AAS_exp%on_w + AAS_exp%on_d)
    do j = 1, nlat
    do i = 1, nlon
       if ( River%celllength(i,j) > 0.0) then
          River%o_coef(i,j) = River%outflowmean(i,j) / &
               ((sinuosity*River%celllength(i,j))*DHG_coef%on_w*DHG_coef%on_d &
                 *(River%outflowmean(i,j)**(DHG_exp%on_w+DHG_exp%on_d)))**River%o_exp(i,j)
       endif
    enddo
    enddo
    River%d_exp  = AAS_exp%on_d
    River%d_coef = DHG_coef%on_d                        &
                    *(River%outflowmean**(DHG_exp%on_d-AAS_exp%on_d))

    module_is_initialized = .TRUE.

  end subroutine river_init

  !#####################################################################
  subroutine update_river ( runoff   , runoff_c,             &
                            discharge2ocean, discharge2ocean_c, &
                            discharge2land, discharge2land_c )
    real, dimension(:,:),   intent(in)  :: runoff
    real, dimension(:,:,:), intent(in)  :: runoff_c
    real, dimension(:,:),   intent(out) :: discharge2ocean
    real, dimension(:,:),   intent(out) :: discharge2land
    real, dimension(:,:,:), intent(out) :: discharge2ocean_c
    real, dimension(:,:,:), intent(out) :: discharge2land_c
    integer, save :: n = 0  ! fast time step with each slow time step

    if (.not.do_rivers) then
      discharge2ocean = 0; discharge2ocean_c = 0
      discharge2land  = 0; discharge2land_c  = 0
      return
    endif

    discharge2ocean   = discharge2ocean_next
    discharge2ocean_c = discharge2ocean_next_c
    discharge2land    = discharge2land_next
    discharge2land_c  = discharge2land_next_c

    !  increment time
    River%Time = increment_time(River%Time, River%dt_fast, 0)
    n = n + 1
    !--- accumulate runoff when read_runoff_forcing is false ---------------------
    if( runoff_source == FROM_LAND_MODEL) then
       if( n == 1) then
          runoff_total = 0.0
          runoff_c_total = 0.0
       endif
       runoff_total   = runoff_total   + runoff
       runoff_c_total = runoff_c_total + runoff_c
    endif
       if (isc_lnd.le.120.and.iec_lnd.ge.120.and.jsc_lnd.le.9.and.jec_lnd.ge.9) then
           write(6,'(2i5,1x,a16,1x,2e10.2)') n,mpp_pe(),'rnf_tot',runoff_total(120,9),gland_mask(120,9)
           write(6,'(2i5,1x,a16,1x,2e10.2)') n,mpp_pe(),'rnfc_tot',runoff_c_total(120,9,1)
           endif
    if(n == num_fast_calls) then
       call mpp_clock_begin(slowclock)
       call update_river_slow(runoff_total(:,:)  /real(num_fast_calls), &
                              runoff_c_total(:,:,:)/real(num_fast_calls)  )
       call mpp_clock_end(slowclock)       
       call mpp_clock_begin(bndslowclock)
       call update_river_bnd_slow
       call mpp_clock_end(bndslowclock)
       n = 0
    endif
    
  end subroutine update_river

  !#####################################################################
  subroutine update_river_slow(runoff, runoff_c)
    real, dimension(:,:), intent(in) :: runoff
    real, dimension(:,:,:), intent(in) :: runoff_c
    integer                             :: travelnow

    slow_step = slow_step + 1

    wrk1 = 0.0
    wrk1_c = 0.0
    select case(runoff_source)
    case(FROM_LAND_MODEL)

       wrk2 = 0.0
       wrk2(isc_lnd:iec_lnd,jsc_lnd:jec_lnd) = runoff(:,:) &
                        *gland_mask(isc_lnd:iec_lnd,jsc_lnd:jec_lnd)
       call mpp_sum(wrk2, nlon_lnd*nlat_lnd)
       call horiz_interp(Interp_land_to_river, wrk2, wrk1 )
       if (isc_lnd.le.120.and.iec_lnd.ge.120.and.jsc_lnd.le.9.and.jec_lnd.ge.9) then
           write(6,'(2i5,1x,a16,1x,2e10.2)') slow_step,mpp_pe(),'rnf_lnd,mask',runoff(120-isc_lnd+1,9-jsc_lnd+1),gland_mask(120,9)
           write(6,'(2i5,1x,a16,1x,2e10.2)') slow_step,mpp_pe(),'rnf_rvr,cellarea',wrk1(598,33),River%cellarea(598,33)
           endif

      do i_species = 1, River%num_phys
       wrk2 = 0.0
       wrk2(isc_lnd:iec_lnd,jsc_lnd:jec_lnd) = runoff_c(:,:, i_species)
       call mpp_sum(wrk2, nlon_lnd*nlat_lnd)
       call horiz_interp(Interp_land_to_river, wrk2*gland_mask, wrk1_c(:,:,i_species) )
       enddo

      if (River%do_age) then
      i_species = 3
       wrk2 = 0.0
       wrk2(isc_lnd:iec_lnd,jsc_lnd:jec_lnd) = runoff_c(:,:, i_species)
       call mpp_sum(wrk2, nlon_lnd*nlat_lnd)
       call horiz_interp(Interp_land_to_river, wrk2*gland_mask, wrk1_c(:,:,i_species) )
       endif

    case(READ_ON_LAND)
       call time_interp_external(id_runoff_forcing, River%Time, wrk2, verbose=.true.)
       where(wrk2 == time_interp_missing) wrk2 = 0.0
       call horiz_interp(Interp_land_to_river, wrk2, wrk1(isc:iec,jsc:jec), mask_in=gland_mask )
       call mpp_sum(wrk1, nlon*nlat)
    case(READ_ON_RIVER)       
       call time_interp_external(id_runoff_forcing, River%Time, wrk4, verbose=.true.)
       where(wrk4 == time_interp_missing) wrk4 = 0.0
       wrk1(is_model:ie_model,js_model:je_model) = wrk4(is_runoff:ie_runoff,js_runoff:je_runoff)
    end select 

    River%infloc   = River%cellarea*wrk1  /DENS_H2O
    River%infloc_c = 0
    do i_species = 1, River%num_phys
      River%infloc_c(:,:,i_species) = &
           River%cellarea*wrk1_c(:,:,i_species)/DENS_H2O
      enddo
    if (River%do_age) then
        i_species = 3
        River%infloc_c(:,:,i_species) = &
           River%cellarea*wrk1_c(:,:,i_species)/DENS_H2O
      endif
    do i_species = num_species-River%num_c+1, num_species  ! create mass flux inputs from c data
      where (River%cellarea.gt.0.)  &
        River%infloc_c(:,:,i_species) = &
           River%infloc*River%source_conc(:,:,i_species)  &
                      + River%source_flux(:,:,i_species)
      enddo
    River%inflow   = 0
    River%inflow_c = 0

    travelnow = maxtravel
    do travelnow = maxtravel, 1, -1
       if( ncells(travelnow).EQ.0 )cycle
       call mpp_clock_begin(physicsclock)
       call river_physics_step (River, petravel, travelnow )
       call mpp_clock_end(physicsclock)
       enddo

call mpp_sum(River%inflow,  nlon*nlat            )
call mpp_sum(River%inflow_c,nlon*nlat*num_species)

    River%disw2o = 0.
    where (River%tocell.eq.0 .and. River%landfrac.lt.1.)
              River%disw2o = River%inflow
      endwhere
    where (River%landfrac.lt.1. .and. River%tocell.eq.0)
              River%disw2o = River%disw2o + River%infloc
      endwhere
    where (River%landfrac.eq.0. .and. River%tocell.ne.0)
              River%disw2o = River%disw2o + River%infloc
      endwhere
    River%disc2o = 0.
    do i_species = 1, num_species
      where (River%tocell.eq.0 .and. River%landfrac.lt.1.)
              River%disc2o(:,:,i_species) = &
               + River%inflow_c(:,:,i_species)
        endwhere
      where (River%landfrac.lt.1. .and. River%tocell.eq.0)
              River%disc2o(:,:,i_species) = &
                 River%disc2o(:,:,i_species) &
                 + River%infloc_c(:,:,i_species)
        endwhere
      where (River%landfrac.eq.0. .and. River%tocell.ne.0)
              River%disc2o(:,:,i_species) = &
                 River%disc2o(:,:,i_species) &
                 + River%infloc_c(:,:,i_species)
        endwhere
      enddo

    River%disw2l = 0.
    where (River%tocell.eq.0 .and. River%landfrac.ge.1.)
              River%disw2l = River%infloc          + River%inflow
      endwhere
    River%disc2l = 0.
    do i_species = 1, num_species
      where (River%tocell.eq.0 .and. River%landfrac.ge.1.)
              River%disc2l(:,:,i_species) = &
                 River%infloc_c(:,:,i_species) &
               + River%inflow_c(:,:,i_species)
        endwhere
      enddo

    call mpp_clock_begin(diagclock)
    if(mod(slow_step, diag_freq) == 0)  call river_diag()
    call mpp_clock_end(diagclock)

  end subroutine update_river_slow

  !#####################################################################

    subroutine update_river_bnd_slow

!!    wrk1 = 0.
!!    where (River%landfrac.lt.1.)&
        wrk1 = (River%disw2o - River%disc2o(:,:,1)) &
                                / River%cellarea
    call horiz_interp(Interp_river_to_land, wrk1, wrk3)
    discharge2ocean_next = DENS_H2O*wrk3

!!    wrk1 = 0.
!!    where (River%landfrac.ge.1.)&
        wrk1 = (River%disw2l - River%disc2l(:,:,1)) &
                                / River%cellarea
    call horiz_interp(Interp_river_to_land, wrk1, wrk3)
    discharge2land_next = DENS_H2O*wrk3

    do i_species = 1, num_species
      wrk1 = 0
      where (River%landfrac.lt.1) &
        wrk1 = River%disc2o(:,:,i_species) / River%cellarea
      call horiz_interp(Interp_river_to_land, wrk1, wrk3)
      discharge2ocean_next_c(:,:,i_species) = DENS_H2O*wrk3
      wrk1 = 0
      where (River%landfrac.ge.1.) &
        wrk1 = River%disc2l(:,:,i_species) / River%cellarea
      call horiz_interp(Interp_river_to_land, wrk1, wrk3)
      discharge2land_next_c(:,:,i_species) = DENS_H2O*wrk3
      enddo

    end subroutine update_river_bnd_slow

  !#####################################################################

  subroutine river_end
    character(len=128)               :: filename
    real,               allocatable  :: tmp(:,:),tmp3(:,:,:)
    !--- write to restart file

    if(.not.do_rivers) return ! do nothing further if rivers are turned off

    allocate(tmp(nlon, nlat) )
    allocate(tmp3(nlon, nlat, num_species) )

    filename = 'RESTART/river.res.nc'

    tmp = River%storage
    call mpp_sum(tmp,nlon*nlat)
    call write_data(filename,'storage', tmp(isc:iec,jsc:jec), domain_river)

    tmp3 = River%storage_c(:,:,:)
    call mpp_sum(tmp3,nlon*nlat*num_species)       ! ????????
    call write_data(filename,'storage_c', tmp3(isc:iec,jsc:jec,:), domain_river)

    !--- write out discharge data
    call write_data(filename,'discharge2ocean'  ,discharge2ocean_next  (isc_lnd:iec_lnd,jsc_lnd:jec_lnd),   domain_lnd)
    call write_data(filename,'discharge2ocean_c',discharge2ocean_next_c(isc_lnd:iec_lnd,jsc_lnd:jec_lnd,:), domain_lnd)
    call write_data(filename,'discharge2land'   ,discharge2land_next   (isc_lnd:iec_lnd,jsc_lnd:jec_lnd),   domain_lnd)
    call write_data(filename,'discharge2land_c' ,discharge2land_next_c (isc_lnd:iec_lnd,jsc_lnd:jec_lnd,:), domain_lnd)

    !--- release memory
    deallocate(tmp, discharge2ocean_next, discharge2ocean_next_c,&
                    discharge2land_next, discharge2land_next_c,&
                    wrk1_c, &
                    wrk1, wrk2, wrk3, runoff_total,  &
                    runoff_c_total)

    deallocate( River%lon, River%lonb)
    deallocate( River%lat, River%latb)
    deallocate(River%cellarea ,     River%basinid        )
    deallocate(River%landfrac )
    deallocate(River%tocell )
    deallocate(River%travel   ,     River%pemask       )
    deallocate(River%outflow  )
    deallocate(River%inflow  )
    deallocate(River%storage        )
    deallocate(River%disw2o        )
    deallocate(River%disc2o        )
    deallocate(River%disw2l        )
!    deallocate(River%diss2l        )
    deallocate(River%disc2l        )
    deallocate(River%infloc   )
    deallocate(River%celllength    )
    deallocate(River%gmask        )
    deallocate(River%So        )
    deallocate(River%depth     )
!    deallocate(River%width     )
!    deallocate(River%vel       )
    deallocate(River%infloc_c ,     River%storage_c    )
    deallocate(River%inflow_c, River%outflow_c )
    deallocate(River%removal_c )
    deallocate(River%vf_ref,River%t_ref,River%q10,River%kinv)
    deallocate(River%d_exp,River%d_coef,River%o_exp,River%o_coef)
    call horiz_interp_del(Interp_land_to_river)
    call horiz_interp_del(Interp_river_to_land)

    module_is_initialized = .FALSE.

  end subroutine river_end

  !#####################################################################
  ! compare ocean grid in grid_spec.nc and river networking. 
  subroutine compare_ocean_grid( )

!!$    integer                           :: ni, nj, i, j, siz(4)
!!$    real, dimension(:,:), allocatable :: lon1, lat1, mask1
!!$    real, dimension(:,:), allocatable :: lon2, lat2, mask2
!!$    character(len=128)                :: grid_file = "INPUT/grid_spec.nc"
!!$
!!$    !--- check if the model land grid is the same as the land grid 
!!$    !--- used to do river_regrid.
!!$
!!$    !--- nullify land domain to aviod size mismatch.
!!$    call field_size(river_src_file, 'mask_ocean', siz)
!!$    ni = siz(1)
!!$    nj = siz(2)
!!$    call field_size(grid_file, "wet", siz)
!!$    if(siz(1) .NE. ni .OR. siz(2) .NE. nj) call mpp_error(FATAL, &
!!$       "river_mod: size mismatch of ocean grid between INPUT/grid_spec.nc and the grid file used in river_regrid.")
!!$
!!$    allocate(lon1(ni,nj), lat1(ni,nj), mask1(ni,nj) )
!!$    allocate(lon2(ni,nj), lat2(ni,nj), mask2(ni,nj) )
!!$    call read_data(river_src_file, 'lon_ocean', lon1, no_domain = .true.)
!!$    call read_data(river_src_file, 'lat_ocean', lat1, no_domain = .true.)
!!$    call read_data(river_src_file, 'mask_ocean', mask1, no_domain = .true.)  
!!$
!!$    call read_data(grid_file, 'wet', mask2, no_domain = .true.)   
!!$    if(field_exist(grid_file, "geolon_t")) then
!!$       call read_data(grid_file, 'geolon_t', lon2, no_domain = .true.)
!!$       call read_data(grid_file, 'geolat_t', lat2, no_domain = .true.)       
!!$    else
!!$       call read_data(grid_file, 'x_T', lon2, no_domain = .true.)
!!$       call read_data(grid_file, 'y_T', lat2, no_domain = .true.)   
!!$    end if
!!$
!!$    do j = 1, nj
!!$       do i = 1, ni
!!$          if(abs(lon1(i,j)-lon2(i,j)) > epsln) call mpp_error(FATAL, &
!!$               "river_mod: ocean longitude mismatch between INPUT/grid_spec.nc and the grid file used in river_regrid.")
!!$          if(abs(lat1(i,j)-lat2(i,j)) > epsln) call mpp_error(FATAL, &
!!$               "river_mod: ocean latitude mismatch between INPUT/grid_spec.nc and the grid file used in river_regrid.")
!!$          if(abs(mask1(i,j)-mask2(i,j)) > epsln) call mpp_error(FATAL, &
!!$               "river_mod: ocean mask mismatch between INPUT/grid_spec.nc and the grid file used in river_regrid.")
!!$       end do
!!$    end do
!!$
!!$    deallocate(lon1, lat1, mask1, lon2, lat2, mask2)

  end subroutine compare_ocean_grid

  !#####################################################################
  subroutine get_river_data(lonb_lnd, latb_lnd)
    real, dimension(:),    intent(in) :: lonb_lnd, latb_lnd

    integer                           :: ni, nj, i, j, siz(4)
    integer                           :: istart, iend, jstart, jend
    real                              :: min_lon, min_lat, max_lon, max_lat
    real, dimension(:,:), allocatable :: tmp(:,:)
    real, dimension(:), allocatable   :: xt, yt, xb, yb
    integer                           :: basin, ii, jj
    real xxx

    call field_size(river_src_file, 'basin', siz)
    ni = siz(1)
    nj = siz(2)
    allocate( xt(ni), yt(nj), xb(ni+1), yb(nj+1) )

    call read_data(river_src_file, 'lon', xt, no_domain = .true.)
    call read_data(river_src_file, 'lat', yt, no_domain = .true.) 
    xb(1) = xt(1) - 0.5*(xt(2)-xt(1))
    if (abs(xb(1)) < epsln) xb(1) = 0.0
    do i = 2, ni
       xb(i) = xt(i-1) + 0.5*(xt(i)-xt(i-1))
    enddo
    xb(ni+1) = xt(ni) + 0.5*(xt(ni)-xt(ni-1))
    if(abs(xb(ni+1) - 360.) < epsln) xb(ni+1) = 360.

    yb(1) = yt(1) - 0.5*(yt(2)-yt(1))
    do j = 2, nj
       yb(j) = yt(j-1) + 0.5*(yt(j)-yt(j-1))
    enddo
    yb(nj+1) = yt(nj) + 0.5*(yt(nj)-yt(nj-1))

    !--- transform to radians, since land model grid use radians.
    xt = xt * D2R
    yt = yt * D2R
    xb = xb * D2R
    yb = yb * D2R

    !--- match the river grid with land grid.
    min_lon = minval(lonb_lnd)
    max_lon = maxval(lonb_lnd)
    min_lat = minval(latb_lnd)
    max_lat = maxval(latb_lnd)
    if(min_lon < xb(1) .OR. min_lat > xb(nlon+1) .OR. min_lat < yb(1) .OR. min_lat > yb(nlon+1) ) then
       call mpp_error(FATAL, "river_mod: land grid should be inside or river network")
    end if
    istart = nearest_index(min_lon, xb)
    iend   = nearest_index(max_lon, xb)
    jstart = nearest_index(min_lat, yb)
    jend   = nearest_index(max_lat, yb)
    if(xb(istart) > min_lon) istart = istart - 1
    if(xb(iend)   < max_lon) iend   = iend   + 1
    if(yb(jstart) > min_lat) jstart = jstart - 1
    if(yb(jend)   < max_lat) jend   = jend   + 1
    iend = iend - 1
    jend = jend - 1

    nlon = iend - istart + 1
    nlat = jend - jstart + 1
    River%nlon = nlon
    River%nlat = nlat

    allocate( River%lon(nlon), River%lonb(nlon+1))
    allocate( River%lat(nlat), River%latb(nlat+1))
    allocate(River%cellarea (nlon,    nlat),     River%basinid      (nlon,nlat)  )
    allocate(River%landfrac (nlon,    nlat) )
    allocate(River%tocell (nlon,    nlat) )
    allocate(River%travel   (nlon,    nlat),     River%pemask       (nlon,nlat))
    allocate(River%inflow   (nlon,    nlat)      )
    allocate(River%outflow  (0:nlon+1,0:nlat+1)  )
    allocate(River%storage  (nlon,    nlat)      )
    allocate(River%disw2o  (nlon,    nlat)      )
    allocate(River%disc2o  (nlon,    nlat, num_species)      )
    allocate(River%disw2l  (nlon,    nlat)      )
 !   allocate(River%diss2l  (nlon,    nlat)      )
    allocate(River%disc2l  (nlon,    nlat, num_species)      )
    allocate(River%infloc   (nlon,    nlat))
    allocate(River%celllength   (nlon,nlat) )
    allocate(River%gmask    (nlon,    nlat),     tmp                (ni, nj)    )
    allocate(River%So       (nlon,    nlat) )
    allocate(River%depth    (nlon,    nlat) )
!    allocate(River%width    (nlon,    nlat) )
!    allocate(River%vel      (nlon,    nlat) )
    allocate(River%infloc_c (nlon,    nlat, num_species),     River%storage_c   (nlon,    nlat, num_species) )
    allocate(River%outflow_c   (0:nlon+1,0:nlat+1, num_species) )
    allocate(River%removal_c (nlon,    nlat, num_species), River%inflow_c (nlon,    nlat, num_species) )
    allocate(River%t_ref(4:num_species),River%vf_ref(4:num_species))
    allocate(River%q10  (4:num_species),River%kinv  (4:num_species))
    allocate(River%d_exp(nlon,nlat))
    allocate(River%d_coef(nlon,nlat))
    allocate(River%o_exp(nlon,nlat))
    allocate(River%o_coef(nlon,nlat))

    River%lon(:)  = xt(istart:iend)
    River%lat(:)  = yt(jstart:jend)
    River%lonb(:) = xb(istart:iend+1)
    River%latb(:) = yb(jstart:jend+1)
    River%infloc = 0.0
    River%infloc_c = 0.0
    River%storage     = 0.0
    River%storage_c   = 0.0
    River%removal_c   = 0.0

    !--- read the data from the source file
    call read_data(river_src_file, 'tocell', tmp, no_domain = .true.) 
    River%tocell(:,:) = tmp(istart:iend, jstart:jend)
    where (River%tocell(:,:).eq.  4) River%tocell(:,:)=3
    where (River%tocell(:,:).eq.  8) River%tocell(:,:)=4
    where (River%tocell(:,:).eq. 16) River%tocell(:,:)=5
    where (River%tocell(:,:).eq. 32) River%tocell(:,:)=6
    where (River%tocell(:,:).eq. 64) River%tocell(:,:)=7
    where (River%tocell(:,:).eq.128) River%tocell(:,:)=8
    call read_data(river_src_file, 'basin', tmp, no_domain = .true.)
    River%basinid(:,:) = tmp(istart:iend, jstart:jend) 
    call read_data(river_src_file, 'travel', tmp, no_domain = .true.) 
    River%travel(:,:) = tmp(istart:iend, jstart:jend)
    call read_data(river_src_file, 'celllength', tmp, no_domain = .true.) 
    River%celllength(:,:) = tmp(istart:iend, jstart:jend)
    call read_data(river_src_file, 'cellarea', tmp, no_domain = .true.) 
    River%cellarea(:,:) = tmp(istart:iend, jstart:jend)
    call read_data(river_src_file, 'mask', tmp, no_domain = .true.) 
    River%landfrac(:,:) = tmp(istart:iend, jstart:jend)
    call read_data(river_src_file, 'So', tmp, no_domain = .true.) 
    River%So(:,:) = tmp(istart:iend, jstart:jend)

    where (River%So .LT. 0.0) River%So = Somin

    !--- search each point on the boundaries to remove truncated river. Only for regional model
    !--- we may need to change the following part in the future.
    if( ni .ne. nlon .or. nj .ne. nlat) then
       write(stdout(),*)' Warning from river_mod: remove truncated river'
       do j = 1, nlat
          do i = 1, nlon
             if(i==1 .or. i==nlon .or. j==1 .or. j==nlat) then
                if(River%travel(i,j) > 0) then
                   basin = River%basinid(i,j) 
                   do jj = 1, nlat
                      do ii = 1, nlon
                         if(River%basinid(ii,jj) == basin) then
                            River%tocell(ii,jj) = missing
                            River%basinid(ii,jj) = missing
                            River%travel(ii,jj) = missing
                         endif
                      enddo
                   enddo
                endif
             endif
          enddo
       enddo
    endif

    !--- define global mask
    where (River%basinid>0)
       River%gmask = .true.
    elsewhere
       River%gmask = .false.
    end where

  end subroutine get_river_data

  !#####################################################################
  !--- the grid in the runoff forcing file does not need to be exact 
  !--- match the river model grid, it will read the runoff data onto 
  !--- the intersection of region of river grid and runoff data grid.
  !--- set runoff to 0 at all other point. This routine will find 
  !--- the intersection region index in river grid, denoted by
  !--- is_model, ie_model, js_model, je_model and the intersection
  !--- region index in runoff forcing data, denoted by 
  !--- is_forcing, ie_forcing, js_forcing, je_forcing.
  subroutine runoff_setup(lonb_lnd, latb_lnd)
    real, dimension(:),    intent(in) :: lonb_lnd, latb_lnd

    integer                      :: ni, nj, i, j, ii, jj
    integer                      ::len1, unit, ndim, nvar, natt, ntime
    character(len=64)            :: runoff_field, name
    logical                      :: found_runoff
    character(len=1)             :: cart
    real                         :: lon_min, lat_min, lon_max, lat_max
    real                         :: lon_lnd, lat_lnd
    real,            allocatable :: lon(:), lat(:)
    type(fieldtype), allocatable :: fields(:)
    type(axistype),  allocatable :: axes(:)

    runoff_field = "runoff"

    if(.not. file_exist(runoff_file)) call mpp_error(FATAL,  &
         "river_mode: nml read_runoff_forcing is set to true, but file " &
         //trim(runoff_file)//" does not exist")

    !--- read the data location.
    call mpp_open(unit, trim(runoff_file),&
         action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)

    allocate(fields(nvar))
    call mpp_get_fields(unit, fields)

    ni=0; nj=0
    found_runoff = .false.
    do i=1,nvar
       call mpp_get_atts(fields(i),name=name)
       if (lowercase(trim(runoff_field)) == lowercase(trim(name))) then
          found_runoff = .true.
          call mpp_get_atts(fields(i),ndim=ndim)
          allocate(axes(ndim))
          call mpp_get_atts(fields(i),axes=axes)
          do j=1,ndim
             call mpp_get_atts(axes(j),len=len1)
             call get_axis_cart(axes(j),cart)
             select case (cart)
             case ('X')
                ni = len1
                allocate(lon(ni))
                call mpp_get_axis_data(axes(j),lon)
             case('Y')
                nj = len1
                allocate(lat(nj))
                call mpp_get_axis_data(axes(j),lat)
             end select
          enddo
          exit
       endif
    end do
 
    if(.not. found_runoff) call mpp_error(FATAL, "river_mode: field "//trim(runoff_field)// &
             " does not exist in file "//trim(runoff_file) )

    if(ni==0 .OR. nj ==0) call mpp_error(FATAL,'river_mod: the cartesian attribute of axis of field ' &
       //trim(runoff_field)//' in file '//trim(runoff_file)//' is not correct' )  

    allocate(wrk4(ni,nj))

    lon = lon*D2R
    lat = lat*D2R

    call time_interp_external_init()
    id_runoff_forcing = init_external_field(runoff_file, runoff_field )

    !--- when runoff_source is READ_ON_LAND, we require the runoff data on the same grid as land model.
    select case(runoff_source)
    case(READ_ON_LAND)
       if( (ni .NE. nlon_lnd) .OR. (nj .NE. nlat_lnd) )  call mpp_error(FATAL, &
           'river_mod: size mismatch between land grid and runoff data')

       do i = 1, nlon_lnd
          lon_lnd = 0.5*(lonb_lnd(i) + lonb_lnd(i+1))
          if(abs(lon(i) - lon_lnd) >epsln )  call mpp_error(FATAL, &
           'river_mod: longitude mismatch between land grid and runoff data')
       enddo

       do j = 1, nlat_lnd
          lat_lnd = 0.5*(latb_lnd(j) + latb_lnd(j+1))
          if(abs(lat(j) - lat_lnd) >epsln )  call mpp_error(FATAL, &
           'river_mod: latitude mismatch between land grid and runoff data')
       enddo

    case(READ_ON_RIVER)
       !--- we don't require the forcing data grid is the same with the river grid, 
       !--- but we require if a runoff grid is contained in the river region, it should 
       !--- match with one river grid, also vice verse. we don't consider about cyclic 
       !--- condition right now.
       lon_min = max(lon(1), River%lon(1))
       lon_max = min(lon(ni), River%lon(River%nlon))
       lat_min = max(lat(1), River%lat(1))
       lat_max = min(lat(nj), River%lat(River%nlat))

       is_runoff = get_index(lon, lon_min)
       ie_runoff = get_index(lon, lon_max)
       js_runoff = get_index(lat, lat_min)
       je_runoff = get_index(lat, lat_max)
       is_model  = get_index(River%lon, lon_min)
       ie_model  = get_index(River%lon, lon_max)
       js_model  = get_index(River%lat, lat_min)
       je_model  = get_index(River%lat, lat_max)

       if(is_runoff == 0 .OR. ie_runoff == 0 .OR. js_runoff == 0 .OR. je_runoff == 0 .OR. &
            is_model  == 0 .OR. ie_model  == 0 .OR. js_model  == 0 .OR. je_model == 0 ) then
          write(stdout(),*)'lon_min = ', lon_min, 'lon_max = ', lon_max, 'lat_min = ', lat_min, 'lat_max = ', lat_max
          write(stdout(),*)'runoff index is ', is_runoff, ie_runoff, js_runoff, je_runoff
          write(stdout(),*)'model index is ', is_model, ie_model, js_model, je_model
          write(stdout(),*)'runoff lon = ', lon
          write(stdout(),*)'river lon = ', River%lon
          write(stdout(),*)'runoff lat = ', lat
          write(stdout(),*)'river lat = ', River%lat
          call mpp_error(FATAL,"river_model: mismatch between runoff data grid and river grid 1")
       endif


       if( (ie_runoff-is_runoff .NE. ie_model-is_model) .OR. (je_runoff-js_runoff .NE. je_model-js_model) ) &
            call mpp_error(FATAL,"river_model: mismatch between runoff data grid and river grid 2")

       !--- make sure the overlapping region is match
       do i = is_runoff, ie_runoff
          ii = i - is_runoff + is_model
          if(abs(lon(i) - River%lon(ii)) > epsln) &
               call mpp_error(FATAL,"river_model: mismatch between runoff data grid and river grid 3")
       enddo

       do j = js_runoff, je_runoff
          jj = j - js_runoff + js_model
          if(abs(lat(j) - River%lat(jj)) > epsln) &
               call mpp_error(FATAL,"river_model: mismatch between runoff data grid and river grid 4")
       enddo

    end select

    contains

    function get_index(array, value)
       real, dimension(:), intent(in) :: array
       real,               intent(in) :: value 
       integer                        :: get_index

       get_index = 0
       do i = 1, size(array)
          if(abs(array(i) - value) < epsln ) then
             get_index = i
             exit
          endif
       enddo       

       return    

    end function get_index

  end subroutine runoff_setup

  !#####################################################################

  subroutine sort_basins

    !   </PUBLICROUTINE>

    ! ---- local vars ----------------------------------------------------------
    integer                            :: i, j, k          ! arbitrary indices
    integer, allocatable, dimension(:) :: numcells_per_pe  ! number of cells per PE
    integer, allocatable, dimension(:) :: pelist           ! pelist
    integer :: maxbasins
    integer, allocatable :: numcells(:), basinpe(:)

    !total number of basins
    maxbasins = maxval(River%basinid)
    !get number of cells per basin
    allocate( numcells(maxbasins) )
    numcells(:) = 0
    do j = 1, nlat
       do i = 1, nlon
          if(River%basinid(i,j)>0) numcells(River%basinid(i,j)) = numcells(River%basinid(i,j)) + 1
       end do
    end do

    !assign to PEs
    allocate( numcells_per_pe(npes) )
    allocate( pelist(npes) )
    call mpp_get_current_pelist(pelist)

    numcells_per_pe(:) = 0
    allocate( basinpe(maxbasins) )
    do i = 1,maxbasins
       !min/maxloc return array of length 1: minval is used to demote that to scalar
       j = minval(minloc(numcells_per_pe)) !PE containing the least cells
       k = minval(maxloc(numcells))       !basin containing the most cells
       basinpe(k) = pelist(j)     !basin k goes to the PE pelist(j)
       numcells_per_pe(j) = numcells_per_pe(j) + numcells(k)
       numcells(k) = -1 !so maxloc(numcells) returns next highest on next pass
    end do

    !set mask to TRUE if basin is computed on this PE
    River%pemask = .false.
    do j = 1, nlat
       do i = 1, nlon
          if(River%basinid(i,j)>0) River%pemask(i,j) = basinpe(River%basinid(i,j)).EQ.pe
       end do
    end do

    !find domain limits at each travel distance
    maxtravel = maxval(River%travel,mask=River%pemask)
    allocate( River%is(0:maxtravel) )
    allocate( River%ie(0:maxtravel) )
    allocate( River%js(0:maxtravel) )
    allocate( River%je(0:maxtravel) )
    do k = 0,maxtravel
       River%is(k) =  huge(1)
       River%ie(k) = -huge(1)
       River%js(k) =  huge(1)
       River%je(k) = -huge(1)
       do j = 1,size(River%travel,2)
          do i = 1,size(River%travel,1)
             if( River%travel(i,j).EQ.k .AND. River%pemask(i,j) )then
                River%is(k) = min( i, River%is(k) )
                River%ie(k) = max( i, River%ie(k) )
                River%js(k) = min( j, River%js(k) )
                River%je(k) = max( j, River%je(k) )
             end if
          end do
       end do
       if(debug_river) then
          write( 6,'(a,6i4,f8.5)' )'PE, travel, is, ie, js, je, cell density=', pe, k, &
               River%is(k), River%ie(k), River%js(k), River%je(k), count(River%travel.EQ.k .AND. River%pemask) &
               /real( (River%ie(k)-River%is(k)+1)*(River%je(k)-River%js(k)+1) )
       endif
    end do

    write(stdout(),*) 'maxbasins=', maxbasins

    do i = 1,npes
       write(stdout(),*) 'pe, numcells=', pelist(i), numcells_per_pe(i)
    end do

end subroutine sort_basins

  !#####################################################################

  subroutine river_diag_init
    integer                       :: id_lon  ! ID of land longitude (X) diag axis
    integer                       :: id_lat  ! ID of land latitude (Y) diag axis
    character(len=11)             :: mod_name = 'river_model'
    real, dimension(isc:iec,jsc:jec) :: tmp
    logical                          :: sent

    !--- diag axis initialization 
    id_lon  = diag_axis_init ( 'riv_lon', River%lon*RADIAN, 'degrees_E', 'X',  &
         'longitude', set_name='land', Domain2=domain_river)
    id_lat = diag_axis_init ( 'riv_lat', River%lat*RADIAN, 'degrees_N', 'Y', &
         'latitude', set_name='land', Domain2=domain_river)

    ! regular diagnostic fields
    id_storage   = register_diag_field ( mod_name, 'storage', (/id_lon, id_lat/), &
         River%Time, 'storage', 'm3', missing_value=missing )
    id_inflow   = register_diag_field ( mod_name, 'inflow', (/id_lon, id_lat/), &
         River%Time, 'inflow', 'm3/s', missing_value=missing )
    id_outflow   = register_diag_field ( mod_name, 'outflow', (/id_lon, id_lat/), &
         River%Time, 'outflow', 'm3/s', missing_value=missing )
    id_infloc    = register_diag_field ( mod_name, 'infloc', (/id_lon, id_lat/), &
         River%Time, 'infloc', 'm3/s', missing_value=missing )
    id_disw2o   = register_diag_field ( mod_name, 'disw2o', (/id_lon, id_lat/), &
         River%Time, 'disw2o', 'm3/s', missing_value=missing )
!    id_diss2o   = register_diag_field ( mod_name, 'diss2o', (/id_lon, id_lat/), &
!         River%Time, 'diss2o', 'm3/s', missing_value=missing )
    id_disw2l   = register_diag_field ( mod_name, 'disw2l', (/id_lon, id_lat/), &
         River%Time, 'disw2l', 'm3/s', missing_value=missing )
!    id_diss2l   = register_diag_field ( mod_name, 'diss2l', (/id_lon, id_lat/), &
!         River%Time, 'diss2l', 'm3/s', missing_value=missing )
    do i_species = 1, num_species
      id_inflow_c(i_species) = register_diag_field ( mod_name, ifname(i_species),     &
         (/id_lon, id_lat/), River%Time, ifname(i_species), ifunits(i_species),    &
         missing_value=missing )
      id_outflow_c(i_species) = register_diag_field ( mod_name, ofname(i_species),     &
         (/id_lon, id_lat/), River%Time, ofname(i_species), ofunits(i_species),    &
         missing_value=missing )
      id_storage_c(i_species) = register_diag_field ( mod_name, stname(i_species),     &
         (/id_lon, id_lat/), River%Time, stname(i_species), stunits(i_species),    &
         missing_value=missing )
      id_infloc_c(i_species) = register_diag_field ( mod_name, rfname(i_species),     &
         (/id_lon, id_lat/), River%Time, rfname(i_species), rfunits(i_species),    &
         missing_value=missing )
      id_removal_c(i_species) = register_diag_field ( mod_name, rmname(i_species),     &
         (/id_lon, id_lat/), River%Time, rmname(i_species), rmunits(i_species),    &
         missing_value=missing )
      id_disc2o(i_species)   = register_diag_field ( mod_name, doname(i_species), &
         (/id_lon, id_lat/), River%Time, doname(i_species), dounits(i_species), &
         missing_value=missing )
      id_disc2l(i_species)   = register_diag_field ( mod_name, dlname(i_species), &
         (/id_lon, id_lat/), River%Time, dlname(i_species), dlunits(i_species), &
         missing_value=missing )
      enddo
    id_depth     = register_diag_field ( mod_name, 'depth', (/id_lon, id_lat/), &
         River%Time, 'depth', 'm', missing_value=missing )
!    id_width     = register_diag_field ( mod_name, 'width', (/id_lon, id_lat/), &
!         River%Time, 'width', 'm', missing_value=missing )
!    id_vel       = register_diag_field ( mod_name, 'vel', (/id_lon, id_lat/), &
!         River%Time, 'velocity', 'm/s', missing_value=missing )

    ! static fields
    id_dx = register_static_field ( mod_name, 'dx', (/id_lon, id_lat/), &
         'cell delta x', 'm', missing_value=missing )
    id_basin = register_static_field ( mod_name, 'basin', (/id_lon, id_lat/), &
         'river basin id', 'none', missing_value=missing )
    id_So = register_static_field ( mod_name, 'So', (/id_lon, id_lat/), &
         'Slope', 'none', missing_value=missing )
    id_travel = register_static_field ( mod_name, 'travel', (/id_lon, id_lat/), &
         'cells left to travel before reaching ocean', 'none', missing_value=missing )
    id_tocell = register_static_field ( mod_name, 'tocell', (/id_lon, id_lat/), &
         'sum of directions from upstream cells', 'none', missing_value=missing )

    if (id_dx>0) then
       sent=send_data(id_dx, River%celllength(isc:iec,jsc:jec), River%Time, &
                          mask=River%gmask(isc:iec,jsc:jec) )
    end if

    if (id_basin>0) then
       tmp = River%basinid(isc:iec,jsc:jec)
       sent=send_data(id_basin, tmp, River%Time, &
                          mask=River%gmask(isc:iec,jsc:jec) )
    end if

    if (id_So>0) then
       sent=send_data(id_So, River%So(isc:iec,jsc:jec), River%Time, &
                          mask=River%gmask(isc:iec,jsc:jec) )
    end if

    if (id_travel>0) then
       tmp = River%travel(isc:iec,jsc:jec)
       sent=send_data(id_travel, tmp, River%Time, &
                          mask=River%gmask(isc:iec,jsc:jec) )
    end if

    if (id_tocell>0) then
       tmp = River%tocell(isc:iec,jsc:jec)
       sent=send_data(id_tocell, tmp, River%Time, &
                          mask=River%gmask(isc:iec,jsc:jec) )
    end if

  end subroutine river_diag_init

  !#####################################################################

  subroutine river_diag
    logical                       :: used   ! logical for send_data
    real, dimension(nlon, nlat)   :: tmp

    if (id_inflow > 0) then
       tmp = River%inflow(1:nlon,1:nlat)
 !      call mpp_sum(tmp,nlon*nlat)
       used = send_data (id_inflow, tmp(isc:iec,jsc:jec), River%Time )
    endif

    if (id_outflow > 0) then
       tmp = River%outflow(1:nlon,1:nlat)
       call mpp_sum(tmp,nlon*nlat)
       used = send_data (id_outflow, tmp(isc:iec,jsc:jec), River%Time) !, &
                       !   mask=River%gmask(isc:iec,jsc:jec) )
    endif

    if (id_storage > 0) then 
       tmp = River%storage
       call mpp_sum(tmp,nlon*nlat)
       used = send_data (id_storage, tmp(isc:iec,jsc:jec), River%Time) !, &
                      !    mask=River%gmask(isc:iec,jsc:jec) )
    endif

    if (id_infloc > 0) then
       tmp = River%infloc
 !      call mpp_sum(tmp,nlon*nlat)
       used = send_data (id_infloc, tmp(isc:iec,jsc:jec), River%Time )
    endif

    if (id_disw2o > 0) then
       tmp = River%disw2o
 !      call mpp_sum(tmp,nlon*nlat)
       used = send_data (id_disw2o, tmp(isc:iec,jsc:jec), River%Time)
    endif

!    if (id_diss2o > 0) then
!       tmp = River%disc2o(:,:,1)
! !      call mpp_sum(tmp,nlon*nlat)
!       used = send_data (id_diss2o, tmp(isc:iec,jsc:jec), River%Time)
!    endif

    if (id_disw2l > 0) then
       tmp = River%disw2l
 !      call mpp_sum(tmp,nlon*nlat)
       used = send_data (id_disw2l, tmp(isc:iec,jsc:jec), River%Time)
    endif

!    if (id_diss2l > 0) then
!!       tmp = River%diss2l
!       tmp = River%disc2l(:,:,1)
! !      call mpp_sum(tmp,nlon*nlat)
!       used = send_data (id_diss2l, tmp(isc:iec,jsc:jec), River%Time)
!    endif

    do i_species = 1, num_species
      if (id_outflow_c(i_species) > 0) then
         tmp = River%outflow_c(1:nlon,1:nlat,i_species)
         call mpp_sum(tmp,nlon*nlat)
         used = send_data (id_outflow_c(i_species), tmp(isc:iec,jsc:jec), River%Time) !, &
                        !  mask=River%gmask(isc:iec,jsc:jec) )
        endif
      if (id_inflow_c(i_species) > 0) then
         tmp = River%inflow_c(1:nlon,1:nlat,i_species)
 !        call mpp_sum(tmp,nlon*nlat)
         used = send_data (id_inflow_c(i_species), tmp(isc:iec,jsc:jec), River%Time) !, &
                        !  mask=River%gmask(isc:iec,jsc:jec) )
        endif
      if (id_storage_c(i_species) > 0) then
         tmp = River%storage_c(:,:,i_species)
         call mpp_sum(tmp,nlon*nlat)
         used = send_data (id_storage_c(i_species), tmp(isc:iec,jsc:jec), River%Time) !, &
                        !  mask=River%gmask(isc:iec,jsc:jec) )
        endif
      if (id_infloc_c(i_species) > 0) then
         tmp = River%infloc_c(:,:,i_species)
 !        call mpp_sum(tmp,nlon*nlat)
         used = send_data (id_infloc_c(i_species), tmp(isc:iec,jsc:jec), River%Time) !, &
                       !   mask=River%gmask(isc:iec,jsc:jec) )
        endif
      if (id_removal_c(i_species) > 0) then
         tmp = River%removal_c(:,:,i_species)
         call mpp_sum(tmp,nlon*nlat)
         used = send_data (id_removal_c(i_species), tmp(isc:iec,jsc:jec), River%Time) !, &
                       !   mask=River%gmask(isc:iec,jsc:jec) )
        endif
      if (id_disc2l(i_species) > 0) then
         tmp = River%disc2l(:,:,i_species)
 !        call mpp_sum(tmp,nlon*nlat)
         used = send_data (id_disc2l(i_species), tmp(isc:iec,jsc:jec), River%Time)
        endif
      if (id_disc2o(i_species) > 0) then
         tmp = River%disc2o(:,:,i_species)
 !        call mpp_sum(tmp,nlon*nlat)
         used = send_data (id_disc2o(i_species), tmp(isc:iec,jsc:jec), River%Time)
        endif

      enddo

!    if (id_width > 0) then
!       tmp = River%width
!       call mpp_sum(tmp,nlon*nlat)
!       used = send_data (id_width, tmp(isc:iec,jsc:jec), River%Time, &
!                          mask=River%gmask(isc:iec,jsc:jec) )
!    endif
!
    if (id_depth > 0) then
       tmp = River%depth
       call mpp_sum(tmp,nlon*nlat)
       used = send_data (id_depth, tmp(isc:iec,jsc:jec), River%Time) !, &
                       !   mask=River%gmask(isc:iec,jsc:jec) )
    endif

!    if (id_vel > 0) then
!       tmp = River%vel
!       call mpp_sum(tmp,nlon*nlat)
!       used = send_data (id_vel, tmp(isc:iec,jsc:jec), River%Time, &
!                          mask=River%gmask(isc:iec,jsc:jec) )
!    endif
!

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

end module river_mod


#ifdef test_river_solo

program river_solo
  use mpp_mod,                  only : mpp_error, mpp_pe, mpp_root_pe, mpp_npes, FATAL
  use mpp_mod,                  only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_domains_mod,          only : mpp_define_layout, mpp_define_domains
  use mpp_domains_mod,          only : mpp_get_compute_domain, domain2d
  use mpp_io_mod,               only : mpp_open, MPP_RDONLY, MPP_NETCDF, MPP_SINGLE
  use mpp_io_mod,               only : MPP_ASCII, MPP_OVERWR, mpp_close
  use mpp_io_mod,               only : mpp_get_atts, mpp_get_axis_data, mpp_get_info
  use mpp_io_mod,               only : mpp_get_axes,  axistype
  use fms_mod,                  only : fms_init, fms_end, stdlog, open_namelist_file
  use fms_mod,                  only : check_nml_error, close_file, file_exist, stdout
  use fms_mod,                  only : field_size, read_data
  use fms_io_mod,               only : fms_io_exit
  use time_manager_mod,         only : time_type, increment_time, set_date, increment_date, set_time
  use time_manager_mod,         only : set_calendar_type, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only : operator(/), operator(-), operator( + ), month_name, get_date
  use diag_manager_mod,         only : diag_manager_init, diag_manager_end
  use river_mod,                only : river_init, river_end, update_river
  use constants_mod,            only : constants_init, PI, radius


  implicit none

  !--- namelist -----------------------------------------

  integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=16)     :: calendar = 'julian'
  integer               :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0
  integer               :: dt_fast     = 0


  namelist /river_solo_nml/ current_date, dt_fast, years, months, days, &
       hours, minutes, seconds, calendar
       
  character(len=128)    :: land_grid_file = 'INPUT/grid_spec.nc'

  !--------------------------------------------------------------------
  type(time_type)      :: Time, Time_start, Time_end, Run_len, Time_step_fast
  integer              :: nr, nf, num_fast_step, unit
  integer              :: yr,mon,day,hr,min,sec, calendar_type=-1
  real, allocatable    :: runoff(:,:), discharge(:,:)
  integer              :: initClock, mainClock, termClock, updateClock 


  call fms_init
  
  initClock = mpp_clock_id( 'Initialization' )
  mainClock = mpp_clock_id( 'Main loop' )
  termClock = mpp_clock_id( 'Termination' )
  updateClock = mpp_clock_id( 'update river')

  call mpp_clock_begin(initClock)
  call river_solo_init
  call mpp_clock_end (initClock) !end initialization

  call mpp_clock_begin(mainClock) !begin main loop
  do nf = 1, num_fast_step
    if(mpp_pe() == mpp_root_pe() ) write(stdout(),*)' at river fast time step ', nf
    call mpp_clock_begin(updateClock)
    call update_river ( runoff, discharge )
    call mpp_clock_end (updateClock)
    Time = Time + Time_step_fast
  enddo
  call mpp_clock_end(mainClock)

  call mpp_clock_begin(termClock)
  call river_end
  call diag_manager_end(Time)
  call get_date(Time,yr,mon,day,hr,min,sec)

  if (mpp_pe() == mpp_root_pe()) then
     call mpp_open(unit, 'RESTART/river_solo.res',form=MPP_ASCII,&
          action=MPP_OVERWR,threading=MPP_SINGLE,fileset=MPP_SINGLE,nohdrs=.true.)
     write(unit,*) yr, mon, day, hr, min, sec
     write(unit,*) calendar_type 
     call mpp_close(unit)
  endif

  call fms_io_exit
  call fms_end


contains

  !#######################################################################
  subroutine river_solo_init

    integer                     :: unit, ndim, nvar, natt, ntime, ierr, io, i
    integer                     :: len, siz(4), layout(2) = (/1,0/)
    integer                     :: ni, nj, npes, isc, iec, jsc, jec
    type(domain2d)              :: Domain
    real,   allocatable         :: lonb_lnd(:), latb_lnd(:)
    real,   allocatable         :: area_lnd(:,:), area_lnd_cell(:,:), gfrac(:,:)
    integer                     :: date(6)
    type(axistype), allocatable :: axes(:)
    character(len=16)           :: name
    character(len=9)            :: month

    call constants_init

    unit = open_namelist_file ()
    ierr=1
    do while (ierr /= 0)
       read  (unit, nml=river_solo_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'river_solo_nml')
    enddo
10  call close_file (unit)

    write(stdlog(), nml= river_solo_nml)

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
    if (file_exist('INPUT/river_solo.res')) then
       call mpp_open(unit,'INPUT/river_solo.res',form=MPP_ASCII,action=MPP_RDONLY)
       read(unit,*) date
       read(unit,*) calendar_type 
       call close_file(unit)
    endif

    call set_calendar_type (calendar_type)

    call diag_manager_init

    if (sum(current_date) <= 0) then
       call mpp_error(FATAL,'==>Error from river_solo_mod: no namelist value for current date')
    else
       Time_start  = set_date(current_date(1),current_date(2), current_date(3), &
            current_date(4),current_date(5),current_date(6))
    endif

    if (file_exist('INPUT/river_solo.res')) then
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

    call mpp_open (unit, 'time_stamp.out', form=MPP_ASCII, action=MPP_OVERWR,threading=MPP_SINGLE)

    month = month_name(current_date(2))
    if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

    call get_date (Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
    month = month_name(date(2))
    if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

    call close_file (unit)  

    !--- get the land grid and set up domain decomposition
    call field_size(land_grid_file, 'AREA_LND', siz)
    ni = siz(1); nj = siz(2)
    allocate( lonb_lnd(ni+1), latb_lnd(nj+1) )
    allocate( area_lnd(ni,nj), area_lnd_cell(ni,nj), gfrac(ni,nj) )
    call read_data(land_grid_file, 'xbl', lonb_lnd)
    call read_data(land_grid_file, 'ybl', latb_lnd)
    call read_data(land_grid_file, 'AREA_LND', area_lnd)
    call read_data(land_grid_file, 'AREA_LND_CELL', area_lnd_cell)
    lonb_lnd = lonb_lnd * PI/180.
    latb_lnd = latb_lnd * PI/180.
    gfrac = area_lnd/area_lnd_cell
    area_lnd = area_lnd*4*pi*radius**2

    npes = mpp_npes()

    !--- define domain ------------------------------------------------
    call mpp_define_layout((/1,ni,1,nj/),npes,layout)
    call mpp_define_domains((/1,ni,1,nj/),layout, Domain)
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)

    allocate(runoff(isc:iec,jsc:jec), discharge(isc:iec,jsc:jec) )

    call river_init( lonb_lnd, latb_lnd, Time_start, Time_step_fast, Domain, area_lnd(isc:iec,jsc:jec), gfrac  )

  end subroutine river_solo_init

  !#####################################################################

end program river_solo

#endif test_river_solo
