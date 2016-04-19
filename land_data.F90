module land_data_mod

use mpp_mod           , only : mpp_get_current_pelist, mpp_pe
use constants_mod     , only : PI
use mpp_domains_mod   , only : domain2d, mpp_get_compute_domain, &
     mpp_define_layout, mpp_define_domains, mpp_define_io_domain, &
     mpp_get_current_ntile, mpp_get_tile_id, CYCLIC_GLOBAL_DOMAIN, &
     mpp_get_io_domain, mpp_get_pelist, mpp_get_layout, mpp_get_domain_npes
use fms_mod           , only : write_version_number, mpp_npes, stdout, &
     file_exist, error_mesg, FATAL
use fms_io_mod        , only : parse_mask_table
use time_manager_mod  , only : time_type
use grid_mod          , only : get_grid_ntiles, get_grid_size, get_grid_cell_vertices, &
     get_grid_cell_centers, get_grid_cell_area, get_grid_comp_area, &
     define_cube_mosaic

implicit none
private

! ==== public interfaces =====================================================
public :: land_data_init
public :: land_data_end
public :: lnd            ! global data 
public :: log_version    ! prints version number

public :: atmos_land_boundary_type ! container for information passed from the 
                         ! atmosphere to land
public :: land_data_type ! container for information passed from land to 
                         ! the atmosphere
public :: land_state_type
! ==== end of public interfaces ==============================================

! ---- module constants ------------------------------------------------------
character(len=*), parameter :: module_name = 'land_data_mod'
#include "shared/version_variable.inc"

! ---- types -----------------------------------------------------------------
type :: atmos_land_boundary_type
   ! data passed from the coupler to the surface
   real, dimension(:,:,:), pointer :: & ! (lon, lat, tile)
        t_flux    => NULL(), &   ! sensible heat flux, W/m2
        lw_flux   => NULL(), &   ! net longwave radiation flux, W/m2
        lwdn_flux => NULL(), &   ! downward longwave radiation flux, W/m2
        sw_flux   => NULL(), &   ! net shortwave radiation flux, W/m2
        swdn_flux => NULL(), &   ! downward shortwave radiation flux, W/m2
        lprec     => NULL(), &   ! liquid precipitation rate, kg/(m2 s)
        fprec     => NULL(), &   ! frozen precipitation rate, kg/(m2 s)
        tprec     => NULL(), &   ! temperature of precipitation, degK
   ! components of downward shortwave flux, W/m2  
        sw_flux_down_vis_dir   => NULL(), & ! visible direct 
        sw_flux_down_total_dir => NULL(), & ! total direct
        sw_flux_down_vis_dif   => NULL(), & ! visible diffuse
        sw_flux_down_total_dif => NULL(), & ! total diffuse
   ! derivatives of the fluxes
        dhdt      => NULL(), &   ! sensible w.r.t. surface temperature
        dhdq      => NULL(), &   ! sensible w.r.t. surface humidity
        drdt      => NULL(), &   ! longwave w.r.t. surface radiative temperature 
   !
        cd_m      => NULL(), &   ! drag coefficient for momentum, dimensionless
        cd_t      => NULL(), &   ! drag coefficient for tracers, dimensionless
        ustar     => NULL(), &   ! turbulent wind scale, m/s
        bstar     => NULL(), &   ! turbulent buoyancy scale, m/s
        wind      => NULL(), &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot     => NULL(), &   ! height of the bottom atmospheric layer above the surface, m
        drag_q    => NULL(), &   ! product of cd_q by wind
        p_surf    => NULL()      ! surface pressure, Pa

   real, dimension(:,:,:,:), pointer :: & ! (lon, lat, tile, tracer)
        tr_flux => NULL(),   &   ! tracer flux, including water vapor flux
        dfdtr   => NULL()        ! derivative of the flux w.r.t. tracer surface value, 
                                 ! including evap over surface specific humidity

   integer :: xtype             !REGRID, REDIST or DIRECT
end type atmos_land_boundary_type


type :: land_data_type
   ! data passed from the surface to the coupler
   logical :: pe ! data presence indicator for stock calculations
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size      => NULL(),  & ! fractional coverage of cell by tile, dimensionless
        t_surf         => NULL(),  & ! ground surface temperature, degK
        t_ca           => NULL(),  & ! canopy air temperature, degK
        albedo         => NULL(),  & ! broadband land albedo [unused?]
        albedo_vis_dir => NULL(),  & ! albedo for direct visible radiation
        albedo_nir_dir => NULL(),  & ! albedo for direct NIR radiation 
        albedo_vis_dif => NULL(),  & ! albedo for diffuse visible radiation 
        albedo_nir_dif => NULL(),  & ! albedo for diffuse NIR radiation
        rough_mom      => NULL(),  & ! surface roughness length for momentum, m
        rough_heat     => NULL(),  & ! roughness length for tracers and heat, m
        rough_scale    => NULL()     ! topographic scaler for momentum drag, m

   real, pointer, dimension(:,:,:,:)   :: &  ! (lon, lat, tile, tracer)
        tr    => NULL()              ! tracers, including canopy air specific humidity

   ! NOTE that in contrast to most of the other fields in this structure, the discharges
   ! hold data per-gridcell, rather than per-tile basis. This, and the order of updates,
   ! have implications for the data reallocation procedure.
   real, pointer, dimension(:,:) :: &  ! (lon, lat)
     discharge           => NULL(),  & ! liquid water flux from land to ocean
     discharge_heat      => NULL(),  & ! sensible heat of discharge (0 C datum)
     discharge_snow      => NULL(),  & ! solid water flux from land to ocean
     discharge_snow_heat => NULL()     ! sensible heat of discharge_snow (0 C datum)

   logical, pointer, dimension(:,:,:):: &
        mask => NULL()               ! true if land

   integer :: axes(2)        ! IDs of diagnostic axes
   type(domain2d) :: domain  ! our computation domain
   integer, pointer, dimension(:) :: pelist
end type land_data_type


! land_state_type combines the general information about state of the land model:
! domain, coordinates, time steps, etc. There is only one variable of this type,
! and it is public in this module.
type :: land_state_type
   integer        :: is,ie,js,je ! compute domain boundaries
   integer        :: nlon,nlat   ! size of global grid
   type(time_type):: dt_fast     ! fast (physical) time step
   type(time_type):: dt_slow     ! slow time step
   
   type(time_type):: time        ! current land model time

   real, allocatable  :: lon (:,:), lat (:,:) ! domain grid center coordinates, radian
   real, allocatable  :: lonb(:,:), latb(:,:) ! domain grid vertices, radian
   real, allocatable  :: area(:,:)  ! land area per grid cell, m2
   real, allocatable  :: cellarea(:,:)  ! grid cell area, m2
   real, allocatable  :: landfrac(:,:)  ! fraction of land in the grid cell
   real, allocatable  :: coord_glon(:), coord_glonb(:) ! longitudes for use in diag axis and such, degrees East
   real, allocatable  :: coord_glat(:), coord_glatb(:) ! latitudes for use in diag axis and such, degrees North
   
   integer :: nfaces ! number of mosaic faces
   integer :: face  ! the current mosaic face
   integer, allocatable :: pelist(:) ! list of processors that run land model
   integer, allocatable :: io_pelist(:) ! list of processors in our io_domain
   ! if io_domain was not defined, then there is just one element in this
   ! array, and it's equal to current PE
   integer :: io_id     ! suffix in the distributed files.
   logical :: append_io_id ! if FALSE, io_id is not appended to the file names
                           ! (for the case io_layout = 1,1)

   type(domain2d) :: domain ! our domain -- should be the last since it simplifies
                            ! debugging in totalview
end type land_state_type

! ---- public module variables -----------------------------------------------
type(land_state_type), save :: lnd

! ---- private module variables ----------------------------------------------
logical :: module_is_initialized = .FALSE.


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
subroutine log_version(version, modname, filename, tag, unit)
  character(len=*), intent(in) :: version
  character(len=*), intent(in), optional :: &
          modname, filename, tag
  integer, intent(in), optional :: unit

  character(512) :: message
  integer :: i
  message=''

  if (present(filename)) then
     ! remove the directory part of the name
     i = scan(filename,'/',back=.true.)

     message = trim(filename(i+1:))
  endif
  if (present(modname)) then
     message = trim(modname)//' ('//trim(message)//')'
  endif
  call write_version_number (trim(message)//': '//trim(version),tag,unit)
end subroutine log_version

! ============================================================================
subroutine land_data_init(layout, io_layout, time, dt_fast, dt_slow, mask_table)
  integer, intent(inout) :: layout(2) ! layout of our domains
  integer, intent(inout) :: io_layout(2) ! layout for land model io
  type(time_type), intent(in) :: &
       time,    & ! current model time
       dt_fast, & ! fast (physical) time step
       dt_slow    ! slow time step
  character(len=*), intent(in) :: mask_table

  ! ---- local vars
  integer :: nlon, nlat ! size of global grid in lon and lat directions
  integer :: ntiles     ! number of tiles in the mosaic grid 
  integer :: ntracers, ndiag ! non-optional output from register_tracers
  integer, allocatable :: tile_ids(:) ! mosaic tile IDs for the current PE
  type(domain2d), pointer :: io_domain ! our io_domain
  integer :: n_io_pes ! number of PEs in our io_domain
  integer :: io_id(1)
  integer :: outunit
  logical :: mask_table_exist
  logical, allocatable :: maskmap(:,:,:)  ! A pointer to an array indicating which
                                          ! logical processors are actually used for
                                          ! the land code. The other logical
                                          ! processors would be all ocean points and
                                          ! are not assigned to actual processors.
                                          ! This need not be assigned if all logical
                                          ! processors are used. 

  ! write the version and tag name to the logfile
  call log_version(version, module_name, &
  __FILE__)

  ! define the processor layout information according to the global grid size 
  call get_grid_ntiles('LND',ntiles)
  lnd%nfaces = ntiles
  call get_grid_size('LND',1,nlon,nlat)
  ! set the size of global grid
  lnd%nlon = nlon; lnd%nlat = nlat

  mask_table_exist = .false.
  outunit = stdout()
  if(file_exist(mask_table)) then
     mask_table_exist = .true.
     write(outunit, *) '==> NOTE from land_data_init:  reading maskmap information from '//trim(mask_table)
     if(layout(1) == 0 .OR. layout(2) == 0 ) call error_mesg('land_model_init', &
        'land_model_nml layout should be set when file '//trim(mask_table)//' exists', FATAL)

     allocate(maskmap(layout(1), layout(2), ntiles))
     call parse_mask_table(mask_table, maskmap, "Land model")
  endif

  if( layout(1)==0 .AND. layout(2)==0 ) &
       call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes()/ntiles, layout )
  if( layout(1)/=0 .AND. layout(2)==0 )layout(2) = mpp_npes()/(layout(1)*ntiles)
  if( layout(1)==0 .AND. layout(2)/=0 )layout(1) = mpp_npes()/(layout(2)*ntiles)

  if( io_layout(1) == 0 .AND. io_layout(2) == 0 ) io_layout = layout

  ! define land model domain
  if (ntiles==1) then
     if( mask_table_exist ) then
        call mpp_define_domains ((/1,nlon, 1, nlat/), layout, lnd%domain, xhalo=1, yhalo=1,&
             xflags = CYCLIC_GLOBAL_DOMAIN, name = 'LAND MODEL', maskmap=maskmap(:,:,1) )
     else
        call mpp_define_domains ((/1,nlon, 1, nlat/), layout, lnd%domain, xhalo=1, yhalo=1,&
             xflags = CYCLIC_GLOBAL_DOMAIN, name = 'LAND MODEL')
     endif
  else
     if( mask_table_exist ) then
        call define_cube_mosaic ('LND', lnd%domain, layout, halo=1, maskmap=maskmap)
     else
        call define_cube_mosaic ('LND', lnd%domain, layout, halo=1)
     endif
  endif
  if(mask_table_exist) deallocate(maskmap)

  ! define io domain
  call mpp_define_io_domain(lnd%domain, io_layout)

  ! set up list of processors for collective io: only the first processor in this
  ! list actually writes data, the rest just send the data to it.
  io_domain=>mpp_get_io_domain(lnd%domain)
  if (associated(io_domain)) then
     n_io_pes = mpp_get_domain_npes(io_domain)
     allocate(lnd%io_pelist(n_io_pes))
     call mpp_get_pelist(io_domain,lnd%io_pelist)
     io_id = mpp_get_tile_id(io_domain)
     lnd%io_id = io_id(1)
  else
     call error_mesg('land_data_init','io_domain is undefined, contact developer', FATAL)
  endif
  lnd%append_io_id = (io_layout(1)/=1.or.io_layout(2)/=1)

  ! get the domain information
  call mpp_get_compute_domain(lnd%domain, lnd%is,lnd%ie,lnd%js,lnd%je)

  ! get the mosaic tile number for this processor: this assumes that there is only one
  ! mosaic tile per PE.
  allocate(tile_ids(mpp_get_current_ntile(lnd%domain)))
  tile_ids = mpp_get_tile_id(lnd%domain)
  lnd%face = tile_ids(1)
  deallocate(tile_ids)

  allocate(lnd%lonb    (lnd%is:lnd%ie+1, lnd%js:lnd%je+1))
  allocate(lnd%latb    (lnd%is:lnd%ie+1, lnd%js:lnd%je+1))
  allocate(lnd%lon     (lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%lat     (lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%area    (lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%cellarea(lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%coord_glon(nlon), lnd%coord_glonb(nlon+1))
  allocate(lnd%coord_glat(nlat), lnd%coord_glatb(nlat+1))

  ! initialize coordinates
  call get_grid_cell_vertices('LND',lnd%face,lnd%coord_glonb,lnd%coord_glatb)
  call get_grid_cell_centers ('LND',lnd%face,lnd%coord_glon, lnd%coord_glat)
  call get_grid_cell_area    ('LND',lnd%face,lnd%cellarea, domain=lnd%domain)
  call get_grid_comp_area    ('LND',lnd%face,lnd%area,     domain=lnd%domain)
  
  ! set local coordinates arrays -- temporary, till such time as the global arrays
  ! are not necessary
  call get_grid_cell_vertices('LND',lnd%face,lnd%lonb,lnd%latb, domain=lnd%domain)
  call get_grid_cell_centers ('LND',lnd%face,lnd%lon, lnd%lat, domain=lnd%domain)
  ! convert coordinates to radian; note that 1D versions stay in degrees
  lnd%lonb = lnd%lonb*pi/180.0 ; lnd%lon = lnd%lon*pi/180.0
  lnd%latb = lnd%latb*pi/180.0 ; lnd%lat = lnd%lat*pi/180.0
  
  ! initialize model's time-related parameters
  lnd%time    = time
  lnd%dt_fast = dt_fast
  lnd%dt_slow = dt_slow

  ! initialize the land model processor list
  allocate(lnd%pelist(0:mpp_npes()-1))
  call mpp_get_current_pelist(lnd%pelist)
end subroutine land_data_init

! ============================================================================
subroutine land_data_end()
  module_is_initialized = .FALSE.
end subroutine land_data_end

end module land_data_mod
