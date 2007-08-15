module land_data_mod

use mpp_domains_mod   , only : domain2d, mpp_get_compute_domain, &
     mpp_define_layout, mpp_define_domains, &
     CYCLIC_GLOBAL_DOMAIN
use fms_mod           , only : write_version_number, mpp_npes
use time_manager_mod  , only : time_type
use tracer_manager_mod, only : NO_TRACER

use land_tile_mod     , only : land_tile_type, land_tile_list_type, &
     land_tile_list_init, land_tile_list_end, nitems

implicit none
private

! ==== public interfaces =====================================================
public :: land_data_init
public :: land_data_end
public :: lnd            ! global data 

public :: atmos_land_boundary_type ! container for information passed from the 
                         ! atmosphere to land
public :: land_data_type ! container for information passed from land to 
                         ! the atmosphere
! both hold information on land grid (that is, after flux exchange translated 
! it from the atmosphere)

public :: dealloc_land2cplr ! deallocates a land_data_type structure
public :: realloc_land2cplr ! allocates a land_data_type members for current 
                            ! number of tiles
public :: dealloc_cplr2land ! deallocates an atmos_land_boundary_type structure
public :: realloc_cplr2land ! allocates an atmos_land_boundary_type members 
                            ! for current number of tiles
! NOTE: realloc_* procedures can be called regardless of the current state
! of the argument data structures, since they deallocate data first.

public :: land_state_type
! ==== end of public interfaces ==============================================

! ---- module constants ------------------------------------------------------
character(len=*),   parameter :: module_name = 'land_data_mod'
character(len=128), parameter :: version = '$Id: land_data.F90,v 15.0 2007/08/14 03:59:32 fms Exp $'
character(len=128), parameter :: tagname = '$Name: omsk $'

! ---- types -----------------------------------------------------------------
type :: atmos_land_boundary_type
   ! data passed from the coupler to the surface
   real, dimension(:,:,:), pointer :: & ! (lon, lat, tile)
        t_flux    => NULL(), &   ! sensible heat flux, W/m2
        lw_flux   => NULL(), &   ! net longwave radiadion flux, W/m2
        lwdn_flux => NULL(), &   ! downward longwave radiation flux, W/m2
        sw_flux   => NULL(), &   ! net shortwave radiation flux, W/m2
        swdn_flux => NULL(), &   ! downward shortwave radiation flux, W/m2
        lprec     => NULL(), &   ! liquid precip, kg/m2/s
        fprec     => NULL(), &   ! frozen precip, kg/m2/s
        tprec     => NULL(), &   !
   ! components of downward shortwave flux, W/m2  
        sw_flux_down_vis_dir   => NULL(), & 
        sw_flux_down_total_dir => NULL(), &
        sw_flux_down_vis_dif   => NULL(), &
        sw_flux_down_total_dif => NULL(), &  
   ! derivatives of the fluxes
        dhdt      => NULL(), &   ! sensible w.r.t. surf temperature
        dhdq      => NULL(), &   ! sensible w.r.t. surf humidity
        drdt      => NULL(), &   ! longwave w.r.t. surf radiative temperature 
   !
        cd_m      => NULL(), &   ! drag coefficient for momentum, dimensionless
        cd_t      => NULL(), &   ! drag coefficient for tracers, dimensionless
        ustar     => NULL(), &   ! turbulent wind scale, m/s
        bstar     => NULL(), &   ! turbulent bouyancy scale, m/s
        wind      => NULL(), &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot     => NULL(), &   ! height of the bottom atmos. layer above surface, m
        drag_q    => NULL(), &   ! product of cd_q by wind
        p_surf    => NULL()      ! surface pressure, Pa

#ifdef LAND_BND_TRACERS
   real, dimension(:,:,:,:), pointer :: & ! (lon, lat, tile,tracer)
        tr_flux => NULL(),   &   ! tracer flux, including water vapor flux
        dfdtr   => NULL()        ! tendency, including evap over surface specific humidity
#else
   real, dimension(:,:,:), pointer :: & ! (lon, lat, tile)
        q_flux => NULL(),    &   ! water vapor flux
        dedt => NULL(),      &   ! evap over surface temperature (assuming saturation @ surf.)
        dedq => NULL()           ! evap over surface specufuc humidity
#endif 

   integer :: xtype             !REGRID, REDIST or DIRECT
end type atmos_land_boundary_type


type :: land_data_type
   logical :: pe ! data presence indicator for stock calculations
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size      => NULL(),  & ! fractional coverage of cell by tile, dimensionless
        t_surf         => NULL(),  & ! ground surface temperature, degK
        t_ca           => NULL(),  & ! canopy air temperature, degK
        albedo         => NULL(),  & ! snow-adjusted land albedo
        albedo_vis_dir => NULL(),  & ! albedo for direct visible-band radiation
        albedo_nir_dir => NULL(),  & ! albedo for direct nir-band radiation 
        albedo_vis_dif => NULL(),  & ! albedo for diffuse visible-band radiation 
        albedo_nir_dif => NULL(),  & ! albedo for diffuse nir-band radiation
        rough_mom      => NULL(),  & ! momentum roughness length, m
        rough_heat     => NULL(),  & ! roughness length for tracers (heat and water), m
        rough_scale    => NULL()     ! topographic scaler for momentum drag, m

#ifdef LAND_BND_TRACERS
   real, pointer, dimension(:,:,:,:)   :: &  ! (lon, lat, tile, tracer)
        tr    => NULL()              ! tracers, including canopy air specific humidity, kg/kg
#else
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        q_ca => NULL()               ! canopy air specific humidity, kg/kg
#endif

   real, pointer, dimension(:,:) :: &  ! (lon, lat)
        discharge      => NULL(),  & ! flux from surface drainage network out of land model
        discharge_snow => NULL()     ! snow analogue of discharge

   logical, pointer, dimension(:,:,:):: &
        mask => NULL()               ! true if land

   integer :: axes(2)        ! IDs of diagnostic axes
   type(domain2d) :: domain  ! our computation domain
   logical, pointer :: maskmap(:,:) 

end type land_data_type


! land_state_type combines the general information about state of the land model:
! domain, coordinates, time steps, etc. There is only one variable of this type,
! and it is public in this module.
type :: land_state_type
   integer        :: is,ie,js,je ! compute domain boundaries
   integer        :: ntprog      ! number of prognostic tracers
   integer        :: isphum      ! index of specific humidity in tracer table
   integer        :: ico2        ! index of carbon dioxide in tracer table
   type(time_type):: dt_fast     ! fast (physical) time step
   type(time_type):: dt_slow     ! slow time step
   type(time_type):: time        ! current time

   real, pointer  :: lon(:), lat(:)     ! local grid coordinates, radian
   real, pointer  :: lonb(:), latb(:)   ! local grid boundaries, radian
   real, pointer  :: lon2(:,:), lat2(:,:) ! 2-d version, useful in some cases 
           ! (e.g diurnal_solar calculations)
   real, pointer  :: glon(:), glat(:)   ! global grid coordinates, radian
   real, pointer  :: glonb(:), glatb(:) ! global grid boundaries, radian
   real, pointer  :: garea(:,:)         ! land area per gridcell, m2

   ! map of tiles
   type(land_tile_list_type), pointer :: tile_map(:,:)
   
   type(domain2d) :: domain ! our domain -- be last since it simplifies
                            ! debugging in totalview
end type land_state_type

! ---- public module variables -----------------------------------------------
type(land_state_type),save :: lnd


! ---- private module variables ----------------------------------------------
logical :: module_is_initialized =.FALSE.


#define __DEALLOC__(x) if (associated(x)) deallocate(x)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



! ============================================================================
subroutine land_data_init(glon, glat, glonb, glatb, garea, layout)
  real   , intent(in) :: glon(:), glat(:), glonb(:), glatb(:)
  real   , intent(in) :: garea(:,:)   ! area of land in grid cells
  integer, intent(inout) :: layout(2) ! layout of our domains

  ! ---- local vars
  integer :: nlon, nlat ! size of global grid in lon and lat directions
  integer :: i,j

  ! write the version and tagname to the logfile
  call write_version_number(version, tagname)

  ! calculate size of the global grid
  nlon = size(glonb) - 1
  nlat = size(glatb) - 1

  ! get the processor layout information, either from upper-layer module
  ! decomposition, or define according to the grid size
  if( layout(1)==0 .AND. layout(2)==0 ) &
       call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes(), layout )
  if( layout(1)/=0 .AND. layout(2)==0 )layout(2) = mpp_npes()/layout(1)
  if( layout(1)==0 .AND. layout(2)/=0 )layout(1) = mpp_npes()/layout(2)

  call mpp_define_domains ( &
       (/1,nlon, 1, nlat/),           &  ! global grid size
       layout,                        &  ! layout for our domains
       lnd%domain,                    &  ! domain to define
       xflags = CYCLIC_GLOBAL_DOMAIN, &
       name = 'LAND MODEL')

  call mpp_get_compute_domain(lnd%domain, lnd%is,lnd%ie,lnd%js,lnd%je)

  allocate(lnd%tile_map (lnd%is:lnd%ie, lnd%js:lnd%je))
  allocate(lnd%glonb    (size(glonb)))
  allocate(lnd%glatb    (size(glatb)))
  allocate(lnd%glon     (size(glon)))
  allocate(lnd%glat     (size(glat)))
  allocate(lnd%garea    (size(glon),size(glat)))
  allocate(lnd%lonb     (lnd%is:lnd%ie+1))
  allocate(lnd%latb     (lnd%js:lnd%je+1))
  allocate(lnd%lon      (lnd%is:lnd%ie))
  allocate(lnd%lat      (lnd%js:lnd%je))
  allocate(lnd%lon2     (lnd%is:lnd%ie, lnd%js:lnd%je))
  allocate(lnd%lat2     (lnd%is:lnd%ie, lnd%js:lnd%je))

  ! initialize coordinates
  lnd%glonb = glonb
  lnd%glatb = glatb
  lnd%lonb  = glonb(lnd%is:lnd%ie+1)
  lnd%latb  = glatb(lnd%js:lnd%je+1)

  lnd%glon  = glon
  lnd%glat  = glat
  lnd%lon   = glon(lnd%is:lnd%ie)
  lnd%lat   = glat(lnd%js:lnd%je)
  do i = lnd%js,lnd%je
     lnd%lon2(:,i) = lnd%lon
  enddo
  do i = lnd%is,lnd%ie
     lnd%lat2(i,:) = lnd%lat
  enddo
  
  ! initialize land tile map
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     call land_tile_list_init(lnd%tile_map(i,j))
  enddo
  enddo

  ! initialize land_area
  lnd%garea = garea

end subroutine land_data_init

! ============================================================================
subroutine land_data_end()

  integer :: i,j
  
  module_is_initialized = .FALSE.

  ! dealocate land tile map here. 
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     call land_tile_list_end(lnd%tile_map(i,j))
  enddo
  enddo

  ! deallocate grid data
  deallocate(&
       lnd%glonb, lnd%glatb, lnd%glon, lnd%glat, &
       lnd%lonb,  lnd%latb,  lnd%lon,  lnd%lat,  &
       lnd%lon2,  lnd%lat2,  lnd%garea           )

end subroutine land_data_end


! ============================================================================
! allocates boundary data for land domain and current number of tiles
subroutine realloc_land2cplr ( bnd )
  type(land_data_type), intent(inout) :: bnd     ! data to allocate

  ! ---- local vars
  integer :: n_tiles

  call dealloc_land2cplr(bnd)

  bnd%domain = lnd%domain
  n_tiles = max_n_tiles()


  ! allocate data according to the domain boundaries
  allocate( bnd%mask(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )

  allocate( bnd%tile_size(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%t_surf(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%t_ca(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
#ifdef LAND_BND_TRACERS
  allocate( bnd%tr(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles,lnd%ntprog) )
#else
  allocate( bnd%q_ca(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
#endif
  allocate( bnd%albedo(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_vis_dir(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_nir_dir(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_vis_dif(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_nir_dif(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%rough_mom(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%rough_heat(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%rough_scale(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )

  allocate( bnd%discharge(lnd%is:lnd%ie,lnd%js:lnd%je) )
  allocate( bnd%discharge_snow(lnd%is:lnd%ie,lnd%js:lnd%je) )

end subroutine realloc_land2cplr


! ============================================================================
! deallocates boundary data memory
subroutine dealloc_land2cplr ( bnd )
  type(land_data_type), intent(inout) :: bnd  ! data to de-allocate

  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%t_surf )
  __DEALLOC__( bnd%t_ca )
#ifdef LAND_BND_TRACERS
  __DEALLOC__( bnd%tr )
#else
  __DEALLOC__( bnd%q_ca )
#endif
  __DEALLOC__( bnd%albedo )
  __DEALLOC__( bnd%albedo_vis_dir )
  __DEALLOC__( bnd%albedo_nir_dir )
  __DEALLOC__( bnd%albedo_vis_dif )
  __DEALLOC__( bnd%albedo_nir_dif )
  __DEALLOC__( bnd%rough_mom )
  __DEALLOC__( bnd%rough_heat )
  __DEALLOC__( bnd%rough_scale )
  __DEALLOC__( bnd%discharge )
  __DEALLOC__( bnd%discharge_snow )
  __DEALLOC__( bnd%mask )

end subroutine dealloc_land2cplr


! ============================================================================
! allocates boundary data for land domain and current number of tiles;
! initializes data for datat override.
! NOTE: previously the body of the procedure was in the flux_exchange_init,
! currently it is called from land_model_init
subroutine realloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  ! ---- local vars
  integer :: kd

  call dealloc_cplr2land(bnd)

  ! allocate data according to the domain boundaries
  kd = max_n_tiles()

  allocate( bnd%t_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%lw_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%lprec(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%fprec(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%tprec(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%dhdt(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%dhdq(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%drdt(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%p_surf(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
#ifdef LAND_BND_TRACERS
  allocate( bnd%tr_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd,lnd%ntprog) )
  allocate( bnd%dfdtr(lnd%is:lnd%ie,lnd%js:lnd%je,kd,lnd%ntprog) )
#else
  allocate( bnd%q_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%dedt(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%dedq(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
#endif
!slm[unwated init?]  bnd%q_flux=0.0
!slm[unwated init?]  bnd%dedt=0.0
!slm[unwated init?]  bnd%dedq=0.0

! initialize boundary values for override experiments (mjh)
!slm[unwated init?]  bnd%t_flux=0.0
!slm[unwated init?]  bnd%lw_flux=0.0
!slm[unwated init?]  bnd%sw_flux=0.0
!slm[unwated init?]  bnd%lprec=0.0
!slm[unwated init?]  bnd%fprec=0.0
!slm[unwated init?]  bnd%dhdt=0.0
!slm[unwated init?]  bnd%drdt=0.0
!slm[unwated init?]  bnd%p_surf=0.0

! + slm May 13 2003 -- fields for lm3 
  allocate( bnd%lwdn_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%swdn_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_vis_dir(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_total_dir(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_vis_dif(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_total_dif(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%cd_t(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%cd_m(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%bstar(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%ustar(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%wind(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%z_bot(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )

!slm[unwated init?]  bnd%lwdn_flux=0.0
!slm[unwated init?]  bnd%swdn_flux=0.0
!slm[unwated init?]  bnd%sw_flux_down_vis_dir=0.0
!slm[unwated init?]  bnd%sw_flux_down_total_dir=0.0
!slm[unwated init?]  bnd%sw_flux_down_vis_dif=0.0
!slm[unwated init?]  bnd%sw_flux_down_total_dif=0.0
!slm[unwated init?]  bnd%cd_t = 1e-3
!slm[unwated init?]  bnd%cd_m = 1e-3
!slm[unwated init?]  bnd%bstar = 0.0
!slm[unwated init?]  bnd%ustar = 0.0
!slm[unwated init?]  bnd%wind = 0.0
!slm[unwated init?]  bnd%z_bot = 30.0
! - slm May 13 2003

  allocate( bnd%drag_q(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )


end subroutine realloc_cplr2land


! ============================================================================
subroutine dealloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  __DEALLOC__( bnd%t_flux )
  __DEALLOC__( bnd%lw_flux )
  __DEALLOC__( bnd%sw_flux )
  __DEALLOC__( bnd%lprec )
  __DEALLOC__( bnd%fprec )
  __DEALLOC__( bnd%dhdt )
  __DEALLOC__( bnd%dhdq )
  __DEALLOC__( bnd%drdt )
  __DEALLOC__( bnd%p_surf )
  __DEALLOC__( bnd%lwdn_flux )
  __DEALLOC__( bnd%swdn_flux )
  __DEALLOC__( bnd%sw_flux_down_vis_dir )
  __DEALLOC__( bnd%sw_flux_down_total_dir )
  __DEALLOC__( bnd%sw_flux_down_vis_dif )
  __DEALLOC__( bnd%sw_flux_down_total_dif )
  __DEALLOC__( bnd%cd_t )
  __DEALLOC__( bnd%cd_m )
  __DEALLOC__( bnd%bstar )
  __DEALLOC__( bnd%ustar )
  __DEALLOC__( bnd%wind )
  __DEALLOC__( bnd%z_bot )
#ifdef LAND_BND_TRACERS
  __DEALLOC__( bnd%tr_flux )
  __DEALLOC__( bnd%dfdtr )
#else
  __DEALLOC__( bnd%q_flux )
  __DEALLOC__( bnd%dedt )
  __DEALLOC__( bnd%dedq )
#endif

end subroutine dealloc_cplr2land

! ============================================================================
! get max number of tiles in the domain
function max_n_tiles() result(n)
  integer :: n
  integer :: i,j

  n=1
  do j=lnd%js,lnd%je
  do i=lnd%is,lnd%ie
     n=max(n, nitems(lnd%tile_map(i,j)))
  enddo
  enddo

end function 

end module land_data_mod
