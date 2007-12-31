! ============================================================================
! top-level core of the Land Dynamics (LaD) model code
! ============================================================================
#define __DEBUG1__(x) write(*,'(a12,99g)')#x,x
#define __DEBUG2__(x1,x2) write(*,'(99(a12,g))')#x1,x1,#x2,x2 
#define __DEBUG3__(x1,x2,x3) write(*,'(99(a12,g))')#x1,x1,#x2,x2,#x3,x3 
#define __DEBUG4__(x1,x2,x3,x4) write(*,'(99(a12,g))')#x1,x1,#x2,x2,#x3,x3,#x4,x4 
#define __DEBUG5__(x1,x2,x3,x4,x5) write(*,'(99(a12,g))')#x1,x1,#x2,x2,#x3,x3,#x4,x4,#x5,x5 

module land_model_mod

use time_manager_mod, only : time_type, get_time, increment_time, time_type_to_real, &
     operator(+)
use mpp_domains_mod, only : domain2d, mpp_get_ntile_count
use mpp_mod, only : mpp_max
use fms_mod, only : write_version_number, error_mesg, FATAL, WARNING, NOTE, mpp_pe, &
     mpp_root_pe, file_exist, open_namelist_file, check_nml_error, close_file, &
     stdlog, get_mosaic_tile_file, mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
     CLOCK_FLAG_DEFAULT, CLOCK_COMPONENT, CLOCK_ROUTINE
use diag_manager_mod, only : diag_axis_init, register_static_field, &
     register_diag_field, send_data
use constants_mod, only : radius, hlf, hlv, hls, tfreeze, pi, rdgas, rvgas, cp_air, &
     stefan
use astronomy_mod, only : diurnal_solar
use sphum_mod, only : qscomp

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR
use glacier_mod, only : read_glac_namelist, glac_init, glac_end, glac_get_sfc_temp, &
     glac_radiation, glac_diffusion, glac_step_1, glac_step_2
use lake_mod, only : read_lake_namelist, lake_init, lake_end, lake_get_sfc_temp, &
     lake_radiation, lake_diffusion, lake_step_1, lake_step_2
use soil_mod, only : read_soil_namelist, soil_init, soil_end, soil_get_sfc_temp, &
     soil_radiation, soil_diffusion, soil_step_1, soil_step_2
use snow_mod, only : read_snow_namelist, snow_init, snow_end, snow_get_sfc_temp, &
     snow_radiation, snow_diffusion, snow_get_depth_area, snow_step_1, snow_step_2
use vegetation_mod, only : read_vegn_namelist, vegn_init, vegn_end, vegn_get_cover, &
     vegn_radiation, vegn_diffusion, vegn_step_1, vegn_step_2, vegn_step_3, &
     update_vegn_slow
use canopy_air_mod, only : read_cana_namelist, cana_init, cana_end, cana_state,&
     cana_step_1, cana_step_2, cana_radiation, cana_roughness
use river_mod, only : river_init, river_end, update_river
use topo_rough_mod, only : topo_rough_init, topo_rough_end, update_topo_rough
use soil_tile_mod, only : soil_cover_cold_start
use vegn_tile_mod, only : vegn_cover_cold_start, vegn_data_rs_min, update_derived_vegn_data
use lake_tile_mod, only : lake_cover_cold_start
use glac_tile_mod, only : glac_pars_type, glac_cover_cold_start
use numerics_mod, only : ludcmp, lubksb
use land_tile_mod, only : land_tile_type, land_tile_list_type, &
     land_tile_enum_type, new_land_tile, insert, nitems, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=), &
     get_elmt_indices, get_tile_tags
use land_data_mod, only : land_data_type, atmos_land_boundary_type, &
     land_state_type, land_data_init, land_data_end, lnd, &
     dealloc_land2cplr, realloc_land2cplr, &
     dealloc_cplr2land, realloc_cplr2land
use nf_utils_mod,  only : nfu_inq_var, nfu_inq_dim, nfu_get_var_int, &
     nfu_get_var_double
use land_io_mod,   only : nearest
use land_tile_io_mod, only : print_netcdf_error, create_tile_out_file, &
    write_tile_data, read_tile_data_r0d_fptr, write_tile_data_r0d_fptr
use land_tile_diag_mod, only : tile_diag_init, tile_diag_end, &
    register_tiled_diag_field, send_tile_data, dump_tile_diag_fields, &
    OP_AVERAGE, OP_SUM
use land_debug_mod, only : land_debug_init, land_debug_end, set_current_point, &
     is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations
public update_land_model_fast   ! time-step integration
public update_land_model_slow   ! time-step integration
public atmos_land_boundary_type ! data from coupler to land
public land_data_type           ! data from land to coupler

public :: Lnd_stock_pe          ! return stocks of conservative quantities
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'land', &
     version     = '$Id: land_model.F90,v 15.1.2.6.2.1 2007/11/20 17:42:17 slm Exp $', &
     tagname     = '$Name: omsk_2007_12 $'

! ==== module variables ======================================================

! ---- namelist --------------------------------------------------------------
logical :: lm2            = .false.
logical :: use_bone_dry   = .false.
real    :: con_fac_large = 1.e6
real    :: con_fac_small = 1.e-6
real    :: cpw = 1952.  ! specific heat of water vapor at constant pressure
real    :: clw = 4218.  ! specific heat of water (liquid)
real    :: csw = 2106.  ! specific heat of water (ice)
real    :: canopy_air_mass = 10.0 ! kg/m2
integer :: num_c = 0
logical :: do_age = .false.
integer :: layout(2) = (/0,0/)
namelist /land_model_nml/ lm2, use_bone_dry, canopy_air_mass, &
                          num_c, do_age, &
                          con_fac_large, con_fac_small, cpw, clw, csw, layout
! ---- end of namelist -------------------------------------------------------

logical  :: module_is_initialized = .FALSE.
logical  :: use_E_max
character(len=256) :: grid_spec_file="INPUT/grid_spec.nc" 
real     :: delta_time ! duration of main land time step (s)
real,     pointer, dimension(:,:,:) ::  zr_lm2
integer,  pointer, dimension(:,:,:) ::  nz_lm2
integer  :: num_species
integer  :: num_phys = 2
real, allocatable :: gfrac    (:,:)    ! fraction of land in cells

! ---- diag field IDs
integer :: id_ntiles, id_frac, id_landarea, id_landfrac, &
  id_lprecv, id_lprecs, id_lprecg, &
  id_fprecv, id_fprecs, &
  id_levapv, id_levaps, id_levapg, &
  id_fevapv, id_fevaps, id_fevapg, &
  id_lrunfs, id_lrunfg, id_frunfs, &
  id_transp, id_meltv, id_melts, &
  id_meltg, id_fswv, id_fsws, id_fswg, id_flwv, id_flws, id_flwg, &
  id_sensv, id_senss, id_sensg, id_gflux, id_hadvecv, id_hadvecs, &
  id_hadvecg, id_LWSv, id_LWSs, id_LWSg, id_FWSv, id_FWSs, id_FWSg, &
  id_HSv, id_HSs, id_HSg, id_Trad, id_Tca, id_qca, &
  id_z0m, id_z0s, id_area, &
  id_lprec_l, id_fprec_l, id_levap, id_fevap, id_lrunf, id_frunf, &
  id_melt, id_fsw, id_flw, id_sens, id_hadvec, id_LWS, id_FWS, id_HS, &
  id_precip, id_wroff, id_sroff, id_water, id_snow, id_groundwater, &
  id_evap, &
  id_swdn_dir, id_swdn_dif, id_swup_dir, id_swup_dif, &
  id_lwdn, &
  id_albedo_dir, id_albedo_dif, &
  id_subs_refl_dif, id_subs_refl_dir, & 
  id_vegn_refl_dif, id_vegn_refl_dir, id_vegn_tran_dif, id_vegn_tran_dir, &
  id_vegn_sctr_dir, id_cosz, id_vegn_tran_lw, id_vegn_refl_lw, &
  id_vegn_cover, id_con_g_h, id_grnd_T, &
  id_e_res_1, id_e_res_2, &
  id_dis_l2o, id_dis_l2l, id_dis_s2o, id_dis_s2l

! ---- global clock IDs
integer :: landClock, landFastClock, landSlowClock


! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains


! ============================================================================
subroutine land_model_init &
     (cplr2land, land2cplr, time_init, time, dt_fast, dt_slow)
! initialize land model using grid description file as an input. This routine
! reads land grid boundaries and area of land from a grid description file

! NOTES: theoretically, the grid description file can specify any regular
! rectangular grid for land, not just lon/lat grid. Therefore the variables
! "xbl" and "ybl" in NetCDF grid spec file are not necessarily lon and lat
! boundaries of the grid.
!   However, at this time the module land_properties assumes that grid _is_
! lon/lat and therefore the entire module also have to assume that the land 
! grid is lon/lat.
!   lon/lat grid is also assumed for the diagnostics, but this is probably not
! so critical. 
  type(atmos_land_boundary_type), intent(inout) :: cplr2land ! boundary data
  type(land_data_type)          , intent(inout) :: land2cplr ! boundary data
  type(time_type), intent(in) :: time_init ! initial time of simulation (?)
  type(time_type), intent(in) :: time      ! current time
  type(time_type), intent(in) :: dt_fast   ! fast time step
  type(time_type), intent(in) :: dt_slow   ! slow time step

  ! ---- local vars ----------------------------------------------------------
  integer :: ncid, varid
  integer :: unit, ierr, io
  integer :: id_lon, id_lat, id_band     ! IDs of land diagnostic axes
  logical :: used                        ! return value of send_data diagnostics routine
  integer :: i,j,k
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce, te
  character(len=256) :: restart_file_name
  ! IDs of local clocks
  integer :: landInitClock

  module_is_initialized = .TRUE.

  ! [1] print out version number
  call write_version_number (version, tagname)

  ! initialize land model clocks
  landClock      = mpp_clock_id('Land'               ,CLOCK_FLAG_DEFAULT,CLOCK_COMPONENT)
  landFastClock  = mpp_clock_id('Update-Land-Fast'   ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)
  landSlowClock  = mpp_clock_id('Update-Land-Slow'   ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)
  landInitClock  = mpp_clock_id('Land init'          ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)

  call mpp_clock_begin(landInitClock)

  ! [ ] initialize land debug output
  call land_debug_init()

  ! [ ] initialize tile-specific diagnostics internals
  call tile_diag_init()

  ! [2] read namelists
  ! [2.1] read land model namelist
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=land_model_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_model_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_model_nml)
     call close_file (unit)
  endif
  ! [2.2] read sub-model namelists: then need to be read before initialization
  ! because they can affect the way cover and tiling is initialized on cold start.
  call read_soil_namelist()
  call read_vegn_namelist()
  call read_lake_namelist()
  call read_glac_namelist()
  call read_snow_namelist()
  call read_cana_namelist()

  ! [ ] initialize land state data, including grid geometry and processor decomposition
  call land_data_init(layout, time, dt_fast, dt_slow)
  delta_time  = time_type_to_real(lnd%dt_fast) ! store in a module variable for convenience

  ! calculate land fraction
  allocate( gfrac(size(lnd%glon,1),size(lnd%glon,2)) )
  gfrac = lnd%garea/lnd%gcellarea

  ! [5] initialize tiling
  call get_mosaic_tile_file('INPUT/land.res.nc',restart_file_name,.FALSE.,lnd%domain)
  if(file_exist(restart_file_name)) then
     ! read map of tiles -- retrieve information from 
     call land_cover_warm_start(lnd)
     ! initialize land model data
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,ncid))
     if (nf_inq_varid(ncid,'lwup',varid)==NF_NOERR) &
          call read_tile_data_r0d_fptr(ncid,'lwup',land_lwup_ptr)
     if (nf_inq_varid(ncid,'e_res_1',varid)==NF_NOERR) &
          call read_tile_data_r0d_fptr(ncid,'e_res_1',land_e_res_1_ptr)
     if (nf_inq_varid(ncid,'e_res_2',varid)==NF_NOERR) &
          call read_tile_data_r0d_fptr(ncid,'e_res_2',land_e_res_2_ptr)
     __NF_ASRT__(nf_close(ncid))
  else
     ! initialize map of tiles -- construct it by combining tiles
     ! from component models
     call land_cover_cold_start(lnd)
  endif

  ! [6] initialize land model diagnostics -- must be before *_data_init so that
  ! *_data_init can write static fields if necessary
  call land_diag_init( lnd%glonb, lnd%glatb, lnd%glon, lnd%glat, time, lnd%domain, &
       id_lon, id_lat, id_band )
  ! set the land diagnostic axes ids for the flux exchange
  land2cplr%axes = (/id_lon,id_lat/)
  ! send some static diagnostic fields to output
  if ( id_landarea > 0 ) used = send_data &
       ( id_landarea, lnd%garea (lnd%is:lnd%ie,lnd%js:lnd%je), lnd%time )
  if ( id_landfrac > 0 ) used = send_data &
       ( id_landfrac,     gfrac (lnd%is:lnd%ie,lnd%js:lnd%je), lnd%time )

  ! [7] initialize individual sub-models
  num_species = num_phys + num_c
  if (do_age) num_species = num_species + 1

  call soil_init ( id_lon, id_lat, id_band, use_E_max )
  call vegn_init ( id_lon, id_lat, id_band )
  call lake_init ( id_lon, id_lat )
  call glac_init ( id_lon, id_lat )
  call snow_init ( id_lon, id_lat )
  call cana_init ( id_lon, id_lat )
  call topo_rough_init( lnd%time, &
       lnd%glonb(lnd%is:lnd%ie+1,lnd%js:lnd%je+1), &
       lnd%glatb(lnd%is:lnd%ie+1,lnd%js:lnd%je+1), &
       lnd%domain, id_lon, id_lat)
  call river_init( lnd%glonb(:,1), lnd%glatb(1,:), lnd%time, lnd%dt_fast, lnd%domain, &
       lnd%garea(lnd%is:lnd%ie,lnd%js:lnd%je), gfrac)
  
  ! [8] initialize boundary data
  ! [8.1] allocate storage for the boundary data 
  call realloc_land2cplr ( land2cplr )
  call realloc_cplr2land ( cplr2land )
  ! [8.2] set the land mask to FALSE everywhere -- update_land_bc_fast
  ! will set it to true where necessary
  land2cplr%mask = .FALSE.
  land2cplr%tile_size = 0.0
  ! [8.2] get the current state of the land boundary for the coupler
  ce = first_elmt(lnd%tile_map,                  &
              is=lbound(cplr2land%t_flux,1), &
              js=lbound(cplr2land%t_flux,2)  )
  te = tail_elmt(lnd%tile_map)
  do while(ce /= te)
     ! calculate indices of the current tile in the input arrays;
     ! assume all the cplr2land components have the same lbounds
     call get_elmt_indices(ce,i,j,k)
     ! set this point coordinates as current for debug output
     call set_current_point(i,j,k)
     ! get pointer to current tile
     tile => current_tile(ce)
     ! advance enumerator to the next tile
     ce=next_elmt(ce)

     ! WARNING: the code below will not reproduce across restarts because
     ! values of DTc, Dqc and p_surf are dummy -- they must be either stored
     ! or excluded from the parameter list
     call update_land_bc_fast (tile, i,j,k, &
        p_surf=1e5, land2cplr=land2cplr, is_init=.true.)
  enddo

  ! [8.2.1] update topographic roughness scaling
  call update_land_bc_slow( land2cplr )

  ! mask error checking
  do j=lnd%js,lnd%je
  do i=lnd%is,lnd%ie
     if(gfrac(i,j)>0.neqv.ANY(land2cplr%mask(i,j,:))) then
        call error_mesg('land_model_init','masks are not equal',FATAL)
     endif
  enddo
  enddo

  call mpp_clock_end(landInitClock)

end subroutine land_model_init


! ============================================================================
subroutine land_model_end (cplr2land, land2cplr)
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars
  integer, allocatable :: idx(:)  ! storage for compressed indices
  real   , allocatable :: frac(:) ! storage for tile fractions
  integer, allocatable :: glac(:), lake(:), soil(:), vegn(:)
     ! storage for tile tags
  integer :: i,j,k,n
  integer :: tile_dim_length ! length of tile dimension in output files == 
                             ! global max of number of tiles per gridcell 
  integer :: unit ! netcdf id of the restart file
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  
  module_is_initialized = .FALSE.

  ! write land model restart
  ! [1] count all land tiles and determine the lenght of tile dimension
  ! sufficient for current domain
  n = 0
  tile_dim_length = 0
  do j = lnd%js, lnd%je
  do i = lnd%is, lnd%ie
     k = nitems(lnd%tile_map(i,j))
     n = n + k
     tile_dim_length = max(tile_dim_length,k)
  enddo
  enddo
  ! write(*,*) 'PE:',mpp_pe(),' TILES:',n,' TILE_DIM_LENGTH:',tile_dim_length

  ! [1.1] calculate the tile dimension length by taking the max across all domains
  call mpp_max(tile_dim_length)  
  
  if (n>0) then
     ! [2] gather relevant information
     allocate(idx(n), frac(n), glac(n), lake(n), soil(n), vegn(n))
     ce = first_elmt(lnd%tile_map, lnd%is, lnd%js)
     te = tail_elmt (lnd%tile_map)
     n = 1
     do while (ce/=te)
        tile=>current_tile(ce)
        call get_elmt_indices(ce,i,j,k)
        idx (n) = (k-1)*size(lnd%garea,1)*size(lnd%garea,2) + (j-1)*size(lnd%garea,1) + (i-1)
        frac(n) = tile%frac
        call get_tile_tags(tile,glac=glac(n),lake=lake(n),soil=soil(n),vegn=vegn(n))
        ce=next_elmt(ce)
        n = n+1
     end do
   
     ! [3] create tile output file
     call error_mesg('land_model_end','writing NetCDF restart',NOTE)
     call create_tile_out_file(unit,'RESTART/land.res.nc', &
          lnd%glon*180.0/PI, lnd%glat*180/PI, idx, tile_dim_length)
     
     ! [4] write out fields
     call write_tile_data(unit,'frac',frac,'fractional area of tile')
     call write_tile_data(unit,'glac',glac,'tag of glacier tiles')
     call write_tile_data(unit,'lake',lake,'tag of lake tiles')
     call write_tile_data(unit,'soil',soil,'tag of soil tiles')
     call write_tile_data(unit,'vegn',vegn,'tag of vegetation tiles')
     ! write the upward long-wave flux 
     call write_tile_data_r0d_fptr(unit,'lwup',land_lwup_ptr,'upward long-wave flux')
     ! write energy residuals
     call write_tile_data_r0d_fptr(unit,'e_res_1',land_e_res_1_ptr,&
          'energy residual in canopy air energy balance equation', 'W/m2')
     call write_tile_data_r0d_fptr(unit,'e_res_2',land_e_res_2_ptr,&
          'energy residual in canopy energy balance equation', 'W/m2')
     

     ! [5] close file
     __NF_ASRT__(nf_close(unit))

     ! [6] clean up
     deallocate(idx,frac,glac, lake, soil, vegn)
  else
     ! if total number of tiles in the current domain is zero, do not
     ! even bother to create tiled restart for it
  endif

  ! we still want to call the *_end procedures for component models, even
  ! if the number of tiles in this domain is zero, in case they are doing 
  ! something else besides saving the restart, of if they want to save
  ! restart anyway
  call glac_end (tile_dim_length)
  call lake_end (tile_dim_length)
  call soil_end (tile_dim_length)
  call snow_end (tile_dim_length)
  call vegn_end (tile_dim_length)
  call cana_end (tile_dim_length)
  call topo_rough_end()
  call river_end()

  deallocate(gfrac)
  call dealloc_land2cplr(land2cplr)
  call dealloc_cplr2land(cplr2land)

  call tile_diag_end()

  ! deallocate tiles
  call land_data_end()

  ! finish up the land debugging diagnostics
  call land_debug_end
  
end subroutine land_model_end


! ============================================================================
subroutine land_cover_cold_start(lnd)
  type(land_state_type), intent(inout) :: lnd

  ! ---- local vars
  real, dimension(:,:,:), pointer :: &
       glac, soil, lake, vegn ! arrays of fractions for respective sub-models
  logical, allocatable :: gmask(:,:), valid_data(:,:)
  integer :: i,j,i1,j1

  ! NOTES: lnd%area must be global
  allocate(gmask(size(lnd%garea,1),size(lnd%garea,2)))
  allocate(valid_data(size(lnd%garea,1),size(lnd%garea,2)))

  ! calculate the global land mask
  gmask = lnd%garea > 0

  ! get the global maps of fractional covers for each of the sub-models
  glac=>glac_cover_cold_start(gmask,lnd%glonb,lnd%glatb)
  lake=>lake_cover_cold_start(gmask,lnd%glonb,lnd%glatb)
  soil=>soil_cover_cold_start(gmask,lnd%glonb,lnd%glatb)
  vegn=>vegn_cover_cold_start(gmask,lnd%glonb,lnd%glatb)

  ! reconcile ground fractions with the land mask within compute domain
  valid_data = gmask.and.(sum(glac,3)+sum(lake,3)+sum(soil,3)>0)
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     if(.not.gmask(i,j)) cycle ! skip ocean points
     if(valid_data(i,j)) cycle ! don't need to do anything with valid points

     call nearest(valid_data,lnd%glon(:,1),lnd%glat(1,:),lnd%glon(i,j),lnd%glat(i,j),i1,j1)
     glac(i,j,:) = glac(i1,j1,:)
     lake(i,j,:) = lake(i1,j1,:)
     soil(i,j,:) = soil(i1,j1,:)
     call set_current_point(i,j,1)
     if(is_watch_point())then
        write(*,*)'###### land_cover_cold_start: reconciling ground with the land mask #####'
        __DEBUG2__(i,j)
        __DEBUG2__(lnd%glon(i,j),lnd%glat(i,j))
        __DEBUG2__(i1,j1)
     endif
  enddo
  enddo

  ! reconcile vegetation fractions with the land mask within compute domain
  valid_data = sum(vegn,3) > 0
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     if(.not.gmask(i,j)) cycle ! skip ocean points
     if(valid_data(i,j)) cycle ! don't need to do anything with valid points
     if(sum(glac(i,j,:))+sum(lake(i,j,:))>=1) &
          cycle                ! skip points fully covered by glaciers or lakes
     call nearest(valid_data,lnd%glon(:,1),lnd%glat(1,:),lnd%glon(i,j),lnd%glat(i,j),i1,j1)
     vegn(i,j,:) = vegn(i1,j1,:)
     call set_current_point(i,j,1)
     if(is_watch_point())then
        write(*,*)'###### land_cover_cold_start: reconciling vegetation with the land mask #####'
        __DEBUG2__(i,j)
        __DEBUG2__(lnd%glon(i,j),lnd%glat(i,j))
        __DEBUG2__(i1,j1)
     endif
  enddo
  enddo
  
  ! create tiles
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     if(.not.gmask(i,j)) cycle ! skip ocean points
     call set_current_point(i,j,1)
     call land_cover_cold_start_0d &
          (lnd%tile_map(i,j),glac(i,j,:),lake(i,j,:),soil(i,j,:),vegn(i,j,:))
     if(nitems(lnd%tile_map(i,j))==0) then
        call error_mesg('land_cover_cold_start',&
             'No tiles were created for a valid land point', FATAL)
     endif
  enddo
  enddo

  deallocate(glac,lake,soil,vegn)
  
end subroutine land_cover_cold_start

! ============================================================================
subroutine land_cover_cold_start_0d (set,glac0,lake0,soil0,vegn0)
  type(land_tile_list_type), intent(inout) :: set 
  real, dimension(:)       , intent(in) :: &
       glac0,lake0,soil0,vegn0 ! fractions of area

  ! ---- local vars
  real :: glac(size(glac0(:))), lake(size(lake0(:))), &
          soil(size(soil0(:))), vegn(size(vegn0(:)))
  type(land_tile_type), pointer :: tile
  integer :: i,j
  real :: factor ! normalizing factor for the tile areas
  real :: frac
  type(land_tile_enum_type) :: first_non_vegn ! position of first non-vegetated tile in the list

  glac = glac0; lake = lake0; soil = soil0; vegn = vegn0
  if (sum(glac)>1) &
       glac=glac/sum(glac)
  if (sum(lake)+sum(glac)>1)&
       lake = lake*(1-sum(glac))/sum(lake)
  if (sum(soil)+sum(glac)+sum(lake)>1)&
       soil = soil*(1-sum(lake)-sum(glac))/sum(soil)
  ! make sure that the sum of the fractions of the soil, lake, and glaciers are 
  ! either one or zero
  factor = sum(soil)+sum(glac)+sum(lake)
  if(factor>0)then
     glac = glac/factor
     lake = lake/factor
     soil = soil/factor
  endif
  if(is_watch_point()) then
     write(*,*)'#### land_cover_cold_start_0d ####'
     __DEBUG1__(glac0)
     __DEBUG1__(lake0)
     __DEBUG1__(soil0)
     __DEBUG1__(vegn0)
     __DEBUG1__(factor)
     __DEBUG1__(glac)
     __DEBUG1__(lake)
     __DEBUG1__(soil)
     __DEBUG1__(vegn)
  endif

  do i = 1,size(glac)
     if (glac(i)>0) then
        tile => new_land_tile(frac=glac(i),glac=i)
        call insert(tile,set)
        if(is_watch_point()) then
           write(*,*)'created glac tile', i, ' frac=', glac(i)
        endif
     endif
  enddo
  do i = 1,size(lake)
     if (lake(i)>0) then
        tile => new_land_tile(frac=lake(i),lake=i)
        call insert(tile,set)
        if(is_watch_point()) then
           write(*,*)'created lake tile',i,' frac=',lake(i)
        endif
     endif
  enddo

  factor = sum(soil)*sum(vegn)
  if (factor/=0) factor = 1/factor
  factor = factor*(1-sum(glac)-sum(lake))
  ! vegetation tiles, if any, are inserted in front of non-vegetated tiles;
  ! this really doesn't matter except for the static vegetation override
  ! case with the data saved by lm3v -- there the vegetation tiles are
  ! in front, so it works more consistently where lad2 has more than 
  ! one tile (e.g. glac/soil or lake/soil), if lad2 vegetation tiles are 
  ! also in front of the list. 
  first_non_vegn=first_elmt(set)
  do i = 1,size(soil)
  do j = 1,size(vegn)
     frac = soil(i)*vegn(j)*factor
     if(frac>0) then
        tile  => new_land_tile(frac=frac,soil=i,vegn=j)
        call insert(tile,first_non_vegn)
        if(is_watch_point()) then
           write(*,*)'created soil tile', i, j, ' frac=', frac
        endif
     endif
  enddo
  enddo

end subroutine land_cover_cold_start_0d

! ============================================================================
! reads the land restart file and restores the tiling structure from this file
subroutine land_cover_warm_start ( lnd )
  type(land_state_type), intent(inout) :: lnd
  
  ! ---- local vars
  integer, allocatable :: idx(:) ! compressed tile index
  integer, allocatable :: glac(:), lake(:), soil(:), snow(:), cana(:), vegn(:) ! tile tags
  real,    allocatable :: frac(:) ! fraction of land covered by tile
  integer :: ncid ! unit number of the input file
  integer :: n    ! total number of tiles in the input size
  integer :: dimids(1) ! id of tile dimension
  character(NF_MAX_NAME) :: tile_dim_name ! name of the tile dimension and respective variable
  integer :: i,j,k,it
  type(land_tile_type), pointer :: tile;
  character(len=256) :: restart_file_name
  
  call get_mosaic_tile_file('INPUT/land.res.nc',restart_file_name,.FALSE.,lnd%domain)
  __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,ncid))
  ! allocate the input data
  __NF_ASRT__(nfu_inq_var(ncid,'frac',varsize=n,dimids=dimids))
  allocate(idx(n))
  allocate(glac(n),lake(n),soil(n),snow(n),cana(n),vegn(n))
  allocate(frac(n))
  ! get the name of the fist (and only) dimension of the variable 'frac' -- this
  ! is supposed to be the compressed dimension, and associated variable will
  ! hold the compressed indices
  __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),name=tile_dim_name))
  ! read the compressed tile indices
  __NF_ASRT__(nfu_get_var_int(ncid,tile_dim_name,idx))
  ! read input data -- fractions and tags
  __NF_ASRT__(nfu_get_var_double(ncid,'frac',frac))
  __NF_ASRT__(nfu_get_var_int(ncid,'glac',glac))
  __NF_ASRT__(nfu_get_var_int(ncid,'lake',lake))
  __NF_ASRT__(nfu_get_var_int(ncid,'soil',soil))
  __NF_ASRT__(nfu_get_var_int(ncid,'vegn',vegn))
  
  ! create tiles
  do it = 1,n
     k = idx(it)
     i = modulo(k,size(lnd%garea,1))+1; k = k/size(lnd%garea,1)
     j = modulo(k,size(lnd%garea,2))+1; k = k/size(lnd%garea,2)
     k = k + 1
     if (i<lnd%is.or.i>lnd%ie) cycle
     if (j<lnd%js.or.j>lnd%je) cycle
     ! the size of the tile set at the point (i,j) must be equal to k
     tile=>new_land_tile(frac=frac(it),&
              glac=glac(it),lake=lake(it),soil=soil(it),vegn=vegn(it))
     call insert(tile,lnd%tile_map(i,j))
  enddo
  __NF_ASRT__(nf_close(ncid))
end subroutine


! ============================================================================
subroutine update_land_model_fast ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(in)    :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local constants
  real,    parameter :: cana_co2 = 350.0e-6 ! co2 concentration in canopy air, mol/mol
  ! indices of variables and equations for implicit time stepping solution :
  integer, parameter :: iqc=1, iTc=2, iTv=3, iwl=4, iwf=5

  ! ---- local vars 
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  
  real :: A(5,5),B1(5),B0(5) ! implicit equation matrix and right-hand side vectors
  real :: A00(5,5),B10(5),B00(5) ! copy of the above, only for debugging
  integer :: indx(5) ! permutation vector
  ! linearization coefficients of various fluxes between components of land
  ! surface scheme
  real :: &
       G0,    DGDTg,  &  ! ground heat flux 
       Ha0,   DHaDTc, &  ! sensible heat flux from the canopy air to the atmosphere 
       Ea0,   DEaDqc, &  ! water vapor flux from canopy air to the atmosphere
       Hv0,   DHvDTv,   DHvDTc, & ! sens heat flux from vegetation
       Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
       Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
       Esi0,  DEsiDTv,  DEsiDqc,  DEsiDwl,  DEsiDwf, & ! sublimation of intercepted snow
       Hg0,   DHgDTg,   DHgDTc, & ! linearization of the sensible heat flux from ground
       Eg0,   DEgDTg,   DEgDqc, & ! linearization of evaporation from ground
       flwv0,  DflwvDTg,  DflwvDTv,& ! linearization of net LW radiation to the canopy
       flwg0,  DflwgDTg,  DflwgDTv,& ! linearization of net LW radiation to the canopy
       vegn_drip_l, vegn_drip_s, & ! drip rate of water and snow, respectively, kg/(m2 s)
       vegn_lai
  
  ! increments of respective variables over time step, results of the implicit
  ! time step:
  real :: delta_qc, delta_Tc, delta_Tv, delta_wl, delta_ws, delta_Tg
  real :: flwg ! updated value of long-wave ground energy balance
  real :: vegn_emis_lw, surf_emis_lw ! emissivities of ground and surface
  real :: vegn_emsn,    surf_emsn    ! emission by vegetation and surface, respectively
  real :: denom ! denominator in the LW radiative balance calculations
  real :: sum0, sum1

  real :: TEMP1, &
       grnd_T, & ! ground temperature
       vegn_T, & ! vegetation (canopy) temperature
       cana_T, & ! canopy air temperature
       soil_uptake_T, & ! average temperatue of water taken up by the vegetation
       vegn_Wl,  vegn_Ws, & ! water and snow mass of the canopy
       vegn_ifrac, & ! intercepted fraction of liquid or frozen precipitation
       vegn_hcap,      & ! vegetation heat capacity, including intercepted water and snow
       hlv_Tv, hlv_Tu, & ! latent heat of vaporization at vegn and uptake temperatures, respectively 
       hls_Tv, &         ! latent heat of sublimation at vegn temperature
       grnd_rh, & 
       grnd_liq, grnd_ice, grnd_subl, &
       grnd_tf, &  ! temperature of freezing on the ground
       grnd_latent, &
       grnd_flux, &
       grnd_E_max, &
       soil_beta, &
       RSv(NBANDS), & ! net short-wave radiation balance of the canopy, W/m2
       con_g_h, con_g_v, & ! turbulent cond. between ground and canopy air, for heat and vapor respectively
       snow_area, &
       cana_q, & ! specific humidity of canopy air
       fswg, evapg, sensg, &
       subs_G, Mg_imp, snow_G_Z, snow_G_TZ, &
       vegn_ovfl_l,  vegn_ovfl_s,  & ! overflow of liquid and solid water from the canopy
       vegn_ovfl_Hl, vegn_ovfl_Hs, & ! heat flux from canopy due to overflow
       delta_fprec, & ! correction of below-canopy solid precip in case it's average T > tfreeze 

       ISa_dn_dir(NBANDS), & ! downward direct sw radiation at the top of the canopy
       ISa_dn_dif(NBANDS), & ! downward diffuse sw radiation at the top of the canopy
       ILa_dn,             & ! downward lw radiation at the to of the canopy
       vegn_flw,vegn_sens,snow_sens,snow_levap,snow_fevap,snow_hadvec,snow_melt,&
       snow_lprec, snow_hlprec,snow_lrunf, precip_s,vegn_levap,vegn_fevap,vegn_uptk,&
       vegn_fsw, vegn_hadvec,vegn_melt,vegn_lprec,vegn_fprec,vegn_hlprec,vegn_hfprec,vegn_LMASS,&
       vegn_FMASS,vegn_HEAT, precip_l,precip_T,snow_fsw,snow_flw,snow_frunf,snow_hlrunf,&
       snow_hfrunf, snow_LMASS,snow_FMASS,snow_HEAT,subs_fsw,subs_flw,subs_sens,&
       subs_DT, subs_M_imp, subs_evap, snow_Tbot, snow_Cbot, subs_levap,&
       subs_fevap,subs_hadvec,subs_melt,subs_lrunf,subs_hlrunf,subs_LMASS,subs_FMASS,&
       subs_HEAT,subs_Ttop,subs_Ctop, subs_subl
  real :: soil_water_supply ! supply of water to roots, per unit active root biomass, kg/m2
  real :: snow_T, snow_rh, snow_liq, snow_ice, snow_subl
  integer :: i, j, k, i_species
  integer :: ii, jj ! indices for debug output
  integer :: ierr
  logical :: bone_dry, conserve_glacier_mass, snow_active
  real :: subs_z0m, subs_z0s, snow_z0m, snow_z0s, grnd_z0s
  real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je) :: &
       heat_pos,      & ! accumulated heat runoff per grid cell
       sol_pos,       & ! accumulated snow runoff per grid cell
       tot_pos          ! accumulated total runoff per grid cell
  real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je,num_species) :: &
       c_flux           ! flux of tracers into the rivers
  real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je) :: &
       discharge2o_l      ! discharge of liquid water to ocean
  real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je) :: &
       discharge2l_l      ! discharge of liquid water to land
  real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je,num_species) :: &
       discharge2o_c      ! discharge of tracers from the rivers to ocean
  real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je,num_species) :: &
       discharge2l_c      ! discharge of tracers from the rivers to land
  real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je) :: &
       discharge_work
  logical :: used         ! return value of send_data diagnostics routine

  ! start clocks
  call mpp_clock_begin(landClock)
  call mpp_clock_begin(landFastClock)

  ! clear the runoff values, for accumulation over the tiles
  heat_pos = 0 ; sol_pos = 0 ; tot_pos = 0 ; c_flux = 0

  ! initialize current tile enumerator
  ce = first_elmt(lnd%tile_map,                  &
              is=lbound(cplr2land%t_flux,1), &
              js=lbound(cplr2land%t_flux,2)  )
  ! get the end marker (end tile enumerator)
  te = tail_elmt(lnd%tile_map)
  ! main tile loop
  do while(ce /= te)
     ! calculate indices of the current tile in the input arrays;
     ! assume all the cplr2land components have the same lbounds
     call get_elmt_indices(ce,i,j,k)
     ! set this point coordinates as current for debug output
     call set_current_point(i,j,k)
     ! get pointer to current tile
     tile => current_tile(ce)
     ! advance enumerator to the next tile
     ce=next_elmt(ce)

     ! get data from atmosphere
     precip_l = cplr2land%lprec(i,j,k)
     precip_s = cplr2land%fprec(i,j,k)
     precip_T = cplr2land%tprec(i,j,k)
     Ha0    =  cplr2land%t_flux(i,j,k)
     DHaDTc =  cplr2land%dhdt  (i,j,k)
#ifdef LAND_BND_TRACERS
     Ea0    = cplr2land%tr_flux(i,j,k, lnd%isphum)
     DEaDqc = cplr2land%dfdtr(i,j,k, lnd%isphum)
#else
     Ea0    = cplr2land%q_flux(i,j,k)
     DEaDqc = cplr2land%dedq(i,j,k)
#endif
     ISa_dn_dir(BAND_VIS) = cplr2land%sw_flux_down_vis_dir(i,j,k)
     ISa_dn_dir(BAND_NIR) = cplr2land%sw_flux_down_total_dir(i,j,k)&
                           -cplr2land%sw_flux_down_vis_dir(i,j,k)
     ISa_dn_dif(BAND_VIS) = cplr2land%sw_flux_down_vis_dif(i,j,k)
     ISa_dn_dif(BAND_NIR) = cplr2land%sw_flux_down_total_dif(i,j,k)&
                           -cplr2land%sw_flux_down_vis_dif(i,j,k)
     ILa_dn               = cplr2land%lwdn_flux(i,j,k)

     ! *** This bone_dry code has _NOT_ been tested successfully against LaD1 code.
     ! It is left here as a starting point in case I decide to pursue it again.
     bone_dry    = .FALSE.
     IF (LM2 .and. use_bone_dry) THEN                  ! (Should get 'temp1' from soil module)
        TEMP1 = max(0.,sum(tile%snow%prog%ws))/3600.0
        if (associated(tile%glac)) then
           TEMP1 = TEMP1 + 1.e16
        else if (associated(tile%lake)) then
           TEMP1 = TEMP1+sum(tile%lake%prog%wl)/3600.0
        else if (associated(tile%soil)) then
           do i = 1, min(size(tile%soil%prog),nz_lm2(i,j,k))
              TEMP1 = TEMP1 + tile%soil%prog(i)%wl/3600.
           enddo
           TEMP1 = TEMP1-150.*zr_lm2(i,j,k)/3600.
        else
           call error_mesg('update_land_model_fast','none of the surface tiles exists',FATAL)
        endif
        TEMP1=max(TEMP1,0.0)
        if (Ea0 > TEMP1) then
           bone_dry     = .true.
           Ea0  = TEMP1
           DEaDqc    = 0.
        end if
     ENDIF

     soil_uptake_T = tfreeze ! just to avoid using un-initialized values
     if (associated(tile%glac)) then
        call glac_step_1 ( tile%glac, &
             grnd_T, grnd_rh, grnd_liq, grnd_ice, grnd_subl, grnd_tf, &
             snow_G_Z, snow_G_TZ, conserve_glacier_mass  )
        grnd_E_max = HUGE(grnd_E_max)
     else if (associated(tile%lake)) then
        call lake_step_1 ( tile%lake, &
             grnd_T, grnd_rh, grnd_liq, grnd_ice, grnd_subl, grnd_tf, &
             snow_G_Z, snow_G_TZ)
        grnd_E_max = HUGE(grnd_E_max)
     else if (associated(tile%soil)) then
        call soil_step_1 ( tile%soil, tile%vegn, &
             grnd_T, soil_uptake_T, soil_beta, soil_water_supply, grnd_E_max, &
             grnd_rh, grnd_liq, grnd_ice, grnd_subl, grnd_tf, &
             snow_G_Z, snow_G_TZ)
     else
        call error_mesg('update_land_model_fast','none of the surface tiles exists',FATAL)
     endif

     subs_subl = grnd_subl

     call snow_step_1 ( tile%snow, snow_G_Z, snow_G_TZ, &
          snow_active, snow_T, snow_rh, snow_liq, snow_ice, &
          snow_subl, snow_area, G0, DGDTg )
     if (snow_active) then
        grnd_T    = snow_T;   grnd_rh   = snow_rh;   grnd_liq  = snow_liq
        grnd_ice  = snow_ice; grnd_subl = snow_subl; grnd_tf   = tfreeze
     endif
     grnd_latent = hlv + hlf*grnd_subl

! since con_v_v is no longer available in this subroutine, the code below is not
! going to work. Possibly introduce another option in the turbulence calculations? 
!!$     IF (LM2.and.associated(tile%soil)) THEN
!!$!        if (snow_active.and.tile%vegn%pars%rs_min>0.0) then   ! *** CLEAN THIS UP
!!$        if (snow_active.and.vegn_data_rs_min(tile%vegn)>0.0) then   ! *** CLEAN THIS UP
!!$           con_g_v = 1.0/vegn_data_rs_min(tile%vegn)
!!$           con_v_v = con_v_v * con_fac_small
!!$        else if (snow_active) then
!!$          con_g_v = con_g_v * con_fac_large
!!$          con_v_v = con_v_v * con_fac_small
!!$        else
!!$          con_g_v = con_g_v * con_fac_small
!!$          con_v_v = con_v_v * con_fac_large
!!$        endif
!!$        con_g_h = con_g_h * con_fac_large
!!$        con_v_h = con_v_h * con_fac_large
!!$     ENDIF

     call cana_state(tile%cana, cana_T, cana_q)

     if (associated(tile%vegn)) then
     ! Calculate net short-wave radiation input to the vegetation
        RSv    = tile%Sv_dir*ISa_dn_dir + tile%Sv_dif*ISa_dn_dif
        call soil_diffusion(tile%soil, subs_z0s, subs_z0m)
        call snow_diffusion(tile%snow, snow_z0s, snow_z0m)
        grnd_z0s = exp( (1-snow_area)*log(subs_z0s) + snow_area*log(snow_z0s))
        
        call vegn_step_1 ( tile%vegn, tile%diag, &
           cplr2land%p_surf(i,j,k), &
           cplr2land%ustar (i,j,k), &
           cplr2land%drag_q(i,j,k), &
           ISa_dn_dir+ISa_dn_dif, RSv, precip_l, precip_s, &
           tile%land_d, tile%land_z0s, tile%land_z0m, grnd_z0s, & 
           soil_beta, soil_water_supply,  soil_uptake_T,&
           cana_T, cana_q, cana_co2, &
           ! output
           con_g_h, con_g_v, &
           vegn_T, vegn_Wl, vegn_Ws, & ! temperature, water and snow mass on the canopy
           vegn_ifrac, vegn_lai, &
           vegn_drip_l, vegn_drip_s,& 
           vegn_hcap, & ! total vegetation heat capacity (including intercepted water/snow)
           Hv0,   DHvDTv,   DHvDTc,            & 
           Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & 
           Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & 
           Esi0,  DEsiDTv,  DEsiDqc,  DEsiDwl,  DEsiDwf  ) 
     else
        RSv    = 0
        con_g_h = con_fac_large ; con_g_v = con_fac_large
        if(associated(tile%glac).and.conserve_glacier_mass.and..not.snow_active) &
             con_g_v = con_fac_small
        vegn_T  = cana_T ; vegn_Wl = 0 ; vegn_Ws = 0
        vegn_ifrac  = 0 ; vegn_lai    = 0
        vegn_drip_l = 0 ; vegn_drip_s = 0
        vegn_hcap = 1.0
        Hv0 =0;  DHvDTv =0;  DHvDTc=0;
        Et0 =0;  DEtDTv =0;  DEtDqc=0;   DEtDwl=0;   DEtDwf=0
        Eli0=0;  DEliDTv=0;  DEliDqc=0;  DEliDwl=0;  DEliDwf=0 
        Esi0=0;  DEsiDTv=0;  DEsiDqc=0;  DEsiDwl=0;  DEsiDwf=0
     endif
     ! calculate net shortwave for ground and canopy
     fswg     = SUM(tile%Sg_dir*ISa_dn_dir + tile%Sg_dif*ISa_dn_dif)
     vegn_fsw = SUM(RSv)
     
     call cana_step_1 (tile%cana, cplr2land%p_surf(i,j,k), con_g_h, con_g_v,   &
          grnd_t, grnd_rh, &
          Hg0,  DHgDTg, DHgDTc, Eg0, DEgDTg, DEgDqc)


! [X.X] using long-wave optical properties, calculate the explicit long-wave 
!       radiative balances and their derivatives w.r.t. temperatures
     vegn_emis_lw = 1 - tile%vegn_refl_lw - tile%vegn_tran_lw
     surf_emis_lw = 1 - tile%surf_refl_lw

     denom = 1-tile%vegn_refl_lw*tile%surf_refl_lw

     vegn_emsn = vegn_emis_lw * stefan * vegn_T**4
     surf_emsn = surf_emis_lw * stefan * grnd_T**4

     flwv0 = ILa_dn * vegn_emis_lw*(1+tile%vegn_refl_lw*tile%surf_refl_lw/denom) &
          + vegn_emsn * (tile%surf_refl_lw*vegn_emis_lw/denom-2) &
          + surf_emsn * vegn_emis_lw/denom
     DflwvDTg = vegn_emis_lw/denom                      * surf_emis_lw * stefan * 4 * grnd_T**3
     DflwvDTv = (tile%surf_refl_lw*vegn_emis_lw/denom-2)* vegn_emis_lw * stefan * 4 * vegn_T**3

     flwg0 = (ILa_dn*tile%vegn_tran_lw + vegn_emsn)*(1-tile%surf_refl_lw)/denom &
          - surf_emsn*(1-tile%vegn_refl_lw)/denom
     DflwgDTg = -(1-tile%vegn_refl_lw)/denom * surf_emis_lw * stefan * 4 * grnd_T**3
     DflwgDTv =  (1-tile%surf_refl_lw)/denom * vegn_emis_lw * stefan * 4 * vegn_T**3

! [X.0] calculate the latent heats of vaporization for appropriate temperatures
     hlv_Tv = hlv       - (cpw-clw)*tfreeze + cpw*vegn_T
     hls_Tv = hlv + hlf - (cpw-csw)*tfreeze + cpw*vegn_T
     hlv_Tu = hlv       - (cpw-clw)*tfreeze + cpw*vegn_T - clw*soil_uptake_T
     ! what about ground? should grnd_latent depend on temperature?
     if(is_watch_point()) then
        write(*,*)'#### input data for the matrix ####'
        __DEBUG1__(delta_time)
        __DEBUG3__(vegn_T,vegn_Wl,vegn_Ws)
        __DEBUG2__(grnd_T,grnd_rh)
        __DEBUG2__(cana_T,cana_q)
        __DEBUG2__(vegn_emis_lw,surf_emis_lw)
        __DEBUG2__(vegn_emsn,surf_emsn)
        __DEBUG3__(precip_l, vegn_drip_l, precip_T)
        __DEBUG2__(precip_s, vegn_drip_s)
        __DEBUG2__(vegn_ifrac, vegn_lai)
        __DEBUG1__(ILa_dn)
        __DEBUG2__(ISa_dn_dir(1),ISa_dn_dir(2))
        __DEBUG2__(ISa_dn_dif(1),ISa_dn_dif(2))
        __DEBUG2__(fswg, vegn_fsw)
        __DEBUG1__(vegn_hcap)
        __DEBUG3__(hlv_Tv, hlv_Tu, hls_Tv)
        __DEBUG2__(G0, DGDTg)
        __DEBUG2__(Ha0, DHaDTc)
        __DEBUG2__(Ea0, DEaDqc)
        __DEBUG3__(Hv0, DHvDTv, DHvDTc)
        __DEBUG5__(Et0,  DEtDTv,  DEtDqc,  DEtDwl,  DEtDwf)
        __DEBUG5__(Eli0, DEliDTv, DEliDqc, DEliDwl, DEliDwf)
        __DEBUG5__(Esi0, DEsiDTv, DEsiDqc, DEsiDwl, DEsiDwf)
        __DEBUG3__(Hg0, DHgDTg, DHgDTc)
        __DEBUG3__(Eg0, DEgDTg, DEgDqc)
        __DEBUG3__(flwv0, DflwvDTg, DflwvDTv)
        __DEBUG3__(flwg0, DflwgDTg, DflwgDTv)
        __DEBUG2__(tile%e_res_1,tile%e_res_2)
     endif

! [X.1] form the system of equations for implicit scheme, such that A*X = B1*delta_Tg+B0
! [X.1.1] equation of canopy air mass balance
     A(iqc,iqc) = canopy_air_mass/delta_time-DEtDqc-DEliDqc-DEsiDqc-DEgDqc+DEaDqc
     A(iqc,iTc) = 0
     A(iqc,iTv) = -DEtDTv-DEliDTv-DEsiDTv
     A(iqc,iwl) = -DEtDwl-DEliDwl-DEsiDwl
     A(iqc,iwf) = -DEtDwf-DEliDwf-DEsiDwf
     B0(iqc)  = Esi0+Eli0+Et0+Eg0-Ea0
     B1(iqc)  = DEgDTg
! [X.1.2] equation of canopy air energy balance
     A(iTc,iqc) = canopy_air_mass*cpw*cana_T/delta_time &
          - cpw*vegn_T*(DEtDqc+DEliDqc+DEsiDqc) - cpw*grnd_T*DEgDqc + cpw*cana_T*DEaDqc
     A(iTc,iTc) = canopy_air_mass*cp_air/delta_time-DHvDTc-DHgDTc+DHaDTc
     A(iTc,iTv) = -DHvDTv-cpw*vegn_T*(DEtDTv+DEliDTv+DEsiDTv)
     A(iTc,iwl) = -cpw*vegn_T*(DEtDwl+DEliDwl+DEsiDwl)
     A(iTc,iwf) = -cpw*vegn_T*(DEtDwf+DEliDwf+DEsiDwf)
     B0(iTc)  = Hv0 + Hg0 - Ha0 + cpw*(vegn_T*(Et0+Eli0+Esi0)+grnd_T*Eg0-cana_T*Ea0) - tile%e_res_1
     B1(iTc)  = DHgDTg + cpw*grnd_T*DEgDTg
! [X.1.3] equation of canopy energy balance
     A(iTv,iqc) = hlv_Tu*DEtDqc + hlv_Tv*DEliDqc + hls_Tv*DEsiDqc
     A(iTv,iTc) = DHvDTc
     A(iTv,iTv) = vegn_hcap/delta_time-DflwvDTv + DHvDTv + &
          hlv_Tu*DEtDTv + hlv_Tv*DEliDTv + hls_Tv*DEsiDTv
     A(iTv,iwl) = clw*vegn_T/delta_time + hlv_Tu*DEtDwl + hlv_Tv*DEliDwl + hls_Tv*DEsiDwl
     A(iTv,iwf) = csw*vegn_T/delta_time + hlv_Tu*DEtDwf + hlv_Tv*DEliDwf + hls_Tv*DEsiDwf
     B0(iTv)  = vegn_fsw + flwv0 - Hv0 - hlv_Tu*Et0 - Hlv_Tv*Eli0 - hls_Tv*Esi0 &
          + clw*precip_l*vegn_ifrac*precip_T + csw*precip_s*vegn_ifrac*precip_T &
          - clw*vegn_drip_l*vegn_T - csw*vegn_drip_s*vegn_T - tile%e_res_2
     B1(iTv)  = DflwvDTg
! [X.1.4] equation of intercepted liquid water mass balance
     A(iwl,iqc) = DEliDqc
     A(iwl,iTc) = 0
     A(iwl,iTv) = DEliDTv
     A(iwl,iwl) = 1.0/delta_time + DEliDwl
     A(iwl,iwf) = DEliDwf
     B0(iwl)  = -Eli0 + precip_l*vegn_ifrac - vegn_drip_l
     B1(iwl)  = 0
! [X.1.5] equation of intercepted frozen water mass balance
     A(iwf,iqc) = DEsiDqc
     A(iwf,iTc) = 0
     A(iwf,iTv) = DEsiDTv
     A(iwf,iwl) = DEsiDwl
     A(iwf,iwf) = 1.0/delta_time + DEsiDwf
     B0(iwf)  = -Esi0 + precip_s*vegn_ifrac - vegn_drip_s
     B1(iwf)  = 0
! [X.1.6] if LAI becomes zero (and, therefore, all fluxes from vegetation and their 
! derivatives must be zero too) we get a degenerate case. Still, the drip may be non-zero
! because some water may remain from before leaf drop, and non-zero energy residual can be
! carried over from the previous time step.
! To prevent temperature from going haywire in those cases, we simply replace the equations 
! of canopy energy and mass balance with the following:
! vegn_T + delta_Tv = cana_T + delta_Tc
! delta_Wl = -vegn_Wl 
! delta_Ws = -vegn_Ws 
     if(vegn_lai==0) then
        ! vegn_T + delta_Tv = cana_T + delta_Tc
        A(iTv,:)   = 0
        A(iTv,iTc) = -1
        A(iTv,iTv) = +1
        B0(iTv) = cana_T - vegn_T
        B1(iTv) = 0
        ! delta_Wl = -vegn_Wl 
        A(iwl,:)   = 0
        A(iwl,iwl) = 1
        B0(iwl) = -vegn_Wl
        B1(iwl) = 0
        ! delta_Ws = -vegn_Ws 
        A(iwf,:)   = 0
        A(iwf,iwf) = 1
        B0(iwf) = -vegn_Ws
        B0(iwf) = 0
     endif



     if(is_watch_point()) then
        write(*,*)'#### A, B0, B1 ####'
        do ii = 1, size(A,1)
           write(*,'(99g)')(A(ii,jj),jj=1,size(A,2)),B0(ii),B1(ii)
        enddo
     endif

     A00 = A
     B00 = B0
     B10 = B1

! [X.2] solve the system for free terms and delta_Tg terms, getting
!       linear equation for delta_Tg
     call ludcmp(A,indx, ierr)
     if (ierr/=0)&
          write(*,*) 'Matrix is singular',i,j,k
     call lubksb(A,indx,B0)
     call lubksb(A,indx,B1)

     if(is_watch_point()) then
        write(*,*)'#### solution: B0, B1 ####'
        do ii = 1, size(A,1)
           __DEBUG2__(B0(ii),B1(ii))
        enddo
!!$        write(*,*)'#### solution check ####'
!!$        do ii = 1, size(A,1)
!!$           sum0 = 0; sum1 = 0;
!!$           do jj = 1, size(A,2)
!!$              sum0 = sum0 + A00(ii,jj)*B0(jj)
!!$              sum1 = sum1 + A00(ii,jj)*B1(jj)
!!$           enddo
!!$           write(*,'(99g)')sum0-B00(ii),sum1-B10(ii)
!!$        enddo
     endif
! the result of this solution is a set of expressions for delta_xx in terms
! of delta_Tg: delta_xx(i) = B0(i) + B1(i)*delta_Tg. Note that A, B0, and B1
! are destroyed in the process: A is replaced with LU-decomposition, and
! B0, B1 are replaced with solutions

     ! solve the non-linear equation for energy balance at the surface.
     if(snow_active.or..not.use_E_max) & ! just set the max exfiltration rate to extremely large value
          grnd_E_max = HUGE(grnd_E_max)
     call land_surface_energy_balance( &
          grnd_T, grnd_liq, grnd_ice, grnd_latent, grnd_Tf, grnd_E_max, &
          fswg, &
          flwg0 + b0(iTv)*DflwgDTv, DflwgDTg + b1(iTv)*DflwgDTv, &
          Hg0   + b0(iTc)*DHgDTc,   DHgDTg   + b1(iTc)*DHgDTc, &
          Eg0   + b0(iqc)*DEgDqc,   DEgDTg   + b1(iqc)*DEgDqc, &
          G0,                       DGDTg, &
          ! output
          delta_Tg, Mg_imp )

! [X.5] calculate final value of other tendencies
     delta_qc = B0(iqc) + B1(iqc)*delta_Tg
     delta_Tc = B0(iTc) + B1(iTc)*delta_Tg
     delta_Tv = B0(iTv) + B1(iTv)*delta_Tg
     delta_wl = B0(iwl) + B1(iwl)*delta_Tg
     delta_ws = B0(iwf) + B1(iwf)*delta_Tg

! [X.6] calculate updated values of energy balance components used in further 
!       calculations
     flwg       = flwg0 + DflwgDTg*delta_Tg + DflwgDTv*delta_Tv
     evapg      = Eg0   + DEgDTg*delta_Tg   + DEgDqc*delta_qc
     sensg      = Hg0   + DHgDTg*delta_Tg   + DHgDTc*delta_Tc
     grnd_flux  = G0    + DGDTg*delta_Tg
     vegn_sens  = Hv0   + DHvDTv*delta_Tv   + DHvDTc*delta_Tc
     vegn_levap = Eli0  + DEliDTv*delta_Tv  + DEliDqc*delta_qc + DEliDwl*delta_wl + DEliDwf*delta_ws
     vegn_fevap = Esi0  + DEsiDTv*delta_Tv  + DEsiDqc*delta_qc + DEsiDwl*delta_wl + DEsiDwf*delta_ws
     vegn_uptk  = Et0   + DEtDTv*delta_Tv   + DEtDqc*delta_qc  + DEtDwl*delta_wl  + DEtDwf*delta_ws
     vegn_flw   = flwv0 + DflwvDTv*delta_Tv + DflwvDTg*delta_Tg
     vegn_hadvec= vegn_uptk*(clw*(soil_uptake_T-tfreeze)-cpw*(vegn_T-tfreeze))
! [X.7] calculate energy residuals due to cross-product of time tendencies
     tile%e_res_1 = canopy_air_mass*cpw*delta_qc*delta_Tc/delta_time
     tile%e_res_2 = delta_Tv*(clw*delta_Wl+csw*delta_Ws)/delta_time
! calculate the final value upward long-wave radiation flux from the land, to be 
! returned to the flux exchange.
     tile%lwup = ILa_dn - vegn_flw - flwg 
 
     if(is_watch_point())then
        write(*,*)'#### ground balance'
        __DEBUG2__(fswg,flwg)
        __DEBUG2__(sensg,evapg*grnd_latent)
        __DEBUG1__(grnd_flux)
        __DEBUG1__(Mg_imp)
        write(*,*)'#### implicit time steps'
        __DEBUG3__(delta_Tg, grnd_T,  grnd_T+delta_Tg )
        __DEBUG3__(delta_qc, cana_q,  cana_q+delta_qc )
        __DEBUG3__(delta_Tc, cana_T,  cana_T+delta_Tc )
        __DEBUG3__(delta_Tv, vegn_T,  vegn_T+delta_Tv )
        __DEBUG3__(delta_wl, vegn_Wl, vegn_Wl+delta_wl)
        __DEBUG3__(delta_ws, vegn_Ws, vegn_Ws+delta_ws)
        __DEBUG2__(tile%e_res_1, tile%e_res_2)
        write(*,*)'#### resulting fluxes'
        __DEBUG4__(flwg, evapg, sensg, grnd_flux)
        __DEBUG3__(vegn_levap,vegn_fevap,vegn_uptk)
        __DEBUG3__(vegn_sens,vegn_flw,vegn_hadvec)
        __DEBUG1__(Ea0+DEaDqc*delta_qc)
        __DEBUG2__(tile%cana%prog%q,cana_q)
     endif

     call cana_step_2 ( tile%cana, delta_Tc, delta_qc )

     if(associated(tile%vegn)) then
        call vegn_step_2 ( tile%vegn, tile%diag, &
             delta_Tv, delta_wl, delta_ws, &
             vegn_melt,  &
             vegn_ovfl_l,   vegn_ovfl_s, &
             vegn_ovfl_Hl, vegn_ovfl_Hs, &
             vegn_LMASS, vegn_FMASS, vegn_HEAT )
        ! calculate total amount of liquid and solid precipitation below the canopy
        vegn_lprec  = (1-vegn_ifrac)*precip_l + vegn_drip_l + vegn_ovfl_l
        vegn_fprec  = (1-vegn_ifrac)*precip_s + vegn_drip_s + vegn_ovfl_s
        ! calculate heat carried by liquid and solid precipitation below the canopy
        vegn_hlprec = clw*((1-vegn_ifrac)*precip_l*(precip_T-tfreeze) &
                         + vegn_drip_l*(vegn_T+delta_Tv-tfreeze)) &
                         + vegn_ovfl_Hl
        vegn_hfprec = csw*((1-vegn_ifrac)*precip_s*(precip_T-tfreeze) &
                         + vegn_drip_s*(vegn_T+delta_Tv-tfreeze)) &
                         + vegn_ovfl_Hs
        ! make sure the temperature of the snow falling below company is below freezing
        ! this correction was introduced in an attempt to fix the problem with ficticious 
        ! heat accumulating in near-zero-mass snow; however it does not seem to make a 
        ! difference.
        if(vegn_hfprec>0)then
           ! solid precipitation from vegetation carries positive energy -- we can't have
           ! that, because that would bring snow T above tfreeze, so convert excess to 
           ! liquid
           delta_fprec = min(vegn_fprec,vegn_hfprec/hlf)
           vegn_fprec = vegn_fprec - delta_fprec
           vegn_lprec = vegn_lprec + delta_fprec
           vegn_hfprec = 0
           ! we don't need to correct the vegn_hlprec since the temperature of additional
           ! liquid precip is tfreeze, and therefore its contribution to vegn_hlprec is
           ! exactly zero
        endif
        ! possibly we need to correct for the opposite situation: negative energy carried
        ! by liquid precipitation.
     else
        vegn_lprec  = precip_l
        vegn_fprec  = precip_s
        vegn_hlprec = precip_l*clw*(precip_T-tfreeze)
        vegn_hfprec = precip_s*csw*(precip_T-tfreeze)
        ! the fields below are only used in diagnostics
        vegn_melt   = 0
        vegn_hadvec = 0
        vegn_LMASS  = 0
        vegn_FMASS  = 0
        vegn_HEAT   = 0
        vegn_fsw    = 0
     endif

     call snow_step_2 ( tile%snow, &
          snow_subl,            &
          vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
          delta_Tg,  Mg_imp,  evapg,  fswg,  flwg,  sensg, &
          ! output:
          subs_DT, subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens, &
          snow_fsw, snow_flw, snow_sens, &
          snow_levap, snow_fevap, snow_hadvec, snow_melt, &
          snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
          snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, &
          snow_LMASS, snow_FMASS, snow_HEAT  )

     if(is_watch_point()) then
        write(*,*) 'subs_M_imp', subs_M_imp
     endif

     if (snow_active) then
        subs_G = snow_G_Z+snow_G_TZ*subs_DT
     else
        subs_G = 0
     endif
     
     if (associated(tile%glac)) then
        call glac_step_2 &
             ( tile%glac, tile%diag, subs_subl, snow_lprec, snow_hlprec, &
             subs_DT, subs_M_imp, subs_evap, &
             subs_levap, subs_fevap, &
             subs_hadvec, subs_melt, &
             subs_lrunf, subs_hlrunf, &
             subs_Ttop, subs_Ctop, &
             subs_LMASS, subs_FMASS, &
             subs_HEAT )
     else if (associated(tile%lake)) then
        call lake_step_2 &
             ( tile%lake, tile%diag, subs_subl, snow_lprec, snow_hlprec, &
             subs_DT, subs_M_imp, subs_evap, &
             subs_levap, subs_fevap, &
             subs_hadvec, subs_melt, &
             subs_lrunf, subs_hlrunf, &
             subs_Ttop, subs_Ctop, &
             subs_LMASS, subs_FMASS, &
             subs_HEAT )
     else if (associated(tile%soil)) then
        call soil_step_2 &
             ( tile%soil, tile%diag, subs_subl, snow_lprec, snow_hlprec, &
             vegn_uptk, &
             subs_DT, subs_M_imp, subs_evap, &
             ! output:
             subs_levap, subs_fevap, &
             subs_hadvec, subs_melt, &
             subs_lrunf, subs_hlrunf, &
             subs_Ttop, subs_Ctop, &
             subs_LMASS, subs_FMASS, &
             subs_HEAT )
     endif
     
! TEMP FIX: MAIN PROG SHOULD NOT TOUCH CONTENTS OF PROG VARS. ******
! ALSO, NEED TO ADJUST GFLUX FOR THIS TRANSFER OF ENERGY ******
! ALSO, DIAGNOSTICS IN COMPONENT MODULES SHOULD _FOLLOW_ THIS ADJUSTMENT******
     IF (LM2) THEN
        tile%snow%prog(1)%T = subs_Ttop
        tile%snow%prog(2)%T = subs_Ttop
        tile%snow%prog(3)%T = subs_Ttop
     ELSE
        if (sum(tile%snow%prog(:)%ws)>0) &
             tile%snow%prog(3)%T = (subs_Ctop*subs_Ttop +snow_Cbot*snow_Tbot) &
                           / (subs_Ctop+snow_Cbot) ! 3 is hardwired****
        if (sum(tile%snow%prog(:)%ws)>0)then
           if(associated(tile%glac)) tile%glac%prog(1)%T = tile%snow%prog(3)%T
           if(associated(tile%lake)) tile%lake%prog(1)%T = tile%snow%prog(3)%T
           if(associated(tile%soil)) tile%soil%prog(1)%T = tile%snow%prog(3)%T
        endif
     ENDIF

     if (associated(tile%vegn)) then
        ! do the calculations that require updated land surface prognostic variables
        call vegn_step_3 (tile%vegn, tile%soil, tile%cana%prog%T, precip_l+precip_s, tile%diag)
     endif
     
     call update_land_bc_fast (tile, i,j,k, cplr2land%p_surf(i,j,k), land2cplr)

     ! TEMPORARY FIX TO PREVENT NEGATIVE MASS LEAKING INTO RIVERS
!!$     if ((snow_frunf.ge.0).and.(subs_lrunf.ge.0).and.(snow_lrunf.ge.0)) then
        heat_pos(i,j) = heat_pos(i,j) + (snow_hfrunf + subs_hlrunf + snow_hlrunf)*tile%frac
        sol_pos (i,j) = sol_pos (i,j) + snow_frunf * tile%frac
        tot_pos (i,j) = tot_pos (i,j) + (snow_frunf + subs_lrunf + snow_lrunf)*tile%frac
!!$     endif

     ! ---- diagnostic section ----------------------------------------------
     call send_tile_data(id_frac, tile%frac, tile%diag)
     call send_tile_data(id_area, tile%frac*lnd%garea(i,j), tile%diag)
     call send_tile_data(id_ntiles, 1.0, tile%diag)     

     call send_tile_data(id_lprecv, cplr2land%lprec(i,j,k)-vegn_lprec, tile%diag)
     call send_tile_data(id_lprecs, vegn_lprec-snow_lprec, tile%diag)
     call send_tile_data(id_lprecg, snow_lprec, tile%diag)
     call send_tile_data(id_fprecv, cplr2land%fprec(i,j,k)-vegn_fprec, tile%diag)
     call send_tile_data(id_fprecs, vegn_fprec, tile%diag)
     call send_tile_data(id_levapv, vegn_levap, tile%diag)
     call send_tile_data(id_levaps, snow_levap, tile%diag)
     call send_tile_data(id_levapg, subs_levap, tile%diag)
     call send_tile_data(id_fevapv, vegn_fevap, tile%diag)
     call send_tile_data(id_fevaps, snow_fevap, tile%diag)
     call send_tile_data(id_fevapg, subs_fevap, tile%diag)

     call send_tile_data(id_lrunfs, snow_lrunf, tile%diag)
     call send_tile_data(id_frunfs, snow_frunf, tile%diag)
     call send_tile_data(id_lrunfg, subs_lrunf, tile%diag)

     call send_tile_data(id_transp, vegn_uptk, tile%diag)
     call send_tile_data(id_meltv, vegn_melt, tile%diag)
     call send_tile_data(id_melts, snow_melt, tile%diag)
     call send_tile_data(id_meltg, subs_melt, tile%diag)
     
     call send_tile_data(id_swdn_dir, ISa_dn_dir, tile%diag)
     call send_tile_data(id_swdn_dif, ISa_dn_dif, tile%diag)
     call send_tile_data(id_swup_dir, ISa_dn_dir*tile%land_refl_dir, tile%diag)
     call send_tile_data(id_swup_dif, ISa_dn_dif*tile%land_refl_dif, tile%diag)
     call send_tile_data(id_lwdn,     ILa_dn, tile%diag)

     call send_tile_data(id_fswv, vegn_fsw, tile%diag)
     call send_tile_data(id_fsws, snow_fsw, tile%diag)
     call send_tile_data(id_fswg, subs_fsw, tile%diag)
     call send_tile_data(id_flwv, vegn_flw, tile%diag)
     call send_tile_data(id_flws, snow_flw, tile%diag)
     call send_tile_data(id_flwg, subs_flw, tile%diag)

     call send_tile_data(id_sensv, vegn_sens, tile%diag)
     call send_tile_data(id_senss, snow_sens, tile%diag)
     call send_tile_data(id_sensg, subs_sens, tile%diag)

     call send_tile_data(id_gflux, subs_G, tile%diag)
     
     call send_tile_data(id_hadvecv, vegn_hadvec, tile%diag)
     call send_tile_data(id_hadvecs, snow_hadvec, tile%diag)
     call send_tile_data(id_hadvecg, subs_hadvec, tile%diag)
     
     call send_tile_data(id_LWSv, vegn_LMASS, tile%diag)
     call send_tile_data(id_LWSs, snow_LMASS, tile%diag)
     call send_tile_data(id_LWSg, subs_LMASS, tile%diag)
     call send_tile_data(id_FWSv, vegn_FMASS, tile%diag)
     call send_tile_data(id_FWSs, snow_FMASS, tile%diag)
     call send_tile_data(id_FWSg, subs_FMASS, tile%diag)
     call send_tile_data(id_HSv, vegn_HEAT, tile%diag)
     call send_tile_data(id_HSs, snow_HEAT, tile%diag)
     call send_tile_data(id_HSg, subs_HEAT, tile%diag)
     
     call send_tile_data(id_Trad, land2cplr%t_surf(i,j,k), tile%diag)
     call send_tile_data(id_Tca, land2cplr%t_ca(i,j,k), tile%diag)
#ifdef LAND_BND_TRACERS
     call send_tile_data(id_qca, land2cplr%tr(i,j,k,lnd%isphum), tile%diag)
#else
     call send_tile_data(id_qca, land2cplr%q_ca(i,j,k), tile%diag)
#endif
     call send_tile_data(id_z0m, land2cplr%rough_mom(i,j,k), tile%diag)
     call send_tile_data(id_z0s, land2cplr%rough_heat(i,j,k), tile%diag)
     
     call send_tile_data(id_lprec_l, cplr2land%lprec(i,j,k), tile%diag)
     call send_tile_data(id_fprec_l, cplr2land%fprec(i,j,k), tile%diag)
     call send_tile_data(id_levap, vegn_levap+snow_levap+subs_levap+vegn_uptk, tile%diag)
     call send_tile_data(id_fevap, vegn_fevap+snow_fevap+subs_fevap, tile%diag)
     
     call send_tile_data(id_lrunf, snow_lrunf+subs_lrunf, tile%diag)
     call send_tile_data(id_frunf, snow_frunf, tile%diag)
     call send_tile_data(id_melt, vegn_melt+snow_melt+subs_melt, tile%diag)
     call send_tile_data(id_fsw, vegn_fsw+snow_fsw+subs_fsw, tile%diag)
     call send_tile_data(id_flw, vegn_flw+snow_flw+subs_flw, tile%diag)
     call send_tile_data(id_sens, vegn_sens+snow_sens+subs_sens, tile%diag)
     call send_tile_data(id_hadvec, vegn_hadvec+snow_hadvec+subs_hadvec, tile%diag)

     call send_tile_data(id_LWS, vegn_LMASS+snow_LMASS+subs_LMASS, tile%diag)
     call send_tile_data(id_FWS, vegn_FMASS+snow_FMASS+subs_FMASS, tile%diag)
     call send_tile_data(id_HS, vegn_HEAT+snow_HEAT+subs_HEAT, tile%diag)
     call send_tile_data(id_precip, cplr2land%lprec(i,j,k)+cplr2land%fprec(i,j,k), tile%diag)
     call send_tile_data(id_wroff, snow_lrunf+subs_lrunf, tile%diag)
     call send_tile_data(id_sroff, snow_frunf, tile%diag)
     call send_tile_data(id_water, subs_LMASS+subs_FMASS, tile%diag)
     call send_tile_data(id_snow, snow_LMASS+snow_FMASS, tile%diag)
     call send_tile_data(id_evap, vegn_uptk+vegn_levap+snow_levap+subs_levap+vegn_fevap+snow_fevap+subs_fevap, tile%diag)

     call send_tile_data(id_con_g_h, con_g_h, tile%diag)

     call send_tile_data(id_e_res_1, tile%e_res_1, tile%diag)
     call send_tile_data(id_e_res_2, tile%e_res_2, tile%diag)

  enddo

  ! set values of tracer fluxes
  c_flux(:,:,1) = sol_pos  ! snow
  c_flux(:,:,2) = heat_pos ! heat
  do i_species = num_phys+1, num_species
    c_flux(:,:,i_species) = 0        ! age, species
    enddo

  ! update river state
  call update_river(tot_pos, c_flux, &
                    discharge2o_l, discharge2o_c, &
                    discharge2l_l, discharge2l_c)
! the following mixing and repartitioning makes the earlier partitioning look useless.
! but i'm allowing for future when river might make a partition without ocean-land leakage.
! at that time, the code below could be eliminated (but adjustment for ocean areas would
! still be needed).
  discharge_work = discharge2o_l + discharge2l_l
  discharge2o_l = 0.
  discharge2l_l = 0.
  where (gfrac(lnd%is:lnd%ie,lnd%js:lnd%je).ge.1.) 
      discharge2l_l = discharge_work
    elsewhere
      discharge2o_l = discharge_work
    endwhere
  do i_species = 1, num_species
    discharge_work = discharge2o_c(:,:,i_species) + discharge2l_c(:,:,i_species)
    discharge2o_c(:,:,i_species) = 0.
    discharge2l_c(:,:,i_species) = 0.
    where (gfrac(lnd%is:lnd%ie,lnd%js:lnd%je).ge.1.) 
        discharge2l_c(:,:,i_species) = discharge_work
      elsewhere
        discharge2o_c(:,:,i_species) = discharge_work
      endwhere
    enddo

  where (gfrac(lnd%is:lnd%ie,lnd%js:lnd%je).lt.1.) 
      land2cplr%discharge      = discharge2o_l        / (1-gfrac(lnd%is:lnd%ie,lnd%js:lnd%je))
      land2cplr%discharge_snow = discharge2o_c(:,:,1) / (1-gfrac(lnd%is:lnd%ie,lnd%js:lnd%je))
    endwhere

  ! advance land model time
  lnd%time = lnd%time + lnd%dt_fast

  ! send the accumulated diagnostics to the output
  call dump_tile_diag_fields(lnd%tile_map, lnd%time)

  if (id_dis_l2o > 0) used = send_data (id_dis_l2o, discharge2o_l, lnd%time) 
  if (id_dis_l2l > 0) used = send_data (id_dis_l2l, discharge2l_l, lnd%time) 
  if (id_dis_s2o > 0) used = send_data (id_dis_s2o, discharge2o_c(:,:,1), lnd%time) 
  if (id_dis_s2l > 0) used = send_data (id_dis_s2l, discharge2l_c(:,:,1), lnd%time) 

  call mpp_clock_end(landFastClock)
  call mpp_clock_end(landClock)
end subroutine update_land_model_fast


! ============================================================================
subroutine update_land_model_slow ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars
  integer :: i,j,k
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce, te

  call mpp_clock_begin(landClock)
  call mpp_clock_begin(landSlowClock)

  call update_vegn_slow( )
  ! send the accumulated diagnostics to the output
  ! call dump_tile_diag_fields(lnd%tile_map, lnd%time)

  ! get the current state of the land boundary for the coupler
  ce = first_elmt(lnd%tile_map,                  &
              is=lbound(cplr2land%t_flux,1), &
              js=lbound(cplr2land%t_flux,2)  )
  te = tail_elmt(lnd%tile_map)
  do while(ce /= te)
     ! calculate indices of the current tile in the input arrays;
     ! assume all the cplr2land components have the same lbounds
     call get_elmt_indices(ce,i,j,k)
     ! set this point coordinates as current for debug output
     call set_current_point(i,j,k)
     ! get pointer to current tile
     tile => current_tile(ce)
     ! advance enumerator to the next tile
     ce=next_elmt(ce)

     ! WARNING: the code below will not reproduce across restarts because
     ! values of DTc, Dqc and p_surf are dummy -- they must be either stored
     ! or excluded from the parameter list
     call update_land_bc_fast (tile, i,j,k, &
        cplr2land%p_surf(i,j,k), land2cplr, is_init=.true.)
  enddo

  call update_land_bc_slow( land2cplr )

  call mpp_clock_end(landClock)
  call mpp_clock_end(landSlowClock)
end subroutine update_land_model_slow


! ============================================================================
! solve for surface temperature. Do not employ any mass-limiting constraints
! on sublimation or evaporation (which can be dealt with by allowing negative
! stores, if necessary), but ensure that melt does not exceed available
! snow or soil ice (which would create energy-balance problems). also ensure
! that melt and temperature are consistent and that evaporation from liquid
! soil water does not exceed exfiltration rate limit, if one is supplied.
subroutine land_surface_energy_balance ( &
     ! surface parameters
     grnd_T,          & ! ground temperature
     grnd_liq, grnd_ice, & ! amount of water available for freeze or melt on the surface, kg/m2
     grnd_latent,     & ! specific heat of vaporization for the ground
     grnd_Tf,         & ! ground freezing temperature
     grnd_E_max,      & ! exfiltration rate limit, kg/(m2 s)
     ! components of the ground energy balance linearization. Note that those
     ! are full derivatives, which include the response of the other surface
     ! cheme parameters to the change in ground temperature.
     fswg,            & ! net short-wave
     flwg0, DflwgDTg, & ! net long-wave
     Hg0,   DHgDTg,   & ! sensible heat
     Eg0,   DEgDTg,   & ! latent heat
     G0,    DGDTg,    & ! sub-surface heat 
     delta_Tg,        & ! surface temperature change for the time step
     Mg_imp          )  ! implicit melt, kg/m2

  real, intent(in) :: & 
     grnd_T,          & ! ground temperature
     grnd_liq, grnd_ice, & ! amount of water available for freeze or melt on the surface
     grnd_latent,     & ! specific heat of vaporization for the ground
     grnd_Tf,         & ! ground freezing temperature
     grnd_E_max,      & ! exfiltration rate limit, kg/(m2 s)
     fswg,            & ! net short-wave
     flwg0, DflwgDTg, & ! net long-wave
     Hg0,   DHgDTg,   & ! sensible heat
     Eg0,   DEgDTg,   & ! latent heat
     G0,    DGDTg       ! sub-surface heat 
  real, intent(out) :: &
     delta_Tg,        & ! chnage in surface temperature
     Mg_imp             ! mass of surface ice melted (or water frozen) during the 
                        ! time step, kg/m2

  real :: grnd_B     ! surface energy balance
  real :: grnd_DBDTg ! full derivative of grnd_B w.r.t. surace temperature

  ! IMPLEMENTATION NOTE: instead of passing "snow_active" and using global variable
  ! "use_E_max", this code simply compares the evaporation value with the maximum
  ! soil effiltration rate passed by the caling subroutine. The trick is that the
  ! input effiltration rate is huge in case when soil_E_max should not be used, so 
  ! that evaroration is guaranteed to be smaller. I feel that this results in cleaner
  ! interface and implementation.

  grnd_B      = fswg + flwg0    - Hg0     - grnd_latent*Eg0    - G0
  grnd_DBDTg  =        DflwgDTg - DHgDTg  - grnd_latent*DEgDTg - DGDTg

  ! make an intial estimate of the ground temperature change
  delta_Tg = - grnd_B/grnd_DBDTg
  ! calculate phase change on the ground, if necessary
  if     (grnd_ice>0.and.grnd_T+delta_Tg>grnd_Tf) then ! melt > 0
     Mg_imp =  min(grnd_ice,  grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
  elseif (grnd_liq>0.and.grnd_T+delta_Tg<grnd_Tf) then ! melt < 0
     Mg_imp = -min(grnd_liq, -grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
  else
     Mg_imp = 0
  endif
  ! adjust temperature change for the ground melt
  delta_Tg = -(grnd_B - Mg_imp*hlf/delta_time)/grnd_DBDTg

  if(is_watch_point())then
     write(*,*)'#### ground balance after melt'
     __DEBUG2__(grnd_B, grnd_DBDTg)
     __DEBUG2__(Mg_imp, delta_Tg)
  endif

  ! Solution above is OK if soil_rh was known. However, if we only know
  ! the effiltration-rate limit from the soil, then we next need to
  ! solve the equation assuming that limit holds. In that case,
  ! compute alternate solution, including melting checks, wherever the
  ! first solution yielded Eg in excess of rate limit.
  if ( Eg0 + DEgDTg*delta_Tg > grnd_E_max ) then
     ! recalculate the ground heat balance with evap fixed to grnd_E_max
     ! and recalculte the balance derivative w.r.t. Tg with fixed ground evap
     grnd_B     = fswg + flwg0    - Hg0  - grnd_latent*grnd_E_max - G0
     grnd_DBDTg =        DflwgDTg - DHgDTg                        - DGDTg
     ! recalculate the ground temperature tendency
     delta_Tg = - grnd_B/grnd_DBDTg
     ! recalculate the ground metl again
     if      ( grnd_ice > 0 .and. grnd_T+delta_Tg > grnd_Tf ) then ! melt>0
        Mg_imp =  min(grnd_ice,  grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
     else if ( grnd_liq > 0 .and. grnd_T+delta_Tg < grnd_Tf ) then ! melt<0
        Mg_imp = -min(grnd_liq, -grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
     else
        Mg_imp = 0
     endif
     ! adjust temperature change for ground melt
     delta_Tg = -(grnd_B - Mg_imp*hlf/delta_time)/grnd_DBDTg
     if(is_watch_point())then
        write(*,*)'#### ground balance after melt inside use_E_max'
        __DEBUG2__(grnd_B,grnd_DBDTg)
        __DEBUG2__(Mg_imp,delta_Tg)
     endif
  endif
  
end subroutine land_surface_energy_balance


! ============================================================================
subroutine update_land_bc_fast (tile, i,j,k, p_surf, land2cplr, is_init)
  type(land_tile_type), intent(inout) :: tile
  integer             , intent(in) :: i,j,k
  real                , intent(in) :: p_surf
  type(land_data_type), intent(inout) :: land2cplr
  logical, optional :: is_init

  ! ---- local vars
  real :: &
         grnd_T, subs_z0m, subs_z0s, &
                 snow_z0s, snow_z0m, &
         snow_area, snow_depth

  real :: subs_refl_dir(NBANDS), subs_refl_dif(NBANDS) ! direct and diffuse albedos
  real :: subs_refl_lw ! reflectance for thermal radiation
  real :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS) ! direct and diffuse albedos of snow
  real :: snow_refl_lw ! snow reflectance for thermal radiation
  real :: snow_emis ! snow emissivity
  real :: snow_area_rad ! "snow area for radiation calculations" -- introduced
                        ! to reproduce lm2 behavior
  real :: vegn_refl_lw, vegn_tran_lw ! reflectance and transmittance of vegetation for thermal radiation
  real :: vegn_refl_dif(NBANDS), vegn_tran_dif(NBANDS) ! reflectance and transmittance of vegetation for diffuse light
  real :: vegn_refl_dir(NBANDS), vegn_tran_dir(NBANDS) ! reflectance and transmittance of vegetation for direct light
  real :: vegn_tran_dir_dir(NBANDS) ! (?)
  real :: &
         vegn_Tv,     &
         vegn_cover,  &
         vegn_height, vegn_lai, vegn_sai, vegn_d_leaf
  logical :: do_update

  real :: cosz    ! cosine of solar zenith angle
  real :: fracday ! daytime fraction of time interval
  real :: rrsun   ! earth-sun distance (r) relative to semi-major axis
                  ! of orbital ellipse (a) : (a/r)**2
  vegn_Tv = 0

  do_update = .not.present(is_init)

  ! on initialization the albedos are calculated for the current time step ( that is, interval
  ! lnd%time, lnd%time+lnd%dt_fast); in the course of the run this subroutine is called
  ! at the end of time step (but before time is advanced) to calculate the radiative properties 
  ! for the _next_ time step
  if (do_update) then
     call diurnal_solar(lnd%glat(i,j), lnd%glon(i,j), lnd%time+lnd%dt_fast, &
          cosz, fracday, rrsun, lnd%dt_fast)
  else
     call diurnal_solar(lnd%glat(i,j), lnd%glon(i,j), lnd%time, &
          cosz, fracday, rrsun, lnd%dt_fast)
  endif
  
  if (associated(tile%glac)) then
     call glac_radiation(tile%glac, subs_refl_dir, subs_refl_dif, subs_refl_lw, tile%grnd_emis)
     call glac_diffusion(tile%glac, subs_z0s, subs_z0m )
  else if (associated(tile%lake)) then
     call lake_radiation(tile%lake, subs_refl_dir, subs_refl_dif, subs_refl_lw, tile%grnd_emis)
     call lake_diffusion(tile%lake, subs_z0s, subs_z0m )
  else if (associated(tile%soil)) then
     call soil_radiation(tile%soil, subs_refl_dir, subs_refl_dif, subs_refl_lw, tile%grnd_emis)
     call soil_diffusion(tile%soil, subs_z0s, subs_z0m )
  else
     call error_mesg('update_land_bc_fast','none of the surface tiles exists',FATAL)
  endif

  call snow_radiation ( tile%snow, snow_refl_dif, snow_refl_dir, snow_refl_lw, snow_emis)
  call snow_get_depth_area ( tile%snow, snow_depth, snow_area )
  call snow_diffusion ( tile%snow, snow_z0s, snow_z0m )

  if (associated(tile%vegn)) then
     call update_derived_vegn_data(tile%vegn)
     ! USE OF SNOWPACK RAD PROPERTIES FOR INTERCEPTED SNOW IS ERRONEOUS,
     ! NEEDS TO BE CHANGED. TEMPORARY.
     call vegn_radiation ( tile%vegn, cosz, snow_depth, snow_refl_dif, snow_emis, &
                   vegn_refl_dif, vegn_tran_dif, &
                   vegn_refl_dir, vegn_tran_dir, vegn_tran_dir_dir, &
                   vegn_refl_lw, vegn_tran_lw)
     ! (later see if we can remove vegn_cover from c-a-radiation...) TEMPORARY
     call vegn_get_cover ( tile%vegn, snow_depth, vegn_cover)
     call vegn_diffusion ( tile%vegn, vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf)
  else
     ! set radiative properties for null vegetation
     vegn_refl_dif     = 0
     vegn_tran_dif     = 1
     vegn_refl_dir     = 0
     vegn_tran_dir     = 0
     vegn_tran_dir_dir = 1
     vegn_refl_lw      = 0 
     vegn_tran_lw      = 1
     ! set cover for null vegetation
     vegn_cover        = 0
     ! set other parameters for null vegetation
     vegn_height       = 0
     vegn_lai          = 0
     vegn_sai          = 0
     vegn_d_leaf       = 0
  endif

  ! store the values of long-wave optical properties to be used in the update_land_model_fast
  tile%surf_refl_lw = subs_refl_lw  + (snow_refl_lw  - subs_refl_lw ) * snow_area
  tile%vegn_refl_lw = vegn_refl_lw
  tile%vegn_tran_lw = vegn_tran_lw

  
  if(is_watch_point()) then
     write(*,*) '#### update_land_bc_fast ### checkpoint 1 ####'
     __DEBUG3__(cosz, fracday, rrsun)
     __DEBUG2__(vegn_lai,vegn_sai)
     __DEBUG1__(subs_refl_dif)
     __DEBUG1__(subs_refl_dir)
     __DEBUG1__(vegn_refl_dif)
     __DEBUG1__(vegn_tran_dif)
     __DEBUG1__(vegn_refl_dir)
     __DEBUG1__(vegn_tran_dir)
     __DEBUG1__(vegn_tran_dir_dir)
     __DEBUG2__(vegn_refl_lw,vegn_tran_lw)
     write(*,*) '#### update_land_bc_fast ### end of checkpoint 1 ####'
  endif

  snow_area_rad = snow_area
  if (lm2) then
     if(associated(tile%glac)                    ) snow_area_rad = 0
     if(associated(tile%soil).and.vegn_cover>0.01) snow_area_rad = 1
  endif

  call cana_radiation( lm2, &
       subs_refl_dir, subs_refl_dif, subs_refl_lw, &
       snow_refl_dir, snow_refl_dif, snow_refl_lw, &
       snow_area_rad,  &
       vegn_refl_dir, vegn_refl_dif, vegn_tran_dir, vegn_tran_dif, &
       vegn_tran_dir_dir, vegn_refl_lw, vegn_tran_lw, &
       vegn_cover, &
       tile%Sg_dir, tile%Sg_dif, tile%Sv_dir, tile%Sv_dif, &
       tile%land_refl_dir, tile%land_refl_dif )

  call cana_roughness( lm2, &
     subs_z0m, subs_z0s, &
     snow_z0m, snow_z0s, snow_area, &
     vegn_cover,  vegn_height, vegn_lai, vegn_sai, &
     tile%land_d, tile%land_z0m, tile%land_z0s)

  if(is_watch_point()) then
     write(*,*) '#### update_land_bc_fast ### checkpoint 2 ####'
     write(*,*) 'Sg_dir', tile%Sg_dir
     write(*,*) 'Sg_dif', tile%Sg_dif
     write(*,*) 'Sv_dir', tile%Sv_dir
     write(*,*) 'Sv_dif', tile%Sv_dif
     write(*,*) 'land_albedo_dir', tile%land_refl_dir
     write(*,*) 'land_albedo_dif', tile%land_refl_dif
     write(*,*) 'land_z0m', tile%land_z0m
     write(*,*) '#### update_land_bc_fast ### end of checkpoint 2 ####'
  endif

  land2cplr%t_surf         (i,j,k) = tfreeze
  land2cplr%t_ca           (i,j,k) = tfreeze
#ifdef LAND_BND_TRACERS
  land2cplr%tr             (i,j,k, lnd%isphum) = 0.0
#else
  land2cplr%q_ca           (i,j,k) = 0.0
#endif
  land2cplr%albedo         (i,j,k) = 0.0
  land2cplr%albedo_vis_dir (i,j,k) = 0.0
  land2cplr%albedo_nir_dir (i,j,k) = 0.0
  land2cplr%albedo_vis_dif (i,j,k) = 0.0
  land2cplr%albedo_nir_dif (i,j,k) = 0.0
  land2cplr%rough_mom      (i,j,k) = 0.1
  land2cplr%rough_heat     (i,j,k) = 0.1
  land2cplr%rough_scale    (i,j,k) = 0.1

  ! Calculate radiative surface temperature. lwup can't be calculated here
  ! based on the available temperatures because it's a result of the implicit 
  ! time step: lwup = lwup0 + DlwupDTg*delta_Tg + ..., so we have to carry it
  ! from the update_land_fast 
  ! Consequence: since update_landbc_fast is called once at the very beginning of 
  ! every run (before update_land_fast is called) lwup from the previous step 
  ! must be stored in the in the restart for reproducibility
!!$  if(do_update) then
     land2cplr.t_surf(i,j,k) = ( tile.lwup/stefan ) ** 0.25
!!$  else
!!$     land2cplr.t_surf(i,j,k) = tile%cana%prog%T
!!$  endif

  if (associated(tile%glac)) call glac_get_sfc_temp(tile%glac, grnd_T)
  if (associated(tile%lake)) call lake_get_sfc_temp(tile%lake, grnd_T)
  if (associated(tile%soil)) call soil_get_sfc_temp(tile%soil, grnd_T)
  if (snow_area > 0)         call snow_get_sfc_temp(tile%snow, grnd_T)

! the piece of code below shouldn't be uncommented, since it would lead to non-conservation
! of water mass.
!!$! *** following was added to get q tracked over glacier and lake....TEMP
!!$! this could help where we want t_ca to jump due to switch of snow presence
!!$  if(do_update) then
!!$     if (associated(tile%glac).or.associated(tile%lake)) then
!!$        tile%cana%prog%T = grnd_T
!!$        call qscomp(grnd_T, p_surf, tile%cana%prog%q)
!!$     endif
!!$  endif

  ! set the boundary conditions for the flux exchange
  land2cplr%mask           (i,j,k) = .TRUE.
  land2cplr%tile_size      (i,j,k) = tile%frac

  land2cplr%t_ca           (i,j,k) = tile%cana%prog%T
#ifdef LAND_BND_TRACERS
  land2cplr%tr             (i,j,k,lnd%isphum) = tile%cana%prog%q
#else
  land2cplr%q_ca           (i,j,k) = tile%cana%prog%q
#endif
  land2cplr%albedo_vis_dir (i,j,k) = tile%land_refl_dir(BAND_VIS)
  land2cplr%albedo_nir_dir (i,j,k) = tile%land_refl_dir(BAND_NIR)
  land2cplr%albedo_vis_dif (i,j,k) = tile%land_refl_dif(BAND_VIS)
  land2cplr%albedo_nir_dif (i,j,k) = tile%land_refl_dif(BAND_NIR)
  land2cplr%albedo         (i,j,k) = SUM(tile%land_refl_dir + tile%land_refl_dif)/4 ! incorrect, replace with proper weighting later
  land2cplr%rough_mom      (i,j,k) = tile%land_z0m           !  TEMPORARY ***
  land2cplr%rough_heat     (i,j,k) = tile%land_z0s           !  TEMPORARY ***
  land2cplr%rough_scale    (i,j,k) = land2cplr%rough_mom (i,j,k)    !  TEMPORARY ***  

  if(is_watch_point()) then
     write(*,*)'#### update_land_bc_fast ### output ####'
     write(*,*)'land2cplr%mask',land2cplr%mask(i,j,k)
     write(*,*)'land2cplr%tile_size',land2cplr%tile_size(i,j,k)
     write(*,*)'land2cplr%t_surf',land2cplr%t_surf(i,j,k)
     write(*,*)'land2cplr%t_ca',land2cplr%t_ca(i,j,k)
     write(*,*)'land2cplr%albedo',land2cplr%albedo(i,j,k)
     write(*,*)'land2cplr%rough_mom',land2cplr%rough_mom(i,j,k)
     write(*,*)'land2cplr%rough_heat',land2cplr%rough_heat(i,j,k)
#ifdef LAND_BND_TRACERS
     write(*,*)'land2cplr%tr',land2cplr%tr(i,j,k,:)
#else
     write(*,*)'land2cplr%q_ca',land2cplr%q_ca(i,j,k)
#endif
     write(*,*)'#### update_land_bc_fast ### end of output ####'
  endif

  ! ---- diagnostic section
  call send_tile_data(id_cosz, cosz, tile%diag)

  call send_tile_data(id_albedo_dir, tile%land_refl_dir, tile%diag)
  call send_tile_data(id_albedo_dif, tile%land_refl_dif, tile%diag)

  call send_tile_data(id_subs_refl_dif, subs_refl_dif, tile%diag)
  call send_tile_data(id_subs_refl_dir, subs_refl_dir, tile%diag)

  call send_tile_data(id_vegn_refl_dif, vegn_refl_dif, tile%diag)
  call send_tile_data(id_vegn_tran_dif, vegn_tran_dif, tile%diag)
  call send_tile_data(id_vegn_refl_dir, vegn_refl_dir,     tile%diag)
  call send_tile_data(id_vegn_tran_dir, vegn_tran_dir_dir, tile%diag)
  call send_tile_data(id_vegn_sctr_dir, vegn_tran_dir,     tile%diag)

  call send_tile_data(id_vegn_refl_lw, vegn_refl_lw, tile%diag)
  call send_tile_data(id_vegn_tran_lw, vegn_tran_lw, tile%diag)

  call send_tile_data(id_vegn_cover, vegn_cover, tile%diag)
  call send_tile_data(id_grnd_T,     grnd_T,     tile%diag)

  ! --- debug section
  if(120.0<land2cplr%t_ca(i,j,k).and.land2cplr%t_ca(i,j,k)<373.0) then
     continue
  else
     write(*,'(a,3i4,g)')'update_land_bc_fast: t_ca out of range',i,j,k,land2cplr%t_ca(i,j,k)
  endif
end subroutine update_land_bc_fast


! ============================================================================
subroutine update_land_bc_slow (land2cplr)
  type(land_data_type), intent(inout) :: land2cplr

  call update_topo_rough(land2cplr%rough_scale)
  where (land2cplr%mask) &
       land2cplr%rough_scale = max(land2cplr%rough_mom,land2cplr%rough_scale)
end subroutine update_land_bc_slow


! ============================================================================
subroutine Lnd_stock_pe(bnd,index,value)

type(land_data_type), intent(in)  :: bnd   ! Not needed by this routine, but needed to satisfy standard interface
integer             , intent(in)  :: index
real                , intent(out) :: value


value = 0.0

end subroutine Lnd_stock_pe


! ============================================================================
! initialize horizontal axes for land grid so that all sub-modules can use them,
! instead of creating their own
subroutine land_diag_init(glonb, glatb, glon, glat, time, domain, &
     id_lon, id_lat, id_band)
  real, intent(in) :: &
       glonb(:,:), glatb(:,:), & ! longitude and latitude boundaries of grid cells,
                                 ! specified for the global grid (not local domain)
       glon(:,:), glat(:,:)      ! lon and lat of respective grid cell centers
  type(time_type), intent(in) :: time ! initial time for diagnostic fields
  type(domain2d), intent(in)  :: domain
  integer, intent(out) :: &
       id_lon, id_lat, id_band   ! IDs of respective diag. manager axes

  ! ---- local vars ----------------------------------------------------------
  integer :: id_lonb, id_latb ! IDs for cell boundaries
  integer :: nlon, nlat       ! sizes of respective axes
  real    :: rad2deg          ! conversion factor radian -> degrees
  integer :: axes(2)          ! array of axes for 2-D fields
  integer :: i

  rad2deg = 180./pi
  nlon = size(glonb,1)-1
  nlat = size(glatb,2)-1

  if(mpp_get_ntile_count(lnd%domain)==1) then
     ! grid has just one tile, so we assume that the grid is regular lat-lon
     ! define longitude axes and its edges
     id_lonb = diag_axis_init ( &
          'lonb', glonb(:,1)*rad2deg, 'degrees_E', 'X', 'longitude edges', &
          set_name='land', domain2=domain )
     id_lon  = diag_axis_init (                                                &
          'lon',  glon(:,1)*rad2deg, 'degrees_E', 'X',  &
          'longitude', set_name='land',  edges=id_lonb, domain2=domain )

     ! define latitude axes and its edges
     id_latb = diag_axis_init ( &
          'latb', glatb(1,:)*rad2deg, 'degrees_N', 'Y', 'latitude edges',  &
          set_name='land',  domain2=domain   )
     id_lat = diag_axis_init (                                                &
          'lat',  glat(1,:)*rad2deg, 'degrees_N', 'Y', &
          'latitude', set_name='land', edges=id_latb, domain2=domain   )
  else
     id_lon = diag_axis_init ( 'grid_xt', (/(real(i),i=1,nlon)/), 'degrees_E', 'X', &
            'T-cell longitude', set_name='land',  domain2=domain )
     id_lat = diag_axis_init ( 'grid_yt', (/(real(i),i=1,nlat)/), 'degrees_N', 'Y', &
             'T-cell latitude', set_name='land',  domain2=domain )
  endif
  id_band = diag_axis_init (                                                &
       'band',  (/1.0,2.0/), 'unitless', 'Z', &
       'spectral band', set_name='land' )

  ! set up an array of axes, for convenience
  axes = (/id_lon, id_lat/)

  ! register static diagnostic fields

  id_landarea = register_static_field ( module_name, 'land_area', axes, &
       'land area in grid cell', 'm2', missing_value=-1.0 )
  id_landfrac = register_static_field ( module_name, 'land_frac', axes, &
       'fraction of land in grid cell','unitless', missing_value=-1.0 ) 

  ! register regular (dynamic) diagnostic fields

  id_ntiles = register_tiled_diag_field(module_name,'ntiles',axes,  &
       time, 'number of tiles', 'unitless', missing_value=-1.0, op=OP_SUM)
  id_frac = register_tiled_diag_field(module_name,'frac', axes,&
       time, 'fraction of land area', 'unitless', missing_value=-1.0, op=OP_SUM )
  id_area = register_tiled_diag_field(module_name,'area', axes,&
       time, 'area in the grid cell', 'm2', missing_value=-1.0, op=OP_SUM )

  id_dis_l2o   = register_diag_field ( module_name, 'disl2o', axes, &
       time, 'disl2o', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_dis_l2l   = register_diag_field ( module_name, 'disl2l', axes, &
       time, 'disl2l', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_dis_s2o   = register_diag_field ( module_name, 'diss2o', axes, &
       time, 'diss2o', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_dis_s2l   = register_diag_field ( module_name, 'diss2l', axes, &
       time, 'diss2l', 'kg/(m2 s)', missing_value=-1.0e+20 )

  id_lprecv = register_tiled_diag_field ( module_name, 'lprecv', axes, time, &
             'net rainfall to vegetation', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_lprecs = register_tiled_diag_field ( module_name, 'lprecs', axes, time, &
             'rainfall to snow, minus drainage', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_lprecg = register_tiled_diag_field ( module_name, 'lprecg', axes, time, &
             'effective rainfall to ground sfc', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fprecv = register_tiled_diag_field ( module_name, 'fprecv', axes, time, &
             'net snowfall to vegetation', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fprecs = register_tiled_diag_field ( module_name, 'fprecs', axes, time, &
             'effective snowfall to snowpack', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_levapv = register_tiled_diag_field ( module_name, 'levapv', axes, time, &
             'vapor flux leaving intercepted liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_levaps = register_tiled_diag_field ( module_name, 'levaps', axes, time, &
             'vapor flux leaving snow liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_levapg = register_tiled_diag_field ( module_name, 'levapg', axes, time, &
             'vapor flux from ground liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fevapv = register_tiled_diag_field ( module_name, 'fevapv', axes, time, &
             'vapor flux leaving vegn ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fevaps = register_tiled_diag_field ( module_name, 'fevaps', axes, time, &
             'vapor flux leaving snow ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fevapg = register_tiled_diag_field ( module_name, 'fevapg', axes, time, &
             'vapor flux leaving ground ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_lrunfs  = register_tiled_diag_field ( module_name, 'lrunfs', axes, time, &
             'rate of liq runoff via calving', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lrunfg  = register_tiled_diag_field ( module_name, 'lrunfg', axes, time, &
             'rate of liq runoff, ground', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_frunfs  = register_tiled_diag_field ( module_name, 'frunfs', axes, time, &
             'rate of solid calving', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_transp  = register_tiled_diag_field ( module_name, 'transp', axes, time, &
             'transpiration; = uptake by roots', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_meltv   = register_tiled_diag_field ( module_name, 'meltv', axes, time, &
             'rate of melt, interception', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_melts   = register_tiled_diag_field ( module_name, 'melts', axes, time, &
             'rate of snow melt', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_meltg   = register_tiled_diag_field ( module_name, 'meltg', axes, time, &
             'rate of ground thaw', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fswv    = register_tiled_diag_field ( module_name, 'fswv', axes, time, &
             'net sw rad to vegetation', 'W/m2', missing_value=-1.0e+20 )
  id_fsws    = register_tiled_diag_field ( module_name, 'fsws', axes, time, &
             'net sw rad to snow', 'W/m2', missing_value=-1.0e+20 )
  id_fswg    = register_tiled_diag_field ( module_name, 'fswg', axes, time, &
             'net sw rad to ground', 'W/m2', missing_value=-1.0e+20 )
  id_flwv    = register_tiled_diag_field ( module_name, 'flwv', axes, time, &
             'net lw rad to vegetation', 'W/m2', missing_value=-1.0e+20 )
  id_flws    = register_tiled_diag_field ( module_name, 'flws', axes, time, &
             'net lw rad to snow', 'W/m2', missing_value=-1.0e+20 )
  id_flwg    = register_tiled_diag_field ( module_name, 'flwg', axes, time, &
             'net lw rad to ground', 'W/m2', missing_value=-1.0e+20 )
  id_sensv   = register_tiled_diag_field ( module_name, 'sensv', axes, time, &
             'sens heat flux from vegn', 'W/m2', missing_value=-1.0e+20 )
  id_senss   = register_tiled_diag_field ( module_name, 'senss', axes, time, &
             'sens heat flux from snow', 'W/m2', missing_value=-1.0e+20 )
  id_sensg   = register_tiled_diag_field ( module_name, 'sensg', axes, time, &
             'sens heat flux from ground', 'W/m2', missing_value=-1.0e+20 )
  id_gflux   = register_tiled_diag_field ( module_name, 'gflux', axes, time, &
             'sens heat into ground from snow', 'W/m2', missing_value=-1.0e+20 )
  id_hadvecv = register_tiled_diag_field ( module_name, 'hadvecv', axes, time, &
             'net advect. of heat to veg', 'W/m2', missing_value=-1.0e+20 )
  id_hadvecs = register_tiled_diag_field ( module_name, 'hadvecs', axes, time, &
             'net advect. of heat to snow', 'W/m2', missing_value=-1.0e+20 )
  id_hadvecg = register_tiled_diag_field ( module_name, 'hadvecg', axes, time, &
             'net advect. of heat to grnd', 'W/m2', missing_value=-1.0e+20 )
  id_LWSv    = register_tiled_diag_field ( module_name, 'LWSv', axes, time, &
             'liquid interception storage', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSs    = register_tiled_diag_field ( module_name, 'LWSs', axes, time, &
             'liquid storage in snowpack', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSg    = register_tiled_diag_field ( module_name, 'LWSg', axes, time, &
             'liquid ground storage', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSv    = register_tiled_diag_field ( module_name, 'FWSv', axes, time, &
             'frozen interception storage', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSs    = register_tiled_diag_field ( module_name, 'FWSs', axes, time, &
             'frozen storage in snowpack', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSg    = register_tiled_diag_field ( module_name, 'FWSg', axes, time, &
             'frozen ground storage', 'kg/m2', missing_value=-1.0e+20 )
  id_HSv     = register_tiled_diag_field ( module_name, 'HSv', axes, time, &
             'interception heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSs     = register_tiled_diag_field ( module_name, 'HSs', axes, time, &
             'heat storage in snowpack', 'J/m2', missing_value=-1.0e+20 )
  id_HSg     = register_tiled_diag_field ( module_name, 'HSg', axes, time, &
             'ground heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_Trad    = register_tiled_diag_field ( module_name, 'Trad', axes, time, &
             'radiative sfc temperature', 'degK', missing_value=-1.0e+20 )
  id_Tca     = register_tiled_diag_field ( module_name, 'Tca', axes, time, &
             'canopy-air temperature', 'degK', missing_value=-1.0e+20 )
  id_qca     = register_tiled_diag_field ( module_name, 'qca', axes, time, &
             'canopy-air specific humidity', 'kg/kg', missing_value=-1.0e+20 )
  id_z0m     = register_tiled_diag_field ( module_name, 'z0m', axes, time, &
             'momentum roughness of land', 'm', missing_value=-1.0e+20 )
  id_z0s     = register_tiled_diag_field ( module_name, 'z0s', axes, time, &
             'scalar roughness of land', 'm', missing_value=-1.0e+20 )

  id_lprec_l = register_tiled_diag_field ( module_name, 'lprec_l', axes, time, &
             'rainfall to land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fprec_l = register_tiled_diag_field ( module_name, 'fprec_l', axes, time, &
             'snowfall to land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_levap   = register_tiled_diag_field ( module_name, 'levap', axes, time, &
             'vapor flux from all liq (inc Tr)', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fevap   = register_tiled_diag_field ( module_name, 'fevap', axes, time, &
             'vapor flux from all ice', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lrunf   = register_tiled_diag_field ( module_name, 'lrunf', axes, time, &
             'total rate of liq runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_frunf   = register_tiled_diag_field ( module_name, 'frunf', axes, time, &
             'total rate of solid runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_melt    = register_tiled_diag_field ( module_name, 'melt', axes, time, &
             'total rate of melt', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fsw     = register_tiled_diag_field ( module_name, 'fsw', axes, time, &
             'net sw rad to land', 'W/m2', missing_value=-1.0e+20 )
  id_flw     = register_tiled_diag_field ( module_name, 'flw', axes, time, &
             'net lw rad to land', 'W/m2', missing_value=-1.0e+20 )
  id_sens    = register_tiled_diag_field ( module_name, 'sens', axes, time, &
             'sens heat flux from land', 'W/m2', missing_value=-1.0e+20 )
  id_hadvec = register_tiled_diag_field ( module_name, 'hadvec', axes, time, &
             'net advect. of heat to land', 'W/m2', missing_value=-1.0e+20 )
  id_LWS = register_tiled_diag_field ( module_name, 'LWS', axes, time, &
             'liquid storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_FWS = register_tiled_diag_field ( module_name, 'FWS', axes, time, &
             'frozen storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_HS = register_tiled_diag_field ( module_name, 'HS', axes, time, &
             'land heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_precip = register_tiled_diag_field ( module_name, 'precip', axes, time, &
             'precipitation rate', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_wroff = register_tiled_diag_field ( module_name, 'wroff', axes, time, &
             'rate of liquid runoff to rivers', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_sroff = register_tiled_diag_field ( module_name, 'sroff', axes, time, &
             'rate of solid runoff to rivers', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_water = register_tiled_diag_field ( module_name, 'water', axes, time, &
             'column-integrated soil water', 'kg/m2', missing_value=-1.0e+20 )
  id_snow = register_tiled_diag_field ( module_name, 'snow', axes, time, &
             'column-integrated snow water', 'kg/m2', missing_value=-1.0e+20 )
!  id_groundwater = register_tiled_diag_field ( module_name, 'groundwater', axes, time, 'column-integrated gw', 'kg/m2', missing_value=-1.0e+20 )
  id_evap = register_tiled_diag_field ( module_name, 'evap', axes, time, &
             'vapor flux up from land', 'kg/(m2 s)', missing_value=-1.0e+20 )

  id_swdn_dir = register_tiled_diag_field ( module_name, 'swdn_dir', (/id_lon,id_lat,id_band/), time, &
       'downward direct short-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
  id_swdn_dif = register_tiled_diag_field ( module_name, 'swdn_dif', (/id_lon,id_lat,id_band/), time, &
       'downward diffuse short-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
  id_swup_dir = register_tiled_diag_field ( module_name, 'swup_dir', (/id_lon,id_lat,id_band/), time, &
       'direct short-wave radiation flux reflected by the land surface', 'W/m2', missing_value=-999.0)
  id_swup_dif = register_tiled_diag_field ( module_name, 'swup_dif', (/id_lon,id_lat,id_band/), time, &
       'diffuse short-wave radiation flux reflected by the land surface', 'W/m2', missing_value=-999.0)
  id_lwdn = register_tiled_diag_field ( module_name, 'lwdn', axes, time, &
       'downward long-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
       
  id_cosz = register_tiled_diag_field ( module_name, 'coszen', axes, time, &
       'cosine of zenith angle', missing_value=-2.0 )
  id_albedo_dir = register_tiled_diag_field ( module_name, 'albedo_dir', &
       (/id_lon,id_lat,id_band/), time, &
       'land surface albedo for direct light', missing_value=-1.0 )
  id_albedo_dif = register_tiled_diag_field ( module_name, 'albedo_dif', &
       (/id_lon,id_lat,id_band/), time, &
       'land surface albedo for diffuse light', missing_value=-1.0 )
  id_subs_refl_dif = register_tiled_diag_field(module_name, 'subs_refl_dif', &
       (/id_lon, id_lat, id_band/), time, &
       'substrate reflectivity for diffuse light',missing_value=-1.0)
  id_subs_refl_dir = register_tiled_diag_field(module_name, 'subs_refl_dir', &
       (/id_lon, id_lat, id_band/), time, &
       'substrate reflectivity for direct light',missing_value=-1.0)
  id_vegn_refl_dif = register_tiled_diag_field(module_name, 'vegn_refl_dif', &
       (/id_lon, id_lat, id_band/), time, &
       'black-background canopy reflectivity for diffuse light',missing_value=-1.0)
  id_vegn_refl_dir = register_tiled_diag_field(module_name, 'vegn_refl_dir', &
       (/id_lon, id_lat, id_band/), time, &
       'black-background canopy reflectivity for direct light',missing_value=-1.0)
  id_vegn_tran_dif = register_tiled_diag_field(module_name, 'vegn_tran_dif', &
       (/id_lon, id_lat, id_band/), time, &
       'black-background canopy transmittance for diffuse light',missing_value=-1.0)
  id_vegn_tran_dir = register_tiled_diag_field(module_name, 'vegn_tran_dir', &
       (/id_lon, id_lat, id_band/), time, &
       'part of direct light that passes through canopy unscattered',missing_value=-1.0)
  id_vegn_sctr_dir = register_tiled_diag_field(module_name, 'vegn_sctr_dir', &
       (/id_lon, id_lat, id_band/), time, &
       'part of direct light scattered downward by canopy',missing_value=-1.0)
  id_vegn_refl_lw = register_tiled_diag_field ( module_name, 'vegn_refl_lw', axes, time, &
       'canopy reflectivity for thermal radiation', missing_value=-1.0)
  id_vegn_tran_lw = register_tiled_diag_field ( module_name, 'vegn_tran_lw', axes, time, &
       'canopy transmittance for thermal radiation', missing_value=-1.0)

  id_vegn_cover = register_tiled_diag_field ( module_name, 'vegn_cover', axes, time, &
             'fraction covered by vegetation', missing_value=-1.0 )

  id_con_g_h = register_tiled_diag_field ( module_name, 'con_g_h', axes, time, &
       'conductance for sensible heat between ground surface and canopy air', &
       'm/s', missing_value=-1.0 )
  id_grnd_T = register_tiled_diag_field ( module_name, 'Tgrnd', axes, time, &
       'ground surface temperature', 'degK', missing_value=-1.0 )

  id_e_res_1 = register_tiled_diag_field ( module_name, 'e_res_1', axes, time, &
       'canopy air energy residual due to nonlinearities', 'W/m2', missing_value=-1e20)
  id_e_res_2 = register_tiled_diag_field ( module_name, 'e_res_2', axes, time, &
       'canopy energy residual due to nonlinearities', 'W/m2', missing_value=-1e20)

end subroutine land_diag_init

#define DEFINE_LAND_ACCESSOR_0D(xtype,x) subroutine land_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))p=>t%x;end subroutine

DEFINE_LAND_ACCESSOR_0D(real,lwup)
DEFINE_LAND_ACCESSOR_0D(real,e_res_1)
DEFINE_LAND_ACCESSOR_0D(real,e_res_2)

end module land_model_mod
   
