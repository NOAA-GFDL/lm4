! ============================================================================
! lake model module
! ============================================================================
module lake_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : error_mesg, file_exist, read_data, check_nml_error, &
     stdlog, close_file, mpp_pe, mpp_root_pe, FATAL, NOTE
use time_manager_mod, only: time_type_to_real
use diag_manager_mod, only: diag_axis_init
use constants_mod, only: tfreeze, hlv, hlf, dens_h2o, grav, vonkarm, rdgas

use land_constants_mod, only : NBANDS
use lake_tile_mod, only : &
     lake_tile_type, read_lake_data_namelist, &
     lake_data_radiation, &
     lake_data_thermodynamics, &
     cpw,clw,csw, lake_width_inside_lake, large_lake_sill_width, &
     lake_specific_width, n_outlet, outlet_face, outlet_i, outlet_j, outlet_width, &
     new_lake_tile, delete_lake_tile, lake_tile_stock_pe, lake_tile_heat
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_tile_diag_mod, only : register_tiled_static_field, &
     register_tiled_diag_field, send_tile_data, diag_buff_type, &
     send_tile_data_r0d_fptr, add_tiled_static_field_alias, &
     set_default_diag_filter
use land_data_mod, only : lnd, log_version
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, get_tile_data, field_exists
use land_debug_mod, only: is_watch_point, check_var_range, land_error_message
use land_utils_mod, only : put_to_tiles_r0d_fptr
use land_io_mod, only : read_field

implicit none
private

! ==== public interfaces =====================================================
public :: read_lake_namelist
public :: lake_init
public :: lake_end
public :: save_lake_restart

public :: lake_get_sfc_temp
public :: lake_radiation
public :: lake_step_1
public :: lake_step_2

public :: large_dyn_small_stat

public :: lake_abstraction
public :: prohibit_shallow_lake
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'lake'
#include "../shared/version_variable.inc"

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
real    :: init_temp            = 288.        ! cold-start lake T
real    :: init_w               = 1000.      ! cold-start w(l)/dz(l)
logical :: use_rh_feedback      = .true.
logical :: make_all_lakes_wide  = .false.
logical :: large_dyn_small_stat = .true.
logical :: relayer_in_step_one  = .false.
logical :: float_ice_to_top     = .false.
logical :: wind_penetrates_ice  = .false.
real    :: min_rat              = 0.4
logical :: do_stratify          = .true.
character(len=16):: albedo_to_use = ''  ! or 'brdf-params'
real    :: K_z_large            = 1.
real    :: K_z_background       = 0.
real    :: K_z_min              = 0.
real    :: K_z_factor           = 1.
real    :: c_drag               = 1.2e-3
real    :: lake_depth_max       = 1.e10
real    :: lake_depth_min       = 1.99
real    :: max_plain_slope      = -1.e10
logical :: do_lake_abstraction  = .false.
logical, public :: use_reservoir        = .false.
real, public    :: ResMin               = 0.1 !public for river_physics
real, public    :: ResMax               = 0.75
logical :: prohibit_shallowlake = .false.

namelist /lake_nml/ init_temp, init_w,       &
                    use_rh_feedback, cpw, clw, csw, &
                    make_all_lakes_wide, large_dyn_small_stat, &
                    relayer_in_step_one, float_ice_to_top, &
                    min_rat, do_stratify, albedo_to_use, K_z_large, &
		    K_z_background, K_z_min, K_z_factor, &
		    lake_depth_max, lake_depth_min, max_plain_slope, &
        do_lake_abstraction, use_reservoir, ResMin, ResMax, prohibit_shallowlake
!---- end of namelist --------------------------------------------------------
real    :: K_z_molec            = 1.4e-7
real    :: tc_molec             = 0.59052 ! dens_h2o*clw*K_z_molec
real    :: tc_molec_ice         = 2.5

logical         :: module_is_initialized =.FALSE.
logical         :: use_brdf
real            :: delta_time

integer         :: num_l              ! # of water layers
real, allocatable:: zfull (:)    ! diag axis, dimensionless layer number
real, allocatable:: zhalf (:)
real            :: max_rat
logical, public :: is_rsv_restart = .false. !public for lake transitions

! ---- diagnostic field IDs
integer :: id_lwc, id_swc, id_temp
integer :: id_evap, id_dz, id_wl, id_ws, id_K_z, id_silld, id_sillw, id_backw
integer :: id_Afrac_rsv, id_Vfrac_rsv, id_rsv_depth
integer :: id_sub_lmass, id_sub_fmass, id_sub_heat, id_sub_cmass
integer :: id_back1
integer :: id_lake_area, id_lake_frac
! ==== end of module variables ===============================================

contains

! ============================================================================
subroutine read_lake_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l

  call read_lake_data_namelist(num_l)

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=lake_nml, iostat=io)
     ierr = check_nml_error(io, 'lake_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=lake_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'lake_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=lake_nml)
  endif

  ! ---- set up vertical discretization
  allocate (zhalf(num_l+1), zfull(num_l))
  zhalf(1) = 0
  do l = 1, num_l
     zhalf(l+1) = zhalf(l) + 1.
     zfull(l) = 0.5*(zhalf(l+1) + zhalf(l))
  enddo

  max_rat = 1. / min_rat
  if (trim(albedo_to_use)=='brdf-params') then
     use_brdf = .true.
  else if (trim(albedo_to_use)=='') then
     use_brdf = .false.
  else
     call error_mesg('lake_init',&
          'option albedo_to_use="'// trim(albedo_to_use)//&
          '" is invalid, use "brdf-params", or nothing ("")',&
          FATAL)
  endif

end subroutine read_lake_namelist


! ============================================================================
! initialize lake model
subroutine lake_init ( id_ug )
  integer,intent(in) :: id_ug !<Unstructured axis id.

  ! ---- local vars
  type(land_tile_enum_type)     :: ce   ! tile list enumerator
  type(land_tile_type), pointer :: tile ! pointer to current tile
  type(land_restart_type) :: restart
  logical :: restart_exists
  real, allocatable :: buffer(:),bufferc(:),buffert(:)
  integer :: i, g, l
  logical :: river_data_exist
  character(*), parameter :: restart_file_name = 'INPUT/lake.res.nc'

  module_is_initialized = .TRUE.
  delta_time = time_type_to_real(lnd%dt_fast)

  allocate(buffer (lnd%ls:lnd%le))
  allocate(bufferc(lnd%ls:lnd%le))
  allocate(buffert(lnd%ls:lnd%le))
  buffer (:) = 0
  bufferc(:) = 0
  buffert(:) = 0

  river_data_exist = file_exist('INPUT/river_data.nc', lnd%sg_domain)
  if (river_data_exist) then
     call error_mesg('lake_init', 'reading lake information from river data file', NOTE)
  else
     call error_mesg('lake_init', 'river data file not present: lake fraction is set to zero', NOTE)
  endif

  IF (LARGE_DYN_SMALL_STAT) THEN

     if (river_data_exist) call read_data('INPUT/river_data.nc', 'connected_to_next', bufferc(:), lnd%sg_domain, lnd%ug_domain)
     call put_to_tiles_r0d_fptr(bufferc, land_tile_map, lake_connected_to_next_ptr)

     if (river_data_exist) call read_data('INPUT/river_data.nc', 'whole_lake_area', buffer(:), lnd%sg_domain, lnd%ug_domain)
     call put_to_tiles_r0d_fptr(buffer, land_tile_map, lake_whole_area_ptr)

     if (river_data_exist) call read_data('INPUT/river_data.nc', 'lake_depth_sill', buffer(:),  lnd%sg_domain, lnd%ug_domain)
     buffer = min(buffer, lake_depth_max)
     buffer = max(buffer, lake_depth_min)
     call put_to_tiles_r0d_fptr(buffer,  land_tile_map, lake_depth_sill_ptr)

     ! lake_tau is just used here as a flag for 'large lakes'
     ! sill width of -1 is a flag saying not to allow transient storage
     if (river_data_exist) call read_data('INPUT/river_data.nc', 'lake_tau', buffert(:),  lnd%sg_domain, lnd%ug_domain)
     buffer = -1.
     !where (bufferc.gt.0.5) buffer = lake_width_inside_lake
     where (bufferc.lt.0.5 .and. buffert.gt.1.) buffer = large_lake_sill_width
     if (lake_specific_width) then
         do i = 1, n_outlet
           g = outlet_j(i)*lnd%nlon + outlet_i(i)
           if(lnd%ug_face.eq.outlet_face(i).and.lnd%gs.le.g.and.lnd%ge.ge.g) then
             l = lnd%l_index(g)
             buffer(l) = outlet_width(i)
           endif
         enddo
     endif
     call put_to_tiles_r0d_fptr(buffer, land_tile_map, lake_width_sill_ptr)

     buffer = 1.e8
     if (river_data_exist .and. max_plain_slope.gt.0.) &
        call read_data('INPUT/river_data.nc', 'max_slope_to_next', buffer(:), lnd%sg_domain, lnd%ug_domain)
     if (river_data_exist) call read_data('INPUT/river_data.nc', 'travel', buffert(:), lnd%sg_domain, lnd%ug_domain)
     bufferc = 0.
     where (buffer.lt.max_plain_slope .and. buffert.gt.1.5) bufferc = 1.
     call put_to_tiles_r0d_fptr(bufferc, land_tile_map, lake_backwater_ptr)
     bufferc = 0
     where (buffer.lt.max_plain_slope .and. buffert.lt.1.5) bufferc = 1.
     call put_to_tiles_r0d_fptr(bufferc, land_tile_map, lake_backwater_1_ptr)

  ELSE
     if (river_data_exist) then
        call read_data('INPUT/river_data.nc', 'whole_lake_area', bufferc(:), lnd%sg_domain, lnd%ug_domain)
        call read_data('INPUT/river_data.nc', 'lake_depth_sill', buffer(:), lnd%sg_domain, lnd%ug_domain)
     endif
     where (bufferc.eq.0.)                      buffer = 0.
     where (bufferc.gt.0..and.bufferc.lt.2.e10) buffer = max(2., 2.5e-4*sqrt(bufferc))
     call put_to_tiles_r0d_fptr(buffer,  land_tile_map, lake_depth_sill_ptr)
     call put_to_tiles_r0d_fptr(bufferc, land_tile_map, lake_whole_area_ptr)

     buffer = 4. * buffer
     where (bufferc.gt.2.e10) buffer = min(buffer, 60.)
     if (river_data_exist) call read_data('INPUT/river_data.nc', 'connected_to_next', bufferc(:), lnd%sg_domain, lnd%ug_domain)
     call put_to_tiles_r0d_fptr(bufferc, land_tile_map, lake_connected_to_next_ptr)

     where (bufferc.gt.0.5) buffer=lake_width_inside_lake
     if (make_all_lakes_wide) buffer = lake_width_inside_lake
     call put_to_tiles_r0d_fptr(buffer, land_tile_map, lake_width_sill_ptr)
  ENDIF

  deallocate (buffer, bufferc, buffert)

  ! -------- initialize lake state --------
  ce = first_elmt(land_tile_map)
  do while (loop_over_tiles(ce,tile))
     if (.not.associated(tile%lake)) cycle

     tile%lake%dz = tile%lake%pars%depth_sill/num_l
     if (init_temp.ge.tfreeze) then
        tile%lake%wl = init_w*tile%lake%dz
        tile%lake%ws = 0
     else
        tile%lake%wl = 0
        tile%lake%ws = init_w*tile%lake%dz
     endif
     tile%lake%T             = init_temp
     tile%lake%sub_lmass     = 0.
     tile%lake%sub_fmass     = 0.
     tile%lake%sub_heat      = 0.
     tile%lake%sub_cmass     = 0.
     !these three vars will be initialized later in lake_transitions_init
     tile%lake%Afrac_rsv     = 0.
     tile%lake%Vfrac_rsv     = 0.
     tile%lake%rsv_depth     = 0.
  enddo

  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then
     call error_mesg('lake_init', 'reading NetCDF restart "'//trim(restart_file_name)//'"', NOTE)
     if (field_exists(restart,'dz')) &
        call get_tile_data(restart, 'dz', 'zfull', lake_dz_ptr)
     call get_tile_data(restart, 'temp', 'zfull', lake_temp_ptr)
     call get_tile_data(restart, 'wl',   'zfull', lake_wl_ptr)
     call get_tile_data(restart, 'ws',   'zfull', lake_ws_ptr)
     if (field_exists(restart,'Afrac_rsv')) &
        call get_tile_data(restart, 'Afrac_rsv', lake_Afrac_rsv_ptr)
     if (field_exists(restart,'Vfrac_rsv')) &
        call get_tile_data(restart, 'Vfrac_rsv', lake_Vfrac_rsv_ptr)
     if (field_exists(restart,'depth_rsv')) &
        call get_tile_data(restart, 'depth_rsv', lake_depth_rsv_ptr)
     if (field_exists(restart,'sub_lmass')) &
        call get_tile_data(restart, 'sub_lmass', lake_sub_lmass_ptr)
     if (field_exists(restart,'sub_fmass')) &
        call get_tile_data(restart, 'sub_fmass', lake_sub_fmass_ptr)
     if (field_exists(restart,'sub_heat')) &
        call get_tile_data(restart, 'sub_heat',  lake_sub_heat_ptr)
     if (field_exists(restart,'sub_cmass')) &
        call get_tile_data(restart, 'sub_cmass', lake_sub_cmass_ptr)
  else
     call error_mesg('lake_init', 'cold-starting lake', NOTE)
  endif

  if(field_exists(restart,'Afrac_rsv').and.field_exists(restart,'Vfrac_rsv').and.field_exists(restart,'depth_rsv'))then
    is_rsv_restart=.true.
  endif

  call free_land_restart(restart)

  call lake_diag_init(id_ug)
  ! ---- static diagnostic section
  call send_tile_data_r0d_fptr(id_sillw, lake_width_sill_ptr)
  call send_tile_data_r0d_fptr(id_silld, lake_depth_sill_ptr)
  call send_tile_data_r0d_fptr(id_backw, lake_backwater_ptr)
  call send_tile_data_r0d_fptr(id_back1, lake_backwater_1_ptr)
end subroutine lake_init


! ============================================================================
subroutine lake_end ()

  deallocate (zfull, zhalf)
  module_is_initialized =.FALSE.

end subroutine lake_end


! ============================================================================
subroutine save_lake_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object

  call error_mesg('lake_end','writing NetCDF restart',NOTE)
! must set domain so that io_domain is available
! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp)//'lake.res.nc'
  call init_land_restart(restart, filename, lake_tile_exists, tile_dim_length)
  call add_restart_axis(restart,'zfull',zfull(1:num_l),'Z','m','full level',sense=-1)

  ! write out fields
  call add_tile_data(restart,'dz',   'zfull', lake_dz_ptr,   'layer thickness','m')
  call add_tile_data(restart,'temp', 'zfull', lake_temp_ptr, 'lake temperature','degrees_K')
  call add_tile_data(restart,'wl',   'zfull', lake_wl_ptr,   'liquid water content','kg/m2')
  call add_tile_data(restart,'ws',   'zfull', lake_ws_ptr,   'solid water content','kg/m2')
  call add_tile_data(restart,'Afrac_rsv', lake_Afrac_rsv_ptr, 'area fraction of reservoir to the lake tile', 'unitless')
  call add_tile_data(restart,'Vfrac_rsv', lake_Vfrac_rsv_ptr, 'volume fraction of reservoir to the lake tile', 'unitless')
  call add_tile_data(restart,'depth_rsv', lake_depth_rsv_ptr, 'reservoir construction depth', 'm')
  call add_tile_data(restart,'sub_lmass', lake_sub_lmass_ptr, 'buried liquid water under lake due to reservoir building', 'kg/m2')
  call add_tile_data(restart,'sub_fmass', lake_sub_fmass_ptr, 'buried frozen water under lake due to reservoir building', 'kg/m2')
  call add_tile_data(restart,'sub_heat',  lake_sub_heat_ptr,  'buried heat under lake due to reservoir building', 'J/m2')
  call add_tile_data(restart,'sub_cmass', lake_sub_cmass_ptr, 'buried carbon under lake due to reservoir building', 'kg C/m2')

  ! save performs io domain aggregation through mpp_io as with regular domain data
  call save_land_restart(restart)
  call free_land_restart(restart)
end subroutine save_lake_restart

! ============================================================================
subroutine lake_get_sfc_temp(lake, lake_T)
  type(lake_tile_type), intent(in) :: lake
  real, intent(out) :: lake_T

  lake_T = lake%T(1)
end subroutine lake_get_sfc_temp


! ============================================================================
! compute lake-only radiation properties
subroutine lake_radiation ( lake, cosz, &
     lake_refl_dir, lake_refl_dif, lake_refl_lw, lake_emis )
  type(lake_tile_type), intent(in) :: lake
  real, intent(in) :: cosz
  real, intent(out) :: lake_refl_dir(NBANDS), lake_refl_dif(NBANDS), lake_refl_lw, lake_emis

  call lake_data_radiation ( lake, cosz, use_brdf, lake_refl_dir, lake_refl_dif, lake_emis )
  lake_refl_lw = 1 - lake_emis
end subroutine lake_radiation


! ============================================================================
! update lake properties explicitly for time step.
! integrate lake-heat conduction equation upward from bottom of lake
! to surface, delivering linearization of surface ground heat flux.
subroutine lake_step_1 ( u_star_a, p_surf, latitude, lake, &
                         lake_T, &
                         lake_rh, lake_liq, lake_ice, lake_subl, lake_tf, lake_G0, &
                         lake_DGDT )

  real, intent(in)   :: u_star_a, p_surf, latitude
  type(lake_tile_type), intent(inout) :: lake
  real, intent(out)  :: &
       lake_T, &
       lake_rh, lake_liq, lake_ice, lake_subl, &
       lake_tf, & ! freezing temperature of lake, degK
       lake_G0, &
       lake_DGDT

  ! ---- local vars
  real                  :: bbb, denom, dt_e, tc_dz_eff
  real                  :: z_cum, z_mid, dz_mid, rho_t_mid, k_neutral
  real, dimension(num_l):: aaa, ccc, thermal_cond, heat_capacity, dz_alt, &
                            z_alt, rho_t
  integer               :: l
  real                  :: k_star, N_sq, Ri, u_star, z_liq, z_ice, rho_a
  real                  :: lake_depth
! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  if(is_watch_point()) then
     write(*,*) 'lake_step_1 checkpoint 1'
     write(*,*) 'mask    ', .true.
     write(*,*) 'T       ', lake_T
     write(*,*) 'rh      ', lake_rh
     write(*,*) 'liq     ', lake_liq
     write(*,*) 'ice     ', lake_ice
     write(*,*) 'subl    ', lake_subl
     write(*,*) 'G0      ', lake_G0
     write(*,*) 'DGDT    ', lake_DGDT
    do l = 1, num_l
      write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
                 ' dz=', lake%dz(l),&
                 ' T =', lake%T(l),&
                 ' wl=', lake%wl(l),&
                 ' ws=', lake%ws(l), &
                 'K_z=', lake%K_z(l)
      enddo
  endif


  if (relayer_in_step_one) call lake_relayer ( lake )

  lake%K_z = 0.
  lake_T = lake%T(1)
  if (use_rh_feedback) then
      lake_depth = (sum(lake%wl(:))+sum(lake%ws(:))) / DENS_H2O
    else
      lake_depth = lake%pars%depth_sill
    endif
  call lake_data_thermodynamics ( lake%pars, lake_depth, lake_rh, &
                                  lake%heat_capacity_dry, thermal_cond )
! Ignore air humidity in converting atmospheric friction velocity to lake value
  rho_a = p_surf/(rdgas*lake_T)
! No momentum transfer through ice cover
  if (lake%ws(1).le.0. .or. wind_penetrates_ice) then
     u_star = u_star_a*sqrt(rho_a/dens_h2o)
     k_star = 2.79e-5*sqrt(sin(abs(latitude)))*u_star**(-1.84)
     k_star = k_star*(c_drag/1.2e-3)**1.84
  else
     u_star = 0.
     k_star = 1.
  endif
! k_star from B. Henderson-Sellers (1985, Appl. Math. Mod., 9)
!  k_star = 2.79e-5*sqrt(sin(abs(latitude)))*u_star**(-1.84)
  z_cum = 0.
  do l = 1, num_l
    heat_capacity(l) = lake%heat_capacity_dry(l) * lake%dz(l) &
            + clw*lake%wl(l) + csw*lake%ws(l)
    dz_alt(l) = (lake%wl(l) + lake%ws(l))/dens_h2o
    z_alt(l) = z_cum + 0.5*dz_alt(l)
    z_cum = z_cum + dz_alt(l)
! rho_t from hostetler and bartlein (1990), citing Heggen (1983)
! is a call available in fms?
    rho_t(l) = 1. - 1.9549e-5*abs(lake%T(l)-277.)**1.68
    enddo

  lake_liq  = max(lake%wl(1), 0.)
  lake_ice  = max(lake%ws(1), 0.)
  if (lake_ice > 0) then
     lake_subl = 1
  else
     lake_subl = 0
  endif

  if(num_l > 1) then
    if (do_stratify) then
        do l = 1, num_l-1
          if (lake%ws(l).le.0..and.lake%ws(l+1).le.0.) then
              dz_mid = z_alt(l+1)-z_alt(l)
              z_mid = 0.5 * (z_alt(l)+z_alt(l+1))
              rho_t_mid = 0.5*(rho_t(l)+rho_t(l+1))
              if (k_star*z_mid .lt. 10.) then
                  k_neutral = vonkarm * u_star * z_mid * exp (-k_star * z_mid)
                else
                  k_neutral = 0.
                endif
              N_sq = (grav / rho_t_mid) * (rho_t(l+1)-rho_t(l)) /dz_mid
              if (N_sq .gt. 0. .and. k_neutral.ne.0.) then
                  ! stability function from B. Henderson-Sellers (1985)
                  Ri = 0.05*(-1. + sqrt(1.+40.*N_sq*(vonkarm*z_mid/u_star)**2 &
                                                   *exp(2.*k_star*z_mid)))
                  lake%K_z(l) = k_neutral / (1. + 37.*Ri*Ri) + K_z_molec
                else if (k_neutral.eq.0.) then
                  lake%K_z(l) = K_z_molec
                else  ! arbitrary constant for unstable mixing
                  lake%K_z(l) = K_z_large
                endif
	      if (lake%pars%depth_sill.gt.2.01) &
	          lake%K_z(l) = K_z_factor &
		   * max(lake%K_z(l) + K_z_background, K_z_min)
              aaa(l+1) = - lake%K_z(l) * delta_time / (dz_alt(l+1)*dz_mid)
              ccc(l)   = - lake%K_z(l) * delta_time / (dz_alt(l  )*dz_mid)
            else
              z_liq = 0.5*(lake%wl(l)+lake%wl(l+1))/dens_h2o
              z_ice = 0.5*(lake%ws(l)+lake%ws(l+1))/dens_h2o
              tc_dz_eff = 1. / (z_liq/tc_molec + z_ice/tc_molec_ice)
              aaa(l+1) = - tc_dz_eff * delta_time / heat_capacity(l+1)
              ccc(l)   = - tc_dz_eff * delta_time / heat_capacity(l)
            endif
          enddo
      else
        do l = 1, num_l-1
          tc_dz_eff = 2 / ( dz_alt(l+1)/thermal_cond(l+1) &
             + dz_alt(l)/thermal_cond(l)   )
          aaa(l+1) = - tc_dz_eff * delta_time / heat_capacity(l+1)
          ccc(l)   = - tc_dz_eff * delta_time / heat_capacity(l)
          enddo
      endif

     bbb = 1.0 - aaa(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(lake%T(num_l) - lake%T(num_l-1)) &
               + lake%geothermal_heat_flux * delta_time / heat_capacity(num_l)
     lake%e(num_l-1) = -aaa(num_l)/denom
     lake%f(num_l-1) = dt_e/denom
     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*lake%e(l)
        dt_e = - ( ccc(l)*(lake%T(l+1) - lake%T(l)  ) &
                  -aaa(l)*(lake%T(l)   - lake%T(l-1)) )
        lake%e(l-1) = -aaa(l)/denom
        lake%f(l-1) = (dt_e - ccc(l)*lake%f(l))/denom
     end do
     denom = delta_time/(heat_capacity(1) )
     lake_G0    = ccc(1)*(lake%T(2)- lake%T(1) &
          + lake%f(1)) / denom
     lake_DGDT  = (1 - ccc(1)*(1-lake%e(1))) / denom
  else  ! one-level case
     denom = delta_time/heat_capacity(1)
     lake_G0    = 0.
     lake_DGDT  = 1. / denom
  end if

  ! set the freezing temperature of the lake
  lake_tf = tfreeze

  if(is_watch_point()) then
     write(*,*) 'lake_step_1 checkpoint 2'
     write(*,*) 'mask    ', .true.
     write(*,*) 'T       ', lake_T
     write(*,*) 'rh      ', lake_rh
     write(*,*) 'liq     ', lake_liq
     write(*,*) 'ice     ', lake_ice
     write(*,*) 'subl    ', lake_subl
     write(*,*) 'G0      ', lake_G0
     write(*,*) 'DGDT    ', lake_DGDT
    do l = 1, num_l
      write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
                 ' dz=', lake%dz(l),&
                 ' T =', lake%T(l),&
                 ' wl=', lake%wl(l),&
                 ' ws=', lake%ws(l), &
                 'K_z=', lake%K_z(l)
      enddo
  endif

end subroutine lake_step_1


! ============================================================================
! apply boundary flows to lake water and move lake water vertically.
  subroutine lake_step_2 ( lake, diag, lake_subl, snow_lprec, snow_hlprec,  &
                           subs_DT, subs_M_imp, subs_evap, &
                           use_tfreeze_in_grnd_latent, &
                           lake_levap, lake_fevap, lake_melt, &
                           lake_Ttop, lake_Ctop )
  type(lake_tile_type), intent(inout) :: lake
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: &
     lake_subl     !
  real, intent(in) :: &
     snow_lprec, &
     snow_hlprec, &
     subs_DT,       &!
     subs_M_imp,       &! rate of phase change of non-evaporated lake water
     subs_evap
  logical, intent(in) :: use_tfreeze_in_grnd_latent
  real, intent(out) :: &
     lake_levap, lake_fevap, lake_melt, &
     lake_Ttop, lake_Ctop

  ! ---- local vars
  real, dimension(num_l) :: del_t, dW_l
  real, dimension(num_l+1) :: flow
  real :: ice_to_move, h_upper, h_lower, h_to_move_up, &
     hcap, dheat, &
     melt_per_deg, melt
  real :: jj
  integer :: l
  real :: v0, v1

  jj = 1.

  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 1 ***** '
    write(*,*) 'mask    ', .true.
    write(*,*) 'subs_evap    ', subs_evap
    write(*,*) 'snow_lprec   ', snow_lprec
    write(*,*) 'subs_M_imp   ', subs_M_imp
    write(*,*) 'theta_s ', lake%pars%w_sat
    do l = 1, num_l
      write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
                 ' dz=', lake%dz(l),&
                 ' T =', lake%T(l),&
                 ' Th=', (lake%ws(l) &
                         +lake%wl(l))/(dens_h2o*lake%dz(l)),&
                 ' wl=', lake%wl(l),&
                 ' ws=', lake%ws(l)
      enddo
  endif

  ! ---- record fluxes ---------
  lake_levap  = subs_evap*(1-lake_subl)
  lake_fevap  = subs_evap*   lake_subl
  lake_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution --------------
  del_t(1) = subs_DT
  lake%T(1) = lake%T(1) + del_t(1)
  if ( num_l > 1) then
    do l = 1, num_l-1
      del_t(l+1) = lake%e(l) * del_t(l) + lake%f(l)
      lake%T(l+1) = lake%T(l+1) + del_t(l+1)
    end do
  end if

  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 2 ***** '
    do l = 1, num_l
       write(*,*) 'level=', l, 'T', lake%T(l)
    enddo
  endif

  v0 = sum(lake%wl+lake%ws)/DENS_H2O !m

  ! ---- extract evap from lake and do implicit melt --------------------
  lake%wl(1) = lake%wl(1) - lake_levap*delta_time
  lake%ws(1) = lake%ws(1) - lake_fevap*delta_time
  hcap = lake%heat_capacity_dry(1)*lake%dz(1) &
                     + clw*lake%wl(1) + csw*lake%ws(1)
  ! T adjustment for nonlinear terms (del_T)*(del_W)
  dheat = delta_time*(clw*lake_levap+csw*lake_fevap)*del_T(1)
  ! take out extra heat not claimed in advance for evaporation
  if (use_tfreeze_in_grnd_latent) dheat = dheat &
          - delta_time*((cpw-clw)*lake_levap+(cpw-csw)*lake_fevap) &
                                 *(lake%T(1)-del_T(1)-tfreeze)
  lake%T(1)  = lake%T(1)  + dheat/hcap
  lake%wl(1) = lake%wl(1) + subs_M_imp
  lake%ws(1) = lake%ws(1) - subs_M_imp
  lake%T(1)  = tfreeze + (hcap*(lake%T(1)-tfreeze) ) &
                            / ( hcap + (clw-csw)*subs_M_imp )

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 2.1 ***** '
     do l = 1, num_l
        write(*,*) 'level=', l, 'T', lake%T(l)
     enddo
  endif

  ! ---- remainder of mass fluxes and associated sensible heat fluxes --------
  ! note that only liquid inputs are received by  lake from snow pack. any
  ! snow fall just creates a snow pack on top  of lake, even if lake is not
  ! frozen. but snow pack on top of unfrozen lake will interact thermally,
  ! so that either lake freezes or snow melts and falls in.
    flow=1
    flow(1)  = snow_lprec *delta_time
    do l = 1, num_l
      flow(l+1) = 0
      dW_l(l) = flow(l) - flow(l+1)
      lake%wl(l) = lake%wl(l) + dW_l(l)
    enddo

  v1 = sum(lake%wl+lake%ws)/DENS_H2O !m
  if(use_reservoir) lake%Vfrac_rsv = (lake%Vfrac_rsv*v0 + (snow_lprec-subs_evap)*lake%Afrac_rsv*delta_time/DENS_H2O) / v1

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 3.3 ***** '
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' wl=', lake%wl(l),&
             'flow=', flow(l)
     enddo
  endif

  hcap = lake%heat_capacity_dry(1)*lake%dz(1) &
                     + clw*(lake%wl(1)-dW_l(1)) + csw*lake%ws(1)
  lake%T(1) = tfreeze + (hcap*(lake%T(1)-tfreeze) +  &
                                 snow_hlprec*delta_time) &
                            / ( hcap + clw*dW_l(1) )

  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 3.4 ***** '
    write(*,*) ' tfreeze', tfreeze
    write(*,*) ' snow_hlprec', snow_hlprec
  endif

  do l = 1, num_l
    ! ---- compute explicit melt/freeze --------------------------------------
    hcap = lake%heat_capacity_dry(l)*lake%dz(l) &
             + clw*lake%wl(l) + csw*lake%ws(l)
    melt_per_deg = hcap/hlf
    if (lake%ws(l)>0 .and. lake%T(l)>tfreeze) then
      melt =  min(lake%ws(l), (lake%T(l)-tfreeze)*melt_per_deg)
    else if (lake%wl(l)>0 .and. lake%T(l)<tfreeze) then
      melt = -min(lake%wl(l), (tfreeze-lake%T(l))*melt_per_deg)
    else
      melt = 0
    endif
    lake%wl(l) = lake%wl(l) + melt
    lake%ws(l) = lake%ws(l) - melt
    lake%T(l) = tfreeze &
       + (hcap*(lake%T(l)-tfreeze) - hlf*melt) &
                            / ( hcap + (clw-csw)*melt )
    lake_melt = lake_melt + melt / delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 5 ***** '
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' dz=', lake%dz(l),&
             ' T =', lake%T(l),&
             ' Th=', (lake%ws(l) +lake%wl(l))/(dens_h2o*lake%dz(l)),&
             ' wl=', lake%wl(l),&
             ' ws=', lake%ws(l)
     enddo
  endif

  if (.not.relayer_in_step_one) call lake_relayer ( lake )

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 6 ***** '
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' dz=', lake%dz(l),&
             ' T =', lake%T(l),&
             ' Th=', (lake%ws(l) +lake%wl(l))/(dens_h2o*lake%dz(l)),&
             ' wl=', lake%wl(l),&
             ' ws=', lake%ws(l)
     enddo
  endif

  if (float_ice_to_top) then
      do l = num_l, 2, -1
        if (lake%ws(l) .gt. 0. .and. lake%wl(l-1) .gt. 0.) then
            ice_to_move = min(lake%ws(l), lake%wl(l-1))
            h_upper = (clw*lake%wl(l-1)+csw*lake%ws(l-1))*lake%T(l-1)
            h_lower = (clw*lake%wl(l  )+csw*lake%ws(l  ))*lake%T(l  )
            lake%wl(l-1) = lake%wl(l-1) - ice_to_move
            lake%ws(l-1) = lake%ws(l-1) + ice_to_move
            lake%wl(l  ) = lake%wl(l  ) + ice_to_move
            lake%ws(l  ) = lake%ws(l  ) - ice_to_move
            h_to_move_up = ice_to_move*(csw*lake%T(l)-clw*lake%T(l-1))
            h_upper  = h_upper + h_to_move_up
            h_lower  = h_lower - h_to_move_up
            lake%T(l-1) = h_upper / (clw*lake%wl(l-1)+csw*lake%ws(l-1))
            lake%T(l  ) = h_lower / (clw*lake%wl(l  )+csw*lake%ws(l  ))
          endif
        enddo
    endif

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 7 ***** '
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' dz=', lake%dz(l),&
             ' T =', lake%T(l),&
             ' Th=', (lake%ws(l) +lake%wl(l))/(dens_h2o*lake%dz(l)),&
             ' wl=', lake%wl(l),&
             ' ws=', lake%ws(l)
     enddo
  endif


  lake_Ttop = lake%T(1)
  lake_Ctop = lake%heat_capacity_dry(1)*lake%dz(1) &
       + clw*lake%wl(1) + csw*lake%ws(1)

! ----------------------------------------------------------------------------
! given solution for surface energy balance, write diagnostic output.
!

  ! ---- diagnostic section
  call send_tile_data (id_dz,   lake%dz,     diag )
  call send_tile_data (id_temp, lake%T,     diag )
  call send_tile_data (id_wl,  lake%wl(1:num_l), diag )
  call send_tile_data (id_ws,  lake%ws(1:num_l), diag )
  call send_tile_data (id_lwc,  lake%wl(1:num_l)/lake%dz(1:num_l), diag )
  call send_tile_data (id_swc,  lake%ws(1:num_l)/lake%dz(1:num_l), diag )
  call send_tile_data (id_K_z,  lake%K_z(1:num_l),        diag )
  call send_tile_data (id_evap, lake_levap+lake_fevap, diag )

  call send_tile_data (id_Afrac_rsv, lake%Afrac_rsv, diag)
  call send_tile_data (id_Vfrac_rsv, lake%Vfrac_rsv, diag)
  call send_tile_data (id_rsv_depth, lake%rsv_depth, diag)
  call send_tile_data (id_sub_lmass, lake%sub_lmass, diag)
  call send_tile_data (id_sub_fmass, lake%sub_fmass, diag)
  call send_tile_data (id_sub_heat,  lake%sub_heat,  diag)
  call send_tile_data (id_sub_cmass, lake%sub_cmass, diag)

end subroutine lake_step_2


! ============================================================================
!
  subroutine lake_relayer ( lake )
  type(lake_tile_type), intent(inout) :: lake

  ! ---- local vars
  integer :: l, l_lowest_thin_layer, l_highest_thick_layer
  real :: new_dz, new_ws, new_wl, new_h, new_T, liq_frac

! now check whether we need to re-layer the lake.
  if ( (lake%wl(1)+lake%ws(1)) &
      /(lake%wl(2)+lake%ws(2)) .gt. max_rat) then
      ! top layer has grown too thick. join two lower layers, and split
      ! top layer into two layers. in special case, just join and
      ! re-split top two layers.
      l_lowest_thin_layer = num_l
      do l = 2, num_l-1
        if (lake%dz(l).lt.0.99*lake%dz(num_l)) l_lowest_thin_layer = l
        enddo
      if (l_lowest_thin_layer.gt.2) then
          new_dz = lake%dz(l_lowest_thin_layer) &
                  + lake%dz(l_lowest_thin_layer-1)
          new_wl = lake%wl(l_lowest_thin_layer) &
                  + lake%wl(l_lowest_thin_layer-1)
          new_ws = lake%ws(l_lowest_thin_layer) &
                  + lake%ws(l_lowest_thin_layer-1)
          new_h  = ( clw*lake%wl(l_lowest_thin_layer) &
                   + csw*lake%ws(l_lowest_thin_layer)) &
                   *     lake%T(l_lowest_thin_layer) &
                 + ( clw*lake%wl(l_lowest_thin_layer-1) &
                   + csw*lake%ws(l_lowest_thin_layer-1)) &
                   *     lake%T(l_lowest_thin_layer-1)
          new_T = new_h / (clw*new_wl+csw*new_ws)
          lake%dz(l_lowest_thin_layer) = new_dz
          lake%wl(l_lowest_thin_layer) = new_wl
          lake%ws(l_lowest_thin_layer) = new_ws
          lake%T(l_lowest_thin_layer)  = new_T
          do l = l_lowest_thin_layer-1, 3, -1
            lake%dz(l) = lake%dz(l-1)
            lake%wl(l) = lake%wl(l-1)
            lake%ws(l) = lake%ws(l-1)
            lake%T(l)  = lake%T(l-1)
            enddo
          liq_frac = lake%wl(1) / (lake%wl(1)+lake%ws(1))
          lake%wl(2) =     liq_frac *DENS_H2O*lake%dz(2)
          lake%ws(2) = (1.-liq_frac)*DENS_H2O*lake%dz(2)
          lake%T(2)  = lake%T(1)
          lake%dz(1) = lake%dz(2)
          lake%wl(1) = lake%wl(1) - lake%wl(2)
          lake%ws(1) = lake%ws(1) - lake%ws(2)
        else
          new_wl = lake%wl(1) + lake%wl(2)
          new_ws = lake%ws(1) + lake%ws(2)
          new_h  = ( clw*lake%wl(1) + csw*lake%ws(1))   &
                                                 * lake%T(1) &
                 + ( clw*lake%wl(2) + csw*lake%ws(2))    &
                                                 * lake%T(2)
          new_T  = new_h / (clw*new_wl+csw*new_ws)
          liq_frac = new_wl / (new_wl+new_ws)
          lake%dz(2) = lake%dz(3)
          lake%wl(2) =     liq_frac *DENS_H2O*lake%dz(2)
          lake%ws(2) = (1.-liq_frac)*DENS_H2O*lake%dz(2)
          lake%T(2)  = new_T
          lake%dz(1) = lake%dz(2)
          lake%wl(1) = new_wl - lake%wl(2)
          lake%ws(1) = new_ws - lake%ws(2)
          lake%T(1)  = new_T
        endif
    else if(  (lake%wl(1)+lake%ws(1)) &
             /(lake%wl(2)+lake%ws(2)) .lt. min_rat) then
      ! top layer has grown too thin. join with next layer down, and split
      ! a lower layer to maintain number of layers.  in special case, just
      ! join and re-split top two layers.
      l_highest_thick_layer = 2
      do l = num_l, 3, -1
        if (lake%dz(l).gt.1.01*lake%dz(2)) l_highest_thick_layer = l
        enddo
      new_wl = lake%wl(1) + lake%wl(2)
      new_ws = lake%ws(1) + lake%ws(2)
      new_h  = ( clw*lake%wl(1) + csw*lake%ws(1))   &
                                             * lake%T(1) &
             + ( clw*lake%wl(2) + csw*lake%ws(2))    &
                                             * lake%T(2)
      new_T  = new_h / (clw*new_wl+csw*new_ws)
      if (l_highest_thick_layer.gt.2) then
          lake%dz(1) = lake%dz(2)
          lake%wl(1) = new_wl
          lake%ws(1) = new_ws
          lake%T(1)  = new_T
          do l = 2, l_highest_thick_layer-2
            lake%dz(l) = lake%dz(l+1)
            lake%wl(l) = lake%wl(l+1)
            lake%ws(l) = lake%ws(l+1)
            lake%T(l)  = lake%T(l+1)
            enddo
          new_dz = lake%dz(l_highest_thick_layer) / 2.
          new_wl = lake%wl(l_highest_thick_layer) / 2.
          new_ws = lake%ws(l_highest_thick_layer) / 2.
          new_T  = lake%T(l_highest_thick_layer)
          lake%dz(l_highest_thick_layer-1) = new_dz
          lake%wl(l_highest_thick_layer-1) = new_wl
          lake%ws(l_highest_thick_layer-1) = new_ws
          lake%T(l_highest_thick_layer-1)  = new_T
          lake%dz(l_highest_thick_layer) = new_dz
          lake%wl(l_highest_thick_layer) = new_wl
          lake%ws(l_highest_thick_layer) = new_ws
          lake%T(l_highest_thick_layer)  = new_T
        else
          liq_frac = new_wl / (new_wl+new_ws)
          lake%dz(2) = lake%dz(3) / 2.
          lake%wl(2) =     liq_frac *DENS_H2O*lake%dz(2)
          lake%ws(2) = (1.-liq_frac)*DENS_H2O*lake%dz(2)
          lake%T(2)  = new_T
          lake%dz(1) = lake%dz(2)
          lake%wl(1) = new_wl - lake%wl(2)
          lake%ws(1) = new_ws - lake%ws(2)
          lake%T(1)  = new_T
        endif
    endif
  end subroutine lake_relayer
! ============================================================================
  subroutine  lake_relayer2 (lake_wl, lake_ws, lake_T, lake_dz)
     real, dimension(num_l), intent(inout) ::  lake_wl, lake_ws, lake_T, lake_dz
     type(lake_tile_type), pointer :: lake_new => NULL()

     lake_new => new_lake_tile(1)
     lake_new%T=lake_T; lake_new%wl=lake_wl; lake_new%ws=lake_ws; lake_new%dz=lake_dz
     call lake_relayer (lake_new)
     lake_T=lake_new%T; lake_wl=lake_new%wl; lake_ws=lake_new%ws; lake_dz=lake_new%dz
     call delete_lake_tile(lake_new)

  end subroutine lake_relayer2

! ============================================================================
  subroutine lake_relayer_converge (lake_wl, lake_ws, lake_T, lake_dz)
     real, dimension(num_l), intent(inout) ::  lake_wl, lake_ws, lake_T, lake_dz

     type(lake_tile_type), pointer :: lake_new => NULL()
     integer :: n
     integer :: n_max = 100000
     real :: r = 0.

     if((lake_wl(2)+lake_ws(2))==0.) &
       call error_mesg('lake_relayer_converge', '(lake_wl(2)+lake_ws(2)) cannot be zero', FATAL)
     n = 0
     r = (lake_wl(1)+lake_ws(1))/(lake_wl(2)+lake_ws(2))
     do while((r.gt.max_rat).or.(r.lt.min_rat))
       lake_new => new_lake_tile(1)
       lake_new%T=lake_T; lake_new%wl=lake_wl; lake_new%ws=lake_ws; lake_new%dz=lake_dz
       call lake_relayer (lake_new)
       lake_T=lake_new%T; lake_wl=lake_new%wl; lake_ws=lake_new%ws; lake_dz=lake_new%dz
       call delete_lake_tile(lake_new)
       n = n+1
       if(n>=n_max) call error_mesg('lake_relayer_converge', 'relayer too many times', FATAL)
       if((lake_wl(2)+lake_ws(2))==0.) &
         call error_mesg('lake_relayer_converge', '(lake_wl(2)+lake_ws(2)) cannot be zero', FATAL)
       r = (lake_wl(1)+lake_ws(1))/(lake_wl(2)+lake_ws(2))
     enddo

  end subroutine lake_relayer_converge
! ============================================================================
  subroutine remove_negative_water(lake_wl, lake_ws, lake_T, lake_dz, lake_dhcap, rsv_zmin)
    real, dimension(num_l), intent(inout) ::  lake_wl, lake_ws, lake_T, lake_dz
    real, dimension(num_l), intent(in)    ::  lake_dhcap
    real, intent(in) :: rsv_zmin  !m

    integer :: n
    integer :: n_max = 100000

    call melt_negative(lake_wl, lake_ws, lake_T, lake_dz, lake_dhcap)
    n = 0
    do while(lake_ws(1)<0.or.lake_wl(1)<0.)
      !if(sum(lake_dz) <= rsv_zmin)then
      !  call error_mesg('remove_negative_water', 'water in tile is too small', NOTE)
      !  exit
      !endif
      call lake_relayer2 (lake_wl, lake_ws, lake_T, lake_dz)
      call melt_negative(lake_wl, lake_ws, lake_T, lake_dz, lake_dhcap)
      n = n+1
      if(n>=n_max) &
        call land_error_message('remove_negative_water :: relayer too many times', FATAL)
    enddo

  end subroutine remove_negative_water

! ============================================================================
  subroutine melt_negative(lake_wl,lake_ws,lake_T, lake_dz, lake_dhcap)
    real, dimension(num_l), intent(inout) ::  lake_wl, lake_ws, lake_T, lake_dz
    real, dimension(num_l), intent(in)    :: lake_dhcap ! dry heat capacity of lakes. Is it another "fictitious heat"?

    integer :: l
    real :: hcap, melt

    do l = 1, num_l
      if(lake_ws(l)<0.or.lake_wl(l)<0.)then
        hcap = lake_dhcap(l)*lake_dz(l) &
             + clw*lake_wl(l) + csw*lake_ws(l)
        if(lake_ws(l)<0..and.lake_wl(l)>0.)then
          melt = lake_ws(l)
        else if(lake_wl(l)<0..and.lake_ws(l)>0.)then
          melt = -lake_wl(l)
        else
          melt = 0.
        endif
        lake_wl(l) = lake_wl(l) + melt
        lake_ws(l) = lake_ws(l) - melt
        lake_T(l) = tfreeze &
                   + (hcap*(lake_T(l)-tfreeze) - hlf*melt) &
                            / ( hcap + (clw-csw)*melt )
        if(l>=2.and.(lake_ws(l)<0.or.lake_wl(l)<0.)) then
           call check_var_range(lake_wl(l),0.0,HUGE(1.0),'melt_negative','lake_wl',FATAL)
           call check_var_range(lake_ws(l),0.0,HUGE(1.0),'melt_negative','lake_ws',FATAL)
        endif
      endif !if(lake_ws(l)<0.or.lake_wl(l)<0.)then
    enddo

  end subroutine melt_negative

! ============================================================================
! conduct lake abstraction for irrigation
! if both use_reservoir and do_lake_abstraction are true, water can only be extracted from reservoir
! if use_reservoir is false and do_lake_abstraction is true, water is extracted from lake
! if use_reservoir is true and do_lake_abstraction is false, no water extraction is allowed, but reservoir will act with lake_abst==0.
subroutine lake_abstraction (is_terminal, &
                             irr_demand, Afrac_rsv, Vfrac_rsv, &
                             influx, influx_c, &
                             tot_area, lake_depth_sill, rsv_depth, env_flow, &
                             lake_T, lake_wl, lake_ws, lake_dz, lake_dhcap, &
                             lake_abst, lake_habst, &
                             rsv_out, rsv_out_s, rsv_out_h, vr1)

  logical, intent(in) :: is_terminal
  real, intent(inout) :: irr_demand !m3
  real, intent(in)    :: Afrac_rsv, Vfrac_rsv
  real, intent(in) :: influx !kg
  real, intent(in) :: influx_c(2) !kg, J
  real, intent(in)    :: tot_area !m2 lake_area
  real, intent(in)    :: lake_depth_sill, rsv_depth !m
  real, intent(in)    :: env_flow !m3
  real, dimension(num_l), intent(inout) ::  lake_T, lake_wl, lake_ws, lake_dz
  real, dimension(num_l), intent(in)    ::  lake_dhcap
  real, intent(out)   :: lake_abst !m3
  real, intent(out)   :: lake_habst !J
  real, intent(inout) :: rsv_out !kg
  real, intent(out)   :: rsv_out_s, rsv_out_h !kg, J
  real, intent(out)   :: vr1 !m3

  real :: res_capacity, lake_avail, lake_abst_vol, lake_abst_tot, &
          lake_this_lev, lake_collected, frac_abst
  real :: v0, vr0 !m3
  real :: v0_liq, vr0_liq !m3
  integer :: n
  integer :: n_max=1000
  real :: abst_thres = 1.e-15 !m3
  real :: lake_ws_thres = 1. !kg/m2, 0.001 m
  real :: r

  lake_abst = 0. !m3
  lake_habst = 0. !J

  !if we use reservoir, all water must be extracted from reservoir, otherwise, all water must be extracted from lake
  if(use_reservoir.and.Afrac_rsv<=0.)then
    rsv_out = 0. ;  vr1 = 0.
    rsv_out_s = 0. ; rsv_out_h = 0.
    return
  endif

  if(use_reservoir)then !We extract water from reservoir only, if there is no reservoir, we don't extract water.
    res_capacity = (Afrac_rsv*tot_area) * rsv_depth  !m3
    frac_abst = Vfrac_rsv
  else ! There is no reservoir at all, and we extract water from lake.
    res_capacity = tot_area * lake_depth_sill !m3
    frac_abst = 1.
  endif


  v0 = (sum(lake_wl+lake_ws)*tot_area - influx)/DENS_H2O !m3
  vr0 = frac_abst*v0 !m3   if v0<0, we must have vr0==0

  if(sum(lake_wl+lake_ws)*tot_area/DENS_H2O <= ResMin*res_capacity)then !m3
   rsv_out = 0. ;  vr1 = 0.
   rsv_out_s = 0. ; rsv_out_h = 0.
   if(use_reservoir) vr1 = vr0 + influx/DENS_H2O !m3
   return
  endif

  if(use_reservoir.or.do_lake_abstraction)then
    call remove_negative_water(lake_wl, lake_ws, lake_T, lake_dz, lake_dhcap, ResMin*res_capacity/tot_area)
    !if(sum(lake_dz)<=ResMin*res_capacity/tot_area)then !m
    !  rsv_out = 0. ;  vr1 = 0.
    !  rsv_out_s = 0. ; rsv_out_h = 0.
    !  if(use_reservoir) vr1 = vr0 + influx/DENS_H2O
    !  return
    !endif
  endif

  v0_liq = (sum(lake_wl)*tot_area - (influx-influx_c(1)))/DENS_H2O !m3
  vr0_liq = frac_abst*v0_liq !m3
  lake_avail = min(vr0_liq + (influx-influx_c(1))/DENS_H2O, & !all liquid water in reservoir
                   vr0 + influx/DENS_H2O - ResMin*res_capacity)   !m3
  lake_avail = min(sum(lake_wl)*tot_area/DENS_H2O, lake_avail) !m3
  lake_avail = min(sum(lake_wl+lake_ws)*tot_area/DENS_H2O-ResMin*res_capacity, lake_avail) !m3
  lake_avail = max(0., lake_avail)


  IF (do_lake_abstraction.and.&
      (lake_ws(1) < lake_ws_thres).and.&
      (irr_demand > abst_thres).and.&
      (lake_avail > 0.)) then

     lake_abst_vol = min(irr_demand, lake_avail) !m3
     lake_abst_tot = lake_abst_vol*DENS_H2O !kg
     lake_collected = 0.
     n=0
     do while( lake_collected<lake_abst_tot-abst_thres*DENS_H2O .and. n<=n_max )
       lake_this_lev = max(0.,min((lake_abst_tot-lake_collected), (tot_area*lake_wl(1)))) !kg
       lake_wl(1) = max(0., lake_wl(1)-lake_this_lev/tot_area) !kg/m2
       lake_habst = lake_habst + clw*(lake_T(1)-tfreeze)*lake_this_lev !J/(kg K) * K * kg = J
       lake_collected = lake_collected + lake_this_lev !kg
       r = (lake_wl(1)+lake_ws(1))/(lake_wl(2)+lake_ws(2))
       if(lake_wl(1)==0..and.r<min_rat)then
         call lake_relayer2(lake_wl, lake_ws, lake_T, lake_dz)
       else
         exit
       endif
       n=n+1
     enddo
     if(n>=n_max) call error_mesg('lake_abstraction', 'relayer too many times', NOTE)
     irr_demand = max(0., irr_demand - lake_collected/DENS_H2O) !kg / kg/m3 = m3
     lake_abst = lake_collected/DENS_H2O !kg / kg/m3 = m3

     if(sum(lake_wl+lake_ws)*tot_area/DENS_H2O < (ResMin-0.01)*res_capacity) &
       call error_mesg('lake_abstraction', 'water in tile is less than ResMin*res_capacity', FATAL)
     if(vr0 + influx/DENS_H2O - lake_abst < (ResMin-0.01)*res_capacity) &
       call error_mesg('lake_abstraction', 'water in reservoir is less than ResMin*res_capacity', FATAL)
  ENDIF !  if (do_lake_abstraction) then

 !calculate reservoir outflow
 if(use_reservoir)then
   if(vr0 + influx/DENS_H2O - lake_abst <= ResMax*res_capacity)then !m3
     rsv_out = 0. !m3
   else
     rsv_out = vr0 + influx/DENS_H2O - lake_abst - ResMax*res_capacity !m3
   endif
   rsv_out = max(env_flow, rsv_out) !m3
   rsv_out = min(vr0 + influx/DENS_H2O - lake_abst - ResMin*res_capacity, rsv_out) !m3
   rsv_out = max(0., rsv_out) !m3
   !if here is river terminal point, and all area is reservoir, then no rsv_out allowed.
   if(is_terminal.and.Afrac_rsv>=1.) rsv_out = 0. !m3
   vr1 = vr0 + influx/DENS_H2O - lake_abst - rsv_out !m3
   !if(vr1<0.) &
     !call error_mesg('lake_abstraction', 'vr1 could not be less than 0', FATAL)
   rsv_out = rsv_out*DENS_H2O !m3 * kg/m3 = kg
   rsv_out_s = 0. ; rsv_out_h = 0.
   if(rsv_out>0..and.Afrac_rsv>=1) then !special case: no lake, only reservoir
     call rsv_outflow_c(lake_wl,lake_ws,lake_T,lake_dz,&
                        tot_area,rsv_out,&
                        rsv_out_s,rsv_out_h) !we need to know rsv_out_s,rsv_out_h only when there is no lake
     vr1 = vr0 + influx/DENS_H2O - lake_abst - rsv_out/DENS_H2O !m3
     !if(vr1<0.) &
       !call error_mesg('lake_abstraction', 'vr1 could not be less than 0', FATAL)
   endif
 else
   rsv_out = 0. ;  vr1 = 0.
   rsv_out_s = 0. ; rsv_out_h = 0.
 endif

end subroutine lake_abstraction

! ============================================================================
subroutine prohibit_shallow_lake(lake)
  type(lake_tile_type), intent(inout) :: lake

  !real :: fac_min = 0.1
  real :: heat0, heat1, lm0, lm1, fm0, fm1
  real :: liq_frac
  real, dimension(num_l) :: dz_frac
  real :: new_z
  integer :: l

  if(.not.prohibit_shallowlake) return
  if(sum(lake%dz)>=ResMin*lake%pars%depth_sill) return

  heat0 = lake_tile_heat(lake)
  call lake_tile_stock_pe(lake, lm0, fm0)

  dz_frac = lake%dz/sum(lake%dz)
  new_z = ResMin*lake%pars%depth_sill
  lake%dz = dz_frac*new_z

  do l = 2, num_l
    liq_frac = lake%wl(l)/(lake%wl(l)+lake%ws(l))
    lake%wl(l)=lake%dz(l)*DENS_H2O*liq_frac
    lake%ws(l)=lake%dz(l)*DENS_H2O*(1.-liq_frac)
  enddo

  heat1 = lake_tile_heat(lake) !J/m2
  call lake_tile_stock_pe(lake, lm1, fm1)

  lake%sub_lmass = lake%sub_lmass - (lm1 - lm0)
  lake%sub_fmass = lake%sub_fmass - (fm1 - fm0)
  lake%sub_heat  = lake%sub_heat  - (heat1 - heat0)


end subroutine prohibit_shallow_lake

! ============================================================================

subroutine rsv_outflow_c(lake_wl,lake_ws,lake_T,lake_dz,&
                         tot_area,qt,&
                         qs,qh)

  real, dimension(num_l), intent(inout) ::  lake_T, lake_wl, lake_ws, lake_dz
  real, intent(in)    :: tot_area !m2 lake_area
  real, intent(inout) :: qt !kg
  real, intent(out)   :: qs !kg
  real, intent(out)   :: qh !J

  integer   :: n
  real      :: ql
  real      :: out_frac = 0.
  real      :: qt_to_flow, qt_this_lev, liq_this_lev, ice_this_lev
  integer   :: n_max=1000

  qt_to_flow = qt !kg
  ql = 0.; qs = 0.; qh = 0.; n = 0
  do while( qt_to_flow>0. .and. n<=n_max )
     if(lake_wl(1)>=0..and.lake_ws(1)>=0.and.(lake_wl(1)+lake_ws(1))>0.)then
       out_frac = lake_wl(1)/(lake_wl(1)+lake_ws(1))
     else
       exit
     endif
     qt_this_lev = max(0.,min(qt_to_flow, (lake_wl(1)+lake_ws(1))*tot_area)) !kg
     liq_this_lev = out_frac*qt_this_lev !kg
     ice_this_lev = (1.-out_frac)*qt_this_lev  !kg
     lake_wl(1) = max(0., lake_wl(1)-liq_this_lev/tot_area) !kg/m2
     lake_ws(1) = max(0., lake_ws(1)-ice_this_lev/tot_area) !kg/m2
     qt_to_flow = qt_to_flow - qt_this_lev !kg
     qh = qh + (clw*liq_this_lev+csw*ice_this_lev)*(lake_T(1)-tfreeze) !J
     ql = ql + liq_this_lev
     qs = qs + ice_this_lev
     if(lake_wl(1)==0.) call lake_relayer2(lake_wl, lake_ws, lake_T, lake_dz)
     n=n+1
  enddo
  if(n>=n_max) call error_mesg('rsv_outflow_c', 'relayer too many times', NOTE)

  qt = ql + qs !kg

end subroutine rsv_outflow_c

! ============================================================================
subroutine lake_diag_init(id_ug)
  integer,intent(in) :: id_ug !<Unstructured axis id.

  ! ---- local vars
  integer :: axes(2)
  integer :: id_zhalf, id_zfull

  ! define vertical axis
  id_zhalf = diag_axis_init ( &
       'zhalf_lake', zhalf(1:num_l+1), 'meters', 'z', 'half level',  -1, set_name='lake' )
  id_zfull = diag_axis_init ( &
       'zfull_lake', zfull(1:num_l),   'meters', 'z', 'full level',  -1, set_name='lake', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/id_ug,id_zfull/)

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('lake')

  ! define static diagnostic fields
  id_sillw = register_tiled_static_field ( module_name, 'lake_width', &
       axes(1:1), 'lake width at outflow', 'm', missing_value=-100.0 )
  id_silld = register_tiled_static_field ( module_name, 'lake_depth', &
       axes(1:1), 'lake depth below sill', 'm', missing_value=-100.0 )
  id_backw = register_tiled_static_field ( module_name, 'backwater', &
       axes(1:1), 'backwater flag', '-', missing_value=-100.0 )
  id_back1 = register_tiled_static_field ( module_name, 'backwater_1', &
       axes(1:1), 'backwater1 flag', '-', missing_value=-100.0 )

  ! define dynamic diagnostic fields
  id_dz  = register_tiled_diag_field ( module_name, 'lake_dz', axes,         &
       lnd%time, 'nominal layer thickness', 'm', missing_value=-100.0 )
  id_wl  = register_tiled_diag_field ( module_name, 'lake_wl', axes,         &
       lnd%time, 'liquid water mass', 'kg/m2', missing_value=-100.0 )
  id_ws  = register_tiled_diag_field ( module_name, 'lake_ws', axes,         &
       lnd%time, 'solid water mass', 'kg/m2', missing_value=-100.0 )
  id_lwc  = register_tiled_diag_field ( module_name, 'lake_liq',  axes,       &
       lnd%time, 'bulk density of liquid water', 'kg/m3',  missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'lake_ice',  axes,       &
       lnd%time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'lake_T',  axes,        &
       lnd%time, 'temperature',            'degK',  missing_value=-100.0 )
  id_K_z  = register_tiled_diag_field ( module_name, 'lake_K_z', axes,         &
       lnd%time, 'vertical diffusivity', 'm2/s', missing_value=-100.0 )
  id_evap  = register_tiled_diag_field ( module_name, 'lake_evap',  axes(1:1),  &
       lnd%time, 'lake evap',            'kg/(m2 s)',  missing_value=-100.0 )

  id_Afrac_rsv = register_tiled_diag_field ( module_name, 'Afrac_rsv',  axes(1:1),  &
       lnd%time, 'reservoir area fraction to the tile', 'none',  missing_value=-100.0 )
  id_Vfrac_rsv = register_tiled_diag_field ( module_name, 'Vfrac_rsv',  axes(1:1),  &
       lnd%time, 'reservoir area fraction to the tile', 'none',  missing_value=-100.0 )
  id_rsv_depth  = register_tiled_diag_field ( module_name, 'rsv_depth',  axes(1:1),  &
       lnd%time, 'reservoir construction depth',            'm',  missing_value=-100.0 )

  id_sub_lmass = register_tiled_diag_field ( module_name, 'sub_lmass',  axes(1:1),  &
       lnd%time, 'buried liquid water under lake due to reservoir building', 'kg/m2',  missing_value=-100.0 )
  id_sub_fmass = register_tiled_diag_field ( module_name, 'sub_fmass',  axes(1:1),  &
       lnd%time, 'buried frozen water under lake due to reservoir building', 'kg/m2',  missing_value=-100.0 )
  id_sub_heat  = register_tiled_diag_field ( module_name, 'sub_heat',   axes(1:1),  &
       lnd%time, 'buried heat under lake due to reservoir building', 'J/m2',  missing_value=-100.0 )
  id_sub_cmass = register_tiled_diag_field ( module_name, 'sub_cmass',  axes(1:1),  &
       lnd%time, 'buried carbon under lake due to reservoir building', 'kgC/m2',  missing_value=-100.0 )

  id_lake_area = register_tiled_diag_field ( module_name, 'lake_area', axes(1:1), &
       lnd%time, 'lake area', 'm2',  missing_value=-100.0 )
  id_lake_frac = register_tiled_diag_field ( module_name, 'lake_frac', axes(1:1), &
       lnd%time, 'lake frac', '-',  missing_value=-100.0 )

  call add_tiled_static_field_alias (id_silld, module_name, 'sill_depth', &
       axes(1:1), 'obsolete, pls use lake_depth (static)','m', &
       missing_value=-100.0 )
  call add_tiled_static_field_alias (id_sillw, module_name, 'sill_width', &
       axes(1:1), 'obsolete, pls use lake_width (static)','m', &
       missing_value=-100.0 )

end subroutine lake_diag_init

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function lake_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   lake_tile_exists = associated(tile%lake)
end function lake_tile_exists


! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine lake_dz_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr => tile%lake%dz(i)
   endif
end subroutine lake_dz_ptr

subroutine lake_temp_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr => tile%lake%T(i)
   endif
end subroutine lake_temp_ptr

subroutine lake_wl_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr => tile%lake%wl(i)
   endif
end subroutine lake_wl_ptr

subroutine lake_ws_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr => tile%lake%ws(i)
   endif
end subroutine lake_ws_ptr

subroutine lake_connected_to_next_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%connected_to_next
   endif
end subroutine lake_connected_to_next_ptr

subroutine lake_depth_sill_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%depth_sill
   endif
end subroutine lake_depth_sill_ptr

subroutine lake_whole_area_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%whole_area
   endif
end subroutine lake_whole_area_ptr

subroutine lake_width_sill_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%width_sill
   endif
end subroutine lake_width_sill_ptr

subroutine lake_backwater_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%backwater
   endif
end subroutine lake_backwater_ptr

subroutine lake_backwater_1_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%pars%backwater_1
   endif
end subroutine lake_backwater_1_ptr

subroutine lake_Afrac_rsv_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%Afrac_rsv
   endif
end subroutine lake_Afrac_rsv_ptr

subroutine lake_Vfrac_rsv_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%Vfrac_rsv
   endif
end subroutine lake_Vfrac_rsv_ptr

subroutine lake_depth_rsv_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%rsv_depth
   endif
end subroutine lake_depth_rsv_ptr

subroutine lake_sub_lmass_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%sub_lmass
   endif
end subroutine lake_sub_lmass_ptr

subroutine lake_sub_fmass_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%sub_fmass
   endif
end subroutine lake_sub_fmass_ptr

subroutine lake_sub_heat_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%sub_heat
   endif
end subroutine lake_sub_heat_ptr

subroutine lake_sub_cmass_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%sub_cmass
   endif
end subroutine lake_sub_cmass_ptr


end module lake_mod



