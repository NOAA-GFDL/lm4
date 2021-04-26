! ============================================================================
! canopy air
! ============================================================================
module canopy_air_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : error_mesg, FATAL, NOTE, file_exist, &
     close_file, check_nml_error, mpp_pe, mpp_root_pe, stdlog, string, lowercase
use constants_mod, only : VONKARM, dens_h2o, pi, grav
use field_manager_mod, only : parse, MODEL_ATMOS, MODEL_LAND
use tracer_manager_mod, only : get_tracer_index, get_tracer_names, &
     query_method, NO_TRACER

use land_constants_mod, only : mol_CO2, mol_air, diffusivity_h2o, kin_visc_air, thermal_diff_air
use land_numerics_mod, only : gamma, gammaln
use land_tracers_mod, only : ntcana, isphum, ico2
use land_debug_mod, only : is_watch_point
use cana_tile_mod, only : cana_tile_type, &
     canopy_air_mass, canopy_air_mass_for_tracers, cpw
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_data_mod, only : lnd, log_version
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, get_tile_data, field_exists
use land_tile_diag_mod, only : register_tiled_diag_field, &
     send_tile_data, diag_buff_type, set_default_diag_filter
use vegn_tile_mod, only: vegn_tile_type, vegn_tile_bwood, vegn_tile_LAI, vegn_tile_SAI
use soil_tile_mod, only: soil_tile_type, soil_data_hydraulic_properties, dz, &
        get_soil_litter_C

implicit none
private

! ==== public interfaces =====================================================
public :: read_cana_namelist
public :: cana_init
public :: cana_end
public :: save_cana_restart
public :: cana_turbulence
public :: cana_roughness
public :: surface_resistances

! there are other public data members, see/search below
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name   = 'canopy_air_mod'
character(len=*), parameter :: diag_mod_name = 'cana'
#include "../shared/version_variable.inc"

! options for turbulence parameter calculations
integer, parameter :: TURB_LM3W = 1, TURB_LM3V = 2

! options of soil surface resistance calculations
integer, parameter :: &
   RESIST_NONE   = 0, & ! no extra soil resistance
   RESIST_HO2013 = 1    ! soil resistance based on Haghighi and Or (2013) and related papers

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
real :: init_T           = 288.
real :: init_T_cold      = 260.
real :: init_q           = 0.
real :: init_co2         = 350.0e-6 ! ppmv = mol co2/mol of dry air
character(len=32) :: turbulence_to_use = 'lm3w' ! or lm3v
logical :: use_SAI_for_heat_exchange = .FALSE. ! if true, con_v_h is calculated for LAI+SAI
   ! traditional treatment (default) is to only use SAI
logical :: save_qco2     = .TRUE.
! resistance-related namelist variables
character(32) :: soil_resistance_to_use = 'none' ! or 'HO2013'
real :: bare_rah_sca      = 0.01 ! bare-ground resistance between ground and canopy air, s/m
  ! resistances in soil upper layer and viscous sublayer
real :: rav_lit_0         = 0.0 ! constant litter resistance to vapor
real :: rav_lit_vi        = 0.0 ! litter resistance to vapor per LAI+SAI
real :: rav_lit_fsc       = 0.0 ! litter resistance to vapor per fsc
real :: rav_lit_ssc       = 0.0 ! litter resistance to vapor per ssc
real :: rav_lit_deadmic   = 0.0 ! litter resistance to vapor per dead microbe C
real :: rav_lit_bwood     = 0.0 ! litter resistance to vapor per bwood
real :: d_visc_max        = 0.1 ! when positive, max thickness of viscous sublayer (m);
                                ! negative or zero turn off limitation
real :: k_over_B          = 2.0 !
! fog-related namelist variables
logical, protected, public :: do_fog   = .FALSE.  ! if true, water vapore in canopy air can form condensate
real,    protected, public :: fog_form_rate = 1.0
real,    protected, public :: fog_diss_time = 0.0   ! e-folding time of fog evaporation in dry environment, s
                                                    ! 0 or below means instant dissipation

namelist /cana_nml/ &
  init_T, init_T_cold, init_q, init_co2, turbulence_to_use, use_SAI_for_heat_exchange, &
  canopy_air_mass, canopy_air_mass_for_tracers, cpw, save_qco2, bare_rah_sca, &
  k_over_B, &
  ! soil resistance parameters
  soil_resistance_to_use, &
  d_visc_max, &
  rav_lit_0, rav_lit_vi, rav_lit_fsc, rav_lit_ssc, rav_lit_deadmic, rav_lit_bwood, &
  ! fog-related namelists
  do_fog, fog_form_rate, fog_diss_time

!---- end of namelist --------------------------------------------------------

logical :: module_is_initialized =.FALSE.
integer :: turbulence_option ! selected option of turbulence parameters
     ! calculations
integer :: soil_resistance_option = -1 ! option of soil resistance parameterization

! ---- diag field IDs
integer :: id_r_litt_evap, id_r_bl_sens, id_r_bl_evap, id_r_sv_evap, &
           id_c_litt_evap, id_c_bl_sens, id_c_bl_evap, id_c_sv_evap, &
           id_d_visc, id_ustar_sfc, id_u_sfc, id_theta_sfc

contains

! ============================================================================
subroutine read_cana_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=cana_nml, iostat=io)
     ierr = check_nml_error(io, 'cana_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=cana_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'cana_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=cana_nml)
  endif

  ! initialize options, to avoid expensive string comparisons during
  ! run-time
  if (trim(turbulence_to_use)=='lm3v') then
     turbulence_option = TURB_LM3V
  else if (trim(turbulence_to_use)=='lm3w') then
     turbulence_option = TURB_LM3W
  else
     call error_mesg('cana_init', 'canopy air turbulence option turbulence_to_use="'// &
          trim(turbulence_to_use)//'" is invalid, use "lm3w" or "lm3v"', FATAL)
  endif

  ! convert symbolic names of surface resistance options into numeric IDs to
  ! avoid expensive string comparisons run-time
  if (trim(lowercase(soil_resistance_to_use))=='none') then
     soil_resistance_option = RESIST_NONE
  else if (trim(lowercase(soil_resistance_to_use))=='ho2013') then
     soil_resistance_option = RESIST_HO2013
  else
     call error_mesg('surface_resistance_init',&
          'soil resistance option soil_resistance_to_use="'//&
          trim(soil_resistance_to_use)//'" is invalid, use "none" or "HO"',&
          FATAL)
  endif
end subroutine read_cana_namelist

! ============================================================================
! initialize canopy air
subroutine cana_init (id_ug)
  integer, intent(in) :: id_ug   !<Unstructured axis id.

  ! ---- local vars ----------------------------------------------------------
  type(land_tile_enum_type)     :: ce ! last and current tile
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  character(*), parameter :: restart_file_name='INPUT/cana.res.nc'
  type(land_restart_type) :: restart
  logical :: restart_exists

  character(32)  :: name  ! name of the tracer
  integer        :: tr, i ! tracer indices
  real           :: init_tr(ntcana) ! initial (cold-start) values of tracers
  real           :: value ! used for parameter parsing
  character(32)  :: scheme
  character(1024) :: parameters

  ! ---- initialize cana state -----------------------------------------------
  ! get the initial conditions for tracers

  ! For  now, we get the initial value from the surface_value parameter of the
  ! *atmospheric* tracer vertical profile, except for co2 and sphum. init_q in
  ! the cana_nml sets the initial value of the specific humidity, and init_co2
  ! sets the initial value of the dry volumetric mixing ration for co2.
  ! If surface_value is not defined in tracer table, then initial condition is zero
  init_tr(:) = 0.0
  do i = 1, ntcana
     call get_tracer_names(MODEL_LAND, i, name=name)
     tr = get_tracer_index(MODEL_ATMOS, name)
     if (tr==NO_TRACER) cycle ! no such tracer in the atmos
     ! TODO: possibly we need to add an ability to read init value from some parameter
     ! in the land tracer table?
     scheme = ''; parameters = ''
     if (query_method('profile_type', MODEL_ATMOS, tr, scheme, parameters)) then
        if (parse(parameters,'surface_value',value)>0) init_tr(i) = value
     endif
  enddo
  init_tr(isphum) = init_q
  init_tr(ico2)   = init_co2*mol_CO2/mol_air*(1-init_tr(isphum)) ! convert to kg CO2/kg wet air

  ! first, set the initial values
  ce = first_elmt(land_tile_map)
  do while(loop_over_tiles(ce, tile))
     if (.not.associated(tile%cana)) cycle

     if (associated(tile%glac)) then
        tile%cana%T = init_T_cold
     else
        tile%cana%T = init_T
     endif
     tile%cana%fog = 0.0
     tile%cana%tr(:) = init_tr(:)
  enddo

  ! then read the restart if it exists
  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then
     call error_mesg('cana_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     call get_tile_data(restart, 'temp', cana_T_ptr)
     if (field_exists(restart,'fog')) &
        call get_tile_data(restart, 'fog', cana_fog_ptr)
     do tr = 1, ntcana
        call get_tracer_names(MODEL_LAND, tr, name=name)
        if (field_exists(restart,trim(name))) then
           call error_mesg('cana_init','reading tracer "'//trim(name)//'"',NOTE)
           call get_tile_data(restart,name,cana_tr_ptr,tr)
        else
           call error_mesg('cana_init', 'tracer "'//trim(name)// &
                '" was set to initial value '//string(init_tr(tr)), NOTE)
        endif
     enddo
  else
     call error_mesg('cana_init',&
          'cold-starting canopy air',&
          NOTE)
  endif
  call free_land_restart(restart)

  module_is_initialized = .TRUE.

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('land')
  id_r_bl_sens = register_tiled_diag_field( diag_mod_name, 'r_bl_sens', &
       (/id_ug/), lnd%time, 'resistance of viscous boundary layer to heat flux', 's/m', missing_value=-1.0 )
  id_r_bl_evap = register_tiled_diag_field( diag_mod_name, 'r_bl_evap', &
       (/id_ug/), lnd%time, 'resistance of viscous boundary layer to water vapor flux', 's/m', missing_value=-1.0 )

  id_c_bl_sens = register_tiled_diag_field( diag_mod_name, 'c_bl_sens', &
       (/id_ug/), lnd%time, 'conductance of viscous boundary layer for heat', 'm/s', missing_value=-1.0 )
  id_c_bl_evap = register_tiled_diag_field( diag_mod_name, 'c_bl_evap', &
       (/id_ug/), lnd%time, 'conductance of viscous boundary layer for water vapor', 'm/s', missing_value=-1.0 )

  id_d_visc = register_tiled_diag_field( diag_mod_name, 'd_visc', &
       (/id_ug/), lnd%time, 'thickness of viscous sublayer', 'm', missing_value=-1.0 )
  id_u_sfc = register_tiled_diag_field( diag_mod_name, 'u_sfc', &
       (/id_ug/), lnd%time, 'near-surface wind velocity', 'm/s', missing_value=-1.0 )
  id_ustar_sfc = register_tiled_diag_field( diag_mod_name, 'ustar_sfc', &
       (/id_ug/), lnd%time, 'friction velocity at the surface', 'm/s', missing_value=-1.0 )

  call set_default_diag_filter('soil')
  id_r_litt_evap = register_tiled_diag_field( diag_mod_name, 'r_litt_evap', &
       (/id_ug/), lnd%time, 'resistance of litter layer to water vapor flux', 's/m', missing_value=-1.0 )
  id_c_litt_evap = register_tiled_diag_field( diag_mod_name, 'c_litt_evap', &
       (/id_ug/), lnd%time, 'conductance of litter layer for water vapor', 'm/s', missing_value=-1.0 )
  id_r_sv_evap = register_tiled_diag_field( diag_mod_name, 'r_sv_evap', &
       (/id_ug/), lnd%time, 'resistance of near-surface soil to for water flux', 's/m', missing_value=-1.0 )
  id_c_sv_evap = register_tiled_diag_field( diag_mod_name, 'c_sv_evap', &
       (/id_ug/), lnd%time, 'conductance of near-surface soil for water flux', 'm/s', missing_value=-1.0 )
  id_theta_sfc = register_tiled_diag_field( diag_mod_name, 'theta_sfc', &
       (/id_ug/), lnd%time, 'relative soil wetness', 'unitless', missing_value=-1.0 )
end subroutine cana_init


! ============================================================================
! release memory
subroutine cana_end ()
  module_is_initialized =.FALSE.
end subroutine cana_end


! ============================================================================
subroutine save_cana_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object
  character(len=32)  :: name,units
  character(len=128) :: longname
  integer :: tr

  call error_mesg('cana_end','writing NetCDF restart',NOTE)
! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp)//'cana.res.nc'
  call init_land_restart(restart, filename, cana_tile_exists, tile_dim_length)

  ! write temperature
  call add_tile_data(restart,'temp',cana_T_ptr,'canopy air temperature','degrees_K')
  call add_tile_data(restart,'fog',cana_fog_ptr,'canopy air condensate mass','kg/m2')
  do tr = 1,ntcana
     call get_tracer_names(MODEL_LAND, tr, name, longname, units)
     if (tr==ico2.and..not.save_qco2) cycle
     call add_tile_data(restart,name,cana_tr_ptr,tr,'canopy air '//trim(longname),trim(units))
  enddo
  call save_land_restart(restart)
  call free_land_restart(restart)
end subroutine save_cana_restart

! ============================================================================
subroutine cana_turbulence (u_star, &
     vegn_cover, vegn_layerfrac, vegn_height, vegn_bottom, vegn_lai, vegn_sai, vegn_d_leaf, &
     land_d, land_z0m, land_z0s, grnd_z0s, &
     con_v_h, con_v_v, con_g_h, con_g_v, u_sfc, ustar_sfc )
  real, intent(in) :: &
       u_star, & ! friction velocity, m/s
       land_d, land_z0m, land_z0s, grnd_z0s, &
       vegn_cover, vegn_height(:), vegn_layerfrac(:), &
       vegn_bottom(:), & ! height of the bottom of the canopy, m
       vegn_lai(:), vegn_sai(:), vegn_d_leaf(:)
  real, intent(out) :: &
       con_v_h(:), con_v_v(:), & ! one-sided foliage-CAS conductance per unit ground area
       con_g_h   , con_g_v,    & ! ground-CAS turbulent conductance per unit ground area
       u_sfc, &                  ! near-surface wind speed, m/s
       ustar_sfc                 ! near-surface friction velocity, m/s

  !---- local constants
  real, parameter :: a_max = 3
  real, parameter :: leaf_co = 0.01 ! meters per second^(1/2)
                                    ! leaf_co = g_b(z)/sqrt(wind(z)/d_leaf)
  real, parameter :: min_thickness = 0.01 ! thickness for switching to thin-canopy approximation, m
  real, parameter :: min_height = 0.1 ! min height of the canopy in TURB_LM3V case, m
  ! ---- local vars
  real :: a        ! parameter of exponential wind profile within canopy:
                   ! u = u(ztop)*exp(-a*(1-z/ztop))
  real :: height   ! height of the current vegetation cohort, m
  real :: ztop     ! height of the tallest vegetation, m
  real :: wind     ! normalized wind on top of canopy, m/s
  real :: Kh_top   ! turbulent exchange coefficient on top of the canopy
  real :: vegn_idx ! total vegetation index = LAI+SAI, sum over cohorts
  real :: rah_sca  ! ground-SCA resistance
  real :: h0       ! height of the canopy bottom, m
  real :: gb       ! aerodynamic resistance per unit leaf (or stem) area

  integer :: i

  ! TODO: check array sizes

  vegn_idx = sum((vegn_lai+vegn_sai)*vegn_layerfrac)  ! total vegetation index

  select case(turbulence_option)
  case(TURB_LM3W)
     a  = max(vegn_cover,0.0)*a_max
     if(vegn_cover > 0) then
        ztop   = maxval(vegn_height(:))

        wind  = u_star/VONKARM*log((ztop-land_d)/land_z0m) ! normalized wind on top of the canopy
        do i = 1,size(vegn_lai)
           height = vegn_height(i) ! effective height of the vegetation
           h0     = vegn_bottom(i) ! height of the bottom of the canopy
           if(height-h0>min_thickness) then
              con_v_h(i) = 2*vegn_lai(i)*leaf_co*sqrt(wind/vegn_d_leaf(i))*ztop/(height-h0)&
                 *(exp(-a/2*(ztop-height)/ztop)-exp(-a/2*(ztop-h0)/ztop))/a
           else
              ! thin cohort canopy limit
              con_v_h(i) = vegn_lai(i)*leaf_co*sqrt(wind/vegn_d_leaf(i))&
                 *exp(-a/2*(ztop-height)/ztop)
           endif
        enddo
        con_g_h = u_star*a*VONKARM*(1-land_d/ztop) &
             / (exp(a*(1-grnd_z0s/ztop)) - exp(a*(1-(land_z0s+land_d)/ztop)))
     else
        con_v_h = 0
        con_g_h = 0
     endif
     con_v_v = con_v_h
  case(TURB_LM3V)
     ztop = max(maxval(vegn_height(:)),min_height)

     a = a_max
     wind=u_star/VONKARM*log((ztop-land_d)/land_z0m) ! normalized wind on top of the canopy

     do i = 1,size(vegn_lai)
        height = max(vegn_height(i),min_height) ! effective height of the vegetation
        h0     = vegn_bottom(i) ! height of the canopy bottom above ground
        if(height-h0>min_thickness) then
           gb = 2*leaf_co*sqrt(wind/vegn_d_leaf(i))*ztop/(height-h0)&
              *(exp(-a/2*(ztop-height)/ztop)-exp(-a/2*(ztop-h0)/ztop))/a
        else
           ! thin cohort canopy limit
           gb = leaf_co*sqrt(wind/vegn_d_leaf(i))&
              *exp(-a/2*(ztop-height)/ztop)
        endif
        con_v_v(i) = vegn_lai(i)*gb
        if (use_SAI_for_heat_exchange) then
           con_v_h(i) = (vegn_lai(i)+vegn_sai(i))*gb
        else
           con_v_h(i) = vegn_lai(i)*gb
        endif
     enddo

     if (land_d > 0.06 .and. vegn_idx > 0.25) then
        Kh_top = VONKARM*u_star*(ztop-land_d)
        rah_sca = ztop/a/Kh_top * &
             (exp(a*(1-grnd_z0s/ztop)) - exp(a*(1-(land_z0m+land_d)/ztop)))
        rah_sca = min(rah_sca,1250.0)
     else
        rah_sca = bare_rah_sca
     endif
     con_g_h = 1.0/rah_sca
  end select
  con_g_v = con_g_h

! u_sfc     = wind * exp(-a)
! ustar_sfc = u_star * exp(-a)
  u_sfc     = wind
  ustar_sfc = u_star/sqrt(2*vegn_idx + 1)

  if (is_watch_point()) then
     __DEBUG3__(vegn_idx,land_d,Kh_top)
     __DEBUG4__(ztop,u_star,wind,rah_sca)
  endif
end subroutine cana_turbulence

! ============================================================================
! update effective surface roughness lengths for CAS-to-atmosphere fluxes
! and conductances for canopy-to-CAS and ground-to-CAS fluxes
!
! Strategy: Always define a canopy present. Non-vegetated situation is simply
! a limit as vegetation density approaches (but isn't allowed to reach) zero.
! Create expressions for the outputs that reduce to the special
! cases of full canopy cover and no canopy. Full canopy solution is that
! from Bonan (NCAR/TN-417+STR, 1996, p. 63). Thus, setting cover=1 in
! recovers Bonan. Letting cover approach 0 makes con_v_coef go to zero,
! preventing exchange with canopy, and makes con_g_coef go infinite,
! removing sub-canopy resistance and putting all resistance above the
! canopy, where it can be affected by stability adjustments.
!
! ** However, there is still a problem with this formulation when a
! canopy is present, because surface flux (I think) is not told to
! subtract out the resistances associated with con_v_coef and con_g_coef,
! which thus seem to be double-counted. For testing LM2, we should set them
! to zero anyway.
subroutine cana_roughness(lm2, &
     subs_z0m, subs_z0s, &
     snow_z0m, snow_z0s, snow_area, &
     vegn_cover, vegn_height, vegn_lai, vegn_sai, &
     land_d, land_z0m, land_z0s, grnd_z0m, grnd_z0s )
  logical, intent(in) :: lm2
  real, intent(in) :: &
       subs_z0m, subs_z0s, snow_z0m, snow_z0s, snow_area, vegn_cover, vegn_height, &
       vegn_lai, vegn_sai
  real, intent(out) :: &
       land_d    ,&
       land_z0m  ,&
       land_z0s  ,&
       grnd_z0m  ,&
       grnd_z0s

  !---- local constants
  real, parameter :: d_h_max = 2./3.
  real, parameter :: z0m_h_max = 1/7.35

  ! ---- local vars
  real :: d_h      ! ratio of displacement height to vegetation height
  real :: z0m_h    ! ratio of roughness length to vegetation height
  real :: z0s_h, z0s_h_max
  real :: vegn_idx ! total vegetation index = LAI+SAI
  real :: height   ! effective vegetation height

  grnd_z0m = exp( (1-snow_area)*log(subs_z0m) + snow_area*log(snow_z0m))
  grnd_z0s = exp( (1-snow_area)*log(subs_z0s) + snow_area*log(snow_z0s))

  select case(turbulence_option)
  case(TURB_LM3W)
     if(vegn_cover > 0) then
        z0s_h_max = z0m_h_max*grnd_z0s/grnd_z0m ! to ensure cover->0 limit works
        d_h = vegn_cover*d_h_max
        if (lm2) then
           if (vegn_lai.gt.1) then  ! TEMP ***
              z0m_h = z0m_h_max
              z0s_h = z0s_h_max
           else
              z0m_h = grnd_z0m/vegn_height
              z0s_h = grnd_z0s/vegn_height
           endif
        else
           z0m_h = exp( vegn_cover*log(z0m_h_max) + (1-vegn_cover)*log(grnd_z0m/vegn_height))
           z0s_h = exp( vegn_cover*log(z0s_h_max) + (1-vegn_cover)*log(grnd_z0s/vegn_height))
        endif
        land_d   = d_h*vegn_height
        land_z0m = z0m_h*vegn_height
        land_z0s = z0s_h*vegn_height
     else
        land_d   = 0
        land_z0m = grnd_z0m
        land_z0s = grnd_z0s
     endif

  case(TURB_LM3V)
     height = max(vegn_height,0.1) ! effective height of the vegetation
     vegn_idx = vegn_lai+vegn_sai  ! total vegetation index
     if(vegn_idx>1e-4) then
        land_d = 1.1*height*log(1+(0.07*vegn_idx)**0.25)
        if(vegn_idx>2.85) then
           land_z0m = 0.3*(height-land_d)
        else
           land_z0m = grnd_z0m + 0.3*height*sqrt(0.07*vegn_idx)
        endif
     else
        land_d   = 0
        land_z0m = grnd_z0m
     endif
     land_z0s = land_z0m*exp(-k_over_B)

  end select

end subroutine cana_roughness

! ============================================================================
! calculate soil surface (laminar) resistances to evaporation and sensible heat
subroutine surface_resistances(tile, T_sfc, u_sfc, ustar_sfc, p, snow_active, &
       r_evap, r_sens)
  type(land_tile_type), intent(inout) :: tile
  real, intent(in) :: T_sfc     ! surface temperature, K
  real, intent(in) :: u_sfc     ! near-surface wind velocity, m/s
  real, intent(in) :: ustar_sfc ! friction velocity at the soil surface, m/s
  real, intent(in) :: p         ! surface pressure, N/m2
  logical, intent(in) :: snow_active ! if TRUE, ground is covered by snow
  ! output
  real, intent(out) :: r_evap ! surface resistance for evaporation, s/m
  real, intent(out) :: r_sens ! surface resistance for sensible heat, s/m

  real :: theta_sfc   ! relative soil wetness at the surface, unitless
  real :: r_litt_evap ! litter resistance, s/m
  real :: r_sv_evap   ! soil surface resistance to evaporation, s/m
  real :: r_bl_evap   ! viscous sublayer resistance to evaporation, s/m
  real :: r_bl_sens   ! viscous sublayer resistance to heat flux, s/m
  real :: d_visc      ! thickness of viscous sublayer, m
  real :: diff_air    ! thermal diffusivity of air, m/s

  r_litt_evap = evap_resistance_litter(tile, snow_active)

  if (is_watch_point()) then
     write(*,*) '#### surface resistance input ####'
     __DEBUG5__(T_sfc,u_sfc,ustar_sfc,p,snow_active)
     write(*,*) '#### end of surface resistance input ####'
  endif

  theta_sfc = 1.0 ! to avoid sending undefined values to diagnotics
  select case(soil_resistance_option)
  case(RESIST_NONE)
      r_sv_evap = 0
      r_bl_evap = 0
      r_bl_sens = 0
      d_visc    = 0
  case(RESIST_HO2013)
      d_visc    = sfc_visc_bl_depth(u_sfc, ustar_sfc, T_sfc, p)
      if (d_visc_max > 0) d_visc = min(d_visc,d_visc_max)
      diff_air  = thermal_diff_air(T_sfc)
      r_bl_sens = d_visc/diff_air

      r_sv_evap = 0.0
      r_bl_evap = d_visc/diffusivity_h2o(T_sfc,p)
      if (associated(tile%soil).and..not.snow_active) then
         r_sv_evap = soil_evap_sv_resistance(tile%soil)
         theta_sfc = max(0.0, tile%soil%wl(1) / (dens_h2o * dz(1)))/tile%soil%pars%vwc_sat
         r_bl_evap = soil_evap_bl_resistance(tile%soil, theta_sfc, T_sfc, p, d_visc)
      endif
  case default
     call error_mesg(module_name, 'invalid surface resistance option', FATAL)
  end select

  r_evap = r_bl_evap + r_sv_evap + r_litt_evap
  r_sens = r_bl_sens
  if (is_watch_point()) then
     __DEBUG5__(r_bl_evap, r_sv_evap, r_bl_sens, d_visc, diff_air)
     write(*,*) '#### surface resistance output ####'
     __DEBUG2__(r_evap, r_sens)
     write(*,*) '#### end of surface resistance output ####'
  endif

  ! diagnostic section
  call send_tile_data(id_r_litt_evap, limited(r_litt_evap),    tile%diag)
  call send_tile_data(id_r_bl_sens,   limited(r_bl_sens),      tile%diag)
  call send_tile_data(id_r_bl_evap,   limited(r_bl_evap),      tile%diag)
  call send_tile_data(id_r_sv_evap,   limited(r_sv_evap),      tile%diag)

  call send_tile_data(id_c_litt_evap, reciprocal(r_litt_evap), tile%diag)
  call send_tile_data(id_c_bl_sens,   reciprocal(r_bl_sens),   tile%diag)
  call send_tile_data(id_c_bl_evap,   reciprocal(r_bl_evap),   tile%diag)
  call send_tile_data(id_c_sv_evap,   reciprocal(r_sv_evap),   tile%diag)

  call send_tile_data(id_d_visc,      d_visc,                  tile%diag)
  call send_tile_data(id_u_sfc,       u_sfc,                   tile%diag)
  call send_tile_data(id_ustar_sfc,   ustar_sfc,               tile%diag)
  call send_tile_data(id_theta_sfc,   theta_sfc,               tile%diag)

  contains

  real elemental function limited(r)
     real, intent(in) :: r
     real, parameter :: max_val = 9999.0

     limited = max(min(r,max_val),-max_val)
  end function limited

  real elemental function reciprocal(r)
     real, intent(in) :: r

     real, parameter :: max_val = 9999.0
     real, parameter :: min_val = 1.0/max_val

     if (abs(r)<min_val) then
        reciprocal = sign(max_val,r) ! sign(a,b) returns value of a with sign of b
     else
        reciprocal = 1.0/r
     endif
  end function reciprocal

end subroutine surface_resistances

! ============================================================================
! additional resistance of litter to the water vapor flux.
! not a good parameterization, but just using for sensitivity analyses now.
! ignores differing biomass and litter turnover rates.
real function evap_resistance_litter(tile, snow_active) result(rav_lit)
  type(land_tile_type), intent(in) :: tile
  logical, intent(in) :: snow_active

  real :: litter_fast_C, litter_slow_C, litter_deadmic_C ! litter carbon pools

  rav_lit = 0.0
  if (.not.associated(tile%soil)) return
  if (snow_active)               return

  call get_soil_litter_C(tile%soil, litter_fast_C, litter_slow_C, litter_deadmic_C)
  rav_lit = rav_lit_0 + rav_lit_vi * (vegn_tile_LAI(tile%vegn)+vegn_tile_SAI(tile%vegn)) &
                      + rav_lit_fsc * litter_fast_C &
                      + rav_lit_ssc * litter_slow_C &
                      + rav_lit_deadmic * litter_deadmic_C &
                      + rav_lit_bwood * vegn_tile_bwood(tile%vegn)
end function evap_resistance_litter

! ============================================================================
! internal soil viscous resistance to evaporation.
! Haghighi et al. (2013): Evaporation rates across a convective air boundary layer
!     are dominated by diffusion. Water Resources Research, 49, No.3, 1602–1610,
!     doi:10.1002/wrcr.20166.
real function soil_evap_sv_resistance(soil) result(r_sv)
  type(soil_tile_type), intent(in) :: soil ! soil properties and parameters

  real, parameter :: gam = 1.73e-5 ! unit conversion constant (Haghighi et al. 2013) eq(13)

  real :: vlc(1), vsc(1) ! volumetric soil and ice content
  real :: psi(1)   ! soil matric potential, unused
  real :: DThDP(1), DKDP(1), DPsi_min, DPsi_max ! unused
  real :: K_z(1)   ! hydraulic conductivity in vertical, kg/(m2 s)
  real :: K_x(1)   ! hydraulic conductivity in horizontal, unused
  ! perhaps we can make soil_data_hydraulic_properties return only requested parameters?

  vlc(1) = max(0.0, soil%wl(1) / (dens_h2o * dz(1)))
  vsc(1) = max(0.0, soil%ws(1) / (dens_h2o * dz(1)))
  call soil_data_hydraulic_properties (soil, vlc, vsc, &
                   psi, DThDP, K_z, K_x, DKDP, DPsi_min, DPsi_max )
  ! calculate resistance to transport within soil upper layer
  r_sv = gam*dens_h2o/(4*K_z(1))
  if(is_watch_point()) then
     __DEBUG3__(psi(1),K_z(1),r_sv)
  endif
end function soil_evap_sv_resistance

! ============================================================================
! soil resistance to evaporation in the viscous sublayer
! Haghighi and Or (2013): Evaporation from porous surfaces into turbulent airflows: Coupling
!      eddy characteristics with pore scale vapor diffusion. Water Resources Research, 49,
!      8432–8442, doi:10.1002/2012WR013324View.
real function soil_evap_bl_resistance(soil, theta_sfc, T_sfc, p, d_bl) result(r_bl)
  type(soil_tile_type), intent(in) :: soil
  real, intent(in) :: theta_sfc ! relative soil wetness at the surface, unitless
  real, intent(in) :: T_sfc     ! surface temperature, K
  real, intent(in) :: p         ! pressure, N/m2
  real, intent(in) :: d_bl      ! thickness of viscous sublayer, m

  real :: f_diff   ! surface wetness dependency factor, unitless
  real :: diff_h2o ! diffusivity of water vapor

  if (theta_sfc > 0) then
     ! surface wetness dependency model, Schlunder (1988):
     f_diff = (sqrt(pi/(4*theta_sfc))-1)/(pi*sqrt(theta_sfc))
     ! to set reasonable value for theta_sfc > pi/4: f_diff = 0 means that water diffuses
     ! like from a completely wet surface
     f_diff = max(f_diff,0.0)
     ! diffusivity of water vapor for current conditions:
     diff_h2o = diffusivity_h2o(T_sfc,p)
     ! finally, calculate resistance
     r_bl = (d_bl + soil%r_pores*sqrt(pi)*f_diff)/diff_h2o
     if (is_watch_point()) then
        __DEBUG4__(f_diff,soil%r_pores,diff_h2o,r_bl)
     endif
  else ! theta_sfc <= 0
     ! returning HUGE value is not incorrect, but it causes overflow in the follow-up
     ! calculations of combined resistances/conductances. Therefore we just return an
     ! arbitrary big value.
     ! r_bl = HUGE(r_bl)
     r_bl = 1e10
     if (is_watch_point()) then
        __DEBUG1__(r_bl)
     endif
  endif
end function soil_evap_bl_resistance

! ============================================================================
! average depth of viscous sublayer
! Haghighi and Or (2013): Evaporation from porous surfaces into turbulent airflows: Coupling
!      eddy characteristics with pore scale vapor diffusion. Water Resources Research, 49,
!      8432–8442, doi:10.1002/2012WR013324View.
real function sfc_visc_bl_depth(u_sfc, ustar_sfc, T, p) result(d_visc)
  real, intent(in) :: u_sfc     ! near-surface wind velocity, m/s
  real, intent(in) :: ustar_sfc ! friction velocity at the soil surface, m/s
  real, intent(in) :: T         ! temperature, degK
  real, intent(in) :: p         ! pressure, N/m2

  real, parameter :: c2 = 2.2   ! parameter of viscous sublayer depth (Haghighi & Or 2013),
                                ! eq (11), (Haghighi & Or 2015) eq (10a)
  real, parameter :: c3 = 112.0 ! parameter of average eddy exposure time (Haghighi & Or,
                                ! 2013) eq (15)

  real :: alpha    ! parameter of eddy exposure time probability distribution, unitless
  real :: visc

  ! parameter of eddy exposure time probability distribution, (Haghighi & Or 2013) eq (A4)
  ! for ustar_sfc > 0.3 u_sfc we set alpha to zero, effectively turning gamma-distribution
  ! to exponential distribution.
  alpha = max(0.3 * u_sfc/ustar_sfc-1.0,0.0)
  ! kinematic viscosity of air
  visc = kin_visc_air(T,p)
  ! it is more robust (and efficient) to calculate the ratio of gamma-functions as the
  ! exponent of logarithm differences: computing gamma function for large arguments can
  ! lead to overflow.
  ! d_visc = visc/ustar_sfc * c2*sqrt(c3)/sqrt(alpha+1) * gamma(alpha+1.5)/gamma(alpha+1)
  d_visc = visc/ustar_sfc * c2*sqrt(c3)/sqrt(alpha+1) * exp(gammaln(alpha+1.5)-gammaln(alpha+1))
  if (is_watch_point()) then
  __DEBUG5__(u_sfc, ustar_sfc, alpha, visc, d_visc)
  endif
end function sfc_visc_bl_depth


! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function cana_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   cana_tile_exists = associated(tile%cana)
end function cana_tile_exists

subroutine cana_T_ptr(t,p)
  type(land_tile_type), pointer :: t
  real,                 pointer :: p

  p=>NULL()
  if(associated(t))then
     if(associated(t%cana))p=>t%cana%T
  endif
end subroutine

subroutine cana_fog_ptr(t,p)
  type(land_tile_type), pointer :: t
  real,                 pointer :: p

  p=>NULL()
  if(associated(t))then
     if(associated(t%cana))p=>t%cana%fog
  endif
end subroutine

subroutine cana_tr_ptr(t,i,p)
  type(land_tile_type), pointer    :: t
  integer,              intent(in) :: i
  real,                 pointer    :: p

  p=>NULL()
  if(associated(t))then
     if(associated(t%cana))p=>t%cana%tr(i)
  endif
end subroutine

end module canopy_air_mod
