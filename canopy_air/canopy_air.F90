!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Land Model 4 (LM4).
!*
!* LM4 is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* LM4 is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with LM4.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
! ============================================================================
! canopy air
! ============================================================================
module canopy_air_mod

#include "../shared/debug.inc"

use mpp_mod, only: input_nml_file
use fms_mod, only : error_mesg, FATAL, NOTE, &
     check_nml_error, mpp_pe, mpp_root_pe, stdlog, string
use constants_mod, only : VONKARM
use field_manager_mod, only : parse, MODEL_ATMOS, MODEL_LAND
use tracer_manager_mod, only : get_tracer_index, get_tracer_names, &
     query_method, NO_TRACER

use land_constants_mod, only : mol_CO2, mol_air
use land_tracers_mod, only : ntcana, isphum, ico2
use cana_tile_mod, only : cana_tile_type, &
     canopy_air_mass, canopy_air_mass_for_tracers, cpw
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_data_mod, only : log_version
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, get_tile_data, field_exists
use land_debug_mod, only : is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public :: read_cana_namelist
public :: cana_init
public :: cana_end
public :: save_cana_restart
public :: cana_turbulence
public :: cana_roughness
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'canopy_air_mod'
#include "../shared/version_variable.inc"

! options for turbulence parameter calculations
integer, parameter :: TURB_LM3W = 1, TURB_LM3V = 2
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
real :: bare_rah_sca     = 0.01     ! bare-ground resistance between ground and canopy air, s/m
namelist /cana_nml/ &
  init_T, init_T_cold, init_q, init_co2, turbulence_to_use, use_SAI_for_heat_exchange, &
  canopy_air_mass, canopy_air_mass_for_tracers, cpw, save_qco2, bare_rah_sca
!---- end of namelist --------------------------------------------------------

logical :: module_is_initialized =.FALSE.
integer :: turbulence_option ! selected option of turbulence parameters
     ! calculations

contains

! ============================================================================
subroutine read_cana_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
     read (input_nml_file, nml=cana_nml, iostat=io)
     ierr = check_nml_error(io, 'cana_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=cana_nml)
  endif
end subroutine read_cana_namelist

! ============================================================================
! initialize canopy air
subroutine cana_init ()

  ! ---- local vars ----------------------------------------------------------
  type(land_tile_enum_type)     :: ce ! last and current tile
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  character(*), parameter :: restart_file_name='INPUT/cana.nc'
  type(land_restart_type) :: restart
  logical :: restart_exists

  character(32)  :: name  ! name of the tracer
  integer        :: tr, i ! tracer indices
  real           :: init_tr(ntcana) ! initial (cold-start) values of tracers
  real           :: value ! used for parameter parsing
  character(32)  :: scheme
  character(1024) :: parameters

  module_is_initialized = .TRUE.

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
     tile%cana%tr(:) = init_tr(:)
  enddo

  ! then read the restart if it exists
  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then
     call error_mesg('cana_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     call get_tile_data(restart, 'temp', cana_T_ptr)
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

  ! initialize options, to avoid expensive character comparisons during
  ! run-time
  if (trim(turbulence_to_use)=='lm3v') then
     turbulence_option = TURB_LM3V
  else if (trim(turbulence_to_use)=='lm3w') then
     turbulence_option = TURB_LM3W
  else
     call error_mesg('cana_init', 'canopy air turbulence option turbulence_to_use="'// &
          trim(turbulence_to_use)//'" is invalid, use "lm3w" or "lm3v"', FATAL)
  endif
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
  filename = 'RESTART/'//trim(timestamp)//'cana.nc'
  call init_land_restart(restart, filename, cana_tile_exists, tile_dim_length)

  ! write temperature
  call add_tile_data(restart,'temp',cana_T_ptr,'canopy air temperature','degrees_K')
  do tr = 1,ntcana
     call get_tracer_names(MODEL_LAND, tr, name, longname, units)
     if (tr==ico2.and..not.save_qco2) cycle
     call add_tile_data(restart,name,cana_tr_ptr,tr,'canopy air '//trim(longname),trim(units))
  enddo
  call save_land_restart(restart)
  call free_land_restart(restart)
end subroutine save_cana_restart

! ============================================================================
subroutine cana_turbulence (u_star,&
     vegn_cover, vegn_layerfrac, vegn_height, vegn_bottom, vegn_lai, vegn_sai, vegn_d_leaf, &
     land_d, land_z0m, land_z0s, grnd_z0s, &
     con_v_h, con_v_v, con_g_h, con_g_v )
  real, intent(in) :: &
       u_star, & ! friction velocity, m/s
       land_d, land_z0m, land_z0s, grnd_z0s, &
       vegn_cover, vegn_height(:), vegn_layerfrac(:), &
       vegn_bottom(:), & ! height of the bottom of the canopy, m
       vegn_lai(:), vegn_sai(:), vegn_d_leaf(:)
  real, intent(out) :: &
       con_v_h(:), con_v_v(:), & ! one-sided foliage-CAS conductance per unit ground area
       con_g_h   , con_g_v       ! ground-CAS conductance per unit ground area

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

  select case(turbulence_option)
  case(TURB_LM3W)
     if(vegn_cover > 0) then
        ztop   = maxval(vegn_height(:))

        wind  = u_star/VONKARM*log((ztop-land_d)/land_z0m) ! normalized wind on top of the canopy
        a     = vegn_cover*a_max
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

     vegn_idx = sum((vegn_lai+vegn_sai)*vegn_layerfrac)  ! total vegetation index
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
     land_d, land_z0m, land_z0s )
  logical, intent(in) :: lm2
  real, intent(in) :: &
       subs_z0m, subs_z0s, snow_z0m, snow_z0s, snow_area, vegn_cover, vegn_height, &
       vegn_lai, vegn_sai
  real, intent(out) :: &
       land_d    ,&
       land_z0m  ,&
       land_z0s

  !---- local constants
  real, parameter :: d_h_max = 2./3.
  real, parameter :: z0m_h_max = 1/7.35

  ! ---- local vars
  real :: d_h      ! ratio of displacement height to vegetation height
  real :: z0m_h    ! ratio of roughness length to vegetation height
  real :: grnd_z0m, grnd_z0s
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
     land_z0s = land_z0m*exp(-2.0)

  end select

end subroutine cana_roughness

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
