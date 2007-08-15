! ============================================================================
! glac model module
! ============================================================================
module glacier_mod

use fms_mod,            only: error_mesg, file_exist,  open_namelist_file,    &
                              check_nml_error, stdlog, write_version_number, &
                              close_file, mpp_pe, mpp_root_pe, FATAL, NOTE
use time_manager_mod,   only: time_type, increment_time, time_type_to_real
use diag_manager_mod,   only: diag_axis_init
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o, PI

use glac_tile_mod,      only: glac_tile_type, glac_pars_type, glac_prog_type, &
     read_glac_data_namelist, glac_data_thermodynamics, glac_data_hydraulics, &
     glac_data_radiation, glac_data_diffusion, max_lev

use land_constants_mod, only : &
     NBANDS
use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=)
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type
use land_data_mod,      only : land_state_type, lnd
use land_io_mod, only : print_netcdf_error
use land_tile_io_mod, only: create_tile_out_file, read_tile_data_r1d_fptr, write_tile_data_r1d_fptr
use nf_utils_mod, only : nfu_def_dim, nfu_put_att_text
use land_debug_mod, only : is_watch_point
implicit none
private

! ==== public interfaces =====================================================
public :: read_glac_namelist
public :: glac_init
public :: glac_end
public :: glac_get_sfc_temp
public :: glac_radiation
public :: glac_diffusion
public :: glac_step_1
public :: glac_step_2
! =====end of public interfaces ==============================================



! ==== module constants ======================================================
character(len=*), parameter :: &
       module_name = 'glacier',&
       version = '',&
       tagname = ''
 
! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2                  = .false.
logical :: use_bucket           = .false.    ! single-layer glac water
logical :: bifurcate            = .false.    ! consider direct evap from bucket
logical :: conserve_glacier_mass = .false.
real    :: init_temp            = 260.       ! cold-start glac T
real    :: init_w               = 150.       ! cold-start w(l)/dz(l)
real    :: init_groundwater     =   0.       ! cold-start gw storage
   real :: cpw = 1952.  ! specific heat of water vapor at constant pressure
   real :: clw = 4218.  ! specific heat of water (liquid)
   real :: csw = 2106.  ! specific heat of water (ice)

namelist /glac_nml/ lm2, use_bucket,             bifurcate,  &
                    conserve_glacier_mass,           &
                    init_temp,      &
                    init_w,       &
                    init_groundwater, cpw, clw, csw
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
type(time_type) :: time
real            :: delta_time       ! fast time step

logical         :: use_beta, use_glac_rh

integer         :: num_l            ! # of water layers
real            :: dz    (max_lev)  ! thicknesses of layers
real            :: zfull (max_lev)
real            :: zhalf (max_lev+1)


! ---- diagnostic field IDs
integer :: id_zhalf, id_zfull
integer :: id_lwc, id_swc, id_temp, id_ie, id_sn, id_bf, id_hie, id_hsn, id_hbf

! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_glac_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l            ! level iterator

  call read_glac_data_namelist(num_l, dz)

  call write_version_number(version, tagname)
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=glac_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'glac_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  if (mpp_pe() == mpp_root_pe()) then
     write (stdlog(), nml=glac_nml)
  endif

  ! ---- set up vertical discretization
  zhalf(1) = 0
  do l = 1, num_l;   
     zhalf(l+1) = zhalf(l) + dz(l)
     zfull(l) = 0.5*(zhalf(l+1) + zhalf(l))
  enddo

end subroutine read_glac_namelist


! ============================================================================
! initialize glacier model
subroutine glac_init ( id_lon, id_lat )
  integer, intent(in)  :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in)  :: id_lat  ! ID of land latitude (Y) axis

  ! ---- local vars
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce ! last and current tile list elements
  type(land_tile_type), pointer :: tile  ! pointer to current tile

  module_is_initialized = .TRUE.
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)

  ! ---- set flags for use in water-availability constraints -----------------
  use_beta    = .false.
  use_glac_rh = .false.
!  if (use_bucket .and. .not.bifurcate) then    !************* TEMPORARY!!!!
!      use_beta = .true.
!    else if (use_bucket .and. bifurcate) then
!    else
!      use_glac_rh = .true.
!    endif
  use_glac_rh = .true.
  use_beta = .true.

  ! -------- initialize glac state --------
  if (file_exist('INPUT/glac.res.nc')) then
     call error_mesg('glac_init','reading NetCDF restart',NOTE)
     __NF_ASRT__(nf_open('INPUT/glac.res.nc',NF_NOWRITE,unit))
     call read_tile_data_r1d_fptr(unit, 'temp'         , glac_temp_ptr  )
     call read_tile_data_r1d_fptr(unit, 'wl'           , glac_wl_ptr )
     call read_tile_data_r1d_fptr(unit, 'ws'           , glac_ws_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater'  , glac_gw_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater_T', glac_gwT_ptr)
     __NF_ASRT__(nf_close(unit))     
  else
     te = tail_elmt (lnd%tile_map)
     ce = first_elmt(lnd%tile_map)
     do while(ce /= te)
        tile=>current_tile(ce) ! get pointer to current tile
        ce=next_elmt(ce)       ! advance position to the next tile
        
        if (.not.associated(tile%glac)) cycle

        if (init_temp.ge.tfreeze.or.lm2) then      ! USE glac TFREEZE HERE
           tile%glac%prog(1:num_l)%wl = init_w*dz(1:num_l)
           tile%glac%prog(1:num_l)%ws = 0
        else
           tile%glac%prog(1:num_l)%wl = 0
           tile%glac%prog(1:num_l)%ws = init_w*dz(1:num_l)
        endif
        tile%glac%prog%T             = init_temp
        tile%glac%prog%groundwater   = init_groundwater
        tile%glac%prog%groundwater_T = init_temp
     enddo
  endif

  call glac_diag_init ( id_lon, id_lat, zfull(1:num_l), zhalf(1:num_l+1) )

end subroutine glac_init


! ============================================================================
subroutine glac_end (tile_dim_length)
  integer, intent(in) :: tile_dim_length ! length of tile dimension in the 
                                         ! output file

  integer :: unit             ! unit for i/o
  logical :: restart_created ! flag indicating that the restart file was created

  module_is_initialized =.FALSE.

  ! ---- write restart file --------------------------------------------------
  
  ! create output file, including internal structure necessary for output
  call error_mesg('glac_end','writing NetCDF restart',NOTE)
  call create_tile_out_file(unit,'RESTART/glac.res.nc', &
          lnd%glon*180.0/PI, lnd%glat*180/PI, glac_tile_exists, tile_dim_length,&
          created=restart_created)

  if (restart_created) then
     ! in addition, define vertical coordinate
     __NF_ASRT__(nfu_def_dim(unit,'zfull',zfull(1:num_l),'full level','m'))
     __NF_ASRT__(nfu_put_att_text(unit,'zfull','positive','down'))
        
     ! write out fields
     call write_tile_data_r1d_fptr(unit,'temp'         ,glac_temp_ptr,'zfull','glacier temperature','degrees_K')
     call write_tile_data_r1d_fptr(unit,'wl'           ,glac_wl_ptr  ,'zfull','liquid water content','kg/m2')
     call write_tile_data_r1d_fptr(unit,'ws'           ,glac_ws_ptr  ,'zfull','solid water content','kg/m2')
     call write_tile_data_r1d_fptr(unit,'groundwater'  ,glac_gw_ptr  ,'zfull')
     call write_tile_data_r1d_fptr(unit,'groundwater_T',glac_gwT_ptr ,'zfull')
   
     ! close file
     __NF_ASRT__(nf_close(unit))
  endif

end subroutine glac_end


! ============================================================================
! returns glacier surface temperature
subroutine glac_get_sfc_temp ( glac, glac_T )
  type(glac_tile_type), intent(in)  :: glac
  real,                 intent(out) :: glac_T

  glac_T = glac%prog(1)%T
end subroutine glac_get_sfc_temp


! ============================================================================
subroutine glac_radiation ( glac, &
     glac_refl_dir, glac_refl_dif, glac_refl_lw, glac_emis )
  type(glac_tile_type), intent(in) :: glac
  real, intent(out) :: &
       glac_refl_dir(NBANDS), glac_refl_dif(NBANDS), & ! glacier albedos for direct and diffuse light
       glac_refl_lw,   &  ! glacier reflectance for longwave (thermal) radiation
       glac_emis          ! glacier emissivity

  call glac_data_radiation ( glac, glac_refl_dir, glac_refl_dif, glac_emis )
  glac_refl_lw = 1 - glac_emis
end subroutine glac_radiation


! ============================================================================
! compute glac-only properties needed to do glac-canopy-atmos energy balance
subroutine glac_diffusion ( glac, glac_z0s, glac_z0m )
  type(glac_tile_type), intent(in) :: glac
  real, intent(out) :: glac_z0s, glac_z0m

  call glac_data_diffusion ( glac, glac_z0s, glac_z0m )
  
end subroutine glac_diffusion


! ============================================================================
! update glac properties explicitly for time step.
! MAY WISH TO INTRODUCE 'UNIVERSAL' SENSITIVITIES FOR SIMPLICITY.
! T-DEPENDENCE OF HYDRAULIC PROPERTIES COULD BE DONE LESS FREQUENTLY.
! integrate glac-heat conduction equation upward from bottom of glac
! to surface, delivering linearization of surface ground heat flux.
subroutine glac_step_1 ( glac, &
                         glac_T, &
                         glac_rh, glac_liq, glac_ice, glac_subl, glac_tf, glac_G0, &
                         glac_DGDT, conserve_glacier_mass_out )
  type(glac_tile_type),intent(inout) :: glac
  real, intent(out) :: &
       glac_T, &
       glac_rh, glac_liq, glac_ice, glac_subl, &
       glac_tf, & ! freezing temperature of glacier, degK
       glac_G0, &
       glac_DGDT
  logical, intent(out) :: conserve_glacier_mass_out

  ! ---- local vars 
  real                   :: bbb, denom, dt_e
  real, dimension(num_l) :: aaa, ccc, thermal_cond, heat_capacity, vlc, vsc
  integer :: l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  conserve_glacier_mass_out = conserve_glacier_mass

  if(is_watch_point()) then
    write(*,*) 'checkpoint gs1 a'
    write(*,*) 'mask    ',  .TRUE.
    write(*,*) 'T       ', glac%prog(1)%T
  endif

  glac_T = glac%prog(1)%T

  if(is_watch_point()) then
     write(*,*) 'checkpoint gs1 b'
     write(*,*) 'mask    ', .TRUE.
     write(*,*) 'glac_T       ', glac_T
  endif

  do l = 1, num_l
     vlc(l) = max(0.0, glac%prog(l)%wl / (dens_h2o * dz(l)))
     vsc(l) = max(0.0, glac%prog(l)%ws / (dens_h2o * dz(l)))
  enddo

  call glac_data_thermodynamics ( glac%pars, vlc(1), vsc(1),  &  
       glac_rh, glac%heat_capacity_dry, thermal_cond )

  do l = 1, num_l
     heat_capacity(l) = glac%heat_capacity_dry(l)*dz(l) &
          + clw*glac%prog(l)%wl + csw*glac%prog(l)%ws
  enddo

  if (.not.use_glac_rh) glac_rh=1.0

  if (lm2) then
     glac_liq = 0
     glac_ice = 1.e6
  else
     glac_liq  = max(glac%prog(1)%wl, 0.0)
     glac_ice  = max(glac%prog(1)%ws, 0.0)
  endif
  if (glac_liq + glac_ice > 0 ) then
     glac_subl = glac_ice / (glac_liq + glac_ice)
  else
     glac_subl = 0
  endif
  
  if(num_l > 1) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz(l+1)/thermal_cond(l+1) &
                + dz(l)/thermal_cond(l)   )
        aaa(l+1) = - dt_e * delta_time / heat_capacity(l+1)
        ccc(l)   = - dt_e * delta_time / heat_capacity(l)
     enddo

     bbb = 1.0 - aaa(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(glac%prog(num_l)%T - glac%prog(num_l-1)%T)
     glac%e(num_l-1) = -aaa(num_l)/denom
     glac%f(num_l-1) = dt_e/denom
     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*glac%e(l)
        dt_e = - ( ccc(l)*(glac%prog(l+1)%T - glac%prog(l)%T  ) &
             -aaa(l)*(glac%prog(l)%T   - glac%prog(l-1)%T) )
        glac%e(l-1) = -aaa(l)/denom
        glac%f(l-1) = (dt_e - ccc(l)*glac%f(l))/denom
     enddo
     denom = delta_time/(heat_capacity(1) )
     glac_G0    = ccc(1)*(glac%prog(2)%T- glac%prog(1)%T &
          + glac%f(1)) / denom
     glac_DGDT  = (1 - ccc(1)*(1-glac%e(1))) / denom   
  else  ! one-level case
     denom = delta_time/heat_capacity(1)
     glac_G0    = 0.
     glac_DGDT  = 1. / denom
  endif

  ! set freezing temperature of glaciers
  glac_tf = glac%pars%tfreeze

  if(is_watch_point())then
     write(*,*) 'checkpoint gs1 c'
     write(*,*) 'mask    ', .TRUE.
     write(*,*) 'T       ', glac_T
     write(*,*) 'rh      ', glac_rh
     write(*,*) 'liq     ', glac_liq
     write(*,*) 'ice     ', glac_ice
     write(*,*) 'subl    ', glac_subl
     write(*,*) 'G0      ', glac_G0
     write(*,*) 'DGDT    ', glac_DGDT
     do l = 1, num_l
        write(*,*) 'T(dbg,l)', glac%prog(l)%T
     enddo

  endif
end subroutine glac_step_1


! ============================================================================
! apply boundary flows to glac water and move glac water vertically.
  subroutine glac_step_2 ( glac, diag, glac_subl, snow_lprec, snow_hlprec,  &
                           subs_DT, subs_M_imp, subs_evap, &
                           glac_levap, glac_fevap, glac_hadvec, glac_melt, &
                           glac_lrunf, glac_hlrunf, glac_Ttop, glac_Ctop, &
                           glac_LMASS, glac_FMASS, glac_HEAT )

  type(glac_tile_type), intent(inout) :: glac
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: &
     glac_subl     !
  real, intent(in) :: &
     snow_lprec, &
     snow_hlprec, &
     subs_DT,       &!
     subs_M_imp,       &! rate of phase change of non-evaporated glac water
     subs_evap
  real, intent(out) :: &
     glac_levap, glac_fevap, glac_hadvec, glac_melt, &
     glac_lrunf, glac_hlrunf, glac_Ttop, glac_Ctop, &
     glac_LMASS, glac_FMASS, glac_HEAT

  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l) :: del_t, eee, fff, &
             psi, DThDP, hyd_cond, DKDP, K, DKDPm, DKDPp, grad, &
             vlc, vsc, dW_l, u_minus, u_plus, DPsi, glac_w_fc
  real, dimension(num_l+1) :: flow
  real, dimension(num_l  ) :: div
  real :: &
     lprec_eff, hlprec_eff, tflow, hcap,cap_flow, &
     melt_per_deg, melt,&
     lrunf_sn,lrunf_ie,lrunf_bf, hlrunf_sn,hlrunf_ie,hlrunf_bf, &
     Qout, DQoutDP,&
     tau_gw, c0, c1, c2, x, aaa, bbb, ccc, ddd, xxx, Dpsi_min, Dpsi_max
  logical :: stiff
  real, dimension(num_l-1) :: del_z
  integer :: l
  real :: jj
  ! --------------------------------------------------------------------------

  jj = 1.
  DPsi = 0.0
  c1   = 0.0
  
  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 1 ***** '
     write(*,*) 'mask    ', .TRUE.
     write(*,*) 'subs_evap    ', subs_evap
     write(*,*) 'snow_lprec   ', snow_lprec
     write(*,*) 'subs_M_imp   ', subs_M_imp
     write(*,*) 'theta_s ', glac%pars%w_sat
     do l = 1, num_l
        write(*,'(i2.2,99(a,g))')l,&
             ' T =', glac%prog(l)%T,&
             ' Th=', (glac%prog(l)%ws+glac%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', glac%prog(l)%wl,&
             ' ws=', glac%prog(l)%ws,&
             ' gw=', glac%prog(l)%groundwater
     enddo
     
  endif

  ! ---- record fluxes ---------
  glac_levap  = subs_evap*(1-glac_subl)
  glac_fevap  = subs_evap*   glac_subl
  glac_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution --------------
  del_t(1) = subs_DT
  glac%prog(1)%T = glac%prog(1)%T + del_t(1)
  if ( num_l > 1) then
    do l = 1, num_l-1
      del_t(l+1) = glac%e(l) * del_t(l) + glac%f(l)
      glac%prog(l+1)%T = glac%prog(l+1)%T + del_t(l+1)
    end do
  end if

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 2 ***** '
     write(*,*) 'levap=',glac_levap
     write(*,*) 'fevap=',glac_fevap
     write(*,*) 'subs_M_imp=',subs_M_imp
     do l = 1, num_l
        write(*,'(i2.2,x,a,g)') l, 'T', glac%prog(l)%T
     enddo
  endif


  glac_hadvec = 0.
  ! ---- extract evap from glac and do implicit melt --------------------
  IF(LM2)THEN
!      glac%prog(1)%wl = glac%prog(1)%wl - glac_levap*delta_time
!      glac%prog(1)%ws = glac%prog(1)%ws - glac_fevap*delta_time
      glac_hadvec = -cpw*subs_evap*(glac%prog(1)%T-tfreeze)
!      hcap = glac%heat_capacity_dry(1)*dz(1) &
!                         + clw*glac%prog(1)%wl + csw*glac%prog(1)%ws
!      glac%prog(1)%T = glac%prog(1)%T + (   &
!                    +((clw-cpw)*glac_levap                              &
!                    + (csw-cpw)*glac_fevap)*(glac%prog(1)%T  -tfreeze) &
!                                                 )*delta_time/ hcap
!      glac%prog(1)%wl = glac%prog(1)%wl + subs_M_imp
!      glac%prog(1)%ws = glac%prog(1)%ws - subs_M_imp
!      glac%prog(1)%T  = tfreeze + (hcap*(glac%prog(1)%T-tfreeze) ) &
!                                / ( hcap + (clw-csw)*subs_M_imp )
  ELSE
    glac%prog(1)%wl = glac%prog(1)%wl - glac_levap*delta_time
    glac%prog(1)%ws = glac%prog(1)%ws - glac_fevap*delta_time
    glac_hadvec = -cpw*subs_evap*(glac%prog(1)%T-tfreeze)
    hcap = glac%heat_capacity_dry(1)*dz(1) &
                       + clw*glac%prog(1)%wl + csw*glac%prog(1)%ws
    glac%prog(1)%T = glac%prog(1)%T + (   &
                  +((clw-cpw)*glac_levap                              &
                  + (csw-cpw)*glac_fevap)*(glac%prog(1)%T  -tfreeze) &
                                               )*delta_time/ hcap
    glac%prog(1)%wl = glac%prog(1)%wl + subs_M_imp
    glac%prog(1)%ws = glac%prog(1)%ws - subs_M_imp
    glac%prog(1)%T  = tfreeze + (hcap*(glac%prog(1)%T-tfreeze) ) &
                              / ( hcap + (clw-csw)*subs_M_imp )
  ENDIF

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3 ***** '
     do l = 1, num_l
        write(*,'(i2.2,99(a,g))') l,&
             ' T =', glac%prog(l)%T,&
             ' wl=', glac%prog(l)%wl,&
             ' ws=', glac%prog(l)%ws
     enddo
  endif

  ! ---- fetch glac hydraulic properties -------------------------------------
  vlc=1;vsc=0
  do l = 1, num_l
     vlc(l) = max(0., glac%prog(l)%wl / (dens_h2o*dz(l)))
     vsc(l) = max(0., glac%prog(l)%ws / (dens_h2o*dz(l)))
  enddo
  call glac_data_hydraulics (glac, vlc, vsc, &
                   psi, DThDP, hyd_cond, DKDP, Dpsi_min, Dpsi_max, tau_gw, &
                   glac_w_fc )
  IF (lm2) THEN ! ********************************
     if(is_watch_point()) then
        write(*,*) ' ***** glac_step_2 checkpoint 3.1 ***** '
        do l = 1, num_l
           write(*,'(i2.2,99(x,a,g))') l, 'vlc', vlc(l),&
                'K  ', hyd_cond(l)
        enddo
     
     endif
     ! ---- remainder of mass fluxes and associated sensible heat fluxes --------
     flow = 1
     flow(1)  = snow_lprec *delta_time + subs_M_imp
     do l = 1, num_l
        flow(l+1) = flow(l)
     enddo
     dW_l = 0
  ELSE   ! ********************************
    div = 0.
    do l = 1, num_l
      div(l) = 0.15*dens_h2o*dz(l)/tau_gw
    enddo
    lrunf_bf = sum(div)

  ! ---- glac-water flow ----------------------------------------------------
    stiff = all(DThDP.eq.0)
    if (snow_lprec/=0 .and. psi(num_l)>0) then
      lrunf_sn = snow_lprec*min((psi(num_l)/zhalf(num_l))**glac%pars%rsa_exp,1.)
      hlrunf_sn = lrunf_sn*snow_hlprec/snow_lprec
    else
      lrunf_sn = 0.
      hlrunf_sn = 0.
    endif
    lprec_eff = snow_lprec - lrunf_sn
    hlprec_eff = snow_hlprec - hlrunf_sn
    flow(1) = delta_time*lprec_eff
    do l = 1, num_l-1
      del_z(l) = zfull(l+1)-zfull(l)
      K(l) = 0.5*(hyd_cond(l)+hyd_cond(l+1))
      DKDPm(l) = 0. !0.5*DKDP(l)
      DKDPp(l) = 0. ! 0.5*DKDP(l+1)
!        K(l) = hyd_cond(l)
!        DKDPm(l) = DKDP(l)
!        DKDPp(l) = 0
      grad(l)  = jj*(psi(l+1)-psi(l))/del_z(l) - 1
    enddo


    if(is_watch_point()) then
       write(*,*) ' ***** glac_step_2 checkpoint 3.1 ***** '
       do l = 1, num_l
          write(*,'(i2.2,x,a,99g)') l, 'DThDP,hyd_cond,psi,DKDP', &
               DThDP(l), hyd_cond(l), psi(l), DKDP(l)
       enddo
       do l = 1, num_l-1
          write(*,'(i2.2,x,a,99g)') l, 'K,DKDPm,DKDPp,grad,del_z', &
               K(l), DKDPm(l), DKDPp(l), grad(l)
      enddo
    endif

    l = num_l
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    aaa =     - ( jj* K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
!      where (stiff)
        bbb = xxx - (- jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
        ddd = - K(l-1) *grad(l-1) - div(l)
!        elsewhere
!          Qout = hyd_cond(l) ! gravity drainage
!          DQoutDP = DKDP(l)  ! gravity drainage
!          Qout = 0.                ! no drainage
!          DQoutDP = 0.             ! no drainage
!          where (psi(l).gt.0.) ! linear baseflow from gw
!              Qout = 0.15*psi(l)/tau_gw
!              DQoutDP = 0.15/tau_gw
!            elsewhere
!              Qout = 0.
!              DQoutDP = 0.
!            endwhere
!          bbb = xxx - (- jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
!                      -DQoutDP )
!          ddd = -Qout - K(l-1) *grad(l-1)
!        endwhere
    eee(l-1) = -aaa/bbb
    fff(l-1) =  ddd/bbb
    
    if(is_watch_point()) then
       write(*,'(a,i,99g)') 'l,a,b, ,d', l,aaa, bbb,ddd       
    endif


    do l = num_l-1, 2, -1
      xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
      aaa = - ( jj*K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
      bbb = xxx-( -jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
                  -jj*K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
      ccc =   - (  jj*K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
      ddd =       K(l)*grad(l) - K(l-1)*grad(l-1) &
                            - div(l)
      eee(l-1) =                    -aaa/(bbb+ccc*eee(l))
      fff(l-1) =  (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
      if(is_watch_point()) then
         write(*,'(a,i,99g)') 'l,a,b,c,d', l,aaa, bbb,ccc,ddd
      endif
    enddo

    l = 1
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    bbb = xxx - ( -jj*K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
    ccc =     - (  jj*K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
    ddd =          flow(1)/delta_time +    K(l)     *grad(l) &
                            - div(l)
    if (stiff) then
      dPsi(l) =  - psi(l)
    else
      dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
      dPsi(l) = min (dPsi(l), Dpsi_max)
      dPsi(l) = max (dPsi(l), Dpsi_min)
    endif
    flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                      - K(l)*grad(l))*delta_time
    lrunf_ie         = lprec_eff - flow(l)/delta_time

    if(is_watch_point()) then
       write(*,'(a,i,99g)') 'l,  b,c,d', l, bbb,ccc,ddd

       write(*,*) ' ***** glac_step_2 checkpoint 3.2 ***** '
       write(*,*) 'ie,sn,bf:', lrunf_ie,lrunf_sn,lrunf_bf
       do l = 1, num_l-1
          write(*,'(a,i,99g)') 'l,eee(l),fff(l)', l,eee(l), fff(l)
       enddo
       write(*,*) 'DThDP(1)', DThDP(1)
       write(*,*) 'ddd(1)', ddd
       write(*,*) 'ccc(1)', ccc
       write(*,*) 'bbb(1)', bbb
       write(*,*) 'dPsi(1)', dPsi(1)
       write(*,*) 'Psi(1)', Psi(1)
    endif

    do l = 2, num_l
      dPsi(l) = eee(l-1)*dPsi(l-1) + fff(l-1)
    enddo

    do l = 1, num_l-1
      flow(l+1) = delta_time*( &
           -K(l)*(grad(l)&
           +jj*(DPsi(l+1)-DPsi(l))/ del_z(l)) &
           -grad(l)*(DKDPp(l)*Dpsi(l+1)+ &
                           DKDPm(l)*Dpsi(l) )  )
      dW_l(l) = flow(l) - flow(l+1) - div(l)*delta_time
      glac%prog(l)%wl = glac%prog(l)%wl + dW_l(l)
    enddo
!  where (stiff)
    flow(num_l+1) = 0.
!    elsewhere
!      flow(num_l+1) = (Qout &
!                         + DQoutDP*DPsi(num_l)) * delta_time
!    endwhere
    dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                          - div(num_l)*delta_time
    glac%prog(num_l)%wl = glac%prog(num_l)%wl + dW_l(num_l)
  
  ENDIF ! ************************************

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',glac%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,'(i2.2,99(a,g))')l,&
             ' Th=', (glac%prog(l)%ws+glac%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', glac%prog(l)%wl,&
             ' ws=', glac%prog(l)%ws,&
             'Dpsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif

  glac_hadvec = glac_hadvec + snow_hlprec
  if  (snow_hlprec.ne.0.) then
    tflow = tfreeze + snow_hlprec/(clw*snow_lprec)
  else
    tflow = tfreeze
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.4 ***** '
     write(*,*) ' tfreeze', tfreeze
     write(*,*) ' snow_hlprec', snow_hlprec
  endif

! For initial testing, use top-down-flow weights to advect heat.
  u_minus = 1.
  u_plus  = 0.
  if (flow(1).lt.0.) u_minus(1) = 0.
  hcap = (glac%heat_capacity_dry(num_l)*dz(num_l) &
                              + csw*glac%prog(num_l)%ws)/clw
  aaa = -flow(num_l) * u_minus(num_l)
  bbb =  hcap + glac%prog(num_l)%wl - dW_l(num_l) - aaa
  eee(num_l-1) = -aaa/bbb
  fff(num_l-1) = aaa*(glac%prog(num_l)%T-glac%prog(num_l-1)%T) / bbb

  do l = num_l-1, 2, -1
    hcap = (glac%heat_capacity_dry(l)*dz(l) &
                              + csw*glac%prog(l)%ws)/clw
    aaa = -flow(l)   * u_minus(l)
    ccc =  flow(l+1) * u_plus (l)
    bbb =  hcap + glac%prog(l)%wl - dW_l(l) - aaa - ccc
    eee(l-1) = -aaa / ( bbb +ccc*eee(l) )
    fff(l-1) = (   aaa*(glac%prog(l)%T-glac%prog(l-1)%T)    &
                       + ccc*(glac%prog(l)%T-glac%prog(l+1)%T)    &
                       - ccc*fff(l) ) / ( bbb +ccc*eee(l) )
  enddo
    
  hcap = (glac%heat_capacity_dry(1)*dz(1) + csw*glac%prog(1)%ws)/clw
  aaa = -flow(1) * u_minus(1)
  ccc =  flow(2) * u_plus (1)
  bbb =  hcap + glac%prog(1)%wl - dW_l(1) - aaa - ccc

  del_t(1) =  (  aaa*(glac%prog(1)%T-tflow          ) &
                     + ccc*(glac%prog(1)%T-glac%prog(2)%T) &
                     - ccc*fff(1) ) / (bbb+ccc*eee(1))
  glac%prog(1)%T = glac%prog(1)%T + del_t(1)

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.4.1 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'aaa', aaa
     write(*,*) 'bbb', bbb
     write(*,*) 'ccc', ccc
     write(*,*) 'del_t(1)', del_t(1)
     write(*,*) ' T(1)', glac%prog(1)%T
  endif

  do l = 1, num_l-1
    del_t(l+1) = eee(l)*del_t(l) + fff(l)
    glac%prog(l+1)%T = glac%prog(l+1)%T + del_t(l+1)
  enddo

  tflow = glac%prog(num_l)%T

!  do l = 1, num_l
!    where (mask)
!        hcap = glac%heat_capacity_dry(l)*dz(l) &
!                 + clw*(glac%prog(l)%wl-dW_l(l)) + csw*glac%prog(l)%ws
!        cap_flow = clw*flow(l)
!        glac%prog(l)%T = (hcap*glac%prog(l)%T + cap_flow*tflow) &
!                         /(hcap                 + cap_flow      )
!        tflow  = glac%prog(l)%T
!      endwhere
!    enddo

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.5 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'cap_flow', cap_flow
     do l = 1, num_l
        write(*,'(i2.2,99(a,g))')l, ' T', glac%prog(l)%T, ' flow ',flow(l)
     enddo
     write(*,*) 'delta_time,tau_gw,c0,c1,c2,x', delta_time,tau_gw,c0,&
          c1,c2,x
     write(*,*) 'level=', num_l+1, ' flow ',flow(num_l+1)
     write(*,*) 'gw(1)',glac%prog(1)%groundwater
  endif

  ! ---- groundwater ---------------------------------------------------------
  ! THIS T AVERAGING IS WRONG, BECAUSE IT NEGLECTS THE MEDIUM  ***
  ! ALSO, FREEZE-THAW IS NEEDED!
  ! PROBABLY THIS SECTION WILL BE DELETED ANYWAY, WITH GW TREATED ABOVE.
  IF (lm2) THEN
    if (conserve_glacier_mass) flow(num_l+1) = snow_lprec*delta_time
    do l = 1, 1      !TEMPORARY LAYER THING !!!!***
      if (glac%prog(l)%groundwater + flow(num_l+1) .ne. 0.) then ! TEMP FIX
         glac%prog(l)%groundwater_T =    &
          + (glac%prog(l)%groundwater*glac%prog(l)%groundwater_T &
             + flow(num_l+1)*tflow) &
           /(glac%prog(l)%groundwater + flow(num_l+1))
      endif
      c0 = delta_time/tau_gw
      c1 = exp(-c0)
      c2 = (1-c1)/c0
      x  = (1-c1)*glac%prog(l)%groundwater/delta_time &
                          + (1-c2)*flow(num_l+1)/delta_time
      glac%prog(l)%groundwater = c1 * glac%prog(l)%groundwater &
                                + c2 * flow(num_l+1)
      glac_lrunf  = x
      glac_hlrunf = x*clw*(glac%prog(l)%groundwater_T-tfreeze)
      glac_hadvec = glac_hadvec - glac_hlrunf
    enddo
    if(is_watch_point()) then
       write(*,*) ' ***** glac_step_2 checkpoint 3.6 ***** '
       write(*,*) 'delta_time,tau_gw,c0,c1,c2,x',&
            delta_time,tau_gw,c0,c1,c2,x
       write(*,*) 'level=', num_l+1, ' flow ',flow(num_l+1)
       write(*,*) 'gw(1)',glac%prog(1)%groundwater
       write(*,*) 'hcap', hcap
       write(*,*) 'cap_flow', cap_flow
       do l = 1, num_l
          write(*,*) 'level=', l, ' T', glac%prog(l)%T
       enddo
    endif

  ELSE
    if (lprec_eff.ne.0. .and. flow(1).ge.0. ) then
      hlrunf_ie = lrunf_ie*hlprec_eff/lprec_eff
    else if (flow(1).lt.0. ) then
      hlrunf_ie = hlprec_eff - (flow(1)/delta_time)*clw &
                         *(glac%prog(1)%T-tfreeze)
    else
      hlrunf_ie = 0.
    endif
    hlrunf_bf = clw*sum(div*(glac%prog%T-tfreeze))
    glac_lrunf  = lrunf_sn + lrunf_ie + lrunf_bf
    glac_hlrunf = hlrunf_sn + hlrunf_bf + hlrunf_ie 
    glac_hadvec = glac_hadvec - glac_hlrunf
  ENDIF

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.7 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'cap_flow', cap_flow
     do l = 1, num_l
        write(*,'(i2.2,99(a,g))')l, ' T', glac%prog(l)%T
     enddo
  endif

  IF(LM2) THEN
  ELSE
    do l = 1, num_l
    ! ---- compute explicit melt/freeze --------------------------------------
      hcap = glac%heat_capacity_dry(l)*dz(l) &
               + clw*glac%prog(l)%wl + csw*glac%prog(l)%ws
      melt_per_deg = hcap/hlf
      if (glac%prog(l)%ws>0 .and. glac%prog(l)%T>glac%pars%tfreeze) then
        melt =  min(glac%prog(l)%ws, (glac%prog(l)%T-glac%pars%tfreeze)*melt_per_deg)
      else if (glac%prog(l)%wl>0 .and. glac%prog(l)%T<glac%pars%tfreeze) then
        melt = -min(glac%prog(l)%wl, (glac%pars%tfreeze-glac%prog(l)%T)*melt_per_deg)
      else
        melt = 0
      endif

      if(is_watch_point()) then
         write(*,'(a,i,99g)') 'l,T,wl(1),ws(1),melt:', l,glac%prog(l)%T, glac%prog(l)%wl, &
              glac%prog(l)%ws, melt
      endif

      glac%prog(l)%wl = glac%prog(l)%wl + melt
      glac%prog(l)%ws = glac%prog(l)%ws - melt
      glac%prog(l)%T = tfreeze &
         + (hcap*(glac%prog(l)%T-tfreeze) - hlf*melt) &
                              / ( hcap + (clw-csw)*melt )
      if(is_watch_point()) then
         write(*,'(a,i,99g)') 'l,T,wl(1),ws(1):', l,glac%prog(l)%T, glac%prog(l)%wl, &
              glac%prog(l)%ws
      endif

      glac_melt = glac_melt + melt / delta_time
    enddo
  ENDIF

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 5 ***** '
     write(*,*) 'i,j,k,melt:',&
          glac_melt*delta_time
     do l = 1, num_l
        write(*,'(i2.2,99(a,g))')l, &
             ' T =', glac%prog(l)%T, &
             ' Th=', (glac%prog(l)%ws+glac%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', glac%prog(l)%wl,&
             ' ws=', glac%prog(l)%ws,&
             ' gw=', glac%prog(l)%groundwater
     enddo
  endif

  glac_LMASS = 0
  glac_FMASS = 0
  glac_HEAT = 0
  do l = 1, num_l
    glac_LMASS = glac_LMASS &
         +      glac%prog(l)%wl
    glac_FMASS = glac_FMASS &
         +      glac%prog(l)%ws
    glac_HEAT = glac_HEAT &
         + (glac%heat_capacity_dry(l)*dz(l) &
             + clw*glac%prog(l)%wl + csw*glac%prog(l)%ws)  &
                        * (glac%prog(l)%T-tfreeze)
  enddo
  glac_LMASS = glac_LMASS +     glac%prog(1)%groundwater
  glac_HEAT = glac_HEAT + clw*glac%prog(1)%groundwater    &
                                  *(glac%prog(1)%groundwater_T-tfreeze)
  glac_Ttop = glac%prog(1)%T
  glac_Ctop = glac%heat_capacity_dry(1)*dz(1) &
    + clw*glac%prog(1)%wl + csw*glac%prog(1)%ws

! ----------------------------------------------------------------------------
! given solution for surface energy balance, write diagnostic output.
!  

  ! ---- increment time
  time = increment_time(time, int(delta_time), 0)
  
  ! ---- diagnostic section
  call send_tile_data (id_temp, glac%prog%T,     diag )
  call send_tile_data (id_lwc,  glac%prog(1:num_l)%wl/dz(1:num_l), diag )
  call send_tile_data (id_swc,  glac%prog(1:num_l)%wl/dz(1:num_l), diag )
  call send_tile_data (id_ie,   lrunf_ie,        diag )
  call send_tile_data (id_sn,   lrunf_sn,        diag )
  call send_tile_data (id_bf,   lrunf_bf,        diag )
  call send_tile_data (id_hie,  hlrunf_ie,       diag )
  call send_tile_data (id_hsn,  hlrunf_sn,       diag )
  call send_tile_data (id_hbf,  hlrunf_bf,       diag )

end subroutine glac_step_2

! ============================================================================
subroutine glac_diag_init ( id_lon, id_lat, zfull, zhalf )
  integer,         intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer,         intent(in) :: id_lat  ! ID of land longitude (X) axis
  real,            intent(in) :: zfull(:)! Full levels, m
  real,            intent(in) :: zhalf(:)! Half levels, m

  ! ---- local vars ----------------------------------------------------------
  integer :: axes(3)

  ! define vertical axis
  id_zhalf = diag_axis_init ( &
       'glac_zhalf', zhalf, 'meters', 'z', 'half level',  -1, set_name='glac' )
  id_zfull = diag_axis_init ( &
       'glac_zfull', zfull, 'meters', 'z', 'full level',  -1, set_name='glac', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/ id_lon, id_lat, id_zfull /)

  ! define diagnostic fields
  id_lwc = register_tiled_diag_field ( module_name, 'glac_liq', axes,        &
       Time, 'bulk density of liquid water', 'kg/m3', missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'glac_ice',  axes,      &
       Time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'glac_T',  axes,       &
       Time, 'temperature',            'degK',  missing_value=-100.0 )
  id_ie  = register_tiled_diag_field ( module_name, 'glac_rie',  axes(1:2),  &
       Time, 'inf exc runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_sn  = register_tiled_diag_field ( module_name, 'glac_rsn',  axes(1:2),  &
       Time, 'satn runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_bf  = register_tiled_diag_field ( module_name, 'glac_rbf',  axes(1:2),  &
       Time, 'baseflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_hie  = register_tiled_diag_field ( module_name, 'glac_hie',  axes(1:2), &
       Time, 'heat ie runf',            'W/m2',  missing_value=-100.0 )
  id_hsn  = register_tiled_diag_field ( module_name, 'glac_hsn',  axes(1:2), &
       Time, 'heat sn runf',            'W/m2',  missing_value=-100.0 )
  id_hbf  = register_tiled_diag_field ( module_name, 'glac_hbf',  axes(1:2), &
       Time, 'heat bf runf',            'W/m2',  missing_value=-100.0 )

  ! define static diagnostic fields

  
end subroutine glac_diag_init

! ============================================================================
! tile existance detector: returns a logical value indicating wether component
! model tile exists or not
logical function glac_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   glac_tile_exists = associated(tile%glac)
end function glac_tile_exists


! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine glac_temp_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%glac)) ptr=>tile%glac%prog%T
   endif
end subroutine glac_temp_ptr

subroutine glac_wl_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%glac)) ptr=>tile%glac%prog%wl
   endif
end subroutine glac_wl_ptr

subroutine glac_ws_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%glac)) ptr=>tile%glac%prog%ws
   endif
end subroutine glac_ws_ptr

subroutine glac_gw_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%glac)) ptr=>tile%glac%prog%groundwater
   endif
end subroutine glac_gw_ptr

subroutine glac_gwT_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%glac)) ptr=>tile%glac%prog%groundwater_T
   endif
end subroutine glac_gwT_ptr

end module glacier_mod



