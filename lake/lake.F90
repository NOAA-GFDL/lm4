! ============================================================================
! lake model module
! ============================================================================
module lake_mod

use fms_mod, only : error_mesg, file_exist,  open_namelist_file, check_nml_error, &
     stdlog, write_version_number, close_file, mpp_pe, mpp_root_pe, FATAL, NOTE, &
     get_mosaic_tile_file
use time_manager_mod,   only: time_type, increment_time, time_type_to_real
use diag_manager_mod,   only: diag_axis_init, register_diag_field,           &
                              register_static_field, send_data
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o, PI

use land_constants_mod, only : &
     NBANDS
use lake_tile_mod, only : &
     lake_tile_type, lake_pars_type, lake_prog_type, read_lake_data_namelist, &
     lake_data_radiation, lake_data_diffusion, &
     lake_data_thermodynamics, lake_data_hydraulics, &
     max_lev, cpw,clw,csw
use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=)
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type
use land_data_mod,      only : land_state_type, lnd
use land_tile_io_mod, only : print_netcdf_error, create_tile_out_file, &
     read_tile_data_r1d_fptr, write_tile_data_r1d_fptr 
use nf_utils_mod, only : nfu_def_dim, nfu_put_att_text
use land_debug_mod, only: is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public :: read_lake_namelist
public :: lake_init
public :: lake_end

public :: lake_get_sfc_temp
public :: lake_radiation
public :: lake_diffusion
public :: lake_step_1
public :: lake_step_2
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'lake',&
    version     = '$Id: lake.F90,v 15.0.2.4 2008/02/17 21:01:29 slm Exp $',&
    tagname     = '$Name: omsk_2008_03 $'

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2                  = .false.
logical :: use_bucket           = .false.    ! single-layer lake water
logical :: bifurcate            = .false.    ! consider direct evap from bucket
real    :: init_temp            = 288.        ! cold-start lake T
real    :: init_w               = 1000.      ! cold-start w(l)/dz(l)
real    :: init_groundwater     =   0.        ! cold-start gw storage

namelist /lake_nml/ lm2, use_bucket,             bifurcate,             &
                    init_temp,      &
                    init_w,       &
                    init_groundwater, cpw, clw, csw
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
type(time_type) :: time
real            :: delta_time

logical         :: use_beta, use_lake_rh

integer         :: num_l              ! # of water layers
real            :: dz    (max_lev)    ! thicknesses of layers
real            :: zfull (max_lev)
real            :: zhalf (max_lev+1)

! ---- diagnostic field IDs
integer :: id_lwc, id_swc, id_temp, id_ie, id_sn, id_bf, id_hie, id_hsn, id_hbf

! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_lake_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l

  call read_lake_data_namelist(num_l,dz)

  call write_version_number(version, tagname)
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
  if (mpp_pe() == mpp_root_pe()) then
     write (stdlog(), nml=lake_nml)
  endif

  ! ---- set up vertical discretization
  zhalf(1) = 0
  do l = 1, num_l;   
     zhalf(l+1) = zhalf(l) + dz(l)
     zfull(l) = 0.5*(zhalf(l+1) + zhalf(l))
  enddo

end subroutine read_lake_namelist


! ============================================================================
! initialize lake model
subroutine lake_init ( id_lon, id_lat )
  integer, intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in) :: id_lat  ! ID of land latitude (Y) axis

  ! ---- local vars 
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce ! last and current tile list elements
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  character(len=256) :: restart_file_name

  module_is_initialized = .TRUE.
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)

  ! ---- set flags for use in water-availability constraints -----------------
  use_beta    = .false.
  use_lake_rh = .false.
!  if (use_bucket .and. .not.bifurcate) then    !************* TEMPORARY!!!!
!      use_beta = .true.
!    else if (use_bucket .and. bifurcate) then
!    else
!      use_lake_rh = .true.
!    endif
  use_lake_rh = .true.
  use_beta = .true.

  ! -------- initialize lake state --------
  call get_mosaic_tile_file('INPUT/lake.res.nc',restart_file_name,.FALSE.,lnd%domain)
  if (file_exist(restart_file_name)) then
     call error_mesg('lake_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
     call read_tile_data_r1d_fptr(unit, 'temp'         , lake_temp_ptr  )
     call read_tile_data_r1d_fptr(unit, 'wl'           , lake_wl_ptr )
     call read_tile_data_r1d_fptr(unit, 'ws'           , lake_ws_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater'  , lake_gw_ptr )
     call read_tile_data_r1d_fptr(unit, 'groundwater_T', lake_gwT_ptr)
     __NF_ASRT__(nf_close(unit))     
  else
     call error_mesg('lake_init',&
          'cold-starting lake',&
          NOTE)
     te = tail_elmt (lnd%tile_map)
     ce = first_elmt(lnd%tile_map)
     do while(ce /= te)
        tile=>current_tile(ce)  ! get pointer to current tile
        ce=next_elmt(ce)        ! advance position to the next tile
        
        if (.not.associated(tile%lake)) cycle
        
        if (init_temp.ge.tfreeze) then      ! USE lake TFREEZE HERE
           tile%lake%prog(1:num_l)%wl = init_w*dz(1:num_l)
           tile%lake%prog(1:num_l)%ws = 0
        else
           tile%lake%prog(1:num_l)%wl = 0
           tile%lake%prog(1:num_l)%ws = init_w*dz(1:num_l)
        endif
        tile%lake%prog%T             = init_temp
        tile%lake%prog%groundwater   = init_groundwater
        tile%lake%prog%groundwater_T = init_temp
     enddo
  endif

  call lake_diag_init ( id_lon, id_lat )

end subroutine lake_init


! ============================================================================
subroutine lake_end (tile_dim_length)
  integer, intent(in) :: tile_dim_length ! length of tile dimension in the 
                                         ! output file

  ! ---- local vars 
  integer :: unit         ! unit for i/o
  logical :: restart_created ! flag indicating that the restart file was created

  module_is_initialized =.FALSE.

  ! ---- write restart file --------------------------------------------------

  ! create output file, including internal structure necessary for output
  call error_mesg('lake_end','writing NetCDF restart',NOTE)
  call create_tile_out_file(unit,'RESTART/lake.res.nc', &
          lnd%glon*180.0/PI, lnd%glat*180/PI, lake_tile_exists, tile_dim_length,&
          created=restart_created)

  if (restart_created) then
     ! in addition, define vertical coordinate
     __NF_ASRT__(nfu_def_dim(unit,'zfull',zfull(1:num_l),'full level','m'))
     __NF_ASRT__(nfu_put_att_text(unit,'zfull','positive','down'))
        
     ! write out fields
     call write_tile_data_r1d_fptr(unit,'temp'         ,lake_temp_ptr,'zfull','glacier temperature','degrees_K')
     call write_tile_data_r1d_fptr(unit,'wl'           ,lake_wl_ptr  ,'zfull','liquid water content','kg/m2')
     call write_tile_data_r1d_fptr(unit,'ws'           ,lake_ws_ptr  ,'zfull','solid water content','kg/m2')
     call write_tile_data_r1d_fptr(unit,'groundwater'  ,lake_gw_ptr  ,'zfull')
     call write_tile_data_r1d_fptr(unit,'groundwater_T',lake_gwT_ptr ,'zfull')
   
     ! close file
     __NF_ASRT__(nf_close(unit))
  endif

end subroutine lake_end


! ============================================================================
subroutine lake_get_sfc_temp(lake, lake_T)
  type(lake_tile_type), intent(in) :: lake
  real, intent(out) :: lake_T

  lake_T = lake%prog(1)%T
end subroutine lake_get_sfc_temp


! ============================================================================
! compute lake-only radiation properties
subroutine lake_radiation ( lake, &
     lake_refl_dir, lake_refl_dif, lake_refl_lw, lake_emis )
  type(lake_tile_type), intent(in) :: lake
  real, intent(out) :: lake_refl_dir(NBANDS), lake_refl_dif(NBANDS), lake_refl_lw, lake_emis

  call lake_data_radiation ( lake, lake_refl_dir, lake_refl_dif, lake_emis )
  lake_refl_lw = 1 - lake_emis
end subroutine lake_radiation


! ============================================================================
! compute lake-only roughness parameters
subroutine lake_diffusion ( lake, lake_z0s, lake_z0m )
  type(lake_tile_type), intent(in) :: lake
  real, intent(out) :: lake_z0s, lake_z0m

  call lake_data_diffusion ( lake, lake_z0s, lake_z0m )
end subroutine lake_diffusion

! ============================================================================
! update lake properties explicitly for time step.
! MAY WISH TO INTRODUCE 'UNIVERSAL' SENSITIVITIES FOR SIMPLICITY.
! T-DEPENDENCE OF HYDRAULIC PROPERTIES COULD BE DONE LESS FREQUENTLY.
! integrate lake-heat conduction equation upward from bottom of lake
! to surface, delivering linearization of surface ground heat flux.
subroutine lake_step_1 ( lake, &
                         lake_T, &
                         lake_rh, lake_liq, lake_ice, lake_subl, lake_tf, lake_G0, &
                         lake_DGDT )

  type(lake_tile_type), intent(inout) :: lake
  real, intent(out)  :: &
       lake_T, &
       lake_rh, lake_liq, lake_ice, lake_subl, &
       lake_tf, & ! freezing temperature of lake, degK
       lake_G0, &
       lake_DGDT

  ! ---- local vars
  real                  :: bbb, denom, dt_e, vlc_sfc, vsc_sfc
  real, dimension(num_l):: aaa, ccc, thermal_cond, heat_capacity, vlc, vsc
  integer               :: l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  lake_T = lake%prog(1)%T
  do l = 1, num_l
     vlc(l) = max(0., lake%prog(l)%wl / (dens_h2o * dz(l)))
     vsc(l) = max(0., lake%prog(l)%ws / (dens_h2o * dz(l)))
  enddo

  vlc_sfc = vlc(1)
  vsc_sfc = vsc(1)

  call lake_data_thermodynamics ( lake%pars, vlc_sfc, vsc_sfc,  &  
                                  lake_rh, &
                                  lake%heat_capacity_dry, thermal_cond )
  do l = 1, num_l
     heat_capacity(l) = lake%heat_capacity_dry(l) *dz(l) &
          + clw*lake%prog(l)%wl + csw*lake%prog(l)%ws
  enddo
  if (.not.use_lake_rh) lake_rh=1

  if (lm2) then
     lake_liq = 0
     lake_ice = 0
  else
     lake_liq  = max(lake%prog(1)%wl, 0.)
     lake_ice  = max(lake%prog(1)%ws, 0.)
  endif
  if (lake_liq + lake_ice > 0 ) then
     lake_subl = lake_ice / (lake_liq + lake_ice)
  else
     lake_subl = 0
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
     dt_e = aaa(num_l)*(lake%prog(num_l)%T - lake%prog(num_l-1)%T)
     lake%e(num_l-1) = -aaa(num_l)/denom
     lake%f(num_l-1) = dt_e/denom
     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*lake%e(l)
        dt_e = - ( ccc(l)*(lake%prog(l+1)%T - lake%prog(l)%T  ) &
                  -aaa(l)*(lake%prog(l)%T   - lake%prog(l-1)%T) )
        lake%e(l-1) = -aaa(l)/denom
        lake%f(l-1) = (dt_e - ccc(l)*lake%f(l))/denom
     end do
     denom = delta_time/(heat_capacity(1) )
     lake_G0    = ccc(1)*(lake%prog(2)%T- lake%prog(1)%T &
          + lake%f(1)) / denom
     lake_DGDT  = (1 - ccc(1)*(1-lake%e(1))) / denom   
  else  ! one-level case
     denom = delta_time/heat_capacity(1)
     lake_G0    = 0.
     lake_DGDT  = 1. / denom
  end if
  
  ! set the freezing temperature of the lake
  lake_tf = lake%pars%tfreeze
  
  if(is_watch_point()) then
     write(*,*) 'lake_step_1 checkpoint 1, lake_fe(dbg)'
     write(*,*) 'mask    ', .true.   
     write(*,*) 'T       ', lake_T
     write(*,*) 'rh      ', lake_rh
     write(*,*) 'liq     ', lake_liq
     write(*,*) 'ice     ', lake_ice
     write(*,*) 'subl    ', lake_subl
     write(*,*) 'G0      ', lake_G0
     write(*,*) 'DGDT    ', lake_DGDT
     do l = 1, num_l
        write(*,*) 'uptake_frac(dbg,l)', 0.0,&
                   'T(dbg,l)          ', lake%prog(l)%T
     enddo
  endif

end subroutine lake_step_1


! ============================================================================
! apply boundary flows to lake water and move lake water vertically.
  subroutine lake_step_2 ( lake, diag, lake_subl, snow_lprec, snow_hlprec,  &
                           subs_DT, subs_M_imp, subs_evap, &
                           lake_levap, lake_fevap, lake_hadvec, lake_melt, &
                           lake_lrunf, lake_hlrunf, lake_Ttop, lake_Ctop, &
                           lake_LMASS, lake_FMASS, lake_HEAT )
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
  real, intent(out) :: &
     lake_levap, lake_fevap, lake_hadvec, lake_melt, &
     lake_lrunf, lake_hlrunf, lake_Ttop, lake_Ctop, &
     lake_LMASS, lake_FMASS, lake_HEAT

  ! ---- local vars
  real, dimension(num_l) :: del_t, eee, fff, &
             psi, DThDP, hyd_cond, DKDP, K, DKDPm, DKDPp, grad, &
             vlc, vsc, dW_l, u_minus, u_plus, DPsi, lake_w_fc
  real, dimension(num_l+1) :: flow
  real, dimension(num_l  ) :: div
  real :: &
     lprec_eff, hlprec_eff, tflow, hcap,cap_flow, &
     melt_per_deg, melt,&
     lrunf_sn,lrunf_ie,lrunf_bf, hlrunf_sn,hlrunf_ie,hlrunf_bf, &
     Qout, DQoutDP,&
     tau_gw, c0, c1, c2, x, aaa, bbb, ccc, ddd, xxx, Dpsi_min, Dpsi_max
  logical  :: stiff
  real, dimension(num_l-1) :: del_z
  integer :: l
  real :: jj

  jj = 1.
  
  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 1 ***** '
    write(*,*) 'mask    ', .true.   
    write(*,*) 'subs_evap    ', subs_evap   
    write(*,*) 'snow_lprec   ', snow_lprec  
    write(*,*) 'subs_M_imp   ', subs_M_imp   
    write(*,*) 'theta_s ', lake%pars%w_sat
    do l = 1, num_l
      write(*,*) ' level=', l,&
                 ' T =', lake%prog(l)%T,&
                 ' Th=', (lake%prog(l)%ws &
                         +lake%prog(l)%wl)/(dens_h2o*dz(l)),&
                 ' wl=', lake%prog(l)%wl,&
                 ' ws=', lake%prog(l)%ws,&
                 ' gw=', lake%prog(l)%groundwater
      enddo
  endif

  ! ---- record fluxes ---------
  lake_levap  = subs_evap*(1-lake_subl)
  lake_fevap  = subs_evap*   lake_subl
  lake_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution --------------
  del_t(1) = subs_DT
  lake%prog(1)%T = lake%prog(1)%T + del_t(1)
  if ( num_l > 1) then
    do l = 1, num_l-1
      del_t(l+1) = lake%e(l) * del_t(l) + lake%f(l)
      lake%prog(l+1)%T = lake%prog(l+1)%T + del_t(l+1)
    end do
  end if

  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 2 ***** '
    do l = 1, num_l
       write(*,*) 'level=', l, 'T', lake%prog(l)%T
    enddo
  endif

  lake_hadvec = 0.
  ! ---- extract evap from lake and do implicit melt --------------------
  lake%prog(1)%wl = lake%prog(1)%wl - lake_levap*delta_time
  lake%prog(1)%ws = lake%prog(1)%ws - lake_fevap*delta_time
  lake_hadvec = -cpw*subs_evap*(lake%prog(1)%T-tfreeze)
  hcap = lake%heat_capacity_dry(1)*dz(1) &
                     + clw*lake%prog(1)%wl + csw*lake%prog(1)%ws
  lake%prog(1)%T = lake%prog(1)%T + (   &
                +((clw-cpw)*lake_levap                              &
                + (csw-cpw)*lake_fevap)*(lake%prog(1)%T  -tfreeze) &
                                             )*delta_time/ hcap
  lake%prog(1)%wl = lake%prog(1)%wl + subs_M_imp
  lake%prog(1)%ws = lake%prog(1)%ws - subs_M_imp
  lake%prog(1)%T  = tfreeze + (hcap*(lake%prog(1)%T-tfreeze) ) &
                            / ( hcap + (clw-csw)*subs_M_imp )

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 2.1 ***** '
     do l = 1, num_l
        write(*,*) 'level=', l, 'T', lake%prog(l)%T
     enddo
  endif

  ! ---- fetch lake hydraulic properties -------------------------------------
  vlc=0;vsc=0
  do l = 1, num_l
    vlc(l) = max(0., lake%prog(l)%wl / (dens_h2o*dz(l)))
    vsc(l) = max(0., lake%prog(l)%ws / (dens_h2o*dz(l)))
  enddo
  call lake_data_hydraulics (lake, vlc, vsc, &
                   psi, DThDP, hyd_cond, DKDP, Dpsi_min, Dpsi_max, tau_gw, &
                   lake_w_fc )

  IF (lm2) THEN ! ********************************

     if(is_watch_point()) then
        write(*,*) ' ***** lake_step_2 checkpoint 3.1 ***** '
        do l = 1, num_l
           write(*,*) 'level=', l, 'vlc', vlc(l), 'K  ', hyd_cond(l)
        enddo
     endif
  ! ---- remainder of mass fluxes and associated sensible heat fluxes --------
    flow=1
    flow(1)  = snow_lprec *delta_time
    do l = 1, num_l
!        flow(l+1) = max(0., lake%prog(l)%wl + flow(l) &
!                       - lake_w_fc(l)*dz(l)*dens_h2o)
!        flow(l+1) = flow(l)
      flow(l+1) = lake%prog(l)%wl + flow(l) &
                     - lake_w_fc(l)*dz(l)*dens_h2o
      dW_l(l) = flow(l) - flow(l+1)
      lake%prog(l)%wl = lake%prog(l)%wl + dW_l(l)
    enddo
  ELSE   ! ********************************
    div = 0.
    do l = 1, num_l
      if (vsc(l).eq.0. .and. psi(l).gt.0.) then
        div(l) = 0.15*dens_h2o*dz(l)/tau_gw
      endif
    enddo
    lrunf_bf = sum(div)

    ! ---- lake-water flow ----------------------------------------------------
    stiff = all(DThDP.eq.0)
    if (snow_lprec.ne.0. .and. psi(num_l).gt.0.) then
      lrunf_sn = snow_lprec*min((psi(num_l)/zhalf(num_l))**lake%pars%rsa_exp,1.)
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
       write(*,*) ' ***** lake_step_2 checkpoint 3.1 ***** '
       do l = 1, num_l
          write(*,*) 'level=', l, 'DThDP,hyd_cond,psi,DKDP', &
               DThDP(l),&
               hyd_cond(l),&
               psi(l),&
               DKDP(l)
       enddo
       do l = 1, num_l-1
          write(*,*) 'interface=', l, 'K,DKDPm,DKDPp,grad,del_z', &
               K(l),&
               DKDPm(l),&
               DKDPp(l),&
               grad(l)
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

    if(is_watch_point()) write(*,*) 'l,a,b, ,d', l,aaa, &
                       bbb,ddd

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
      if(is_watch_point()) write(*,*) 'l,a,b,c,d', l,aaa, &
                       bbb,ccc,ddd
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

    if(is_watch_point()) write(*,*) 'l,  b,c,d', l, &
                       bbb,ccc,ddd

    if(is_watch_point()) then
       write(*,*) ' ***** lake_step_2 checkpoint 3.2 ***** '
       write(*,*) 'ie,sn,bf:', lrunf_ie,lrunf_sn,lrunf_bf
       do l = 1, num_l-1
          write(*,*) 'l,eee(l),fff(l)', l,eee(l), fff(l)
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
      lake%prog(l)%wl = lake%prog(l)%wl + dW_l(l)
    enddo
!  where (stiff)
    flow(num_l+1) = 0.
!    elsewhere
!      flow(num_l+1) = (Qout &
!                         + DQoutDP*DPsi(num_l)) * delta_time
!    endwhere
    dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                            - div(num_l)*delta_time
    lake%prog(num_l)%wl = lake%prog(num_l)%wl + dW_l(num_l)
  
  ENDIF ! ************************************

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',lake%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,*) ' level=', l,&
             ' Th=', (lake%prog(l)%ws +lake%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', lake%prog(l)%wl,&
             'Dpsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif
  
  lake_hadvec = lake_hadvec + snow_hlprec
  if (snow_lprec.ne.0.) then
    tflow = tfreeze + snow_hlprec/(clw*snow_lprec)
  else
    tflow = tfreeze
  endif

  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 3.4 ***** '
    write(*,*) ' tfreeze', tfreeze
    write(*,*) ' snow_hlprec', snow_hlprec
  endif

! For initial testing, use top-down-flow weights to advect heat.
  u_minus = 1.
  u_plus  = 0.
  if (flow(1).lt.0.) u_minus(1) = 0.
  hcap = (lake%heat_capacity_dry(num_l)*dz(num_l) &
                              + csw*lake%prog(num_l)%ws)/clw
  aaa = -flow(num_l) * u_minus(num_l)
  bbb =  hcap + lake%prog(num_l)%wl - dW_l(num_l) - aaa
  eee(num_l-1) = -aaa/bbb
  fff(num_l-1) = aaa*(lake%prog(num_l)%T-lake%prog(num_l-1)%T) / bbb

  do l = num_l-1, 2, -1
    hcap = (lake%heat_capacity_dry(l)*dz(l) &
                              + csw*lake%prog(l)%ws)/clw
    aaa = -flow(l)   * u_minus(l)
    ccc =  flow(l+1) * u_plus (l)
    bbb =  hcap + lake%prog(l)%wl - dW_l(l) - aaa - ccc
    eee(l-1) = -aaa / ( bbb +ccc*eee(l) )
    fff(l-1) = (   aaa*(lake%prog(l)%T-lake%prog(l-1)%T)    &
                       + ccc*(lake%prog(l)%T-lake%prog(l+1)%T)    &
                       - ccc*fff(l) ) / ( bbb +ccc*eee(l) )
  enddo
    
  hcap = (lake%heat_capacity_dry(1)*dz(1) + csw*lake%prog(1)%ws)/clw
  aaa = -flow(1) * u_minus(1)
  ccc =  flow(2) * u_plus (1)
  bbb =  hcap + lake%prog(1)%wl - dW_l(1) - aaa - ccc

  del_t(1) =  (  aaa*(lake%prog(1)%T-tflow          ) &
                     + ccc*(lake%prog(1)%T-lake%prog(2)%T) &
                     - ccc*fff(1) ) / (bbb+ccc*eee(1))
  lake%prog(1)%T = lake%prog(1)%T + del_t(1)

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 3.4.1 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'aaa', aaa
     write(*,*) 'bbb', bbb
     write(*,*) 'ccc', ccc
     write(*,*) 'del_t(1)', del_t(1)
     write(*,*) ' T(1)', lake%prog(1)%T
  endif

  do l = 1, num_l-1
    del_t(l+1) = eee(l)*del_t(l) + fff(l)
    lake%prog(l+1)%T = lake%prog(l+1)%T + del_t(l+1)
  enddo

  tflow = lake%prog(num_l)%T

!  do l = 1, num_l
!    where (mask)
!        hcap = lake%heat_capacity_dry(l)*dz(l) &
!                 + clw*(lake%prog(l)%wl-dW_l(l)) + csw*lake%prog(l)%ws
!        cap_flow = clw*flow(l)
!        lake%prog(l)%T = (hcap*lake%prog(l)%T + cap_flow*tflow) &
!                         /(hcap                 + cap_flow      )
!        tflow  = lake%prog(l)%T
!      endwhere
!    enddo


  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 3.5 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'cap_flow', cap_flow
     do l = 1, num_l
        write(*,*) 'level=', l, ' T', lake%prog(l)%T
     enddo
  endif

  ! ---- groundwater ---------------------------------------------------------
  ! THIS T AVERAGING IS WRONG, BECAUSE IT NEGLECTS THE MEDIUM  ***
  ! ALSO, FREEZE-THAW IS NEEDED!
  ! PROBABLY THIS SECTION WILL BE DELETED ANYWAY, WITH GW TREATED ABOVE.
  IF (lm2) THEN
    do l = 1, 1      !TEMPORARY LAYER THING !!!!***
      if (lake%prog(l)%groundwater + flow(num_l+1) .ne. 0.) then ! TEMP FIX
          lake%prog(l)%groundwater_T =    &
           + (lake%prog(l)%groundwater*lake%prog(l)%groundwater_T &
              + flow(num_l+1)*tflow) &
            /(lake%prog(l)%groundwater + flow(num_l+1))
      endif
      c0 = delta_time/tau_gw
      c1 = exp(-c0)
      c2 = (1-c1)/c0
      x  = (1-c1)*lake%prog(l)%groundwater/delta_time &
                          + (1-c2)*flow(num_l+1)/delta_time
      lake%prog(l)%groundwater = c1 * lake%prog(l)%groundwater &
                                + c2 * flow(num_l+1)
      lake_lrunf  = x
      lake_hlrunf = x*clw*(lake%prog(l)%groundwater_T-tfreeze)
      lake_hadvec = lake_hadvec - lake_hlrunf
    enddo
  ELSE
    if (lprec_eff.ne.0. .and. flow(1).ge.0. ) then
      hlrunf_ie = lrunf_ie*hlprec_eff/lprec_eff
    else if (flow(1).lt.0. ) then
      hlrunf_ie = hlprec_eff - (flow(1)/delta_time)*clw &
                             *(lake%prog(1)%T-tfreeze)
    else
      hlrunf_ie = 0.
    endif
    hlrunf_bf = clw*sum(div*(lake%prog%T-tfreeze))
    lake_lrunf  = lrunf_sn + lrunf_ie + lrunf_bf
    lake_hlrunf = hlrunf_sn + hlrunf_bf + hlrunf_ie 
    lake_hadvec = lake_hadvec - lake_hlrunf
  ENDIF

  do l = 1, num_l
    ! ---- compute explicit melt/freeze --------------------------------------
    hcap = lake%heat_capacity_dry(l)*dz(l) &
             + clw*lake%prog(l)%wl + csw*lake%prog(l)%ws
    melt_per_deg = hcap/hlf
    if (lake%prog(l)%ws>0 .and. lake%prog(l)%T>lake%pars%tfreeze) then
      melt =  min(lake%prog(l)%ws, (lake%prog(l)%T-lake%pars%tfreeze)*melt_per_deg)
    else if (lake%prog(l)%wl>0 .and. lake%prog(l)%T<lake%pars%tfreeze) then
      melt = -min(lake%prog(l)%wl, (lake%pars%tfreeze-lake%prog(l)%T)*melt_per_deg)
    else
      melt = 0
    endif
    lake%prog(l)%wl = lake%prog(l)%wl + melt
    lake%prog(l)%ws = lake%prog(l)%ws - melt
    lake%prog(l)%T = tfreeze &
       + (hcap*(lake%prog(l)%T-tfreeze) - hlf*melt) &
                            / ( hcap + (clw-csw)*melt )
    lake_melt = lake_melt + melt / delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 5 ***** '
     do l = 1, num_l
        write(*,*) ' level=', l,&
             ' T =', lake%prog(l)%T,&
             ' Th=', (lake%prog(l)%ws +lake%prog(l)%wl)/(dens_h2o*dz(l)),&
             ' wl=', lake%prog(l)%wl,&
             ' ws=', lake%prog(l)%ws,&
             ' gw=', lake%prog(l)%groundwater
     enddo
  endif

  lake_LMASS = 0
  lake_FMASS = 0
  lake_HEAT = 0
  do l = 1, num_l
    lake_LMASS = lake_LMASS &
         +      lake%prog(l)%wl
    lake_FMASS = lake_FMASS &
         +      lake%prog(l)%ws
    lake_HEAT = lake_HEAT &
         + (lake%heat_capacity_dry(l)*dz(l) &
             + clw*lake%prog(l)%wl + csw*lake%prog(l)%ws)  &
                        * (lake%prog(l)%T-tfreeze)
  enddo
  lake_LMASS = lake_LMASS +     lake%prog(1)%groundwater
  lake_HEAT = lake_HEAT + clw*lake%prog(1)%groundwater    &
                                  *(lake%prog(1)%groundwater_T-tfreeze)
  lake_Ttop = lake%prog(1)%T
  lake_Ctop = lake%heat_capacity_dry(1)*dz(1) &
       + clw*lake%prog(1)%wl + csw*lake%prog(1)%ws

! ----------------------------------------------------------------------------
! given solution for surface energy balance, write diagnostic output.
!  

  ! ---- increment time
  time = increment_time(time, int(delta_time), 0)

  ! ---- diagnostic section
  call send_tile_data (id_temp, lake%prog%T,     diag )
  call send_tile_data (id_lwc,  lake%prog(1:num_l)%wl/dz(1:num_l), diag )
  call send_tile_data (id_swc,  lake%prog(1:num_l)%wl/dz(1:num_l), diag )
  call send_tile_data (id_ie,   lrunf_ie,        diag )
  call send_tile_data (id_sn,   lrunf_sn,        diag )
  call send_tile_data (id_bf,   lrunf_bf,        diag )
  call send_tile_data (id_hie,  hlrunf_ie,       diag )
  call send_tile_data (id_hsn,  hlrunf_sn,       diag )
  call send_tile_data (id_hbf,  hlrunf_bf,       diag )

end subroutine lake_step_2


! ============================================================================
subroutine lake_diag_init ( id_lon, id_lat )
  integer,         intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer,         intent(in) :: id_lat  ! ID of land longitude (X) axis

  ! ---- local vars
  integer :: axes(3)
  integer :: id_zhalf, id_zfull

  ! define vertical axis
  id_zhalf = diag_axis_init ( &
       'zhalf_lake', zhalf(1:num_l+1), 'meters', 'z', 'half level',  -1, set_name='lake' )
  id_zfull = diag_axis_init ( &
       'zfull_lake', zfull(1:num_l),   'meters', 'z', 'full level',  -1, set_name='lake', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/ id_lon, id_lat, id_zfull /)

  ! define diagnostic fields
  id_lwc = register_tiled_diag_field ( module_name, 'lake_liq', axes,         &
       Time, 'bulk density of liquid water', 'kg/m3', missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'lake_ice',  axes,       &
       Time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'lake_T',  axes,        &
       Time, 'temperature',            'degK',  missing_value=-100.0 )
  id_ie  = register_tiled_diag_field ( module_name, 'lake_rie',  axes(1:2),   &
       Time, 'inf exc runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_sn  = register_tiled_diag_field ( module_name, 'lake_rsn',  axes(1:2),   &
       Time, 'satn runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_bf  = register_tiled_diag_field ( module_name, 'lake_rbf',  axes(1:2),   &
       Time, 'baseflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_hie  = register_tiled_diag_field ( module_name, 'lake_hie',  axes(1:2),  &
       Time, 'heat ie runf',            'W/m2',  missing_value=-100.0 )
  id_hsn  = register_tiled_diag_field ( module_name, 'lake_hsn',  axes(1:2),  &
       Time, 'heat sn runf',            'W/m2',  missing_value=-100.0 )
  id_hbf  = register_tiled_diag_field ( module_name, 'lake_hbf',  axes(1:2),  &
       Time, 'heat bf runf',            'W/m2',  missing_value=-100.0 )

  ! define static diagnostic fields
  
end subroutine lake_diag_init

! ============================================================================
! tile existance detector: returns a logical value indicating wether component
! model tile exists or not
logical function lake_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   lake_tile_exists = associated(tile%lake)
end function lake_tile_exists


! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine lake_temp_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%prog%T
   endif
end subroutine lake_temp_ptr

subroutine lake_wl_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%prog%wl
   endif
end subroutine lake_wl_ptr

subroutine lake_ws_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%prog%ws
   endif
end subroutine lake_ws_ptr

subroutine lake_gw_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%prog%groundwater
   endif
end subroutine lake_gw_ptr

subroutine lake_gwT_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%lake)) ptr=>tile%lake%prog%groundwater_T
   endif
end subroutine lake_gwT_ptr

end module lake_mod



