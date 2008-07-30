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
     lake_data_thermodynamics, &
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
    version     = '$Id: lake.F90,v 16.0 2008/07/30 22:12:49 fms Exp $',&
    tagname     = '$Name: perth $'

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
real    :: init_temp            = 288.        ! cold-start lake T
real    :: init_w               = 1000.      ! cold-start w(l)/dz(l)
real    :: init_groundwater     =   0.        ! cold-start gw storage
logical :: use_rh_feedback      = .true.

namelist /lake_nml/ init_temp,      &
                    init_w,       &
                    init_groundwater, use_rh_feedback, cpw, clw, csw
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
type(time_type) :: time
real            :: delta_time

integer         :: num_l              ! # of water layers
real            :: dz    (max_lev)    ! thicknesses of layers
real            :: zfull (max_lev)
real            :: zhalf (max_lev+1)

! ---- diagnostic field IDs
integer :: id_lwc, id_swc, id_temp, id_ie, id_sn, id_bf, id_hie, id_hsn, id_hbf
integer :: id_evap
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
        
        if (init_temp.ge.tfreeze) then
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
  real, dimension(num_l):: aaa, ccc, thermal_cond, heat_capacity, dz_alt
  integer               :: l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  lake_T = lake%prog(1)%T
  if (use_rh_feedback) then
      vlc_sfc = max(0., lake%prog(1)%wl / (dens_h2o * dz(1)))
      vsc_sfc = max(0., lake%prog(1)%ws / (dens_h2o * dz(1)))
    else
      vlc_sfc = 1
      vsc_sfc = 0
    endif
  call lake_data_thermodynamics ( lake%pars, vlc_sfc, vsc_sfc,  &  
                                  lake_rh, &
                                  lake%heat_capacity_dry, thermal_cond )
  do l = 1, num_l
    heat_capacity(l) = lake%heat_capacity_dry(l) * dz(l) &
            + clw*lake%prog(l)%wl + csw*lake%prog(l)%ws
    dz_alt(l) = (lake%prog(l)%wl + lake%prog(l)%ws)/dens_h2o
    enddo

  lake_liq  = max(lake%prog(1)%wl, 0.)
  lake_ice  = max(lake%prog(1)%ws, 0.)
  if (lake_liq + lake_ice > 0 ) then
     lake_subl = lake_ice / (lake_liq + lake_ice)
  else
     lake_subl = 0
  endif

  if(num_l > 1) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz_alt(l+1)/thermal_cond(l+1) &
             + dz_alt(l)/thermal_cond(l)   )
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
  lake_tf = tfreeze
  
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
             dW_l, u_minus, u_plus, DPsi, lake_w_fc
  real, dimension(num_l+1) :: flow
  real, dimension(num_l  ) :: div
  real :: &
     lprec_eff, hlprec_eff, hcap, &
     melt_per_deg, melt,&
     lrunf_sn,lrunf_ie,lrunf_bf, hlrunf_sn,hlrunf_ie,hlrunf_bf
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

  ! ---- remainder of mass fluxes and associated sensible heat fluxes --------
    flow=1
    flow(1)  = snow_lprec *delta_time
    do l = 1, num_l
      flow(l+1) = 0
      dW_l(l) = flow(l) - flow(l+1)
      lake%prog(l)%wl = lake%prog(l)%wl + dW_l(l)
    enddo

  if(is_watch_point()) then
     write(*,*) ' ***** lake_step_2 checkpoint 3.3 ***** '
     do l = 1, num_l
        write(*,*) ' level=', l,&
             ' wl=', lake%prog(l)%wl,&
             'flow=', flow(l)
     enddo
  endif
  
  lake_hadvec = lake_hadvec + snow_hlprec
  hcap = lake%heat_capacity_dry(1)*dz(1) &
                     + clw*(lake%prog(1)%wl-dW_l(1)) + csw*lake%prog(1)%ws
  lake%prog(1)%T = tfreeze + (hcap*(lake%prog(1)%T-tfreeze) +  &
                                 snow_hlprec*delta_time) &
                            / ( hcap + clw*dW_l(1) )

  if(is_watch_point()) then
    write(*,*) ' ***** lake_step_2 checkpoint 3.4 ***** '
    write(*,*) ' tfreeze', tfreeze
    write(*,*) ' snow_hlprec', snow_hlprec
  endif

    lrunf_sn = 0
    lrunf_ie = 0
    lrunf_bf = 0
    hlrunf_sn = 0
    hlrunf_ie = 0
    hlrunf_bf = 0
    lake_lrunf  = lrunf_sn + lrunf_ie + lrunf_bf
    lake_hlrunf = hlrunf_sn + hlrunf_bf + hlrunf_ie 

  do l = 1, num_l
    ! ---- compute explicit melt/freeze --------------------------------------
    hcap = lake%heat_capacity_dry(l)*dz(l) &
             + clw*lake%prog(l)%wl + csw*lake%prog(l)%ws
    melt_per_deg = hcap/hlf
    if (lake%prog(l)%ws>0 .and. lake%prog(l)%T>tfreeze) then
      melt =  min(lake%prog(l)%ws, (lake%prog(l)%T-tfreeze)*melt_per_deg)
    else if (lake%prog(l)%wl>0 .and. lake%prog(l)%T<tfreeze) then
      melt = -min(lake%prog(l)%wl, (tfreeze-lake%prog(l)%T)*melt_per_deg)
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
    lake_LMASS = lake_LMASS + lake%prog(l)%wl
    lake_FMASS = lake_FMASS + lake%prog(l)%ws
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
  call send_tile_data (id_swc,  lake%prog(1:num_l)%ws/dz(1:num_l), diag )
  call send_tile_data (id_ie,   lrunf_ie,        diag )
  call send_tile_data (id_sn,   lrunf_sn,        diag )
  call send_tile_data (id_bf,   lrunf_bf,        diag )
  call send_tile_data (id_hie,  hlrunf_ie,       diag )
  call send_tile_data (id_hsn,  hlrunf_sn,       diag )
  call send_tile_data (id_hbf,  hlrunf_bf,       diag )
  call send_tile_data (id_evap, lake_levap+lake_fevap, diag )

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
  id_evap  = register_tiled_diag_field ( module_name, 'lake_evap',  axes(1:2),  &
       Time, 'lake evap',            'kg/(m2 s)',  missing_value=-100.0 )

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



