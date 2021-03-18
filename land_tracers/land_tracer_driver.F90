module land_tracer_driver_mod

#include "../shared/debug.inc"

  use constants_mod, only : rdgas,wtmair
  use time_manager_mod, only : time_type, time_type_to_real
  use fms_mod, only : lowercase, stdout, stdlog
  use field_manager_mod , only : MODEL_ATMOS, MODEL_LAND, parse
  use tracer_manager_mod, only : NO_TRACER, get_tracer_index, get_tracer_names, query_method
  use table_printer_mod

  use land_constants_mod, only : d608, kBoltz
  use land_debug_mod, only : is_watch_point
  use land_data_mod, only : lnd, log_version
  use land_tracers_mod, only : ntcana, isphum, ico2
  use land_tile_mod, only : land_tile_type, land_tile_grnd_T
  use land_tile_diag_mod, only : set_default_diag_filter, &
       register_tiled_diag_field, send_tile_data

  use cana_tile_mod, only : canopy_air_mass_for_tracers
  use vegn_data_mod, only : spdata
  use vegn_tile_mod, only : vegn_tile_LAI, vegn_tile_SAI  
  use vegn_cohort_mod, only : get_vegn_wet_frac

  use soil_tile_mod, only : num_l, soil_theta
  use snow_mod,      only : snow_get_depth_area

  ! import interfaces from non-generic tracer modules, e.g.:
  use land_dust_mod, only : land_dust_init, land_dust_end, update_land_dust

  implicit none
  private

  ! ==== public interfaces =====================================================
  public :: land_tracer_driver_init
  public :: land_tracer_driver_end
  public :: update_cana_tracers
  ! ==== end of public interfaces ==============================================

  ! ---- module constants ------------------------------------------------------
  character(len=*), parameter :: module_name = 'land_tracer_driver_mod'
#include "../shared/version_variable.inc"
  character(len=*), parameter :: diag_name   = 'land_tracers'

  real :: diffusivity_h2o = 0.282e-4 ! diffusivity of water vapor m2/s,
  ! Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems (2nd ed.).
  ! New York: Cambridge University Press. ISBN 0-521-45078-0.


  !for dry deposition
  real, parameter  :: mw_air = WTMAIR/1000
  
  real   :: desert_biomass    = 0.25 !kg/m2                                                                                                                   
  real   :: ws_min            = 1. !1mm  (threshold to decide whether lake is frozen)                                                                         
  real   :: theta_wetland_thr = 0.9 !above this call it wetland                                                                                               

  real :: r_gs_lake         = 20.    
  real :: r_gs_swamp        = 50.    
  real :: r_gs_desert       = 700.
  real :: r_gs_glac         = 70.    !this will get overwritten by r_snow                                                                                     
  real :: r_snows           = 70.

  real :: r_go_lake         = 500.   
  real :: r_go_swamp        = 500.   
  real :: r_go_desert       = 500.   
  real :: r_go_glac         = 2000.
  real :: r_snowo           = 2000.

  real :: A_aer_lake         = -999    
  real :: A_aer_swamp        = 10.e-3
  real :: A_aer_desert       = -999 

  real :: gamma_aer_lake         = 0.50    
  real :: gamma_aer_swamp        = 0.54       
  real :: gamma_aer_desert       = 0.54
  real :: gamma_aer_glac         = 0.54

  real :: alpha_aer_lake         = 100.    
  real :: alpha_aer_swamp        = 50.        
  real :: alpha_aer_desert       = 50.
  real :: alpha_aer_snow         = 50.

  ! ---- data types -----------------------------------------------------------
  type :: tracer_data_type
     character(32) :: name = ''  ! tracer name                                                                                
     integer :: tr_atm  = NO_TRACER ! index of this tracer in atmos tracer array                                              
     logical :: is_generic    = .TRUE. ! flag of generic tracer; initialization of non-generic tracers should turn it to FALSE
     logical :: do_deposition = .TRUE. ! if true, generic dry deposition is used                                              
     ! dry deposition parameters. The default values are set as O3 parameters from (Wesely, 1989)                             
     real    :: reactivity    = 1.0    ! normalized reactivity factor                                                         
     real    :: alpha         = -1     ! scaling factor relative to SO2                                                       
     real    :: r_mx          = -999    ! no resistance                                                                       
     real    :: mw            = -9999.9  ! kg/mol                                                                             
     logical :: coldTc        = .false. ! cold t increases resistance                                                         
     real    :: diff_ratio    = 1.6    ! ratio of water vapor molecular diffusivity in the air to that of the tracer, unitless

     !for aerosol
     real           :: radius = -999., density = 2e3                                                                                                         
     logical        :: is_aerosol = .false.
     
     real           :: conv_flux

     integer :: & ! diag field IDs
          id_emis,      id_ddep,  &
          id_flux_atm,  id_dfdtr, &
          id_con_v_lam, id_con_g_lam, &
          id_con_v,     id_con_g, &
          id_conc
  end type tracer_data_type

  ! ---- private module variables ----------------------------------------------
  logical :: module_is_initialized = .FALSE.
  real, save :: dt ! fast time step, s
  type(tracer_data_type), allocatable :: trdata(:)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! ============================================================================
  subroutine land_tracer_driver_init(id_ug)
    integer,intent(in) :: id_ug !<Unstructured axis id.

    integer :: tr ! tracer index
    real    :: value ! temporary storage for parsing input
    character(32)  :: name, units, funits ! name and units of the tracer and flux
    character(128) :: longname ! long name of the tracer
    character(32)  :: method
    character(1024) :: parameters
    type(table_printer_type) :: table

    ! write the version and tag name to the logfile
    call log_version(version, module_name, &
         __FILE__)

    ! calculate time step
    dt  = time_type_to_real(lnd%dt_fast) ! store in a module variable for convenience

    ! allocate storage for tracer informations
    allocate(trdata(ntcana))
    trdata(isphum)%is_generic = .FALSE.
    trdata(ico2)%is_generic   = .FALSE.

    ! initialize non-generic tracers, e.g.:
    call land_dust_init(id_ug, trdata(:)%is_generic)
    ! NOTE that (1) the non-generic tracer init must skip all non-generic tracers
    ! that have already been initialized (in case there is a conflict), and
    ! (2) it must set trdata(:)%is_generic to FALSE for the tracers it claims

    ! initialize generic tracer parameters
    do tr = 1, ntcana
       if (.not.trdata(tr)%is_generic) cycle ! skip all non-generic tracers

       call get_tracer_names(MODEL_LAND, tr, trdata(tr)%name)
       trdata(tr)%tr_atm = get_tracer_index (MODEL_ATMOS, trdata(tr)%name)

       ! set up deposition flag
       trdata(tr)%do_deposition = .FALSE.
       method = ''; parameters = ''
       if (trdata(tr)%tr_atm >0) then
          if (query_method('dry_deposition', MODEL_ATMOS, trdata(tr)%tr_atm, method)) then
             trdata(tr)%do_deposition = (index(lowercase(method),'land:lm3')>0)
          endif
       endif
       ! set up deposition parameters
       if(query_method('dry_deposition', MODEL_LAND, tr, method, parameters)) then
          !ratio of tracer diffusivity to h2o diffusivity
          if ( parse(parameters, 'diff_ratio',  value) > 0 ) then
             trdata(tr)%diff_ratio  = value
          else
             if (trdata(tr)%mw.gt.0.) then
                !graham law                                                                    
                trdata(tr)%diff_ratio  = sqrt(trdata(tr)%mw/18e-3)
             end if
          end if
          if ( parse(parameters, 'reactivity',  value) > 0 ) trdata(tr)%reactivity  = value
          if ( parse(parameters, 'alpha',       value) > 0 ) trdata(tr)%alpha       = value        
       endif
    enddo

    call init_with_headers(table, trdata(:)%name)
    call add_row(table, 'do_deposition',    trdata(:)%do_deposition)
    call add_row(table, 'diff_ratio',       trdata(:)%diff_ratio)
    call add_row(table, 'alpha',            trdata(:)%alpha)    
    call add_row(table, 'reactivity',       trdata(:)%reactivity)

    call print(table,stdlog())
    call print(table,stdout())

    ! register diag fields for generic tracers
    call set_default_diag_filter('land')
    do tr = 1, ntcana
       call get_tracer_names(MODEL_LAND, tr, name, longname, units)
       call flux_units(units,funits,trdata(tr)%conv_flux)
       
       trdata(tr)%id_flux_atm = &
            register_tiled_diag_field(diag_name, trim(name)//'_flux_atm', &
            (/id_ug/),  lnd%time, trim(name)//' flux to the atmosphere', &
            trim(funits), missing_value=-1.0)
       ! TODO: verify units of dfdtr
       trdata(tr)%id_dfdtr = &
            register_tiled_diag_field(diag_name, trim(name)//'_dfdtr', &
            (/id_ug/),  lnd%time,'derivative of '//trim(name)//' flux to the atmosphere', &
            trim(funits), missing_value=-1.0)
       if (trdata(tr)%do_deposition) then
          ! TODO: initialize parameters of generic dry deposition here

          trdata(tr)%id_ddep = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep', &
               (/id_ug/),  lnd%time, trim(name)//' dry deposition', 'kg/(m2 s)', &
               missing_value=-1.0)
          trdata(tr)%id_con_v_lam = &
               register_tiled_diag_field(diag_name, trim(name)//'_con_v_lam', &
               (/id_ug/),  lnd%time, 'quasi-laminar conductance between canopy and canopy air for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_con_g_lam = &
               register_tiled_diag_field(diag_name, trim(name)//'_con_g_lam', &
               (/id_ug/),  lnd%time, 'quasi-laminar conductance between ground and canopy air for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_con_v = &
               register_tiled_diag_field(diag_name, trim(name)//'_con_v', &
               (/id_ug/),  lnd%time, 'total conductance between canopy and canopy air for'//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_con_g = &
               register_tiled_diag_field(diag_name, trim(name)//'_con_g', &
               (/id_ug/),  lnd%time, 'total conductance between ground and canopy air for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_conc = &
               register_tiled_diag_field(diag_name, trim(name), &
               (/id_ug/),  lnd%time, 'concentration or '//trim(name)//' in canopy air', &
               units, missing_value=-1.0)
       endif
    enddo
    module_is_initialized = .TRUE.
  end subroutine land_tracer_driver_init


  ! ============================================================================
  subroutine land_tracer_driver_end()
    ! call finalizers for all non-generic tracers, e.g.:
    call land_dust_end()

    deallocate(trdata)
    module_is_initialized = .FALSE.
  end subroutine land_tracer_driver_end

  ! ============================================================================
  ! updates concentration of tracers in the canopy air, taking into account dry
  ! deposition and exchange with the atmosphere
  subroutine update_cana_tracers(tile, l, tr_flux, dfdtr, &
       precip_l, precip_s, pressure, ustar, con_g, con_v_v, con_v_stem, stomatal_cond, r_bl_h2o, con_atm )
    type(land_tile_type), intent(inout) :: tile
    integer :: l ! grid cell indices (global)
    real, intent(in) :: tr_flux(:) ! fluxes of tracers
    real, intent(in) :: dfdtr  (:) ! derivatives of tracer fluxes w.r.t. concentrations
    real, intent(in) :: pressure   ! atmospheric pressure, N/m2
    real, intent(in) :: precip_l   ! liquid precipitation, kg/(m2 s)
    real, intent(in) :: precip_s   ! solid precipitation, kg/(m2 s)
    real, intent(in) :: ustar ! friction velocity, m/s
    real, intent(in) :: con_g ! aerodynamic conductance between canopy air and ground for tracers
    real, intent(in) :: con_v_v(:) ! aerodynamic conductance between canopy air and canopy
    ! for tracers, by cohort, per unit area occupied by cohort
    real, intent(in) :: con_v_stem(:)
    real, intent(in) :: stomatal_cond(:) ! integral stomatal conductance of each cohort
    ! canopy (that is, multiplied by LAI), for water vapor, m/s
    real, intent(in) :: r_bl_h2o
    real, intent(in) :: con_atm !Ra from flux exchange

    integer :: tr      ! tracer index
    integer :: k       ! cohort index
    real    :: rho     ! density of canopy air
    real    :: dq      ! canopy air tracer tendency per time step
    real    :: con_v   ! laminar conductance to leaves
    real    :: con_st  ! total stomatal conductance, scaled by dry leaf area, m/s
    real    :: con_stem ! laminar conductance to stem
    real    :: con_mx  ! "mesophyll conductance", m/s
    real    :: con_mx_st !mesophyll+stomatal
    real    :: con_cu  ! total cuticular conductance, including dry and wet areas, m/s
    real    :: con_gr  ! "ground conductance", m/s
    real    :: cv0, cv1, cv2, cg0 ! intermediate values for total conductance calculations, m/s
    real    :: cv, cg  ! total conductances for vegetation and ground surface, m/s
    real    :: f_atm   ! flux of the tracer to the atmosphere, kg/(m2 s)
    real    :: ddep    ! dry deposition of the tracer, kg/(m2 s)
    real    :: LAI     ! leaf area index, m2/m2
    real    :: SAI     ! stem area index, m2/m2
    real    :: kvis,dvis
    real    :: emis(ntcana) ! tracer sources
    real    ::  ft, & ! fraction of canopy not covered by intercepted water/snow
         fw, & ! fraction of canopy covered by intercepted water
         fs    ! fraction of canopy covered by intercepted snow

    !fraction of ground that is frozen, wet, dry
    real    :: gfrac_frz,gfrac_wet,gfrac_dry

    real    :: con_cu_dry, con_cu_wet, con_cu_frz, con_bk                                           
    real    :: con_gr_dry, con_gr_wet, con_gr_frz
    real    :: r_gs, frac_desert
    real    :: con_v_v_tr, con_v_stem_tr, con_bl_tr, con_g_tr

    if (is_watch_point()) then
       write(*,*) 'update_cana_tracers input'
       __DEBUG1__(tr_flux)
       __DEBUG1__(dfdtr)
       __DEBUG1__(pressure)
       __DEBUG2__(precip_l, precip_s)
       __DEBUG1__(ustar)
       write(*,*) 'end of update_cana_tracers input'
    endif

    ! update non-generic tracers, e.g.:
    call update_land_dust(tile, l, tr_flux, dfdtr, &
         precip_l, precip_s, pressure, ustar, con_g, con_v_v )
    ! wind10 is not passed to this subroutine yet

    ! update generic tracers
    ! calculate tracers sources
    emis(:) = 0.0
    ! TODO: add non-zero sources for generic tracers

    if (associated(tile%vegn)) then
       LAI = vegn_tile_LAI(tile%vegn)
       SAI = vegn_tile_SAI(tile%vegn)   
    else
       LAI = 0.0
       SAI = 0.0
    endif

    call get_tile_property(tile,gfrac_frz,gfrac_wet,frac_desert)
    gfrac_dry = max(1.-gfrac_wet-gfrac_frz,0.)  

    ! loop for generic tracers only
    do tr = 1, ntcana

       if (.not.trdata(tr)%is_generic) cycle

       if (trdata(tr)%do_deposition) then

          if (.not. trdata(tr)%is_aerosol) then

             !conductance to the vegetation
             if (associated(tile%vegn)) then

                !get fraction of ground covered with snow

                do k = 1, tile%vegn%n_cohorts
                   associate(c=>tile%vegn%cohorts(k),sp=>spdata(tile%vegn%cohorts(k)%species))

                     call get_vegn_wet_frac ( c, fw=fw, fs=fs ); ft = 1-fw-fs
                     !need to implement rh dependence for con_cu_dry (f1p)
                     !need to implement rt dependence
                     !need to implement snow modulation
                     con_cu_dry  = get_conductance_tracer(trdata(tr),sp%r_cus,sp%r_cuo)
                     con_cu_wet  = get_conductance_tracer(trdata(tr),sp%r_cus_wet,sp%r_cuo_wet)
                     con_cu_frz  = get_conductance_tracer(trdata(tr),r_snows,r_snowo)                   

                     !here we use the bulk leaf property for the cohort. This is different from the LM3 implementation.

                     con_cu = c%lai * &
                          ( ft * con_cu_dry &
                          + fw * con_cu_wet &
                          + fs * con_cu_frz )

                     con_bk = c%sai * get_conductance_tracer(trdata(tr),sp%r_bks,sp%r_bko)

                     if (frac_desert.gt.0.) then
                        r_gs = r_gs_desert*frac_desert+sp%r_gs*(1.-frac_desert)
                     end if

                     con_gr_dry  = get_conductance_tracer(trdata(tr),r_gs,sp%r_go)
                     con_gr_wet  = get_conductance_tracer(trdata(tr),r_gs_swamp,r_go_swamp)
                     con_gr_frz  = get_conductance_tracer(trdata(tr),r_snows,r_snowo)

                     !note that for the ground we are calculating the tile average
                     con_gr      = con_gr +               &
                          c%layerfrac  *                  &
                          ( gfrac_dry  * con_gr_dry       &
                          + gfrac_wet  * con_gr_wet       &
                          + gfrac_frz  * con_gr_frz )

                     con_st      = stomatal_cond(k) / trdata(tr)%diff_ratio

                     if (trdata(tr)%r_mx > 0) then
                        con_mx    = 1./trdata(tr)%r_mx
                     else
                        con_mx    = get_conductance_tracer(trdata(tr),1.,100.)
                     end if

                     ! combined mesophyll + stomatal conductance
                     con_mx_st = con_mx*con_st/(con_mx+con_st)

                     ! calculate contribution of this cohort to the overall vegetation conductance
                     !con_v_v and con_bk are for H2O, we need to scale by (Di/Dw)**(2./3.)

                     con_v_v_tr = con_v_v(k)*1./trdata(tr)%diff_ratio**(2./3.)
                     con_v_stem_tr  = con_v_stem(k)*1./trdata(tr)%diff_ratio**(2./3.)

                     cv = cv &
                          + c%layerfrac*(con_v_v_tr*(con_mx_st+con_cu)/(con_v_v_tr+con_mx_st+con_cu) &
                          + con_v_stem_tr*con_bk/(con_bk+con_v_stem_tr))

                   end associate

                end do

             elseif (associated(tile%lake)) then
                con_gr = trdata(tr)%alpha/r_gs_lake + trdata(tr)%reactivity/r_go_lake
                !f1p need to deal with frozen lake
             elseif (associated(tile%glac)) then
                con_gr = trdata(tr)%alpha/r_gs_glac + trdata(tr)%reactivity/r_go_glac
             endif

             con_bl_tr  = 1./r_bl_h2o   * 1./trdata(tr)%diff_ratio**(2./3.)
             con_g_tr   = con_g*con_bl_tr/(con_g+con_bl_tr)

             cg = con_gr*con_g_tr/(con_g_tr+con_gr)

          endif
       else
          cv=0.
          cg=0.
       end if


       rho = pressure/(rdgas*tile%cana%T *(1+d608*tile%cana%tr(isphum)))

       ! update tracer concentration -- this needs to be done even when do_deposition
       ! is FALSE, because source might be non-zero
       dq = (-tr_flux(tr)-rho*(cv+cg)*tile%cana%tr(tr)+emis(tr)) &
            / (canopy_air_mass_for_tracers/dt + dfdtr(tr) + rho*(cv+cg))
       if (is_watch_point()) then
          write(*,*) 'update_cana_tracers : ', trim(trdata(tr)%name)
          __DEBUG4__(tile%cana%tr(tr), tr_flux(tr), dfdtr(tr), emis(tr))
          __DEBUG3__(rho,cv,cg)
          __DEBUG2__(dq,tile%cana%tr(tr) + dq)
       endif
       tile%cana%tr(tr) = tile%cana%tr(tr) + dq
       ! ---- final values of the fluxes, for diagnostics
       ddep  = rho*(cv+cg)*tile%cana%tr(tr)
       f_atm = tr_flux(tr)+dfdtr(tr)*dq
       ! ---- diagnostic section
       !call send_tile_data(trdata(tr)%id_con_v_lam,  con_v_lam,  tile%diag)
       !call send_tile_data(trdata(tr)%id_con_g_lam,  con_g_lam,  tile%diag)
       call send_tile_data(trdata(tr)%id_con_v,      cv,         tile%diag)
       call send_tile_data(trdata(tr)%id_con_g,      cg,         tile%diag)
       call send_tile_data(trdata(tr)%id_emis,       emis(tr),   tile%diag)
       call send_tile_data(trdata(tr)%id_ddep,       ddep*trdata(tr)%conv_flux,       tile%diag)
       call send_tile_data(trdata(tr)%id_flux_atm,   ddep*trdata(tr)%conv_flux,       tile%diag)
    enddo
    ! send concentrations for all tracers, generic or not
    do tr = 1, ntcana
       call send_tile_data(trdata(tr)%id_conc,       tile%cana%tr(tr), tile%diag)
       call send_tile_data(trdata(tr)%id_dfdtr,      dfdtr(tr),  tile%diag)
    enddo

  end subroutine update_cana_tracers

  subroutine get_tile_property(tile,gfrac_frz,gfrac_wet,frac_desert)

    type(land_tile_type), intent(in)   :: tile

    real,    intent(out)             :: gfrac_frz,gfrac_wet,frac_desert

    real                             :: depth_to_wt_2b, biomass
    real,dimension(num_l)            :: theta
    integer                          :: k
    real                             :: snow_depth, snow_area


    frac_desert = 0.
    gfrac_frz   = 0.
    gfrac_wet   = 0.

    if (associated(tile%snow)) then
       call snow_get_depth_area ( tile%snow, snow_depth, snow_area )
       gfrac_frz = snow_area
    end if

    if (associated(tile%glac)) then
       gfrac_frz = 1.
    end if

    if (associated(tile%lake)) then
       gfrac_wet = 1.
       if ( tile%lake%ws(1) .gt. ws_min ) then
          gfrac_frz  = 1.
          gfrac_wet  = 0.
       end if
    end if

    if (associated(tile%vegn)) then
       theta = soil_theta(tile%soil) !top layer
       if ( theta(1) .gt. theta_wetland_thr ) then
          gfrac_wet = 1.
       else
          gfrac_wet = 0.
       end if

       biomass = 0.
       associate(cc=>tile%vegn%cohorts)      
         do k = 1, tile%vegn%n_cohorts
            biomass = &
                 cc(k)%bl  + cc(k)%blv + &
                 cc(k)%br  + cc(k)%bwood + &
                 cc(k)%bsw + cc(k)%bseed + &
                 cc(k)%nsc 
         end do
       end associate


       if ( biomass .lt. desert_biomass ) then
          frac_desert = min(max(1.-biomass/desert_biomass,0.),1.)
       end if
    end if

  end subroutine get_tile_property

  elemental real function snow_scale(T) result(s)

    !erisman (1994) showed that resistance increases as temperature decreases from 70 to 500
    real, intent(in) :: T
    real, parameter  :: max_scale = 500./70.

    s = max(min(max_scale,275.15-T),1.)

  end function snow_scale

  elemental real function calc_rt(tk) result(r_t)
    !scaling factor for resistance at cold temperature
    real, intent(in) :: tk

    r_t   = max(min(2.,exp(0.2*(-1-(tk-273.15)))),1.)
  end function calc_rt

  elemental real function get_conductance_tracer(trdata,r_s,r_o) result(con)
    real, intent(in)                   :: r_s, r_o
    type(tracer_data_type), intent(in) :: trdata

    con = trdata%alpha/r_s + trdata%reactivity/r_o

  end function get_conductance_tracer



  ! ============================================================================
  subroutine flux_units(tracer_units,units,conv)
    character(*), intent(in)   :: tracer_units
    character(32), intent(out) :: units
    real, intent(out)          :: conv

    select case (trim(lowercase(tracer_units)))
    case ('mmr')
       units = 'kg/m2/s'
       conv  = 1.
    case ('kg/kg')
       units = 'kg/m2/s'
       conv  = 1.
    case ('vmr')
       units = 'mole/m2/s'
       conv  = 1./mw_air
    case ('mol/mol')
       units = 'mole/m2/s'
       conv  = 1./mw_air
    case ('mole/mole')
       units = 'mole/m2/s'
       conv = 1./mw_air
    case default
       units = trim(tracer_units)//' kg/(m2 s)'
       conv = 1.
    end select
  end subroutine flux_units

end module land_tracer_driver_mod
