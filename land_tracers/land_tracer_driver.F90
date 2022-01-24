module land_tracer_driver_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file
#endif

  use ieee_arithmetic

  use constants_mod, only : rdgas,wtmair,grav,pi,pstd_mks,avogno,DENS_H2O,epsln
  use time_manager_mod, only : time_type, time_type_to_real
  use fms_mod, only : lowercase, stdout, stdlog, mpp_pe, mpp_root_pe, file_exist, close_file, &
       check_nml_error
  use field_manager_mod , only : MODEL_ATMOS, MODEL_LAND, parse
  use tracer_manager_mod, only : NO_TRACER, get_tracer_index, get_tracer_names, query_method
  use table_printer_mod

  use land_constants_mod, only : d608,kin_visc_air,dyn_visc_air
  use land_debug_mod, only : is_watch_point, check_var_range
  use land_data_mod, only : lnd, log_version
  use land_tracers_mod, only : ntcana, isphum, ico2
  use land_tile_mod, only : land_tile_type, land_tile_grnd_T, loop_over_tiles, &
       first_elmt, land_tile_enum_type, land_tile_map
  use land_tile_diag_mod, only : set_default_diag_filter, &
       register_tiled_diag_field, send_tile_data

  use cana_tile_mod, only : canopy_air_mass_for_tracers
  use vegn_data_mod, only : spdata
  use vegn_cohort_mod, only : get_vegn_wet_frac
  use vegn_tile_mod, only : vegn_tile_fw_fs

  use soil_tile_mod, only : num_l, soil_theta, soil_ice_porosity, zhalf, n_dim_soil_types
  use snow_mod,      only : snow_get_depth_area

  use sat_vapor_pres_mod, only: compute_qs
  
  ! import interfaces from non-generic tracer modules, e.g.:
  use land_dust_mod, only : land_dust_init, land_dust_end, update_land_dust

  implicit none
  private

  ! ==== public interfaces =====================================================
  public :: land_tracer_driver_init
  public :: land_tracer_driver_end
  public :: update_cana_tracers
  ! ==== end of public interfaces ==============================================

  !maximum increase of snow resistance with lower temperature
  real :: max_scale_snow_T=500./70.
  !maximum increase associated with cold T
  real :: max_scale_cold_T=2.
  !maximum increase associated with low biomass
  real :: max_scale_desert=2.5
  real :: cg_aer_frz = 0.03e-2   !from fisher (2011) vd(aer) = 0.03 cm/s, conductance of aerosol over snow.

  real   :: desert_biomass    = 0.25 !kg/m2
  real   :: ws_min            = 1. !1mm  (threshold to decide whether lake is frozen)
  real   :: theta_wetland_thr = 0.9 !above this call it wetland

  real :: r_gs_lake         = 20.     !ground resistance, lake, SO2 (s/m)
  real :: r_gs_wet          = 100.    !ground resistance, wet surface, SO2 (s/m)
  real :: r_gs_dry          = 200.    !ground reistance, dry surface, SO2 (s/m)
  real :: r_snows           = 70.     !ground resitance, snow surface, SO2 (s/m)

  real :: r_go_lake         = 500.    !ground resistance, lake, O3 (s/m)
  real :: r_go_wet          = 500.    !ground resistance, wet surface, O3 (s/m)
  real :: r_go_dry          = 200.    !ground resistance, dry surface, O3 (s/m)
  real :: r_snowo           = 7000.   !ground resistance, snow surface, O3 (s/m). Default value updated to match Clifton (2020)

  real :: A_aer_lake         = -999   !characteristic aerosol radius for deposition, lake, m
  real :: A_aer_swamp        = 10.e-3 !characteristic aerosol radius for deposition, wet, m

  !Sc power for aerosol deposition, unitless
  real :: gamma_aer_lake         = 0.50
  real :: gamma_aer_swamp        = 0.54
  real :: gamma_aer_desert       = 0.54
  real :: gamma_aer_frz          = 0.54

  !alpha is used for E_im in aerosol deposition
  real :: alpha_aer_lake         = 100.
  real :: alpha_aer_swamp        = 50.
  real :: alpha_aer_desert       = 50.
  real :: alpha_aer_frz          = 50.

  !in LM4p1, the following soil categories are used
  !'clay(heavy)', 'silty clay', 'clay (light)', 'silty clay loam', 'clay loam', 'silt', 'silt loam', 'sandy clay', 'loam', 'sandy clay loam', 'sandy loam', 'loamy sand', 'sand', 'undefined'

  real :: h2_b(n_dim_soil_types)     = (/  4.05, 10.40, 11.40, 7.75, 8.52, 5.3, 5.3, 10.4, 5.4, 7.12, 4.9, 4.38, 4.05,  4.05/),    &
          h2_st(n_dim_soil_types)    = (/  0.11, 0.480, 0.530, 0.37, 0.45, 0.31, 0.31, 0.44, 0.24, 0.32, 0.18, 0.11, 0.11,  0.11/),    &
          h2_N(n_dim_soil_types)     = (/0.1452, 0.131, 0.118,  0.152, 0.163, 0.215, 0.215, 0.112, 0.180, 0.147, 0.170, 0.145, 0.157,0.157/),    &
          h2_betab(n_dim_soil_types) = (/5.3693, 2.293, 2.232, 2.541, 2.152, 2.270, 2.270, 2.641, 2.942, 2.85, 3.596, 5.3693, 4.977, 4.977/),        &
          h2_a   = 1.4

  real :: c_snow=0.025, c_dry=0.1, c_wet=0.9 !strength of R increase with decreasing T under <5C (Clifton 2020)


  namelist /land_tracer_nml/ &
       max_scale_snow_T,  max_scale_cold_T, max_scale_desert, cg_aer_frz, &
       desert_biomass,ws_min,theta_wetland_thr, &
       r_gs_lake,r_gs_wet,r_gs_dry,r_snows,     &
       r_go_lake,r_go_wet,r_go_dry,r_snowo,     &
       A_aer_lake,A_aer_swamp,                  &
       gamma_aer_lake,gamma_aer_swamp,gamma_aer_desert,gamma_aer_frz, &
       alpha_aer_lake,alpha_aer_swamp,alpha_aer_desert,alpha_aer_frz, &
       h2_b, h2_st, h2_N, h2_a, h2_betab, &
       c_snow, c_dry, c_wet


  ! ---- module constants ------------------------------------------------------
  character(len=*), parameter :: module_name = 'land_tracer_driver_mod'
#include "../shared/version_variable.inc"
  character(len=*), parameter :: diag_name   = 'land_tracers'

  !for dry deposition
  real, parameter  :: mw_air = WTMAIR/1000
  real, parameter  :: rgas  = (rdgas*mw_air)
  real, parameter  :: kb     = rgas/avogno

  !for wetness diag
  integer, parameter :: nwet_diag = 11
  real      :: wet_diag_thr(nwet_diag) = (/ 0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95/)
  character(len=5) :: wet_str(nwet_diag)
  data wet_str/'wet05','wet10','wet20','wet30','wet40','wet50','wet60','wet70','wet80','wet90','wet95'/
  

  ! ---- data types -----------------------------------------------------------
  type :: tracer_data_type
     character(32) :: name = ''          ! tracer name
     integer :: tr_atm  = NO_TRACER      ! index of this tracer in atmos tracer array
     logical :: is_generic    = .TRUE.   ! flag of generic tracer; initialization of non-generic tracers should turn it to FALSE
     logical :: do_deposition = .FALSE.  ! if true, generic dry deposition is used
     ! dry deposition parameters. The default values are set as O3 parameters from (Wesely, 1989)
     real    :: reactivity    = 1.0      ! normalized reactivity factor
     real    :: alpha         = -1       ! scaling factor relative to SO2
     real    :: r_mx          = 1e-5     ! negligible resistance
     real    :: mw            = -9999.9  ! kg/mol
     real    :: diff_ratio    = 1.6      ! ratio of water vapor molecular diffusivity in the air to that of the tracer, unitless
     real    :: scale_stom    = 1.       ! additional
     !for aerosol
     real           :: radius = 0.25e-6, rho = 1500.
     logical        :: is_aerosol = .false.

     real           :: conv_flux

     logical        :: is_h2=.FALSE.

     real           :: gamma=1., km=0.03, depth=0.1, con_gr = -999  !for h2

     integer :: & ! diag field IDs
          id_emis,      id_ddep,  &
          id_flux_atm,  id_dfdtr, &
          id_con_v,     id_con_g, &
          id_con_mx_st, id_con_cu, id_con_stem, id_con_gr, &
          id_conc,      id_tcond, id_tcond_wet(nwet_diag), &
          id_econ_v,    id_econ_g, &
          id_ddep_v,    id_ddep_g, &
          id_econ_g_dry,id_econ_g_wet, id_econ_g_frz, &
          id_ddep_g_dry,id_ddep_g_wet, id_ddep_g_frz, &
          id_econ_stem, id_econ_stom, id_econ_cu, &
          id_ddep_stem, id_ddep_stom, id_ddep_cu, &
          id_econ_cu_dry,id_econ_cu_wet, id_econ_cu_frz, &
          id_ddep_cu_dry,id_ddep_cu_wet, id_ddep_cu_frz, &
          id_con_v_v, id_con_v_stem, id_con_v_g, &
          id_Eb, id_Ein, id_Eim
  end type tracer_data_type

  integer :: id_fw_avg, id_fs_avg, id_fd_avg
  integer :: id_fw_wet(nwet_diag)
  integer :: id_con_atm
  integer :: id_gfrac_dry, id_gfrac_wet, id_gfrac_frz, id_frac_desert
  integer :: id_h2_fm, id_h2_ft, id_h2_sdiff
  integer :: id_h2_n, id_h2_st, id_h2_betab, id_h2_b

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
    character(32)  :: value_str
    character(32)  :: name, units, funits ! name and units of the tracer and flux
    character(128) :: longname ! long name of the tracer
    character(32)  :: method
    character(1024) :: parameters
    type(table_printer_type) :: table

    integer :: unit         ! unit for namelist i/o
    integer :: io           ! i/o status for the namelist
    integer :: ierr         ! error code, returned by i/o routines

    integer :: soil_tag, iw
    type(land_tile_enum_type)     :: ce   ! tile list enumerator
    type(land_tile_type), pointer :: tile ! pointer to current tile
    
    
    ! write the version and tag name to the logfile
    call log_version(version, module_name, &
         __FILE__)

    !read namelist
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=land_tracer_nml, iostat=io)
    ierr = check_nml_error(io, 'land_tracer_nml')
#else
    if (file_exist('input.nml')) then
       unit = open_namelist_file()
       ierr = 1;
       do while (ierr /= 0)
          read (unit, nml=land_tracer_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'land_tracer_nml')
       enddo
10     continue
       call close_file (unit)
    endif
#endif
    if (mpp_pe() == mpp_root_pe()) then
       unit=stdlog()
       write(unit, nml=land_tracer_nml)
    endif


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
       call get_tracer_names(MODEL_LAND, tr, trdata(tr)%name)
       if (.not.trdata(tr)%is_generic) cycle ! skip all non-generic tracers
       trdata(tr)%tr_atm = get_tracer_index (MODEL_ATMOS, trdata(tr)%name)

       if (lowercase(trim(trdata(tr)%name)).eq."h2") trdata(tr)%is_h2=.true.


       ! set up deposition flag
       trdata(tr)%do_deposition = .FALSE.
       method = ''; parameters = ''
       if (trdata(tr)%tr_atm >0) then
          if (query_method('dry_deposition', MODEL_ATMOS, trdata(tr)%tr_atm, method,parameters)) then
             if (parse(parameters,'land', value_str)>0) then
                if (trim(value_str)=="interactive") trdata(tr)%do_deposition=.TRUE.
                if (trim(value_str)=="interactive_aer") then
                   trdata(tr)%do_deposition = .TRUE.
                   trdata(tr)%is_aerosol    = .TRUE.
                end if
             end if
          endif
       endif
       ! set up deposition parameters
       if(query_method('dry_deposition', MODEL_LAND, tr, method, parameters)) then

          if ( parse(parameters, 'reactivity',  value) > 0 ) trdata(tr)%reactivity  = value
          if ( parse(parameters, 'alpha',       value) > 0 ) trdata(tr)%alpha       = value
          if ( parse(parameters, 'mw',  value) > 0 )         trdata(tr)%mw          = value

          if ( parse(parameters, 'gamma',  value) > 0 )    trdata(tr)%gamma  = value
          if ( parse(parameters, 'km',  value) > 0 )       trdata(tr)%km     = value
          if ( parse(parameters, 'depth',  value) > 0 )    trdata(tr)%depth  = value
          if ( parse(parameters, 'con_gr',  value) > 0 )   trdata(tr)%con_gr = value

          !ratio of h2o to tracer diffusivity
          if ( parse(parameters, 'diff_ratio',  value) > 0 ) then
             trdata(tr)%diff_ratio  = value
          else
             if (trdata(tr)%mw.gt.0.) then
                !graham law
                trdata(tr)%diff_ratio  = sqrt(trdata(tr)%mw/18e-3)
             end if
          end if

          if (parse(parameters, 'scale_stom', value) > 0) then
             trdata(tr)%scale_stom = value
          else
             trdata(tr)%scale_stom = 1./trdata(tr)%diff_ratio
          end if
          if ( parse(parameters, 'r_mx',  value) > 0 )   trdata(tr)%r_mx  = max(value,epsln)
          if ( parse(parameters, 'radius', value) > 0 )  trdata(tr)%radius  = value
          if ( parse(parameters, 'rho', value) > 0 )     trdata(tr)%rho = value


       endif
    enddo

    call init_with_headers(table, trdata(:)%name)
    call add_row(table, 'do_deposition',    trdata(:)%do_deposition)
    call add_row(table, 'diff_ratio',       trdata(:)%diff_ratio)
    call add_row(table, 'alpha',            trdata(:)%alpha)
    call add_row(table, 'reactivity',       trdata(:)%reactivity)
    call add_row(table, 'scale_stom',       trdata(:)%scale_stom)
    call add_row(table, 'is_aerosol',       trdata(:)%is_aerosol)
    call add_row(table, 'radius',           trdata(:)%radius)
    call add_row(table, 'rho',              trdata(:)%rho)

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
               (/id_ug/),  lnd%time, trim(name)//' dry deposition', trim(funits), &
               missing_value=-1.0)

          trdata(tr)%id_con_v = &
               register_tiled_diag_field(diag_name, trim(name)//'_con_v', &
               (/id_ug/),  lnd%time, 'total conductance between canopy and canopy air for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_con_v_v = &
               register_tiled_diag_field(diag_name, trim(name)//'_con_v_v', &
               (/id_ug/),  lnd%time, 'conductance between canopy and leaf for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_con_v_stem = &
               register_tiled_diag_field(diag_name, trim(name)//'_con_v_stem', &
               (/id_ug/),  lnd%time, 'conductance between canopy and stem for '//trim(name), &
               'm/s', missing_value=-1.0)

          trdata(tr)%id_con_v_g = &
               register_tiled_diag_field(diag_name, trim(name)//'_con_v_g', &
               (/id_ug/),  lnd%time, 'conductance between canopy and ground for '//trim(name), &
               'm/s', missing_value=-1.0)


          trdata(tr)%id_Eb =                register_tiled_diag_field(diag_name, trim(name)//'_Eb', &
               (/id_ug/),  lnd%time, 'Eb collection efficiency from Brownian diffusion for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_Ein =                register_tiled_diag_field(diag_name, trim(name)//'_Ein', &
               (/id_ug/),  lnd%time, 'Ein collection efficiency from interception for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_Eim =                register_tiled_diag_field(diag_name, trim(name)//'_Eim', &
               (/id_ug/),  lnd%time, 'Eim collection efficiency from impaction for '//trim(name), &
               'm/s', missing_value=-1.0)


          trdata(tr)%id_con_g = &
               register_tiled_diag_field(diag_name, trim(name)//'_con_g', &
               (/id_ug/),  lnd%time, 'total conductance between ground and canopy air for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_conc = &
               register_tiled_diag_field(diag_name, trim(name), &
               (/id_ug/),  lnd%time, 'concentration or '//trim(name)//' in canopy air', &
               units, missing_value=-1.0)
          trdata(tr)%id_tcond = &
               register_tiled_diag_field(diag_name, trim(trdata(tr)%name)//'_tot_con', &
               (/id_ug/),  lnd%time,'total conductance of '//trim(trdata(tr)%name), &
               "m/s", missing_value=-1.0)
          trdata(tr)%id_con_mx_st = &
               register_tiled_diag_field(diag_name, trim(trdata(tr)%name)//'_con_mx_st', &
               (/id_ug/),  lnd%time, 'mesophyl + stomatal conductance for '//trim(trdata(tr)%name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_con_cu = &
               register_tiled_diag_field(diag_name, trim(trdata(tr)%name)//'_con_cu', &
               (/id_ug/),  lnd%time, 'cuticular conductance for '//trim(trdata(tr)%name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_con_stem = &
               register_tiled_diag_field(diag_name, trim(trdata(tr)%name)//'_con_stem', &
               (/id_ug/),  lnd%time, 'stem conductance for '//trim(trdata(tr)%name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_con_gr = &
               register_tiled_diag_field(diag_name, trim(trdata(tr)%name)//'_con_gr', &
               (/id_ug/),  lnd%time, 'ground conductance for '//trim(trdata(tr)%name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_econ_v = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_v', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the vegetation/stem for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_econ_g = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_g', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the ground for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_econ_stom = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_stom', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the stomata for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_econ_stem = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_stem', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the stem for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_econ_cu = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_cu', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the cuticles for '//trim(name), &
               'm/s', missing_value=-1.0)

          trdata(tr)%id_econ_cu_dry = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_cu_dry', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the dry cuticles for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_econ_cu_wet = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_cu_wet', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the wet cuticles for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_econ_cu_frz = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_cu_frz', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the frozen cuticles for '//trim(name), &
               'm/s', missing_value=-1.0)

          trdata(tr)%id_econ_g_dry = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_g_dry', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the dry cuticles for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_econ_g_wet = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_g_wet', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the wet cuticles for '//trim(name), &
               'm/s', missing_value=-1.0)
          trdata(tr)%id_econ_g_frz = &
               register_tiled_diag_field(diag_name, trim(name)//'_econ_g_frz', &
               (/id_ug/),  lnd%time, 'effective deposition velocity to the frozen cuticles for '//trim(name), &
               'm/s', missing_value=-1.0)


          trdata(tr)%id_ddep_v = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_v', &
               (/id_ug/),  lnd%time, 'deposition to the vegetation/stem for '//trim(name), &
               trim(funits), missing_value=-1.0)
          trdata(tr)%id_ddep_g = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_g', &
               (/id_ug/),  lnd%time, 'deposition to the ground for '//trim(name), &
               trim(funits), missing_value=-1.0)
          trdata(tr)%id_ddep_stom = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_stom', &
               (/id_ug/),  lnd%time, 'deposition to the stomata for '//trim(name), &
               trim(funits), missing_value=-1.0)
          trdata(tr)%id_ddep_stem = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_stem', &
               (/id_ug/),  lnd%time, 'deposition to the stem for '//trim(name), &
               trim(funits), missing_value=-1.0)
          trdata(tr)%id_ddep_cu = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_cu', &
               (/id_ug/),  lnd%time, 'deposition to the cuticles for '//trim(name), &
               trim(funits), missing_value=-1.0)
          trdata(tr)%id_ddep_cu_wet = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_cu_wet', &
               (/id_ug/),  lnd%time, 'deposition to the dry cuticles for '//trim(name), &
               trim(funits), missing_value=-1.0)
          trdata(tr)%id_ddep_cu_dry = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_cu_dry', &
               (/id_ug/),  lnd%time, 'deposition to the wet cuticles for '//trim(name), &
               trim(funits), missing_value=-1.0)
          trdata(tr)%id_ddep_cu_frz = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_cu_frz', &
               (/id_ug/),  lnd%time, 'deposition to the frozen cuticles for '//trim(name), &
               trim(funits), missing_value=-1.0)

          trdata(tr)%id_ddep_g_wet = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_g_wet', &
               (/id_ug/),  lnd%time, 'deposition to the dry cuticles for '//trim(name), &
               trim(funits), missing_value=-1.0)
          trdata(tr)%id_ddep_g_dry = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_g_dry', &
               (/id_ug/),  lnd%time, 'deposition to the wet cuticles for '//trim(name), &
               trim(funits), missing_value=-1.0)
          trdata(tr)%id_ddep_g_frz = &
               register_tiled_diag_field(diag_name, trim(name)//'_ddep_g_frz', &
               (/id_ug/),  lnd%time, 'deposition to the frozen cuticles for '//trim(name), &
               trim(funits), missing_value=-1.0)

       endif
    enddo

    
    id_con_atm = &
         register_tiled_diag_field(diag_name, 'con_atm', &
         (/id_ug/),  lnd%time,'1/Ra', &
         "m/s", missing_value=-1.0)

    id_gfrac_dry = register_tiled_diag_field(diag_name, 'gfrac_dry', &
         (/id_ug/),  lnd%time, 'ground dry fraction', &
         'unitless', missing_value=-1.0)
    id_gfrac_frz  = register_tiled_diag_field(diag_name, 'gfrac_frz', &
         (/id_ug/),  lnd%time, 'ground froz fraction', &
         'unitless', missing_value=-1.0)
    id_gfrac_wet = register_tiled_diag_field(diag_name, 'gfrac_wet', &
         (/id_ug/),  lnd%time, 'ground wet fraction', &
         'unitless', missing_value=-1.0)
    id_frac_desert = register_tiled_diag_field(diag_name, 'frac_desert', &
         (/id_ug/),  lnd%time, 'ground desert fraction', &
         'unitless', missing_value=-1.0)

    call set_default_diag_filter('soil')    
    id_h2_fm  = register_tiled_diag_field(diag_name, 'h2_fm', &
         (/id_ug/),  lnd%time, 'h2_fm', &
         'unitless', missing_value=-1.0)
    id_h2_ft  = register_tiled_diag_field(diag_name, 'h2_ft', &
         (/id_ug/),  lnd%time, 'h2_ft', &
         'unitless', missing_value=-1.0)

    id_h2_sdiff = register_tiled_diag_field(diag_name, 'h2_sdiff', &
         (/id_ug/),  lnd%time, 'h2 soil diffusivity', &
         'm2/s', missing_value=-1.0)

    id_h2_n = register_tiled_diag_field(diag_name, 'h2_n', &
         (/id_ug/),  lnd%time, 'normalization constant for the modified beta distribution of the biological sink', &
         'unitless', missing_value=-1.0)
    id_h2_st = register_tiled_diag_field(diag_name, 'h2_st', &
         (/id_ug/),  lnd%time, 'soil moisture threshold for bacterial activity', &
         'unitless', missing_value=-1.0)
    id_h2_betab = register_tiled_diag_field(diag_name, 'h2_betab', &
         (/id_ug/),  lnd%time, 'second exponent of the beta distribution', &
         'unitless', missing_value=-1.0)
    id_h2_b = register_tiled_diag_field(diag_name, 'h2_b', &
         (/id_ug/),  lnd%time, 'Cambpell exponent', &
         'unitless', missing_value=-1.0)

    do iw=1,nwet_diag
       id_fw_wet(iw) =  register_tiled_diag_field(diag_name,'f_wet_'//trim(wet_str(iw)), &
            (/id_ug/),  lnd%time, 'fraction of the time when canopy is more than '//trim(wet_str(iw))//' wet', &
            'unitless', missing_value= -1.0)

       do tr = 1, ntcana
             trdata(tr)%id_tcond_wet(iw) = register_tiled_diag_field(diag_name, trim(trdata(tr)%name)//'_tot_con_'//trim(wet_str(iw)), &
                  (/id_ug/),  lnd%time, 'total conductance of '//trim(trdata(tr)%name)//' with fwet>'//trim(wet_str(iw)), &
                  'm/s', missing_value=-1.0)
          end do
       
    end do

    id_fw_avg = &
         register_tiled_diag_field(diag_name, 'fw_avg', &
         (/id_ug/),  lnd%time,'canopy avg wet fraction', &
         "unitless", missing_value=-1.0)
    id_fs_avg = &
         register_tiled_diag_field(diag_name, 'fs_avg', &
         (/id_ug/),  lnd%time,'canopy avg snow fraction', &
         "unitless", missing_value=-1.0)
    id_fd_avg = &
         register_tiled_diag_field(diag_name, 'fd_avg', &
         (/id_ug/),  lnd%time,'canopy avg dry fraction', &
         "unitless", missing_value=-1.0)

    
    !save soil properties used for H2 soil removal
    ce = first_elmt(land_tile_map)    
    do while (loop_over_tiles(ce,tile))
       if (associated(tile%soil)) then
          soil_tag = tile%soil%tag          
          if (id_h2_n.gt.0)     call send_tile_data(id_h2_n,     h2_n(soil_tag),     tile%diag)
          if (id_h2_st.gt.0)    call send_tile_data(id_h2_st,    h2_st(soil_tag),    tile%diag)
          if (id_h2_betab.gt.0) call send_tile_data(id_h2_betab, h2_betab(soil_tag), tile%diag)
          if (id_h2_b.gt.0)     call send_tile_data(id_h2_b,     h2_b(soil_tag),     tile%diag)          
       end if
    end do
        
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
    integer :: iw      ! wet threshold index
    real    :: rho     ! density of canopy air
    real    :: dq      ! canopy air tracer tendency per time step
    real    :: con_v   ! laminar conductance to leaves
    real    :: con_st_tr ! total stomatal conductance, scaled by dry leaf area, m/s
    real    :: con_mx  ! "mesophyll conductance", m/s
    real    :: con_mx_st !mesophyll+stomatal
    real    :: con_cu  ! total cuticular conductance, including dry and wet areas, m/s
    real    :: con_gr  ! "ground conductance", m/s
    real    :: cv, cg  ! total conductances for vegetation and ground surface, m/s
    real    :: f_atm   ! flux of the tracer to the atmosphere, kg/(m2 s)
    real    :: ddep    ! dry deposition of the tracer, kg/(m2 s)
    real    :: emis(ntcana) ! tracer sources
    real    ::  ft, & ! fraction of canopy not covered by intercepted water/snow
         fw, & ! fraction of canopy covered by intercepted water
         fs    ! fraction of canopy covered by intercepted snow

    !fraction of ground that is frozen, wet, dry
    real    :: gfrac_frz,gfrac_wet,gfrac_dry

    real    :: con_cu_dry, con_cu_wet, con_cu_frz, con_stem
    real    :: con_gr_dry, con_gr_wet, con_gr_frz
    real    :: r_gs, frac_desert
    real    :: con_v_v_tr, con_v_stem_tr, con_bl_tr, con_g_tr

    real    :: con_cu_diag, con_mx_st_diag, con_stem_diag, fdiag
    real    :: con_v_v_tr_diag, con_v_stem_tr_diag
    real    :: econ_cu, econ_stem, econ_mx_st
    real    :: econ_cu_dry, econ_cu_wet, econ_cu_frz
    real    :: dvel
    real    :: tmp

    real    :: alpha_aere, gamma_aere, A_aere
    real    :: Eb, Eim, Ein
    real    :: fw_avg, fs_avg, rh


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

    call get_tile_property(tile,gfrac_frz,gfrac_wet,frac_desert)
    gfrac_dry = max(1.-gfrac_wet-gfrac_frz,0.)

    !calculate rh
    call compute_qs (tile%cana%T, pressure, rh, q=tile%cana%tr(isphum))
    RH = tile%cana%tr(isphum)/RH
    !cap RH
    RH = max(min(RH,0.995),0.)    

    !calculate the vegn fw and fs
    if (associated(tile%vegn)) then
       call vegn_tile_fw_fs(tile%vegn,fw_avg,fs_avg)
       call send_tile_data(id_fw_avg,fw_avg, tile%diag)
       call send_tile_data(id_fs_avg,fs_avg, tile%diag) 
       call send_tile_data(id_fd_avg,1.-fw_avg-fs_avg, tile%diag)
       do iw=1,nwet_diag       
          if ( fw_avg .gt. wet_diag_thr(iw) ) then
             call send_tile_data(id_fw_wet(iw), 1., tile%diag)
          else
             call send_tile_data(id_fw_wet(iw), 0., tile%diag)
          end if
       end do       
    end if

    
    ! loop for generic tracers only
    do tr = 1, ntcana

       cv             = 0.
       cg             = 0.

       if (.not.trdata(tr)%is_generic) cycle

       if (trdata(tr)%do_deposition) then

          con_cu_diag    = 0.
          con_stem_diag  = 0.
          con_mx_st_diag = 0.
          con_gr_dry     = 0.
          con_gr_wet     = 0.
          con_gr_frz     = 0.

          con_v_v_tr_diag = 0.
          con_v_stem_tr_diag = 0.
          con_g_tr       = 0.

          econ_cu        = 0.
          econ_cu_dry    = 0.
          econ_cu_wet    = 0.
          econ_cu_frz    = 0.
          econ_stem      = 0.
          econ_mx_st     = 0.

          if (.not. trdata(tr)%is_aerosol) then

             !conductance to the vegetation
             if (associated(tile%vegn)) then

                !get fraction of ground covered with snow
                do k = 1, tile%vegn%n_cohorts
                   associate(c=>tile%vegn%cohorts(k),sp=>spdata(tile%vegn%cohorts(k)%species))

                     call get_vegn_wet_frac ( c, fw=fw, fs=fs ); ft = 1-fw-fs

                     con_cu_dry  = ft * c%lai * get_conductance_tracer(trdata(tr),sp%r_cus,sp%r_cuo) / scale_r_T(c%Tv,c_dry) * exp(RH)
                     con_cu_wet  = fw * c%lai * get_conductance_tracer(trdata(tr),sp%r_cus_wet,sp%r_cuo_wet) / scale_r_T(c%Tv,c_wet)
                     con_cu_frz  = fs * c%lai * get_conductance_tracer(trdata(tr),r_snows,r_snowo)

                     !here we use the bulk leaf property for the cohort. This is different from the LM3 implementation.
                     con_cu   = con_cu_dry+con_cu_wet+con_cu_frz

                     !comment-out temperature dependence based on Clifton (2020)
                     !con_stem = c%sai * get_conductance_tracer(trdata(tr),sp%r_stems,sp%r_stemo) / scale_r_T(c%Tv)
                     con_stem = c%sai * get_conductance_tracer(trdata(tr),sp%r_stems,sp%r_stemo)

                     con_st_tr    = stomatal_cond(k) * trdata(tr)%scale_stom

                     if (trdata(tr)%r_mx .gt. 0) then
                        con_mx    = 1./trdata(tr)%r_mx
                     else
                        con_mx    = c%lai * get_conductance_tracer(trdata(tr),1.,100.)
                     end if

                     con_mx_st = conductance_series(con_mx,con_st_tr)

                     !calculate contribution of this cohort to the overall vegetation conductance
                     !con_v_v and con_stem are for H2O, we need to scale by (Di/Dw)**(2./3.)
                     con_v_v_tr     = con_v_v(k)*1./trdata(tr)%diff_ratio**(2./3.)
                     con_v_stem_tr  = con_v_stem(k)*1./trdata(tr)%diff_ratio**(2./3.)

                     econ_mx_st       = econ_mx_st + c%layerfrac*con_mx_st/(con_mx_st+con_cu+epsln)*conductance_series(con_v_v_tr,con_mx_st+con_cu)

                     tmp              = c%layerfrac * con_cu/(con_mx_st+con_cu+epsln)*conductance_series(con_v_v_tr,con_mx_st+con_cu)

                     econ_cu          = econ_cu     + tmp
                     econ_cu_wet      = econ_cu_wet + tmp * con_cu_wet/(con_cu+epsln)
                     econ_cu_frz      = econ_cu_frz + tmp * con_cu_frz/(con_cu+epsln)
                     econ_cu_dry      = econ_cu_dry + tmp * con_cu_dry/(con_cu+epsln)

                     econ_stem        = econ_stem   + c%layerfrac*conductance_series(con_v_stem_tr,con_stem)

                     !for diagnostics
                     con_mx_st_diag   = con_mx_st_diag + c%layerfrac*con_mx_st
                     con_cu_diag      = con_cu_diag    + c%layerfrac*con_cu
                     con_stem_diag    = con_stem_diag  + c%layerfrac*con_stem

                     con_v_v_tr_diag     = con_v_v_tr_diag + c%layerfrac*con_v_v_tr
                     con_v_stem_tr_diag  = con_v_stem_tr_diag + c%layerfrac*con_v_stem_tr

                   end associate
                end do

                cv = econ_mx_st + econ_cu + econ_stem
             end if

             !note that for the ground we are calculating the tile average
             if (trdata(tr)%is_h2) then
                !for now set constant conductance
                con_gr_dry = con_h2(trdata(tr),tile,pressure)
             else
                con_gr_dry = gfrac_dry     * get_conductance_tracer(trdata(tr),r_gs_dry,r_go_dry) * 1./scale_r_T(land_tile_grnd_T(tile),c_dry) * 1./scale_biomass(frac_desert)
             end if

             con_gr_frz =    gfrac_frz     * get_conductance_tracer(trdata(tr),r_snows,r_snowo) * 1./scale_r_T(land_tile_grnd_T(tile),c_snow)

             con_gr_wet =    gfrac_wet     * get_conductance_tracer(trdata(tr),r_gs_wet,r_go_wet)  * 1./scale_r_T(land_tile_grnd_T(tile),c_wet)

             if (associated(tile%lake)) then
                if (tile%lake%ws(1).le.ws_min) then
                   con_gr_wet = get_conductance_tracer(trdata(tr),r_gs_lake,r_go_lake)
                end if
             end if


             con_bl_tr  = 1./(r_bl_h2o+epsln)   * 1./trdata(tr)%diff_ratio**(2./3.)

             con_g_tr = conductance_series(con_g,con_bl_tr)

             con_gr   = con_gr_dry+con_gr_wet+con_gr_frz
             cg       = conductance_series(con_gr,con_g_tr)

             call send_tile_data(trdata(tr)%id_con_mx_st,con_mx_st_diag, tile%diag)
             call send_tile_data(trdata(tr)%id_con_cu,con_cu_diag, tile%diag)
             call send_tile_data(trdata(tr)%id_con_stem,con_stem_diag, tile%diag)
             call send_tile_data(trdata(tr)%id_con_gr,con_gr, tile%diag)
             call send_tile_data(trdata(tr)%id_con_v_v, con_v_v_tr_diag, tile%diag)
             call send_tile_data(trdata(tr)%id_con_v_stem, con_v_stem_tr_diag,  tile%diag)
             call send_tile_data(trdata(tr)%id_con_v_g, con_g_tr, tile%diag)

          else
             !aerosol
             alpha_aere   = 0.
             gamma_aere   = 0.
             A_aere       = 0.

             if (associated(tile%vegn)) then
                tmp =  0.
                do k = 1, tile%vegn%n_cohorts
                   associate(c=>tile%vegn%cohorts(k),sp=>spdata(tile%vegn%cohorts(k)%species))
                     A_aere     = A_aere     + c%layerfrac*sp%A_aer
                     gamma_aere = gamma_aere + c%layerfrac*sp%gamma_aer
                     alpha_aere = alpha_aere + c%layerfrac*sp%alpha_aer
                     tmp        = tmp        + c%layerfrac
                   end associate
                end do

                A_aere     = A_aere/(tmp+epsln)
                gamma_aere = gamma_aere/(tmp+epsln)
                alpha_aere = alpha_aere/(tmp+epsln)

                gamma_aere = gamma_aere*(1-frac_desert)+gamma_aer_desert*frac_desert
                alpha_aere = alpha_aere**(1.-frac_desert)+alpha_aer_desert*frac_desert
             end if

             if (associated(tile%glac)) then
                A_aere       = -999.
                gamma_aere   = gamma_aer_frz
                alpha_aere   = alpha_aer_frz
             end if

             if (associated(tile%lake)) then
                A_aere       = -999.
                gamma_aere   = gamma_aer_lake
                alpha_aere   = alpha_aer_lake
             end if

             cv = 0.
             call cg_aer(trdata(tr),tile%cana%T,ustar,pressure,alpha_aere,gamma_aere,A_aere,gfrac_wet,frac_desert,cg,Eb,Eim,Ein)

             if (cg_aer_frz.gt.0.) &
                  cg = cg*(1.-gfrac_frz)+cg_aer_frz*gfrac_frz

             call send_tile_data(trdata(tr)%id_Eb, Eb,  tile%diag)
             call send_tile_data(trdata(tr)%id_Eim,Eim, tile%diag)
             call send_tile_data(trdata(tr)%id_Ein,Ein, tile%diag)


          endif
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
       ddep  = rho*(cv+cg)*tile%cana%tr(tr)*trdata(tr)%conv_flux
       f_atm = tr_flux(tr)+dfdtr(tr)*dq*trdata(tr)%conv_flux
       ! ---- diagnostic section
       dvel = con_atm*(cv+cg)/(con_atm+cv+cg)
       call send_tile_data(trdata(tr)%id_con_v,      cv,         tile%diag)
       call send_tile_data(trdata(tr)%id_con_g,      cg,         tile%diag)
       call send_tile_data(trdata(tr)%id_emis,       emis(tr),   tile%diag)
       call send_tile_data(trdata(tr)%id_ddep,       ddep,       tile%diag)
       call send_tile_data(trdata(tr)%id_flux_atm,   f_atm,      tile%diag)
       call send_tile_data(trdata(tr)%id_tcond,dvel, tile%diag)

       if (associated(tile%vegn)) then
          do iw=1,nwet_diag
             if (trdata(tr)%id_tcond_wet(iw).gt.0) then
                if (fw_avg.gt.wet_diag_thr(iw)) then
                   call send_tile_data(trdata(tr)%id_tcond_wet(iw),dvel, tile%diag)
                else
                   call send_tile_data(trdata(tr)%id_tcond_wet(iw),0., tile%diag)
                end if
             end if
          end do
       end if
       
       !save deposition to the vegetation and ground
       fdiag = min(max(cv/(cv+cg+epsln),0.),1.)
       call send_tile_data(trdata(tr)%id_econ_g, (1.-fdiag)*dvel,tile%diag)
       call send_tile_data(trdata(tr)%id_econ_v, fdiag*dvel,     tile%diag)
       call send_tile_data(trdata(tr)%id_ddep_g, (1.-fdiag)*ddep,tile%diag)
       call send_tile_data(trdata(tr)%id_ddep_v, fdiag*ddep,     tile%diag)

       call send_tile_data(trdata(tr)%id_econ_g_wet,   con_gr_wet/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*dvel,     tile%diag)
       call send_tile_data(trdata(tr)%id_econ_g_dry,   con_gr_dry/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*dvel,     tile%diag)
       call send_tile_data(trdata(tr)%id_econ_g_frz,   con_gr_frz/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*dvel,     tile%diag)

       call send_tile_data(trdata(tr)%id_ddep_g_wet,   con_gr_wet/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*ddep,     tile%diag)
       call send_tile_data(trdata(tr)%id_ddep_g_dry,   con_gr_dry/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*ddep,     tile%diag)
       call send_tile_data(trdata(tr)%id_ddep_g_frz,   con_gr_frz/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*ddep,     tile%diag)


       call send_tile_data(trdata(tr)%id_econ_cu,       fdiag*dvel*econ_cu/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
       call send_tile_data(trdata(tr)%id_econ_cu_wet,   fdiag*dvel*econ_cu_wet/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
       call send_tile_data(trdata(tr)%id_econ_cu_dry,   fdiag*dvel*econ_cu_dry/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
       call send_tile_data(trdata(tr)%id_econ_cu_frz,   fdiag*dvel*econ_cu_frz/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)

       call send_tile_data(trdata(tr)%id_econ_stem, fdiag*dvel*econ_stem/(econ_cu+econ_stem+econ_mx_st+epsln),   tile%diag)
       call send_tile_data(trdata(tr)%id_econ_stom, fdiag*dvel*econ_mx_st/(econ_cu+econ_stem+econ_mx_st+epsln),  tile%diag)

       call send_tile_data(trdata(tr)%id_ddep_cu,       fdiag*ddep*econ_cu/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
       call send_tile_data(trdata(tr)%id_ddep_cu_wet,   fdiag*ddep*econ_cu_wet/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
       call send_tile_data(trdata(tr)%id_ddep_cu_dry,   fdiag*ddep*econ_cu_dry/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
       call send_tile_data(trdata(tr)%id_ddep_cu_frz,   fdiag*ddep*econ_cu_frz/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)

       call send_tile_data(trdata(tr)%id_ddep_stem, fdiag*ddep*econ_stem/(econ_cu+econ_stem+econ_mx_st+epsln),   tile%diag)
       call send_tile_data(trdata(tr)%id_ddep_stom, fdiag*ddep*econ_mx_st/(econ_cu+econ_stem+econ_mx_st+epsln),  tile%diag)


    enddo
        
    ! send concentrations for all tracers, generic or not
    do tr = 1, ntcana
       call send_tile_data(trdata(tr)%id_conc,       tile%cana%tr(tr), tile%diag)
       call send_tile_data(trdata(tr)%id_dfdtr,      dfdtr(tr),  tile%diag)
    enddo

    call send_tile_data(id_con_atm,      con_atm,      tile%diag)
    call send_tile_data(id_gfrac_dry,    gfrac_dry,    tile%diag)
    call send_tile_data(id_gfrac_wet,    gfrac_wet,    tile%diag)
    call send_tile_data(id_gfrac_frz,    gfrac_frz,   tile%diag)
    call send_tile_data(id_frac_desert,  frac_desert,  tile%diag)

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

    if (associated(tile%lake)) then
       if (tile%lake%ws(1).gt.ws_min) then
          !lake is assumed to be frozen
          gfrac_frz = 1.
          gfrac_wet = 0.
       else
          gfrac_wet = 1.
          gfrac_frz = 0.
       end if
    end if

    if (associated(tile%glac)) then
       gfrac_frz = 1.
       gfrac_wet = 0.
    endif

    if (associated(tile%snow)) then
       call snow_get_depth_area ( tile%snow, snow_depth, snow_area )
       gfrac_frz = snow_area
    end if

    if (associated(tile%vegn)) then
       theta = soil_theta(tile%soil) !top layer
       if ( theta(1) .gt. theta_wetland_thr ) then
          gfrac_wet = max(1. - gfrac_frz,0.)
       else
          gfrac_wet = 0.
       end if

       biomass = 0.
       associate(cc=>tile%vegn%cohorts)
         do k = 1, tile%vegn%n_cohorts
            biomass = biomass +            &
                 cc(k)%layerfrac * (       &
                 cc(k)%bl  + cc(k)%blv   + &
                 cc(k)%br  + cc(k)%bwood + &
                 cc(k)%bsw + cc(k)%bseed + &
                 cc(k)%nsc )
         end do
       end associate

       if ( biomass .lt. desert_biomass ) then
          frac_desert = min(max(1.-biomass/desert_biomass,0.),1.)
       end if


    end if

  end subroutine get_tile_property

  elemental real function calc_rt(tk) result(r_t)
    !scaling factor for resistance at cold temperature
    real, intent(in) :: tk

    r_t   = max(min(2.,exp(0.2*(-1-(tk-273.15)))),1.)
  end function calc_rt

  elemental real function get_conductance_tracer(tr_data,r_s,r_o) result(con)
    real, intent(in)                   :: r_s, r_o
    type(tracer_data_type), intent(in) :: tr_data

    con = tr_data%alpha/r_s + tr_data%reactivity/r_o
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

  elemental real function conductance_series(con1,con2) result(con)
    real, intent(in) :: con1, con2

    con = (con1*con2)/(con1+con2+epsln)

  end function conductance_series

  elemental real function scale_snow_T(T) result(s)
    !erisman (1994) showed that resistance increases as temperature decreases from 70 to 500
    real, intent(in) :: T
    s = max(min(max_scale_snow_T,275.15-T),1.)
  end function scale_snow_T


  elemental real function scale_r_T(T,c) result(s)
    !zhang (2003) equation 10a
    !updated based on Clifton (2020)
    real, intent(in) :: T,c
    !    s = min(max(exp(0.2*(-1-T+273.15)),max_scale_cold_T),1.)
    s = max(exp(-c*(T-273.15-5)),1.)
  end function scale_r_T

  elemental real function scale_biomass(frac_desert) result(s)
    real, intent(in) :: frac_desert

    s = max(frac_desert*max_scale_desert+(1-frac_desert),1.)

  end function scale_biomass

  subroutine  cg_aer(tr_data,T,ustar,pressure,alpha,gamma,A,frac_wet,frac_desert,con, Eb, Eim, Ein)

    !from Zhang 2001 as presented by Seinfeld (19.27)
    !rb = 1/ (3*ustar * (Sc^-gamma + (St/(alpha+St))^2 + 1/2 * (Dp/A)^2) * R1 )
    type(tracer_data_type), intent(in) :: tr_data

    real, intent(in) :: T !temperature
    real, intent(in) :: alpha, gamma, A !depend on land type
    real, intent(in) :: pressure
    real, intent(in) :: ustar
    real, intent(in) :: frac_wet, frac_desert

    real, intent(out):: con, Eb, Eim, Ein
    !Eb:  collection efficiency from Brownian motion
    !Eim: collection efficiency from impaction
    !Ein: collection efficiency from interception

    real             :: Sc, diff_aer
    real             :: R1, St
    !R1: correction factor representing the fractionof particles that stick to the surface
    real             :: free_path, C_c
    real             :: kvis,dvis
    real             :: rwet, ratio_r, rho_wet, vts

    real, parameter  :: e0   = 3.
    real, parameter  :: beta = 2.

    real             :: rho_p, rp

    rho_p = tr_data%rho
    rp    = tr_data%radius

    kvis = kin_visc_air(T,pressure)
    dvis = dyn_visc_air(T)

    !calculate settling velocity
    rwet          = rp !TODO: hygroscopic growth
    ratio_r       = (rp/rwet)**3  ! Ratio dry over wet radius cubic power
    rho_wet       = ratio_r*rho_p+(1.-ratio_r)*DENS_H2O ! Density of wet aerosol [kg/m3]

    !from sedimentation_velocity [land_dust]
    free_path = 6.6e-8*T/293.15*(PSTD_MKS/pressure)
    C_c        = 1.0 + free_path/rp * (1.257+0.4*exp(-1.1*rp/free_path)) ! slip correction
    vts       = 2./9.*C_c*GRAV*rho_p*rp**2/dvis  ! Settling velocity [m/s]

    if  ( A .gt. epsln ) then
       St       = vts*ustar/(grav*A)*(1.-frac_desert) + vts*ustar**2/(grav*kvis)*frac_desert
       Ein      = (0.5*(2*rp/A)**2)*(1.-frac_desert) !Ein is 0. for desert
    else
       St       = vts*ustar**2/(grav*kvis)
       Ein = 0.
    end if

    R1      = frac_wet + (1-frac_wet)*exp(-St**0.5)

    !brownian diffusion
    diff_aer= kb*T*C_c / ( 6.*pi*dvis*rp )
    Sc      = kvis/diff_aer

    Eb      = 1./Sc**gamma
    Eim     = (St/(alpha+St))**beta

    con     = e0 * ustar * R1 * (Eb + Eim + Ein)

  end subroutine cg_aer

  real function con_h2(tr_data,tile,p) result(con)

    type(tracer_data_type), intent(in) :: tr_data
    type(land_tile_type),   intent(inout) :: tile
    real,                   intent(in) :: p

    real    :: delta !inactive layer
    integer :: isoil, soil_tag
    real    :: frac_water_pores_avg, frac_ice_pores_avg, frac_air_pores_avg, porosity
    real,dimension(num_l) :: frac_water_pores, frac_ice_pores
    real    :: dz, T_avg, T_avgC
    real    :: diff_h2_air, diff_h2
    real    :: f_T, f_M
    real    :: h2_gamma, h2_km, h2_depth
    real    :: gdelta, snow_depth, snow_area
    real, parameter :: snow_porosity = 0.64

    con = 0.

    if (tr_data%con_gr.gt.0.) then
       con = tr_data%con_gr

    else
       delta = 0. !no inactive layer for now

       if (associated(tile%soil)) then

          h2_gamma = tr_data%gamma
          h2_km    = tr_data%km
          h2_depth = tr_data%depth

          soil_tag = tile%soil%tag

          frac_water_pores  = soil_theta(tile%soil)
          frac_ice_pores    = soil_ice_porosity(tile%soil)

          porosity          = tile%soil%pars%vwc_sat

          isoil = 1
          T_avg = 0.
          frac_water_pores_avg = 0.
          frac_ice_pores_avg   = 0.

          do while (zhalf(isoil).lt.h2_depth)
             dz                   = min(zhalf(isoil+1),h2_depth)-zhalf(isoil)
             T_avg                = T_avg + dz*tile%soil%T(isoil)
             frac_water_pores_avg = frac_water_pores_avg + dz*frac_water_pores(isoil)
             frac_ice_pores_avg   = frac_ice_pores_avg + dz*frac_ice_pores(isoil)
             isoil                = isoil+1
          end do

          T_avg                 = T_avg/(h2_depth-delta)
          T_avgC                = T_avg - 273.15

          frac_water_pores_avg  = max(frac_water_pores_avg/(h2_depth-delta),0.)
          frac_ice_pores_avg    = max(frac_ice_pores_avg/(h2_depth-delta),0.)
          frac_air_pores_avg    = max(1.-frac_water_pores_avg-frac_ice_pores_avg,0.)

          diff_h2_air    = 0.668e-4*(101325./p)*(T_avg/273)**1.75
          diff_h2        = diff_h2_air * porosity**2.*frac_air_pores_avg**(2+3./h2_b(soil_tag))

          f_T = 1/(1+exp(-(T_avgC-3.8)/6.7)) + 1./(1.+exp((T_avgC - 62.2)/7.7)) - 1.

          if (frac_water_pores_avg.gt.1.) then
             f_M = 0.
          elseif (frac_water_pores_avg.lt.h2_st(soil_tag)) then
             f_M = 0.
          else
             f_M = 1/h2_N(soil_tag)*(frac_water_pores_avg-h2_st(soil_tag))**(h2_a-1)*(1.-frac_water_pores_avg)**(h2_betab(soil_tag)-1)
          end if

          con = sqrt(h2_gamma*f_T*f_M*h2_km*diff_h2)

          if (associated(tile%snow)) then
             call snow_get_depth_area ( tile%snow, snow_depth, snow_area )
             if (snow_depth.gt.epsln) then
                snow_depth = snow_depth/(snow_area+epsln)
                gdelta = diff_h2_air*snow_porosity/(snow_depth+epsln)

                con = (1.-snow_area)*con + snow_area*(con*gdelta)/(con+gdelta)
             end if
          end if

          call send_tile_data(id_h2_fm,f_M,                       tile%diag)
          call send_tile_data(id_h2_ft,f_T,                       tile%diag)
          call send_tile_data(id_h2_sdiff,diff_h2,                tile%diag)
                    
       end if

    end if

  end function con_h2

end module land_tracer_driver_mod
