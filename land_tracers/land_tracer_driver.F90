module land_tracer_driver_mod

#include "../shared/debug.inc"

use constants_mod,      only: rdgas,wtmair,grav,pi,pstd_mks,avogno,DENS_H2O,epsln, WTMH2O
use field_manager_mod , only: MODEL_ATMOS, MODEL_LAND, parse
use fms_mod,            only: lowercase, stdout, stdlog, mpp_pe, mpp_root_pe, check_nml_error, error_mesg, WARNING, FATAL, NOTE
use fms_mod,            only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_MODULE
use ieee_arithmetic
use mpp_mod,            only: input_nml_file
use table_printer_mod

use cana_tile_mod,      only: canopy_air_mass_for_tracers
use snow_tile_mod,      only: N_SNOW_TRACERS, SNOW_TR_BC, SNOW_TR_MD, SNOW_TR_OM
use land_constants_mod, only: d608,kin_visc_air,dyn_visc_air,N_LITTER_POOLS,LITT_LEAF
use land_data_mod,      only: lnd, log_version
use land_debug_mod,     only: is_watch_point, check_var_range, check_temp_range
use land_dust_mod,      only: land_dust_init, land_dust_end, update_land_dust
use land_tracers_mod,   only: ntcana, isphum, ico2
use land_tile_mod,      only: land_tile_type, land_tile_grnd_T, loop_over_tiles, first_elmt, land_tile_enum_type, land_tile_map
use land_tile_diag_mod, only: diag_buff_type, set_default_diag_filter, register_tiled_diag_field, send_tile_data
use sat_vapor_pres_mod, only: compute_qs
use soil_tile_mod,      only: num_l, soil_theta, soil_ice_porosity, zhalf, n_dim_soil_types
use soil_mod,           only: soil_get_sfc_temp
use time_manager_mod,   only: time_type, time_type_to_real
use tracer_manager_mod, only: NO_TRACER, get_tracer_index, get_tracer_names, query_method
use vegn_data_mod,      only: spdata
use vegn_cohort_mod,    only: get_vegn_wet_frac
use vegn_tile_mod,      only: vegn_tile_fw_fs

! import interfaces from non-generic tracer modules, e.g.:

implicit none
private

! ==== public interfaces =====================================================
public :: land_tracer_driver_init
public :: land_tracer_driver_end
public :: update_cana_tracers
! ==== end of public interfaces ==============================================


real, parameter  :: litter_leaf_density_C = 0.03/2.*1000. !https://doi.org/10.1093/sjaf/33.1.29
real, parameter  :: litter_leaf_porosity = 0.955 !leaf: 0.955; Twig: 0.8; decomposing matter 0.45-0.55
!NOTE: This is not correct. Wood in lm4p2 means "fallen branches" not duff.
!real, parameter  :: litter_wood_density_C = 0.1/2. *1000. !https://doi.org/10.1093/sjaf/33.1.29
!real, parameter  :: litter_wood_porosity = 0.7   !10 .3389/fmech.2019.00053/full

!min/max S snow resistances
real :: r_snows_min     = 100
real :: r_snows_max     = 500
!maximum increase associated with low biomass
real :: max_scale_desert=2.5
real :: cg_aer_frz = 0.03e-2   !from fisher (2011) vd(aer) = 0.03 cm/s, conductance of aerosol over snow.

real :: desert_biomass    = 0.25 !kg/m2
real :: ws_min            = 1. !1mm  (threshold to decide whether lake is frozen)
real :: theta_wetland_thr = 0.9 !above this call it wetland

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

real :: b_lai_aer              = 0 !there is some confusion in the litterature on whether cg_aer should (b_lai_aer=1) or shoundn't (b_lai_aer=0) be multiplied by lai.

!Note about default resistances (specified in field table)
!From Zhang et al. "A size-segregated particle dry deposition scheme for an atmospheric aerosol module" Atmospheric Environment 35 (2001) 549-560. We do not account for seasonal variations in A (use midsummer)
!From Zhang et al. "A revised parameterization for gaseous dry deposition in air-quality models" Atmos. Chem. Phys., 3, 2067–2082, 2003

!LUC refers to the entry in Table 3 of Zhang et al. (2001)

!prioria (LUC=2) - Evergreen broadleaf trees
!r_cus=2500, r_cuo=6000, r_cuo_wet = 400, A_aer=5.e-3, gamma_aer=0.58, alpha_aer=0.6
!
!picea (LUC=1) - Evergreen needleleaf trees
!r_cus=2000, r_cuo=4000, r_cuo_wet=200, A_aer=2.e-3, gamma_aer=0.56, alpha_aer=1.
!
!larix (LUC=3) - Deciduous needleleaf trees
!r_cus=2000, r_cuo=4000, r_cuo_wet=200, A_aer=2.e-3, gamma_aer=0.56, alpha_aer=1.1
!
!acer (LUC=4) - Deciduous broadleaf trees
!r_cus=2500, r_cuo=6000, r_cuo_wet=400, A_aer=5.e-3, gamma_aer=0.56, alpha_aer=0.8
!
!c4grass (LUC=6) - Grass
!r_cus=1000, r_cuo=4000, r_cuo_wet=200, A_aer=2.e-3, gamma_aer=0.54, alpha_aer=1.2
!
!c3grass (LUC=6) - Grass
!r_cus=1000, r_cuo=4000, r_cuo_wet=200, A_aer=2.e-3, gamma_aer=0.54, alpha_aer=1.2

real :: c_snow=0.025, c_dry=0.1, c_wet=0.9 !strength of R increase with decreasing T under <5C (Clifton 2020)

real :: e_lai_dry=1.0, e_lai_frz=1.0, e_lai_wet=1.0 ! exponent for LAI dependence of cuticle
! conductance (lai**e_lai). Zhang et al. (2003, doi:10.5194/acp-3-2067-2003) proposes
! 0.5 for dry condition, 0.25 for wet condition, and 0 for snow covered leaves
real :: e_ustar=0.0 ! exponent for ustar dependence of cuticle conductance (u_star**e_ustar).
! Zhang proposed 1.

!for h2
real    :: h2_km           =  0.03    !prefactor (s-1)
real    :: h2_depth        =  0.1     !depth over which soil water and soil carbon are averaged for h2 calculation, m
real    :: h2_psi_ws       = -100e2   !minimum psi for HA-HOB activation (m) - 1MPa is 100m
real    :: h2_psi_opt      = -0.5e2   !optimal psi for HA-HOB (m)
real    :: h2_beta1        = 1.       !exponent (see Bertagni (GBC, 2021)

real    :: h2_precip_min   = -1       !if precip_ann is less than h2_precip_min, reduce h2 uptake


character(32) :: h2_soilC_mod      = "NONE"         !name of the soilC parameterization
real          :: h2_soilC_param(2) = (/-1,-1/)!Definition depends on soilC modulation
                                              !Similar to Paulot (2021) but h2_soilC(1) is in in kgC/m2 (h2_soilC(2) is not used)
                                              !For Reji (2024): h2_soilC(1) and h2_soilC(2) are unitless (a*soilC+b)
real    :: h2_litterC_mod  = 1.       !scale litter depth (turned off if negative)

logical :: pmod_lai_frz, pmod_lai_wet, pmod_lai_dry


namelist /land_tracer_driver_nml/ &
            r_snows_min,  r_snows_max, max_scale_desert, cg_aer_frz, &
            desert_biomass,ws_min,theta_wetland_thr, &
            r_gs_lake,r_gs_wet,r_gs_dry,r_snows,     &
            r_go_lake,r_go_wet,r_go_dry,r_snowo,     &
            A_aer_lake,A_aer_swamp,                  &
            gamma_aer_lake,gamma_aer_swamp,gamma_aer_desert,gamma_aer_frz,    &
            alpha_aer_lake,alpha_aer_swamp,alpha_aer_desert,alpha_aer_frz,    &
            h2_psi_ws, h2_psi_opt, h2_beta1, h2_km, h2_depth, h2_soilC_mod, h2_soilC_param, h2_litterC_mod, h2_precip_min, &
            c_snow, c_dry, c_wet, e_lai_dry,e_lai_wet,e_lai_frz, e_ustar, &
            r_snows_max, r_snows_max, b_lai_aer

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

integer        :: h2_soilC_mod_id    !identifier for h2_soilC_mod

!parameterization
integer, parameter :: AEROSOL_DEFAULT = -1
integer, parameter :: GAS_PARAM       = 0
integer, parameter :: GAS_DEFAULT     = 1
integer, parameter :: GAS_BERTAGNI    = 2
integer, parameter :: GAS_CODEP       = 3

integer, parameter :: H2_SOILC_NO_MOD   = -1
integer, parameter :: H2_SOILC_PAULOT21 = 1
integer, parameter :: H2_SOILC_REJI25   = 2

! ---- data types -----------------------------------------------------------
type :: tracer_data_type
   character(32)  :: name = ''          ! tracer name
   integer        :: tr_atm  = NO_TRACER      ! index of this tracer in atmos tracer array
   logical        :: is_generic    = .TRUE.   ! flag of generic tracer; initialization of non-generic tracers should turn it to FALSE
   logical        :: do_deposition = .FALSE.  ! if true, generic dry deposition is used
   ! dry deposition parameters. The default values are set as O3 parameters from (Wesely, 1989)
   real           :: reactivity    = 0.0      ! normalized reactivity factor
   real           :: alpha         = 0.0      ! scaling factor relative to SO2
   real           :: r_mx          = 1e-5     ! negligible resistance
   real           :: mw            = -9999.9  ! kg/mol
   real           :: diff_ratio    = 0      ! ratio of water vapor molecular diffusivity in the air to that of the tracer, unitless
   real           :: scale_stom    = 1.       ! additional
   real           :: diff_ratio23  = 1
   integer        :: nb_n_ox  = 0
   integer        :: nb_n_red = 0
   !for aerosol
   real           :: radius = 0.25e-6, rho = 1500.
   integer        :: parameterization = 0
   character(32)  :: param_name = "N/A"

   real           :: a_RH = 3.

   integer        :: km_mod=-1

   logical        :: is_vmr

   integer        :: & ! diag field IDs
                     id_emis,      id_ddep,  &
                     id_flux_atm,  id_dfdtr, &
                     id_con_v,     id_con_g, &
                     id_con_mx_st, id_con_cu, id_con_stem, id_con_gr, &
                     id_conc,      id_tcond, id_tcond_wet(nwet_diag), id_tcond_new, &
                     id_econ_v,    id_econ_g, &
                     id_ddep_v,    id_ddep_g, &
                     id_codep,                &
                     id_econ_g_dry,id_econ_g_wet, id_econ_g_frz, &
                     id_ddep_g_dry,id_ddep_g_wet, id_ddep_g_frz, &
                     id_econ_stem, id_econ_stom, id_econ_cu, &
                     id_ddep_stem, id_ddep_stom, id_ddep_cu, &
                     id_econ_cu_dry,id_econ_cu_wet, id_econ_cu_frz, &
                     id_ddep_cu_dry,id_ddep_cu_wet, id_ddep_cu_frz, &
                     id_con_v_v, id_con_v_stem, id_con_v_g, &
                     id_Eb, id_Ein, id_Eim

   !for some compounds, we use the same exact parameterization. This a
   !integer        :: map_to_index
   !character*32   :: map_to
   !real, allocatable:: con_cu_dry(:), con_cu_wet(:), con_cu_frz(:),  con_stem(:), con_mx(:)
   !real             :: con_gr_lake, con_gr_frz, con_gr_dry, con_gr_wet

end type tracer_data_type

integer :: id_fw_avg, id_fs_avg, id_fd_avg
integer :: id_fw_wet(nwet_diag)
integer :: id_con_atm
integer :: id_gfrac_dry, id_gfrac_wet, id_gfrac_frz, id_frac_desert
integer :: id_h2_fm, id_h2_ft, id_h2_sdiff, id_h2_ilayer, id_h2_km
integer :: id_h2_frac_water_pores_avg, id_h2_frac_ice_pores_avg
integer :: id_h2_R_bact, id_h2_R_inactive, id_h2_R_snow, id_h2_R_litter
integer :: id_h2_sws, id_h2_sopt, id_h2_sup, id_h2_moist_r1, id_h2_moist_r2
integer :: id_h2_depth_litter, id_h2_soilC
integer :: id_con_h2_no_snow, id_con_h2_no_litter

integer :: id_ddep_noy, id_ddep_nhx
integer :: id_ddep_bc, id_ddep_bc_neg, id_ddep_bc_neg_freq
integer :: id_ddep_oa, id_ddep_oa_neg, id_ddep_oa_neg_freq
integer :: id_ddep_md, id_ddep_md_neg, id_ddep_md_neg_freq
integer :: id_acid_ratio

integer :: nomphilic, nbcphilic, nomphobic, nbcphobic,nsoa, nh2
integer :: nso2, nhno3, nnh3, nh2so4

! ---- private module variables ----------------------------------------------
logical :: module_is_initialized = .FALSE.
logical :: do_generic_ddep_emis = .FALSE.  ! if true, do dust calculations
real, save :: dt ! fast time step, s
type(tracer_data_type), allocatable :: trdata(:)

integer :: land_tracer_clock, land_tracer_ddep_clock, land_tracer_ddep_aerosol_clock, land_tracer_ddep_gas_clock, land_tracer_ddep_gas_h2_clock
integer :: land_tracer_ddep_gas_vegn_clock, land_tracer_ddep_gas_grnd_clock

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine land_tracer_driver_init(id_ug,id_zfull)
   integer,intent(in) :: id_ug, id_zfull !<Unstructured axis id.
   integer                  :: tr  ! tracer index
   real                     :: value ! temporary storage for parsing input
   character(32)            :: name, units, funits ! name and units of the tracer and flux
   character(128)           :: longname ! long name of the tracer
   character(64)            :: method
   character(1024)          :: parameters
   type(table_printer_type) :: table

   integer :: unit         ! unit for namelist i/o
   integer :: io           ! i/o status for the namelist
   integer :: ierr         ! error code, returned by i/o routines

   integer :: iw

   ! write the version and tag name to the logfile
   call log_version(version, module_name, &
   __FILE__)

   !read namelist
   read (input_nml_file, nml=land_tracer_driver_nml, iostat=io)
   ierr = check_nml_error(io, 'land_tracer_driver_nml')
   if (mpp_pe() == mpp_root_pe()) then
      unit=stdlog()
      write(unit, nml=land_tracer_driver_nml)
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

   if (ntcana==2) then
      do_generic_ddep_emis = .FALSE.
      return
   else
      do_generic_ddep_emis = .TRUE.
   end if

   h2_soilc_mod_id = -1
   if (trim(lowercase(h2_soilc_mod)).eq."paulot21") then
      h2_soilC_mod_id = H2_SOILC_PAULOT21
      if (h2_soilC_param(1) < 0) then
         call error_mesg("land_tracer_driver","using default parameter for H2_PAULOT", NOTE)
         h2_soilC_param(1) = 7 !kgC/m3
      end if
   elseif (trim(lowercase(h2_soilc_mod)).eq."reji25") then
      h2_soilC_mod_id = H2_SOILC_REJI25
      !soilC in %

      if (h2_soilC_param(1) < 0) then
         call error_mesg("land_tracer_driver","using default parameter for H2_REJI [1]", NOTE)
         h2_soilC_param(1) = 0.0024 !slope
      end if
      if (h2_soilC_param(2) < 0) then
            call error_mesg("land_tracer_driver","using default parameter for H2_REJI [2]", NOTE)
            h2_soilC_param(2) = 0.00387 !minimum (soilC = 0)
      end if
   elseif (trim(lowercase(h2_soilc_mod)).eq."none") then
      h2_soilC_mod_id = H2_SOILC_NO_MOD
   else
      call error_mesg("land_tracer_driver","SoilC H2 modulation not recognized", FATAL)
   end if

   !hard-coded deposition fields for cmip -  need to be defined early for sanity checks later
   id_ddep_bc  = register_tiled_diag_field(diag_name, 'bc_ddep', (/id_ug/),  lnd%time, &
            'dry deposition of black carbon', 'kg/m2/s', missing_value=-1.0e20)
   id_ddep_bc_neg  = register_tiled_diag_field(diag_name, 'bc_ddep_neg', (/id_ug/),  lnd%time, &
            'negative dry deposition of black carbon', 'kg/m2/s', missing_value=-1.0e20)
   id_ddep_bc_neg_freq  = register_tiled_diag_field(diag_name, 'bc_ddep_neg_freq', (/id_ug/),  lnd%time, &
            'frequency of negative dry deposition of black carbon', 'kg/m2/s', missing_value=-1.0e20)

   id_ddep_oa  = register_tiled_diag_field(diag_name, 'oa_ddep', (/id_ug/),  lnd%time, &
            'dry deposition of organic aerosols', 'kg/m2/s', missing_value=-1.0e20)
   id_ddep_oa_neg  = register_tiled_diag_field(diag_name, 'oa_ddep_neg', (/id_ug/),  lnd%time, &
            'negative dry deposition of organic aerosols', 'kg/m2/s', missing_value=-1.0e20)
   id_ddep_oa_neg_freq  = register_tiled_diag_field(diag_name, 'oa_ddep_neg_freq', (/id_ug/),  lnd%time, &
            'frequency of negative dry deposition of organic aerosols', 'kg/m2/s', missing_value=-1.0e20 )

   id_ddep_md  = register_tiled_diag_field(diag_name, 'md_ddep', (/id_ug/),  lnd%time, &
            'dry deposition of mineral dust', 'kg/m2/s', missing_value=-1.0e20)
   id_ddep_md_neg  = register_tiled_diag_field(diag_name, 'md_ddep_neg', (/id_ug/),  lnd%time, &
            'negative dry deposition of mineral dust', 'kg/m2/s', missing_value=-1.0e20)
   id_ddep_md_neg_freq  = register_tiled_diag_field(diag_name, 'md_ddep_neg_freq', (/id_ug/),  lnd%time, &
            'frequency of negative dry deposition of mineral dust', 'kg/m2/s', missing_value=-1.0e20 )

   id_ddep_noy = register_tiled_diag_field(diag_name, 'noy_ddep', &
                                          (/id_ug/),  lnd%time, 'noy dry deposition', 'mole/m2/s', &
                                          missing_value=-1.0)
   id_ddep_nhx = register_tiled_diag_field(diag_name, 'nhx_ddep', &
                                          (/id_ug/),  lnd%time, 'nhx dry deposition', 'mole/m2/s', &
                                          missing_value=-1.0)

   nbcphobic = get_tracer_index(MODEL_LAND,'bcphob')
   nbcphilic = get_tracer_index(MODEL_LAND,'bcphil')
   nomphobic = get_tracer_index(MODEL_LAND,'omphob')
   nomphilic = get_tracer_index(MODEL_LAND,'omphil')
   nsoa      = get_tracer_index(MODEL_LAND,'soa')
   nh2       = get_tracer_index(MODEL_LAND,'h2')
   nnh3      = get_tracer_index(MODEL_LAND,'nh3')
   nhno3     = get_tracer_index(MODEL_LAND,'hno3')
   nso2      = get_tracer_index(MODEL_LAND,'so2')
   nh2so4    = get_tracer_index(MODEL_LAND,'h2so4')

   ! initialize generic tracer parameters
   do tr = 1, ntcana
      call get_tracer_names(MODEL_LAND, tr, trdata(tr)%name)
      if (.not.trdata(tr)%is_generic) cycle ! skip all non-generic tracers
      trdata(tr)%tr_atm = get_tracer_index (MODEL_ATMOS, trdata(tr)%name)

      ! set up deposition flag
      trdata(tr)%do_deposition = .FALSE.
      method = ''; parameters = ''
      if (trdata(tr)%tr_atm >0) then
         if (query_method('dry_deposition', MODEL_ATMOS, trdata(tr)%tr_atm, method,parameters)) then
            trdata(tr)%do_deposition=(index(lowercase(method),'land:')>0 )
            if (trdata(tr)%do_deposition) then
               if (query_method('chem_param', MODEL_ATMOS, trdata(tr)%tr_atm, method,parameters)) then
                  if ( parse(parameters, 'nb_N_ox',  value) > 0 ) trdata(tr)%nb_N_ox    = value
                  if ( parse(parameters, 'nb_N_red',  value) > 0 ) trdata(tr)%nb_N_red  = value
                  if ( parse(parameters, 'mw',  value) > 0 )       trdata(tr)%mw        = value * 1e-3 !convert to kg/mol
               endif
            endif
         end if

         ! set up deposition parameters
         if(query_method('dry_deposition', MODEL_LAND, tr, method, parameters)) then
            if (.not. trdata(tr)%do_deposition) call error_mesg("land_tracer_driver","missmatch between atm and land configuration for "//trim(trdata(tr)%name), FATAL)
            method = trim(lowercase(method))
            trdata(tr)%param_name = method
            if (method.eq."gas") then
               trdata(tr)%parameterization = GAS_DEFAULT
               if (trdata(tr)%mw<0) call error_mesg("land_tracer_driver","mw is not defined for "//trim(trdata(tr)%name), FATAL)
            elseif (method.eq."aerosol") then
               trdata(tr)%parameterization = AEROSOL_DEFAULT
            elseif (method.eq."gas_bertagni") then
               trdata(tr)%parameterization = GAS_BERTAGNI
            elseif (method.eq."gas_codep") then
               trdata(tr)%parameterization = GAS_CODEP
               if (trdata(tr)%name.ne."nh3" .and. trdata(tr)%name.ne."so2") call error_mesg("land_tracer_driver_init", "codeposition is not implemented for "//trim(trdata(tr)%name), FATAL)
            else
               call error_mesg("land_tracer_driver_init", 'unrecognized ddep scheme for '//trim(trdata(tr)%name), FATAL)
            end if

            if ( parse(parameters, 'reactivity',  value) > 0 ) trdata(tr)%reactivity  = value
            if ( parse(parameters, 'alpha',       value) > 0 ) trdata(tr)%alpha       = value
            if ( parse(parameters, 'a_RH',        value) > 0 ) trdata(tr)%a_RH        = value

            !ratio of tracer diffusivity to h2o
            if ( parse(parameters, 'diff_ratio',  value) > 0 ) then
               trdata(tr)%diff_ratio  = value
            else
               if (trdata(tr)%mw.gt.0.) then
                  !graham law
                  trdata(tr)%diff_ratio  = sqrt(18e-3/trdata(tr)%mw)
               end if
            end if

            if (trdata(tr)%diff_ratio>0) trdata(tr)%diff_ratio23 = trdata(tr)%diff_ratio ** (2./3.)

            if (parse(parameters, 'scale_stom', value) > 0) then
               trdata(tr)%scale_stom = value
            else
               trdata(tr)%scale_stom = trdata(tr)%diff_ratio
            end if
            if ( parse(parameters, 'r_mx',  value) > 0 )   trdata(tr)%r_mx    = max(value,epsln)
            if ( parse(parameters, 'radius', value) > 0 )  trdata(tr)%radius  = value
            if ( parse(parameters, 'rho', value) > 0 )     trdata(tr)%rho     = value
         endif
      end if
   enddo

   call init_with_headers(table, trdata(:)%name)
   call add_row(table, 'do_deposition',    trdata(:)%do_deposition)
   call add_row(table, 'parameterization', trdata(:)%param_name)
   !call add_row(table, 'parameterization', trdata(:)%parameterization)
   call add_row(table, 'mw',               trdata(:)%mw)
   call add_row(table, 'alpha',            trdata(:)%alpha)
   call add_row(table, 'reactivity',       trdata(:)%reactivity)
   call add_row(table, 'radius',           trdata(:)%radius)
   call add_row(table, 'rho',              trdata(:)%rho)
   call add_row(table, 'a_RH',              trdata(:)%a_RH)
   !         call add_row(table, 'map_to',           trdata(:)%map_to)
   !         call add_row(table, 'map_to_index',     trdata(:)%map_to_index)

   call print(table,stdlog(),transposed=.TRUE.)
   call print(table,stdout(),transposed=.TRUE.)

   !to speed up calculations
   pmod_lai_dry=.TRUE.;pmod_lai_frz=.TRUE.;pmod_lai_wet=.TRUE.
   if (abs(e_lai_dry-1.0) <= 1e-5) pmod_lai_dry=.FALSE.
   if (abs(e_lai_wet-1.0) <= 1e-5) pmod_lai_wet=.FALSE.
   if (abs(e_lai_frz-1.0) <= 1e-5) pmod_lai_frz=.FALSE.

   ! register diag fields for generic tracers
   call set_default_diag_filter('land')

   id_acid_ratio = &
      register_tiled_diag_field(diag_name, 'acid_ratio', &
         (/id_ug/),  lnd%time, 'acid ratio used for codeposition', &
         'unitless', missing_value=-1.0)

   do tr = 1, ntcana
      call get_tracer_names(MODEL_LAND, tr, name, longname, units)
      call flux_units(units,funits,trdata(tr)%is_vmr)


      if ((trdata(tr)%nb_n_ox .gt. 0) .and. (trim(funits).ne.'mole/m2/s') .and. (id_ddep_nhx .gt. 0)) then
         call error_mesg("land_tracer_driver_init", trim(funits) // ' for '//trim(trdata(tr)%name) // 'incompatible with noy ddep', FATAL)
      endif
      if ((trdata(tr)%nb_n_red .gt. 0) .and. (trim(funits).ne.'mole/m2/s') .and. (id_ddep_noy .gt. 0 )) then
         call error_mesg("land_tracer_driver_init", trim(funits) // ' for '//trim(trdata(tr)%name) // 'incompatible with nhx ddep', FATAL)
      endif

      trdata(tr)%id_flux_atm = &
         register_tiled_diag_field(diag_name, trim(name)//'_flux_atm', &
         (/id_ug/),  lnd%time, trim(name)//' flux to the atmosphere', &
         trim(funits), missing_value=-1.0)
      trdata(tr)%id_dfdtr = &
         register_tiled_diag_field(diag_name, trim(name)//'_dfdtr', &
         (/id_ug/),  lnd%time,'derivative of '//trim(name)//' flux to the atmosphere', &
         trim(funits), missing_value=-1.0)
      trdata(tr)%id_conc = &
           register_tiled_diag_field(diag_name, trim(name), &
           (/id_ug/),  lnd%time, 'concentration or '//trim(name)//' in canopy air', &
           units, missing_value=-1.0)

      if (trdata(tr)%do_deposition) then
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

         trdata(tr)%id_codep = &
            register_tiled_diag_field(diag_name, trim(name)//'_codep', &
               (/id_ug/),  lnd%time, 'codeposition factor for '//trim(name), &
               'unitless', missing_value=-1.0)

         !only available for aerosols
         if (trdata(tr)%parameterization .lt. GAS_PARAM) then
            trdata(tr)%id_Eb =                register_tiled_diag_field(diag_name, trim(name)//'_Eb',     &
               (/id_ug/),  lnd%time, 'Eb collection efficiency from Brownian diffusion for '//trim(name), &
               'm/s', missing_value=-1.0)
            trdata(tr)%id_Ein =                register_tiled_diag_field(diag_name, trim(name)//'_Ein', &
               (/id_ug/),  lnd%time, 'Ein collection efficiency from interception for '//trim(name),    &
               'm/s', missing_value=-1.0)
            trdata(tr)%id_Eim =                register_tiled_diag_field(diag_name, trim(name)//'_Eim', &
               (/id_ug/),  lnd%time, 'Eim collection efficiency from impaction for '//trim(name),       &
               'm/s', missing_value=-1.0)
         end if

         trdata(tr)%id_con_g = &
         register_tiled_diag_field(diag_name, trim(name)//'_con_g', &
            (/id_ug/),  lnd%time, 'total conductance between ground and canopy air for '//trim(name), &
            'm/s', missing_value=-1.0)
         trdata(tr)%id_tcond = &
            register_tiled_diag_field(diag_name, trim(trdata(tr)%name)//'_tot_con', &
            (/id_ug/),  lnd%time,'total conductance of '//trim(trdata(tr)%name), &
            "m/s", missing_value=-1.0)
         trdata(tr)%id_tcond_new = &
            register_tiled_diag_field(diag_name, trim(trdata(tr)%name)//'_tot_con_new', &
            (/id_ug/),  lnd%time,'total conductance of '//trim(trdata(tr)%name)//' new', &
            "m/s", missing_value=-1.0)
         !these diagnostics are not relevant for aerosol tracers
         !if (.not. trdata(tr)%is_aerosol) then
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

   call set_default_diag_filter('soil')

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

   do iw=1,nwet_diag
      id_fw_wet(iw) =  register_tiled_diag_field(diag_name,'f_'//wet_str(iw), &
         (/id_ug/),  lnd%time, 'fraction of the time when canopy is more than '//wet_str(iw)(4:5)//' wet', &
         'unitless', missing_value= -1.0)

      do tr = 1, ntcana
         if (trdata(tr)%do_deposition) then
            trdata(tr)%id_tcond_wet(iw) = register_tiled_diag_field(diag_name, trim(trdata(tr)%name)//'_tot_con_'//wet_str(iw), &
               (/id_ug/),  lnd%time, 'total conductance of '//trim(trdata(tr)%name)//' with fwet>'//wet_str(iw)(4:5), &
               'm/s', missing_value=-1.0)

            if (trdata(tr)%id_tcond_wet(iw).gt.0 .and. id_fw_wet(iw).lt.0) &
               call error_mesg("land_tracer_driver",trim(trdata(tr)%name)//"_tot_con_"//wet_str(iw)//" requires"//"f_"//wet_str(iw),FATAL)
         end if
      end do
   end do

   id_fw_avg = register_tiled_diag_field(diag_name, 'fw_avg', &
      (/id_ug/),  lnd%time,'canopy avg wet fraction', &
      "unitless", missing_value=-1.0)
   id_fs_avg = register_tiled_diag_field(diag_name, 'fs_avg', &
      (/id_ug/),  lnd%time,'canopy avg snow fraction', &
      "unitless", missing_value=-1.0)
   id_fd_avg = register_tiled_diag_field(diag_name, 'fd_avg', &
      (/id_ug/),  lnd%time,'canopy avg dry fraction', &
      "unitless", missing_value=-1.0)

   id_h2_frac_water_pores_avg  = register_tiled_diag_field(diag_name, 'h2_frac_lw_pores_avg', &
      (/id_ug/),  lnd%time, 'liquid water fraction used for H2 soil removal', &
      'unitless', missing_value=-1.0)
   id_h2_frac_ice_pores_avg  = register_tiled_diag_field(diag_name, 'h2_frac_iw_pores_avg', &
      (/id_ug/),  lnd%time, 'ice water fraction used for H2 soil removal', &
      'unitless', missing_value=-1.0)
   id_h2_ilayer  = register_tiled_diag_field(diag_name, 'h2_ilayer', &
      (/id_ug/),  lnd%time, 'h2_ilayer', &
      'm', missing_value=-1.0)
   id_h2_depth_litter  = register_tiled_diag_field(diag_name, 'h2_depth_litter', &
      (/id_ug/),  lnd%time, 'h2_depth_litter', &
      'm', missing_value=-1.0)
   id_h2_soilC  = register_tiled_diag_field(diag_name, 'h2_soilC', &
      (/id_ug/),  lnd%time, 'soilC used for H2 modulation', &
      'kgC/m2', missing_value=-1.0)
   id_h2_R_bact  = register_tiled_diag_field(diag_name, 'h2_R_bact', &
      (/id_ug/),  lnd%time, 'h2_R_bact', &
      's/m', missing_value=-1.0)
   id_h2_R_inactive  = register_tiled_diag_field(diag_name, 'h2_R_inactive', &
      (/id_ug/),  lnd%time, 'h2_R_inactive', &
      's/m', missing_value=-1.0)
   id_h2_R_litter  = register_tiled_diag_field(diag_name, 'h2_R_litter', &
      (/id_ug/),  lnd%time, 'h2_R_litter', &
      's/m', missing_value=-1.0)
   id_h2_R_snow  = register_tiled_diag_field(diag_name, 'h2_R_snow', &
      (/id_ug/),  lnd%time, 'h2_R_snow', &
      's/m', missing_value=-1.0)

   id_h2_km  = register_tiled_diag_field(diag_name, 'h2_km', &
      (/id_ug/),  lnd%time, 'Km (H2)', &
      '1/s', missing_value=-1.0)
   id_h2_fm  = register_tiled_diag_field(diag_name, 'h2_fm', &
      (/id_ug/),  lnd%time, 'Moisture factor (H2)', &
      'unitless', missing_value=-1.0)
   id_h2_ft  = register_tiled_diag_field(diag_name, 'h2_ft', &
      (/id_ug/),  lnd%time, 'Temperature factor (H2)', &
      'unitless', missing_value=-1.0)

   id_h2_sdiff = register_tiled_diag_field(diag_name, 'h2_sdiff', &
      (/id_ug/),  lnd%time, 'H2 soil diffusivity', &
      'm2/s', missing_value=-1.0)

   id_h2_sopt = register_tiled_diag_field(diag_name, 'h2_s_opt', &
      (/id_ug,id_zfull/),  lnd%time,  'liquid soil moisture optimum for H2-HOB', &
      'unitless', missing_value=-1.0)
   id_h2_sws = register_tiled_diag_field(diag_name, 'h2_s_ws', &
      (/id_ug,id_zfull/),  lnd%time, 'liquid soil moisture threshold for H2-HOB', &
      'unitless', missing_value=-1.0)
   id_h2_sup = register_tiled_diag_field(diag_name, 'h2_s_up', &
      (/id_ug,id_zfull/),  lnd%time, 'liquid soil moisture maximum for H2-HOB', &
      'unitless', missing_value=-1.0)
   id_h2_moist_r1 = register_tiled_diag_field(diag_name, 'h2_moist_r1', &
      (/id_ug,id_zfull/),  lnd%time, 'soil moisture below s_ws (fraction)', &
      'unitless', missing_value=-1.0)
   id_h2_moist_r2 = register_tiled_diag_field(diag_name, 'h2_moist_r2', &
      (/id_ug,id_zfull/),  lnd%time, 'soil moisture above s_ws but below s_opt (fraction)', &
      'unitless', missing_value=-1.0)

   id_con_h2_no_litter =  register_tiled_diag_field(diag_name, 'h2_con_gr_no_litter', &
      (/id_ug/),  lnd%time, 'H2 surface conductance w/o litter', &
      'm/s', missing_value=-1.0)
   id_con_h2_no_snow =  register_tiled_diag_field(diag_name, 'h2_con_gr_no_snow', &
      (/id_ug/),  lnd%time, 'H2 surface conductance w/o snow', &
      'm/s', missing_value=-1.0)

   land_tracer_clock = mpp_clock_id( 'land_tracer', &
                           grain=CLOCK_MODULE )
   land_tracer_ddep_clock = mpp_clock_id( 'land_tracer:ddep', &
                            grain=CLOCK_MODULE )
   land_tracer_ddep_aerosol_clock = mpp_clock_id( 'land_tracer:ddep_aerosol', &
                            grain=CLOCK_MODULE )
   land_tracer_ddep_gas_clock = mpp_clock_id( 'land_tracer:ddep_gas', &
                           grain=CLOCK_MODULE )
   land_tracer_ddep_gas_vegn_clock = mpp_clock_id( 'land_tracer:ddep_gas_vegn', &
                           grain=CLOCK_MODULE )
   land_tracer_ddep_gas_grnd_clock = mpp_clock_id( 'land_tracer:ddep_gas_grnd', &
                           grain=CLOCK_MODULE )
   land_tracer_ddep_gas_h2_clock = mpp_clock_id( 'land_tracer:ddep_gas_h2', &
                           grain=CLOCK_MODULE )

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
   precip_l, precip_s, pressure, ustar, con_g, con_v_v, con_v_stem, stomatal_cond, r_bl_h2o, con_atm, &
   ! output
   dep_to_snow )

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
   ! output
   real, intent(out):: dep_to_snow(N_SNOW_TRACERS) ! dry deposition of snow tracers, kg/m2/s

   integer :: tr      ! tracer index
   integer :: k       ! cohort index
   integer :: iw      ! wet threshold index
   real    :: rho     ! density of canopy air
   real    :: dq      ! canopy air tracer tendency per time step
   real    :: con_st_tr ! total stomatal conductance, scaled by dry leaf area, m/s
   real    :: con_mx  ! "mesophyll conductance", m/s
   real    :: con_mx_st !mesophyll+stomatal
   real    :: con_cu  ! total cuticular conductance, including dry and wet areas, m/s
   real    :: con_gr  ! "ground conductance", m/s
   real    :: cv, cg  ! total conductances for vegetation and ground surface, m/s
   real    :: f_atm   ! flux of the tracer to the atmosphere, kg/(m2 s)
   real    :: ddep    ! dry deposition of the tracer, kg/(m2 s)
   real    :: emis(ntcana) ! tracer sources
   real    :: gamma_codep(ntcana) !codeposition modulation
   real    ::  ft, & ! fraction of canopy not covered by intercepted water/snow
               fw, & ! fraction of canopy covered by intercepted water
               fs    ! fraction of canopy covered by intercepted snow

   !fraction of ground that is frozen, wet, dry
   real    :: gfrac_frz,gfrac_wet,gfrac_dry

   real    :: con_cu_dry, con_cu_wet, con_cu_frz, con_stem
   real    :: con_gr_dry, con_gr_wet, con_gr_frz
   real    :: frac_desert
   real    :: con_v_v_tr, con_v_stem_tr, con_bl_tr, con_g_tr

   real    :: con_cu_diag, con_mx_st_diag, con_stem_diag, fdiag
   real    :: con_v_v_tr_diag, con_v_stem_tr_diag
   real    :: econ_cu, econ_stem, econ_mx_st
   real    :: econ_cu_dry, econ_cu_wet, econ_cu_frz
   real    :: dvel, dvel_new
   real    :: tmp

   real    :: alpha_aer, gamma_aer, A_aer, cg_aer_v, cg_aer_g
   real    :: fw_avg, fs_avg, rh

   real    :: ddep_oa, ddep_bc, ddep_md, ddep_noy, ddep_nhx
   real    :: ustar_mod, tcond

   real    :: e_RH, acid_ratio, acid, base, ustar_s

   call mpp_clock_begin (land_tracer_clock)

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
              precip_l, precip_s, pressure, ustar, con_g, con_v_v, &
              ddep_md )
   ! wind10 is not passed to this subroutine yet

   call mpp_clock_begin (land_tracer_ddep_clock)
   if (do_generic_ddep_emis) then

      ! update generic tracers
      ! calculate tracers sources
      emis(:) = 0.0
      ! TODO: add non-zero sources for generic tracers

      ddep_oa = 0. ; ddep_bc = 0. ; ddep_nhx = 0. ; ddep_noy = 0.

      call get_tile_property(tile,gfrac_frz,gfrac_wet,frac_desert)
      gfrac_dry = max(1.-gfrac_wet-gfrac_frz,0.)

      !calculate rh
      call check_temp_range(tile%cana%T, 'update_cana_tracers', 'cana_T')
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

      !precalculate a few quantities
      ustar_mod = ustar**e_ustar

      !calculate acid ratio. Here I follow the definition of Masad (2010) 2*SO2  and I added SO4 and NH4. Note that EMEP used [SO2]/[NH3] i.e., 1/2 the acid/base ratio
      acid = 0.; base = 0.
      if (nhno3.gt.0)   acid = acid + tile%cana%tr(nhno3)
      if (nh2so4.gt.0)  acid = acid + 2.*tile%cana%tr(nh2so4)
      if (nso2.gt.0)    acid = acid + tile%cana%tr(nso2)  !note that Massad x2
      if (nnh3.gt.0)    base = base + tile%cana%tr(nnh3)

      !
      acid_ratio         = max(acid,epsln)/max(base,epsln)
      !to avoid very large changes cap acid ratio
      acid_ratio         = min(max(0.25,acid_ratio),4.)

      gamma_codep(:)     = 1.
      !from Simpson (2012). Note that Simpson defines the acid ratio as SO2/NH3, while we use 2.*SO2/NH3.
      !Force gamma_code(nso2) to be 1 when acid_ratio = 1 (near neutral conditions)
      if (nso2.gt.0.) gamma_codep(nso2) = exp(-1.1*acid_ratio)*3.0

      !from Massad (2010)
      if (nnh3.gt.0.) gamma_codep(nnh3)  = acid_ratio !Massad has no observations above 3.
      !Note that there are large differences in the treatment of RH
      !Zhang (2003) suggests that Rcut = Rcut_ref * exp(-0.03*RH)
      !Massad (2010)  suggests that Rcut = Rcut_ref *  exp(0.18*(100-RH)) = R_cut_ref_p * exp(-[0.03--0.18]*RH)
      !Based on a 0.03, R_cut_ref_p = R_cut_ref*20 or 31*.5*20 = 630. For forest, Zhang proposed R_cut ~= 2000.
      !The two expression will thus be equal for an AR of 0.315


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

            if (trdata(tr)%parameterization.gt.GAS_PARAM) then

               call mpp_clock_begin (land_tracer_ddep_gas_clock)

               e_RH      = exp(RH*trdata(tr)%a_RH)

               !conductance to the vegetation
               call mpp_clock_begin (land_tracer_ddep_gas_vegn_clock)
               if (associated(tile%vegn)) then
                  do k = 1, tile%vegn%n_cohorts
                     associate(c=>tile%vegn%cohorts(k),sp=>spdata(tile%vegn%cohorts(k)%species))
                        call get_vegn_wet_frac ( c, fw=fw, fs=fs ); ft = 1-fw-fs

                        con_cu_dry  = ft * ustar_mod * get_conductance_tracer(trdata(tr),sp%r_cus,sp%r_cuo)  / scale_r_T(c%Tv,c_dry) * e_RH
                        if (pmod_lai_dry) then
                           con_cu_dry  = con_cu_dry* c%lai**e_lai_dry
                        else
                           con_cu_dry = con_cu_dry * c%lai
                        end if

                        con_cu_wet  = fw * ustar_mod * get_conductance_tracer(trdata(tr),sp%r_cus_wet,sp%r_cuo_wet) / scale_r_T(c%Tv,c_wet)
                        if (pmod_lai_wet) then
                           con_cu_wet  = con_cu_wet* c%lai**e_lai_wet
                        else
                           con_cu_wet = con_cu_wet * c%lai
                        end if

                        if (trdata(tr)%parameterization.eq.GAS_CODEP) then
                           con_cu_dry = con_cu_dry * gamma_codep(tr)
                           con_cu_wet = con_cu_wet * gamma_codep(tr)
                        end if

                        con_cu_frz  = fs * get_conductance_tracer(trdata(tr),get_snows(c%Tv),r_snowo)

                        if (pmod_lai_frz) then
                           con_cu_frz  = con_cu_frz* c%lai**e_lai_frz
                        else
                           con_cu_frz = con_cu_frz * c%lai
                        end if

                        !here we use the bulk leaf property for the cohort. This is different from the LM3 implementation.
                        con_cu   = con_cu_dry+con_cu_wet+con_cu_frz

                        !comment-out temperature dependence based on Clifton (2020)
                        !con_stem = c%sai * get_conductance_tracer(trdata(tr),sp%r_stems,sp%r_stemo) / scale_r_T(c%Tv)
                        con_stem =  c%sai * get_conductance_tracer(trdata(tr),sp%r_stems,sp%r_stemo)

                        if (trdata(tr)%parameterization.eq.GAS_CODEP) then
                           con_stem = con_stem * gamma_codep(tr)
                        end if

                        con_st_tr    = stomatal_cond(k) * trdata(tr)%scale_stom

                        if (trdata(tr)%r_mx .gt. 0) then
                           con_mx    = 1./trdata(tr)%r_mx
                        else
                           con_mx    = c%lai * get_conductance_tracer(trdata(tr),1.,100.)
                        end if

                        con_mx_st = conductance_series(con_mx,con_st_tr)

                        !calculate contribution of this cohort to the overall vegetation conductance
                        !con_v_v and con_stem are for H2O, we need to scale by (Di/Dw)**(2./3.)
                        con_v_v_tr     = con_v_v(k)*trdata(tr)%diff_ratio23
                        con_v_stem_tr  = con_v_stem(k)*trdata(tr)%diff_ratio23

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

                        con_v_v_tr_diag     = con_v_v_tr_diag    + c%layerfrac*con_v_v_tr
                        con_v_stem_tr_diag  = con_v_stem_tr_diag + c%layerfrac*con_v_stem_tr
                     end associate
                  end do
                  cv = econ_mx_st + econ_cu + econ_stem
               end if
               call mpp_clock_end (land_tracer_ddep_gas_vegn_clock)

               !note that for the ground we are calculating the tile average
               call mpp_clock_begin (land_tracer_ddep_gas_grnd_clock)
               if (tr==nh2 .and. trdata(tr)%parameterization.eq.GAS_BERTAGNI) then
                  !for now set constant conductance
                  call mpp_clock_begin (land_tracer_ddep_gas_h2_clock)
                  con_gr_dry = con_h2(tile,pressure)
                  con_gr_wet = 0.
                  con_gr_frz = 0.
                  call mpp_clock_end (land_tracer_ddep_gas_h2_clock)
               else
                  con_gr_dry =    gfrac_dry * get_conductance_tracer(trdata(tr),r_gs_dry,r_go_dry) * 1./scale_r_T(land_tile_grnd_T(tile),c_dry) * 1./scale_biomass(frac_desert)
                  con_gr_wet =    gfrac_wet * get_conductance_tracer(trdata(tr),r_gs_wet,r_go_wet) * 1./scale_r_T(land_tile_grnd_T(tile),c_wet)
                  con_gr_frz =    gfrac_frz * get_conductance_tracer(trdata(tr),get_snows(land_tile_grnd_T(tile)),r_snowo) * 1./scale_r_T(land_tile_grnd_T(tile),c_snow)

               !if (trdata(tr)%parameterization.eq.GAS_CODEP) then
               !   con_gr_dry = con_gr_dry !* gamma_codep(tr) - do not apply correction to ground
               !   con_gr_wet = con_gr_wet !* gamma_codep(tr)
               !end if

                  if (associated(tile%lake)) then
                     if (tile%lake%ws(1).le.ws_min) then
                        con_gr_wet = get_conductance_tracer(trdata(tr),r_gs_lake,r_go_lake)
                     end if
                  end if
               end if
               call mpp_clock_end( land_tracer_ddep_gas_grnd_clock )

               con_bl_tr  = 1./(r_bl_h2o+epsln) * trdata(tr)%diff_ratio23

               con_g_tr   = conductance_series(con_g,con_bl_tr)

               con_gr     = con_gr_dry+con_gr_wet+con_gr_frz
               cg         = conductance_series(con_gr,con_g_tr)

               call send_tile_data(trdata(tr)%id_con_mx_st,con_mx_st_diag, tile%diag)
               call send_tile_data(trdata(tr)%id_con_cu,con_cu_diag, tile%diag)
               call send_tile_data(trdata(tr)%id_con_stem,con_stem_diag, tile%diag)
               call send_tile_data(trdata(tr)%id_con_gr,con_gr, tile%diag)
               call send_tile_data(trdata(tr)%id_con_v_v, con_v_v_tr_diag, tile%diag)
               call send_tile_data(trdata(tr)%id_con_v_stem, con_v_stem_tr_diag,  tile%diag)
               call send_tile_data(trdata(tr)%id_con_v_g, con_g_tr, tile%diag)

               call mpp_clock_end(land_tracer_ddep_gas_clock)

            elseif (trdata(tr)%parameterization.lt.GAS_PARAM) then

               call mpp_clock_begin(land_tracer_ddep_aerosol_clock)
               !aerosol
               if (associated(tile%vegn)) then
                  cv = 0
                  do k = 1, tile%vegn%n_cohorts
                     associate(c=>tile%vegn%cohorts(k),sp=>spdata(tile%vegn%cohorts(k)%species))

                        call get_vegn_wet_frac ( c, fw=fw, fs=fs ); ft = 1-fw-fs
                        cg_aer_v = cg_aer(trdata(tr),                 &
                           tile%cana%T,ustar,pressure,              &
                           sp%alpha_aer,                            &
                           sp%gamma_aer,                            &
                           sp%A_aer,                                &
                           ft,fw, fs)

                        if (b_lai_aer.gt.epsln) then
                           cg_aer_v = cg_aer_v*c%lai**b_lai_aer
                        end if

                        cv = cv + c%layerfrac*conductance_series(con_v_v(k),cg_aer_v)

                     end associate
                  end do


               else
                  cv = 0.
               end if

               if (associated(tile%glac)) then
                  A_aer       = -999.
                  gamma_aer   = gamma_aer_frz
                  alpha_aer   = alpha_aer_frz
               elseif (associated(tile%lake)) then
                  A_aer       = -999.
                  gamma_aer   = gamma_aer_lake
                  alpha_aer   = alpha_aer_lake
               else
                  A_aer       = -999.
                  gamma_aer   = gamma_aer_desert
                  alpha_aer   = alpha_aer_desert
               end if

               ustar_s=ustar*exp(-10.*tile%land_d) !same as dust
               cg_aer_g =  cg_aer(trdata(tr),land_tile_grnd_T(tile),ustar_s,pressure, &
                                 alpha_aer,gamma_aer,A_aer,gfrac_dry,gfrac_wet,gfrac_frz)

               cg = conductance_series(con_g,cg_aer_g)

               call mpp_clock_end(land_tracer_ddep_aerosol_clock)

            endif
            end if

            rho = pressure/(rdgas*tile%cana%T *(1+d608*tile%cana%tr(isphum)))

            if (tile%cana%tr(tr).lt.0.) then
               tcond = 0
            else
               tcond = (cv+cg)
            end if

            ! update tracer concentration -- this needs to be done even when do_deposition
            ! is FALSE, because source might be non-zero
            dq = (-tr_flux(tr)-rho*tcond*tile%cana%tr(tr)+emis(tr)) &
               / (canopy_air_mass_for_tracers/dt + dfdtr(tr) + rho*tcond)

            if (is_watch_point()) then
               write(*,*) 'update_cana_tracers : ', trim(trdata(tr)%name)
               __DEBUG4__(tile%cana%tr(tr), tr_flux(tr), dfdtr(tr), emis(tr))
               __DEBUG4__(rho,cv,cg,tcond)
               __DEBUG2__(dq,tile%cana%tr(tr) + dq)
            endif

            tile%cana%tr(tr) = tile%cana%tr(tr) + dq
            ! ---- final values of the fluxes, for diagnostics
            ddep  = rho*tcond*tile%cana%tr(tr)
            f_atm = (tr_flux(tr)+dfdtr(tr)*dq)

            if (trdata(tr)%is_vmr) then
               tmp   = 1e-3*WTMAIR*WTMH2O/((1.-tile%cana%tr(isphum))*WTMH2O+tile%cana%tr(isphum)*WTMAIR) !(kg(air)/mol(air))
               ddep  = ddep/tmp
               f_atm = f_atm/tmp
            end if



            if (trdata(tr)%do_deposition) then
            ! ---- diagnostic section
            if (con_atm.gt.epsln) then
               dvel = con_atm*tcond/(con_atm+tcond)
            else
               dvel = 0.
            end if
            if (trdata(tr)%id_tcond_new>0) then
               if (dfdtr(tr).gt.epsln) then
                  dvel_new =  dfdtr(tr)*tcond/(dfdtr(tr)+rho*tcond)
               else
                  dvel_new = 0.
               end if
               call send_tile_data(trdata(tr)%id_tcond_new,dvel_new, tile%diag)
            end if

            call send_tile_data(trdata(tr)%id_con_v,      cv,                tile%diag)
            call send_tile_data(trdata(tr)%id_con_g,      cg,                tile%diag)
            call send_tile_data(trdata(tr)%id_emis,       emis(tr),          tile%diag)
            call send_tile_data(trdata(tr)%id_ddep,       ddep,              tile%diag)
            call send_tile_data(trdata(tr)%id_flux_atm,   f_atm,             tile%diag)
            call send_tile_data(trdata(tr)%id_tcond,      dvel,              tile%diag)
            call send_tile_data(trdata(tr)%id_codep,      gamma_codep(tr),   tile%diag)


            !save temporary array for
            if (id_ddep_noy.gt.0 .and. trdata(tr)%nb_n_ox.gt.0)  ddep_noy = ddep_noy + ddep*trdata(tr)%nb_n_ox
            if (id_ddep_nhx.gt.0 .and. trdata(tr)%nb_n_red.gt.0) ddep_nhx = ddep_nhx + ddep*trdata(tr)%nb_n_red

            if (tr .eq. nomphilic) ddep_oa  = ddep_oa + ddep
            if (tr .eq. nomphobic) ddep_oa  = ddep_oa + ddep
            if (tr .eq. nsoa)      ddep_oa  = ddep_oa + ddep

            if (tr .eq. nbcphilic) ddep_bc  = ddep_bc + ddep
            if (tr .eq. nbcphobic) ddep_bc  = ddep_bc + ddep

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
            if (trdata(tr)%id_econ_g>0)  call send_tile_data(trdata(tr)%id_econ_g, (1.-fdiag)*dvel,tile%diag)
            if (trdata(tr)%id_econ_v>0)  call send_tile_data(trdata(tr)%id_econ_v, fdiag*dvel,     tile%diag)
            if (trdata(tr)%id_ddep_g>0)  call send_tile_data(trdata(tr)%id_ddep_g, (1.-fdiag)*ddep,tile%diag)
            if (trdata(tr)%id_ddep_v>0)  call send_tile_data(trdata(tr)%id_ddep_v, fdiag*ddep,     tile%diag)


            !the following diagnostics are not (yet) defined if the tracer is an aerosol species
            if (trdata(tr)%parameterization .gt. GAS_PARAM) then
               if (trdata(tr)%id_econ_g_wet>0)  call send_tile_data(trdata(tr)%id_econ_g_wet,   con_gr_wet/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*dvel,     tile%diag)
               if (trdata(tr)%id_econ_g_dry>0)  call send_tile_data(trdata(tr)%id_econ_g_dry,   con_gr_dry/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*dvel,     tile%diag)
               if (trdata(tr)%id_econ_g_frz>0)  call send_tile_data(trdata(tr)%id_econ_g_frz,   con_gr_frz/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*dvel,     tile%diag)

               if (trdata(tr)%id_ddep_g_wet>0)  call send_tile_data(trdata(tr)%id_ddep_g_wet,   con_gr_wet/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*ddep,     tile%diag)
               if (trdata(tr)%id_ddep_g_dry>0)  call send_tile_data(trdata(tr)%id_ddep_g_dry,   con_gr_dry/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*ddep,     tile%diag)
               if (trdata(tr)%id_ddep_g_frz>0)  call send_tile_data(trdata(tr)%id_ddep_g_frz,   con_gr_frz/(con_gr_wet+con_gr_dry+con_gr_frz+epsln)*(1.-fdiag)*ddep,     tile%diag)

               if (trdata(tr)%id_econ_cu>0)     call send_tile_data(trdata(tr)%id_econ_cu,       fdiag*dvel*econ_cu/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
               if (trdata(tr)%id_econ_cu_wet>0) call send_tile_data(trdata(tr)%id_econ_cu_wet,   fdiag*dvel*econ_cu_wet/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
               if (trdata(tr)%id_econ_cu_dry>0) call send_tile_data(trdata(tr)%id_econ_cu_dry,   fdiag*dvel*econ_cu_dry/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
               if (trdata(tr)%id_econ_cu_frz>0) call send_tile_data(trdata(tr)%id_econ_cu_frz,   fdiag*dvel*econ_cu_frz/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)

               if (trdata(tr)%id_econ_stem>0)   call send_tile_data(trdata(tr)%id_econ_stem, fdiag*dvel*econ_stem/(econ_cu+econ_stem+econ_mx_st+epsln),   tile%diag)
               if (trdata(tr)%id_econ_stom>0)   call send_tile_data(trdata(tr)%id_econ_stom, fdiag*dvel*econ_mx_st/(econ_cu+econ_stem+econ_mx_st+epsln),  tile%diag)

               if (trdata(tr)%id_ddep_cu>0)     call send_tile_data(trdata(tr)%id_ddep_cu,       fdiag*ddep*econ_cu/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
               if (trdata(tr)%id_ddep_cu_wet>0) call send_tile_data(trdata(tr)%id_ddep_cu_wet,   fdiag*ddep*econ_cu_wet/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
               if (trdata(tr)%id_ddep_cu_dry>0) call send_tile_data(trdata(tr)%id_ddep_cu_dry,   fdiag*ddep*econ_cu_dry/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)
               if (trdata(tr)%id_ddep_cu_frz>0) call send_tile_data(trdata(tr)%id_ddep_cu_frz,   fdiag*ddep*econ_cu_frz/(econ_cu+econ_stem+econ_mx_st+epsln),     tile%diag)

               if (trdata(tr)%id_ddep_stem>0) call send_tile_data(trdata(tr)%id_ddep_stem, fdiag*ddep*econ_stem/(econ_cu+econ_stem+econ_mx_st+epsln),   tile%diag)
               if (trdata(tr)%id_ddep_stom>0) call send_tile_data(trdata(tr)%id_ddep_stom, fdiag*ddep*econ_mx_st/(econ_cu+econ_stem+econ_mx_st+epsln),  tile%diag)
            end if
         end if

      end do

      ! send the deposition to the diagnostics and remove negatives
      call diag_ddep(ddep_bc, id_ddep_bc, id_ddep_bc_neg, id_ddep_bc_neg_freq, tile%diag)
      call diag_ddep(ddep_oa, id_ddep_oa, id_ddep_oa_neg, id_ddep_oa_neg_freq, tile%diag)
      call diag_ddep(ddep_md, id_ddep_md, id_ddep_md_neg, id_ddep_md_neg_freq, tile%diag)

      ! set up dry deposition of light-absorbing particles to snow
      dep_to_snow(:)          = 0.0
      dep_to_snow(SNOW_TR_BC) = ddep_bc
      dep_to_snow(SNOW_TR_OM) = ddep_oa
      dep_to_snow(SNOW_TR_MD) = ddep_md

   !    call check_var_range(ddep_bc, 0.0, HUGE(1.0), 'update_cana_tracers', 'ddep_bc', WARNING)
   !    call check_var_range(ddep_oa, 0.0, HUGE(1.0), 'update_cana_tracers', 'ddep_oa', WARNING)
   !    call check_var_range(ddep_md, 0.0, HUGE(1.0), 'update_cana_tracers', 'ddep_md', WARNING)

   !    call send_tile_data(id_ddep_bc,  ddep_bc,  tile%diag)
   !    call send_tile_data(id_ddep_oa,  ddep_oa,  tile%diag)
      call send_tile_data(id_ddep_noy, ddep_noy, tile%diag)
      call send_tile_data(id_ddep_nhx, ddep_nhx, tile%diag)

      call send_tile_data(id_con_atm,      con_atm,      tile%diag)
      call send_tile_data(id_gfrac_dry,    gfrac_dry,    tile%diag)
      call send_tile_data(id_gfrac_wet,    gfrac_wet,    tile%diag)
      call send_tile_data(id_gfrac_frz,    gfrac_frz,    tile%diag)
      call send_tile_data(id_frac_desert,  frac_desert,  tile%diag)

      call send_tile_data(id_acid_ratio, acid_ratio, tile%diag)

   end if

   ! send concentrations for all tracers, generic or not
   do tr = 1, ntcana
      call send_tile_data(trdata(tr)%id_conc,       tile%cana%tr(tr), tile%diag)
      call send_tile_data(trdata(tr)%id_dfdtr,      dfdtr(tr),  tile%diag)
   enddo

   call mpp_clock_end(land_tracer_ddep_clock)
   call mpp_clock_end(land_tracer_clock)

contains
   subroutine diag_ddep(ddep_diag, id_ddep, id_negatives, id_negative_freq, diag)
      real,    intent(in) :: ddep_diag        ! dry deposition, updated to be positive
      integer, intent(in) :: id_ddep          ! diag ID of dry deposition, prior to update (that is, including negatives)
      integer, intent(in) :: id_negatives     ! diag ID of the negative deposit, to keep track of average negative deposition
      integer, intent(in) :: id_negative_freq ! diag ID of the frequency of negatives
      type(diag_buff_type), intent(inout) :: diag ! diagnostic buffer of the tile

      real :: freq

      call send_tile_data(id_ddep, ddep_diag, diag)
      call send_tile_data(id_negatives,  min(ddep_diag,0.0), diag)
      if (ddep_diag<0.0) then
         freq=1.0
      else
         freq=0.0
      endif
      call send_tile_data(id_negative_freq, freq, diag)
   end subroutine

end subroutine update_cana_tracers

subroutine get_tile_property(tile,gfrac_frz,gfrac_wet,frac_desert)

   type(land_tile_type), intent(in)   :: tile

   real,    intent(out)             :: gfrac_frz,gfrac_wet,frac_desert

   real                             :: biomass
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
      call tile%snow%get_depth_area(snow_depth, snow_area)
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
subroutine flux_units(tracer_units,units,is_vmr)
   character(*), intent(in)   :: tracer_units
   character(32), intent(out) :: units
   logical, intent(out)       :: is_vmr

   select case (trim(lowercase(tracer_units)))
   case ('mmr')
      units = 'kg/m2/s'
      is_vmr = .FALSE.
   case ('kg/kg')
      units = 'kg/m2/s'
      is_vmr = .FALSE.
   case ('vmr')
      units = 'mole/m2/s'
      is_vmr = .TRUE.
   case ('mol/mol')
      units = 'mole/m2/s'
      is_vmr = .TRUE.
   case ('mole/mole')
      units = 'mole/m2/s'
      is_vmr = .TRUE.
   case default
      units = trim(tracer_units)//' kg/(m2 s)'
      is_vmr = .FALSE.
   end select
end subroutine flux_units

elemental real function conductance_series(con1,con2) result(con)
   real, intent(in) :: con1, con2

   con = (con1*con2)/(con1+con2+epsln)

end function conductance_series

elemental real function get_snows(T) result(s)
   !erisman (1994) showed that resistance increases as temperature decreases from 70 to 500
   real, intent(in) :: T

   s = max(min(r_snows*(275.15-T),r_snows_max),r_snows_min)

end function get_snows

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

function  cg_aer(tr_data,T,ustar,pressure,alpha,gamma,A,frac_dry,frac_wet,frac_snow) result(con)

   !from Zhang 2001 as presented by Seinfeld (19.27)
   !rb = 1/ (3*ustar * (Sc^-gamma + (St/(alpha+St))^2 + 1/2 * (Dp/A)^2) * R1 )

   type(tracer_data_type), intent(in) :: tr_data

   real, intent(in) :: T !temperature
   real, intent(in) :: alpha, gamma, A !depend on land type
   real, intent(in) :: pressure
   real, intent(in) :: ustar
   real, intent(in) :: frac_wet, frac_dry, frac_snow


   real :: con
   real :: Ein,Eim,Eb

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
      St       = vts*ustar/(grav*A)
      Ein      = (0.5*(2*rp/A)**2)
   else
      St       = vts*ustar**2/(grav*kvis)
      Ein      = 0.
   end if

   R1      = frac_wet + frac_dry*exp(-St**0.5)

   !brownian diffusion
   diff_aer= kb*T*C_c / ( 6.*pi*dvis*rp )
   Sc      = kvis/diff_aer

   Eb      = 1./Sc**gamma
   Eim     = (St/(alpha+St))**beta

   con = (e0 * ustar * R1 * (Eb + Eim + Ein))
   if (cg_aer_frz.gt.0.)  then
     con     = (frac_dry+frac_wet) * con + frac_snow*cg_aer_frz
                end if

end function cg_aer

real function con_h2(tile,p) result(con)

   type(land_tile_type),   intent(inout) :: tile
   real,                   intent(in) :: p

   integer :: isoil
   real    :: frac_water_pores_avg, frac_ice_pores_avg
   real,dimension(num_l) :: frac_water_pores, frac_ice_pores, soil_C_layer, s_ws, s_opt, s_upc, psi_sat_ref, h2_moist_r1, h2_moist_r2
   real    :: dz, T_avg, T_avgC, dz_tot
   real    :: diff_h2
   real    :: f_T, f_M
   real    :: snow_depth, snow_area, Tsnow

   real    :: R_inactive, R_snow, R_bact, inactive_layer, R_litter, R_litter_leaf

   real    :: beta1, beta2, norm
   real    :: litterC_leaf
   real    :: depth_litter_leaf, depth_litter
   real    :: grnd_T

   real    :: h2_km_eff
   real    :: s_opt_avg, s_upc_avg, s_ws_avg
   real    :: soil_C, litt_C(N_LITTER_POOLS)

   real, parameter :: s_up = 1. !no cap on h2 activity

   logical :: TOP_LAYER,FLOODED

   con = 0.

   if (associated(tile%soil)) then

      frac_water_pores  = max(min(soil_theta(tile%soil),1.),0.)
      frac_ice_pores    = max(min(soil_ice_porosity(tile%soil),1.),0.)

      isoil           = 1

      h2_km_eff       = h2_km
      f_M             = 0.
      f_T             = 0.
      diff_H2         = 0.
      R_bact          = 1.e20
      R_litter        = 0.
      R_litter_leaf   = 0.
      R_snow          = 0.
      R_inactive      = 0.
      inactive_layer  = 0.
      depth_litter    = 0.
      depth_litter_leaf = 0.

      s_opt_avg       = 0.
      s_upc_avg       = 0.
      s_ws_avg        = 0.

      T_avg                = 0.
      frac_water_pores_avg = 0.
      frac_ice_pores_avg   = 0.

      dz_tot          = 0.
      soil_C          = 0.

      s_ws(:)         = -1.
      s_opt(:)        = -1.
      psi_sat_ref(:)  = -1.

      h2_moist_r1(:) = -1.
      h2_moist_r2(:) = -1.

      s_ws(:)  = -1
      s_opt(:) = -1
      s_upc(:) = -1

      beta1       = h2_beta1

      if (h2_soilC_mod_id .gt. 0) then
         !get C for modulation do not include litter
         soil_C = tile%soilc%total_soil_C_to_depth(h2_depth)
      end if

      TOP_LAYER = .TRUE.
      isoil     = 1

      !loop over soil layers
      do while (zhalf(isoil).lt.h2_depth)
         dz = min(zhalf(isoil+1),h2_depth)-zhalf(isoil)

         !calculate soil properties
         !could be optimized as it only needs to be calculated once
         psi_sat_ref(isoil) = tile%soil%pars%psi_sat_ref/tile%soil%alpha(isoil) !m

         if (frac_ice_pores(isoil).gt.0) psi_sat_ref(isoil) = psi_sat_ref(isoil) /2.2

         !soil activation threshold
         s_ws(isoil)      = min(max((psi_sat_ref(isoil)/h2_psi_ws)**(1./tile%soil%pars%chb),0.),1.)
         s_opt(isoil)     = min(max((psi_sat_ref(isoil)/h2_psi_opt)**(1./tile%soil%pars%chb),0.),1.)
         s_upc(isoil)     = 1. !min(max(s_up - frac_ice_pores(isoil),0.),1.)

         !NOTE that this equation is only valid between if Xl_eff>psi_min and <Xsat. This is ok as long as psi_h2 is >-100e2 m

         h2_moist_r1(isoil) = 0.
         h2_moist_r2(isoil) = 0.
         if (frac_water_pores(isoil)  .lt. s_ws(isoil)) h2_moist_r1(isoil) = 1.
         if ((frac_water_pores(isoil) .ge. s_ws(isoil)) .and. (frac_water_pores(isoil).lt.s_opt(isoil))) h2_moist_r2(isoil) = 1.

         flooded = .FALSE.

         if (frac_water_pores(isoil).lt.s_ws(isoil) .and. TOP_LAYER) then
            !we have yet to encounter a wet enough layer
            R_inactive   = R_inactive +                                                   &
                           dz/max(diff_H2_soil(tile%soil%T(isoil),p,                      &
                                  tile%soil%pars%vwc_sat,                                 &
                                  frac_water_pores(isoil)+frac_ice_pores(isoil),          &
                                  tile%soil%pars%chb),1.e-20)
            inactive_layer = inactive_layer + dz
         elseif (((1.-frac_ice_pores(isoil)+frac_water_pores(isoil)).lt.epsln) .and. TOP_LAYER) then
            R_inactive = 1.e20 !H2 won't diffusive to active sites
            flooded    = .TRUE.
         else
            !this layer exceeds minimum s_ws. Add all layers below until h2_depth
            T_avg                = T_avg                + dz*tile%soil%T(isoil)
            frac_water_pores_avg = frac_water_pores_avg + dz*frac_water_pores(isoil)
            frac_ice_pores_avg   = frac_ice_pores_avg   + dz*frac_ice_pores(isoil)

            s_opt_avg            = s_opt_avg + dz*s_opt(isoil)
            s_upc_avg            = s_upc_avg + dz*s_upc(isoil)
            s_ws_avg             = s_ws_avg  + dz*s_ws(isoil)

            dz_tot               = dz_tot               + dz

            TOP_LAYER = .FALSE.
         end if
         isoil                = isoil+1
      end do

      R_snow = 0.
      if (associated(tile%snow)) then
         call tile%snow%get_depth_area(snow_depth, snow_area)
         if (snow_depth.gt.epsln) then
            Tsnow  = tile%snow%sfc_temp()
            R_snow = snow_depth/diff_H2_snow(Tsnow,p,tile%snow%porosity())
         end if
      end if

      !dz_tot>0. H2 uptake is possible
      if (dz_tot.gt.epsln) then
         T_avg                 = T_avg/dz_tot
         T_avgC                = T_avg - 273.15
         frac_water_pores_avg  = frac_water_pores_avg/dz_tot
         frac_ice_pores_avg    = frac_ice_pores_avg/dz_tot

         s_ws_avg              = s_ws_avg/dz_tot
         s_opt_avg             = s_opt_avg/dz_tot
         s_upc_avg             = s_upc_avg/dz_tot

         !s = min(max((psi_sat_ref/psi)**(1./b),0.),1.)
         !derivative : beta1*(s-s_ws)**(beta1-1)* (s_up-s)**beta2 - beta2* (s-s_ws)**beta1 * (s_up-s)**(beta2-1)
         !=> beta1*(s_up-s_opt) - beta2*(s_opt-s_ws) = 0
         !=> beta2 = beta1*(s_up-s_opt)/(s_opt-s_ws)
         beta2     = beta1 * (s_upc_avg-s_opt_avg)/(s_opt_avg-s_ws_avg)
         norm      = (s_opt_avg-s_ws_avg)**beta1*(s_upc_avg-s_opt_avg)**beta2

         if (frac_water_pores_avg.gt.s_upc_avg) then
            f_M = 0.
         elseif (frac_water_pores_avg.lt.s_ws_avg) then
            f_M = 0.
         else
            f_M = 1/norm*(frac_water_pores_avg-s_ws_avg)**beta1*(s_upc_avg-frac_water_pores_avg)**beta2
         end if

         if (f_M.lt.0.) then
            write(*,*) norm,frac_water_pores_avg,s_ws_avg,beta1,s_upc_avg,frac_water_pores_avg,beta2
            call error_mesg("land_tracer_driver","f_M<0",FATAL)
         end if

         if (f_M.gt.1.01) then
            write(*,*) 'f_M,norm,frac_water_pores_avg,s_ws,beta1,s_opt,beta2',f_M,norm,frac_water_pores_avg,s_ws_avg,beta1,s_opt_avg,beta2
            call error_mesg("land_tracer_driver","f_M>1",FATAL)
         end if

         !From Ehhalt (2011)
         f_T = 1/(1+exp(-(T_avgC-3.8)/6.7)) + 1./(1.+exp((T_avgC - 62.2)/7.7)) - 1.

         diff_H2 = diff_H2_soil( T_avg,                                            &
                                 p,                                                &
                                 tile%soil%pars%vwc_sat,                           &
                                 frac_ice_pores_avg+frac_water_pores_avg,          &
                                 tile%soil%pars%chb)

         if (h2_soilC_mod_id.eq.H2_SOILC_PAULOT21) then
            h2_km_eff = (h2_km*soil_C/h2_depth)/(soil_C/h2_depth+h2_soilC_param(1))
         elseif (h2_soilC_mod_id.eq.H2_SOILC_REJI25) then
            !assume soil density of 2650 kg/m3 -> convert to %
            h2_km_eff = h2_km*(h2_soilC_param(1)*soil_C/h2_depth*1./2650.*100.+h2_soilC_param(2))
         elseif (h2_soilC_mod_id.eq.H2_SOILC_NO_MOD) then
            h2_km_eff = h2_km
         end if

         if (h2_precip_min.gt.0. .and. associated(tile%vegn)) then
            !this is a short-term solution to reduce h2 uptake in deserts.
            !scale h2_km_eff linearly from 0 to 1. between 0 and h2_precip_min (kg/m2)
            if (tile%vegn%p_ann.lt.h2_precip_min) then
               h2_km_eff = h2_km_eff * max(tile%vegn%p_ann,0.)/h2_precip_min
            end if
         end if

         !regardless, we turn off h2 uptake if there is no soilC
         if (soil_C.le.epsln) h2_km_eff=0.

         if (h2_km_eff.lt.0.) then
            write(*,*) 'h2_km',h2_km,'soil_C',soil_C,'h2_depth',h2_depth,'mod (Reji)', &
                 (h2_soilC_param(1)*soil_C/h2_depth*1./2650.*100.+h2_soilC_param(2))
            call error_mesg("land_tracer_driver","h2_km_eff<0",FATAL)
         end if

         R_bact = 1./max(sqrt(f_T*f_M*h2_km_eff*diff_H2),1.e-20)
      else
         !too dry for anything
         !I am setting diff_H2 to 0., f_T to 0., so that we can normalize by (1-h2_moist_r1)
         f_M         = 0.
         f_T         = 0.
         diff_H2     = 0.
         h2_km_eff   = 0.
         soil_C      = 0.
      end if

      !litter
      if (h2_litterC_mod.gt.0) then
         call tile%soilc%get_littC(litt_C)
         litterC_leaf = litt_C(LITT_LEAF)

         depth_litter_leaf = max(litterC_leaf/litter_leaf_density_C,0.) * h2_litterC_mod !m
!         depth_litter_wood = max(litterC_wood/litter_wood_density_C,0.) * h2_litterC_mod !m
            call soil_get_sfc_temp(tile%soil, grnd_T)
         if (depth_litter_leaf .gt. 0.) then
            R_litter_leaf = depth_litter_leaf/(diff_H2_air(grnd_T,p)*litter_leaf_porosity**2)
         end if
!         if (depth_litter_wood .gt. 0.) then
!            R_litter_wood = depth_litter_wood/(diff_H2_air(grnd_T,p)*litter_wood_porosity**2)
!         end if

!         R_litter = R_litter_wood + R_litter_leaf
         R_litter = R_litter_leaf
         depth_litter = depth_litter_leaf
      end if


      con = 1./(R_litter+R_inactive+R_snow+R_bact)

      if (id_con_h2_no_litter.gt.0) &
           call send_tile_data(id_con_h2_no_litter,  1./(R_inactive+R_snow+R_bact),      tile%diag)
      if (id_con_h2_no_snow.gt.0) &
           call send_tile_data(id_con_h2_no_snow,    1./(R_litter+R_inactive+R_bact),    tile%diag)

      call send_tile_data(id_h2_moist_r1, h2_moist_r1,        tile%diag)
      call send_tile_data(id_h2_moist_r2, h2_moist_r2,        tile%diag)
      call send_tile_data(id_h2_km,h2_km_eff,                 tile%diag)
      call send_tile_data(id_h2_fm,f_M,                       tile%diag)
      call send_tile_data(id_h2_ft,f_T,                       tile%diag)
      call send_tile_data(id_h2_sdiff,diff_h2,                tile%diag)
      call send_tile_data(id_h2_R_bact,R_bact,                tile%diag)
      call send_tile_data(id_h2_R_litter,R_litter,            tile%diag)
      call send_tile_data(id_h2_R_snow,R_snow,                tile%diag)
      call send_tile_data(id_h2_R_inactive,R_inactive,        tile%diag)
      call send_tile_data(id_h2_ilayer,inactive_layer,        tile%diag)
      call send_tile_data(id_h2_depth_litter,depth_litter,    tile%diag)
      call send_tile_data(id_h2_soilC,soil_C,                 tile%diag)


      call send_tile_data(id_h2_sws,    s_ws,  tile%diag)
      call send_tile_data(id_h2_sopt,   s_opt, tile%diag)
      call send_tile_data(id_h2_sup,    s_upc, tile%diag)

      call send_tile_data(id_h2_frac_water_pores_avg,frac_water_pores_avg, tile%diag)
      call send_tile_data(id_h2_frac_ice_pores_avg,frac_ice_pores_avg, tile%diag)

   end if

end function con_h2

elemental real function diff_H2_air(T,p) result(D)

   real, intent(in) :: T,p

   D = 0.668e-4*(101325./p)*(T/273.)**1.75

end function diff_H2_air

elemental real function diff_H2_snow(T,p,snow_porosity) result(D)
   real, intent(in) :: T,p,snow_porosity

   !https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009JD012459#jgrd15785-tbl-0001
   !D_eff = D * psi/tau^2
   !Assuming snow density of 200 kg/m3 and ice density of 917 kg m−3
   !psi = 1 - snow_density/ice_density = 1 - 200/917 = 0.78
   !tau^2 = 1.16 - 1.33 -> pick 1.2
   !This yields an effective snow porosity of 0.65

   !https://link.springer.com/article/10.1007/s10533-009-9302-3
   !proposes a different relationship
   !D_eff = D * psi * tau where tau is psi**(1/3)
   !D_eff = D * psi**4/3. This would yield an effective porosity of 0.72
   !When using eq. 4, this becomes 0.635

   !https://tc.copernicus.org/articles/16/967/2022/
   !porosity     = 1. - 300./917.
   !porosity_off = 0.078
   !porosity_res = max((porosity - porosity_off)/(1.-porosity_off),0.)
   !eff_snow_porosity   = porosity_res**1.61   => 0.645

   D = diff_H2_air(T,p) *  snow_porosity**(4./3.)

end function diff_H2_snow

elemental real function diff_H2_soil(T,p,n,st,b) result(D)
   !T: temperature (K)
   !p: pressure
   !n: porosity
   !st: water+ice fraction
   !b
   real, intent(in) :: T,p,n,st,b

   D = diff_H2_air(T,p) * n**2 * (max(1.-st,0.))**(2.+3./b)

 end function diff_H2_soil
end module land_tracer_driver_mod
