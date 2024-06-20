module snow_evolution_mod


#include <fms_platform.h>
#include "../../shared/debug.inc"

use mpp_mod, only: input_nml_file
use fms_mod, only : check_nml_error, stdlog, mpp_pe, mpp_root_pe, lowercase, &
       FATAL, WARNING, NOTE
use time_manager_mod, only: time_type_to_real
use constants_mod, only : GRAV, HLF, HLV, TFREEZE, PI
use land_constants_mod, only : NBANDS, BAND_NIR, BAND_VIS, &
! MODIS BRDF model parameters
    g_iso, g0_iso, g1_iso, g2_iso, &
    g_vol, g0_vol, g1_vol, g2_vol, &
    g_geo, g0_geo, g1_geo, g2_geo
use land_data_mod, only : lnd, log_version
use land_debug_mod, only : is_watch_point, land_error_message

use snicar_mod, only: compute_snicar_albedo
use snowpack_mod, only : snowpack_t, snow_layer_type, rho_water, rho_ice, LAI_ext, LAI_ssa, eps, &
    add_liquid_to_layer, compute_snow_grain_shape, merge_layers
use snow_tile_mod, only : NTRACERS, distinct_snow_on_glacier, cpw, clw, csw


implicit none
private

! public new_snow_density
! public snow_history_type
public :: snow_evolution_init
public :: read_F06_data
public :: gl_snow_step_2
public :: gl_sweep_tiny_snow
public :: gl_sweep_huge_snow
public :: gl_compute_snow_albedo
public :: read_snow_evolution_namelist
public :: use_internal_sources
! public :: min_snow_depth
public :: do_mgimplicit
public :: albedo_to_use
public :: thresh_snow_depth_swheat
public :: assign_substrate_sw_to_surface



! ! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snow_evolution_mod'
#include "../../shared/version_variable.inc"

!< scavenging coefficients for the tracters
!                                        (BC,  MD,  OM)
real, parameter :: SCAVENG(NTRACERS) = (/ 0.2, 0.0, 0.0 /)


! structure to store Flanner and Zender 2006 data
type :: data_F06_type
    REAL, ALLOCATABLE :: xT(:)
    REAL, ALLOCATABLE :: xDT(:)
    REAL, ALLOCATABLE :: xRHO(:)
    REAL, ALLOCATABLE :: drdt0s(:,:,:)
    REAL, ALLOCATABLE :: taus(:,:,:)
    REAL, ALLOCATABLE :: kappas(:,:,:)
    INTEGER nT
    INTEGER nDT
    INTEGER nRHO
end type

type(data_F06_type) :: dF06 ! create global module structure to store F06 data


! variable to store Snicar optical data
! type :: data_snicar_type
! end type

! type(data_snicar_type) :: dFSNI ! create global module structure to store F06 dat


! type :: site_data
! REAL lon ! site longitude
! REAL lat ! site latitude
! REAL elev ! site elevation [m msl]
! REAL z_wind ! height of wind speed measurements [m]
! REAL z0 ! momentum roughness heights [m]
! end type


! type(site_data) :: site ! information on current experimental site




! ! type variable to store simulation results for each model time step
! type :: snow_history_step_type
!     real doyt  !day of year + fractional time [days], use it with stored initial date
!     real snow_depth ! [m]
!     real runoff  ! [kg m^-2 s^-1]
!     real T_surface  ! [K]
!     real T_bottom  ! [K]
!     real levap   ! [kg m^-2 s^-1]
!     real fevap   ! [kg m^-2 s^-1]
!     real swe     ! [kg m^-2]
!     real liq     ! [kg m^-2]
!     real ice     ! [kg m^-2]
!     real heat    ! [J m^-2]
!     integer nlayers ! number of layers in the snowpack
!     real H ! sensible heat flux ! [W m^-2]
!     real LE ! latent heat flux [W m^-2]
!     real refl_dir_vis ! vis albedo, direct light[-]
!     real refl_dir_nir ! nir albedo, direct light[-]
!     real refl_dif_vis ! vis albedo, diffuse light[-]
!     real refl_dif_nir ! nir albedo, diffuse light[-]
!     real beta_vis ! vis light penetration length scale [m]
!     real beta_nir ! nir light penetration length scale [m]
!     real density ! [Kg m^-3]
!     real avrg_age ! [days]
!     real avrg_sph ! average sphericity [number in [0,1]]
!     real avrg_optd ! average optical diameter [m]
!     real avrg_T ! average temperature  [K]
!     real nrsf_age ! near surface layer snow age [days]
!     real nrsf_density ! near surface layer snow density [kg m^-3]
!     real nrsf_sph ! near surface layer snow sphericity [number in [0,1]]
!     real nrsf_optd ! near surface layer snow optical diameter [m]
!     real nrsf_bceq_im ! near surface layer equiv BC concentration, internally mixed [ppm]
!     real nrsf_bceq_em ! near surface layer equiv BC concentration, externally mixed [ppm]
! end type


! type :: snow_history_type
!     integer nsteps
!     type(snow_history_step_type), ALLOCATABLE :: hist(:) ! history
! contains
!     procedure :: print => snow_history_print
!     procedure :: write_timestep => snow_history_write_timestep
! end type



! ---- namelist
logical :: do_compaction = .true.
logical :: do_metamorph = .true.
logical :: do_wind_drift = .true.
logical :: do_mgimplicit = .true.
logical :: use_internal_sources = .false.
character(len=8) :: wlmax_to_use = 'CROCUS'  ! CROCUS, ANDERSON
character(len=12) :: albedo_to_use = 'CROCUS'  ! CROCUS, BRDF, HE  choices avail
character(len=12) :: albedo_correction_to_use = 'HE'  ! albedo correction due to impurities
character(len=12) :: metamor_model = 'f06'  ! available f06, c13
character(len=100) :: file_data_F06 = "INPUT/Snow_Flanner_drdt_bst_fit_60.nc"
logical :: do_split = .true. ! do relayering - split
logical :: do_merge = .true. ! do relayering - merge
logical :: do_snow_check_cons = .true.
real :: min_snow_mass = 0.1 ! to sweep tiny snow below this threshold [kg m^-3]
real :: min_snow_depth = 0.0 ! to sweep tiny snow below this threshold [m]
real :: max_snow = 1000.0 ! to sweep huge snow above this threshold [kg m^-3]
logical :: prevent_tiny_snow = .true.
logical :: correct_surface_T = .false.
real :: depth_surface_T_corr = 0.2
real :: thresh_snow_depth_swheat = 0.05 ! snow depth threshold [m] above which internal sw heat sources are computed
logical :: assign_substrate_sw_to_surface = .FALSE.
real :: min_fresh_density = 50.0 ! [kg/m3] minimum density for newly formed snow layers

namelist /snow_evolution_nml/ &
         do_compaction, do_metamorph, do_wind_drift, do_split, do_merge, &
         use_internal_sources, do_snow_check_cons, &
         min_snow_mass, min_snow_depth, max_snow, prevent_tiny_snow, do_mgimplicit, &
         metamor_model, file_data_F06, wlmax_to_use, albedo_to_use, &
         albedo_correction_to_use, correct_surface_T, depth_surface_T_corr, &
         thresh_snow_depth_swheat, assign_substrate_sw_to_surface, min_fresh_density
! ---- end of namelist

! ---- module data
real, public, protected :: delta_time ! model physics time step, s

contains  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



! ============================================================================
subroutine read_snow_evolution_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l            ! layer iterator

  call log_version(version, module_name, &
  __FILE__)
  read (input_nml_file, nml=snow_evolution_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_evolution_nml')
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=snow_evolution_nml)
  endif

  delta_time = time_type_to_real(lnd%dt_fast) ! [s]
end subroutine read_snow_evolution_namelist

! !> \Write current timestep variables to output
! subroutine snow_history_write_timestep(history, snowpack, step_index, &
!                                        runoff, fevap, levap, doyt, H, LE)


!   class(snow_history_type), intent(inout) :: history
!   class(snowpack_t), intent(in) :: snowpack
!   integer step_index
!   real runoff, fevap, levap, H, LE, doyt
!   real T_surface, T_bottom

!   if (snowpack%nlayers>0) then
!     T_surface = snowpack%snow(1)%T
!     T_bottom = snowpack%snow(snowpack%nlayers)%T
!   else
!     T_surface = -9999.9
!     T_bottom =  -9999.9
!   endif

!   ! save result to history file
!   history%hist(step_index)%nlayers = snowpack%nlayers
!   history%hist(step_index)%doyt = doyt
!   history%hist(step_index)%H = H
!   history%hist(step_index)%LE = LE
!   history%hist(step_index)%refl_dir_vis = snowpack%snow_refl_dir(1)
!   history%hist(step_index)%refl_dir_nir = snowpack%snow_refl_dir(2)
!   history%hist(step_index)%refl_dif_vis = snowpack%snow_refl_dif(1)
!   history%hist(step_index)%refl_dif_nir = snowpack%snow_refl_dif(2)
!   history%hist(step_index)%beta_vis = snowpack%beta_rad(1)
!   history%hist(step_index)%beta_nir = snowpack%beta_rad(2)
!   history%hist(step_index)%levap = levap
!   history%hist(step_index)%fevap = fevap
!   history%hist(step_index)%runoff = runoff
!   history%hist(step_index)%T_surface = T_surface
!   history%hist(step_index)%T_bottom = T_bottom
!   history%hist(step_index)%snow_depth = snowpack%depth()
!   history%hist(step_index)%swe = snowpack%swe()
!   history%hist(step_index)%liq = snowpack%liq()
!   history%hist(step_index)%ice = snowpack%ice()
!   history%hist(step_index)%heat = snowpack%heat()
!   history%hist(step_index)%density = snowpack%density()
!   history%hist(step_index)%avrg_age = snowpack%avrg_age()
!   history%hist(step_index)%avrg_sph = snowpack%avrg_sph()
!   history%hist(step_index)%avrg_optd = snowpack%avrg_optd()
!   history%hist(step_index)%avrg_T = snowpack%avrg_T()
!   history%hist(step_index)%nrsf_age = snowpack%nearsurf_age
!   history%hist(step_index)%nrsf_density = snowpack%nearsurf_rho
!   history%hist(step_index)%nrsf_sph = snowpack%nearsurf_sph
!   history%hist(step_index)%nrsf_optd = snowpack%nearsurf_optd
!   history%hist(step_index)%nrsf_bceq_im = snowpack%nearsurf_bceq_im
!   history%hist(step_index)%nrsf_bceq_em = snowpack%nearsurf_bceq_em

! end subroutine snow_history_write_timestep


! !> \Print state of snow history
! subroutine snow_history_print(h, filename)
!   class(snow_history_type), intent(in) :: h
!   ! real :: z
!   integer :: k
!   character(len=50) :: filename
!   ! write(*,*) "the number of history time steps is", h%nsteps
!   ! open(19, file=filename, status='REPLACE', access='SEQUENTIAL')
!   open(19, file=filename, status='REPLACE')
!   ! write(*,'(a2,99(",",a9,:))') "k","top","dz","T","ws","wl"
!   write(19,'(a6, a7, 99(",",a15,:))') "step", "nlayers", "doyt", "depth", "runoff", &
!                 "T_surface", "T_bottom", "levap", "fevap", "swe", "liq", "ice", "heat", &
!                 "sens_hf", "latent_hf", &
!                 "refl_dir_vis", "refl_dir_nir", &
!                 "refl_dif_vis", "refl_dif_nir",  &
!                 "beta_vis", "beta_nir",  &
!                 "density", "avrg_age", "avrg_sph", "avrg_optd", "avrg_T", &
!                 "nrsf_age", "nrsf_density", "nrsf_sph", "nrsf_optd", &
!                 "nrsf_bceq_im", "nrsf_bceq_em"
!   do k = 1, h%nsteps
!     !  write(*,'(i2.2,99(",",f9.4,:))') k, h%hist(k)%snow_depth, h%hist(k)%runoff, &
!             ! h%hist(k)%levap, h%hist(k)%fevap, &
!             ! h%hist(k)%swe, h%hist(k)%liq, h%hist(k)%ice, h%hist(k)%heat
!     write(19,'(i6, i7, 99(",",f23.8,:))') k, h%hist(k)%nlayers, h%hist(k)%doyt, h%hist(k)%snow_depth, &
!             h%hist(k)%runoff, h%hist(k)%T_surface, h%hist(k)%T_bottom, &
!             h%hist(k)%levap*1000000, h%hist(k)%fevap*1000000, &
!             h%hist(k)%swe, h%hist(k)%liq, h%hist(k)%ice, h%hist(k)%heat, &
!             h%hist(k)%H, h%hist(k)%LE, &
!             h%hist(k)%refl_dir_vis, h%hist(k)%refl_dir_nir, &
!             h%hist(k)%refl_dif_vis, h%hist(k)%refl_dif_nir, &
!             h%hist(k)%beta_vis, h%hist(k)%beta_nir, &
!             h%hist(k)%density, h%hist(k)%avrg_age, h%hist(k)%avrg_sph, &
!             h%hist(k)%avrg_optd, h%hist(k)%avrg_T, &
!             h%hist(k)%nrsf_age, h%hist(k)%nrsf_density, &
!             h%hist(k)%nrsf_sph, h%hist(k)%nrsf_optd, &
!             h%hist(k)%nrsf_bceq_im, h%hist(k)%nrsf_bceq_em
!   enddo
!   close(19)
! end subroutine snow_history_print


!> \initialize snowpack module, in particular read namelist parameters
subroutine snow_evolution_init()
  integer :: io, k, n
  real    :: dz ! layer thickness, for initialization of optimal vertical discretization, m

!   open (701, file='nml/input_snow_evolution.nml')
  open (701, file='nml/input.nml')
  read (701, snow_evolution_nml, iostat=io)
  if (io /= 0) stop 'Error reading input namelist "snow_evolution_nml"'
  close (701)
  write(*,snow_evolution_nml)

  ! initialize forcing

! call read_snicar_optics_data()

  ! initialize snowpack

  ! initialize snow history

end subroutine snow_evolution_init


!> \Check netcdf data
SUBROUTINE check(istatus)
  use netcdf
  INTEGER, INTENT (IN) :: istatus
  IF (istatus /= nf90_noerr) THEN
   write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
  END IF
END SUBROUTINE check


!> \Read netcdf F06 data fields
SUBROUTINE readgrid_F06(ncid,xpos,ypos, zpos,vdata1,vdata2, vdata3, NX,NY,NZ)
  USE netcdf
  REAL, DIMENSION(NX), INTENT(OUT) :: xpos
  REAL, DIMENSION(NY), INTENT(OUT) :: ypos
  REAL, DIMENSION(NZ), INTENT(OUT) :: zpos
  ! EAL, DIMENSION(NX,NY,NZ), INTENT(OUT) :: vdata1
  REAL, DIMENSION(NZ,NY,NX), INTENT(OUT) :: vdata1, vdata2, vdata3
  INTEGER, INTENT(IN) :: NX, NY, NZ
  INTEGER, DIMENSION(3) :: dimids
  INTEGER :: xtype, ndims, varid
  INTEGER, INTENT(IN) :: ncid
  CHARACTER(LEN=50) :: xname, yname, vname
  CALL check(nf90_inquire_variable(ncid,1,vname,xtype,ndims,dimids))
  CALL check(nf90_inq_varid(ncid,vname,varid))
  CALL check(nf90_get_var(ncid,varid,xpos))
!   write(*,*) "xpos var 1 = ", varid, vname
  CALL check(nf90_inquire_variable(ncid,2,vname,xtype,ndims,dimids))
  CALL check(nf90_inq_varid(ncid,vname,varid))
  CALL check(nf90_get_var(ncid,varid,ypos))
!   write(*,*) "ypos var 2 = ", varid, vname
  CALL check(nf90_inquire_variable(ncid,3,vname,xtype,ndims,dimids))
  CALL check(nf90_inq_varid(ncid,vname,varid))
  CALL check(nf90_get_var(ncid,varid,zpos))
!   write(*,*) "ypos var 3 = ", varid, vname
  CALL check(nf90_inquire_variable(ncid,4,vname,xtype,ndims,dimids))
  CALL check(nf90_inq_varid(ncid,vname,varid))
  CALL check(nf90_get_var(ncid,varid,vdata1))
!   write(*,*) "dimids, vname var 4 (1st data)= ", dimids, vname
  CALL check(nf90_inquire_variable(ncid,5,vname,xtype,ndims,dimids))
  CALL check(nf90_inq_varid(ncid,vname,varid))
  CALL check(nf90_get_var(ncid,varid,vdata2))
!   write(*,*) "dimids, vname var 5 (2nd data) = ", dimids, vname
  CALL check(nf90_inquire_variable(ncid,6,vname,xtype,ndims,dimids))
  CALL check(nf90_inq_varid(ncid,vname,varid))
  CALL check(nf90_get_var(ncid,varid,vdata3))
!   write(*,*) "dimids, vname var 6 (3rd data) = ", dimids, vname
END SUBROUTINE readgrid_F06


!> \Read netcdf F06 data
SUBROUTINE read_F06_data()
  USE netcdf
  INTEGER :: ncid
  INTEGER, DIMENSION(3) :: dimids
  CHARACTER(LEN=100) :: infile, ndims, varid
  CHARACTER(LEN=100) :: xname, yname, zname, v1name, v2name, v3name
!   write(*,*) "-------------------------------------------------------------------------"
!   write(*,*) "Reading the F06 netcdf data file"
!   write(*,*) "file: ", file_data_F06
  CALL check(nf90_open(file_data_F06, nf90_nowrite, ncid))
  CALL check(nf90_inquire_dimension(ncid,1,xname,dF06%nT))
  CALL check(nf90_inquire_dimension(ncid,2,yname,dF06%nDT))
  CALL check(nf90_inquire_dimension(ncid,3,zname,dF06%nRHO))
  ! allocate arrays
  ALLOCATE(dF06%drdt0s(dF06%nRHO, dF06%nDT, dF06%nT))
  ALLOCATE(dF06%taus  (dF06%nRHO, dF06%nDT, dF06%nT))
  ALLOCATE(dF06%kappas(dF06%nRHO, dF06%nDT, dF06%nT))
  ALLOCATE(dF06%xT(df06%nT))
  ALLOCATE(dF06%xDT(df06%nDT))
  ALLOCATE(dF06%xRHO(df06%nRHO))
  call readgrid_F06(ncid, dF06%xT, dF06%xDT, dF06%xRHO, dF06%taus, dF06%kappas,  &
                  dF06%drdt0s, dF06%nT, dF06%nDT, dF06%nRHO)
!   write(*,*) "nT = ", dF06%nT
!   write(*,*) "nDT = ", dF06%nDT
!   write(*,*) "nRHO = ", dF06%nRHO
!   write(*,*) "xT = ", dF06%xT
!   write(*,*) "xDT = ", dF06%xDT
!   write(*,*) "xRHO = ", dF06%xRHO
  CALL check(nf90_close(ncid))
!   write(*,*) "-------------------------------------------------------------------------"
END SUBROUTINE read_F06_data


!> \If snow amount is below specified limit, sweeps it into runoff [adapt. from lm4p2]
subroutine gl_sweep_tiny_snow(snowpack, lrunf, frunf, hlrunf, hfrunf,lost_wc_em, lost_wc_im )
  type(snowpack_t), intent(inout) :: snowpack
  real, intent(out) :: lrunf, frunf, hlrunf, hfrunf
  real, intent(out), dimension(NTRACERS) :: lost_wc_em, lost_wc_im
  real :: snow_mass
  integer :: il, it
  real snow_depth

   lost_wc_em = (/0.0, 0.0, 0.0 /)
   lost_wc_im = (/0.0, 0.0, 0.0 /)
  lrunf=0 ; frunf=0 ; hlrunf=0 ; hfrunf=0
  if (.not.prevent_tiny_snow) return ! do nothing, return zeros
  if (snowpack%nlayers==0) return ! no snowpack,return zeros
  ! does not change top water and top snow deficit if any
!    lost_wc_em = (/0.0, 0.0, 0.0 /)
!    lost_wc_im = (/0.0, 0.0, 0.0 /)

  snow_mass  = snowpack%ice() - snowpack%topsnowdeficit ! do not sweep aweay the deficit if any
  snow_depth = snowpack%depth()
!   min_snow_depth = 1E-4 ! 1/2 mm
!   min_snow_depth = 1E-3 ! 1 mm
!   min_snow_depth = 1E-2 ! 1 cm
  ! check if the snow is small enough to warrant sweeping
  if ( snow_mass<0 .or. snow_mass >= min_snow_mass ) return
!   if ( snow_mass<0 .or. ((snow_mass >= min_snow_mass).and.(snow_depth >= min_snow_depth) )) return

  lrunf  = snowpack%liq() ! do also sweep away topwater if any
  frunf  = snow_mass
  hfrunf = 0.0 ! do not sweep aweay the top snow deficit if any
!   hlrunf = snowpack%topwheat + lrunf*HLF
!   hlrunf = snowpack%topwheat - lrunf*HLF ! WAS
!   hlrunf = snowpack%topwheat - snowpack%topwater*HLF ! IS
  hlrunf = snowpack%topwheat ! IS
  do il = 1, snowpack%nlayers
     hlrunf = hlrunf + clw*snowpack%snow(il)%wl*(snowpack%snow(il)%T-tfreeze) + snowpack%snow(il)%wl * HLF
    !  hfrunf = hfrunf + csw*snowpack%snow(il)%ws*(snowpack%snow(il)%T-tfreeze)
    !  hlrunf = hlrunf + clw*snowpack%snow(il)%wl*(snowpack%snow(il)%T-tfreeze)
    !  hfrunf = hfrunf + csw*snowpack%snow(il)%ws*(snowpack%snow(il)%T-tfreeze) - snowpack%snow(il)%ws * HLF ! WAS
     hfrunf = hfrunf + csw*snowpack%snow(il)%ws*(snowpack%snow(il)%T-tfreeze) ! IS
     do it = 1, NTRACERS
       lost_wc_em(it) =lost_wc_em(it)+ snowpack%snow(il)%wc_em(it)
       lost_wc_im(it) =lost_wc_im(it)+ snowpack%snow(il)%wc_im(it)
     enddo
  enddo
!   write(*,*) "sweeping tiny snow ..."
   deallocate(snowpack%snow)
   snowpack%nlayers = 0.0
   snowpack%topwater = 0.0 ! sweep away liquid only in this case
   snowpack%topwheat = 0.0 ! sweep away liquid only in this case
end subroutine gl_sweep_tiny_snow


!> \If snow amount is above specified limit, sweeps it into runoff [adapt. from lm4p2]
subroutine gl_sweep_huge_snow(snowpack, lrunf, frunf, hlrunf, hfrunf, lost_wc_em, lost_wc_im )
  type(snowpack_t), intent(inout) :: snowpack
  real, intent(out) :: lrunf, frunf, hlrunf, hfrunf
  real, intent(out), dimension(NTRACERS) :: lost_wc_em, lost_wc_im
  real :: snow_mass
  integer :: il, it
  real snow_depth
  real cumsum_snow
  integer il_max_to_keep

    lost_wc_em = (/0.0, 0.0, 0.0 /)
    lost_wc_im = (/0.0, 0.0, 0.0 /)
    lrunf=0 ; frunf=0 ; hlrunf=0 ; hfrunf=0
    snow_mass = snowpack%ice() + snowpack%liq()
    ! write(*,*) "Checking whether Sweeping snow above max_snow = ", max_snow
    if (( snow_mass .le. max_snow ).or.(snowpack%nlayers==0)) return
    ! ELSE : REMOVE ALL SNOW IN EXCESS OF max_snow
    ! write(*,*) "Sweeping snow above max_snow = ", max_snow

    ! first, compute cumulative snow mass for each layer starting from top
    ! Then, remove bottom layers
    cumsum_snow = 0.0
    il_max_to_keep = -1
    ! do il=1,snowpack%nlayers
    !     cumsum_snow = cumsum_snow + snowpack%snow(il)%ws + snowpack%snow(il)%wl
    !     if (cumsum_snow>max_snow) then
    !         il_max_to_keep = il
    !     endif
    ! enddo
    il = 1
    do while ((il .le. snowpack%nlayers).and.(il_max_to_keep<0))
        cumsum_snow = cumsum_snow + snowpack%snow(il)%ws + snowpack%snow(il)%wl
        if (cumsum_snow > max_snow) then
            il_max_to_keep = il
        endif
        il = il + 1
    enddo


    ! NOW SINCE WE ARE IN THE CASE OF SNOW MASS > MAX_SNOW THIS SHOULD HOLD:
    ! (UNLESS TOPWATER IS UNPHYSICALLY LARGE)
    if (il_max_to_keep < 0) then
        write(*,*) "snowpack%nlayers = ", snowpack%nlayers
        write(*,*) "snowpack%ice() = ", snowpack%ice()
        call land_error_message("gl_sweep_huge_snow in snow_evolution_mod :: il_max_to_keep should be positive here!", FATAL)
    endif

    ! FIRST GET RID OF ALL ADDITIONAL LAYERS IF ANY
    if (il_max_to_keep < snowpack%nlayers) then
        do il=il_max_to_keep+1,snowpack%nlayers
            lrunf = lrunf + snowpack%snow(il)%wl
            frunf = frunf + snowpack%snow(il)%ws
            hlrunf = hlrunf + clw*snowpack%snow(il)%wl*(snowpack%snow(il)%T-tfreeze) + snowpack%snow(il)%wl * HLF
            hfrunf = hfrunf + csw*snowpack%snow(il)%ws*(snowpack%snow(il)%T-tfreeze)
            do it = 1, NTRACERS
              lost_wc_em(it) =lost_wc_em(it)+ snowpack%snow(il)%wc_em(it)
              lost_wc_im(it) =lost_wc_im(it)+ snowpack%snow(il)%wc_im(it)
            enddo
        enddo
        snowpack%nlayers = il_max_to_keep
        snowpack%snow = snowpack%snow(1:snowpack%nlayers) ! not really needed
    endif
    ! DO NOT MODIFY TOPWATER AND TOP SNOW DEFICIT
end subroutine gl_sweep_huge_snow


!> \Density of freshly fallen snow, from CROCUS (Vionnet et al., 2012)
subroutine new_snow_density(rho_new, Tatm, Ubar)
   real, INTENT(OUT) :: rho_new ! density of fresh fallen snow [kg m^-3]
   real, INTENT(IN) :: Tatm ! atm temp [K]
   real, INTENT(IN) :: Ubar ! wind speed [m/s]
   real, PARAMETER :: ar = 109.0 ! [Kg m^-3]
   real, PARAMETER :: br = 6.0 ! [Kg m^-3 K-1]
   real, PARAMETER :: cr = 26.0 ! Kg m^-7/2 s^-1/2
!    real, PARAMETER :: rho_min = 50.0 ! Kg m-3 ! // FIXME was 50 in CROCUS
   rho_new = ar + br * (Tatm - TFREEZE) + cr * sqrt(Ubar)
   rho_new = max(rho_new, min_fresh_density)
end subroutine new_snow_density


!> \compute the maximum liquid water content for a snow layer
real function compute_wlmax(depth_snow_layer, rho_snow_layer) result(wlmax)
real, intent(in) :: depth_snow_layer ! depth of the snow layer [m]
real, intent(in) :: rho_snow_layer ! density of the snow layer [kg m^-3]
real theta_crocus, theta_englesson, theta
real Crmax, Crmin, gamma_e
! Fraction of pore spaces occupied by water
! Or, from Englesson 1971 / from Dingman's book:
! not used for now
! theta_englesson = -0.0735*(rho_snow_layer/rho_water)+2.67*10.0**(-4)*(rho_snow_layer**2/rho_water)
! write(*,*) "theta englesson", theta_englesson
    if (trim(lowercase(wlmax_to_use)) == "crocus") then
        ! MODEL BY VIONNET ET AL., 2012
        theta_crocus = 0.05
        theta = theta_crocus
        wlmax = theta * rho_water * depth_snow_layer * (1.0 - rho_snow_layer/rho_ice)
        ! wlmax = theta * rho_water * depth_snow_layer
        ! write(*,*) "wlmax2", wlmax
    else if (trim(lowercase(wlmax_to_use)) == "anderson") then
        ! MODEL BY ANDERSON 1976 [Used e.g., by Shresta et al., 2006 and Arduini et al., 2019]
        Crmin = 0.03
        Crmax = 0.1
        gamma_e = 200.0 ! [Kg/m^3]
        if (rho_snow_layer - gamma_e > 0.0) then
            wlmax = Crmin
        else
            wlmax = Crmin + (Crmax-Crmin) * (1.0 - rho_snow_layer/gamma_e)
        endif
        wlmax = wlmax * rho_water * depth_snow_layer
        ! print("psurf, Tsurf, wlmax = ")

    else
        call land_error_message("Error in compute_wlmax in snow_evolution_mod: Must specify a valid snow liquid holding capacity model!", FATAL)
        wlmax = 0.0
    endif

end function compute_wlmax


!> \compute dendricity based on sphericity and dopt (Eq. (2) in Carmagnola et al., 2013)
real function den_from_dopt(sph, dopt) result(den)
real, intent(in) :: sph ! snow sphericity
real, intent(in) :: dopt ! snow optical diameter
real alpha
alpha = 1E-4
den = (dopt/alpha - 4.0 + sph)/(sph - 3.0)
end function den_from_dopt



!> \Compute snowpack vertical temperature grandient
subroutine vert_temperature_gradient(s, ilayer, G)
   ! Compute the (abs. value of the) vertical temperature gradient
   ! around layer i of snowpack s
   ! If i=1, estimate with forward differences
   ! If i=nlayers, estimate with backward differences
   ! Else, use differences centered in the interval
   ! approximate, due to the irregular width of layers
   class(snowpack_t), intent(inout) :: s !< state of snowpack
   integer, intent(in) :: ilayer ! index for snow layers
   real, intent(out) :: G ! temperature gradient [K m^-1] for selected snow layer i
   real dz1, dz2, dz, dT, dz3

    if (s%nlayers > 1) then
        if (ilayer==1) then
            dz1 = s%snow(ilayer)%dz
            dz2 = s%snow(ilayer+1)%dz
            dz = dz1/2.0 + dz2/2.0
            dT = s%snow(ilayer+1)%T - s%snow(ilayer)%T
        else if (ilayer==s%nlayers) then
            dz1 = s%snow(ilayer-1)%dz
            dz2 = s%snow(ilayer)%dz
            dz = dz1/2.0 + dz2/2.0
            dT = s%snow(ilayer)%T - s%snow(ilayer-1)%T
        else
            dz1 = s%snow(ilayer-1)%dz
            dz2 = s%snow(ilayer)%dz
            dz3 = s%snow(ilayer+1)%dz
            dz = dz1/2.0 + dz2 + dz3/2.0
            dT = s%snow(ilayer+1)%T - s%snow(ilayer-1)%T
        endif
    else if (s%nlayers == 1) then
        ! Only one snow layer, gradient undefined: just pass a small value
        dT = 1.0
        dz = 1.0
    else ! no snow
        ! error stop "attempting to compute temperature gradients but there is no snow"
        dT = 1.0
        dz = 1.0
    endif
    ! write(*,*) "Temp grad: dT, dz = ", dT, dz
   G = abs(dT/dz)
end subroutine vert_temperature_gradient


!> \Compute snow metamorphism
subroutine snow_metamorph(snowpack, dt, verbose)

   class(snowpack_t), intent(inout) :: snowpack !< state of snowpack
   real, intent(in) :: dt ! model time step [s]
   logical, intent(in) :: verbose ! print to screen some additional info
   integer il ! index for snow layers
   real Gi, Ti ! vertical temperature grandient (abs. val.) around layer i [K m^-1]
   real theta_i ! liquid fraction = fraction of snow layer mass which is liquid [percent]
   real rho_i ! density of current snow layer [kg m^-3] liquid + soild
   real dsph, ddopt ! current increments in sphericity [dimless] and optical diameter [m]
   real ddendr ! current increment in dendricity [dimless]
   real dt_days ! model time step in days
   real dt_hours ! model time step in hours
   logical is_wet ! switch for selecting wet or dry metamorphism law
   real max_optd, min_optd

   dt_days = dt/86400.0

   if (snowpack%nlayers > 0) then
    do il=1,snowpack%nlayers ! starts from top
        call vert_temperature_gradient(snowpack, il, Gi)
        Ti = snowpack%snow(il)%T
        ! liquid water content of snow layer: [theta = 100 wl/(wl+ws) = 100 wl/rho/depth ]
        theta_i = 100.0 * snowpack%snow(il)%wl / (snowpack%snow(il)%wl + snowpack%snow(il)%ws) ! percent !
        ! snow layer density [kg m^-3] liq + solid
        rho_i = (snowpack%snow(il)%ws + snowpack%snow(il)%wl)/snowpack%snow(il)%dz
        if(verbose) write(*,*) "snow metamorph: snow layer density rho_i and water content theta_i = ", rho_i, theta_i

        ! In any case, use C13 for wet metamorphism, F06 for dry only if selected
        ! if ((metamor_model=='C13').or.((metamor_model == 'F06').and.(theta_i > eps))) then
        ! Updated: F06 routine includes wet snow aging from Brun now
        ! Note in the C13 case is selected, dendriticy is not a prognostic model variable
        ! It is recomputed each time based on (sph, dendr values)
        if (trim(lowercase(metamor_model)) == 'c13') then ! use Carmagnola 2013 formulation

            if (theta_i > eps) then ! wet snow metamorphism
                is_wet = .TRUE.
                if(verbose) write(*,*) "Computing Carmagnola 2013 [c13] wet metamorphism ..."
                call wet_metamorph_carmagnola(dsph, ddopt, &
                            snowpack%snow(il)%sph, snowpack%snow(il)%optd,&
                            theta_i, Ti,  dt_days, verbose=verbose)
            ! if (ddopt < 0.0) then
            ! error stop "wet metamo - Found decrease in optical diameter!"
            ! endif
            else
                is_wet = .FALSE.
                if(verbose) write(*,*) "Computing Carmagnola 2013 [c13] dry metamorphism ..."
                call dry_metamorph_carmagnola(dsph, ddopt, &
                            snowpack%snow(il)%sph, snowpack%snow(il)%optd,&
                            Ti, Gi, rho_i,  dt_days, verbose=verbose)
            ! if (ddopt < 0.0) then
            ! error stop "dry metamo - Found decrease in optical diameter!"
            endif
        else if (trim(lowercase(metamor_model)) == 'f06') then
            if (verbose) write(*,*) "Computing dry snow metamorphism according to Flanner and Zender, 2006 [f06]"
            dt_hours = dt_days*24.0
            call metamorph_FlannerZender2006( &
                ddopt, snowpack%snow(il)%optd,&
                snowpack%snow(il)%ws, snowpack%snow(il)%wl,&
                Ti, Gi, rho_i,  dt_hours, verbose=verbose)
            ! F06 does not compute grain sphericity and dendriticy
            ! Use Brun 1992 Journal of Glaciology formulation:
            if (theta_i>eps) then
                call wet_metamorph_brun(dsph, ddendr, snowpack%snow(il)%sph,snowpack%snow(il)%dendr, theta_i, dt_days)
            else
                call dry_metamorph_brun(dsph, ddendr, snowpack%snow(il)%sph,snowpack%snow(il)%dendr, Ti, Gi, rho_i, dt_days)
            endif

        else
            ! error stop "ERROR snow_metamorph in snow_evolution module: specify a valid snow metamorphism model!"
            call land_error_message("ERROR snow_metamorph in snow_evolution module: specify a valid snow metamorphism model!", FATAL)
        endif

        ! update sphericity, optical diameter and dendriticy for current time step
        snowpack%snow(il)%sph = snowpack%snow(il)%sph + dsph
        snowpack%snow(il)%sph = max(0.0, min(1.0, snowpack%snow(il)%sph))
        snowpack%snow(il)%optd = snowpack%snow(il)%optd + ddopt
        snowpack%snow(il)%optd = max(1E-4, snowpack%snow(il)%optd)

        ! update dendriticy here because this is done only for F06 formaulation
        ! in C13 it is a derived quantity from optical diameter and sphericity
        ! else in F06 case evolve it dynamically using Brun 1992 laws
        if (trim(lowercase(metamor_model)) == 'c13') then
            snowpack%snow(il)%dendr = den_from_dopt(snowpack%snow(il)%sph, snowpack%snow(il)%optd)
        else if (trim(lowercase(metamor_model)) == 'f06') then
            snowpack%snow(il)%dendr = snowpack%snow(il)%dendr + ddendr
        else
            ! error stop "ERROR snow_metamorph in snow_evolution module: specify a valid snow metamorphism model!"
            call land_error_message( "ERROR snow_metamorph in snow_evolution module: specify a valid snow metamorphism model!", FATAL)
        endif
        snowpack%snow(il)%dendr = max(0.0, min(1.0, snowpack%snow(il)%dendr))


            ! max_optd = 1500*10**(-6) ! meters, diameter
            ! min_optd = 30*10**(-6) ! meters, diameter
            ! snowpack%snow(il)%optd = max(min( snowpack%snow(il)%optd , max_optd), min_optd)


        if ((snowpack%snow(il)%optd <0.0).or.(snowpack%snow(il)%optd > 1.0)) then !
            write(*,*) "ddopt = ", ddopt
            write(*,*) "s%dopt = ", snowpack%snow(il)%optd
            write(*,*) "dsph = ", dsph
            write(*,*) "s%sph = ", snowpack%snow(il)%sph
            call snowpack%print()
            if (is_wet) then
                ! error stop "ERROR snow_metamorph in snow_evolution module: optical diameter out of bounds after snow wet metamorph calculation"
                call land_error_message("ERROR snow_metamorph in snow_evolution module: optical diameter out of bounds after snow wet metamorph calculation", FATAL)
            else
                write(*,*) "Ti = ", Ti
                write(*,*) "Gi = ", Gi
                write(*,*) "rho_i = ", rho_i
                ! error stop "ERROR snow_metamorph in snow_evolution module: optical diameter out of bounds after snow dry metamorph calculation"
                call land_error_message("ERROR snow_metamorph in snow_evolution module: optical diameter out of bounds after snow dry metamorph calculation", FATAL)
            endif
        endif
        if ((snowpack%snow(il)%sph <0.0).or.(snowpack%snow(il)%sph > 1.0)) then ! sph bounds
            if (is_wet) then
                ! error stop "ERROR snow_metamorph in snow_evolution module: sphericity out of bounds after snow wet metamorph calculation"
                call land_error_message("ERROR snow_metamorph in snow_evolution module: sphericity out of bounds after snow wet metamorph calculation", FATAL)
            else
                ! error stop "ERROR snow_metamorph in snow_evolution module: sphericity out of bounds after snow dry metamorph calculation"
                call land_error_message("ERROR snow_metamorph in snow_evolution module: sphericity out of bounds after snow dry metamorph calculation", FATAL)
            endif
        endif
        enddo
    endif
end subroutine snow_metamorph


!> \compute snow metamorphism following Flanner and Zender, 2013, and Olson et al., 2010
! includes dry and wet metamorph component
subroutine metamorph_FlannerZender2006(ddopt, dopt, ws, wl, Ti, Gi, rho_i, &
                                       dt_hours, verbose)
    real, intent(OUT) ::ddopt ! predicted change in opt diameter [m]
    real, intent(IN) :: dopt ! current optical diameter [m]
    real, intent(IN) :: dt_hours ! time step [hours]
    logical, intent(IN) :: verbose
    real, INTENT(IN) :: ws, wl ! solid ice and liquid water [kg/m2]
    real, intent(IN) :: Ti, Gi, rho_i ! temp [K], temp grad [K/m] and density [kg m^-3] of current layer
    real term1, term2, term3
    real TiC
    real fT, hR, gG, Phi ! Marbouty's empirical functions
    real drdt, drdt0, re0
    real tau, k ! parameters from lookup table
    real old_re, new_re, fliq
    real kap, eta
    real f_new, f_refr, f_old ! fractions of new, refrozen and old snow used in CLM, not used here
    real re_0, re_refr
    real dre_dry, dre_wet
    real rho_i2, Ti2, Gi2
    integer irho, iGG, iTT
    integer irho2, iGG2, iTT2
    real re_abs_diff

    ! make sure vars are within bounds

    ! make sure vars are within bounds
    rho_i2 = min(max(50.0, rho_i), 400.0)
    Ti2 = max(223.0, Ti)
    Gi2 = Gi
    ! irho =MAX( MIN( ABS( INT( (rho_i2 - 25.0) / 50.0       ) + 1 ), 8  ), 1)
    ! iGG = MAX( MIN( ABS( INT( (Gi2 - 5.0   ) / 10.0 + 2.0  )     ), 31 ), 1)
    ! iTT = MAX( MIN( ABS( INT( (Ti2-225.65   ) / 5.0 + 2.0  )     ), 11 ), 1)
    irho =MAX( MIN( NINT( (rho_i2 - 50.0) / 50.0      )  + 1, 8  ), 1)
    iGG = MAX( MIN( NINT( (Gi2          ) / 10.0      )  + 1, 31 ), 1)
    iTT = MAX( MIN( NINT( (Ti2-223.0    ) / 5.0       )  + 1, 11 ), 1)

    ! get indices for lookup table
    ! irho2 = int( minloc( abs(dF06%xRHO - rho_i2), dim=1))
    ! iGG2 =  int( minloc( abs(dF06%xDT - Gi2), dim=1))
    ! iTT2 =  int( minloc( abs(dF06%xT - Ti2), dim=1))

    ! irho2 = min(max(int( minloc( abs(dF06%xRHO - rho_i2), dim=1)), 1), 8)
    ! iGG2 =  min(max(int( minloc( abs(dF06%xDT - Gi2), dim=1)), 1), 31)
    ! iTT2 =  min(max(int( minloc( abs(dF06%xT - Ti2), dim=1)), 1), 11)


    ! if (irho .ne. irho2) then
    !     write(*,*) "irho different index: irho, irho2 = ", irho, irho2
    !     error stop "irho different index."
    ! endif
    !     if (iGG .ne. iGG2) then
    !     write(*,*) "iGG different index: iGG, iGG2 = ", iGG, iGG2
    !     error stop "iGG different index."
    ! endif
    !     if (iTT .ne. iTT) then
    !     write(*,*) "iTT different index: iTT, iTT2 = ", iTT, iTT2
    !     error stop "iTT different index."
    ! endif

    ! write(*,*) "size taus = ", size(dF06%taus)
    tau = dF06%taus(irho, iGG, iTT)
    kap = dF06%kappas(irho, iGG, iTT)
    drdt0 = dF06%drdt0s(irho, iGG, iTT)

    ! write(*,*) "-------------------------------------"
    ! write(*,*) "Flanner F06 - check table values:"
    ! write(*,*) "rho, irho = ", rho_i2, irho
    ! write(*,*) "T, iT = ", Ti2, iTT
    ! write(*,*) "Gi, iGi = ", Gi2, iGG


    ! write(*,*) "Flanner tau = ", tau
    ! write(*,*) "Flanner kappa = ", kap
    ! write(*,*) "Flanner drdt0 = ", drdt0
    ! write(*,*) "-------------------------------------"

! // TODO: add -> track fresh and refrozen water as in CLM?
! instead, we treat wet metamorph separately as done in CROCUS (C13)
    f_new = 0.0
    f_old = 1.0
    f_refr = 0.0

    old_re = dopt/2.0 * 1E6 ! switch to a RADIUS in [\mu m] only for this parameterization
    ! paramaters fixed in the model:
    re_0 = 50.0 ! [\mu m] similar to 54.5 [\mu m] recommended, equivalent to my optical diameter
               ! of freshly fallen show which is set to 1E-4 [m] following Carmagnola, 2013
    re_refr = 1000.0 ! [\mu m] = 1 mm ! radius of refrozen snow  [not used here]

    ! As in CLM (Olson et al., 2010) use Wet snow aging from Brun (1989)
    ! after subtracting offset due to dry aging, here already accounted for
    ! here we sum the two contributions from wet and dry aging
    ! here expressed for radius change dre_wet in [\mu m]
    fliq = wl/(wl+ws)
    dre_wet = (dt_hours*3600.0) * (10.0**18 * 4.22 * 10.0**(-13) * (fliq)**3 )/(4.0*PI*old_re**2)

    re_abs_diff =old_re-re_0
    ! if (re_abs_diff < -1E-6) then
    !     write(*,*) "Flanner 2006: Re, Re0, Re - Re0 = ", old_re, re_0, re_abs_diff
    !     call land_error_message("Error in metamorph_FlannerZender2006 in snow_evolution_mod: re-re0 < 0 found!", FATAL)
    ! else
    !     re_abs_diff = max(1E-8, re_abs_diff)
    ! endif
    ! since the parameterization for snow drift can lead to optical diameters
    ! smaller than minimum value here, we set the difference to a positive value
    ! when this happens if the wind drift routine is used in the model.
    re_abs_diff = max(1E-8, re_abs_diff)

    dre_dry = drdt0 * (tau/(re_abs_diff + tau))**(1.0/kap) * dt_hours
    ! new_re = ( old_re + dre_dry + dre_wet )*f_old + re_0 * f_new + re_refr * f_refr

    ! return increment in optical diameter (=2 * radius... ), going back to [m]
    ddopt = ( dre_dry + dre_wet ) * 2.0 / 1E6


end subroutine metamorph_FlannerZender2006



!> \compute wet snow metamorphism according to Carmagnola, 2013
subroutine wet_metamorph_carmagnola(ds, ddopt, s, dopt, theta_i, Ti, dt_days, verbose)
    real, intent(OUT) :: ds, ddopt ! predicted changes in sphericity [adim.] and opt diameter [m]
    real, intent(IN)  :: s, dopt ! current sphericity [adim.] and optical diameter [m] of current layer.
    real, intent(IN) :: dt_days ! time step [days]
    real, intent(IN) :: theta_i ! water volumetric content [percent] of current layer
    real, intent(IN) :: Ti
    logical, intent(IN) :: verbose
    real dv ! Change in volume of grains - as in Brun 1989 [mm^3]
    real new_volume ! [mm**3]
    real ZVDENT1
    real dd, dgs
    real gs
    real dendr
    real fliq

    logical is_dendritic
    logical sphericity_is_1

    ! TODO: mult by 10**-9 to get res in [m^3/s]
    real, PARAMETER :: v0p = 1.28E-8 ! [mm^3 s^-1] From Brun et al., 1989
    real, PARAMETER :: v1p = 4.22E-10 ! [mm^3 s^-1] From Brun et al., 1989

    ! is_dendritic = (d > epsilon) ! Vionnet 2012
    sphericity_is_1 = (s > 1.0 - eps)

    is_dendritic = dopt < 1E-4 * (4.0 - s) ! Carmagnola 2013
    ! is_dendritic = dopt < 1E-4 * (4.0 - s) - eps ! Carmagnola 2013
    ! if(verbose) write(*,*) "snow id dendritic =", is_dendritic

    ! FROM SURFEX CODE :: evol of dendricity, humid case ::
    ! ZVDENT1 = MAX( XVDENT2 * ZTELM**NVDENT1, XVDENT1 * EXP(XVVAP1/XTT) )
    ! coeff is different because dt in seconds
    ! ZVDENT1 = MAX( 7.2338E-7 * theta_i**3, 2314.81481 * EXP(-6000.0 / Ti ) )
    !                  1/16                     2 * 10 ^ 8
    ! dd = MAX( -1.0/16.0*theta_i**3 * dt_days, -2.0*10**8 * EXP(-6000.0/Ti) * dt_days  )
    dd = -1.0/16.0*theta_i**3 * dt_days


    ! double check the boundaries of each case - which toll? does it matter?
    ! Following Vionnet et al., 2012 here.
    !  Updated with values from Brun et al., 1989
    if ( (.not. sphericity_is_1) .and. (.not. is_dendritic)) then
        if(verbose) write(*,*) "not spheric (sph < 1) and not dendritic snow"

        ds = -dd
        ! ddopt =  1E-4 * 2.0 * s * ds ! INCREASING, OK

        ddopt = 0.0 ! WRONG

        ! write(*,*) "WET C13 - NOT sphericity=1, NOT dendritic :: dopt, ddopt = ", dopt, ddopt

    else if ( (.not. sphericity_is_1) .and. is_dendritic) then
        if(verbose) write(*,*) "not spheric (sph < 1) and dendritic snow"

        ds = -dd
        ddopt = 1E-4  * (dd * (s-3.0) + ds * (dopt/1E-4 -1.0)/(s-3.0))

        ! write(*,*) "WET C13 - NOT sphericity=1, dendritic :: dopt, ddopt = ", dopt, ddopt

    else if ( (sphericity_is_1) .and. (.not. is_dendritic)) then
        if(verbose) write(*,*) "spheric (sph = 1) and not dendritic snow"
        ! change only in grain size in this case
        ! Use time in seconds only here for Brun parameterization of wet accretion
        ! theta_i in [percent]
        ! write(*,*) "WET C13 - Using Brun 1992 Snow grain growth"
        gs = dopt
    !  dv = (v0p + v1p*theta_i**3) * dt_days * 86400.0/10E9 !!!! From Brun 1989 - volume change in [mm^3]
    !   dv = (v0p + v1p*(theta_i/100.0)**3) * dt_days * 86400.0/1E9 !!!! From Brun 1989 - volume change in [mm^3]
    !   ! how can I convert this to a change in grain size?
    ! !   new_volume = gs**3 + dv*1E-9 ! EZDEV - volume back to m^3
    !   ! AS DONE IN SURFEX - SPHERICAL GRAIN WITH DIAMETER gs
    !   new_volume = 4.0 * PI / 3.0 * (gs/2.0)**3 + dv
    !   ! dgs = 2.0 * 3.0/(4.0*PI)*(new_volume)**(1.0/3.0) - gs ! increase in linear size
    !   dgs = 2.0 * (3.0/(4.0*PI)*new_volume)**(1.0/3.0) - gs ! increase in linear size
    !   ds = 0.0
    ! !   ddopt = 0.0
    !   ddopt = dgs ! correspond in this case
        ! write(*,*) "ddopt (3) = ", ddopt

        fliq = theta_i/100.0
        ! fliq = theta_i
        ! dre_wet = (dt_hours*3600.0) * (10.0**18 * 4.22 * 10.0**(-13) * (fliq)**3 )/(4.0*PI*old_re**2)
        ! \mu m
        ! ddopt = (dt_days*86400.0) * (10.0**18 * 4.22 * 10.0**(-13) * (fliq)**3 )/(4.0*PI*(dopt/2.0*1E6)**2)
        ! SAME AS IN CLM, BUT HERE ADD v0p COEFF -> TO INCLUDE EFFECT OF DRY METAMORPH,
        ddopt = (dt_days*86400.0) * (1.28E-8*1E9 + 10.0**18 * 4.22 * 10.0**(-13) * (fliq)**3 )/(4.0*PI*(dopt/2.0*1E6)**2)
        ddopt = ddopt * 2.0 / 1E6

        ! write(*,*) "WET C13 - Using Brun 1992 (SPHERICITY=1, NOT DENDRITIC) : dopt, ddopt = ", dopt, ddopt

    else if ( (sphericity_is_1) .and. is_dendritic) then
        if(verbose) write(*,*) "spheric (sph = 1) and dendritic snow"
        ! ds = 1.0/16.0*theta_i**3 * dt_days ! OK
        ! ddopt =  2.0 * 1E-4 * s * dd /dt_days
        ! ds = - dd
        ! ddopt =  - 2.0 * 1E-4 * s * dd /dt_days


        ds = -dd
        ddopt = 1E-4  * (dd * (s-3.0) + ds * (dopt/1E-4 -1.0)/(s-3.0))

        ! write(*,*) "WET C13 - sphericity=1, dendritic :: dopt, ddopt = ", dopt, ddopt

    else
        ! error stop "ERROR: wet_metamorph_carmagnola in snow_evolution module: this should not happen!"
        call land_error_message("ERROR: wet_metamorph_carmagnola in snow_evolution module: this should not happen!", FATAL)
    endif


    ! if (is_dendritic) then
    !   ds = 1.0/16.0*theta_i**3 * dt_days
    !   dd = - 1.0/16.0*theta_i**3 * dt_days
    !   dgs = 0.0
    ! else ! reached sphericity
    !     ds = 0.0
    !     dd = 0.0
    !     dgs = 0.0
    ! endif


end subroutine wet_metamorph_carmagnola


!> \compute dry snow metamorphism according to Carmagnola, 2013
subroutine dry_metamorph_carmagnola(ds, ddopt, s, dopt, Ti, Gi, rho_i, dt_days, verbose)
    real, intent(OUT)::  ds, ddopt ! predicted changes in sphericity [adim] and opt diameter [m]
    real, intent(IN) :: s, dopt ! current sphericity of current layer [adim.] and opt diameter [m]
    real, intent(IN) :: dt_days ! time step [days]
    real, intent(IN) :: Ti, Gi, rho_i ! temp [K], temp grad [K/m] and density for current layer
    logical, intent(IN) :: verbose
    real term1, term2, term3 !
    real TiC
    real fT, hR, gG, Phi ! empirical functions
    logical is_dendritic


    TiC = Ti-TFREEZE
    ! empirical functions from Marbouty (1980)
    if (TiC <= - 40.0) then
        fT = 0.0
    else if ((TiC > -40.0).and.(TiC <= -22.0)) then
        fT = 0.011 * ( TiC + 40.0 )
    else if ((Tic > -22.0).and.(TiC <= -6)) then
        fT = 0.2 + 0.05 * (TiC + 22.0)
    else
        fT = 1.0 - 0.05 * (TiC)
    endif


    if (rho_i < 150.0) then
    else if ((rho_i >= 150.0).and.(rho_i < 400.0)) then ! kg m^-3
        hR = 1.0 - 0.004 * (rho_i - 150.0)
    else
        hR = 0.0
    endif

    if (Gi < 15.0) then
        gG = 0.0
    else if ((Gi >= 15.0).and.(Gi<25.0)) then
        gG = 0.01 * (Gi - 15.0)
    else if ((Gi >= 25.0).and.(Gi<40.0)) then
        gG = 0.1 + 0.037 * (Gi - 25.0)
    else if ((Gi >= 40.0).and.(Gi<50.0)) then
        gG = 0.65 + 0.02 * (Gi - 40.0)
    else if ((Gi >= 50.0).and.(Gi<70.0)) then
        gG = 0.85 + 0.0075 * (Gi - 50.0)
    else
        gG = 1.0
    endif

    Phi = 1.0417*10.0**(-9) ! [m s^-1]


    term1 = 1E9 * exp( - 6000.0 / Ti )
    term2 =  (-2.0 * 10**8) * exp( -6000.0 / Ti)
    term3 = (dopt/1E-4 - 1.0)/(s-3.0)



    ! is_dendritic = (d > 0.0001)
    is_dendritic = dopt < 1E-4 * (4.0 - s) ! Carmagnola 2013
    ! is_dendritic = dopt < 1E-4 * (4.0 - s) -eps ! Carmagnola 2013
    ! if(verbose) write(*,*) "snow id dendritic =", is_dendritic

    if (.not.is_dendritic) then
        if (Gi <= 5.0) then
            if(verbose) write(*,*) " not dendritic snow, mild gradient dT/dz"
            ds = dt_days * term1           ! OK, INCREASE
            ddopt = - 2.0 * 1E-4 * s * ds  ! OK, DECREASE
        else if ((Gi > 5.0) .and. (Gi <= 15.0)) then
            if(verbose) write(*,*) " not dendritic snow, intermediate gradient dT/dz"
            ds = dt_days * term2 * Gi**0.4 ! OK, DECREASE
            ddopt = - 2.0 * 1E-4 * s * ds  ! OK, INCREASE
        else
            if (s > 1E-5) then
                if(verbose) write(*,*) " not dendritic snow, steep gradient dT/dz, sphericity > 0"
                ds = dt_days * term2 * Gi**0.4 ! DECREASE
                ddopt = - 2.0 * 1E-4 * s * ds ! INCREASE
            else ! case s=0
                if(verbose) write(*,*) " not dendritic snow, steep gradient dT/dz, sphericity == 0"
                ds = 0.0     ! OK
                ddopt = 0.5 * dt_days * fT * hR * gG * Phi ! OK, INCREASE
            endif
        endif
    else ! case of dendritic snow
        if (Gi <= 5.0) then
        if(verbose) write(*,*) " dendritic snow, mild gradient dT/dz"
            ds = dt_days * term1 ! OK, INCREASE
            ddopt = dt_days * 1E-4 * (term2*(s-3.0) + ds/dt_days * term3 ) ! OK - First term -, Second + if dopt > alpha as it should be
            ! write(*,*)  'term A = ', dt_days * 1E-4 * term2*(s-3.0)
            ! write(*,*)  'term B = ', dt_days * 1E-4 * ds/dt_days * term3
            ! write(*,*)  'ds/dt = ', ds/dt_days
            ! write(*,*)  's = ', s
            ! write(*,*) "ddopt (5) = ", ddopt
            ! write(*,*) "term3 = ", dopt/1E-4
        else  ! medium or high gradient & dendritic snow [treated to gether in B92 and C13]
        if(verbose) write(*,*) " dendritic snow, intermediate gradient dT/dz"
            ds = dt_days * term2 * Gi**0.4 ! OK, DECREASE
            ddopt = dt_days * 1E-4 * (term2 * Gi**0.4*(s-3.0) + ds/dt_days * term3 ) ! OK - Both terms negative
            ! write(*,*)  'term A = ', dt_days * 1E-4 * term2*(s-3.0)
            ! write(*,*)  'term B = ', dt_days * 1E-4 * ds/dt_days * term3
            ! write(*,*) "ddopt (6) = ", ddopt
        endif
    endif
end subroutine dry_metamorph_carmagnola


!> \compute wet snow metamorphism according to CROCUS
! only for sphericity and dendricity - do opt diam separately
! Following Vionnet et al., 2012 and Brun et al., 1992
subroutine wet_metamorph_brun(ds, dd, s, d, theta_i, dt_days)
    real, intent(OUT)::  ds, dd ! predicted changes in sphericity and dendricity [adim.]
    real, intent(IN) :: s, d ! current sphericity and dendricity of current layer [adim.]
    real, intent(IN) :: dt_days ! time step [days]
    real, intent(IN) :: theta_i ! water volumetric content [percent] of current layer
    logical is_dendritic
    logical sphericity_is_1
    is_dendritic = (d > 1E-8)
    sphericity_is_1 = (s > 1E-8)
    if ( sphericity_is_1) then
        ds = 0.0
    else
        ds = 1.0/16.0*theta_i**3 * dt_days
    endif
    if ( is_dendritic) then
        dd = -1.0/16.0*theta_i**3 * dt_days
    else
        dd = 0.0
    endif
end subroutine wet_metamorph_brun


!> \compute dry snow metamorphism according to CROCUS
! only for sphericity and dendricity - do opt diam separately
! Following Vionnet et al., 2012 and Brun et al., 1992
subroutine dry_metamorph_brun(ds, dd, s, d, Ti, Gi, rho_i, dt_days)
    real, intent(OUT) :: ds, dd ! predicted changes in sphericity and dendricity [adim.]
    real, intent(IN)  :: s, d ! current sphericity and dendricity of current layer [adim.]
    real, intent(IN)  :: dt_days ! time step [days]
    real, intent(IN)  :: Ti, Gi, rho_i ! temp [K], temp grad [K/m] and density for current layer
    real term1, term2 !
    real TiC
    real fT, hR, gG, Phi ! empirical functions
    logical is_dendritic

    TiC = Ti-TFREEZE
    ! empirical functions from Marbouty (1980)
    if (TiC <= - 40.0) then
        fT = 0.0
    else if ((TiC > -40.0).and.(TiC <= -22.0)) then
        fT = 0.011 * ( TiC + 40.0 )
    else if ((Tic > -22.0).and.(TiC <= -6)) then
        fT = 0.2 + 0.05 * (TiC + 22.0)
    else
        fT = 1.0 - 0.05 * (TiC)
    endif

    if (rho_i < 150.0) then
    else if ((rho_i >= 150.0).and.(rho_i < 400.0)) then ! kg m^-3
        hR = 1.0 - 0.004 * (rho_i - 150.0)
    else
        hR = 0.0
    endif

    if (Gi < 15.0) then
        gG = 0.0
    else if ((Gi >= 15.0).and.(Gi<25.0)) then
        gG = 0.01 * (Gi - 15.0)
    else if ((Gi >= 25.0).and.(Gi<40.0)) then
        gG = 0.1 + 0.037 * (Gi - 25.0)
    else if ((Gi >= 40.0).and.(Gi<50.0)) then
        gG = 0.65 + 0.02 * (Gi - 40.0)
    else if ((Gi >= 50.0).and.(Gi<70.0)) then
        gG = 0.85 + 0.0075 * (Gi - 50.0)
    else
        gG = 1.0
    endif

    Phi = 1.0417*10.0**(-9) ! [m s^-1]


    term1 = 1E9 * exp( - 6000.0 / Ti )
    term2 =  -2.0 * 10**8 * exp( -6000.0 / Ti)
    is_dendritic = (d > 1E-8)

    if (.not.is_dendritic) then
        dd = 0.0
        if (Gi <= 5.0) then
            ds = dt_days * term1
        else if (Gi <= 15) then
            ds = dt_days * term2 * Gi**0.4
        else
            if (s > 0.0) then
                ds = dt_days * term2 * Gi**0.4
            else
                ds = 0.0
            endif
        endif
    else ! case of dendritic snow
        if (Gi <= 5.0) then ! mild temperature gradient
            dd = dt_days * term2
            ds = dt_days * term1
        else ! medium/high temperature gradient
            dd = dt_days * term2 * Gi**0.4
            ds = dt_days * term2 * Gi**0.4
        endif
    endif
end subroutine dry_metamorph_brun


!> \Effect of wind drift on snow, following Brun et al., 1997 and Vionnet et al., 2012
subroutine snow_wind_drift_C13(snowpack, dt, Ubar, verbose)

    class(snowpack_t), intent(inout) :: snowpack !< state of snowpack
    real, INTENT(IN) :: dt ! time step [s]
    real, INTENT(IN) :: Ubar ! wind speed ! [m s^-1]
    logical, intent(IN) :: verbose
    integer il
    real, PARAMETER :: rho_min = 50.0 ! kg m^-3
    real, PARAMETER :: rho_max = 350.0 ! kg m^-3
    real, PARAMETER :: tau_48h = 48.0 ! 48 hours in [hours]
    real M0i ! mobility index for layer il = potntial for snow erosion
    real SLi ! driftability index
    real rho_i
    real new_rho_i
    real tau_i
    real Gamma_i_drift
    real pseudo_zi
    real dt_hours
    logical is_dendritic
    real term_a, term_b
    real ds, dd, dgs, drho ! variations in snow properties due to wind drift
    real den, sph
    real Frho
    real gs, dendr ! computed from s and optical diameter
    real dopt, ddopt
    real ddendr

    dt_hours = dt/3600.0

    ! depth_below_surface = 0.0
    pseudo_zi = 0.0
    if(verbose) write(*,*) "Running wind drift routine!"

    if (snowpack%nlayers > 0) then
        do il = 1, snowpack%nlayers

            ! CARMAGNOLA 2013 - alpha = 1E-4
            sph = snowpack%snow(il)%sph
            ! is_dendritic = snowpack%snow(il)%optd < 1E-4 * (4.0 - sph)
            is_dendritic = snowpack%snow(il)%optd < 1E-4 * (4.0 - sph) - eps
            ! endc

            gs = 1E-4 * (4.0 - snowpack%snow(il)%sph)
            ! TODO: use a function insetad
            dendr = (snowpack%snow(il)%optd / 1E-4 - 4.0 + snowpack%snow(il)%sph )/(snowpack%snow(il)%sph - 3.0)

            rho_i = (snowpack%snow(il)%ws + snowpack%snow(il)%wl)/snowpack%snow(il)%dz ! layer density [kg m^-3] liq + solid
            ! rho_i = (snowpack%snow(il)%ws)/snowpack%snow(il)%dz ! layer density [kg m^-3] solid snow only
            ! is_dendritic = (snowpack%snow(il)%d > 0.0001)
            Frho = 1.25 - 0.0042 * (max(rho_min, rho_i) - rho_min)
            if (is_dendritic) then
                M0i = 0.34 * ( 0.75*dendr - 0.5*snowpack%snow(il)%sph + 0.5) + 0.66 * Frho
            else
                M0i = 0.34 * (-0.583*gs - 0.833*snowpack%snow(il)%sph + 0.833) + 0.66 * Frho
            endif

            ! compute the driftabiliy index (Guyormac'h and Merindol, 1998)
            ! SLi > 0 : snowdrifting occurs
            ! SLi = 0: threshold for wind transport of snow
            SLi = -2.868 * exp ( - 0.085 * Ubar) + 1.0 + M0i

            ! wind drift changes the properties of saltating particles
            ! it makes them more round, more fragmented
            ! For snow layer i
            ! characteristic time for snow grain change due to wind drift
            Gamma_i_drift = max(0.0, SLi * exp ( - pseudo_zi / 0.1))
            tau_i = tau_48h/Gamma_i_drift

            ! VIONNET 2012 - dgs and dd deprecated
            if (.not.is_dendritic) then
                ds = dt_hours * (1.0-snowpack%snow(il)%sph)/tau_i
                ! dgs = dt_hours * 5.0*10.0**(-4)/tau_i
                ! dd = 0.0
            else
                ! dd = dt_hours * snowpack%snow(il)%d / 2.0 / tau_i
                ds = dt_hours * (1.0 - snowpack%snow(il)%sph)/tau_i
            endif

            ! ! CARMAGNOLA 2013 - alpha = 1E-4
            ! ! sph = snowpack%snow(il)%s
            ! ! is_dendritic = snowpack%snow(il)%dopt < 1E-4 * (4.0 - sph)
            ! if (.not.is_dendritic) then
            !     ddopt = -2.0 * 1E-4 * sph * dt_hours *( 1.0 - sph ) / tau_i
            ! else
            !     ! dendricity from dopt and s
            !     den = den_from_dopt(sph, dopt)
            !     term_a = den * (sph - 3.0) / 2.0 / tau_i
            !     term_b = (1.0 - sph) / tau_i * (den - 1.0)
            !     ddopt = 1E-4 * dt_hours * ( term_a + term_b )
            ! endif

            ! CARMAGNOLA 2013 - alpha = 1E-4 - REVISED
            ! sph = snowpack%snow(il)%s
            ! is_dendritic = snowpack%snow(il)%dopt < 1E-4 * (4.0 - sph)
            ds = dt_hours * (1.0 - snowpack%snow(il)%sph)/tau_i
            if (.not.is_dendritic) then
                ddopt = -2.0 * 1E-4 * sph * dt_hours *( 1.0 - sph ) / tau_i ! unchanged, was ok
            else
                ! dendricity from dopt and s
                den = den_from_dopt(sph, dopt)
                term_a = den * (sph - 3.0) / 2.0 / tau_i ! was ok
                term_b = (1.0 - sph) / tau_i * (den - 1.0)
                ddopt = 1E-4 * dt_hours * ( term_a + term_b )
            endif

            ! update now density (i.e. vertical dim. now)
            ! note: this modifies the vertical z profile of the snowpack
            drho = dt_hours * (rho_max - rho_i) / tau_i
            new_rho_i = rho_i + drho
            new_rho_i = max(min(rho_max, new_rho_i), rho_min)
            ! apply constraints - soild snow density
            ! snowpack%snow(il)%dz = snowpack%snow(il)%ws/drho
            snowpack%snow(il)%dz = snowpack%snow(il)%ws/new_rho_i

            ! update snow properties
            ! Updated based on Carmagnola 2013
            ! snowpack%snow(il)%sph = min( snowpack%snow(il)%sph + ds, 1.0)
            ! snowpack%snow(il)%d = snowpack%snow(il)%d + dd
            ! snowpack%snow(il)%gs = snowpack%snow(il)%gs + dgs

            snowpack%snow(il)%optd = snowpack%snow(il)%optd + ddopt


            snowpack%snow(il)%sph = max(min( snowpack%snow(il)%sph + ds, 1.0), 0.0)

            ! add checks
            ! if (snowpack%snow(il)%optd < 1E-4) then
            !     error stop "wind drift: computed optical diameter smaller than alpha = 1E-4"
            ! endif
            ! if (ddopt < 0.0) then
            !     write(*,*) "dendritic snow? ", is_dendritic
            !     error stop "wind drift: computed negative increment in optical diameter!"
            ! endif



            if ((snowpack%snow(il)%optd <0.0).or.(snowpack%snow(il)%optd > 1.0)) then ! 1m upper bound ...
                write(*,*) "current optical diameter value: layer , optd = ", il, snowpack%snow(il)%optd
                ! error stop "ERROR snow_wind_drift_C13 in snow_evolution module: dopt out of bounds after snow wind drift calculation"
                call land_error_message("ERROR snow_wind_drift_C13 in snow_evolution module: dopt out of bounds after snow wind drift calculation", FATAL)
            endif
            if ((snowpack%snow(il)%sph <0.0).or.(snowpack%snow(il)%sph > 1.0)) then ! sph bounds
                write(*,*) "current sphericity value: layer , sph = ", il, snowpack%snow(il)%sph
                ! error stop "ERROR snow_wind_drift_C13 in snow_evolution module: sphericity out of bounds after snow wind drift calculation"
                call land_error_message("ERROR snow_wind_drift_C13 in snow_evolution module: sphericity out of bounds after snow wind drift calculation", FATAL)
            endif

            ! update psuedo-depth above next layer
            pseudo_zi = pseudo_zi + snowpack%snow(il)%dz * (3.25 - SLi)

        enddo
    endif

end subroutine snow_wind_drift_C13



!> \Effect of wind drift on snow, following Brun et al., 1997 and Vionnet et al., 2012
subroutine snow_wind_drift(snowpack, dt, Ubar, verbose)

    class(snowpack_t), intent(inout) :: snowpack !< state of snowpack
    real, INTENT(IN) :: dt ! time step [s]
    real, INTENT(IN) :: Ubar ! wind speed ! [m s^-1]
    logical, intent(IN) :: verbose
    integer il
    real, PARAMETER :: rho_min = 50.0 ! kg m^-3
    ! real, PARAMETER :: rho_max = 350.0 ! kg m^-3
    real, PARAMETER :: rho_max = 400.0 ! kg m^-3
    real, PARAMETER :: tau_48h = 48.0 ! 48 hours in [hours]
    real M0i ! mobility index for layer il = potntial for snow erosion
    real SLi ! driftability index
    real rho_i
    real new_rho_i
    real tau_i
    real Gamma_i_drift
    real pseudo_zi
    real dt_hours
    logical is_dendritic
    real term_a, term_b
    real ds, dd, dgs, drho ! variations in snow properties due to wind drift
    real sph
    ! real den, sph
    real Frho
    real gs, dendr ! computed from s and optical diameter
    real dopt, ddopt
    real ddendr
    real max_optd, min_optd

    dt_hours = dt/3600.0

    ! depth_below_surface = 0.0
    pseudo_zi = 0.0
    if(verbose) write(*,*) "Running wind drift routine!"

    if (snowpack%nlayers > 0) then
        do il = 1, snowpack%nlayers

            sph = snowpack%snow(il)%sph
            dendr =  snowpack%snow(il)%dendr
            is_dendritic = dendr > 1E-7


            gs = 1E-4 * (4.0 - snowpack%snow(il)%sph)
            ! // TODO: use a function insetad
            ! dendr = (snowpack%snow(il)%optd / 1E-4 - 4.0 + snowpack%snow(il)%sph )/(snowpack%snow(il)%sph - 3.0)

            rho_i = (snowpack%snow(il)%ws + snowpack%snow(il)%wl)/snowpack%snow(il)%dz ! layer density [kg m^-3] liq + solid
            ! rho_i = (snowpack%snow(il)%ws)/snowpack%snow(il)%dz ! layer density [kg m^-3] solid snow only
            ! is_dendritic = (snowpack%snow(il)%d > 0.0001)
            Frho = 1.25 - 0.0042 * (max(rho_min, rho_i) - rho_min)
            if (is_dendritic) then
                M0i = 0.34 * ( 0.75*dendr - 0.5*sph + 0.5) + 0.66 * Frho
            else
                M0i = 0.34 * (-0.583*gs - 0.833*sph + 0.833) + 0.66 * Frho
            endif

            ! compute the driftabiliy index (Guyormac'h and Merindol, 1998)
            ! SLi > 0 : snowdrifting occurs
            ! SLi = 0: threshold for wind transport of snow
            SLi = -2.868 * exp ( - 0.085 * Ubar) + 1.0 + M0i

            ! wind drift changes the properties of saltating particles
            ! it makes them more round, more fragmented
            ! For snow layer i
            ! characteristic time for snow grain change due to wind drift
            Gamma_i_drift = max(0.0, SLi * exp ( - pseudo_zi / 0.1))
            tau_i = tau_48h/Gamma_i_drift

            ds = dt_hours * (1.0 - sph)/tau_i ! POSITIVE
            if (.not.is_dendritic) then
                ddopt = -2.0 * 1E-4 * sph * dt_hours * (1.0 - sph)/tau_i ! NEGATIVE
            else
                ! dendricity from dopt and s
                ! den = den_from_dopt(sph, dopt)
                term_a = - dendr * (sph - 3.0) / 2.0 / tau_i
                term_b = (1.0 - sph) / tau_i * (dendr - 1.0)
                ddopt = 1E-4 * dt_hours * ( term_a + term_b )
            endif

            ! if actually dendritic, update dendriticy:
            if (snowpack%snow(il)%dendr > 1E-7 ) then
                ddendr = - dendr/2.0/tau_i * dt_hours
            else
                ddendr = 0.0
            endif



            ! update now density (i.e. vertical dim. now)
            ! note: this modifies the vertical z profile of the snowpack
            drho = dt_hours * (rho_max - rho_i) / tau_i
            drho = max(0.0, drho)
            new_rho_i = rho_i + drho
            new_rho_i = max(min(rho_max, new_rho_i), rho_min)
            ! apply constraints - soild snow density
            ! snowpack%snow(il)%dz = snowpack%snow(il)%ws/drho
            snowpack%snow(il)%dz = (snowpack%snow(il)%ws+snowpack%snow(il)%wl)/new_rho_i

            ! update snow properties
            ! Updated based on Carmagnola 2013
            ! snowpack%snow(il)%sph = min( snowpack%snow(il)%sph + ds, 1.0)
            ! snowpack%snow(il)%d = snowpack%snow(il)%d + dd
            ! snowpack%snow(il)%gs = snowpack%snow(il)%gs + dgs

            snowpack%snow(il)%optd = snowpack%snow(il)%optd + ddopt
            ! snowpack%snow(il)%dendr = snowpack%snow(il)%dendr + ddendr

            ! max_optd = 1500*10**(-6) ! meters, diameter
            ! min_optd = 30*10**(-6) ! meters, diameter
            ! snowpack%snow(il)%optd = max(min( snowpack%snow(il)%optd, max_optd), min_optd)


            snowpack%snow(il)%sph = max(min( snowpack%snow(il)%sph + ds, 1.0), 0.0)
            snowpack%snow(il)%dendr = max(min( snowpack%snow(il)%dendr + ddendr, 1.0), 0.0)

            ! add checks
            ! if (snowpack%snow(il)%optd < 1E-4) then
            !     error stop "wind drift: computed optical diameter smaller than alpha = 1E-4"
            ! endif
            ! if (ddopt < 0.0) then
            !     write(*,*) "dendritic snow? ", is_dendritic
            !     error stop "wind drift: computed negative increment in optical diameter!"
            ! endif



            if ((snowpack%snow(il)%optd <0.0).or.(snowpack%snow(il)%optd > 1.0)) then ! 1m upper bound ...
                write(*,*) "current optical diameter value: layer , optd = ", il, snowpack%snow(il)%optd
                ! error stop "ERROR snow_wind_drift in snow_evolution module: dopt out of bounds after snow wind drift calculation"
                call land_error_message("ERROR snow_wind_drift in snow_evolution module: dopt out of bounds after snow wind drift calculation", FATAL)
            endif
            if ((snowpack%snow(il)%sph <0.0).or.(snowpack%snow(il)%sph > 1.0)) then ! sph bounds
                write(*,*) "current sphericity value: layer , sph = ", il, snowpack%snow(il)%sph
                ! error stop "ERROR snow_wind_drift in snow_evolution module: sphericity out of bounds after snow wind drift calculation"
                call land_error_message("ERROR snow_wind_drift in snow_evolution module: sphericity out of bounds after snow wind drift calculation", FATAL)
            endif
            if ((snowpack%snow(il)%dendr <0.0).or.(snowpack%snow(il)%dendr > 1.0)) then ! dendr bounds
                write(*,*) "current dendricity value: layer , dendr = ", il, snowpack%snow(il)%dendr
                ! error stop "ERROR snow_wind_drift in snow_evolution module: dendricity out of bounds after snow wind drift calculation"
                call land_error_message("ERROR snow_wind_drift in snow_evolution module: dendricity out of bounds after snow wind drift calculation", FATAL)
            endif

             ! if (ddopt<0.0) then
             !      error stop "Wind drift: ddopt should be >0!"
             ! endif



            ! update psuedo-depth above next layer
            pseudo_zi = pseudo_zi + snowpack%snow(il)%dz * (3.25 - SLi)

        enddo
    endif

end subroutine snow_wind_drift


!> \Remove solid snow due to evaporation - updated version compatible with lm4p2
!> \ Includes temporary snow deficit on top of snowpack in case sublim exceedes top layer
subroutine snow_sublimation(s, dt, snow_levap, snow_fevap, hfevap, hlevap, dheat_fevap, &
            use_tfreeze_in_grnd_latent, del_T_toplayer, &
            Mg_imp, snow_melt, &
            lswept1, fswept1, hlswept1, hfswept1, &
            subs_m_imp, lost_wc_em, lost_wc_im, thick_enough_for_evap, verbose)
    class(snowpack_t), intent(inout) :: s !< state of snowpack
    real, intent(out) :: lswept1, fswept1, hlswept1, hfswept1 ![kg m^-2]
    real, intent(in) :: snow_levap ! liquid evaporation rate [kg m^-2 s^-1]
    real, intent(in) :: snow_fevap ! solid sublimation rate [kg m^-2 s^-1]
    real, intent(out) :: hfevap, hlevap ! heat released by sublim [and evap], rate  [J m^-2 s^-1]
    real, intent(out) :: dheat_fevap ! corr in heat released = Dc * DT [J m^-2]
    real, intent(in) :: dt ! model time step [s]
    logical, intent(in), optional :: verbose
    real, intent(IN) :: Mg_imp
    logical, intent(IN) :: thick_enough_for_evap
    real, intent(OUT) :: snow_melt,subs_m_imp
    real, intent(out), dimension(ntracers) :: lost_wc_em, lost_wc_im ! mass of tracers lost from the system [mg/m^2]
    real mass_to_subl, current_mass, rho1
    integer il, it
    real mc_fict, del_T_toplayer, temptop, cap0, dheat, initial_snow_depth
    integer new_upper_layer
    real, ALLOCATABLE :: M_layer(:) ! local variable needed for implicit melt
    real init_ws, init_wl, init_T
    real Told_check, old_ws_check, old_density_1, old_density_il
    logical use_tfreeze_in_grnd_latent
    class(snow_layer_type), ALLOCATABLE :: snow_top, snow_temp
    type(snow_layer_type), allocatable :: snow1(:) ! new snow array
    real old_density, old_heat, new_heat, zerot_heat, excess_heat
    real total_mass, dheat_over_cap0
    real init_heat
    logical stay_in_da_loop
    real trial_new_T, trial_old_T
    logical try_to_merge_snow_deficit
    real Cap1
    real addf, wdef, hdef
    real borrowed_ws
    real DT_max, DT_try, excess_e
    real max2add, ener2add
    real excess_e2, excess_e2_cond

    addf = 1.0
    try_to_merge_snow_deficit = .TRUE.
    excess_e2 = 0.0

    initial_snow_depth = s%depth()
    init_heat = s%heat()
    hfevap = 0.0
    lost_wc_em = 0.0
    lost_wc_im = 0.0
    dheat_fevap = 0.0
    lswept1 = 0
    fswept1 = 0
    hlswept1 = 0
    hfswept1 = 0

    ! ---- evaporation and sublimation -----------------------------------------
    if (initial_snow_depth>0) then

         if(is_watch_point()) then
            write(*,*) '#### gl_snow_step_2 - snow_sublimation ### checkpoint 1 ####'
            write(*,*) "Before sublim, nlayers = ", s%nlayers
            __DEBUG4__(hlevap, hfevap, dheat, dheat_fevap)
            write(*,*) "SUBL CHECKPOINT #1 T[1] = ", s%snow(1)%T
            write(*,*) "SUBL CHECKPOINT #1 - SWE, nlayers, snowdef, heatdef = ", s%SWE(), s%nlayers, s%topsnowdeficit, s%topsnowheatdeficit
            if(verbose) write(*,*) "check 1: ws, dz, rho = ", s%snow(1)%ws,s%snow(1)%dz, old_density_1
        endif

        old_ws_check =s%snow(1)%ws
        old_density_1 = s%snow(1)%ws/s%snow(1)%dz
        s%snow(1)%wl = s%snow(1)%wl - snow_levap*dt
        s%snow(1)%ws = s%snow(1)%ws - snow_fevap*dt

        if (s%snow(1)%ws<0) then
            s%topsnowheatdeficit = s%topsnowheatdeficit + CSW*(s%snow(1)%ws-1E-7)*(s%snow(1)%T-TFREEZE)
            s%topsnowdeficit = s%topsnowdeficit + s%snow(1)%ws - 1E-7
            s%snow(1)%ws = 1E-7
        endif

        s%snow(1)%dz = s%snow(1)%ws/old_density_1 ! note can be < 0 if ws < 0 here
        cap0 = clw*s%snow(1)%wl + csw*s%snow(1)%ws
        ! T adjustment for nonlinear terms (del_T)*(del_W)
        !   dheat = delta_time*(clw*snow_levap+csw*snow_fevap)*del_T(1)
        dheat = dt*(clw*snow_levap+csw*snow_fevap)*del_T_toplayer !
        dheat_over_cap0 = dheat/cap0

        ! take out extra heat not claimed in advance for evaporation
        if (use_tfreeze_in_grnd_latent) dheat = dheat &
            - dt*((cpw-clw)*snow_levap+(cpw-csw)*snow_fevap) &
                               *(s%snow(1)%T-del_T_toplayer-tfreeze)
        hfevap = snow_fevap*CSW*(s%snow(1)%T-TFREEZE) ! for output /check energy balance only
        hlevap = snow_levap*CLW*(s%snow(1)%T-TFREEZE) ! for output /check energy balance only
        Told_check = s%snow(1)%T
        ! s%snow(1)%T  = s%snow(1)%T  + dheat_over_cap0

        !-------------------
        ! STORE EXCESS ENERGY TO AVOID HIGH/LOO T IN TOP LAYER DUE TO SUBL
        DT_max = 0.0 ! max T change in a single time step
        DT_try = dheat_over_cap0
        if (DT_try > DT_max) then
            ! write(*,*) "CASE A"
            s%snow(1)%T  = s%snow(1)%T  + DT_max
            excess_e = (DT_try - DT_max)*cap0 ! assign excess energy to topsnowdeficit, switch sign
            ! s%topsnowheatdeficit = s%topsnowheatdeficit + excess_e
            excess_e2 = excess_e
        else if (DT_try < -DT_max) then
            ! write(*,*) "CASE B"
            s%snow(1)%T  = s%snow(1)%T  - DT_max
            excess_e = (DT_try + DT_max)*cap0 ! assign excess energy to topsnowdeficit, switch sign
            ! s%topsnowheatdeficit = s%topsnowheatdeficit + excess_e
            excess_e2 = excess_e
        else
            ! write(*,*) "CASE C"
            s%snow(1)%T  = s%snow(1)%T  + dheat_over_cap0
            excess_e2 = 0.0
        endif
        ! write(*,*) "cap0, dheat, Told, Tnew, excess_e2 = ", cap0 ,dheat, Told_check, s%snow(1)%T, excess_e2
        ! write(*,*) "Heat after = ", s%heat()
        !-------------------

          dheat_fevap = dheat


        !!!!! -------- NOW DO IMPLICIT MELT OR FREEZE ------------
        if(is_watch_point()) then
            write(*,*) '#### gl_snow_step_2 - snow_sublimation ### checkpoint 2 ####'
            if(s%nlayers>0) write(*,*) "Cap0, dheat, dheat_over_Cap0 = ",cap0, dheat, dheat_over_cap0
            if(s%nlayers>0) write(*,*) "SUBL CHECKPOINT #2 T[1] = ", s%snow(1)%T
            if(s%nlayers>0)  write(*,*) "SUBL CHECKPOINT #3 ws[1], wl[1] = ", s%snow(1)%ws, s%snow(1)%wl
            write(*,*) "SUBL CHECKPOINT #2 - SWE, nlayers, snowdef, heatdef = ", s%SWE(), s%nlayers, s%topsnowdeficit, s%topsnowheatdeficit
        endif
        allocate(M_layer(s%nlayers))
        if (initial_snow_depth>0) then  ! // TODO remove if, already in this case here surely
            snow_melt = Mg_imp/dt
        else
            snow_melt = 0.0
        endif
        M_layer = 0.0
        subs_M_imp = Mg_imp
        do il = 1, s%nlayers
            if (initial_snow_depth>0 .and. subs_M_imp.gt.0) then ! MELT case, subs_M_imp > 0
                ! M_layer(il) =  min( subs_M_imp, max(0.0,s%snow(il)%ws) )
                M_layer(il) =  min( subs_M_imp, max(0.0,s%snow(il)%ws - 1E-9) ) ! EZSNOW
                subs_M_imp = subs_M_imp - M_layer(il)
            endif
        enddo
        if (initial_snow_depth>0) then ! Case of Freeze, or remaining ! EZEVAP - ASSIGN THIS TO SUBS INSTEAD
            M_layer(1) = M_layer(1) + subs_M_imp
            subs_M_imp = 0
            if ((M_layer(1)>0).and.(s%snow(1)%ws<M_layer(1))) then ! BORROW THE MISSING ICE
                borrowed_ws = M_layer(1)-s%snow(1)%ws + 1E-7 ! positive ice mass
                s%topsnowheatdeficit = s%topsnowheatdeficit + CSW*(-borrowed_ws)*(s%snow(1)%T-TFREEZE)
                s%topsnowdeficit = s%topsnowdeficit - borrowed_ws
                s%snow(1)%ws = s%snow(1)%ws + borrowed_ws
            endif
        endif
        if(is_watch_point()) then
            write(*,*) "Start -> Case of neagtive ws: T, wl, ws, MELT[1] = ", s%snow(1)%T , s%snow(1)%wl, s%snow(1)%ws, M_layer(1)
            write(*,*) "Start -> Case of negative ws: ws*Cs, wl*Cl, wl*Cl - abs(ws)*Cs = ", s%snow(1)%ws*CSW, s%snow(1)%wl*CLW, s%snow(1)%wl*CLW-abs(s%snow(1)%ws*CSW)
        endif
    do il = 1, s%nlayers
        if (initial_snow_depth>0) then
            old_density_il = s%snow(il)%ws/s%snow(il)%dz ! original layer density
            cap0 = s%snow(il)%hCap() ! original heat capacity of layer
            init_wl =s%snow(il)%wl
            init_ws =s%snow(il)%ws
            init_T =s%snow(il)%T
            s%snow(il)%wl = s%snow(il)%wl + M_layer(il) ! melt if positive, freeze if negative
            s%snow(il)%ws = s%snow(il)%ws - M_layer(il)
            s%snow(il)%dz = s%snow(il)%ws/old_density_il ! maintain original density of the layer -> shrink layer thickness
            s%snow(il)%T  = TFREEZE + (cap0*(s%snow(il)%T-TFREEZE) ) &
                                                    / ( cap0 + (CLW-CSW)*M_layer(il) )
            if(is_watch_point() .and.(il==1)) then
                write(*,*) "End Case of negative ws: ws, wl, -abs(ws)+wl = ", s%snow(1)%ws, s%snow(1)%wl, -abs(s%snow(1)%ws)+s%snow(1)%wl
                write(*,*) "End Case of negative ws: ws*Cs, wl*Cl, wl*Cl - abs(ws)*Cs = ", s%snow(1)%ws*CSW, s%snow(1)%wl*CLW, s%snow(1)%wl*CLW-abs(s%snow(1)%ws*CSW)
                write(*,*) "End Case of neagtive ws: new T, Cap0, Cap1 = ", s%snow(1)%T , Cap0, Cap1
                write(*,*) "Tnew = ", s%snow(1)%T
            endif
        endif
    enddo
    DEALLOCATE(M_layer)
    !!!!! ------ END IMPLICIT MELT --------

    if(is_watch_point()) then
        write(*,*) '#### gl_snow_step_2 - snow_sublimation ### checkpoint 3 ####'
        if(s%nlayers>0)  write(*,*) "SUBL CHECKPOINT #3 T[1] = ", s%snow(1)%T
        if(s%nlayers>0)  write(*,*) "SUBL CHECKPOINT #3 ws[1], wl[1] = ", s%snow(1)%ws, s%snow(1)%wl
        write(*,*) "SUBL CHECKPOINT #3 - SWE, nlayers, snowdef, heatdef = ", s%SWE(), s%nlayers, s%topsnowdeficit, s%topsnowheatdeficit
    endif
    ! if top layer has negative ice, remove it and create top snow layer
    ! do it also if the first layer is too small
    ! so as to avoid issue that if too much sublimation, too low T1 due to dheat budget
    if (s%snow(1)%ws < 0.0) then ! changed_eps
        if (s%snow(1)%wl + s%snow(1)%ws > 0) then ! ASSIGN ALL TO TOPWATER
            s%topwater = s%topwater + s%snow(1)%wl + s%snow(1)%ws
            s%topwheat = s%topwheat + CLW*s%snow(1)%wl*(s%snow(1)%T-TFREEZE) + HLF*s%snow(1)%wl + CSW*s%snow(1)%ws*(s%snow(1)%T-TFREEZE)
        else ! ASSIGN ALL TO TOPSNOWDEFICIT
            s%topsnowdeficit = s%topsnowdeficit + s%snow(1)%wl + s%snow(1)%ws
            s%topsnowheatdeficit = s%topsnowheatdeficit + CLW*s%snow(1)%wl*(s%snow(1)%T-TFREEZE) + HLF*s%snow(1)%wl + CSW*s%snow(1)%ws*(s%snow(1)%T-TFREEZE)
        endif

        ! remove layer now
        if (s%nlayers > 1) then
            ! first pass any tracers to layer below, then remove 1st layer
            do it = 1, NTRACERS
                s%snow(2)%wc_im(it) = s%snow(2)%wc_im(it) + s%snow(1)%wc_im(it)
                s%snow(2)%wc_em(it) = s%snow(2)%wc_em(it) + s%snow(1)%wc_em(it)
            enddo
            s%snow(1:s%nlayers-1) = s%snow(2:s%nlayers)
            s%nlayers = s%nlayers - 1 ! in both cases
        else if (s%nlayers == 1) then ! case only one layer  since snowpack must be activve
            do it = 1, NTRACERS
                lost_wc_em(it) = lost_wc_em(it) + s%snow(1)%wc_em(it)
                lost_wc_im(it) = lost_wc_im(it) + s%snow(1)%wc_im(it)
            enddo
            deallocate(s%snow)
            s%nlayers = 0
        else
             call land_error_message("ERROR snow_sublimation in snow_evolution module: The number of layers should not be zero here!", FATAL)
        endif

    ! EZDEV - RESTART HERE
    if(is_watch_point()) then
        write(*,*) '#### gl_snow_step_2 - snow_sublimation ### checkpoint 4 ####'
        if (s%nlayers>0)  write(*,*) "SUBL CHECKPOINT #4 T[1] = ", s%snow(1)%T
        if (s%nlayers>0)  write(*,*) "SUBL CHECKPOINT #4 ws[1], wl[1] = ", s%snow(1)%ws, s%snow(1)%wl
         write(*,*) "SUBL CHECKPOINT #4 - SWE, nlayers, snowdef, heatdef = ", s%SWE(), s%nlayers, s%topsnowdeficit, s%topsnowheatdeficit
    endif
    ! //TODO: This was added to increase numerical stability and reduce low T occurrence
    else if ((s%snow(1)%ws < 1E-2)) then
        ! add it to the topwater pool instead
        ! get rid of 1st layer and add it to the second layer
        ! call merge_layers(s%snow(1), s%snow(2))

        s%topwater = s%topwater + s%snow(1)%ws + s%snow(1)%wl
        s%topwheat = s%topwheat + s%snow(1)%ws*CSW*(s%snow(1)%T-TFREEZE)
        s%topwheat = s%topwheat + s%snow(1)%wl*CLW*(s%snow(1)%T-TFREEZE) + s%snow(1)%wl*HLF

        if (s%nlayers > 1) then
            ! first pass any tracers to layer below, then remove 1st layer
            do it = 1, NTRACERS
                s%snow(2)%wc_im(it) = s%snow(2)%wc_im(it) + s%snow(1)%wc_im(it)
                s%snow(2)%wc_em(it) = s%snow(2)%wc_em(it) + s%snow(1)%wc_em(it)
            enddo
            s%snow(1:s%nlayers-1) = s%snow(2:s%nlayers)
            s%nlayers = s%nlayers - 1
        else if (s%nlayers == 1) then
            do it = 1, NTRACERS
                lost_wc_em(it) = lost_wc_em(it) + s%snow(1)%wc_em(it)
                lost_wc_im(it) = lost_wc_im(it) + s%snow(1)%wc_im(it)
            enddo
            s%nlayers = 0
            deallocate(s%snow)
        endif
    endif


    excess_e2 = excess_e2 + s%topsnowheatdeficit
    s%topsnowheatdeficit = 0.0
    excess_e2_cond = excess_e2
    if (.not.(s%nlayers>0)) then
        s%topsnowheatdeficit = s%topsnowheatdeficit + excess_e2
        excess_e2 = 0.0
    else
        if (excess_e2_cond > 0.0) then
            do il=1,s%nlayers
                if ((excess_e2 > 0.0) .and. (s%snow(il)%T-TFREEZE < 0.0) .and. (s%snow(il)%hCap()>0) ) then ! case Tdef > TF, Ti < TF
                    max2add =  (s%snow(il)%T-TFREEZE) * s%snow(il)%hCap()  ! NEG
                    ener2add = min(-max2add, excess_e2)  ! POS
                    excess_e2 = excess_e2 - ener2add
                    s%snow(il)%T = s%snow(il)%T  + ener2add / s%snow(il)%hCap()
                endif
            enddo
            s%topsnowheatdeficit = s%topsnowheatdeficit + excess_e2
            excess_e2 = 0.0
        else
            do il=1,s%nlayers
                if ((excess_e2 < 0.0) .and. (s%snow(il)%T-TFREEZE > 0.0) .and. (s%snow(il)%hCap()>0) ) then ! case Tdef > TF, Ti < TF
                    max2add =  (s%snow(il)%T-TFREEZE) * s%snow(il)%hCap() ! POS
                    ener2add = min(max2add, -excess_e2) ! POS
                    excess_e2 = excess_e2 + ener2add
                    s%snow(il)%T = s%snow(il)%T - ener2add / s%snow(il)%hCap()
                endif
            enddo
            ! add any remaining excess_e2 to top layer, regardless of resulting T
            s%topsnowheatdeficit = s%topsnowheatdeficit + excess_e2
            excess_e2 = 0.0
        endif
    endif

    if (try_to_merge_snow_deficit) then

    if(is_watch_point()) then
        write(*,*) '#### gl_snow_step_2 - snow_sublimation [before merge] ### checkpoint 5 ####'
        if (s%nlayers>0)  write(*,*) "SUBL CHECKPOINT #5 T[1] = ", s%snow(1)%T
        write(*,*) "SUBL CHECKPOINT #5 - SWE, nlayers, snowdef, heatdef = ", s%SWE(), s%nlayers, s%topsnowdeficit, s%topsnowheatdeficit
        call s%print()
    endif
    ! if there is negative ice on top of snow,
    ! attemp to merge it with underlying layers until deficit is filled
    ! if deficit exceedes total ice in the snowpack, retain negative mass on toplayer
    ! to be filled later on by fresh snowfall
    stay_in_da_loop = .true.
    do while (s%nlayers > 0 .and. s%topsnowdeficit < 0.0 .and. stay_in_da_loop)

        if (s%snow(1)%ws > - s%topsnowdeficit*addf ) then ! enough mass to fill deficit

            old_density = (s%snow(1)%ws)/s%snow(1)%dz
            old_heat = s%snow(1)%heat()  ! without top snow, original snowpack layer
            new_heat = old_heat + s%topsnowheatdeficit*addf! new snowpack, with top snow deficit added to it
            zerot_heat = HLF*(s%snow(1)%ws + s%snow(1)%wl + s%topsnowdeficit*addf)
            excess_heat = new_heat - zerot_heat
            if (excess_heat > 0) then
                ! all melted:
                s%snow(1)%wl = s%snow(1)%wl  + s%snow(1)%ws + s%topsnowdeficit*addf
                s%snow(1)%ws = 0.0
                s%snow(1)%T = TFREEZE + excess_heat/(CLW*s%snow(1)%wl)
                s%topsnowdeficit = s%topsnowdeficit*(1-addf)
                s%topsnowheatdeficit = s%topsnowheatdeficit*(1-addf)
            else if (new_heat > 0) then ! mixed phases
                total_mass = s%snow(1)%ws  + s%snow(1)%wl + s%topsnowdeficit*addf
                s%snow(1)%wl = new_heat/HLF
                s%snow(1)%ws = total_mass - s%snow(1)%wl
                s%snow(1)%T = TFREEZE
                s%topsnowdeficit = s%topsnowdeficit*(1-addf)
                s%topsnowheatdeficit = s%topsnowheatdeficit*(1-addf)
            else ! energy < 0, all solid
                trial_old_T = s%snow(1)%T
                trial_new_T = TFREEZE + new_heat/(CSW*(s%snow(1)%ws  + s%snow(1)%wl + s%topsnowdeficit))
                if (trial_new_T < 200.0) then ! do not do the merge
                    ! don't merge
                    stay_in_da_loop = .false.
                else
                    ! do the merging
                    s%snow(1)%ws = s%snow(1)%ws  + s%snow(1)%wl + s%topsnowdeficit*addf
                    s%snow(1)%wl = 0.0
                    s%snow(1)%T = TFREEZE + new_heat/(CSW*s%snow(1)%ws)
                    s%topsnowdeficit = s%topsnowdeficit*(1-addf)
                    s%topsnowheatdeficit = s%topsnowheatdeficit*(1-addf)
                endif
            endif
            ! preserve original layer density, and other properties
            s%snow(1)%dz = s%snow(1)%ws/old_density
        else ! not enough mass in current layer to make up for deficit
            ! add existing layer's mass and heat to toplayer and proceed to next layer
            s%topsnowdeficit = s%topsnowdeficit + s%snow(1)%ws ! this must still be negative here
            if (s%topsnowdeficit>0) call land_error_message("ERROR snow_sublimation in snow_evolution module: topsnowdeficit should still be negative here!!", FATAL)
            s%topsnowheatdeficit = s%topsnowheatdeficit + s%snow(1)%ws*CSW*(s%snow(1)%T-TFREEZE) ! is it ok summing energy to energy deficit?
            s%topwater = s%topwater + s%snow(1)%wl
            s%topwheat = s%topwheat + s%snow(1)%wl*CLW*(s%snow(1)%T-TFREEZE) + HLF*s%snow(1)%wl
            ! remove current layer from stack and pass any tracers to layer below
            if (s%nlayers>1) then
                do it = 1, NTRACERS
                    s%snow(2)%wc_im(it) = s%snow(2)%wc_im(it) + s%snow(1)%wc_im(it)
                    s%snow(2)%wc_em(it) = s%snow(2)%wc_em(it) + s%snow(1)%wc_em(it)
                enddo
                s%snow(1:s%nlayers-1) = s%snow(2:s%nlayers)
                s%nlayers = s%nlayers - 1
            else
                do it = 1, NTRACERS
                    lost_wc_em(it) = lost_wc_em(it) + s%snow(1)%wc_em(it)
                    lost_wc_im(it) = lost_wc_im(it) + s%snow(1)%wc_im(it)
                enddo
                s%nlayers = 0
                deallocate(s%snow)
            endif
        endif
    enddo

    endif


    else ! case of no snow layers
        subs_M_imp = Mg_imp
        snow_melt = 0.0
    endif !! end case of nlayers >0

    if(is_watch_point()) then
        write(*,*) '#### gl_snow_step_2 - snow_sublimation [after merge] ### checkpoint 6 ####'
        if (s%nlayers>0)  write(*,*) "SUBL CHECKPOINT #6 T[1] = ", s%snow(1)%T
        write(*,*) "SUBL CHECKPOINT #6 - SWE, nlayers, snowdef, heatdef = ", s%SWE(), s%nlayers, s%topsnowdeficit, s%topsnowheatdeficit
        call s%print()
    endif

end subroutine snow_sublimation

!> \add solid precipitation to the snowpack
subroutine snow_solid_balance(s, fprec, fevap, lprec, levap, tprec, wetdep, drydep, &
                              Ubar, Tatm, lost_wc_em, lost_wc_im, dt, verbose_in)

    class(snowpack_t), intent(inout) :: s !< state of snowpack
    real, intent(in) :: fevap ! snow sublimation rate [kg m^-2 s^-1]
    real, intent(in) :: fprec ! solid precipitation rate [kg m^-2 s^-1]
    real, intent(in) :: levap ! liquid evaporation
    real, INTENT(IN) :: tprec ! temperature of precipitation [K]
    real, INTENT(IN) :: lprec ! temperature of precipitation [K]
    real, INTENT(IN) :: Ubar ! average wind speed [m s^-1]
    real, INTENT(IN) :: Tatm ! atmos temperature [K]
    real, intent(IN) :: dt ! time step [s]
    logical, intent(in), OPTIONAL :: verbose_in
    real, intent(out), dimension(NTRACERS) :: lost_wc_em, lost_wc_im ! [mg/m2]
    real, intent(IN) :: wetdep(NTRACERS) ! wet deposition of tracers from atmosphere [ppm]
    real, intent(IN) :: drydep(NTRACERS) ! wet deposition of tracers from atmosphere [mg m^-2 s^-1]
    real :: new_snow_depth
    real :: T_new_snow
    integer :: il, it ! counters (layers, tracers)
    real wetdepf(NTRACERS) ! wet deposition IN SNOW ONLY of tracers from atmosphere [mg m^-2 s^-1]
    real mass_to_subl, current_mass
    type(snow_layer_type) snow0 ! new snow instance of size 1
    type(snow_layer_type), allocatable :: snow1(:) ! new snow array
    real s_fall ! sphericity of new snow [adim, in (0,1)]
    real d_fall ! dendricity of new snow [adim, in (0,1)]
    real dopt_fall ! optical diameter of new snow [m]
    real rho_fall ! density of fresh snow [kg m^-3]
    integer n_new_layers
    integer min_nlayers_single_event
    integer nlpnl
    real rho1
    real topdz
    real evap_left, my_evap
    real fprec2 ! precipitation left after filling top snow deficit
    logical verbose

    if (.not.PRESENT(verbose_in))  then
        verbose = .FALSE. ! default argument
    else
        verbose = verbose_in
    endif

    lost_wc_em = 0.0
    lost_wc_im = 0.0


    ! add dry deposition of impurities
    if (s%nlayers > 0) then
        s%snow(1)%wc_em = s%snow(1)%wc_em + drydep * dt
    else
        do it =1, NTRACERS
            lost_wc_em(it) = lost_wc_em(it) + drydep(it) * dt ! UNITS [mg/m2]
        enddo
    endif


    ! compute the wet deposition due to solid snow only
    do it =1, NTRACERS
        ! if ((fprec) > 0.0) then
        if ((fprec) > 1E-9) then
            ! wetdepf(it) = wetdep(it) * fprec/(fprec + lprec)
            ! UNITS: [wetdepf] = mg/m2/s
            ! UNITS: [wetdep] = ppm
            ! UNITS: [fprec] = kg/m2/s
            wetdepf(it) = wetdep(it) * fprec
        else
            wetdepf(it) = 0.0
        endif
    enddo

    ! properties of fresh snow ::
    ! use Carmagnola 2013 approach instead of Vionnet et al., 2012
    d_fall = min(max(1.29-0.17*Ubar, 0.20), 1.0) ! Vionnet et al., 2012
    s_fall = min(max(0.08*Ubar+0.38, 0.5), 0.9) ! Vionnet et al., 2012
    dopt_fall = 1E-4
    !    gs_fall = 3.5 * 10.0**(-4) ! 3.5 mm diam of freshly fallen snow
    !    h_fall = .false. ! not yet melt / refrozen
    call new_snow_density(rho_fall, Tatm, Ubar)
    if(verbose) write(*,*) "snowfall: new snow density [kg m^-3] = ", rho_fall

    ! first, try to fill any snow deficit on top with new snow
    ! if all snow was used to fill deficit and there is snow, add wet deposition below
    fprec2 = fprec
    if (fprec>0.0 .and. s%topsnowdeficit < 0.0) then
        if (fprec*dt > - s%topsnowdeficit) then
            ! note: the heat deficit can remain non zero, and in general it will!
            ! water enters at tprec temperature...
            s%topsnowheatdeficit = s%topsnowheatdeficit - CSW*s%topsnowdeficit*(tprec-TFREEZE) ! subtract a negative mass -> neg heat capacity
            fprec2 = fprec +s%topsnowdeficit/dt ! this is the snow rate remaining
            s%topsnowdeficit = 0.0
            if(verbose) write(*,*) "add solid, cover deficit + some left: precip, precip2 = ", fprec, fprec2
        else
            if(verbose) write(*,*) "add solid, cover only deficit"
            s%topsnowheatdeficit = s%topsnowheatdeficit + CSW*fprec*dt*(tprec-TFREEZE)
            fprec2 = 0.0
            s%topsnowdeficit = s%topsnowdeficit + fprec*dt
            ! used it all up, but still deliver impurities to snow
                if (s%nlayers > 0.0) then
                    do it = 1, NTRACERS
                        s%snow(1)%wc_im(it) = s%snow(1)%wc_im(it)  + wetdepf(it) * dt
                    enddo
                else ! if all snow is used up for filling deficit, and no other snow layers
                    do it = 1,NTRACERS
                    lost_wc_im(it) = lost_wc_im(it) + wetdepf(it) * dt
                    enddo
                endif
        endif
    endif
    !  address case of negative solid precipitation
    if (fprec<0.0) then
        fprec2 = 0.0
        s%topsnowdeficit = s%topsnowdeficit + fprec*dt
        s%topsnowheatdeficit = s%topsnowheatdeficit + CSW*fprec*dt*(tprec-TFREEZE)
    endif

    ! now add the new fresh snow to the top of the remaining snowpack
    if(verbose) write(*,*) "fprec left, ", fprec2
    if(verbose) write(*,*) "rho_fall, ", rho_fall
    if(verbose) write(*,*) "dt, ", dt
    new_snow_depth = fprec2 * dt / rho_fall

    if(verbose) write(*,*) "snowfall: old snow depth = ", s%depth()
    if(verbose) write(*,*) "snowfall: new snow depth = ", new_snow_depth
    if(verbose) write(*,*) "snowfall: original snow mass (before deficit) = ", fprec*dt
    if(verbose) write(*,*) "snowfall: new snow mass = ", fprec2*dt

    ! do nothing if fprec = 0 and fevap = 0
    if (s%nlayers > 0) then
        topdz = s%snow(1)%dz
    else
        topdz = 1.0
    endif

    if (.not.new_snow_depth >0) then
        if (verbose) write(*,*) "no new snow left to deposit, leave now"
    else

    !  changed_eps
    if ((s%nlayers > 0).and.(new_snow_depth < 0.5 * topdz).and.(new_snow_depth > 0.0)) then ! changed_eps
        if(verbose) write(*,*) "adding [little] snow, merging it to existing snowpack with nlayers = ", s%nlayers
        !create a new layer, and then merge layers
        snow0%T = tprec
        snow0%dz = new_snow_depth
        snow0%wl = 0.0
        snow0%ws = fprec2 * dt
        do it = 1, NTRACERS
            snow0%wc_im(it) = wetdepf(it) * dt
            snow0%wc_em(it) = 0.0 ! added
        enddo
        snow0%age = 0.0
        snow0%dendr = d_fall
        snow0%optd = dopt_fall
        snow0%sph = s_fall
        call merge_layers(snow0, s%snow(1)) ! merge the new layer snow0 into the top snowpack layer
    else ! snow does not exist on the ground, or a lot of snow
        if (verbose) write(*,*) "either no old snow, or a lot of new snow: create new layers"
        ! in this case add the new snow in a series of new layers
        ! all layers will be created equal : same mass of snow and of tracers
        ! determine number of layers to add to the snowpack
        ! if there are already layers, the minimum number can be reduced up to 1
        ! we want at least 3 layers to solve diffusion eqn (Vionnet et al., 2012)
        if (verbose) write(*,*) "new snow depth = ", new_snow_depth
        if (new_snow_depth > 0.0) then
            ! if there are less than 3 layers of snow, make sure we go to three
            ! else set the minimum number of additional layers to 1.0
            min_nlayers_single_event = int( max(3.0 - real(s%nlayers), 1.0) )
            ! then as in Vionnet et al., 2012, new number of layers bewteen 1/MINVAL and 5
            n_new_layers =  max(min_nlayers_single_event,  min( 5, ceiling(100.0*new_snow_depth)))
            ! allocate new snowpack only if it snows and not already snow on the ground
            ! size of the new array of snow layers - after adding new snow
            nlpnl = n_new_layers + s%nlayers
            allocate(snow1( nlpnl))
            if(verbose) write(*,*) "min nlayers single event", min_nlayers_single_event
            if(verbose) write(*,*) 'n_new_layers =', n_new_layers
            if(verbose) write(*,*) 'current value of s%nlayers =', s%nlayers
            if(verbose) write(*,*) 'current size of snow1 =', size(snow1)
            if(verbose) write(*,*) 'current size of s%snow =', size(s%snow)

            if (n_new_layers == 0) call land_error_message("ERROR snow_solid_balance in snow_evolution module: n_new_layers should not be zero!", FATAL)
            if (verbose) write(*,*) "case no old snow -> create fresh snow layers"
            do il = 1, n_new_layers
                do it =1, NTRACERS ! subdivide equally - all new layers are equal
                    snow1(il)%wc_im(it) = wetdepf(it) / real(n_new_layers) * dt ! if new layers, all are created equal
                    snow1(il)%wc_em(it) = 0.0 ! added
                enddo
                snow1(il)%age = 0.0 ! in days, real value
                snow1(il)%sph = s_fall
                snow1(il)%dendr = d_fall
                snow1(il)%optd = dopt_fall ! fresh snow value in Carmagnola et al., 2013
                snow1(il)%dz = new_snow_depth/real(n_new_layers) ! new layers all equal sized
                if(verbose) write(*,*) "new dz = ",new_snow_depth/real(n_new_layers)
                snow1(il)%wl = 0.0 ! no water raining in the snow for now - add later
                snow1(il)%ws = rho_fall * snow1(il)%dz  ! store mass in [kg /m^2]
                snow1(il)%T = tprec
            enddo
            ! if there is no snow, allocate new snow layers now:
            if (s%nlayers==0) then
                if (verbose) write(*,*) "no initial snowpack, create one"
                if (ALLOCATED(s%snow)) deallocate(s%snow)
                ALLOCATE(s%snow(n_new_layers))
                s%snow = snow1
                s%nlayers = n_new_layers
            else
                if (verbose) write(*,*) "existing snow present; adding new layers to existing snowpack"
                ! copy the existing snow layers in the bottom part of the new snow column
                snow1(n_new_layers+1:n_new_layers+s%nlayers) = s%snow(1:s%nlayers)
                s%snow = snow1
                s%nlayers = s%nlayers + n_new_layers
            endif
            DEALLOCATE(snow1)
        endif
    endif ! end case in which we add a number of new snow layer due to snowfall
    endif ! case of precip > 0 to deposit

end subroutine snow_solid_balance


!> \add liquid precipitation to the snowpack, do melt and vertical liquid water flow
subroutine snow_liquid_balance(s, lprec, levap, fprec, tprec, wetdep, snow_lprec, &
                                    snow_hlprec, lost_wc_em, lost_wc_im,  dt, verbose_in)

    ! water flushed down leaves instantaneously the snowpack as runoff
    ! the scavenged tracers are also lost
    class(snowpack_t), intent(inout) :: s !< state of snowpack
    real, intent(in) :: lprec ! liquid precipitation rate [kg m^-2 s^-1]
    real, intent(in) :: levap ! liquid evaporation rate [kg m^-2 s^-1]
    real, intent(in) :: fprec ! solid precipitation rate [kg m^-2 s^-1]
    real, intent(in) :: tprec ! temperature of precipitation [K]
    real, intent(in) :: dt ! time step [s]
    real, intent(out), dimension(NTRACERS) :: lost_wc_em, lost_wc_im ! [mg/m2]
    real, intent(in) :: wetdep(NTRACERS) ! wet deposition of tracers from atmosphere [ppm]
    real, intent(out) :: snow_lprec, snow_hlprec ! rates, heat wrt liquid at TF in LM4p2
    logical, intent(in), optional :: verbose_in
    real :: wl_excess
    real :: zflux_wl, zflux_T ! verical mass of liquid water [Kg m^-2] moved down the snowpack
    real :: zflux_wc_em(NTRACERS), zflux_wc_im(NTRACERS) ! flux scavenged for each im or em
    integer :: il, it ! counter
    real SWE_il ! snow water equivalent of layer il [kg m^-2]
    integer n_melt_layers, new_layer_counter
    real wetdepl(NTRACERS) ! wet deposition of tracers from atmosphere [mg m^-2 s^-1]
    real wl_max
    real delta_wl
    real rho_snow_il
    real wl_excess_heat
    integer orig_n_layers
    type(snow_layer_type), allocatable :: snow1(:) ! new snow array
    real zflux_ws, old_rho1, zflux_Tsol
    real eps_kill_layer
    real eps_water ! set it to this val to avoid negative values
    logical verbose

    if (.not.PRESENT(verbose_in)) then
        verbose = .FALSE. ! default argument
    else
        verbose = verbose_in
    endif
    eps_kill_layer = 1E-8 ! solid ice mass threhold to remove a layer from the pack

    ! init fluxes of liquid water, heat and impurities percolating down to substrate
    snow_lprec = 0.0
    snow_hlprec = 0.0
    lost_wc_em = 0.0
    lost_wc_im = 0.0

    if (allocated(snow1)) DEALLOCATE(snow1)
    n_melt_layers = 0 ! init number of completely melted layers

    ! compute the fraction of wet deposition carried by liquid precipitation only
    do it = 1, NTRACERS
        if ((lprec) > 1E-9) then
            ! wetdepl(it) = wetdep(it) * lprec / (lprec + fprec)
            ! UNITS: [wetdepf] = mg/m2/s
            ! UNITS: [wetdep] = ppm
            ! UNITS: [fprec] = kg/m2/s
            wetdepl(it) = wetdep(it) * lprec
        else
            wetdepl(it) = 0.0
        endif
    enddo

    ! write(*,*) "liquid balance, test initial LAI content = ", sum(s%lai_em() + s%lai_im() + lost_wc_em + lost_wc_im  )
    ! write(*,*) "liquid balance, test initial LAI content = ", sum(s%lai_em() + s%lai_im() + lost_wc_em + lost_wc_im  + wetdepl*dt)
    ! write(*,*) "Before liquid balance, nlayers = ", s%nlayers

    if (s%nlayers==0) then ! if no snow, get rid immediately of all wet deposited LAIs
        do it = 1, NTRACERS
        lost_wc_im(it) = lost_wc_im(it) + wetdepl(it)*dt ! if there is no snow, flush tracers away
        enddo
    endif

    if(is_watch_point()) then
        write(*,*)'#### Snow step 2 : snow_liquid_balance, initial checkpoint [1] ####'
        write(*,*) "liquid balance, test initial LAI content = ", sum(s%lai_em() + s%lai_im() + lost_wc_em + lost_wc_im  )
        write(*,*) "liquid balance, test initial LAI content = ", sum(s%lai_em() + s%lai_im() + lost_wc_em + lost_wc_im  + wetdepl*dt)
        write(*,*) "Before liquid balance, nlayers = ", s%nlayers
        call s%print()
    endif

    if (s%nlayers > 0) then
        do il = 1, s%nlayers ! loop on layers, do liquid water balance in each
            if (il == 1) then ! water balance for top snow layer
                if (lprec > 0.0) then
                    call add_liquid_to_layer(s%snow(il), lprec*dt, tprec) ! add precip to 1st layer
                endif
                s%snow(il)%wc_im = s%snow(il)%wc_im + wetdepl*dt   ! add IM tracers to 1st layer
                if(is_watch_point()) then
                    write(*,*)'#### Snow step 2 : snow_liquid_balance, inter checkpoint [1.5] ####'
                    call s%print()
                endif

                if (s%topwater > 0.0) then ! also add existing snow topwater to 1st layer if any
                    call add_liquid_to_layer(s%snow(il), s%topwater, TFREEZE+(s%topwheat-s%topwater*HLF)/s%topwater/CLW)
                    s%topwater = 0.0 ! this can be nonzero only between evaporation and liquid balance
                    s%topwheat = 0.0 ! this can be nonzero only between evaporation and liquid balance
                endif

                if(is_watch_point()) then
                    write(*,*)'#### Snow step 2 : snow_liquid_balance, second checkpoint [2] ####'
                    write(*,*) "1. snowpack heat now is = ", s%heat() - s%topwater*HLF
                    write(*,*) "1. snowpack heat now is = ", s%heat()
                    write(*,*) "snow liquid balance: state of snowpack after adding liq precip:"
                    write(*,*) "SWE after adding precip to 1st layer: lprec*dt = ", s%SWE()
                    call s%print()
                endif

            else ! layer other than the top layer
                ! add fluxes of impurities from layers above
                s%snow(il)%wc_im = s%snow(il)%wc_im + zflux_wc_im   ! vector dim=3
                s%snow(il)%wc_em = s%snow(il)%wc_em + zflux_wc_em   ! vector dim=3

                if (zflux_wl > 0.0) then
                    call add_liquid_to_layer(s%snow(il), zflux_wl, zflux_T) ! add water flushed from above
                    ! add any impurities scavenged from above [for now, keep externally mixed]
                    ! s%snow(il)%wc_im = s%snow(il)%wc_im + zflux_wc_im   ! vector dim=3
                    ! s%snow(il)%wc_em = s%snow(il)%wc_em + zflux_wc_em   ! vector dim=3

                endif
                if (zflux_ws > 0.0) then
                    ! add also solid if any
                    ! if (zflux_ws>0.0) then
                    !     snowl1 = s%snow(il)
                    !     old_rho1 = snowl1%ws/snowl1%dz
                    !     snowl1%ws = zflux_ws
                    !     snowl1%T = zflux_T
                    !     snowl1%wl = 0.0
                    !     snowl1%dz = snowl1%ws/old_rho1
                    !     call merge_layers(snowl1, s%snow(il))
                    ! endif

                    ! add solid from above too as cooled liquid with same total heat content
                    ! zflux_Tsol = TFREEZE + (CSW*(zflux_T-TFREEZE)+HLF)/CLW
                    zflux_Tsol = TFREEZE + (CSW*(zflux_T-TFREEZE)-HLF)/CLW
                    call add_liquid_to_layer(s%snow(il), zflux_ws, zflux_Tsol)
                    ! // note in this case layer above was depleted
                    ! so we should add any remaining LAIs
                    ! but it is done at the end
                endif
            endif

            ! if a layer is completely melted after adding liquid, remove it from the stack
            if (s%snow(il)%ws < eps_kill_layer) then ! changed_eps
                n_melt_layers = n_melt_layers + 1
                zflux_wl = s%snow(il)%wl
                zflux_T = s%snow(il)%T
                zflux_wc_im = s%snow(il)%wc_im   ! vector dim=3
                zflux_wc_em = s%snow(il)%wc_em   ! vector dim=3
                s%snow(il)%wc_im = 0.0
                s%snow(il)%wc_em = 0.0
                s%snow(il)%wl = 0.0
                zflux_ws = s%snow(il)%ws
                s%snow(il)%ws = 0.0
                ! I need to remove the il-th layer from allocated snow array. Done at the end
                ! what happens if it is the bottom layer to melt completely?
                ! nothing, just get rid of water as runoff, and tracers will disappear
            else !case layer is still there:
                ! compute pore space only now after updating the ws mass due to any melt / freeze
                rho_snow_il = s%snow(il)%ws/s%snow(il)%dz ! density of solid snow in current layer !
                wl_max = compute_wlmax(s%snow(il)%dz, rho_snow_il)
                wl_excess = s%snow(il)%wl - wl_max   ! excess above available liquid water storage
                zflux_ws = 0.0
                zflux_T = s%snow(il)%T

                if (wl_excess > 0.0) then ! now any excess water is flushed down the snowpack.
                    zflux_wl = wl_excess
                    s%snow(il)%wl = wl_max ! retain pore space capacity full of water (*)
                else
                    zflux_wl = 0.0 ! snow layer pore space is not saturated - nothing to flush
                endif

                ! compute amounts of tracers flushed down with the water
                do it = 1, NTRACERS
                    SWE_il = s%snow(il)%ws + s%snow(il)%wl ! layer snow water equivalent
                    ! note the denominator, since as this point zfluz_wl was already removed in (*)
                    ! scavenge proportionally from internally and externally mixed impurities
                    zflux_wc_em(it) = s%snow(il)%wc_em(it) * SCAVENG(it) * zflux_wl / (SWE_il + zflux_wl)
                    zflux_wc_im(it) = s%snow(il)%wc_im(it) * SCAVENG(it) * zflux_wl / (SWE_il + zflux_wl)
                    s%snow(il)%wc_im(it) = s%snow(il)%wc_im(it) - zflux_wc_im(it) ! lose tracers headed down (snow or ground)
                    s%snow(il)%wc_em(it) = s%snow(il)%wc_em(it) - zflux_wc_em(it) ! lose tracers headed down (snow or ground)
                enddo
            endif

            ! now we reached bottom of snowpack - we remove excess liquid water as runoff
            if (il==s%nlayers) then
                ! // TODO no need to init and incerement here, done only for bottom layer
                snow_lprec = snow_lprec + zflux_wl/dt ! rate
                snow_hlprec = snow_hlprec + zflux_wl/dt*CLW*(zflux_T-TFREEZE) + zflux_wl/dt*HLF ! // with resp to liquid, rate
                snow_lprec = snow_lprec + zflux_ws/dt ! rate
                snow_hlprec = snow_hlprec + zflux_ws/dt*CSW*(zflux_T-TFREEZE)
                ! DO FLUSH THE TRACERS HERE ::
                do it = 1, NTRACERS
                    SWE_il = s%snow(il)%ws + s%snow(il)%wl ! layer snow water equivalent
                    lost_wc_em(it) = lost_wc_em(it) + zflux_wc_em(it)
                    lost_wc_im(it) = lost_wc_im(it) + zflux_wc_im(it)
                enddo
            endif
        enddo ! end loop on layers

        if(is_watch_point()) then
            write(*,*)'#### Snow step 2 : snow_liquid_balance, third checkpoint [3] ####'
            write(*,*) "intermediate check: SWE = ", s%SWE()
            write(*,*) "intermediate check: heat = ", s%heat()
            write(*,*) "state of snowpack before removing empty layers ::"
            write(*,*) "number of layers, numbe of layers to melt = ", s%nlayers, n_melt_layers
            if (s%nlayers > 0) write(*,*) "First layer: ws[1], wl[1] = ", s%snow(1)%ws, s%snow(1)%wl
            call s%print()
        endif

        if (n_melt_layers == s%nlayers) then ! melting all snow layers
            lost_wc_em = lost_wc_em + s%lai_em()
            lost_wc_im = lost_wc_im + s%lai_im()
            s%nlayers = 0
            DEALLOCATE(s%snow)
        else if ((n_melt_layers > 0).and.(n_melt_layers < s%nlayers)) then
            new_layer_counter = 1
            orig_n_layers = s%nlayers
            s%nlayers = orig_n_layers - n_melt_layers
            ALLOCATE(snow1(orig_n_layers - n_melt_layers))
            do il=1, orig_n_layers
                if (s%snow(il)%ws >= eps_kill_layer) then
                    snow1(new_layer_counter) = s%snow(il)
                    new_layer_counter = new_layer_counter + 1
                endif
            enddo

            s%snow(1:s%nlayers) = snow1
            deallocate(snow1)

        else if (n_melt_layers > s%nlayers) then
            write(*,*) "orig_n_layers, n_melt_layers = ", orig_n_layers, n_melt_layers
            call land_error_message( "ERROR snow_liquid_balance in snow_evolution module: Too many layers to melt!", FATAL)
        else if (n_melt_layers==0) then ! ok, no layers to melt
        else
            write(*,*) "orig_n_layers, n_melt_layers = ", orig_n_layers, n_melt_layers
            call land_error_message("ERROR snow_liquid_balance in snow_evolution module: wrong number of layers to melt!", FATAL)
        endif

        if(is_watch_point()) then
            write(*,*)'#### Snow step 2 : snow_liquid_balance, fourth checkpoint [4] ####'
            write(*,*) "state of snowpack after removing empty layers ::"
            write(*,*) "after removing empty layer, check: SWE = ", s%SWE()
            call s%print()
        endif

    else ! case of no snowpack (nlayers=0): Just send any rainfall directly to runoff, and any reamining topwater from sublimated layers
        snow_lprec = snow_lprec + lprec + s%topwater/dt ! rate
        snow_hlprec = snow_hlprec + s%topwheat/dt + lprec*CLW*(tprec-TFREEZE) + lprec*HLF ! // snowpack definition
        s%topwater = 0.0
        s%topwheat = 0.0
    endif

    ! additionally, if prec < 0 must also be passed down as it was not added to the snowpack
    !  but it's already included in snow_lprec if there are zero snow layers
    if ((s%nlayers>0).and.(lprec<0)) then
        snow_lprec = snow_lprec + lprec
        snow_hlprec = snow_hlprec + lprec*CLW*(tprec-TFREEZE) + lprec*HLF
    endif

    if (s%nlayers > 0) then
        if (s%snow(1)%dz<0.0) then
            write(*,*) "liquid balance ends with negative dz = ",s%snow(1)%dz
            call land_error_message("ERROR snow_liquid_balance in snow_evolution module: liquid balance ends with negative layer thickness", FATAL)
        endif
    endif

end subroutine snow_liquid_balance



!> \compute snow compaction for a snow layer
subroutine snow_compaction(s, dt, verbose)

   class(snowpack_t), intent(inout) :: s !< state of snowpack
   real, intent(in) :: dt ! time step
   logical, intent(in) :: verbose
   real :: delta_depth ! change in depth of current snow layer [m]
   real rho ! snow density [kg m^-3]
   real gs ! snow grain size [m]
   real eta ! snow viscosity [kg s^-1 m^-1]
   real :: sigma ! normal; stress from weight of snow in above layers
   real :: f1, f2 ! correction factors to adjust viscosity based on snow micro structure
   real :: mass_on_top
   real :: current_layer_mass
   integer il
   real dz_old
   real min_ws

   real, PARAMETER :: eta0 = 7.62237*10**6 ! [kg s^-1 m^-1]
   real, PARAMETER :: a_eta = 0.1 ! [K^-1]
   real, PARAMETER :: b_eta = 0.023 ! [m^3 kg^-1]
   real, PARAMETER :: c_eta = 250.0 ! [kg m^-3]
   real, PARAMETER :: g1 = 0.4 ! [mm]
   real, PARAMETER :: g2 = 0.2 ! [mm]
   real, PARAMETER :: g3 = 0.1 ! [mm]



  ! modification:
  ! to avoid numerical problems, do not compact snow if the layer is very very small
  ! it will be merged when relayering  the snowpack
  min_ws = 1E-3 ! relayer only if more mass than this low threshold [kg/m^2]


   mass_on_top = 0.0 ! weight above snow layer il
   if (s%nlayers > 0) then
    do il = 1, s%nlayers

            gs = 1E-4 * (4.0 - s%snow(il)%sph) ! using Carmagnola 2013
            ! gs = s%snow(il)%optd
            current_layer_mass = s%snow(il)%ws + s%snow(il)%wl
            sigma = GRAV * (mass_on_top + 0.5 * current_layer_mass)
            mass_on_top = mass_on_top + current_layer_mass

            f1 = (1.0 + 60.0* s%snow(il)%wl / rho_water / s%snow(il)%dz)**(-1)
            f2 = min(4.0, exp(min(g1, max(0.0, gs*1000.0-g2))/g3)) ! gs*1000 in [mm] here

            ! snow viscosity
            rho = (s%snow(il)%ws + s%snow(il)%wl) / s%snow(il)%dz ! density of solid + liquid components
            eta = f1*f2*eta0*rho/c_eta*exp(a_eta*(TFREEZE-s%snow(il)%T) + b_eta*rho)
            delta_depth = - sigma * dt / eta ! [dimensionless] snow layer deformation DL/L
            if(verbose) write(*,*) "snow compaction :: delta_depth = ", delta_depth
            dz_old = s%snow(il)%dz

           if(is_watch_point()) then
              write(*,*)'#### Snow step 2 : snow_compaction, layer il = ', il
              __DEBUG4__(rho, eta, delta_depth, dz_old)
          endif

            if (s%snow(il)%ws > min_ws) then
                s%snow(il)%dz = s%snow(il)%dz * (1 + delta_depth)
                ! only in this case do check that density makes sense
                if ((rho < 0.0).or.(rho > rho_ice)) then
                    write(*,*) "rho = ", rho
                    write(*,*) "eta = ", eta
                    write(*,*) "sph = ", s%snow(il)%sph
                    write(*,*) "sigma = ", sigma
                    write(*,*) "T = ", s%snow(il)%T
                    write(*,*) "dz_old = ", dz_old
                    write(*,*) "dz_new = ", s%snow(il)%dz
                    write(*,*) "ws = ", s%snow(il)%ws
                    write(*,*) "wl = ", s%snow(il)%wl
                    write(*,*) "delta_depth = ", delta_depth
                    call s%print()
                    call land_error_message("ERROR snow_compaction in snow_evolution module: snow density out of bounds after snow compaction calculation!", FATAL)
                endif
            endif

        enddo
    endif

end subroutine snow_compaction



!> \given new temp profile, compute solid melt and / or liquid freeze
subroutine snow_melt_and_freeze(s, dt, snow_lprec, snow_hlprec, lost_wc_em, lost_wc_im, verbose_in)

    class(snowpack_t), intent(inout) :: s !< state of snowpack
    real, intent(in) :: dt ! time step [s]
    real, intent(out),  DIMENSION(NTRACERS) :: lost_wc_im, lost_wc_em ! [mg/m2]
    real, intent(out) :: snow_lprec, snow_hlprec ! rates, heat wrt liquid at TF in LM4p2
    logical, intent(in), optional :: verbose_in
    integer il
    real DTold, DTnew ! temperature (wrt TFREEZE) computed by diff (DT, [K])
                 ! and corrected after accounting for malt/freeze (DTp, [K])
    real melt, freeze, Qsink, Qsource
    logical verbose
    type(snow_layer_type), allocatable :: snow1(:) ! new snow array
    integer n_melt_layers, new_layer_counter, origin_n_layers
    real zflux_wl, zflux_T
    real, DIMENSION(NTRACERS) :: zflux_wc_im, zflux_wc_em
    real hCap0, heat0, hCap1, heat1
    real original_ws
    real rho_start, rho_ends
    real max_freeze, wl_max, rho_snow_il
    logical previous_layer_melted
    integer it

    if (.not.PRESENT(verbose_in)) then
        verbose = .FALSE. ! default argument
    else
        verbose = verbose_in
    endif

    snow_lprec = 0.0
    snow_hlprec = 0.0
    lost_wc_im = 0.0
    lost_wc_em = 0.0
    snow_lprec = 0.0
    snow_hlprec = 0.0

    if (allocated(snow1)) DEALLOCATE(snow1)
    if(verbose) write(*,*) "before melt-freeze: heat, SWE, LIQ, ICE, nlayers= ", s%heat(), s%SWE(), s%liq(), s%ice(), s%nlayers
    n_melt_layers = 0 ! counter for the layers melting
    rho_start = s%density()


    if (s%nlayers > 0) then
    do il = 1, s%nlayers
        DTold = s%snow(il)%T - TFREEZE

        ! condition for snow melt
        if ( ( DTold > 0.0) .and. (s%snow(il)%ws > 0.0)) then
            heat0 = s%snow(il)%heat()
            hCap0 = s%snow(il)%hCap()
            melt = min( s%snow(il)%ws, s%snow(il)%hCap()*DTold/HLF )
            Qsink = melt * HLF ! energy lost by melt latent heat
            original_ws =s%snow(il)%ws
            s%snow(il)%ws = s%snow(il)%ws - melt
            s%snow(il)%wl = s%snow(il)%wl + melt
            if (s%snow(il)%ws < 0.0) call land_error_message("Error in snow_melt_and_freeze in snow_evolution_mod:: Found negative ws value!", FATAL)
            s%snow(il)%dz = s%snow(il)%dz * (s%snow(il)%ws)/original_ws
            ! note :: now wl could exceede the available pore storage
            ! In that case, excess liquid will be flushed now in the next time step
            hCap1 = s%snow(il)%hCap()
            ! compute new temperature of the layer
            ! account for specific heats of both solid and liquid component
            ! if water is there must keep in equil with ice !
            DTnew = (heat0 - s%snow(il)%wl*HLF) / hCap1
            s%snow(il)%T = DTnew + TFREEZE ! new temp in [K]

        else if ( ( DTold < 0.0) .and. (s%snow(il)%wl > 0.0)) then
            max_freeze =s%snow(il)%wl
            freeze = min( max_freeze, -s%snow(il)%hCap()*DTold/HLF )
            Qsource = freeze * HLF
            heat0 = s%snow(il)%heat()
            hCap0 = s%snow(il)%hCap()

            s%snow(il)%ws = s%snow(il)%ws + freeze
            s%snow(il)%wl = s%snow(il)%wl - freeze
            if (s%snow(il)%wl < 0.0) call land_error_message("Error in snow_melt_and_freeze in snow_evolution_mod: wl < 0 value found!", FATAL)
            ! set max layer density after freezing
            if ((s%snow(il)%ws + s%snow(il)%wl)/s%snow(il)%dz > rho_ice) then
                s%snow(il)%dz = (s%snow(il)%ws + s%snow(il)%wl) / (0.95 * rho_ice)
            endif
            hCap1 = s%snow(il)%hCap()
            DTnew = (heat0  - s%snow(il)%wl*HLF) / hCap1
            s%snow(il)%T = DTnew + TFREEZE ! new temp in [K]
        endif

        if ((il > 1).and.(previous_layer_melted)) then ! add any fluxes from above if complete melt happened above
            call add_liquid_to_layer(s%snow(il), zflux_wl, zflux_T)
            s%snow(il)%wc_im = s%snow(il)%wc_im + zflux_wc_im
            s%snow(il)%wc_em = s%snow(il)%wc_em + zflux_wc_em
        endif

        if (s%snow(il)%ws <= eps) then ! was eps
            previous_layer_melted = .TRUE.
            n_melt_layers = n_melt_layers + 1
            zflux_wc_im = s%snow(il)%wc_im ! add impurities to downward flux
            zflux_wc_em = s%snow(il)%wc_em
            s%snow(il)%wc_im = 0.0
            s%snow(il)%wc_em = 0.0
            zflux_wl = s%snow(il)%wl  ! add water to downward flux
            zflux_T = s%snow(il)%T
            s%snow(il)%wl = 0.0
            ! pass excess water and tracers to the layer below
            if(il==s%nlayers) then ! last layer
                if (verbose) write(*,*) "melt runoff updated"
                snow_lprec = snow_lprec + zflux_wl/dt !
                snow_hlprec = snow_hlprec + zflux_wl/dt*CLW*(zflux_T-TFREEZE) + zflux_wl/dt*HLF
                lost_wc_em = lost_wc_em + zflux_wc_em ! added here
                lost_wc_im = lost_wc_im + zflux_wc_im
            endif
        else
            previous_layer_melted = .FALSE.
        endif
    enddo


    ! write(*,*) "M & F: Intermediate bounds check"
    ! call s%check_bounds("check bounds - intermediate melt and freeze ....")
    ! write(*,*) "melt_and_freeze: intermediate check: nlayers, SWE, heat = ", s%nlayers, s%SWE(), s%heat()
    ! call s%print()
    ! now allocate new snow array with only the non-zero layers
    ! write(*,*) "nlayers, n_melt_layers =", s%nlayers, n_melt_layers
    ! write(*,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
    ! call s%print()
    if (n_melt_layers == s%nlayers) then
        s%nlayers = 0
        DEALLOCATE(s%snow)
    else if ((n_melt_layers > 0).and.(n_melt_layers < s%nlayers)) then
        new_layer_counter = 0
        origin_n_layers = s%nlayers
        ALLOCATE(snow1(origin_n_layers - n_melt_layers))
        do il=1, origin_n_layers
            if (s%snow(il)%ws >= eps) then ! changed_eps
                new_layer_counter = new_layer_counter + 1
                snow1(new_layer_counter) = s%snow(il)
            endif
        enddo
        s%nlayers = origin_n_layers - n_melt_layers
        s%snow(1:s%nlayers) = snow1
        deallocate(snow1)
    else if (n_melt_layers==0) then ! no layers to melt, pass
    else
        write(*,*) "nlayers, n_melt_layers =", s%nlayers, n_melt_layers
        call land_error_message("ERROR snow_melt_and_freeze in snow_evolution module: Something wrong with the number of layers to remove!", FATAL)
    endif
    rho_ends = s%density()
    endif ! end case of nlayers >0

    ! write(*,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
    ! call s%check_bounds("check bounds - final melt and freeze ....")
    ! write(*,*) "M & F :: final melt-freeze: heat, SWE, LIQ, ICE, nlayers= ", s%heat(), s%SWE(), s%liq(), s%ice(), s%nlayers
    ! if (verbose)

end subroutine snow_melt_and_freeze


!> \compute snow albedo in lm4p2 using one of the available models
! subroutine gl_compute_snow_albedo(s, snow_T, cosz, on_glacier, p_atm, & ! input
                ! snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis) ! output
subroutine gl_compute_snow_albedo(s, snow_T, cosz, on_glacier, p_atm, subs_refl_dif, & ! input
                                    snow_refl_dir, snow_refl_dif) ! output
  class(snowpack_t), intent(inout) :: s !< state of snowpack
  real, intent(in) :: p_atm  ! ! atm pressure [Pa] from forcing
  real, intent(in) :: snow_T  ! snow temperature, deg K [get it from s instead?]
  real, intent(in) :: cosz ! cosine of zenith angle
  logical, intent(in) :: on_glacier ! TRUE if snow is on glacier
  real, dimension(NBANDS), intent(IN) :: subs_refl_dif
  real, dimension(NBANDS), intent(OUT) :: snow_refl_dir
  real, dimension(NBANDS), intent(OUT) :: snow_refl_dif
!   real, intent(OUT) :: snow_refl_lw, snow_emis

    call s%nearsurf_properties()

    ! snow_emis = 0.0
    ! snow_refl_lw  = 1 - snow_emis

    ! write(*,*) "computing snow albedo:"
    ! write(*,*) "albedo to use = ", albedo_to_use
    ! write(*,*) "albedo correction to use = ", albedo_correction_to_use
    ! write(*,*) "use internal sources = ", use_internal_sources
    ! write(*,*) "checking snowpack nml parameters"
    ! write(*,*) "opt_layer_N = ", opt_layer_N
    ! write(*,*) "opt_layer_R = ", opt_layer_R
    ! write(*,*) "opt_layer_max = ", opt_layer_max

    if (.not. s%nlayers > 0) then
        snow_refl_dir = (/ -9999.9, -9999.9 /)
        snow_refl_dif = (/ -9999.9, -9999.9 /)
        ! snow_emis = -9999.9
        ! snow_refl_lw = -9999.9
        s%beta_rad(BAND_VIS) = -9999.9
        s%beta_rad(BAND_NIR) = -9999.9
    else


    ! write(*,*) "compute snow albedo:"
    ! write(*,*) "albedo_to_use = ", albedo_to_use
    ! write(*,*) "albedo_correction_to_use = ", albedo_correction_to_use

        ! if (albedo_correction_to_use == 'HE') then
        !     call land_error_message("ERROR compute_snow_albedo in snow_evolution module: LAI Albedo correction still needs to be implemented!", FATAL)
        ! endif

        call compute_beta_rad_crocus(s, p_atm) ! first call is used only for the light penetration depth
        ! call compute_albedo_lm4p2(s, cosz, on_glacier) ! first call only for longwave snow propertie
        if (trim(lowercase(albedo_to_use))=='brdf') then
            call compute_albedo_lm4p2(s, snow_T, cosz, on_glacier)
        else if (trim(lowercase(albedo_to_use)) == 'he') then
            call compute_albedo_he(s, cosz) ! only for the penetration depth
        else if (trim(lowercase(albedo_to_use)) == 'crocus') then
            call compute_albedo_crocus(s, p_atm) ! add to it cos dependence through modificed snow grain?
        else if (trim(lowercase(albedo_to_use)) == 'snicar') then
            call compute_snicar_albedo(s, cosz, subs_refl_dif) ! add to it cos dependence through modified snow grain?
        else
            ! error stop "ERROR compute_snow_albedo in snow_evolution module: Must specify a valid albedo model!"
            call land_error_message( "ERROR compute_snow_albedo in snow_evolution module: Must specify a valid albedo model!", FATAL)
        endif
           ! TODO: compute these from crocus regardless of the albedo model chosen
           snow_refl_dif = s%snow_refl_dif ! arrays of size 2 = (VIS, NIR)
           snow_refl_dir = s%snow_refl_dir ! arrays of size 2 = (VIS, NIR)

    if (is_watch_point()) then
        write(*,*) "snow  - gl_compute_snow_albedo:: computed albedo values:"
        write(*,*) "albedo_to_use = ", albedo_to_use
        write(*,*) "albedo_correction_to_use = ", albedo_correction_to_use
        write(*,*) snow_refl_dif
        write(*,*) snow_refl_dir
        write(*,*) s%beta_rad
    endif
    endif

end subroutine gl_compute_snow_albedo



subroutine compute_beta_rad_crocus(s, patm)
    class(snowpack_t), intent(inout) :: s !< state of snowpack
    real, intent(in) :: patm ! atm pressure [Pa] needed for effect of elevation on snow aging
    real :: beta_vis
    real :: beta_nir
    real :: beta_mir
    real :: beta_vis_sh, beta_nir_sh
    real :: rho ! snow density [kg m^-3]
    real :: age ! snow age [days]
    real :: optd ! snow layer optical diameter [m]
    real :: dprime, alphai, delta_alpha_age
    real :: P_over_PCDP

   ! use properties of top snow layer for now
!    rho = s%snow(1)%ws / s%snow(1)%dz
!    age = s%snow(1)%age
   age = s%nearsurf_age
   rho = s%nearsurf_rho
   optd = s%nearsurf_optd

   ! real :: dopt_dendritic, dopt_non_dendritic ! Not needed in Carmagnola 2013 Formulation
   ! real, intent(in) :: d  ! grain variables - snow dendricity [-]
   ! real, intent(in) :: s  ! grain variables - snow sphericity [-]
   ! real, intent(in) :: gs ! grain variables - snow grain size [m]

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! snow optical diameter:
   ! Pre -factor 0.5 added as in S, Morin et al., 2013, " Measurments and modeling of the
   ! vertical profile of specific surface area of an alpine snowpack"
   ! dopt_dendritic =  0.5 * 1E-4*(d + (1.0-d)*(4.0-s))
   ! dopt_non_dendritic = 0.5 * gs * s + (1-s)*max(4E-4, 0.5*gs)
   ! ! snow remains dendritic until d reaches zero (Vionnet et al., 2012):
   ! if (d > 1E-6) then
   !    dopt = dopt_dendritic
   ! else ! then, rounded crystals:
   !    dopt = dopt_non_dendritic
   ! endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! USE THE APPROACH BY CARMAGNOLA 2013 INSTEAD:
   ! IT IS BETTER TO USE DIRECTLY dopt AS A PROGNOSTIC VARIABLE IN THE SNOW MODEL
   ! snow albedo:
   P_over_PCDP = patm/87000.0 ! wrt 870 hPa. Need patm in [Pa] here
   delta_alpha_age = min(1.0, max(P_over_PCDP, 0.5))*0.2*age/60.0 ! age in [days]
   alphai = min(0.92, 0.96-1.58*sqrt(optd))
   dprime = min(optd, 0.0023)
!    snow_refl_vis = max(0.6, alphai-delta_alpha_age)
!    snow_refl_nir = max(0.3, 0.9-15.4*sqrt(optd))
!    snow_refl_mir = 346.3*dprime -32.31*sqrt(dprime) + 0.88

   ! penetration coefficient beta [m^-1] as in Vionnet et al., 2012
!    beta_vis = max(40.0, 0.00192*rho/sqrt(optd))
!    beta_nir = max(100.0, 0.01098*rho/sqrt(optd))
!    beta_mir = 1E6 ! (~infinite! - little to no light penetration in snowpack)
! Values from Brun 1992 [without the lower bound added in Vionnet et al., 2012]
   beta_vis = 0.00192*rho/sqrt(optd)
   beta_nir = 0.01098*rho/sqrt(optd)
   beta_mir = 1E4 ! (~infinite! - little to no light penetration in snowpack)

   ! assign to snowpack
   ! for now neglect third band from CROCUS ..
!    s%snow_refl_dir(BAND_VIS) = snow_refl_vis
!    s%snow_refl_dir(BAND_NIR) = (snow_refl_nir + snow_refl_mir)/2 ! todo: use correct weights
!    s%snow_refl_dir(BAND_NIR) = snow_refl_nir ! todo: use correct weights

!    s%snow_refl_dif(BAND_VIS) = snow_refl_vis
!    s%snow_refl_dif(BAND_NIR) = (snow_refl_nir + snow_refl_mir)/2 ! todo: use correct weights
!    s%snow_refl_dif(BAND_NIR) = snow_refl_nir  ! todo: use correct weights

    ! beta_vis = 1.0
    ! beta_nir = 1.0

    ! COMPUTE VALUES FOR GFDL BANDS BASED ON CROCUS BANDS
    ! CROCUS BANDS HAVE WEIGHTS [Vionnet et al., 2012]:
    ! wc1 = 0.71 ! [0.3 - 0.8]
    ! wc2 = 0.21 ! [0.8 - 1.5]
    ! wc3 = 0.08 ! [1.5 - 4.0]

!    s%beta_rad(BAND_VIS) = beta_vis
! !    s%beta_rad(BAND_NIR) = (beta_nir + beta_mir)/2 ! todo: use correct weights
!    s%beta_rad(BAND_NIR) = (0.21 * beta_nir + 0.08 * beta_mir)/(0.29)! todo: use correct weights

   ! INSTEAD USE JORNAN 1991 VALUES, AS IN SHRESTA ET AL 2006
   beta_vis_sh = 0.003759*rho/sqrt(optd)
   beta_nir_sh = 400.0
!    beta_vis_sh = 5.0
!    beta_nir_sh = 5.0
   s%beta_rad(BAND_VIS) = beta_vis_sh
   s%beta_rad(BAND_NIR) = beta_nir_sh

end subroutine compute_beta_rad_crocus

subroutine compute_albedo_crocus(s, patm)
    class(snowpack_t), intent(inout) :: s !< state of snowpack
    real, intent(in) :: patm ! atm pressure [Pa] needed for effect of elevation on snow aging
    real :: snow_refl_vis
    real :: snow_refl_nir
    real :: snow_refl_mir
    ! real :: beta_vis
    ! real :: beta_nir
    ! real :: beta_mir
    real :: rho ! snow density [kg m^-3]
    real :: age ! snow age [days]
    real :: optd ! snow layer optical diameter [m]
    real :: dprime, alphai, delta_alpha_age
    real :: P_over_PCDP

   ! use properties of top snow layer for now
!    rho = s%snow(1)%ws / s%snow(1)%dz
!    age = s%snow(1)%age
   age = s%nearsurf_age
   rho = s%nearsurf_rho
   optd = s%nearsurf_optd

   ! real :: dopt_dendritic, dopt_non_dendritic ! Not needed in Carmagnola 2013 Formulation
   ! real, intent(in) :: d  ! grain variables - snow dendricity [-]
   ! real, intent(in) :: s  ! grain variables - snow sphericity [-]
   ! real, intent(in) :: gs ! grain variables - snow grain size [m]

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! snow optical diameter:
   ! Pre -factor 0.5 added as in S, Morin et al., 2013, " Measurments and modeling of the
   ! vertical profile of specific surface area of an alpine snowpack"
   ! dopt_dendritic =  0.5 * 1E-4*(d + (1.0-d)*(4.0-s))
   ! dopt_non_dendritic = 0.5 * gs * s + (1-s)*max(4E-4, 0.5*gs)
   ! ! snow remains dendritic until d reaches zero (Vionnet et al., 2012):
   ! if (d > 1E-6) then
   !    dopt = dopt_dendritic
   ! else ! then, rounded crystals:
   !    dopt = dopt_non_dendritic
   ! endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! USE THE APPROACH BY CARMAGNOLA 2013 INSTEAD:
   ! IT IS BETTER TO USE DIRECTLY dopt AS A PROGNOSTIC VARIABLE IN THE SNOW MODEL
   ! snow albedo:
   P_over_PCDP = patm/87000.0 ! wrt 870 hPa. Need patm in [Pa] here
   delta_alpha_age = min(1.0, max(P_over_PCDP, 0.5))*0.2*age/60.0 ! age in [days]
   alphai = min(0.92, 0.96-1.58*sqrt(optd))
   dprime = min(optd, 0.0023)
   snow_refl_vis = max(0.6, alphai-delta_alpha_age)
   snow_refl_nir = max(0.3, 0.9-15.4*sqrt(optd))
   snow_refl_mir = 346.3*dprime -32.31*sqrt(dprime) + 0.88

   ! penetration coefficient beta [m^-1]
!    beta_vis = max(40.0, 0.00192*rho/sqrt(optd))
!    beta_nir = max(100.0, 0.01098*rho/sqrt(optd))
!    beta_mir = 1E6 ! (~infinite! - little to no light penetration in snowpack)

   ! assign to snowpack
   ! for now neglect third band from CROCUS ..
   s%snow_refl_dir(BAND_VIS) = snow_refl_vis
!    s%snow_refl_dir(BAND_NIR) = (snow_refl_nir + snow_refl_mir)/2 ! todo: use correct weights
!    s%snow_refl_dir(BAND_NIR) = snow_refl_nir ! todo: use correct weights
   s%snow_refl_dir(BAND_NIR) = (0.21 * snow_refl_nir + 0.08 * snow_refl_mir)/(0.29)

   s%snow_refl_dif(BAND_VIS) = snow_refl_vis
!    s%snow_refl_dif(BAND_NIR) = (snow_refl_nir + snow_refl_mir)/2 ! todo: use correct weights
!    s%snow_refl_dif(BAND_NIR) = snow_refl_nir  ! todo: use correct weights
   s%snow_refl_dif(BAND_NIR) = (0.21 * snow_refl_nir + 0.08 * snow_refl_mir)/(0.29)

!    s%beta_rad(BAND_VIS) = beta_vis
!    s%beta_rad(BAND_NIR) = (beta_nir + beta_mir)/2 ! todo: use correct weights
!    s%beta_rad(BAND_NIR) = beta_nir! todo: use correct weights

end subroutine compute_albedo_crocus




! compute snow properties needed to do soil-canopy-atmos energy balance
subroutine compute_albedo_lm4p2(s, snow_T, cosz, on_glacier)

   class(snowpack_t), intent(inout) :: s !< state of snowpack
!    real snow_T  ! snow temperature, deg K
   real, intent(in) :: cosz ! cosine of zenith angle
   real, intent(in) :: snow_T ! snow surf T [K]
   logical, intent(in) :: on_glacier ! TRUE if snow is on glacier
!    logical, intent(in) :: brdf_distinct_snow_on_glacier
   ! real, intent(out) :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_refl_lw, snow_emis
   real :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_refl_lw, snow_emis



   real :: f_iso_cold(NBANDS) = (/ 0.92, 0.58  /) ! VERONICA'S PAPER VALUES - VIS
   real :: f_vol_cold(NBANDS) = (/ 0.06, 0.08 /) ! VERONICA'S PAPER VALUES - VIS
   real :: f_geo_cold(NBANDS) = (/ 0.0,  0.0 /) ! VERONICA'S PAPER VALUES - VIS
   real :: f_iso_warm(NBANDS) = (/ 0.77, 0.43 /) ! VERONICA'S PAPER VALUES - VIS
   real :: f_vol_warm(NBANDS) = (/ 0.06, 0.08 /) ! VERONICA'S PAPER VALUES - VIS
   real :: f_geo_warm(NBANDS) = (/ 0.0,  0.0 /) ! VERONICA'S PAPER VALUES - VIS


            ! reflectance of snow on glaciers, otherwise snow reflectance does not depend
            ! on the underlying surface (except overlap).
   real :: f_iso_cold_on_glacier(NBANDS) = (/ 0.92, 0.73 /)
   real :: f_vol_cold_on_glacier(NBANDS) = (/ 0.06, 0.08 /)
   real :: f_geo_cold_on_glacier(NBANDS) = (/ 0.0, 0.0 /)
   real :: f_iso_warm_on_glacier(NBANDS) = (/ 0.77, 0.580 /)
   real :: f_vol_warm_on_glacier(NBANDS) = (/ 0.06, 0.08 /)
   real :: f_geo_warm_on_glacier(NBANDS) = (/ 0.0, 0.0 /)

   real    :: refl_snow_max_dir(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
   real    :: refl_snow_max_dif(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
   real    :: refl_snow_min_dir(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
   real    :: refl_snow_min_dif(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM

   real :: refl_snow_max_dir_on_glacier(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
   real :: refl_snow_max_dif_on_glacier(NBANDS) = (/ 0.8,  0.8  /) ! reset to 0.6 for MCM
   real :: refl_snow_min_dir_on_glacier(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
   real :: refl_snow_min_dif_on_glacier(NBANDS) = (/ 0.65, 0.65 /) ! reset to 0.45 for MCM
   ! real    :: emis_snow_max         = 0.95      ! reset to 1 for MCM
   ! real    :: emis_snow_min         = 0.90      ! reset to 1 for MCM
   ! ########### END PARAMETERS FOR THE ALBEDO MODEL CURRENTLY USED IN LM4P2 #############


!    if (s%nlayers == 0) then
!         ! write(*,*) "albedo lm4p2 warning: there is no snow on the ground! can't compute albedo"
!         ! snow_T = 273.15 - 10.0
!         call land_error_message("ERROR in compute_albedo_lm4p2 in snow_evolution module: There is no snow on the ground when routine was called!", FATAL)
!    else
!         ! snow_T = 273.15 - 10.0
!        snow_T = s%snow(1)%T
!    endif



if (on_glacier.and.distinct_snow_on_glacier) then
   call snow_rad_calculations_lm4p2 ( snow_T, cosz, &
      f_iso_warm_on_glacier, f_vol_warm_on_glacier, f_geo_warm_on_glacier, &
      f_iso_cold_on_glacier, f_vol_cold_on_glacier, f_geo_cold_on_glacier, &
      refl_snow_min_dir_on_glacier, refl_snow_max_dir_on_glacier, &
      refl_snow_min_dif_on_glacier, refl_snow_max_dif_on_glacier, &
      snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis )
else
   call snow_rad_calculations_lm4p2 ( snow_T, cosz, &
      f_iso_warm, f_vol_warm, f_geo_warm, &
      f_iso_cold, f_vol_cold, f_geo_cold, &
      refl_snow_min_dir, refl_snow_max_dir, &
      refl_snow_min_dif, refl_snow_max_dif, &
      snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis )


endif

      s%snow_refl_dif = snow_refl_dif
      s%snow_refl_dir = snow_refl_dir
   ! if need snow_refl_lw or snow_emis just add them fields the the snowpack structure
      ! update using T average ove thickness?
      ! add if needed the penetration length computed as in CROCUS?

end subroutine compute_albedo_lm4p2

! ============================================================================
subroutine snow_rad_calculations_lm4p2 ( snow_T, cosz, &
   f_iso_warm, f_vol_warm, f_geo_warm, &
   f_iso_cold, f_vol_cold, f_geo_cold, &
   refl_snow_min_dir, refl_snow_max_dir, &
   refl_snow_min_dif, refl_snow_max_dif, &
   snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis )
real, intent(in) :: snow_T  ! snow temperature, deg K
real, intent(in) :: cosz ! cosine of zenith angle
real, intent(in), dimension(NBANDS) :: &
   f_iso_warm, f_vol_warm, f_geo_warm, &
   f_iso_cold, f_vol_cold, f_geo_cold, &
   refl_snow_min_dir, refl_snow_max_dir, refl_snow_min_dif, refl_snow_max_dif
real, intent(out) :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_refl_lw, snow_emis

! ---- local vars
real :: blend
real :: warm_value_dir(NBANDS), cold_value_dir(NBANDS)
real :: warm_value_dif(NBANDS), cold_value_dif(NBANDS)
real :: zenith_angle, zsq, zcu

   logical, parameter :: use_brdf = .true. ! in lm4p2 this is set in nml, can change
!    logical, parameter :: use_brdf = .false. ! in lm4p2 this is set in nml, can change
   real, parameter :: t_range = 10.0 ! degK ! range of temperatures for ramp between "warm" and "cold" albedo

!    real    :: emis_snow_max         = 0.95      ! reset to 1 for MCM
!    real    :: emis_snow_min         = 0.90      ! reset to 1 for M
   real    :: emis_snow_max         = 1.0      ! reset to 1 for MCM
   real    :: emis_snow_min         = 1.0      ! reset to 1 for M

blend = max(0.,min(1.,1.-(tfreeze-snow_T)/t_range))
if (use_brdf) then
   zenith_angle = acos(cosz)
   zsq = zenith_angle*zenith_angle
   zcu = zenith_angle*zsq
   warm_value_dir = f_iso_warm*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                  + f_vol_warm*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                  + f_geo_warm*(g0_geo+g1_geo*zsq+g2_geo*zcu)
   cold_value_dir = f_iso_cold*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                  + f_vol_cold*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                  + f_geo_cold*(g0_geo+g1_geo*zsq+g2_geo*zcu)
   cold_value_dif = g_iso*f_iso_cold + g_vol*f_vol_cold + g_geo*f_geo_cold
   warm_value_dif = g_iso*f_iso_warm + g_vol*f_vol_warm + g_geo*f_geo_warm
else
   warm_value_dir = refl_snow_min_dir
   cold_value_dir = refl_snow_max_dir
   warm_value_dif = refl_snow_min_dif
   cold_value_dif = refl_snow_max_dif
endif
! write(*,*) "cold value dir = ", cold_value_dir
! write(*,*) "cold value dif = ", cold_value_dif
snow_refl_dir = cold_value_dir + blend*(warm_value_dir-cold_value_dir)
snow_refl_dif = cold_value_dif + blend*(warm_value_dif-cold_value_dif)
snow_emis     = emis_snow_max + blend*(emis_snow_min-emis_snow_max  )
snow_refl_lw  = 1 - snow_emis
! write(*,*) "snow_emis  = ", snow_emis
! write(*,*) "snow_refl_lw  = ", snow_refl_lw
end subroutine snow_rad_calculations_lm4p2


subroutine compute_albedo_he(s, cosz)
    ! ------------------------------------------------------------------------------------
    ! Snow albedo parameterization following "Impact of Grain Shape and Multiple
    ! Black Carbon Internal Mixing on Snow Albedo: Parameterization and Radiative
    ! Effect Analysis" By Celin He et al., 2018, JGR
    ! ------------------------------------------------------------------------------------

   class(snowpack_t), intent(inout) :: s !< state of snowpack
   real, intent(in) :: cosz ! cosine of solar zenith angle
   real rho_snow, snow_depth
   real snow_refl_vis_dir, snow_refl_nir_dir
   real snow_refl_vis_dif, snow_refl_nir_dif
   real :: grain_radius ! radius
   real :: ceqns, ceqns_im, ceqns_em ! black carbon equaivalent concentration in the snow surface layer [ppm]
   real :: Dcosz, phiVIS, phiNIR, effR_DcosVIS, effR_DcosNIR, effR
   real, PARAMETER :: aVIS = 0.781 ! FOR VIS BAND
   real, PARAMETER :: aNIR = 0.791 ! FOR NIR BAND
   real, PARAMETER :: aBROAD = 0.786 ! FOR BROADBAND
   real, PARAMETER :: b0 = 9.80508E-1 ! used by Veronica [Spherical snow, VIS band]
   real, PARAMETER :: b1 = -1.36104E-2
   real, PARAMETER :: b2 = -1.95416E-2
   real, PARAMETER :: d0 = 3.01470E-3
   real, PARAMETER :: d1 = 4.68312E-1
   real, PARAMETER :: d2 = 1.27961E-1

    integer idxshp
    real Re
    real delta_albedo_vis_im, delta_albedo_nir_im
    real delta_albedo_vis_em, delta_albedo_nir_em
    real delta_albedo_vis, delta_albedo_nir


    ! Select bands VIS = 0.3 - 0.7
    ! Select bands NIR = 0.7 - 4.0
    ! coeffs of the albedo parameterization
    real :: BB0(4,2) ! dimension (grain shape, band VIS=1, NIR=2)
    real :: BB1(4,2) ! dimension (grain shape, band VIS=1, NIR=2)
    real :: BB2(4,2) ! dimension (grain shape, band VIS=1, NIR=2)
    ! coeffs of the albedo reduction for internal mixing
    real :: D0_IM(4,2) ! dimension (grain shape, band VIS=1, NIR=2)
    real :: D1_IM(4,2) ! dimension (grain shape, band VIS=1, NIR=2)
    real :: D2_IM(4,2) ! dimension (grain shape, band VIS=1, NIR=2)
    ! coeffs of the albedo reduction for external mixing
    real :: D0_EM(4,2) ! dimension (grain shape, band VIS=1, NIR=2)
    real :: D1_EM(4,2) ! dimension (grain shape, band VIS=1, NIR=2)
    real :: D2_EM(4,2) ! dimension (grain shape, band VIS_IM=1, NIR=2)

    ! Albedo parameters
    data BB0(:,1) /9.80508E-01, 9.81658E-01, 9.84617E-01, 9.87015E-01  /  ! VIS
    data BB0(:,2) /6.44881E-01, 6.41406E-01, 6.66866E-01, 6.87695E-01  /  ! NIR

    data BB1(:,1) /-1.36104E-02, -1.56109E-02, -1.44812E-02, -1.34250E-02  /  ! VIS
    data BB1(:,2) /-1.84202E-01, -1.88446E-01, -1.77890E-01, -1.72401E-01  /  ! NIR

    data BB2(:,1) /-1.95416E-02, -1.93986E-02, -1.56859E-02, -1.33287E-02 /  ! VIS
    data BB2(:,2) /-2.64706E-02, -2.47061E-02, -2.75902E-02, -2.65446E-02  /  ! NIR



    ! In this case, hexagon values as the same as Koch as in He2018
    ! Internally mixed
    !                 SHPERE 1  - SPHEROID 2  - HEXAGONAL 3  - KOCH 4
    data D0_IM(:,1) /4.00225E-03, 4.00213E-03, 2.80587E-03, 2.80587E-03  /  ! VIS
    data D0_IM(:,2) /1.59399E-04, 1.48671E-04, 1.29313E-03, 1.29313E-03  /  ! NIR

    data D1_IM(:,1) /4.47263E-01, 4.57966E-01, 4.72359E-01, 4.72359E-01  /  ! VIS
    data D1_IM(:,2) /7.38319E-01, 7.46812E-01, 4.28793E-01, 4.28793E-01  /  ! NIR

    data D2_IM(:,1) /1.34864E-01, 1.25910E-01, 1.33016E-01, 1.33016E-01  /  ! VIS
    data D2_IM(:,2) /5.60124E-02, 5.15241E-02, 1.13120E-01, 1.13120E-01  /  ! NIR

    ! Externally mixed
    data D0_EM(:,1) /3.01470E-03, 2.39483E-03, 1.57041E-03, 1.57041E-03  /  ! VIS
    data D0_EM(:,2) /1.10568E-04, 7.33963E-05, 5.03177E-05, 5.03177E-05  /  ! NIR

    data D1_EM(:,1) /4.68312E-01, 4.82428E-01, 5.02587E-01, 5.02587E-01  /  ! VIS
    data D1_EM(:,2) /7.53732E-01, 7.78737E-01, 7.92972E-01, 7.92972E-01  /  ! NIR

    data D2_EM(:,1) /1.27961E-01, 1.25756E-01, 1.27417E-01, 1.27417E-01  /  ! VIS
    data D2_EM(:,2) /6.62209E-02, 6.55470E-02, 7.16926E-02, 7.16926E-02  /  ! NIR



    ! This parameterization holds for optically thick snow
    ! Here we use snow average properties average in the snowpack surface layer
    grain_radius = s%nearsurf_optd/2.0 !
    ! correct based on snow shape
    ! equivalent concentration of black carbon in surface layer [ppm]
    ! sum internally and externally mixed components

    ! ceqns = s%nearsurf_bceq_im + s%nearsurf_bceq_em
    ceqns_im = s%nearsurf_bceq_im
    ceqns_em = s%nearsurf_bceq_em
    rho_snow = s%nearsurf_rho
    snow_depth = s%depth()

    ! TODO: Add sphericity from near surface, and use to determine which parameterization
    ! must be use



    ! use near-surface snow properties to select which parameterization to use
    ! if (s%nearsurf_dendr > 0.5) then
    !     idxshp = 4 ! KOCH SNOWFLAKE [Case of fresh edndritic snow]
    ! else if (s%nearsurf_sph > 0.8) then
    !     idxshp = 1 ! SPHERE
    ! else if (s%nearsurf_sph < 0.2) then
    !     idxshp = 3 ! HEXAGONAL [Not spherical at all..]
    ! else ! remaining case of snow not very dendritic, and of intermediate sphericity
    !     idxshp = 2 ! SPHEROID
    ! endif

call compute_snow_grain_shape(s%nearsurf_dendr, s%nearsurf_sph, idxshp)




    ! idxshp = 1 ! do spherical snow only for now
    ! add external option to force shape?

    if (idxshp == 4) then
    ! CASE OF KOCH SNOWFLAKE
    ! NEED TO CORRECT EFFECTIVE RADIUS BECUASE THIS SHAPE IS NOT CONVEX
        Re = grain_radius / 0.544
    else ! NO NEED TO CORRECT EFF RADIUS FOR SPHERE, SPHERIOD OR HEXAGON
        Re = grain_radius
    endif


! correction of effective grain size based on Marshall (1989)
! and Wiscombe and Warren (1980) to account for cosz in snow albedo estimate
! should be done only for direct light [?]

   Dcosz = cosz - cos(49.5/180.0*PI) ! = cosz - 0.65
   phiVIS = (1.0 + aVIS *Dcosz)**2
   phiNIR = (1.0 + aNIR *Dcosz)**2
   ! effR = r * phi/100.0 ! /mu m
   effR = Re / 1E-4 ! No cosz correction for diffuse component
   ! note that the He parameterization was derived for cosz = 49.5 deg, which
   ! represents the insolation-weighted mean solar zenith cosine for sunlit Earth hemisphere
   effR_DcosVIS = Re * phiVIS / 1E-4 ! [den=100 \mu m, effective grain radius is in m]
   effR_DcosNIR = Re * phiNIR / 1E-4 ! [den=100 \mu m, effective grain radius is in m]
!    write(*,*) "HE EFF RADIUS = ", effR, phi, r, Dcosz, cosz, PI
!    snow_refl_vis = (b0 + b1*log(effR) + b2*log(effR**2)) - d0*(ceqns)**(d1*effR**d2)
   ! NOTE :: passing ceqns to PPB (see He at al., 2018)

! compute the albedo reduction due to impurities later on
!    snow_refl_vis = (b0 + b1*log10(effR) + b2*log10(effR)**2)
!    delta_albedo = d0*(ceqns*1E9)**(d1*effR**d2)

    ! direct flux albedo, VIS and NIR
    snow_refl_vis_dir = (BB0(idxshp, 1) + BB1(idxshp, 1)*log10(effR_DcosVIS) + BB2(idxshp, 1)*log10(effR_DcosVIS)**2)
    snow_refl_nir_dir = (BB0(idxshp, 2) + BB1(idxshp, 2)*log10(effR_DcosNIR) + BB2(idxshp, 2)*log10(effR_DcosNIR)**2)

    ! do not use cosz-based correction of grain size for diffuse rad -
    ! Thus, the snow effective grain radius is the same for VIS and NIR bands in this case
    snow_refl_vis_dif = (BB0(idxshp, 1) + BB1(idxshp, 1)*log10(effR) + BB2(idxshp, 1)*log10(effR)**2)
    snow_refl_nir_dif = (BB0(idxshp, 2) + BB1(idxshp, 2)*log10(effR) + BB2(idxshp, 2)*log10(effR)**2)

    ! write(*,*) "ceqnc_im [ppb] = ", ceqns_im
    ! write(*,*) "ceqnc_em [ppb] = ", ceqns_em

    ! pass concetration from [ppm] to [ppb]
    if (trim(lowercase(albedo_correction_to_use))=='he') then
        delta_albedo_vis_im = D0_IM(idxshp,1)*(ceqns_im*1E3)**(D1_IM(idxshp, 1)*effR**D2_IM(idxshp, 1))
        delta_albedo_nir_im = D0_IM(idxshp,2)*(ceqns_im*1E3)**(D1_IM(idxshp, 2)*effR**D2_IM(idxshp, 2))
        delta_albedo_vis_em = D0_EM(idxshp,1)*(ceqns_em*1E3)**(D1_EM(idxshp, 1)*effR**D2_EM(idxshp, 1))
        delta_albedo_nir_em = D0_EM(idxshp,2)*(ceqns_em*1E3)**(D1_EM(idxshp, 2)*effR**D2_EM(idxshp, 2))

    else if (trim(lowercase(albedo_correction_to_use))=='crocus') then
        delta_albedo_vis_em = min(0.2, 0.2 * s%nearsurf_age / 60.0)
        delta_albedo_nir_em = 0.0
        delta_albedo_vis_im = 0.0
        delta_albedo_nir_im = 0.0
    else if (trim(lowercase(albedo_correction_to_use))=='none') then

    ! delta_albedo_vis_em = D0_IM(idxshp,1)*((ceqns_im+ceqns_em)*1E3)**(D1_IM(idxshp, 1)*effR**D2_IM(idxshp, 1))
    ! delta_albedo_nir_em = D0_IM(idxshp,2)*((ceqns_im+ceqns_em)*1E3)**(D1_IM(idxshp, 2)*effR**D2_IM(idxshp, 2))
    delta_albedo_vis_em = 0.0
    delta_albedo_nir_em = 0.0
    delta_albedo_vis_im = 0.0
    delta_albedo_nir_im = 0.0

    else
        call land_error_message("Error in compute_albedo_he in snow_evolution_mod :: must specify a valid albedo_correction_to_use!", FATAL)
    endif

    ! apply albedo reduction due to impurities
    snow_refl_vis_dir = snow_refl_vis_dir - delta_albedo_vis_im - delta_albedo_vis_em
    snow_refl_vis_dif = snow_refl_vis_dif - delta_albedo_vis_im - delta_albedo_vis_em

    snow_refl_nir_dir = snow_refl_nir_dir - delta_albedo_nir_im - delta_albedo_nir_em
    snow_refl_nir_dif = snow_refl_nir_dif - delta_albedo_nir_im - delta_albedo_nir_em

    snow_refl_vis_dir = max(0.2, snow_refl_vis_dir )
    snow_refl_vis_dif = max(0.2, snow_refl_vis_dif )
    snow_refl_nir_dir = max(0.2, snow_refl_nir_dir )
    snow_refl_nir_dif = max(0.2, snow_refl_nir_dif )

    ! write(*,*) "albebo He :: effR , nearsurf optd= ", effR, s%nearsurf_optd
    ! write(*,*) "effR_DcosVIS, effR_DcosNIR", effR_DcosVIS, effR_DcosNIR
    s%snow_refl_dir = (/ snow_refl_vis_dir, snow_refl_nir_dir /)
    s%snow_refl_dif = (/ snow_refl_vis_dif, snow_refl_nir_dif /)
    ! write(*,*) "snow properties"
    ! call s%print()
    ! write(*,*) "snow surface properties"
    ! write(*,*) s%nearsurf_age
    ! write(*,*) s%nearsurf_optd
    ! write(*,*) s%nearsurf_rho
    ! write(*,*) s%nearsurf_sph
    ! write(*,*) s%nearsurf_bceq_em
    ! write(*,*) s%nearsurf_bceq_im
    ! write(*,*) s%nearsurf_bceq_tot

end subroutine compute_albedo_he


! subroutine compute_albedo_malinka(cosz, snow_refl_vis, D, snow_depth, ceqns, omega_bc, S_ext_bc, rho_snow, r)
! //FIXME deprecated
subroutine compute_albedo_malinka(s, cosz)


   class(snowpack_t), intent(inout) :: s !< state of snowpack
   real, intent(in) :: cosz
   real term1, term2, tau
   real :: snow_depth, rho_snow , grain_radius
   real x
   real r
   real S_ext_ice, S_ext_tot, omega_ice
   real :: ceqns ! SSA and extinction coeff of black carbon (from AM4???)

   real, PARAMETER :: g = 0.895 ! asymmetry parameter
   real, PARAMETER :: S_abs_ice = 5.795E-14 ! ice absorption cross section [m^2 kg^-1] ! approx!
   real, PARAMETER :: n = 1.313 ! real part of the refractive index ! approx!
   real, PARAMETER :: Tdiff = 0.9368 ! for visible wavelength! approx!

    grain_radius = s%nearsurf_optd / 2.0 ! METERS
    ! ceqns = s%nearsurf_bceq
    ceqns = s%nearsurf_bceq_im + s%nearsurf_bceq_em
    rho_snow = s%nearsurf_rho
    snow_depth = s%depth()
    ! D = 0.5 ! ratio diffuse to total light (instead, return both values!)

  ! //
   ! z = snow_depth
!    x = n**2 * r * S_abs_ice * rho_ice
   x = n**2 * grain_radius * S_abs_ice * rho_ice
   omega_ice = 1.0 - (x*Tdiff)/(x + Tdiff)
   S_ext_ice = S_abs_ice/(1-omega_ice)
   ! S_ext_tot = S_ext_ice + ceqns * S_ext_bc
   S_ext_tot = S_ext_ice + ceqns * LAI_ext(1) ! saved in constants
   tau = snow_depth * S_ext_tot * rho_snow

!    write(*,*) "tau Malinka = ", tau
   term1 = 1.0 - 12.0/7.0*(1.0 + 2.0*cosz)/(tau + 4.0)
   term2 = tau/(tau+4.0)
!    snow_refl_vis = (1.0-D)*term1 + D*term2

   ! model for VIS only
   ! USE HE always for NIR
    ! s%snow_refl_dir = (/ term1, -9999.0 /)
    ! s%snow_refl_dif = (/ term2, -9999.0 /)
   ! // FIXME do not pass the same for NIR if it not computed here ..
    s%snow_refl_dir = (/ term1, term1 /)
    s%snow_refl_dif = (/ term2, term2 /)

end subroutine compute_albedo_malinka


subroutine compute_albedo_rozenberg(s, cosz)
   ! # D = Fraction of diffuse to total SW radiation
   ! //TODO: separate albedo for DIR and DIF components
    ! //FIXME deprecated

   class(snowpack_t), intent(inout) :: s !< state of snowpack
   real, intent(in) :: cosz
   real term1, term2, tau
   real snow_depth, rho_snow
   real gamma ! asymptotic attenuation coefficient
   real ceqns ! conc of black carbon equivalent impurities
   real grain_radius
   real y !
   real omega0 ! single scattering albedo
   real x
   real r ! grain size radius?
   real omega_ice, S_ext_ice, S_ext_tot
   ! real, intent(in) :: omega_bc, S_ext_bc ! SSA and extinction coeff of black carbon (from AM4???)
   ! real, intent(in) :: r ! snow grain radius

   real, PARAMETER :: g = 0.895 ! asymmetry parameter
   real, PARAMETER :: S_abs_ice = 5.795E-14 ! ice absorption cross section [m^2 kg^-1] ! approx!
   real, PARAMETER :: n = 1.313 ! real part of the refractive index ! approx!
   real, PARAMETER :: Tdiff = 0.9368 ! for visible wavelength! approx!

    grain_radius = s%nearsurf_optd/2.0 ! METERS
    ! ceqns = s%nearsurf_bceq
    ceqns = s%nearsurf_bceq_im + s%nearsurf_bceq_em
    rho_snow = s%nearsurf_rho
    snow_depth = s%depth()

   x = n**2 * grain_radius * S_abs_ice * rho_ice
!    x = n**2 * r * S_abs_ice * rho_snow
!    omega_ice = 1.0 - (x*Tdiff)/(x + Tdiff)
   write(*,*) "*", x*Tdiff
   write(*,*) "+", x+Tdiff
   omega_ice = 1.0 - (x*Tdiff)/(x + Tdiff)
   S_ext_ice = S_abs_ice/(1-omega_ice)
   ! S_ext_tot = S_ext_ice + ceqns * S_ext_bc
   S_ext_tot = S_ext_ice + ceqns * LAI_ext(1)
   tau = snow_depth * S_ext_tot * rho_snow
   ! omega0 = (S_ext_ice * omega_ice + S_ext_bc * ceqns * omega_bc)/S_ext_tot
   omega0 = (S_ext_ice * omega_ice + ceqns * LAI_ext(1) * LAI_ssa(1))/S_ext_tot
   gamma = sqrt(3.0*(1-omega0)*(1-omega0*g))
   y = 4.0*sqrt( (1.0-omega0)/(3.0*(1.0-omega0*g)))
!    write(*,*) "Rozenberg : omega_ice, x, Tdiff, n, Sabsice, rhoice", omega_ice, x, Tdiff, n, S_abs_ice, rho_ice
!    write(*,*) "Rozenberg : snowe_depth, rho, S_ext_tot, S_ext_ice", snow_depth, rho_snow, S_ext_tot, S_ext_ice
!    write(*,*) "Rozenberg : gamma, tau, y, y2", gamma, tau,y,  y*(1.0-3.0/7.0*(1.0+2.0*cosz))
   term1 = sinh( gamma*tau + y*(1.0-3.0/7.0*(1.0+2.0*cosz))    ) / sinh(gamma*tau + y)
   term2 = sinh(gamma*tau)/sinh(gamma*tau+ y)


!    write(*,*) "Rozenberg : den", gamma*tau + y
!    write(*,*) "Rozenberg : den", sinh(gamma*tau )
!    write(*,*) "Rozenberg : term1, term2", term1, term2
!    snow_refl_vis = (1-D)*term1 + D*term2

   ! // FIXME do not pass the same for NIR if it not computed here ..
    s%snow_refl_dir = (/ term1, term1 /)
    s%snow_refl_dif = (/ term2, term1 /)
end subroutine compute_albedo_rozenberg


!> \Brief Snow step 2
! Perform Water and Energy balance for snowpack and update snow properties
subroutine gl_snow_step_2 ( s, snow_subl,                     &
                        vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
                        DTg,  Mg_imp,  evapg,  fswg,  flwg,  sensg,  &
                        use_tfreeze_in_grnd_latent, &
                        ! output
                        subs_DT, &
                        subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens,  &
                        snow_fsw, snow_flw, snow_sens, &
                        snow_levap, snow_fevap, snow_melt, &
                        snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
                        snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, &
                        snow_avrg_T , &
                        ! additional input/output added by Enrico for standalone model only
                        !    snow_rho, snow_age, snow_sph, snow_optd, & ! average snow properties
                        ! heat1, verbose, hfevap, dt, wind_atm, t_atm, &
                        dt, wind_atm, t_atm, p_surf, &
                        wetdep, drydep, grnd_T_preprec, &
                        ! for conservation checks only :
                        begw_check, begh_check, &
                        G0, DGDTg, snow_G_Z, snow_G_TZ, &
                        mass_lai_em_1, mass_lai_im_1, &
                        lost_wc_em_st, lost_wc_im_st, &
                        lost_wc_em, lost_wc_im)
                        ! delta_heat_DTg)
    type(snowpack_t), intent(inout) :: s ! snowpack = instance of snow tile object
    type(snowpack_t) :: s_check ! for debug only
    real, intent(in) :: snow_subl ! fraction of sublimation (equal to 1 when snow is there)
    real, intent(in) :: vegn_lprec ! precip below canopy [Kg m^-2 s^-1]
    real, intent(in) :: vegn_fprec ! precip below canopy [Kg m^-2 s^-1]
    real, intent(in) :: vegn_hlprec ! heat carried by precip [J m^-2 s^-1]
    real, intent(in) :: vegn_hfprec ! heat carried by precip [J m^-2 s^-1]
    real, intent(in), DIMENSION(NTRACERS) :: mass_lai_em_1, mass_lai_im_1 ! mass of LAIs at beginning of step, for mass cons checks
    real, intent(in), DIMENSION(NTRACERS) :: lost_wc_em_st, lost_wc_im_st ! from sweep tiny snow, to check mass balance
    real, intent(out), DIMENSION(NTRACERS) :: lost_wc_em, lost_wc_im
    real, intent(in) :: DTg, Mg_imp, evapg, fswg, flwg, sensg
    logical, intent(in) :: use_tfreeze_in_grnd_latent
    real, intent(out) :: &
         subs_DT, subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens, &
         snow_fsw, snow_flw, snow_sens, &
         snow_levap, snow_fevap, snow_melt, &
         snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
         snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, snow_avrg_T
    ! added in ezsnow
    real, intent(out) :: grnd_T_preprec
    real, intent(in) :: dt ! delta time step
    real, intent(in) :: wind_atm, t_atm
    real, intent(in) :: p_surf ! surface atm pressure in [Pa]
    real, intent(in) :: wetdep(NTRACERS) ! wet deposition of tracers from atmosphere [ppm]
    real, intent(in) :: drydep(NTRACERS) ! dry deposition of tracers from atmosphere [mg m^-2 s^-1]
    real, intent(in) :: begw_check, begh_check
    real, intent(in) :: G0, DGDTg, snow_G_Z, snow_G_TZ
    !  local variables
    real hfevap
    logical :: verbose
    real check_heat0, check_heat1, ftprec, ltprec
    real snow_lprec1, snow_hlprec1, snow_lprec2, snow_hlprec2
    real heat1a, heat1b, heat1c, heat1d, heat1e, heat1f, heat1g
    real netmass1, netmass2, netmassdiff
    real netheat1, netheat2, netheatdiff
    real lswept1, fswept1, hlswept1, hfswept1 ! from sublimation
    real lswept2, fswept2, hlswept2, hfswept2
    integer il
    real, dimension(NTRACERS) :: lost_wc_em1, lost_wc_im1, lost_wc_em2, lost_wc_im2  ! [mg/m2]
    real, dimension(NTRACERS) :: lost_wc_em3, lost_wc_im3, lost_wc_em4, lost_wc_im4, lost_wc_em5, lost_wc_im5 ! [mg/m2]
    real laimass1 , laimass2 , netlaimass ! [mg/m2]
    real laimass_wetdep_rainf, laimass_wetdep_snowf
    real dheat_fevap
    real hlevap
    real vegn_hlprec_r ! lm4p2 conv, used only for energy balance checks
    real vegn_hfprec_ch ! used only for energy balance checks
    real endh_check, neth_check, endw_check, netw_check
    real, dimension(NTRACERS) :: total_wetdep_check, net_delta_lai, mass_lai_em_2, mass_lai_im_2
    real sum_swheat
    real delta_time
    logical thick_enough_for_evap
    real frunf_from_deficit, hfrunf_from_deficit, frac_of_deficit
    real total_depth

    if(is_watch_point()) then
        write(*,*)'###### Beginning GLASS Snow step 2 ######'
        write(*,*) "vegn_lprec, vegn_hlprec = ", vegn_lprec, vegn_hlprec
    endif

    delta_time = dt
    verbose = .False.

    call s%update_age(dt) ! update age of existing snow layers

    heat1a = s%heat()
    if(verbose) write(*,*) "STEP2: heat check A = ", heat1a

    if (s%nlayers>0) then
        snow_fsw   = fswg
        snow_flw   = flwg
        snow_sens  = sensg
        snow_levap = evapg *(1-snow_subl)
        snow_fevap = evapg *   snow_subl
        subs_fsw = fswg - snow_fsw ! EZSNOW
        subs_flw = flwg - snow_flw
        subs_evap = evapg - snow_levap - snow_fevap
        subs_sens = sensg - snow_sens
    else
        snow_fsw    = 0.0
        snow_flw    = 0.0
        snow_sens   = 0.0
        snow_levap  = 0.0
        snow_fevap  = 0.0
        subs_fsw = fswg
        subs_flw = flwg
        subs_evap = evapg
        subs_sens = sensg
    endif

    ! ! EZSNOW: If snow too thin, sublimate from substrate insetad
    ! if(s%depth() < min_snow_depth) then
    !     snow_levap = 0.0
    !     snow_fevap = 0.0
    !     subs_evap = evapg
    !     ! subs_evap = 0.0
    !     thick_enough_for_evap = .FALSE.
    ! else
    !     thick_enough_for_evap = .TRUE.
    ! endif

!   if(is_watch_point()) then
!      write(*,*) '#### gl_snow_step_2 ### checkpoint 1 ####'
!      __DEBUG3__(thick_enough_for_evap, subs_evap, evapg)
!   endif

    snow_lrunf = 0 ! init water and heat runoff terms
    snow_frunf = 0
    snow_hlrunf = 0
    snow_hfrunf = 0

    netheat1 = s%heat() ! heat cons check

    ! check this only outside
    ! update snowpack temperature vertical profile using coeffs computed in step 1 (e,f)
    call s%step2(DTg, subs_DT) ! return same dT if there is no snow
    ! heat1 = s%heat() ! energy balance from here?
    ! delta_heat_DTg = heat1 - heat1a ! change in snowpack energy due to heat conduction
    ! write(*,*) "STEP-2-HEAT CHECKPOINT 2 [After update T profile]:", s%heat()

        if (s%nlayers>0) then
            ! write(*,*) "T1 after updating T progile = ", s%snow(1)%T
            s%preprec_surfT = s%snow(1)%T ! get surface T before subl and snowfall
            ! grnd_T_preprec = s%snow(1)%T ! get surface T before subl and snowfall
        else
            s%preprec_surfT = -1.0 ! start with fill value in case there is no snow
            ! grnd_T_preprec = -1.0 ! start with fill value in case there is no snow
        endif

    netheat2 = s%heat() ! heat cons check
    netheatdiff = netheat2 - netheat1
    ! if (snow_check_cons .and. abs(netheatdiff )>1E-6) then
    ! write(*,*) "start ----- DE due to update temperature profile: ........."
    !     write(*,*) "snow nlayers = ", s%nlayers
    !     write(*,*) "initial heat = ", netheat1
    !     write(*,*) "final heat = ", netheat2
    !     write(*,*) "heat net difference = ", netheatdiff
    !     write(*,*) "DTg, subs_DT", DTg, subs_DT
    !     ! error stop "heat balance violation after snow step 2: Updating temperature profile!"
    !     ! write(*,*) "end ---- DE due to update temperature profile: ........."
    ! endif

    ! check that some snow related quantities are within physical boundaries
    call s%check_bounds("check bounds before sublimation ....")


    ! if(verbose) write(*,*) "SNOW STEP 2 : Sweep tiny snow"
    ! call sweep_tiny_snow(s,lswept2, fswept2, hlswept2, hfswept2, lost_wc_em5, lost_wc_im5)
    ! ! lswept2 = 0; fswept2 = 0; hlswept2 = 0; hfswept2 = 0
    ! snow_lrunf = snow_lrunf + lswept2/dt
    ! snow_frunf = snow_frunf + fswept2/dt
    ! snow_hlrunf = snow_hlrunf + hlswept2/dt
    ! snow_hfrunf = snow_hfrunf + hfswept2/dt

        !    if(is_watch_point()) then
        !       write(*,*)'#### Snow step 2 : before snow sublimation ####'
        !       call s%print()
        !       write(*,*) "evapg, Mg_imp = ", evapg, Mg_imp
        !   endif

    ! write(*,*) "Snow step 2: Snowf = ", vegn_fprec


    netmass1 = s%SWE() ! mass cons check
    netheat1 = s%heat() ! heat cons check
    laimass1 = sum(s%lai_em() + s%lai_im()) ! sum across NTRACERS dimension for purposes of mass cons check

    heat1b = s%heat()
    if(verbose) write(*,*) "STEP2: heat check B = ", heat1b
    if(verbose) write(*,*) "STEP2: heat check B - A = ", heat1b - heat1a

    if(verbose) write(*,*) "SNOW STEP 2 : Do snow sublimation"
    call snow_sublimation(s, dt, snow_levap, snow_fevap, hfevap, hlevap, dheat_fevap, &
                use_tfreeze_in_grnd_latent, DTg, Mg_imp, snow_melt, &
                lswept1, fswept1, hlswept1, hfswept1, &
                subs_m_imp, lost_wc_em1, lost_wc_im1, thick_enough_for_evap, verbose=.FALSE.)


    if(is_watch_point()) then
        write(*,*)'#### Snow step 2 : after snow sublimation ####'
        call s%print()
        write(*,*) "snow_levap, snow_fevap, hfevap, hlevap, snow_melt, subs_m_imp",snow_levap, snow_fevap, hfevap, hlevap, snow_melt, subs_m_imp
    endif


    netheat2 = s%heat() ! heat cons check
    netheatdiff = netheat2 - netheat1 - (Mg_imp-subs_M_imp)*HLF + hfevap*dt - dheat_fevap + hlswept1 + hfswept1
    if (do_snow_check_cons .and. (abs(netheatdiff )>1E-2)) then
        write(*,*) "DEN = ", - Mg_imp*HLF + hfevap*dt - dheat_fevap
        write(*,*) "snow nlayers = ", s%nlayers
        write(*,*) "initial heat = ", netheat1
        write(*,*) "final heat = ", netheat2
        write(*,*) "heat net difference = ", netheatdiff
        write(*,*) "heat due to Explicit [percomputed in SEB] melt/freeze:", Mg_imp*HLF
        write(*,*) "fevap*dt*HLF, levap*dt*HLF = ",snow_fevap*dt*HLF, snow_levap*dt*HLF
        write(*,*) "hfevap*dt = ",hfevap*dt
        write(*,*) "hlevap*dt = ",hlevap*dt
        write(*,*) "dheat_fevap = ",dheat_fevap
        write(*,*) "Mg_imp = ",Mg_imp
        write(*,*) "Mg_imp*HLF = ",Mg_imp*HLF
        write(*,*) "snow_melt = ",snow_melt
        write(*,*) "subs_M_imp = ",subs_M_imp
        write(*,*) "snow_melt * HLF = ",snow_melt*HLF
        call land_error_message("ERROR gl_snow_step_2 in snow_evolution module: heat balance violation after snow_sublimation!", FATAL)
    endif
    ! delta_heat_DTg = dheat_fevap ! export this quantity for global energy conservation checks


    laimass2 = sum(s%lai_em() + s%lai_im())
    netlaimass = laimass2 - laimass1 + sum(lost_wc_em1) + sum(lost_wc_im1)
    if (do_snow_check_cons .and.(abs(netlaimass)>1E-6)) then
        write(*,*) "checkpoint after sublimation: net lai mass difference = ", netlaimass
        write(*,*) "checkpoint after sublimation: initial LAI content = ", laimass1
        write(*,*) "checkpoint after sublimation: final LAI content = ", laimass2
        write(*,*) "checkpoint after sublimation: lost LAI content = ",sum(lost_wc_em1) + sum(lost_wc_im1)
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: LAI balance violation after snow_sublimation!", FATAL)
    endif

    netmass2 = s%SWE() + (snow_levap + snow_fevap)*dt  +lswept1 + fswept1 ! mass cons check
    if (do_snow_check_cons .and.(abs(netmass2 - netmass1)>1E-6)) then
        write(*,*) "checkpoint after sublimation: SWE = ", s%SWE()
        call land_error_message("ERROR gl_snow_step_2 in snow_evolution module: mass balance violation after snow_sublimation!", FATAL)
    endif

    if (s%nlayers > 0) then
        if (s%snow(1)%T < 100.0) then
            call s%print()
            write(*,*) "After sublimation T(1) = ", s%snow(1)%T
            call land_error_message("ERROR gl_snow_step_2 in snow_evolution module: MIN TEMPERATURE violation after snow_sublimation!", FATAL)
        endif
    endif


    call s%check_bounds("check bounds after sublimation ....")


    if(verbose) write(*,*) "SNOW STEP 2 : Do snow melt and freeze"
    if(is_watch_point()) then
        write(*,*)'#### Snow step 2 : before snow melt and freeze ####'
        call s%print()
    endif
    netmass1 = s%SWE() ! init mass cons check
    netheat1 = s%heat()  ! init heat cons check
    laimass1 = sum(s%lai_em() + s%lai_im()) ! sum across NTRACERS dimension for purposes of LAI mass cons check
    call snow_melt_and_freeze( &
            s, dt, snow_lprec1, &
            snow_hlprec1, lost_wc_em2, lost_wc_im2,  verbose_in=.FALSE.)
    netmass2 = s%SWE() + (snow_lprec1)*dt ! mass cons check
    if (do_snow_check_cons .and.(abs(netmass2 - netmass1)>1E-6)) then
        write(*,*) "netmass1 = ", netmass1
        write(*,*) "netmass2 = ", netmass2
        write(*,*) "SWE2, snow_lprec1*dt = ", s%SWE(), snow_lprec1*dt
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: mass balance violation after melt_and_freeze!", FATAL)
    endif
    netheat2 = s%heat()
    netheatdiff = netheat2 - (netheat1 - snow_hlprec1*dt )
    if (do_snow_check_cons .and.(abs(netheatdiff )>1E-2)) then
        write(*,*) "snow nlayers = ", s%nlayers
        write(*,*) "initial heat = ", netheat1
        write(*,*) "final heat = ", netheat2
        write(*,*) "heat net difference = ", netheatdiff
        write(*,*) "sink percolation heat = ", snow_hlprec1*dt
        write(*,*) "sink percolation, only HLF = ", snow_lprec1*dt*HLF
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: heat balance violation after snow_melt_freeze!", FATAL)
    endif
    laimass2 = sum(s%lai_em() + s%lai_im())
    netlaimass = laimass2 - laimass1 + sum(lost_wc_em2) + sum(lost_wc_im2)
    if (do_snow_check_cons .and.(abs(netlaimass)>1E-6)) then
        call s%print()
        write(*,*) "checkpoint after snow_melt_and_freeze: net lai mass difference = ", netlaimass
        write(*,*) "checkpoint after snow_melt_and_freeze: initial LAI content = ", laimass1
        write(*,*) "checkpoint after snow_melt_and_freeze: final LAI content = ", laimass2
        write(*,*) "checkpoint after snow_melt_and_freeze: lost LAI content = ",sum(lost_wc_em2) + sum(lost_wc_im2)
        write(*,*) "current snow nlayers = ", s%nlayers
        write(*,*) "water before: netmass1 = ", netmass1
        write(*,*) "water after: netmass2 = ", netmass2
        write(*,*) "current SWE, water lost: snow_lprec1*dt = ", s%SWE(), snow_lprec1*dt
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: LAI balance violation after snow_melt_and_freeze!", FATAL)
    endif
    if(is_watch_point()) then
        write(*,*)'#### Snow step 2 : after snow melt and freeze ####'
        call s%print()
    endif
    call s%check_bounds("check bounds after melt and freeze ....")


    netmass1 = s%SWE() ! mass cons check
    netheat1 = s%heat() ! heat cons check

    ! //FIXME added 2nd relayering step here
    if(verbose) write(*,*) "SNOW STEP 2 : Do snowpack relayering"
    if (s%nlayers > 0) then
        ! write(*,*) "numbers of snow layers before relayering = ", s%nlayers
        if (do_merge) call s%attempt_merge_layers()
        ! write(*,*) "numbers of snow layers during relayering (before split, after merge) = ", s%nlayers
        if (do_split) call s%attempt_split_layers()
        ! write(*,*) "numbers of snow layers after relayering = ", s%nlayers
    endif

    netmass2 = s%SWE() ! mass cons check
    if (do_snow_check_cons .and.(abs(netmass2 - netmass1)>1E-6)) then
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: mass balance violation after snowpack relayering!", FATAL)
    endif
    netheat2 = s%heat()
    if (do_snow_check_cons .and.(abs(netheat2 - netheat1)>1E-2)) then
        write(*,*) "before: heat = ", netheat1
        write(*,*) "difference: dheat = ", netheat2 - netheat1
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: heat balance violation during snowpack relayering!", FATAL)
    endif


    call s%check_bounds("check bounds after relayering ....")



    ! SNOW SOLID BALANCE
    netmass1 = s%SWE() ! mass cons check
    netheat1 = s%heat() ! heat cons check
    laimass1 = sum(s%lai_em() + s%lai_im())
    if (abs(vegn_fprec)>0.0) then
        ftprec = TFREEZE + vegn_hfprec/CSW/vegn_fprec ! get T from heat content - LM4p2 does not include latent
    else
        ftprec = 273.15
    endif
    if(verbose) write(*,*) "SNOW STEP 2 : Do snow solid balance"
    call snow_solid_balance(s, vegn_fprec, snow_fevap, vegn_lprec, snow_levap, ftprec, &
                            wetdep, drydep, wind_atm, t_atm, lost_wc_em3, lost_wc_im3, dt, verbose_in=.FALSE.)
                            ! //TODO: remove evap from here
    call s%check_bounds("check bounds after solid balance ....")

    laimass_wetdep_rainf = 0.0
    laimass_wetdep_snowf = 0.0
    if (vegn_fprec > 1E-9) then
        laimass_wetdep_snowf = sum(wetdep)*vegn_fprec
    endif
    if (vegn_lprec > 1E-9) then
        laimass_wetdep_rainf = sum(wetdep)*vegn_lprec
    endif
    laimass2 = sum(s%lai_em() + s%lai_im())
    netlaimass = laimass2 - laimass1 + sum(lost_wc_em3) + sum(lost_wc_im3) - laimass_wetdep_snowf*dt - sum(drydep)*dt
    if (do_snow_check_cons .and.(abs(netlaimass )>1E-6)) then
        write(*,*) "checkpoint after snow_solid_balance: net lai mass difference = ", netlaimass
        write(*,*) "checkpoint after snow_solid_balance: initial LAI content = ", laimass1
        write(*,*) "checkpoint after snow_solid_balance: final LAI content = ", laimass2
        write(*,*) "checkpoint after snow_solid_balance: total dry deposition = ",sum(drydep)*dt
        write(*,*) "checkpoint after snow_solid_balance: total wet deposition due to snowfall = ",laimass_wetdep_snowf*dt
        write(*,*) "checkpoint after snow_solid_balance: lost LAI content = ",sum(lost_wc_em3) + sum(lost_wc_im3)
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: LAI balance violation after snow_solid_balance!", FATAL)
    endif
    netmass2 = s%SWE() - (vegn_fprec)*dt ! mass cons check
    netheat2 = s%heat()
    if (do_snow_check_cons .and.(abs(netmass2 - netmass1)>1E-6)) then
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: mass balance violation after snow_solid_balance!", FATAL)
    endif
    ! added if to make sure that model is not fed vegn_fprec = 0 with vegn_hfprec != 0
    if (vegn_fprec>0.) then
        vegn_hfprec_ch = vegn_hfprec
        ! netheatdiff = netheat2 - (netheat1 + vegn_hfprec*dt)
    else
        vegn_hfprec_ch = 0.
        ! netheatdiff = netheat2 - (netheat1)
    endif
    netheatdiff = netheat2 - (netheat1 + vegn_hfprec_ch*dt)
    if (do_snow_check_cons .and.(abs(netheatdiff )>1E-2)) then
        write(*,*) "snow nlayers = ", s%nlayers
        write(*,*) "initial heat = ", netheat1
        write(*,*) "final heat = ", netheat2
        write(*,*) "vegn_hfprec*dt = ", vegn_hfprec_ch*dt
        write(*,*) "[temp. of frozen precip]ftprec = ", ftprec
        write(*,*) "[frozen precip added] vegn_fprec * dt = ", vegn_fprec*dt
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: heat balance violation after snow_solid_balance!", FATAL)
    endif
    ! END SNOW SOLID BALANCE

    if (abs(vegn_lprec)>0.0) then
        ! ltprec = TFREEZE + vegn_hlprec/CLW/abs(vegn_lprec) ! get T from heat content - LM4p2 does not include latent
        ltprec = TFREEZE + vegn_hlprec/CLW/vegn_lprec ! get T from heat content - LM4p2 does not include latent
    else
        ltprec = 273.15  ! no energy added, but needs to be defined
    endif
    if(is_watch_point()) then
        write(*,*)'#### Snow step 2 : before snow liquid balance ####'
        write(*,*) "vegn_lprec, vegn_hlprec, ltprec = ", vegn_lprec, vegn_hlprec, ltprec
    endif
    call s%check_bounds("check bounds before liquid balance ....")
    netmass1 = s%SWE() ! mass cons check
    netheat1 = s%heat() ! heat cons check
    laimass1 = sum(s%lai_em() + s%lai_im())  ! LAIs cons check
    if(verbose) write(*,*) "SNOW STEP 2 : Do snow liquid water balance"
    call snow_liquid_balance(s, vegn_lprec, snow_levap, &
            vegn_fprec, ltprec, wetdep, snow_lprec2, snow_hlprec2, &
            lost_wc_em4, lost_wc_im4, dt, verbose_in=.FALSE.)

    laimass2 = sum(s%lai_em() + s%lai_im())
    netlaimass = laimass2 - laimass1 + sum(lost_wc_em4) + sum(lost_wc_im4) - laimass_wetdep_rainf*dt
    if (do_snow_check_cons .and.(abs(netlaimass )>1E-6)) then
        write(*,*) "checkpoint after snow_liquid_balance: net lai mass difference = ", netlaimass
        write(*,*) "checkpoint after snow_liquid_balance: initial LAI content = ", laimass1
        write(*,*) "checkpoint after snow_liquid_balance: final LAI content = ", laimass2
        write(*,*) "checkpoint after snow_liquid_balance: total wet deposition due to rainfall = ",laimass_wetdep_rainf*dt
        write(*,*) "checkpoint after snow_liquid_balance: lost LAI content = ",sum(lost_wc_em4) + sum(lost_wc_im4)
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: LAI balance violation after snow_liquid_balance!", FATAL)
    endif

    netmass2 = s%SWE() - (vegn_lprec)*dt + (snow_lprec2)*dt ! mass cons check
    if (do_snow_check_cons .and.(abs(netmass2 - netmass1)>1E-6)) then
        write(*,*) "checkpoint after snow liquid balance: SWE = ", s%SWE()
        write(*,*) "before: SWE = ", netmass1
        write(*,*) "after: SWE, vegn_lprec*dt, snow_lprec*dt",  s%SWE(), (vegn_lprec)*dt, (snow_lprec2)*dt
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: mass balance violation after snow_liquid_balance!", FATAL)
    endif

    ! note input heat from lm4p2 is wrt liquid at TF
    netheat2 = s%heat()
    ! netheatdiff = netheat2 - (netheat1 + vegn_hlprec*dt - snow_hlprec2*dt) ! ENEROLD
    netheatdiff = netheat2 - (netheat1 + vegn_hlprec*dt + vegn_lprec*HLF*dt - snow_hlprec2*dt)  ! vegn_hlprec does not include HLF
    ! if (do_snow_check_cons .and.(abs(netheatdiff )>1E-6)) then
    if (do_snow_check_cons .and.(abs(netheatdiff )>1E-2)) then
        write(*,*) "snow nlayers = ", s%nlayers
        write(*,*) "initial heat = ", netheat1
        write(*,*) "final heat = ", netheat2
        write(*,*) "added precip heat = ", vegn_hlprec*dt
        write(*,*) "sink percolation heat = ", snow_hlprec*dt
        write(*,*) "net heat difference: netheatdiff  = ",netheatdiff
        write(*,*) "END - START  = ",netheat2 - (netheat1 )
        write(*,*) "IN - OUT  = ",vegn_hlprec*dt - snow_hlprec2*dt
        write(*,*) "latent heat in / out: vegn_lprec*dt*HLF , snow_hlprec2*dt*HLF", vegn_lprec*dt*HLF, snow_lprec2*dt*HLF
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: heat balance violation after snow_liquid_balance!", FATAL)
    endif

    if(is_watch_point()) then
        write(*,*)'#### Snow step 2 : after snow liquid balance ####'
        write(*,*) "vegn_lprec, vegn_hlprec, ltprec = ", vegn_lprec, vegn_hlprec, ltprec
        write(*,*) "snow_lprec * HLF, snow_hlprec1, snow_hlprec2, snow_hlprec1+snow_hlprec2 = ",snow_lprec * HLF, snow_hlprec1, snow_hlprec2, snow_hlprec1+snow_hlprec2 ! should include HLF here
    !   call s%print()
    endif


    netmass1 = s%SWE() ! mass cons check
    if(verbose) write(*,*) "SNOW STEP 2 : Sweep tiny snow"
    ! FIXME: removed sweep tiny snow
    call gl_sweep_tiny_snow(s,lswept2, fswept2, hlswept2, hfswept2, lost_wc_em5, lost_wc_im5)
    ! lswept2 = 0; fswept2 = 0; hlswept2 = 0; hfswept2 = 0
    ! lost_wc_em5=0; lost_wc_im5=0
    snow_lrunf = snow_lrunf + lswept1/dt +  lswept2/dt
    snow_frunf = snow_frunf + fswept1/dt + fswept2/dt
    snow_hlrunf = snow_hlrunf + hlswept1/dt + hlswept2/dt
    snow_hfrunf = snow_hfrunf + hfswept1/dt + hfswept2/dt
    netmass2 = s%SWE() + lswept2 + fswept2 ! mass cons check
    if (do_snow_check_cons .and.(abs(netmass2 - netmass1)>1E-6)) then
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: mass balance violation after sweep_tiny_snow!", FATAL)
    endif

    ! sum contributions from melt_and_freeze and from liquid_balance routines
    snow_lprec = snow_lprec1 + snow_lprec2
    snow_hlprec = snow_hlprec1 + snow_hlprec2

    ! sum contributions to lost LAIs for LAI mass balance
    lost_wc_em = lost_wc_em1 + lost_wc_em2 + lost_wc_em3 + lost_wc_em4 + lost_wc_em5
    lost_wc_im = lost_wc_im1 + lost_wc_im2 + lost_wc_im3 + lost_wc_im4 + lost_wc_im5

    netmass1 = s%SWE() ! mass cons check
    netheat1 = s%heat() ! heat cons check


    if(verbose) write(*,*) "SNOW STEP 2 : Do Compaction"
    if (do_compaction) call snow_compaction(s, dt, verbose=.FALSE.) ! it modifies the z - levels


    if(verbose) write(*,*) "SNOW STEP 2 : Do Metamorph"
    if (do_metamorph) call snow_metamorph(s, dt, verbose=.FALSE.)


    if(verbose) write(*,*) "SNOW STEP 2 : Do wind drift"
    if (do_wind_drift) call snow_wind_drift(s, dt, wind_atm, verbose=.FALSE.)

    if(verbose) write(*,*) "SNOW STEP 2 : Do snowpack relayering"
    if (s%nlayers > 0) then
        ! write(*,*) "numbers of snow layers before relayering = ", s%nlayers
        if (do_merge) call s%attempt_merge_layers()
        ! write(*,*) "numbers of snow layers during relayering (before split, after merge) = ", s%nlayers
        if (do_split) call s%attempt_split_layers()
        ! write(*,*) "numbers of snow layers after relayering = ", s%nlayers
    endif


    netmass2 = s%SWE() ! mass cons check
    if (do_snow_check_cons .and.(abs(netmass2 - netmass1)>1E-6)) then
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: mass balance violation after compaction, metamorph, and wind drift!", FATAL)
    endif

    netheat2 = s%heat()
    if (do_snow_check_cons .and.(abs(netheat2 - netheat1)>1E-2)) then
        write(*,*) "before: heat = ", netheat1
        write(*,*) "difference: dheat = ", netheat2 - netheat1
        write(*,*) "norm difference: dheat = ", (netheat2 - netheat1)/netheat1
        call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: heat balance violation after compaction, metamorph, and wind drift!", FATAL)
    endif


   ! here convert between snowpack and lm4p2 energy reference - remove liq*HLF
   snow_hlrunf = snow_hlrunf - snow_lrunf * HLF
   snow_hlprec = snow_hlprec - snow_lprec * HLF
!    vegn_hlprec_r = vegn_hlprec - vegn_lprec * HLF ! ENEROLD
   vegn_hlprec_r = vegn_hlprec ! //TODO cleanup unneccessary var

    ! update snowpack near-surface properties to be saved in diag fields
    call s%nearsurf_properties()

    if(is_watch_point()) then
        write(*,*) "#### snow step 2, final checkpoint"
        call s%print()
        write(*,*) "vegn_lprec, vegn_fprec, ltprec, ftprec = ",vegn_lprec, vegn_fprec, ltprec, ftprec
        write(*,*) "snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, snow_hlrunf, snow_hfrunf",snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, snow_hlrunf, snow_hfrunf
        write(*,*) "snow_levap, snow_fevap, snow_melt = ", snow_levap, snow_fevap, snow_melt
    endif


   ! mass conservation balance for light absorbing impurities
   ! sum LAI lost in sweep tiny snow and snow step 2
   ! Note: quantity not needed, if not for checking LAI mass balance
   lost_wc_em = lost_wc_em1 + lost_wc_em2
   lost_wc_im = lost_wc_im1 + lost_wc_im2
   mass_lai_im_2 = s%lai_im()
   mass_lai_em_2 = s%lai_em()
   total_wetdep_check = 0.0
   if (vegn_lprec > 1E-9) then
      ! total_wetdep_check = total_wetdep_check + wetdep * vegn_lprec / (vegn_lprec + vegn_fprec)
      total_wetdep_check = total_wetdep_check + wetdep * vegn_lprec
      endif
   if (vegn_fprec > 1E-9) then
      ! total_wetdep_check = total_wetdep_check + wetdep * vegn_fprec / (vegn_lprec + vegn_fprec)
      total_wetdep_check = total_wetdep_check + wetdep * vegn_fprec
   endif

   net_delta_lai = mass_lai_im_2 + mass_lai_em_2 - mass_lai_im_1 -  mass_lai_em_1 &
              - drydep*delta_time - total_wetdep_check*delta_time + lost_wc_em + lost_wc_im
   ! if (net_delta_lai(1) > 1E-3) then
   if (do_snow_check_cons .and. ((net_delta_lai(1) > 1E-6).or.(net_delta_lai(2) > 1E-6).or.(net_delta_lai(3) > 1E-6))) then
      write(*,*) "TOTAL DRYDEP = ", sum(drydep)*delta_time
      write(*,*) "TOTAL WETDEP = ", sum(total_wetdep_check)*delta_time
      write(*,*) "TOTAL LOST (IM) = ", sum(lost_wc_im)
      write(*,*) "TOTAL LOST (EM) = ", sum(lost_wc_em)
      write(*,*) "TOTAL LOST = ", sum(lost_wc_em + lost_wc_im)
      write(*,*) "INIT STORAGE = ", sum(mass_lai_em_1 + mass_lai_im_1)
      write(*,*) "FINAL STORAGE = ", sum(mass_lai_em_2 + mass_lai_im_2)
      write(*,*) "Outer rainf, snowf = ", vegn_fprec, vegn_lprec
      write(*,*) "total wetdep due to snowfall + rainfall = ", sum(total_wetdep_check) * delta_time
      write(*,*) "total drydep * dt = ", sum(drydep) * delta_time
      write(*,*) "snowpack nlayers, total ice = ", s%nlayers, s%ice()
      write(*,*) "initial lai mass = ", mass_lai_im_1 + mass_lai_em_1
      write(*,*) "final lai mass = ", mass_lai_im_2 + mass_lai_em_2
      write(*,*) "total deposition: = ", drydep*delta_time + total_wetdep_check*delta_time
      write(*,*) "total LAIs lost (1): = ", lost_wc_em1 + lost_wc_im1
      write(*,*) "total LAIs lost (2): = ", lost_wc_em2 + lost_wc_im2
      write(*,*) "total LAIs lost (1) + (2): = ", lost_wc_em + lost_wc_im
      write(*,*) "net_delta_lai = ", net_delta_lai
      call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: snowpack: light absorbing impurities not conserved after snow step 2", FATAL)
   endif


    endw_check =  s%SWE()
    netw_check = endw_check -  begw_check - (vegn_lprec + vegn_fprec)*delta_time &
                  + (snow_lrunf + snow_frunf + snow_levap + snow_fevap+ snow_lprec)*delta_time
    if (do_snow_check_cons .and. (abs(netw_check) > 1E-6)) then
      write(*,*) "Intial mass = ", begw_check
      write(*,*) "Final mass = ", endw_check
      write(*,*) "Precip = ", (vegn_lprec + vegn_fprec)*delta_time
      write(*,*) "Evap ", (snow_levap + snow_fevap )*delta_time
      write(*,*) "Runoff ", (snow_lrunf + snow_frunf  )*delta_time
      write(*,*) "Infiltration ", (snow_lprec )*delta_time
      call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: Mass balance violated after snow step 2", FATAL)
    endif

    endh_check =  s%heat()
    sum_swheat = sum(s%swheat)
    ! USE ENERGY CONVENTION INTERNAL TO SNOWPACK MODULE FOR LATENT HEAT

    neth_check = endh_check -  heat1b  &
      - (vegn_hlprec_r + vegn_hfprec_ch )*delta_time - vegn_lprec*delta_time*HLF &
      + (snow_hlrunf + snow_hfrunf + snow_lrunf*HLF + snow_hlprec + snow_lprec*HLF )*delta_time &
    !   - Mg_imp*HLF + hfevap*delta_time - delta_heat_DTg &
      - (Mg_imp-subs_M_imp)*HLF + hfevap*delta_time - dheat_fevap
    !   + (snow_G_Z+snow_G_TZ*subs_DT)*delta_time - (  G0  + DGDTg*DTg)*delta_time



    if (do_snow_check_cons .and. (abs(neth_check) > 1E-2)) then
      write(*,*) "-----------------------------------------------------------------------"
      write(*,*) "snow nlayers = ", s%nlayers
      write(*,*) "snow total ice = ", s%ice()
      write(*,*) "snow total liquid = ", s%liq()
      write(*,*) "vegn_fprec * dt =  = ", vegn_fprec * dt
      write(*,*) "vegn_lprec * dt =  = ", vegn_lprec * dt
      write(*,*) "-----------------------------------------------------------------------"
      write(*,*) "Net Energy difference = ", neth_check
      write(*,*) "Sum of internal SW sources * dt = ", sum_swheat*delta_time
      write(*,*) "Intial Energy = ", begh_check
      write(*,*) "Final Energy = ", endh_check
      write(*,*) "(vegn_lprec)*delta_time*HLF = ", (vegn_lprec)*delta_time*HLF
      write(*,*) "(vegn_hlprec)*delta_time = ", (vegn_hlprec_r)*delta_time
      write(*,*) "(vegn_hfprec)*delta_time = ", (vegn_hfprec_ch)*delta_time
      write(*,*) "dheat_fevap", dheat_fevap
      write(*,*) "DE due to implicit melt Mg_imp ", Mg_imp*HLF
      write(*,*) "Mg_imp * HLF = ", Mg_imp*HLF
      write(*,*) "subs_M_imp * HLF = ", subs_M_imp*HLF
      write(*,*) "ground flux", (G0 + DGDTg*DTg)*delta_time
      write(*,*) "soil flux ", (snow_G_Z+snow_G_TZ*subs_DT)*delta_time
      write(*,*) "soil flux - ground flux ", (    snow_G_Z+snow_G_TZ*subs_DT - ( G0 + DGDTg*DTg)  )*delta_time
      write(*,*) "(snow_hlrunf  )*delta_time ", (snow_hlrunf)*delta_time
      write(*,*) "(snow_hfrunf  )*delta_time ", (snow_hfrunf)*delta_time
      write(*,*) "(snow_lrunf*HLF)*delta_time ", (snow_lrunf*HLF)*delta_time
      write(*,*) "(snow_hlprec )*delta_time ", (snow_hlprec )*delta_time
      write(*,*) "(snow_lprec * HLF )*delta_time ", (snow_lprec * HLF )*delta_time
      write(*,*) "-----------------------------------------------------------------------"
      call land_error_message( "ERROR gl_snow_step_2 in snow_evolution module: updated_land_model_fast_0d :: Energy balance violated after snow step 2", FATAL)
    endif

    if(allocated(s%e)) deallocate(s%e)
    if(allocated(s%f)) deallocate(s%f)
    if(allocated(s%swheat)) deallocate(s%swheat)

    ! FIX TO CORRECT SURFACE TEMPERATURE:
    ! SET T=Taverage up to a certain depth depth_taves
    if (correct_surface_T) then
        write(*,*) "Correcting surface T up to depth = ", depth_surface_T_corr
        if (s%nlayers > 0) then
            snow_C = 0.0
            snow_avrg_T = 0.0
            total_depth = 0.0
            do il = 1, s%nlayers
                if (total_depth<depth_surface_T_corr) then
                    total_depth = total_depth + s%snow(il)%dz ! bottom of current layer
                    snow_C = snow_C + s%snow(il)%hCap()
                    snow_avrg_T = snow_avrg_T + (s%snow(il)%T - TFREEZE) * s%snow(il)%hCap()
                endif
            enddo
            snow_avrg_T = TFREEZE + snow_avrg_T / snow_C
            do il = 1, s%nlayers
                if (total_depth<depth_surface_T_corr) then
                    s%snow(il)%T = snow_avrg_T
                endif
            enddo
        endif
    endif
    ! END FIX TO CORRECT NEAR-SURF TEMP

    ! save additional snow variables needed in land model:
    if (s%nlayers > 0) then
        snow_Tbot = s%snow(s%nlayers)%T
        snow_Cbot = s%snow(s%nlayers)%hCap()
        snow_C = 0.0
        snow_avrg_T = 0.0
        do il = 1, s%nlayers
            snow_C = snow_C + s%snow(il)%hCap()
            snow_avrg_T = snow_avrg_T + (s%snow(il)%T - TFREEZE) * s%snow(il)%hCap()
        enddo
        snow_avrg_T = TFREEZE + snow_avrg_T / snow_C
    !    do il=1, s%nlayers
    !         s%snow(il)%T = snow_avrg_T ! // TODO clean up
    !     enddo
    else
        snow_Tbot = TFREEZE
        snow_Cbot = 0.0
        snow_C = 0.0
        snow_avrg_T = TFREEZE
    endif

    !//TODO clean up
    ! s%preprec_surfT = snow_avrg_T

    if(is_watch_point()) then
        write(*,*)'#### Snow step 2 : final check ####'
        call s%print()
    endif


end subroutine gl_snow_step_2


end module snow_evolution_mod
