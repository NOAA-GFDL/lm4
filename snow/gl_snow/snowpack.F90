module snowpack_mod

!----- added for lm4p2 -----
#include <fms_platform.h>
#include "../../shared/debug.inc"

use mpp_mod, only: input_nml_file
use fms_mod, only: error_mesg, check_nml_error, stdlog, mpp_pe, mpp_root_pe, lowercase, &
       string, FATAL, WARNING, NOTE
use land_data_mod,  only : lnd, log_version
use land_debug_mod, only : is_watch_point, land_error_message, check_var_range
use land_constants_mod, only : NBANDS
use constants_mod,  only : tfreeze, hlv, hlf, PI, dens_h2o

use snow_tile_mod, only : N_SNOW_TRACERS, csw, clw, snow_data_area

implicit none
private

public :: snowpack_t
public :: dzopt_t ! it is public only to test optimal thicknesses
public :: snow_layer_type ! need public to update layers with new snowfall
public :: merge_layers
! public :: merge_phases
public :: add_liquid_to_layer
public :: snowpack_init
public :: compute_snow_grain_shape
public :: lap_albedo_include_bc
public :: lap_albedo_include_md
public :: lap_albedo_include_om

public :: LAI_ext, LAI_ssa, eps, rho_water, rho_ice

integer, parameter :: TR_BC = 1      !< Index of black carbon - tracer 1
integer, parameter :: TR_MD = 2      !< Index of mineral dust - tracer 2
integer, parameter :: TR_OM = 3      !< Index of organic carbon - tracer 3

real, parameter :: rho_ice = 917.0 ! ice density [kg / m^3]
! real, parameter :: rho_water = 997.0 ! water density [kg / m^3]
real, parameter :: rho_water = dens_h2o ! water density [kg / m^3]
real, parameter :: rho_refrozen = 300.0 ! refrozen water assumed density [kg / m^3]
real, parameter :: thickness_for_surface_optical_props = 0.03 ! [m] 3cm as in Vionnet et al., 2012 - updated to 5cm
! optical properties od BC, MD and OC (respectively) from Veronica's paper
! using the default value for Dust here - see paper for additional values
! single scattering albedos
real, parameter :: LAI_ssa(N_SNOW_TRACERS) = (/ 0.209, 0.857, 0.963   /) ! single scattering albedo [adim.]
real, parameter :: LAI_ext(N_SNOW_TRACERS) = (/ 9267.0, 474.0, 3289.0 /) ! extinction cross section [m^2 kg^-1]
real, parameter :: LAI_sca(N_SNOW_TRACERS) = (/ 1937.0, 406.0, 3167.0 /) ! scattering cross section [m^2 kg^-1]
real, parameter :: LAI_abs(N_SNOW_TRACERS) = (/ 7330.0, 67.8, 122.0   /) ! absorption cross section [m^2 kg^-1]
real, parameter :: eps = 1E-8 ! a small number


!> \brief state of snow layer
type :: snow_layer_type
    real :: T  !< temperature, K
    real :: ws !< solid phase water, kg/m2 # EZ: changed from density to mass/area
    real :: wl !< liquid phase water, kg/m2 # EZ: changed from density to mass/area
    real :: dz !< layer thickness, m
    real :: optd  !< snow optical diameter [m] [Carmagnola et al., 2013, Flanner and Zender 2006]
    real :: dendr !< snow layer densdricity [dim.less number in [0,1] with 0 = Not dendritic]
    real :: age  !< age of snow layer, [days]
    real :: sph  !< snow grain sphericity [number in [0,1] with 1 = spherical grains]
    real :: wc_em(N_SNOW_TRACERS) !< mass of impurities of each type (array, dim=N_SNOW_TRACERS) - externally mixed only (em) [mg/m2]
    real :: wc_im(N_SNOW_TRACERS) !< mass of impurities of each type (array, dim=N_SNOW_TRACERS) - internally mixed only (im) [mg/m2]
contains
    procedure :: hCap => snow_heat_capacity    !< heat capacity of the layer, J/m2/K
    procedure :: hCon => snow_heat_conductance !< heat conductance of snow, W/m/K
    procedure :: heat => snow_heat_content     ! heat content of the layer, J/m2
    procedure :: density => snow_density     ! snow density of the layer, Kg/m3
end type snow_layer_type

!> \brief State of the snowpack
type :: snowpack_t
    integer :: nlayers = 0 !< number of snow layers
    type(snow_layer_type), allocatable :: snow(:) !< layers of snow
    real, allocatable :: e(:), f(:) !< coefficients for backsubstitution, heat diffusion
    real, allocatable :: swheat(:) ! internal heat source due to shortwave radiation [W m^-2]
    real, allocatable :: sw_frac_dir(:,:) ! frac of sw down absorbed by each layer in (VIS, NIR) bands , direct
    real, allocatable :: sw_frac_dif(:,:) ! frac of sw down absorbed by each layer in (VIS, NIR) bands , diffuse
    real, DIMENSION(NBANDS) :: snow_refl_dir ! direct albedo of snowpack in (VIS, NIR) bands
    real, DIMENSION(NBANDS) :: snow_refl_dif ! direct albedo of snowpack in (VIS, NIR) bands
    real, DIMENSION(NBANDS) :: beta_rad ! penetration length of radiation in snowpack
    real :: topwater ! liquid water temporarily stored on top of snow during soild water balance [kg m^-2]
    real :: topwheat ! heat of topliquid water [J m^-2] with resp. to soild ice at TFREEZE [J m^-2]
    real :: topsnowdeficit ! temporary deficit of snow at the top, due to sublimation [Kg m^-2] - it is a negative mass of soild ice!
    real :: topsnowheatdeficit ! heat deficit connected with  "topsnowdeficit" - energy with opposite sign due to negative mass (~ ws * cs * (T-TF) with ws < 0, so usually [but not if T>0] a positive quantity)
    real :: nearsurf_bceq_tot ! near surface property for albedo calculation - black carbon equivalent concentration - total [ppm]
    real :: nearsurf_bceq_em ! near surface property for albedo calculation - black carbon equivalent concentration [ppm] - externally mixed
    real :: nearsurf_bceq_im ! near surface property for albedo calculation - black carbon equivalent concentration [ppm] - internally mixed
    real :: nearsurf_optd ! near surface property for albedo calculation - snow optical diameter [m]
    real :: nearsurf_sph ! near surface property for albedo calculation - snow grain sphericity [real number in [0,1]]
    real :: nearsurf_dendr ! near surface property for albedo calculation - dendricity [real number in [0,1]]
    real :: nearsurf_rho ! near surface property for albedo calculation - snow density [kg/m^3]
    real :: nearsurf_age ! near surface property for albedo calculation - snow age [days]
    real :: nearsurf_T ! near surface snow temperature [K]
    real :: tag ! tag from lm4p2
    real :: preprec_surfT ! snow surface temperature after heat diffusion, but before evap/subl and new snowfall
contains
!     procedure :: empty  => snowpack_empty !< empty snowpack (in case of complete melt / sublimation)
    procedure :: update_age  => snowpack_update_age !< update the age [in days] of existing snow layers
    procedure :: nearsurf_properties  => snowpack_nearsurf_properties !< compute some near-surface properties
    procedure :: sw_sources  => snowpack_sw_sources !< compute albedo of the snowpack
    procedure :: heat  => snowpack_heat !< heat content of the snowpack
    procedure :: ice   => snowpack_ice  !< solid phase water content of the snowpack
    procedure :: liq   => snowpack_liq  !< liquid phase water content of the snowpack
    procedure :: density   => snowpack_density  !< average density of the snow
    procedure :: avrg_age   => snowpack_avrg_age  !< average age of the snow
    procedure :: avrg_sph   => snowpack_avrg_sph  !< average sphericity of the snow
    procedure :: avrg_optd   => snowpack_avrg_optd  !< average opt. diameter of the snow
    procedure :: avrg_dendr   => snowpack_avrg_dendr  !< average dendricity of the snow
    procedure :: avrg_bceq_im   => snowpack_avrg_bceq_im  !< average conc. of LAIS internally mixed (IM)
    procedure :: avrg_bceq_em   => snowpack_avrg_bceq_em  !< average conc. of LAIS externally mixed (EM)
    procedure :: avrg_bceq_tot  => snowpack_avrg_bceq_tot  !< average conc. of LAIS (IM + EM)
    procedure :: avrg_bc_tot   => snowpack_avrg_bc_tot  !< average conc. of black carbon (IM + EM)
    procedure :: avrg_md_tot   => snowpack_avrg_md_tot  !< average conc. of mineral dust (IM + EM)
    procedure :: avrg_om_tot   => snowpack_avrg_om_tot  !< average conc. of organic carbon (IM + EM)
    procedure :: lai_im   => snowpack_lai_im  !< content of internally mixed LAIs of the snowpack [mg/m2]
    procedure :: lai_em   => snowpack_lai_em  !< content of externally mixed LAIs of the snowpack [mg/m2]
    procedure :: SWE   => snowpack_SWE  !< total water content of the snowpack
    procedure :: depth => snowpack_depth  !< total depth of the snowpack
    procedure :: area => snowpack_area  !< fractional area covered by snow (used only for albedo purposes)
!     procedure :: step1 => snowpack_step_1 !< forward elimination of tridiagonal solver
    procedure :: step1a => snowpack_step_1a !< forward elimination of tridiagonal solver
    procedure :: step1b => snowpack_step_1b !< forward elimination of tridiagonal solver
    procedure :: step2 => snowpack_step_2 !< back-substitution part of tridiagonal solver
    procedure :: check_bounds => snowpack_check_bounds !< sanity check for snowpack variables
    procedure :: attempt_split_layers => attempt_split_layers
    procedure :: attempt_merge_layers => attempt_merge_layers
    procedure :: print => snowpack_print
end type snowpack_t

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snowpack_mod' ! lm4p2
#include "../../shared/version_variable.inc"



!> optimal snowpack layer distribution and associated calculations
integer, parameter :: MAX_OPT_LAYERS = 32
! EZSNOW: made this allocatable
real, ALLOCATABLE :: opt_layer(:)   !< prescribed layer thicknesses
real, ALLOCATABLE :: opt_layer_z(:) !< lower boundary of optimal layers, m

type :: dzopt_t
    integer :: n ! number of elements of z, l
    ! real :: z(MAX_OPT_LAYERS+1) !< boundaries of the layers, m
    ! real :: l(MAX_OPT_LAYERS+1) !< layer number corresponding to layer boundaries
    real, ALLOCATABLE :: z(:) !< boundaries of the layers, m
    real, ALLOCATABLE :: l(:) !< layer number corresponding to layer boundaries
contains
    procedure :: init  => dzopt_init  !< set up optimal layer calculations for specific snow depth
    procedure :: dz    => dzopt_dz    !< calculate optimal layer thickness for given layer depth, m
    procedure :: layer => dzopt_layer !< given depth, calculate (fractional) layer number
    procedure :: depth => dzopt_depth !< given (fractional) layer, calculate depth, m
    procedure :: penalty => dzopt_penalty !< calculate distance from the optimal vertical discretization
    procedure :: print => dzopt_print
end type dzopt_t

! ---- namelist
! integer i
real :: opt_layer_N   = 0.03 !< thickness of the bottom layer, m
real :: opt_layer_max = 1.0  !< maximum optimum layer thickness, m
real :: opt_layer_R   = 1.5  !< factor of increase for the layers in the middle of the snowpack, unitless
! real :: opt_layer(MAX_OPT_LAYERS) = [0.01, (-1.0,i=2,MAX_OPT_LAYERS)] !< prescribed layer thicknesses
logical, protected :: lap_albedo_include_bc = .TRUE.
logical, protected :: lap_albedo_include_md = .TRUE.
logical, protected :: lap_albedo_include_om = .TRUE.
character(len=12) :: heat_cond_to_use = 'yen'  ! available: yen, vapor

namelist /snowpack_nml/ &
         opt_layer_R, opt_layer_N, opt_layer_max, heat_cond_to_use, &
         lap_albedo_include_bc, lap_albedo_include_md, lap_albedo_include_om

! ---- end of namelist

integer :: heat_cond_option = -1
integer, parameter :: &
    HEAT_COND_CAL   = 1, &
    HEAT_COND_VAPOR = 2, &
    HEAT_COND_YEN   = 3

! real :: opt_layer_z(MAX_OPT_LAYERS+1) ! lower boundary of optimal layers, m

contains  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!> initialize optimal layer thickness calculations
! Initialization sets up data arrays for calculation of z(layer) and its inverse layer(z),
! where "layer" is a real number, not an integer. these functions are used to calculate
! optimal thickness of layers for any given depth within the snowpack.
!
! NOTE: the thickness of the lowest layer can be increased with total snowpack depth,
! since we do not expect it to matter when the snow is very deep (and therefore
! likely to be in equilibrium with underlying substrate)
subroutine dzopt_init(dzopt, depth)
  class(dzopt_t), intent(inout) :: dzopt
  real,           intent(in)    :: depth !< snow depth

  real :: dz, d1, scale
  integer :: k

  if (.not. allocated(dzopt%z)) allocate(dzopt%z(size(opt_layer_z)))
  if (.not. allocated(dzopt%l)) allocate(dzopt%l(size(opt_layer_z)))

  call update_dzopt_size(depth)
  ! write(*,*) "size(opt_layer_z), size(dzopt%z))", size(opt_layer_z), size(dzopt%z)
  if(size(opt_layer_z)>size(dzopt%z)) then
    deallocate(dzopt%z)
    deallocate(dzopt%l)
    allocate(dzopt%z( size(opt_layer_z) )) !< boundaries of the layers, m
    allocate(dzopt%l( size(opt_layer_z) )) !< layer number corresponding to layer boundaries
  endif

  if (depth<=0.0) then
     ! for zero-depth snow, just create a dummy distribution of optimal layers so
     ! that dzopt%depth(L) and dzopt%layer(Z) return something
     dzopt%n = 2
     dzopt%z(1) = 0.0; dzopt%z(2) = 1.0e-6 ! slm: this small value sets up a very steep L(Z) slope; is this a problem?
     dzopt%l(1) = 0.0; dzopt%l(2) = 1.0
     if (is_watch_point()) then
        write(*,*) "#### dzopt_init: zero snow depth ####"
        __DEBUG1__(depth)
        call dzopt%print()
     endif
     return
  endif

  dzopt%z(:) = opt_layer_z(:)
  dzopt%l(:) = [(float(k-1), k=1,size(opt_layer_z))]
  dzopt%n    = size(opt_layer_z)

  d1 = depth - opt_layer_N  ! depth to the near-soil layer
  d1 = max(d1,opt_layer(1)) ! to avoid division by zero for snow thinner
                            ! than opt_layer(1)/2-opt_layer_N
  k = bisect(dzopt%z(:), d1)

  if (is_watch_point()) then
     write(*,*) "#### dzopt_init 1 ####"
     __DEBUG4__(depth,d1,opt_layer_N,k)
!      call dzopt%print()
  endif

  if (opt_layer(k)<=opt_layer_N) then
      ! bottom layer is thin enough as it is
      dzopt%n = k+2
  else
      ! scale layers to fit integer number in depth
      if (d1 < opt_layer_z(k)+opt_layer(k)/2) then
         scale = d1/opt_layer_z(k)
      else
         scale = d1/opt_layer_z(k+1)
         k = k+1
      endif
      dzopt%z(:) = dzopt%z(:)*scale

      if(depth > dzopt%z(k)) then
         ! add a thin layer at the bottom
         dzopt%z(k+1) = depth
         dzopt%n      = k+1
      else
         dzopt%n      = k
      endif
  endif

  if (is_watch_point()) then
     write(*,*) "#### dzopt_init 2 ####"
     __DEBUG1__(scale)
     __DEBUG4__(depth,d1,opt_layer_N,k)
     call dzopt%print()
  endif
end subroutine dzopt_init

!>\brief  Finds a position of point in array of bounds. Returns i, such that x is
!!  between xx(i) and xx(i+1).
real function bisect(xx, x)
  real, intent(in)              :: xx(:)     !< array of boundaries
  real, intent(in)              :: x         !< point to locate

  ! ---- local vars
  integer :: low, high, mid
  integer :: n              ! size of the input array
  logical :: ascending      ! if true, the coordinates are in ascending order

  n = size(xx)

  ! find the coordinates
  if (x >= xx(1).and.x<=xx(n)) then
     low = 1; high = n
     ascending = xx(n) > xx(1)
     do while (high-low > 1)
        mid = (low+high)/2
        if (ascending.eqv.xx(mid) <= x) then
           low = mid
        else
           high = mid
        endif
     enddo
     bisect = low
  else if (x < xx(1)) then
     bisect = 1
  else if (x > xx(n)) then
     bisect = n-1
  endif
end function bisect

!\> given depth within snowpack, calculate corresponding (fractional) optimal layer number
real function dzopt_layer(dzopt, d) result (layer)
  class(dzopt_t), intent(in) :: dzopt
  real,           intent(in) :: d

  integer :: i
  ! TODO: check that dzopt is initialized
  associate(z=>dzopt%z, l=>dzopt%l, n=>dzopt%n)
     i = bisect(dzopt%z(1:n),d)
     layer = l(i) + (d-z(i)) * (l(i+1)-l(i))/(z(i+1)-z(i))
  end associate
end function dzopt_layer

!> given layer (perhaps fractional), calculate corresponding depth
real function dzopt_depth(dzopt, layer) result(depth)
  class(dzopt_t), intent(in) :: dzopt
  real,           intent(in) :: layer !< layer number, perhaps fractional

  integer :: i
  ! TODO: check that dzopt is initialized
  associate(z=>dzopt%z, l=>dzopt%l, n=>dzopt%n)
     i = bisect(dzopt%l(1:n),layer)
     depth = z(i) + (layer-l(i)) * (z(i+1)-z(i))/(l(i+1)-l(i))
  end associate
end function dzopt_depth

!> calculate optimal thickness of layer with the top at z
real function dzopt_dz(dzopt,z) result(dz)
  class(dzopt_t), intent(in) :: dzopt
  real,           intent(in) :: z

  dz = dzopt%depth(dzopt%layer(z)+1.0) - z
end function dzopt_dz

!> \brief calculate penalty value for given layer distribution
real function dzopt_penalty(dzopt, z) result(p)
  class(dzopt_t), intent(in) :: dzopt
  real,           intent(in) :: z(:) !< boundaries of the snow layers

  integer :: k
  real :: dz_opt, dz

  p = 0
  do k = 1, size(z)-1
     dz_opt = dzopt%dz(z(k))
     dz     = z(k+1)-z(k)
     p = p + (dz-dz_opt)**2
  enddo
end function dzopt_penalty

!> \brief print optimal vertical discretization
subroutine dzopt_print(dzopt)
  class(dzopt_t), intent(in) :: dzopt

  integer :: k
  write(*,'(a2,99(",",a9,:))') "k", "z","l","dz_opt"
  do k = 1, dzopt%n
     write(*,'(i2.2,99(",",f9.4,:))') k, dzopt%z(k), dzopt%l(k), dzopt%dz(dzopt%z(k))
  enddo
end subroutine dzopt_print

!> initialize snowpack module, in particular read namelist parameters
! version for lm4p2
subroutine snowpack_init()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  character(256) :: msg   ! namelist error message text
  integer :: k, n
  real    :: dz ! layer thickness, for initialization of optimal vertical discretization, m

  call log_version(version, module_name, &
  __FILE__)
  read (input_nml_file, nml=snowpack_nml, iostat=io, iomsg=msg)
  ierr = check_nml_error(io, 'snowpack_nml :: '//trim(msg))
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=snowpack_nml)
  endif

  if (trim(lowercase(heat_cond_to_use))=='cal') then
     heat_cond_option = HEAT_COND_CAL
  else if (trim(lowercase(heat_cond_to_use))=='vapor') then
     heat_cond_option = HEAT_COND_VAPOR
  else if (trim(lowercase(heat_cond_to_use))=='yen') then
     heat_cond_option = HEAT_COND_YEN
  else
     call error_mesg('read_snowpack_namelist', &
        'heat_cond_to_use='//trim(heat_cond_to_use)//' in snowpack_nml in incorrect: valid options are "Cal", "vapor", or "Yen"', FATAL)
  endif

  !!! ------ EZSNOW : made these allocatable ------ !!!
  ! allocate(z(MAX_OPT_LAYERS+1)) !< boundaries of the layers, m
  ! allocate(l(MAX_OPT_LAYERS+1)) !< layer number corresponding to layer boundaries
  allocate(opt_layer(MAX_OPT_LAYERS))  !< prescribed layer thicknesses
  allocate(opt_layer_z(MAX_OPT_LAYERS+1)) !< lower boundary of optimal layers, m
  opt_layer = [0.05, (-1.0,k=2,MAX_OPT_LAYERS)] !< start assigning prescribed layer thicknesses
  !!! ---------------------------------------------- !!!

  ! initialize optimal layer distribution for infinite lower bound
  ! using prescribed thicknesses for the shallow snow depths, and
  ! limited exponential increase of each sequential layer thickness below

  ! check that there are positive values in layer thickness array
  n = count(opt_layer(:)>0)
  if (n == 0) &
     call land_error_message( 'No positive values in layer thickness array "opt_layer"', FATAL)
  if (count(opt_layer(1:n)>0) < n) &
     call land_error_message( 'Positive layer thickness values "opt_layer" are intermingled with negatives', FATAL)



  opt_layer_z(1) = 0.0
  do k = 1, size(opt_layer)
     if (opt_layer(k) > 0) then
        dz = opt_layer(k)
     else
        dz = min(dz*opt_layer_R, opt_layer_max)
        opt_layer(k) = dz ! store for future use
     endif
     opt_layer_z(k+1) = opt_layer_z(k) + dz
  enddo

  ! Now, add layers in case snowpack is too thick and requires a larger number of optimal layers

  ! write(*,*) """

!   write(*,'(99(i8.2,:,","))') (k, k=1,size(opt_layer))
!   write(*,'(99(f8.3,:,","))') opt_layer
end subroutine snowpack_init

! if snowpack is too thick, increase number of layers in dzopt
subroutine update_dzopt_size(snow_depth)
  real snow_depth
  integer i
  integer :: io, k, n
  real    :: dz ! layer thickness, for initialization of optimal vertical discretization, m
  integer :: NEW_OPT_LAYERS
  integer cuml
  real cumz


  ! write(*,*) "update dzopt size: depth, opt_layer_z = ", snow_depth, opt_layer_z
  if (snow_depth > opt_layer_z(size(opt_layer_z)-1)) then

  ! write(*,*) "updating the size of opt_layer_z array ..."


    !!! ------ Determine new number of layers needed for dzopt
    cumz = opt_layer_z(size(opt_layer_z)) ! max depth in opt layer z, bottom of opt snowpack
    cuml = size(opt_layer) ! current maximum number of layers in optimum snowpack
    dz = opt_layer(size(opt_layer)) ! current optimal thickness of bottom opt layer
    do while(cumz<snow_depth)
      dz = min(dz*opt_layer_R, opt_layer_max)
      cumz = cumz + dz
      cuml = cuml + 1
    enddo

    !!! ------ Now cuml is the new number of layers needed
    !!! add 5 to reduce the number of times the arrays need to be allocated
    NEW_OPT_LAYERS = cuml + 5

    !!! ------ EZSNOW : made these allocatable ------ !!!
    ! deallocate(dzopt%z)
    ! deallocate(dzopt%l)
    deallocate(opt_layer)
    deallocate(opt_layer_z)
    ! allocate(dzopt%z(NEW_OPT_LAYERS+1)) !< boundaries of the layers, m
    ! allocate(dzopt%l(NEW_OPT_LAYERS+1)) !< layer number corresponding to layer boundaries
    allocate(opt_layer(NEW_OPT_LAYERS))  !< prescribed layer thicknesses
    allocate(opt_layer_z(NEW_OPT_LAYERS+1)) !< lower boundary of optimal layers, m
    opt_layer = [0.05, (-1.0,i=2,NEW_OPT_LAYERS)] !< start assigning prescribed layer thicknesses
    !!! ---------------------------------------------- !!!

    ! initialize optimal layer distribution for infinite lower bound
    ! using prescribed thicknesses for the shallow snow depths, and
    ! limited exponential increase of each sequential layer thickness below

    ! check that there are positive values in layer thickness array
    n = count(opt_layer(:)>0)
    if (n == 0) &
      call land_error_message( 'No positive values in layer thickness array "opt_layer"', FATAL)
    if (count(opt_layer(1:n)>0) < n) &
      call land_error_message( 'Positive layer thickness values "opt_layer" are intermingled with negatives', FATAL)

    opt_layer_z(1) = 0.0
    do k = 1, size(opt_layer)
      if (opt_layer(k) > 0) then
          dz = opt_layer(k)
      else
          dz = min(dz*opt_layer_R, opt_layer_max)
          opt_layer(k) = dz ! store for future use
      endif
      opt_layer_z(k+1) = opt_layer_z(k) + dz
    enddo

  endif

  if (opt_layer_z(size(opt_layer_z))< snow_depth) then
      write(*, * ) "Optimal layer thickness distribution error: snow depth = ", snow_depth
      write(*, * ) "Optimal layer thickness distribution error: NEW_OPT_LAYERS  = ", NEW_OPT_LAYERS
      call land_error_message("Optimal layer thickness distribution is not thick enough for given snow depth", FATAL)
  endif

end subroutine update_dzopt_size


!> update age of existing snow layers [in days]
subroutine snowpack_update_age(s, dt)
  ! add to snow age [in days] the current model time step [dt in seconds!]
  class(snowpack_t), intent(inout) :: s
  integer il
  real dt ! time incremenet [in seconds]
  if (s%nlayers > 0) then
    do il=1,s%nlayers
      s%snow(il)%age = s%snow(il)%age + dt/86400.0
    enddo
  endif
end subroutine snowpack_update_age


!> empty snowpack (to use in case of complete melt or sublimation)
subroutine snowpack_empty(s)
  class(snowpack_t) :: s
  deallocate(s%snow)
  s%nlayers = 0
  ! do nothing else
end subroutine snowpack_empty

!> \brief Calculate snow layer heat capacity, J/m2/K
real function snow_heat_capacity(snow) result(res)
  class(snow_layer_type), intent(in) :: snow !< snow state
!   res = (CSW*snow%ws + CLW*snow%wl)*snow%dz
  res = CSW*snow%ws + CLW*snow%wl
end function snow_heat_capacity

!> \brief Calculate snow layer heat content, J/m2
!!
!! snow heat content is calculated relative to the heat content of solid snow at
!! freezing temperature, with the same mass as total water.
real function snow_heat_content(snow) result(res)
  class(snow_layer_type), intent(in) :: snow !< snow state
  res = snow%hCap()*(snow%t-TFREEZE) + snow%wl * HLF
end function snow_heat_content

!> \brief Calculate heat conductance of snow W/m/K
real function snow_heat_conductance(snow, surfT, surfP) result(res)
  class(snow_layer_type), intent(in) :: snow
  real, INTENT(IN) :: surfT, surfP ! atm temperature [K] and pressure [Pa]
  real, parameter :: al = 2.22   ! [ W m^-1 K^-1]
  real, parameter :: lmin = 4E-2 ! [W m^-2 K^-2]
  ! real, parameter :: lmin = 9E-2 ! [W m^-2 K^-2] ! CM value
  real, parameter :: expon = 1.88
  real rho_snow , temp_snow, press_srfc
  ! should I also consider the density of liquid water in the snow layer
  ! rho_snow = snow%ws/snow%dz ! [kg m^-3]
      real a1, a2, a3, b1, b2, b3, P0, kc, kwv
      ! real surfT, surfP ! move to input

! //TODO: remove surfT as input, not needed

  ! if (use_cm_conductance) then
  select case(heat_cond_option)
  case (HEAT_COND_CAL)
    rho_snow = (snow%ws + snow%wl)/snow%dz ! [kg m^-3]
    a1 = 2.5E-6
    a2 = 1.23E-4
    a3 = 0.024
    res = a1*rho_snow**2 - a2*rho_snow + a3
  case (HEAT_COND_VAPOR)
      ! res = 0.0005 ! as in CM model
      ! res = 50.0 ! as in CM model
      ! res = 0.0003 ! as in CM model
      ! res = 0.09 ! as in CM model

      ! USE MODEL BY CALONNE 2011 + WATER VAPOR BY SUN 1991 [CM = Calonne Modified with water vapor]
      a1 = 2.5E-6
      a2 = 1.23E-4
      a3 = 0.024
      P0 = 1000.0 ! reference pressure [hPa]
      ! surfT = 273.15 ! move to input
      ! surfP = 1000.0 ! move to input
      b1 = -0.06023 ! [W m^-1 K^-1]
      b2 = 2.5425 ! [W m^-1]
      b3 = 289.99 ! [K]
      rho_snow = (snow%ws + snow%wl)/snow%dz ! [kg m^-3]
      temp_snow = max(240.0, min(snow%T, 280.0))
      press_srfc = max(50000.0, min(surfP, 110000.0) )
      kc = a1*rho_snow**2 - a2*rho_snow + a3
      ! surfP input here is given in [Pa], convert to [hPa]
      kwv = 100000.0/surfP * (b1-b2/( temp_snow-b3))
      res = kc + kwv

      ! if (temp_snow>285.0) then
      !   call land_error_message( "ERROR in snow_heat_conductance in snowpack_mod :: snow temperature is too high!", FATAL)
      ! endif

      if (is_watch_point()) then
        write(*,*) "#### snow heat conductance - water vapor effect calculations :: ####"
        write(*,*) "surfT, surfP = ", surfT, surfP
        write(*,*) "kc, kwv, kc+kwv = ", kc, kwv, res
      endif
  case (HEAT_COND_YEN)
      rho_snow = (snow%ws + snow%wl)/snow%dz ! [kg m^-3]
      res = max(lmin, al*(rho_snow/rho_water)**expon) ! YEN 1981, CROCUS
      ! res = 0.023 + (7.75*1E-5 * rho_snow + 1.105 * 1E-6 * rho_snow**2) * (2.29 - 0.023) ! JORDAN 1991, SHRESTA 2006

  case default
      call land_error_message("Error in snow_heat_conductance in snowpack_mod: Must specify a valid snow_heat_cond_to_use!", FATAL)
  end select
  ! write(*,*) "snow heat conductance = ", res
  ! parameterization by Yen (1981), used by default in CROCUS
  ! see Lafaysse et al., 2017
end function snow_heat_conductance

!> \brief Heat content of the snowpack, J/m2
!! note: sum the heat content of each layer
!! add the heat of any "topwater" (water on top of snowpack)
!! add the heat deficit of any snow mass deficit temporarily present
real function snowpack_heat(snowpack) result(heat)
  class(snowpack_t), intent(in) :: snowpack
  integer :: k
  heat = 0.0
  heat = heat + snowpack%topwheat + snowpack%topsnowheatdeficit
  if (snowpack%nlayers>0) then
    do k = 1, snowpack%nlayers
        heat = heat+snowpack%snow(k)%heat()
    enddo
  endif
end function snowpack_heat

!> \brief Solid phase water content of the snowpack, kg/m2
!! add any snow deficit present (negative solid mass)
real function snowpack_ice(s) result(sum)
  class(snowpack_t), intent(in) :: s
  integer :: k
  sum = 0
  sum = sum + s%topsnowdeficit
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%ws
    enddo
  endif
end function snowpack_ice

!> \brief Liquid phase water content of the snowpack, kg/m2
!! add any topwater, if present (liquid water temporarily "on top" of the snowpack)
real function snowpack_liq(s) result(sum)
  class(snowpack_t), intent(in) :: s
  integer :: k
  sum = 0
  sum = sum + s%topwater
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%wl
    enddo
  endif
end function snowpack_liq

!> \brief Total water content of the snowpack, kg/m2
real function snowpack_SWE(s)
  class(snowpack_t), intent(in) :: s
  ! topwater and topsnow deficit are now already accounted for in ice and liq!
  snowpack_SWE = snowpack_ice(s)+snowpack_liq(s)
end function snowpack_SWE



!> \brief density of the snowpack, kg/m3
real function snowpack_density(s) result(dens)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real mass, depth
  mass = s%ice()+s%liq()
  ! mass = s%ice()
  depth = s%depth()
  if (depth>0) then
    dens = mass/depth
  else
    dens = 0.0
  endif
end function snowpack_density


!> \brief density of a snow layer, kg/m3
real function snow_density(s) result(dens)
  class(snow_layer_type), intent(in) :: s
  dens = (s%ws + s%wl)/s%dz
end function snow_density


!> \brief average age of the snowpack, days
! consider only the age of the solid phase in each layer
! and do weighted average
real function snowpack_avrg_age(s) result(age)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real sum,sum_ws
  sum = 0.0
  sum_ws = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%age * s%snow(k)%ws
      sum_ws = sum_ws + s%snow(k)%ws
    enddo
  endif
  if (sum_ws>0) then
    age = sum/sum_ws
  else
    age = 0.0
  endif
end function snowpack_avrg_age


!> \brief average sphericity of the snowpack, [dimless]
! consider only the weight of the solid phase in each layer
! and do weighted average
real function snowpack_avrg_sph(s) result(sph)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real sum,sum_ws
  sum = 0.0
  sum_ws = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%sph * s%snow(k)%ws
      sum_ws = sum_ws + s%snow(k)%ws
    enddo
  endif
  if (sum_ws>0) then
    sph = sum/sum_ws
  else
    sph = 0.0
  endif
end function snowpack_avrg_sph


!> \brief average optical diameter of the snowpack, [m]
! consider only the weight of the solid phase in each layer
! and do weighted average
real function snowpack_avrg_optd(s) result(optd)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real sum,sum_ws
  sum = 0.0
  sum_ws = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%optd * s%snow(k)%ws
      sum_ws = sum_ws + s%snow(k)%ws
    enddo
  endif
  if (sum_ws>0) then
    optd = sum/sum_ws
  else
    optd = 0.0
  endif
end function snowpack_avrg_optd


!> \brief average dendricity of the snowpack
! consider only the weight of the solid phase in each layer
! and do weighted average
real function snowpack_avrg_dendr(s) result(dendr)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real sum,sum_ws
  sum = 0.0
  sum_ws = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%dendr * s%snow(k)%ws
      sum_ws = sum_ws + s%snow(k)%ws
    enddo
  endif
  if (sum_ws>0) then
    dendr = sum/sum_ws
  else
    dendr = 0.0
  endif
end function snowpack_avrg_dendr


!> \brief average concentration of black carbon (IM + EM) in the snowpack [ppm]
real function snowpack_avrg_bc_tot(s) result(res)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real sum,sum_wsl
  sum = 0.0
  sum_wsl = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%wc_im(1) + s%snow(k)%wc_em(1)
      sum_wsl = sum_wsl + s%snow(k)%ws + s%snow(k)%wl
    enddo
  endif
  if (sum_wsl>0) then
    res = sum/sum_wsl
  else
    res = 0.0
  endif
end function snowpack_avrg_bc_tot


!> \brief average concentration of mineral dust (IM + EM) in the snowpack [ppm]
real function snowpack_avrg_md_tot(s) result(res)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real sum,sum_wsl
  sum = 0.0
  sum_wsl = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%wc_im(2) + s%snow(k)%wc_em(2)
      sum_wsl = sum_wsl + s%snow(k)%ws + s%snow(k)%wl
    enddo
  endif
  if (sum_wsl>0) then
    res = sum/sum_wsl
  else
    res = 0.0
  endif
end function snowpack_avrg_md_tot


!> \brief average concentration of organic carbon (IM + EM) in the snowpack [ppm]
real function snowpack_avrg_om_tot(s) result(res)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real sum,sum_wsl
  sum = 0.0
  sum_wsl = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%wc_im(3) + s%snow(k)%wc_em(3)
      sum_wsl = sum_wsl + s%snow(k)%ws + s%snow(k)%wl
    enddo
  endif
  if (sum_wsl>0) then
    res = sum/sum_wsl
  else
    res = 0.0
  endif
end function snowpack_avrg_om_tot


!> \brief average concentration of LAIs (IM + EM) in the snowpack [ppm, equivalent black carbon]
real function snowpack_avrg_bceq_tot(s) result(res)
  class(snowpack_t), intent(in) :: s
  res = s%avrg_bceq_em() + s%avrg_bceq_im()
end function snowpack_avrg_bceq_tot


!> \brief average concentration of LAIs internally mixed (IM) in the snowpack [ppm, equivalent black carbon]
real function snowpack_avrg_bceq_im(s) result(res)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real sum,sum_wsl
  sum = 0.0
  sum_wsl = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%wc_im(1) + s%snow(k)%wc_im(2)*(LAI_abs(2)/LAI_abs(1)) + &
                                       s%snow(k)%wc_im(3)*(LAI_abs(3)/LAI_abs(1))
      sum_wsl = sum_wsl + s%snow(k)%ws + s%snow(k)%wl
      ! sum_wsl = sum_wsl + s%snow(k)%ws
    enddo
  endif
  if (sum_wsl>0) then
    res = sum/sum_wsl
  else
    res = 0.0
  endif
end function snowpack_avrg_bceq_im


!> \brief average concentration of LAIs externally mixed (EM) in the snowpack [ppm, equivalent black carbon]
real function snowpack_avrg_bceq_em(s) result(res)
  class(snowpack_t), intent(in) :: s
  integer :: k
  real sum,sum_wsl
  sum = 0.0
  sum_wsl = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      sum = sum + s%snow(k)%wc_em(1) + s%snow(k)%wc_em(2)*(LAI_abs(2)/LAI_abs(1)) + &
                                       s%snow(k)%wc_em(3)*(LAI_abs(3)/LAI_abs(1))
      sum_wsl = sum_wsl + s%snow(k)%ws + s%snow(k)%wl
      ! sum_wsl = sum_wsl + s%snow(k)%ws
    enddo
  endif
  if (sum_wsl>0) then
    res = sum/sum_wsl
  else
    res = 0.0
  endif
end function snowpack_avrg_bceq_em


!> \Internally mixed LAIs content of the snowpack [mg/m2]
function snowpack_lai_im(s)
  class(snowpack_t), intent(in) :: s
  real, DIMENSION(N_SNOW_TRACERS) :: snowpack_lai_im
  integer :: k
  snowpack_lai_im = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      snowpack_lai_im = snowpack_lai_im + s%snow(k)%wc_im(:)
    enddo
  endif
end function snowpack_lai_im


!> \Externally mixed LAIs content of the snowpack [mg/m2]
function snowpack_lai_em(s)
  class(snowpack_t), intent(in) :: s
  real, DIMENSION(N_SNOW_TRACERS) :: snowpack_lai_em
  integer :: k
  snowpack_lai_em = 0.0
  if (s%nlayers>0) then
    do k = 1,s%nlayers
      snowpack_lai_em = snowpack_lai_em + s%snow(k)%wc_em(:)
    enddo
  endif
end function snowpack_lai_em


!> \Depth of the snowpack, m
real function snowpack_depth(s)
  class(snowpack_t), intent(in) :: s
  integer :: il
  snowpack_depth = 0
  if (s%nlayers > 0) then
    do il = 1,s%nlayers
      snowpack_depth = snowpack_depth + s%snow(il)%dz
    enddo
  endif
end function snowpack_depth

!> \compute snow area - used only for albedo purposes for now
real function snowpack_area(snowpack)
  class(snowpack_t), intent(in) :: snowpack !< state of snowpack

  call snow_data_area (snowpack%depth(), snowpack_area)
end function snowpack_area


!> \Print some summary information about the state of the snowpack
subroutine snowpack_print(s)
  class(snowpack_t), intent(in) :: s
  real :: z
  integer :: k
write(*,*) "___________<< state of snowpack >>______________"
  if (s%nlayers > 0) then
  write(*,'(a2,99(",",a14,:))') "k","top","dz","T - TF","ws","wl", "rho", "wc_BC_em", "wc_MD_em", "wc_OM_em", "wc_BC_im", "wc_MD_im", "wc_OM_im", "dendr", "sph", "dopt", "age"
  z = 0
  do k = 1, s%nlayers
     write(*,'(i2.2,99(",",f14.4,:))') k, z, s%snow(k)%dz, s%snow(k)%T-TFREEZE, &
        s%snow(k)%ws, s%snow(k)%wl, s%snow(k)%ws/s%snow(k)%dz, &
        s%snow(k)%wc_em(TR_BC), s%snow(k)%wc_em(TR_MD), s%snow(k)%wc_em(TR_OM), &
        s%snow(k)%wc_im(TR_BC), s%snow(k)%wc_im(TR_MD), s%snow(k)%wc_im(TR_OM), &
        s%snow(k)%dendr, s%snow(k)%sph, s%snow(k)%optd, s%snow(k)%age
     z = z+s%snow(k)%dz
  enddo
  write(*,'("nlayers = ",i2.2)') s%nlayers
  write(*,'("depth = ",f9.4)') z
  write(*,'("SWE = ",f9.4)') s%SWE()
  write(*,'("size of s%snow = ",i2.2)') size(s%snow)
write(*,'("bands    : " 99(a15,:))') "VIS", "NIR"
write(*,'("refl dir = ", 99(f15.4,:))') s%snow_refl_dir(1), s%snow_refl_dir(2)
write(*,'("refl dif = ", 99(f15.4,:))') s%snow_refl_dif(1), s%snow_refl_dif(2)
write(*,'("beta rad = ", 99(f15.4,:))') s%beta_rad(1), s%beta_rad(2)
else
  write(*,*) "There is no snow here at this time -> nlayers = 0"
endif
! These variables are relevant even if there are no snow layers
  write(*,'("topwater, topwheat = ", 99(f15.4,:))') s%topwater, s%topwheat
  write(*,'("topsnow def, topsnow heat def = ", 99(f15.4,:))') s%topsnowdeficit, s%topsnowheatdeficit
write(*,*) "___________<< end state of snowpack >>______________"
end subroutine snowpack_print

!> \Check that relevant variables are within physical bounds, if not throw an error
subroutine snowpack_check_bounds(s, message)
  class(snowpack_t), intent(inout) :: s !< state of snowpack
  integer il
  real rho_snow_il
  character(*), optional :: message
  character(200) :: def_message

  logical bound_exceeded
  bound_exceeded = .FALSE.

  def_message = 'end ---> checking bounds of some relevant variables'

  if (present(message)) then
    def_message = message
  endif

  ! customize message to identify check

  ! check bounds of a few variables

  if (s%nlayers > 0) then
  do il = 1, s%nlayers
    if (s%snow(il)%ws > 1E-3) then
      ! rho_snow_il = (s%snow(il)%ws +s%snow(il)%wl)/s%snow(il)%dz
      rho_snow_il = (s%snow(il)%ws)/s%snow(il)%dz ! there can be temporary water during the time step
      if (rho_snow_il<0.0 .or. rho_snow_il > rho_water) then
        bound_exceeded = .TRUE.
        write(*,*) "Detected snow density out of bounds: rho snow = ", rho_snow_il
        write(*,*) "ws > 0:", s%snow(il)%ws
        write(*,*) "dz > 0:", s%snow(il)%dz
        ! call s%print()
        ! error stop "snow density out of bounds!"
      endif
    endif
  enddo
  endif

  if (s%nlayers > 0) then
  do il = 1, s%nlayers
    if (s%snow(il)%ws > 1E-3) then
      rho_snow_il = (s%snow(il)%ws +s%snow(il)%wl)/s%snow(il)%dz
      if ( s%snow(il)%T < 0.0 .or. s%snow(il)%T > 400.0  ) then
        bound_exceeded = .TRUE.
        write(*,*) "Detected T out of bounds: s%snow(il)%T = ", s%snow(il)%T
        write(*,*) "ws > 0:", s%snow(il)%ws
        write(*,*) "dz > 0:", s%snow(il)%dz
        write(*,*) "rho snow = ", rho_snow_il
        ! call s%print()
        ! error stop "snow temperaure out of bounds!"
      endif
    endif
  enddo
  endif

if (bound_exceeded) then
  ! def_message = 'end ---> checking bounds of some relevant variables'
  call s%print()
  write(*,*) def_message
  ! error stop "ERROR in snwoapck_check_bounds in snowpack module"
  call land_error_message( "ERROR in snowpack_check_bounds in snowpack module", FATAL)
endif

end subroutine snowpack_check_bounds


!> \brief Forward elimination pass of tridiagonal solver
! Deprecated in lm4p2 to account for internal sw sources, split in 2 parts (step1 a and b)
subroutine snowpack_step_1(s, G0, DGDT, &
                          ! output -- additional output  for lm4p2
                          snow_active, snow_T1, snow_rh, snow_liq, snow_ice, snow_subl, snow_area, F0, DFDT, &
                          t_atm, p_surf ,&
                           dt) ! dt need to be passed only in standalone model, global var in lm4p2 [delta_time]
  class(snowpack_t), intent(inout) :: s !< state of snowpack
  real, intent(in)  :: t_atm, p_surf      !< atm temp [K] and pressure [Pa]
  real, intent(in)  :: dt      !< time step, s
  real, intent(in)  :: G0      !< ground heat flux, positive downward, W/m2
  real, intent(in)  :: DGDT    !< derivative of G0 w.r.t. surface temperature, W/m2/K
  ! EZDEV: now they are stored in the snowpack structure
  real, intent(out) :: F0      !< downward heat flux from surface, W/m2
  real, intent(out) :: DFDT    !< derivative of snow surface flux w.r.t. temperature, W/m2/K
  ! other output variables needed in lm4p2 snow step 1
  real, intent(out) :: snow_T1 ! temperature of snow top layer [K]
  real, intent(out) :: snow_rh ! relative humidity, set to 1 currently in lm4p2
  real, intent(out) :: snow_ice ! total solid mass of snowpack (avail. for implicit melt)
  real, intent(out) :: snow_liq ! total liquid of top layer only (avail. for implicit freeze)
  real, intent(out) :: snow_subl ! sublimation (=1) or liq evap (=0) or maybe fractional?
  real, intent(out) :: snow_area ! as in lm4p2, for albedo only
  ! real, intent(out) :: snow_E_max ! max sublim flux from snowpack [kg/(m2 s)]

  integer :: N ! shorthand for number of layers
  logical, intent(out) :: snow_active ! True if snowpack is there (nlayers >= 0)

  N = s%nlayers
  if (allocated(s%e)) deallocate(s%e)
  if (allocated(s%f)) deallocate(s%f)

  if (N > 0) then
    snow_active = .TRUE.
    allocate(s%e(N-1))
    allocate(s%f(N-1))
    ! call step1(s%snow(1:N), N, dt, G0, DGDT, s%swheat, F0, DFDT, s%e(1:N-1), s%f(1:N-1))
    call step1(s%snow(1:N), N, dt, G0, DGDT, s%swheat, F0, DFDT, s%e(1:N-1), s%f(1:N-1), t_atm, p_surf)

    ! compute additional variables needed in LM4P2 snow_step1
      snow_T1 = s%snow(1)%T
      snow_liq = s%snow(1)%wl ! set this to zero to turn off implicit melt in lm4p2
      snow_rh = 1.0
!      snow_ice = s%ice() ! set this to zero to turn off implicit melt in lm4p2
      snow_ice = s%ice() ! set this to zero to turn off implicit melt in lm4p2
      snow_subl = 1.0
      snow_area = s%area()
      ! snow_E_max = s%snow(1)%ws/dt

  else ! no snow
  snow_active = .FALSE.
    snow_T1 = TFREEZE ! value not defined - do not pass it?
    snow_liq = 0.0
    snow_ice = 0.0
    snow_rh = 1.0 ! value not defined - do not pass it?
    snow_area = 0.0
    snow_subl = 0.0
    ! snow_E_max = 0.0
    F0 = G0 ! no snow, same flux passed to the underlying soil/lake/glac layer directly
    DFDT = DGDT ! no snow, same flux passed to the underlying soil/lake/glac layer directly
  endif
end subroutine snowpack_step_1


!> \brief Forward elimination pass of tridiagonal solver
subroutine snowpack_step_1a(s, &
                            ! G0, DGDT, &
                          ! output -- additional output  for lm4p2
                          snow_active, snow_T1, snow_rh, snow_liq, snow_ice, snow_subl, snow_area, snow_E_max, &
                          ! F0, DFDT, &
                           dt, do_mgimplicit, grnd_T) ! dt need to be passed only in standalone model, global var in lm4p2 [delta_time]
  class(snowpack_t), intent(inout) :: s !< state of snowpack
  real, intent(in)  :: dt      !< time step, s
  logical, intent(in)  :: do_mgimplicit      !< if true, do the "implicit" snow melt
  real, intent(in)  :: grnd_T ! temperature of underlying soil surface [K]
  ! real, intent(in)  :: G0      !< ground heat flux, positive downward, W/m2
  ! real, intent(in)  :: DGDT    !< derivative of G0 w.r.t. surface temperature, W/m2/K
  ! EZDEV: now they are stored in the snowpack structure
  ! real, intent(out) :: F0      !< downward heat flux from surface, W/m2
  ! real, intent(out) :: DFDT    !< derivative of snow surface flux w.r.t. temperature, W/m2/K
  ! other output variables needed in lm4p2 snow step 1
  real, intent(out) :: snow_T1 ! temperature of snow top layer [K]
  real, intent(out) :: snow_rh ! relative humidity, set to 1 currently in lm4p2
  real, intent(out) :: snow_ice ! total solid mass of snowpack (avail. for implicit melt)
  real, intent(out) :: snow_liq ! total liquid of top layer only (avail. for implicit freeze)
  real, intent(out) :: snow_subl ! sublimation (=1) or liq evap (=0) or maybe fractional?
  real, intent(out) :: snow_area ! as in lm4p2, for albedo only

  ! integer :: N ! shorthand for number of layers
  logical, intent(out) :: snow_active ! True if snowpack is there (nlayers >= 0)
  real, intent(out) :: snow_E_max ! max sublim flux from snowpack [kg/(m2 s)]

      integer il

  ! N = s%nlayers
  ! if (allocated(s%e)) deallocate(s%e)
  ! if (allocated(s%f)) deallocate(s%f)

  ! call s%force_split_layers(3) ! if thers is snow, make sure it is split in at least (3) layers

  if (s%nlayers > 0) then
    snow_active = .TRUE.
    ! allocate(s%e(N-1))
    ! allocate(s%f(N-1))
    ! call step1(s%snow(1:N), N, dt, G0, DGDT, s%swheat, F0, DFDT, s%e(1:N-1), s%f(1:N-1))

    ! compute additional variables needed in LM4P2 snow_step1
      snow_T1 = s%snow(1)%T
      ! snow_T1 = s%snow(s%nlayers)%T ! // TODO, not currently used, clean
      ! instead of 1st layer, use average temperature of top 3 cm?
    ! call s%nearsurf_properties()
      ! snow_T1 = s%nearsurf_T


  !  write(*,*) "-------------------------------------------"
  !  write(*,*) "AVRG T BEFORE STEP 1 = ", s%avrg_T() ! // TODO not currently used, clean


      snow_liq = s%snow(1)%wl
      snow_rh = 1.0
      if (do_mgimplicit) then
        snow_ice = s%ice() ! set this to zero to turn off implicit melt in lm4p2
      else
        snow_ice = 0.0
      endif
      snow_subl = 1.0
      snow_area = s%area()
      snow_E_max = s%snow(1)%ws/dt

  else ! no snow
  snow_active = .FALSE.
    snow_T1 = TFREEZE ! value not defined - do not pass it?
    snow_liq = 0.0
    snow_ice = 0.0
    snow_rh = 1.0 ! value not defined - do not pass it?
    snow_area = 0.0
    snow_subl = 0.0
    snow_E_max = 0.0
    ! F0 = G0 ! no snow, same flux passed to the underlying soil/lake/glac layer directly
    ! DFDT = DGDT ! no snow, same flux passed to the underlying soil/lake/glac layer directly
  endif
end subroutine snowpack_step_1a



!> \brief Forward elimination pass of tridiagonal solver
subroutine snowpack_step_1b(s, G0, DGDT, &
                          ! output -- additional output  for lm4p2
                          ! snow_active, snow_T1, snow_rh, snow_liq, snow_ice, snow_subl, snow_area,
                          F0, DFDT, &
                          t_atm, p_surf,&
                           dt) ! dt need to be passed only in standalone model, global var in lm4p2 [delta_time]
  class(snowpack_t), intent(inout) :: s !< state of snowpack
  real, intent(in)  :: t_atm, p_surf      !< atm temp [K] and pressure [Pa]
  real, intent(in)  :: dt      !< time step, s
  real, intent(in)  :: G0      !< ground heat flux, positive downward, W/m2
  real, intent(in)  :: DGDT    !< derivative of G0 w.r.t. surface temperature, W/m2/K
  ! EZDEV: now they are stored in the snowpack structure
  real, intent(out) :: F0      !< downward heat flux from surface, W/m2
  real, intent(out) :: DFDT    !< derivative of snow surface flux w.r.t. temperature, W/m2/K
  ! other output variables needed in lm4p2 snow step 1
  ! real, intent(out) :: snow_T1 ! temperature of snow top layer [K]
  ! real, intent(out) :: snow_rh ! relative humidity, set to 1 currently in lm4p2
  ! real, intent(out) :: snow_ice ! total solid mass of snowpack (avail. for implicit melt)
  ! real, intent(out) :: snow_liq ! total liquid of top layer only (avail. for implicit freeze)
  ! real, intent(out) :: snow_subl ! sublimation (=1) or liq evap (=0) or maybe fractional?
  ! real, intent(out) :: snow_area ! as in lm4p2, for albedo only

  integer :: N ! shorthand for number of layers
  ! logical, intent(out) :: snow_active ! True if snowpack is there (nlayers >= 0)

  ! MAKE SURE THERE ARE AT LEAST 2 SNOW LAYERS FOR HEAT DIFFUSION
  ! call s%force_split_layers(3)

  if(is_watch_point())then
    write(*,*)'#### glass :: snowpack_step_1b checkpoint 1 [input]'
    __DEBUG2__(t_atm,p_surf)
    write( *, *) "heat flux between soil and snow [positive downward]"
    __DEBUG2__(G0, DGDT)
    ! __DEBUG2__(F0, DFDT)
    write(*,*) "s%swheat", s%swheat
    call s%print()
  endif

  N = s%nlayers
  if (allocated(s%e)) deallocate(s%e)
  if (allocated(s%f)) deallocate(s%f)

  if (N > 0) then
    ! snow_active = .TRUE.
    allocate(s%e(N-1))
    allocate(s%f(N-1))
    call step1(s%snow(1:N), N, dt, G0, DGDT, s%swheat, F0, DFDT, s%e(1:N-1), s%f(1:N-1), t_atm, p_surf)

    ! ! compute additional variables needed in LM4P2 snow_step1
    !   snow_T1 = s%snow(1)%T
    !   snow_liq = s%snow(1)%wl ! set this to zero to turn off implicit melt in lm4p2
    !   snow_rh = 1.0
    !   snow_ice = s%ice() ! set this to zero to turn off implicit melt in lm4p2
    !   snow_subl = 1.0
    !   snow_area = s%area()

  else ! no snow
  ! snow_active = .FALSE.
  !   snow_T1 = TFREEZE ! value not defined - do not pass it?
  !   snow_liq = 0.0
  !   snow_ice = 0.0
  !   snow_rh = 1.0 ! value not defined - do not pass it?
  !   snow_area = 0.0
  !   snow_subl = 0.0
    F0 = G0 ! no snow, same flux passed to the underlying soil/lake/glac layer directly
    DFDT = DGDT ! no snow, same flux passed to the underlying soil/lake/glac layer directly
  endif
  if(is_watch_point())then
    write(*,*)'#### glass :: snowpack_step_1b checkpoint 2 [output]'
    write(*,*) "heat flux between surface and snow [positive downward]"
    __DEBUG2__(F0, DFDT)
  endif
end subroutine snowpack_step_1b

!> \brief Forward elimination pass of tridiagonal solver
!!
!! Given the state of snow layers, time step, expression for the downward ground
!! heat flux to the soil G0 + DGDT*\Delta T,  and external heat sources for every
!! layer (e.g. due to absorption of short-wave radiation by the snow), calculates
!! expression for the downward snow heat flux at the surface F0 + DFDT*\Delta T
!! and back-substitution coefficients e and f
subroutine step1(snow, N, dt, G0, DGDT, S, F0, DFDT, e, f, t_atm, p_surf)
  type(snow_layer_type), intent(in)  :: snow(:) !< state of snow layers
  integer,      intent(in)  :: N       !< number of snow levels
  real,         intent(in)  :: dt      !< time step, s
  real,         intent(in)  :: G0      !< ground heat flux, positive downward, W/m2
  real,         intent(in)  :: DGDT    !< derivative of G0 w.r.t. surface temperature, W/m2/K
  real,         intent(in)  :: S(N)    !< energy source for each layer, W/m2
  real,         intent(out) :: F0      !< downward heat flux from surface, W/m2
  real,         intent(out) :: DFDT    !< derivative of snow surface flux w.r.t. temperature, W/m2/K
  real,         intent(out) :: e(N-1)  !< back-substitution coefficient
  real,         intent(out) :: f(N-1)  !< back-substitution coefficient
  real,         intent(in)  :: t_atm, p_surf !< atm temp [K] and pressure [Pa]

  integer :: k
  real :: lambda(N) ! heat conductance of snow layer, W/m/K
  real :: cond(N)   ! heat conductance between layers k and k+1, W/m2/K
  real :: H0(N)     ! downward heat flux from layer k to k+1, W/m2
  real :: gamma     ! coefficient for upward pass
  real :: Q, DQDT   ! full flux and its derivative
  real :: resist

  ! real, dimension(N) :: aaa
  ! real, dimension(N-1) :: ccc
  ! real :: denom, dt_e, bbb
  ! integer l

  do k = 1, N
     ! heat conductance of half-layer (from center to interface)
    !  lambda(k) = snow(k)%hCon()
     lambda(k) = snow(k)%hCon(t_atm, p_surf)
  enddo

  do k = 1, N-1
     ! conductance between layers k and k+1
     ! conductance used to be set to constant 0.3 in lm4p2
    !  cond(k) = 2.0/(snow(k)%dz/lambda(k)+snow(k+1)%dz/lambda(k+1))
    !  if (k==1) then
    !    resist = snow(k)%dz/lambda(k) + snow(k+1)%dz/lambda(k+1)/2.0
    !  else if (k==N-1) then
    !    resist = snow(k)%dz/lambda(k)/2.0 + snow(k+1)%dz/lambda(k+1)
    !  else
    !    resist = snow(k)%dz/lambda(k)/2.0 + snow(k+1)%dz/lambda(k+1)/2.0
    !  endif
     resist = snow(k)%dz/lambda(k)/2.0 + snow(k+1)%dz/lambda(k+1)/2.0
     cond(k) = 1.0/resist
     ! explicit estimate of downward heat flux from layer k to k+1
     H0(k)   = cond(k) * (snow(k)%t - snow(k+1)%t)
  enddo

  Q = G0 ; DQDT = DGDT
  do k = N, 2, -1
     gamma  = 1.0/(snow(k)%hCap()/dt+cond(k-1)+DQDT)
    !  gamma  = 1.0/(snow(k)%hCap()/dt-cond(k-1)+DQDT)
     e(k-1) = gamma*cond(k-1)
    !  f(k-1) = gamma*(H0(k-1) - Q)
     f(k-1) = gamma*(H0(k-1) + S(k) - Q)
!     write(*,'(i2.2, 99(x,a,g12.4))')k, 'Q', Q, 'DQDT', DQDT
     ! for the next level
     Q      = H0(k-1)   - f(k-1)*cond(k-1) ! full flux downward, W/m2
     DQDT   = cond(k-1) - e(k-1)*cond(k-1) ! its full derivative wrt T(K)
  enddo
  F0   = Q - S(1)
  ! F0   = Q
  DFDT = snow(1)%hCap()/dt + DQDT
!  write(*,'(100(x,a,g12.4))')'F0', F0, 'DFDT', DFDT

end subroutine step1

!> \brief Back-substitution pass of tridiagonal solver
!!
!! Given snow layers, back-substitution coefficients calculated in step1, and surface
!! temperature tendency, updates temperatures of snow layers
subroutine snowpack_step_2(snowpack, dT0, dTN)
  class(snowpack_t), intent(inout) :: snowpack
  real,    intent(in)  :: dT0        !< first layer temperature change
  real,    intent(out) :: dTN        !< last layer temperature change, for ground

  real :: dT
  integer :: k
  if (snowpack%nlayers > 0) then
    dT = dT0
    do k = 1, snowpack%nlayers
      snowpack%snow(k)%T = snowpack%snow(k)%T + dT
      if (k<snowpack%nlayers) dT = snowpack%e(k)*dT + snowpack%f(k)
    enddo
    dTN = dT
  else
    dTN = dT0 ! no snow
  endif
end subroutine snowpack_step_2

!> \brief insert layer into snowpack at a given index k
!! the properties of the layer i remain unchanged, but duplicate is inserted as i+1
subroutine snowpack_duplicate_layer(s,i)
  class(snowpack_t), intent(inout) :: s
  integer,           intent(in)    :: i ! index to insert new element to

  integer :: k
  type(snow_layer_type), allocatable :: snow(:)

  ! real heat0, heat1
  ! DUPLICATE LAYER - CHECK ENERGY CONS
  ! heat0 = s%heat()

  ! reallocate array of snow layers as necessary
  if (size(s%snow) <= s%nlayers) then
     ! number of new snow layers is larger than 1 to minimize the number of
     ! allocations
     allocate(snow(int(s%nlayers*1.5+3)))
     snow(1:s%nlayers) = s%snow(1:s%nlayers)
     call move_alloc(snow,s%snow)
  endif
  do k = s%nlayers,i,-1
     s%snow(k+1) = s%snow(k)
  enddo
  s%nlayers = s%nlayers +1

  ! CHECK CONS
end subroutine snowpack_duplicate_layer



! compute average properties for a layer thick thickn at the surface of snowpack
! if snow is thinner than that, get average properties over entire snowapck depth
subroutine snowpack_nearsurf_properties(s)

  class(snowpack_t), intent(inout) :: s
  real :: thickn
  real, DIMENSION(s%nlayers) :: bceq_im, bceq_em
  real, DIMENSION(s%nlayers) :: rho_snow
  integer il, it, il_final
  real cumz ! cumulative depth
  real cumm ! cumulative mass
  real cuml ! cumulative mass liquid + soild
  real cumden ! cumulative snow dendricity (multiplied - weighted by mass)
  real cumd ! cumulative optical diameter (multiplied - weighted by mass)
  real cumc_im ! cumulative mass black carbon equivalent impurities - im
  real cumc_em ! cumulative mass black carbon equivalent impurities - em
  real cuma ! cumulative age [days]
  real cums ! cumulative sphericity
  real zleft, layer_frac_in
  real cumhc, cum_T

  ! properties to average:
  ! optical diameter
  ! density [must be computed] [kg m^-3]
  ! LAI concentration [must be computed] [ppm] - see params
  ! BC equivalent concentration as needed for albedo calculations

  ! first compute the properties needed for each layer and not yet available
  if (s%nlayers > 0) then
  do il = 1, s%nlayers
    rho_snow(il) = s%snow(il)%ws / s%snow(il)%dz
    ! Note, this is the equivalent MASS of black carbon - bring it to concentration here [ppm = mug/g]

    ! bceq_im(il) = s%snow(il)%wc_im(1) + s%snow(il)%wc_im(2)*(LAI_abs(2)/LAI_abs(1)) + &
    !                                     s%snow(il)%wc_im(3)*(LAI_abs(3)/LAI_abs(1))
    ! bceq_em(il) = s%snow(il)%wc_em(1) + s%snow(il)%wc_em(2)*(LAI_abs(2)/LAI_abs(1)) + &
    !                                     s%snow(il)%wc_em(3)*(LAI_abs(3)/LAI_abs(1))

    ! added option to selectively turn off one of the specie in the computation
    ! of near surface LAP concentration used in albedo calculations
    bceq_im(il) = 0.0
    bceq_em(il) = 0.0
    if (lap_albedo_include_bc) then
      bceq_im(il) = bceq_im(il) + s%snow(il)%wc_im(1)
      bceq_em(il) = bceq_em(il) + s%snow(il)%wc_em(1)
    endif
    if (lap_albedo_include_md) then
      bceq_im(il) = bceq_im(il) + s%snow(il)%wc_im(2)*(LAI_abs(2)/LAI_abs(1))
      bceq_em(il) = bceq_em(il) + s%snow(il)%wc_em(2)*(LAI_abs(2)/LAI_abs(1))
    endif
    if (lap_albedo_include_om) then
      bceq_im(il) = bceq_im(il) + s%snow(il)%wc_im(3)*(LAI_abs(3)/LAI_abs(1))
      bceq_em(il) = bceq_em(il) + s%snow(il)%wc_em(3)*(LAI_abs(3)/LAI_abs(1))
    endif

    ! bceq_im(il) = bceq_im(il) / (s%snow(il)%ws) ! concentration in ppm = mg/kg
    ! bceq_em(il) = bceq_em(il) / (s%snow(il)%ws) ! concentration in ppm = mg/kg
    bceq_im(il) = bceq_im(il) / (s%snow(il)%ws + s%snow(il)%wl ) ! concentration in ppm = mg/kg
    bceq_em(il) = bceq_em(il) / (s%snow(il)%ws + s%snow(il)%wl ) ! concentration in ppm = mg/kg

  enddo

  ! thickn = 0.03 ! set layer thickness to 3cm as in Vionnet et al., 2012
  thickn = thickness_for_surface_optical_props
  cumz = 0.0
  cumm = 0.0
  cuml = 0.0
  cumd = 0.0
  cumden = 0.0
  cumc_im = 0.0
  cumc_em = 0.0
  cuma = 0.0
  cums = 0.0
  cumhc = 0.0
  cum_T = 0.0

  if (thickn < s%snow(1)%dz) then ! use only first snow layer
    s%nearsurf_bceq_im = bceq_im(1)
    s%nearsurf_bceq_em = bceq_em(1)
    s%nearsurf_rho = rho_snow(1)
    s%nearsurf_optd = s%snow(1)%optd
    s%nearsurf_dendr = s%snow(1)%dendr
    s%nearsurf_age = s%snow(1)%age
    s%nearsurf_sph = s%snow(1)%sph
    s%nearsurf_T = s%snow(1)%T


  else if (thickn >= s%depth()) then ! compute average properties over entire snowpack
    do il=1,s%nlayers
      cumz = cumz + s%snow(il)%dz
      cumm = cumm + s%snow(il)%ws
      cuml = cuml + s%snow(il)%ws + s%snow(il)%wl
      cumhc = cumhc + s%snow(il)%hCap() ! sum layers heat capacity
      cumden = cumden + s%snow(il)%dendr * s%snow(il)%ws
      cumd = cumd + s%snow(il)%optd * s%snow(il)%ws
      cuma = cuma + s%snow(il)%age * s%snow(il)%ws
      cums = cums + s%snow(il)%sph * s%snow(il)%ws
      cumc_im = cumc_im + bceq_im(il) * ( s%snow(il)%ws + s%snow(il)%wl )
      cumc_em = cumc_em + bceq_em(il) * ( s%snow(il)%ws + s%snow(il)%wl )
      cum_T = cum_T + s%snow(il)%T * s%snow(il)%hCap() ! T weighted by hCap
      ! cumc = cumc + bceq(il)
    enddo
    ! now average properties over entire snowpack
    ! s%nearsurf_bceq = cumc/cumm ! concentration     [ppm]
    s%nearsurf_bceq_im = cumc_im/cuml ! concentration [ppm]
    s%nearsurf_bceq_em = cumc_em/cuml ! concentration [ppm]
    s%nearsurf_rho = cumm/cumz ! density [total mass / total depth] [kg m^-3]
    s%nearsurf_optd = cumd/cumm ! mass weighted average optical diameter [m]
    s%nearsurf_dendr = cumden/cumm ! mass weighted average snow dendricity [number in [0,1]]
    s%nearsurf_age = cuma/cumm ! mass weighted snow age [days]
    s%nearsurf_sph = cums/cumm ! mass weighted snow sphericity [number in [0,1]]
    s%nearsurf_T = cum_T/cumhc ! hCap weighted snow temperature [K]


  else ! ! compute average properties over layer <thickn> only
    il = 1 ! start from top layer
    ! do while (cumz < thickn)
    do while (cumz + s%snow(il)%dz < thickn)
      cumhc = cumhc + s%snow(il)%hCap() ! sum layers heat capacity
      cuma = cuma + s%snow(il)%age * s%snow(il)%ws
      cumz = cumz + s%snow(il)%dz
      cumm = cumm + s%snow(il)%ws
      cuml = cuml + s%snow(il)%ws + s%snow(il)%wl
      cumden = cumden + s%snow(il)%dendr * s%snow(il)%ws
      cumd = cumd + s%snow(il)%optd * s%snow(il)%ws
      cums = cums + s%snow(il)%sph * s%snow(il)%ws
      ! cumc = cumc + bceq(il)
      cumc_im = cumc_im + bceq_im(il) * ( s%snow(il)%ws + s%snow(il)%wl  )
      cumc_em = cumc_em + bceq_em(il) * ( s%snow(il)%ws + s%snow(il)%wl  )
      cum_T = cum_T + s%snow(il)%T * s%snow(il)%hCap() ! T weighted by hCap
      il = il + 1
    enddo
    ! then add the missing piece from the last layer, and average results
    ! il_final = il-1
    ! zleft = cumz + s%snow(il)%dz - thickn
    zleft = thickn - cumz
    layer_frac_in = zleft/s%snow(il)%dz
    ! if (zleft < 0) call land_error_message("error snowpack_nearsurf_properties in snowpack_mod: zleft must be > 0!", FATAL)
    cumhc = cumhc + layer_frac_in * s%snow(il)%hCap() ! sum layers heat capacity
    cuma = cuma + layer_frac_in * s%snow(il)%age * s%snow(il)%ws
    cumz = cumz + layer_frac_in * s%snow(il)%dz
    cumm = cumm + layer_frac_in * s%snow(il)%ws
    cuml = cuml + layer_frac_in * ( s%snow(il)%ws + s%snow(il)%wl  )
    cumd = cumd + layer_frac_in * s%snow(il)%optd * s%snow(il)%ws
    cumden = cumden + layer_frac_in * s%snow(il)%dendr * s%snow(il)%ws
    cums = cums + layer_frac_in * s%snow(il)%sph * s%snow(il)%ws
    ! cumc = cumc + layer_frac_in * bceq(il)
    cumc_em = cumc_em + layer_frac_in * bceq_em(il) * (s%snow(il)%ws + s%snow(il)%wl )
    cumc_im = cumc_im + layer_frac_in * bceq_im(il) * (s%snow(il)%ws + s%snow(il)%wl )
    cum_T = cum_T + layer_frac_in * s%snow(il)%T * s%snow(il)%hCap() ! T weighted by hCap

    if ((zleft < 0.0) .or. (zleft > s%snow(il)%dz )) then
      write(*,*) "ERROR COMPUTING NEAR SURFACE PROPERTIES - CHECK LAYERING!"
      ! error stop "ERROR in snowpack_nearsurf_properties in snowpack module: layering error"
      call land_error_message( "ERROR in snowpack_nearsurf_properties in snowpack module: layering error", FATAL)
    endif

    ! now average properties
    ! s%nearsurf_bceq = cumc/cumm ! concentration     [ppm]
    s%nearsurf_bceq_im = cumc_im/cuml ! concentration [ppm]
    s%nearsurf_bceq_em = cumc_em/cuml ! concentration [ppm]
    s%nearsurf_rho = cumm/cumz ! snow density [total mass / total depth] [kg m^-3]
    s%nearsurf_optd = cumd/cumm ! mass weighted average optical diameter [m]
    s%nearsurf_dendr = cumden/cumm ! mass weighted average optical diameter [m]
    s%nearsurf_age = cuma/cumm ! mass weighted average snow age [days]
    s%nearsurf_sph = cums/cumm  ! mass weighted average snow grain sphericity [number in [0,1]]
    s%nearsurf_T = cum_T/cumhc ! hCap weighted snow temperature [K]

  endif

  else ! case of no snow on the ground

    ! s%nearsurf_bceq_im = -9999.9 ! concentration [ppm]
    ! s%nearsurf_bceq_em = -9999.9 ! concentration [ppm]
    ! s%nearsurf_rho = -9999.9  ! snow density [total mass / total depth] [kg m^-3]
    ! s%nearsurf_optd = -9999.9  ! mass weighted average optical diameter [m]
    ! s%nearsurf_dendr = -9999.9  ! mass weighted average snow grain dendricity [number in [0,1]]
    ! s%nearsurf_age = -9999.9  ! mass weighted average snow age [days]
    ! s%nearsurf_sph = -9999.9  ! mass weighted average snow grain sphericity [number in [0,1]]
    ! s%nearsurf_T = -9999.9  ! hCap weighted average snow temperature [K]

    ! updated to zero when no snow because these are saved as diagnostics weighted by snow area fraction.
    s%nearsurf_bceq_im = 0.0 ! concentration [ppm]
    s%nearsurf_bceq_em = 0.0 ! concentration [ppm]
    s%nearsurf_rho = 0.0  ! snow density [total mass / total depth] [kg m^-3]
    s%nearsurf_optd = 0.0  ! mass weighted average optical diameter [m]
    s%nearsurf_dendr = 0.0  ! mass weighted average snow grain dendricity [number in [0,1]]
    s%nearsurf_age = 0.0  ! mass weighted average snow age [days]
    s%nearsurf_sph = 0.0  ! mass weighted average snow grain sphericity [number in [0,1]]
    s%nearsurf_T = 0.0  ! hCap weighted average snow temperature [K]

  endif
  s%nearsurf_bceq_tot = s%nearsurf_bceq_im + s%nearsurf_bceq_em ! concentration [ppm]

end subroutine snowpack_nearsurf_properties


! STANDALONE MODEL VERSION
!> \given net SW radiation, distribute absorption through snow layers
subroutine snowpack_sw_sources(s, swnet_dir, swnet_dif, swdn_ground)
  class(snowpack_t), intent(inout) :: s
  real, intent(out) :: swdn_ground ! sw radiation exiting the snowpack and passed to the ground [W m^-2]
  ! real swdn_ground ! sw radiation exiting the snowpack and passed to the ground [W m^-2]
  real, INTENT(IN), DIMENSION(NBANDS) :: swnet_dir, swnet_dif ! net sw rad to snowpack
  integer il
  real Q_vis_dir, Q_nir_dir, Q_vis_dif, Q_nir_dif
  real zztop_dir, zzbottom_dir, zztop_dif, zzbottom_dif
  real total0, total1, swnet_in_total

  swnet_in_total = swnet_dif(1) + swnet_dif(2) + swnet_dir(1) + swnet_dir(2)

  ! to compute the radiation absorbed by each snow layer, use the expression for the
  ! radiation flux at depth z below surface : Qz = Qsurface * exp(-beta * z)
  ! See e.g., CROCUS papers  -  Brun 1992 and Vionnet et al., 2012
  ! assume exponential distribution of radiation absorbed within snowpack
  ! e-folding depth is not the same in general for VIS and NIR band
  if (s%nlayers > 0) then
    ALLOCATE(s%swheat(s%nlayers))
    zztop_dir = 0.0
    zzbottom_dir = 0.0
    zztop_dif = 0.0
    zzbottom_dif = 0.0
    do il = 1, s%nlayers
      zzbottom_dir = zzbottom_dir + s%snow(il)%dz
      zzbottom_dif = zzbottom_dif + s%snow(il)%dz
      Q_vis_dir = swnet_dir(1) * ( exp(-s%beta_rad(1)*zztop_dir)  - exp(-s%beta_rad(1)*zzbottom_dir))
      Q_nir_dir = swnet_dir(2) * ( exp(-s%beta_rad(2)*zztop_dir)  - exp(-s%beta_rad(2)*zzbottom_dir))
      Q_vis_dif = swnet_dif(1) * ( exp(-s%beta_rad(1)*zztop_dif)  - exp(-s%beta_rad(1)*zzbottom_dif))
      Q_nir_dif = swnet_dif(2) * ( exp(-s%beta_rad(2)*zztop_dif)  - exp(-s%beta_rad(2)*zzbottom_dif))
      s%swheat(il) = Q_vis_dir + Q_nir_dir + Q_vis_dif + Q_nir_dif
      zztop_dir = zztop_dir + s%snow(il)%dz
      zztop_dif = zztop_dif + s%snow(il)%dz
    enddo
    ! whatwever survives is passed to the ground
    ! radiation passed down to the ground [not done for now, all absorbed by snowpack here]
    swdn_ground = swnet_dir(1)*exp(-s%beta_rad(1)*s%depth()) + &
                  swnet_dir(2)*exp(-s%beta_rad(2)*s%depth()) + &
                  swnet_dif(1)*exp(-s%beta_rad(1)*s%depth()) + &
                  swnet_dif(2)*exp(-s%beta_rad(2)*s%depth())

    !  rescale SW absorbed by snow to match total - assume no penetration to underlying soil
    ! if (swnet_in_total-swdn_ground > 0.0) then
    ! do il = 1, s%nlayers
    !   s%swheat(il) = s%swheat(il) * swnet_in_total/(swnet_in_total-swdn_ground)
    ! enddo
    ! swdn_ground = 0.0
    ! endif

    ! check conservation
    total0 = swnet_in_total
    total1 = swdn_ground + sum(s%swheat(:))

    if ( abs(total0-total1) > eps) then
      write(*,*) "swnet_in_total = ", swnet_in_total
      write(*,*) "total0 = ", total0
      write(*,*) "total1 = ", total1
      write(*,*) "swdn_grnd = ", swdn_ground
      call land_error_message( "ERROR in snowpack_nearsurf_properties in snowpack module: Computing SW sources: energy not conserved!", FATAL)
    endif

  endif


end subroutine snowpack_sw_sources


!> \brief Attempt to split one layer to better match optimal layer thickness distribution
subroutine attempt_split_layers(s)
  class(snowpack_t), intent(inout) :: s

  type(dzopt_t) :: dzopt ! optimal vertical discretization
  real :: z      ! current depth, m
  real :: dz_opt ! optimal layer thickness, m
  real :: dz0    ! original layer thickness, m
  real :: dz1    ! thickness of layer that is split off, m
  real :: f1     ! fraction of layer that is split off
  real :: penalty0, penalty1 ! "distances" from optimal vertical discretization
  integer :: k   ! layer iterator
  type(snow_layer_type), allocatable :: snow1(:) ! new snow array

  ! initialize optimal thickness calculations
  ! write(*,*) "snowpack depth = ", s%depth()
  call dzopt%init(s%depth())

  ! write(*,'(/a)') 'Optimal vertical discretization:'
  ! call dzopt%print()

  z = 0; k = 1
  do while ( k <= s%nlayers ) ! while loop because s%nlayers may change inside
     dz0      = s%snow(k)%dz
     dz_opt   = dzopt%dz(z)
     if (dz0 > dz_opt+epsilon(dz_opt)) then ! do nothing if the layer is thinner than optimum for this depth
        ! calculate thickness of the layer that we might split off
        dz1 = dz_opt
        ! calculate distance of the new layer distribution from the optimal
        penalty0 = dzopt%penalty([z,z+dz0])
        penalty1 = dzopt%penalty([z,z+dz1,z+dz0])
!        write(*,*) k, 'dz(k)',dz(k),'dz_opt',dz_opt,'dz1',dz1,'p0',penalty0, 'p1',penalty1
        if (penalty1 < penalty0) then
           call snowpack_duplicate_layer(s,k)
           f1 = dz1/dz0
           s%snow(k)  %dz = dz0*f1 ! dz0
           s%snow(k+1)%dz = dz0*(1-f1)  ! dz1 - dz0
           ! Q: copy the layers or attempt to match gradients?
           ! in the latter case, what are the assumptions about vertical distribution of water/snow?
           ! in any case, matching gradients may result in overcooled water or overheated snow

           !!!!!!! Additional extensive properties to split:
           !!!!!!! because we are defining these as masses [kg/m2], not densities
           s%snow(k)%ws = s%snow(k)%ws * f1
           s%snow(k)%wl = s%snow(k)%wl * f1
           s%snow(k)%wc_im = s%snow(k)%wc_im * f1 ! vector of size N_SNOW_TRACERS
           s%snow(k)%wc_em = s%snow(k)%wc_em * f1 ! vector of size N_SNOW_TRACERS

           s%snow(k+1)%ws = s%snow(k+1)%ws * (1.0 - f1)
           s%snow(k+1)%wl = s%snow(k+1)%wl * (1.0 - f1)
           s%snow(k+1)%wc_im = s%snow(k+1)%wc_im * (1.0 - f1) ! vector of size N_SNOW_TRACERS
           s%snow(k+1)%wc_em = s%snow(k+1)%wc_em * (1.0 - f1) ! vector of size N_SNOW_TRACERS
           !!!!!!!
        endif
     endif
     z = z + s%snow(k)%dz ! not that it may be different from previous thickness dz0
     k = k+1
  enddo
end subroutine attempt_split_layers

!> \brief checks if two layers can be merged, based on their properties
logical function layers_can_be_merged(s1,s2) result(ans)
  class(snow_layer_type), intent(in) :: s1, s2
  ! TODO: add checks of various snow properties
  logical cond_sph ! condition on sphericity
  logical cond_dopt ! condition on optical diameter
!   logical cond_rho ! condition on density ! avoid for now
  logical cond_ceq ! condition on equiv conc of impurities
  logical cond_dens ! condition snow density



  ! arbitrary thresholds for now
  cond_sph = abs(s1%sph - s2%sph) < 0.2
  cond_dopt = abs(s1%optd - s2%optd) < 1E-4
  cond_dens = abs(s1%density() - s2%density()) < 30.0
  ! cond_dopt = abs(s1%optd - s2%optd)/s2%optd < 0.5 ! see characteristic values. 30% diff for now
  ! cond_ceq = abs(sum(s1%wc) - sum(s2%wc))/sum(s2%wc) < 0.5 ! maybe weight by radiative effects???
  ! cond_sph = .true.
  ! cond_dopt = .true.
  cond_ceq = .true. ! no condition on impurities content for now

  if (cond_dopt .and. cond_sph .and. cond_ceq .and. cond_dens) then
    ans = .true.
  else
    ans = .false.
  endif
end function layers_can_be_merged


subroutine add_liquid_to_layer(s, wl2, T2)


  ! compute thermal equilibrium between a snow layer and some additional liquid
  ! and update the properties of the layer
  ! note: no pore storage is enforced, all mass remains in the layer here
  type(snow_layer_type), intent(inout) :: s !< snow layer to add liquid into
  ! real, intent(in) :: ws1, ws2, wl1, wl2, T1, T2
  real :: ws3, wl3, T3, dz3 ! placeholder for result
  real initws1, initwl1
  real, intent (in) :: wl2, T2
  real :: heat, dz, heat2melt, heatleft, heat2freeze, rho_s, delta_solid
  real :: original_total_liquid, original_total_solid, mass, initT1, initT2
  real final_heat
  integer it
  ! real rho_refrozen
  ! rho_refrozen = 300.0 ! kg/m^2 assume this density for liquid re-freezing

  ! energy conservation
  initwl1 = s%wl
  initws1 = s%ws
  initT1 = s%T
  initT2 = T2
  original_total_liquid = s%wl + wl2
  original_total_solid = s%ws
  mass = original_total_liquid + original_total_solid
  ! heat  = (s%ws*CSW+s%wl*CLW)*(s%T-TFREEZE) + (wl2*CLW)*(T2-TFREEZE) + (s%wl + wl2) * HLF
  heat  = s%hCap()*(s%T-TFREEZE) + (wl2*CLW)*(T2-TFREEZE) + (s%wl + wl2) * HLF
  ! heat2melt = s%ws * HLF ! // old energy balance
  heat2melt = (s%ws + s%wl + wl2) * HLF
  heat2freeze = 0.0 ! all heat is computed wrt solid at freezing termperature
  heatleft = heat - heat2melt

  if(is_watch_point()) then
     write(*,*)'#### add_liquid_to_layer ::: input'
     __DEBUG4__(s%dz,s%ws,s%wl,s%T)
     __DEBUG2__(wl2,T2)
  endif

  if (heatleft > 0.0) then
    ! enough energy to melt everything
    ws3 = 0.0
    wl3 = mass
    T3 = (heatleft)/(ws3*CSW + wl3*CLW) + TFREEZE

    if(is_watch_point()) then
       write(*,*)'#### add_liquid_to_layer ::: case warm'
    endif

  else if (heat < heat2freeze) then
    ! resulting temperature will be <= 0.0
    ! all matter will be in soild phase
    wl3 = 0.0
    ws3 = mass
    ! negative heat will determine the negative temperature of single solid phase
    ! T3 = (heatleft)/(ws3*CSW + wl3*CLW) + TFREEZE
    T3 = (heat)/(ws3*CSW + wl3*CLW) + TFREEZE

    if(is_watch_point()) then
       write(*,*)'#### add_liquid_to_layer ::: case cold, heat < heat2freeze '
    endif
  else
    ! intermediate case: heat2freeze < heat < heat2melt
    ! layer will be at freezing temperature
    T3 = TFREEZE
    ! excess heat will be used to melt as much water as possible
    wl3 = heat/HLF
    ws3 = mass - wl3
    if(is_watch_point()) then
       write(*,*)'#### add_liquid_to_layer ::: case intermediate, heat2freeze < heat < heat2melt'
    endif
  endif

  ! now update thickness of the new layer
  ! if there is a net melt, assume the density of remaining soild stays constant
  ! instead, if there is a net freeze (unlikely, but say we add supercooled water)
  ! assign to the net newly formed solid the density of old snow (350 kg/m3), and do weighted average
  delta_solid = ws3 - original_total_solid
  if(is_watch_point()) then
     __DEBUG4__(delta_solid, ws3, original_total_solid,rho_refrozen)
  endif

!   if (.not.(s%dz>0.0)) then
!      call land_error_message('add_liquid_to_layer: s%dz = '//string(s%dz)//' < 0',FATAL)
!   endif
  ! rho_s = s%ws / max(s%dz, 1E-9)
  if (s%dz > 0.0) then
     rho_s = max(s%ws / s%dz, 10.0)
  else
     ! slm: This is an arbitrary choice in case the layer thickness is zero.
     rho_s = rho_refrozen
  endif
  if (delta_solid > 0.0) then ! net freeze
    dz3 = s%ws / rho_s + delta_solid / rho_refrozen ! sum of layer thickness due to original and newly frozen solid
  else ! net melt
    dz3 = ws3 / rho_s ! preserve density of initial solid phase
  endif
  ! dz3 = max(1E-9, dz3)

  if(is_watch_point()) then
     write(*,*)'#### add_liquid_to_layer'
     __DEBUG3__(delta_solid, ws3, original_total_solid)
     __DEBUG3__(dz3, rho_s, s%ws)
     __DEBUG2__(s%dz, rho_refrozen)
     ! write(*,*) "rho_s, s"
  endif
  ! finally assign new values to snow layer structure
  s%ws = ws3
  s%wl = wl3
  s%T = T3
  s%dz = dz3

  ! final_heat =
  ! final_heat  = (s%ws*CSW+s%wl*CLW)*(s%T-TFREEZE) + (s%wl) * HLF
  final_heat  = s%hCap()*(s%T-TFREEZE) + (s%wl) * HLF

  ! if (  abs(final_heat - heat )> eps ) then
  if (  abs(final_heat - heat )> 1E-4 ) then
     write(*,*) "-------add liquid to layer:: energy not conserved!------"
     write(*,*) "heat = ", heat
     write(*,*) "heat2freeze = ", heat2freeze
     write(*,*) "heat2melt = ", heat2melt
     write(*,*) "heatleft = ", heatleft
     write(*,*) "initial layer: ws, wl, T = ", initws1, initwl1, initT1
     write(*,*) "liquid to add:, wl_add, T_add = ", wl2, initT2
     write(*,*) "final values: ws3, wl3, T3:", ws3, wl3, T3
     write(*,*) "Delta heat = ", final_heat - heat
     write(*,*) "heat = ", heat
     write(*,*) "final_heat = ", final_heat
     call land_error_message( "ERROR in add_liquid_to_layer in snwopack module: energy not conserved!", FATAL)
  endif

  ! additionally, one could change grain properties due to freeze or melt : optd, sph
  !  not done  here
end subroutine add_liquid_to_layer

! // TODO merge this function with add liquid  - they are partially a duplicate of each other
subroutine merge_phases(ws1, wl1, T1, ws2, wl2, T2, ws3, wl3, T3)
  ! compute equilibrium between two two-phases mixtures of water and ice
  real, intent(in) :: ws1, ws2, wl1, wl2, T1, T2
  real, intent(out) :: ws3, wl3, T3
  real :: heat, dz, heat2melt, heatleft, heat2freeze
  real :: original_total_liquid, original_total_solid, mass, initT1, initT2
  integer it

  ! energy conservation
  initT1 = T1
  initT2 = T2
  original_total_liquid = wl1 + wl2
  original_total_solid = ws1 + ws2
  mass = original_total_liquid + original_total_solid
  heat  = (ws1*CSW+wl1*CLW)*(T1-TFREEZE) + (ws2*CSW+wl2*CLW)*(T2-TFREEZE) + (wl1 + wl2) * HLF
  ! heat2melt = (ws1 + ws2) * HLF
  heat2melt = (ws1 + ws2 + wl1 + wl2) * HLF
  heat2freeze = 0.0 ! heat is computed wrt solid at freezing termperature
  heatleft = heat - heat2melt
  if (heatleft > 0.0) then
    ! enough energy to melt everything
    ws3 = 0.0
    wl3 = mass
    T3 = (heatleft)/(ws3*CSW + wl3*CLW) + TFREEZE !`all is liquid here, can remove soild...

  else if (heat < heat2freeze) then
    ! resulting temperature will be <= 0.0
    ! all matter will be in soild phase
    wl3 = 0.0
    ws3 = mass
    ! negative heat will determine the negative temperature of single solid phase
    T3 = heat/(ws3*CSW + wl3*CLW) + TFREEZE ! all is solid here, can remove liquid ...
  else
    ! intermediate case: heat2freeze < heat < heat2melt
    ! layer will be at freezing temperature
    T3 = TFREEZE
    ! all heat will be used to melt as much water as possible
    wl3 = heat/HLF
    ws3 = mass - wl3

  endif
end subroutine merge_phases


!> \brief merge s1 into s2
subroutine merge_layers(s1,s2)
  type(snow_layer_type), intent(in)    :: s1 !< snow layer to merge
  type(snow_layer_type), intent(inout) :: s2 !< snow layer to merge into

  ! real :: heat, dz, heat2melt, heatleft, heat2freeze
  ! real :: original_total_liquid, original_total_solid, mass, initT1, initT2
  real ws3, wl3, T3, dz3, old_rho
  integer it

  ! ezdev added properties
  ! weighted average between the two old layers
  ! so this before changing ws
  s2%sph = (s2%ws * s2%sph  + s1%ws * s1%sph)/(s1%ws + s2%ws)
  s2%optd = (s2%ws * s2%optd + s1%ws * s1%optd)/(s1%ws + s2%ws)
  s2%dendr = (s2%ws * s2%dendr + s1%ws * s1%dendr)/(s1%ws + s2%ws)
  s2%age = (s2%ws * s2%age + s1%ws * s1%age)/(s1%ws + s2%ws)
  old_rho = (s1%ws + s2%ws)/(s1%dz + s2%dz)
  do it = 1, N_SNOW_TRACERS ! sum impurities masses
    ! s2%wc(it) = s2%wc(it) + s1%wc(it)
    s2%wc_im(it) = s2%wc_im(it) + s1%wc_im(it)
    s2%wc_em(it) = s2%wc_em(it) + s1%wc_em(it)
  enddo

  call merge_phases(s1%ws, s1%wl, s1%T, s2%ws, s2%wl, s2%T, ws3, wl3, T3)

   ! assign density to new layer
   ! if overall melting, use old densities
   ! if more frozen mass, assign the density of ice to it
   ! sum        dz_old                          dz_new [if any]

   if ( ws3-(s1%ws + s2%ws) > 0.0) then ! overall freezing
    ! dz3 = (s1%ws + s2%ws)/old_rho + max(0.0, ws3-(s1%ws + s2%ws))/rho_ice
    dz3 = (s1%ws + s2%ws)/old_rho + (ws3-(s1%ws + s2%ws))/rho_refrozen
   else ! overall melting, maintain old density
    dz3 = ws3/old_rho
   endif
   ! if overall freezing, assign to the additional mass the density of ice
   s2%ws = ws3
   s2%wl = wl3
   s2%T = T3
   s2%dz = dz3
end subroutine

!> \brief Attempt to merge two of snowpack layers to better match optimal layer thickness
!! distribution.
subroutine attempt_merge_layers(s)
  class(snowpack_t), intent(inout) :: s

  real :: z, z1 ! current depth
  real :: dz_opt
  type(dzopt_t) :: dzopt
  real :: penalty0, penalty1
  integer :: k, k1, i

  if(is_watch_point()) then
     write(*,*) '#### attempt_merge_layers input'
     do i = 1,s%nlayers
        write(*,'(i2.2)', advance='NO') i
        call dpri('dz',s%snow(i)%dz)
        call dpri('sph',s%snow(i)%sph)
        call dpri('optd',s%snow(i)%optd)
        call dpri('density',s%snow(i)%density())
        write(*,*)
     enddo
  endif

  call dzopt%init(s%depth())

  if(is_watch_point()) then
     write(*,*) '#### attempt_merge_layers loop'
  endif
  z = 0; k = 1
  do while (k < s%nlayers) ! while loop because s%nlayers changes inside
     dz_opt = dzopt%dz(z)
     k1 = k+1 ; z1 = z+s%snow(k)%dz ! index and depth for the next step
     if (is_watch_point()) then
        write(*,'(i2.2)', advance='NO') k
        call dpri('dz',s%snow(k)%dz)
        __DEBUG___(z1)
        __DEBUG___(dz_opt)
     endif
     if (s%snow(k)%dz < dz_opt .and. layers_can_be_merged(s%snow(k), s%snow(k+1))) then
        penalty0 = dzopt%penalty([z, z+s%snow(k)%dz, z+s%snow(k)%dz+s%snow(k+1)%dz])
        penalty1 = dzopt%penalty([z,                 z+s%snow(k)%dz+s%snow(k+1)%dz])
        if (is_watch_point()) then
           __DEBUG___(penalty0)
           __DEBUG___(penalty1)
        endif
        if (penalty1 < penalty0) then
           call merge_layers(s%snow(k+1), s%snow(k))
           do i = k+1, s%nlayers-1
              s%snow(i) = s%snow(i+1)
           enddo
           s%nlayers = s%nlayers-1
           k1 = k ; z1 = z ! do next step with the same layer, except with increased thickness
        endif
     endif
     if (is_watch_point()) then
        write(*,*)
     endif
     k = k1; z = z1
  enddo

  if(is_watch_point()) then
     write(*,*) '#### attempt_merge_layers end'
     do i = 1,s%nlayers
        write(*,'(i2.2)', advance='NO') i
        call dpri('dz',s%snow(i)%dz)
        call dpri('sph',s%snow(i)%sph)
        call dpri('optd',s%snow(i)%optd)
        call dpri('density',s%snow(i)%density())
        write(*,*)
     enddo
  endif
end subroutine attempt_merge_layers


subroutine compute_snow_grain_shape(dendr, sph, idxshp)
    ! based on snow dendricity and sphericity, compute snow shape
    ! the index idxshp identifies the shape
    ! idxshp = 1 => SPHERE
    ! idxshp = 1 => SPHEROID
    ! idxshp = 3 => HEXAGONAL
    ! idxshp = 4 => KOCH
    real, intent(in) :: dendr ! snow dendricity
    real, intent(in) :: sph ! snow sphericity
    integer, intent(out) :: idxshp
    if (dendr > 0.5) then
        idxshp = 4 ! KOCH SNOWFLAKE [Case of fresh dendritic snow]
    else if (sph > 0.8) then
        idxshp = 1 ! SPHERE
    else if (sph < 0.2) then
        idxshp = 3 ! HEXAGONAL [Not spherical at all..]
    else ! remaining case of snow not very dendritic, and of intermediate sphericity
        idxshp = 2 ! SPHEROID
    endif
end subroutine compute_snow_grain_shape


end module snowpack_mod
