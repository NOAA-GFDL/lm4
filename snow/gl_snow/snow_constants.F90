module snow_constants_mod

implicit none

! // TODO for some of these variable, use the one already defined in lm4p2

! real, public, parameter :: PI     = 3.14159265358979
! real, public, parameter :: GGRAV  = 9.80665           !< Terrestrial gravitational constant [m s^-2]
! real, public, parameter :: RDGAS  = 287.04           !< Gas constant for dry air [J/kg/deg]
! real, public, parameter :: KAPPA  = 2.0/7.0  !< RDGAS / CP_AIR [dimensionless]
! real, public, parameter :: CP_AIR = RDGAS/KAPPA              !< Specific heat capacity of dry air at constant pressure [J/kg/deg]
! real, public, parameter :: TFREEZE = 273.15    !< Freezing temperature of fresh water [K]
real, public, parameter :: CLW     = 4218.0    !< specific heat of liquid water [J/kg/K]
real, public, parameter :: CSW     = 2106.0   !< specific heat of ice [J/kg/K] 
real, public, parameter :: CPW     = 1952.0    !< specific heat of water vapor at constant pressure [J/kg/K]
! real, public, parameter :: HLF     = 334000.0 !< Latent heat of fusion [J/kg]
! real, public, parameter :: HLV     = 2.257e6  !< Latent heat of vaporization [J/kg]
! real, public, parameter :: STEFAN  = 5.6734e-8 !< Stefan-Boltzmann constant [W/m^2/deg^4]
! real, public, parameter :: VONKARM = 0.40      !< Von Karman constant [dimensionless]

integer, public, parameter :: NBANDS = 2      !< Number of shortwave bands
integer, public, parameter :: BAND_VIS = 1      !< Index of visible band
integer, public, parameter :: BAND_NIR = 2      !< Index of NIR band
integer, public, parameter :: NTRACERS = 3      !< Number of tracers tracked
integer, public, parameter :: TR_BC = 1      !< Index of black carbon - tracer 1
integer, public, parameter :: TR_MD = 2      !< Index of mineral dust - tracer 2
integer, public, parameter :: TR_OM = 3      !< Index of organic carbon - tracer 3

 !< scavenging coefficients for the tracters                  (BC,  MD,  OM)
real, dimension(NTRACERS), public, parameter :: SCAVENG = (/ 0.2, 0.0, 0.0 /)     
! real, dimension(NTRACERS), public, parameter :: SCAVENG = (/ 0.2, 0.2, 0.2 /)     



real, public, PARAMETER :: rho_refrozen = 300.0 ! refrozen water assumed density [kg / m^3]
real, public, PARAMETER :: rho_ice = 917.0 ! ice density [kg / m^3]
real, public, PARAMETER :: rho_water = 997.0 ! water density [kg / m^3]
real, public, PARAMETER :: eps = 1E-8 ! a small number

! optical properties od BC, MD and OC (respectively) from Veronica's paper
! using the default value for Dust here - see paper for additional values

real, parameter :: thickness_for_surface_optical_props = 0.03 ! [m] 3cm as in Vionnet et al., 2012 - updated to 5cm
real, dimension(NTRACERS), parameter :: LAI_ssa = (/ 0.209, 0.857, 0.963   /) ! single scattering albedo [adim.]
real, dimension(NTRACERS), parameter :: LAI_ext = (/ 9267.0, 474.0, 3289.0 /) ! extinction cross section [m^2 kg^-1]
real, dimension(NTRACERS), parameter :: LAI_sca = (/ 1937.0, 406.0, 3167.0 /) ! scattering cross section [m^2 kg^-1]
real, dimension(NTRACERS), parameter :: LAI_abs = (/ 7330.0, 67.8, 122.0   /) ! absorption cross section [m^2 kg^-1]
! single scattering albedos


end module snow_constants_mod