module land_constants_mod

use constants_mod, only : rdgas, rvgas, wtmair, dens_h2o, grav, cp_air

implicit none
private

! ==== public interfaces =====================================================
integer, public, parameter :: &
     NBANDS   = 2, & ! number of spectral bands for short-wave radiation calculations
     BAND_VIS = 1, & ! visible radiation (wavelength range?)
     BAND_NIR = 2    ! near infra-red radiation (wavelength range?)

integer, public, parameter :: MAX_SOIL_LEV = 100 ! max number of soil layers (max dimension of arrays)
! litter pool constants
integer, parameter, public :: &
     N_LITTER_POOLS    = 2, &
     LITT_LEAF         = 1, & ! leaf litter
     LITT_CWOOD        = 2    ! coarse wood litter
! names of the litter pools, for i/o
character(16), parameter, public :: &
     l_shortname(N_LITTER_POOLS) = [ 'leaf            ', 'cwood           '  ], & ! for restart field names
     l_longname (N_LITTER_POOLS) = [ 'leaf            ', 'coarse wood     '  ], & ! for long names
     l_diagname (N_LITTER_POOLS) = [ 'lf              ', 'cw              '  ]    ! for diag field names

integer, public, parameter :: N_C_TYPES = 3  ! Carbon chemical species (Cellulose, lignin, microbial products)
integer, public, parameter :: & ! indices of carbon chemical species
     C_FAST = 1, & ! cellulose (fast)
     C_SLOW = 2, & ! lignin (slow)
     C_MIC  = 3    ! microbial products
! names of the carbon types, for i/o
character(len=12), public, parameter :: &
     c_shortname(N_C_TYPES) = [ 'fast        ', 'slow        ', 'deadmic     ' ], & ! for restart field names
     c_longname (N_C_TYPES) = [ 'fast        ', 'slow        ', 'dead microbe' ], & ! for long names
     c_diagname (N_C_TYPES) = [ 'fast        ', 'slow        ', 'dmic        ' ]    ! for diag field names

real, public, parameter :: d622 = rdgas/rvgas
real, public, parameter :: d378 = 1.0-d622
real, public, parameter :: d608 = d378/d622

real, public, parameter :: Rugas = 8.314472 ! universal gas constant, J K-1 mol-1
real, public, parameter :: kBoltz= 1.3807e-23 ! Boltzmann's constant, J K-1 Rugas/Avogadro number

real, public, parameter :: days_per_year = 365.0
real, public, parameter :: seconds_per_year = 86400.0*days_per_year
real, public, parameter :: mol_C = 12.0e-3 ! molar mass of carbon, kg
real, public, parameter :: mol_air = wtmair/1000.0 ! molar mass of air, kg
real, public, parameter :: mol_CO2 = 44.00995e-3 ! molar mass of CO2,kg
real, public, parameter :: mol_h2o = 18.0e-3 ! molar mass of water, kg

real, public, parameter :: dens_ice = 917.0 ! density of ice [kg m-3]

real, public, parameter :: MPa_per_m = dens_h2o*grav*1.0e-6 ! pressure of one meter of water, Mega Pascal

! BRDF parameters from the MODIS brdf/albedo product user guide:
! C. B. Schaaf et al. (2002): First operational BRDF, albedo nadir reflectance products from
! MODIS. Remote Sensing of Environment, 83, No.1-2, 135–148, doi:10.1016/s0034-4257(02)00091-3.
real, public, parameter :: g_iso  = 1.
real, public, parameter :: g_vol  = 0.189184
real, public, parameter :: g_geo  = -1.377622
real, public, parameter :: g0_iso = 1.0
real, public, parameter :: g1_iso = 0.0
real, public, parameter :: g2_iso = 0.0
real, public, parameter :: g0_vol = -0.007574
real, public, parameter :: g1_vol = -0.070987
real, public, parameter :: g2_vol =  0.307588
real, public, parameter :: g0_geo = -1.284909
real, public, parameter :: g1_geo = -0.166314
real, public, parameter :: g2_geo =  0.041840


!real, public, parameter :: kin_visc_air = 1.568e-5 ! kinematic viscosity of air, m2/s

public :: diffusivity_h2o ! (T,p) diffusivity of H2O in air, m2/s
public :: dyn_visc_air    ! (T)   dynamic viscosity of air, kg/(m s)
public :: kin_visc_air    ! (T,p) kinematic viscosity of dry air, m2/s
public :: thermal_diff_air! (T)   thermal diffusivity of air, m2/s
contains

! ---------------------------------------------------------------------------------------
! diffusivity of H2O in air, m2/s
! W. J. Massman (1998): A review of the molecular diffusivities of H2O, CO2, CH4, CO, O3,
! SO2, NH3, N2O, NO, and NO2 in air, O2 and N2 near STP. Atmospheric Environment, 32, No.6,
! 1111–1127, doi:10.1016/s1352-2310(97)00391-9.
real function diffusivity_h2o(T,p) result(D)
    real, intent(in) :: T ! temperature, degK
    real, intent(in) :: p ! pressure, N/m2

    real, parameter :: T0 = 273.15    ! reference temperature, degK
    real, parameter :: p0 = 101325.0  ! reference pressure, N/m2
    D = 0.2178e-4*(p0/p)*(T/T0)**1.81
end function diffusivity_h2o

! ---------------------------------------------------------------------------------------
! dynamic viscosity of air, kg/(m s)
! Smithsonian tables
real function dyn_visc_air(T) result(mu)
    real, intent(in) :: T ! temperature, degK

    real, parameter :: mu0 = 1.8325e-5 ! reference dynamic viscosity, kg/(m s)
    real, parameter :: T0  = 296.6     ! reference temperature, degK
    real, parameter :: C   = 120.0     ! Sutherland model constant, degK

    mu = mu0*(T0+C)/(T+C)*sqrt((T/T0)**3)
end function dyn_visc_air

! ---------------------------------------------------------------------------------------
! kinematic viscosity of dry air, m2/s
! We neglect air specific humidity effect on density here. Even if we do take it into
! account here, the dynamic viscosity calculation neglects it.
real function kin_visc_air(T, p) result(nu)
    real, intent(in) :: T ! temperature, degK
    real, intent(in) :: p ! pressure, N/m2

    nu = dyn_visc_air(T)*rdgas*T/p
end function kin_visc_air

! ---------------------------------------------------------------------------------------
! thermal diffusivity of air, m2/s
! https://en.wikipedia.org/wiki/Prandtl_number
real function thermal_diff_air(T) result(k)
    real, intent(in) :: T ! temperature, degK

    real, parameter :: Pr = 0.71 ! Prandtl number for air
    k = dyn_visc_air(T) / Pr
end function thermal_diff_air


end module
