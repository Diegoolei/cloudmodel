
!> Description: Defines parameters related to grid dimensions and intervals.
!!              These parameters are used in numerical simulations.
module dimensions
   !> Number of points in the x-direction
   integer, parameter :: nx1 = 50
   !> Number of points in the z-direction (first grid)
   integer, parameter :: nz1 = 45
   !> Number of points in the z-direction (possibly a different grid)
   integer, parameter :: nz = 64
   !> Spatial interval (grid spacing) in the x-direction
   real, parameter :: dx1 = 300.0
   !> Time intervals for various purposes
   real, parameter :: dt1 = 2.0, dt2 = 1.0, dt3 = 0.2
 end module dimensions

!> Module: constants
!! Description: Defines mathematical and physical constants relevant to precipitation.
!!              These constants are used in numerical simulations and atmospheric modeling.
 module constants
   !> Mathematical constants
   real, parameter :: pi = 3.1415926
   !> Constants related to gamma functions
   real, parameter :: gam1p8 = 0.9134, gam2p8 = 1.6765, gam3p8 = 4.6941742, gam4p8 = 17.837862
   !> Acceleration due to gravity (m/s^2)
   real, parameter :: G = 9.8
   !> Specific gas constants (dry air and water vapor)
   real, parameter :: Rd = 287.04, Rv = 461.05
   !> Ratio of specific heats (dry air)
   real, parameter :: Kapa = 0.2857
   !> Ratio of molecular weights (water vapor to dry air)
   real, parameter :: Eps = 0.622646
   !> Reference temperature (Kelvin)
   real, parameter :: T0 = 273.15
   !> Reference pressure (Pa)
   real, parameter :: P00 = 101300
   !> Kinematic viscosity of air (m^2/s)
   real, parameter :: Kair = 2.40e-2
   !> Specific heat capacities (liquid water, water vapor, ice)
   real, parameter :: Cwl = 4218., Cwv = 1839., Cwi = 2106.
   !> Latent heats (liquid water, sublimation, vaporization)
   real, parameter :: Lvl0 = 2.500e6, Lsl0 = 79.7, Lvs0 = 2.835e6
   !> Specific heat capacities (constant pressure and constant volume)
   real, parameter :: Cp = 1003., Cv = 716.
   !> Saturation vapor pressures (liquid water and ice)
   real, parameter :: elvs0 = 610.78, esvs0 = 610.918
   !> Diffusivity of water vapor in air (m^2/s)
   real, parameter :: Dv0 = 2.11e-5
   !> Kinematic viscosity of air (m^2/s)
   real, parameter :: Vis0 = 1.718e-5
   !> Density of water (kg/m^3)
   real, parameter :: rhow = 1000.
   !> Densities of various substances (kg/m^3)
   real, parameter :: rhocri = 900., rhonie = 100., rhogra = 500.
   !> Number concentrations (cloud droplets, raindrops, ice crystals, graupel)
   real, parameter :: N0got = 2.9e24, N0llu = 4912189., N0nie = 1.66e5, N0gra = 310.
   !> Terminal fall velocity of ice crystals (m/s)
   real, parameter :: Av0 = 1455.
   !> Threshold velocity for ice crystal nucleation (m/s)
   real, parameter :: Vtnie0 = 0.5
   !> Efficiency factors for collision and coalescence
   real, parameter :: Efcol = 0.8, Efcolgn = 0.7
 
   !> Temperature profiles (for specific ranges)
   real, dimension(210:320) :: Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs, Eautcn, Eacrcn
 end module constants

