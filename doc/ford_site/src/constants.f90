
module dimensions
!!Defines parameters related to grid dimensions and intervals.
!!              These parameters are used in numerical simulations.

   integer, parameter :: nx1 = 50 !! Number of points in the x-direction

   integer, parameter :: nz1 = 45 !! Number of points in the z-direction (first grid)

   integer, parameter :: nz = 64 !! Number of points in the z-direction (possibly a different grid)

   real, parameter :: dx1 = 300.0 !! Spatial interval (grid spacing) in the x-direction

   real, parameter :: dt1 = 2.0 !! Time interval for the main simulation
   real, parameter :: dt2 = 1.0 !! Time interval for microphysics
   real, parameter :: dt3 = 0.2 !! Time interval for speed-pressure
end module dimensions

module constants
!! constants
!! Defines mathematical and physical constants relevant to precipitation.
!! These constants are used in numerical simulations and atmospheric modeling.

   real, parameter :: pi = 3.1415926  !! Mathematical constants

   real, parameter :: gam1p8 = 0.9134 !! Constant related to gamma function
   real, parameter :: gam2p8 = 1.6765 !! Constant related to gamma function
   real, parameter :: gam3p8 = 4.6941742 !! Constant related to gamma function
   real, parameter :: gam4p8 = 17.837862 !! Constant related to gamma function

   real, parameter :: G = 9.8 !! Acceleration due to gravity (m/s^2)

   real, parameter :: Rd = 287.04 !! Specific gas constants (dry air and water vapor)
   real, parameter :: Rv = 461.05 !! Specific gas constants (dry air and water vapor)

   real, parameter :: Kapa = 0.2857 !! Ratio of specific heats (dry air)

   real, parameter :: Eps = 0.622646 !! Ratio of molecular weights (water vapor to dry air)

   real, parameter :: T0 = 273.15 !! Reference temperature (Kelvin)

   real, parameter :: P00 = 101300 !! Reference pressure (Pa)

   real, parameter :: Kair = 2.40e-2 !! Kinematic viscosity of air (m^2/s)

   real, parameter :: Cwl = 4218. !! Specific heat capacities (liquid water)
   real, parameter :: Cwv = 1839. !! Specific heat capacities (water vapor)
   real, parameter :: Cwi = 2106. !! Specific heat capacities (ice)

   real, parameter :: Lvl0 = 2.500e6 !! Latent heats (vapor-liquid)
   real, parameter :: Lsl0 = 79.7 !! Latent heats (solid-liquid)
   real, parameter :: Lvs0 = 2.835e6  !! Latent heats (vapor-solid)

   real, parameter :: Cp = 1003. !! Specific heat capacities (constant pressure)
   real, parameter :: Cv = 716. !! Specific heat capacities (constant volume)

   real, parameter :: elvs0 = 610.78 !! Saturation vapor pressures (liquid vapor solid)
   real, parameter :: esvs0 = 610.918 !! Saturation vapor pressures (solid vapor solid)

   real, parameter :: Dv0 = 2.11e-5 !! Diffusivity of water vapor in air (m^2/s)

   real, parameter :: Vis0 = 1.718e-5 !! Kinematic viscosity of air (m^2/s)

   real, parameter :: rhow = 1000. !! Density of water (kg/m^3)

   real, parameter :: rhocri = 900. !! Densities of crystals (kg/m^3)
   real, parameter :: rhonie = 100. !! Densities of snow (kg/m^3)
   real, parameter :: rhogra = 500. !! Densities of hail (kg/m^3)

   real, parameter :: N0got = 2.9e24 !! Number concentrations (cloud droplets)
   real, parameter :: N0llu = 4912189. !! Number concentrations (raindrops)
   real, parameter :: N0nie = 1.66e5 !! Number concentrations (ice crystals)
   real, parameter :: N0gra = 310. !! Number concentrations (graupel)

   real, parameter :: Av0 = 1455. !! Terminal fall velocity of ice crystals (m/s)

   real, parameter :: Vtnie0 = 0.5 !! Threshold velocity for ice crystal nucleation (m/s)

   real, parameter :: Efcol = 0.8 !! Efficiency factor for collision
   real, parameter :: Efcolgn = 0.7 !! Efficiency factor for collision (graupel)

   real :: Tvis(210:320) !! Temperature profile (viscous)
   real :: Tlvl(210:320) !! Temperature profile (liquid)
   real :: Tlsl(210:320) !! Temperature profile (solid-liquid)
   real :: Tlvs(210:320) !! Temperature profile (liquid-vapor-solid)
   real :: Telvs(210:320) !! Temperature profile (equilibrium liquid-vapor-solid)
   real :: Tesvs(210:320) !! Temperature profile (equilibrium solid-vapor-solid)
   real :: Eautcn(210:320) !! Equilibrium constants (liquid-vapor)
   real :: Eacrcn(210:320) !! Equilibrium constants (solid-liquid)
end module constants

