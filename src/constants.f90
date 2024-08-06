
module dimensions
!!Defines parameters related to grid dimensions and intervals.
!!              These parameters are used in numerical simulations.

   integer, parameter :: nx1 = 50 !! Number of points in the x-direction
   integer, parameter :: nz = 64 !! Number of points in the z-direction (possibly a different grid)
   real, parameter :: dt1 = 2.0 !! Time interval for the main simulation
   real, parameter :: dt2 = 1.0 !! Time interval for microphysics
   real, parameter :: dt3 = 0.2 !! Time interval for speed-pressure

   real :: dx1 !! Spatial interval (grid spacing) in the x-direction
   integer :: nz1 !! Number of points in the z-direction (first grid)
contains
   subroutine set_dimensions(dx1_in, nz1_in)
      !! Sets the dimensions and intervals for the grid.
      !! Parameters:
      !! dx1: Spatial interval (grid spacing) in the x-direction
      !! nz1: Number of points in the z-direction (first grid)
      implicit none
      integer, intent(in) :: dx1_in, nz1_in
      dx1 = dx1_in
      nz1 = nz1_in
   end subroutine set_dimensions
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

   real, parameter :: Eps = 0.622646 !! Ratio of molecular weights (water vapor to dry air)

   real, parameter :: Kair = 2.40e-2 !! Kinematic viscosity of air (m^2/s)

   real, parameter :: Cwl = 4218. !! Specific heat capacities (liquid water)
   real, parameter :: Cwv = 1839. !! Specific heat capacities (water vapor)
   real, parameter :: Cwi = 2106. !! Specific heat capacities (ice)

   real, parameter :: Lvs0 = 2.835e6  !! Latent heats (vapor-solid)

   real, parameter :: Cp = 1003. !! Specific heat capacities (constant pressure)
   real, parameter :: Cv = 716. !! Specific heat capacities (constant volume)

   real, parameter :: elvs0 = 610.78 !! Saturation vapor pressures (liquid vapor solid)
   real, parameter :: esvs0 = 610.918 !! Saturation vapor pressures (solid vapor solid)

   real, parameter :: Dv0 = 2.11e-5 !! Diffusivity of water vapor in air (m^2/s)

   real, parameter :: rhow = 1000. !! Density of water (kg/m^3)

   real, parameter :: rhocri = 900. !! Densities of crystals (kg/m^3)
   real, parameter :: rhonie = 100. !! Densities of snow (kg/m^3)

   real, parameter :: N0got = 2.9e24 !! Number concentrations (cloud droplets)
   real, parameter :: N0llu = 4912189. !! Number concentrations (raindrops)
   real, parameter :: N0nie = 1.66e5 !! Number concentrations (ice crystals)
   real, parameter :: N0gra = 310. !! Number concentrations (graupel)

   real, parameter :: Efcol = 0.8 !! Efficiency factor for collision
   real, parameter :: Efcolgn = 0.7 !! Efficiency factor for collision (graupel)

   integer, parameter :: mod_nz1 = 45

   real :: G !! Acceleration due to gravity (m/s^2)
   real :: Rd !! Specific gas constants (dry air and water vapor)
   real :: Rv !! Specific gas constants (dry air and water vapor)
   real :: Kapa !! Ratio of specific heats (dry air)
   real :: T0 !! Reference temperature (Kelvin)
   real :: P00 !! Reference pressure (Pa)
   real :: Lvl0 !! Latent heats (vapor-liquid)
   real :: Lsl0 !! Latent heats (solid-liquid)
   real :: Vis0 !! Kinematic viscosity of air (m^2/s)
   real :: rhogra !! Densities of hail (kg/m^3)
   real :: Av0 !! Terminal fall velocity of ice crystals (m/s)
   real :: Vtnie0 !! Threshold velocity for ice crystal nucleation (m/s)

   real :: Tlvl(210:320) !! Temperature profile (liquid)
   real :: Tlsl(210:320) !! Temperature profile (solid-liquid)
   real :: Tlvs(210:320) !! Temperature profile (liquid-vapor-solid)
   real :: Telvs(210:320) !! Temperature profile (equilibrium liquid-vapor-solid)
   real :: Tesvs(210:320) !! Temperature profile (equilibrium solid-vapor-solid)
   real :: Tvis(210:320) !! Temperature profile (viscous)
   real :: Eautcn(210:320) !! Equilibrium constants (liquid-vapor)
   real :: Eacrcn(210:320) !! Equilibrium constants (solid-liquid)
contains
   subroutine set_constants(G_in, Rd_in, Rv_in, Kapa_in, T0_in, P00_in,&
      Lvl0_in, Lsl0_in, Vis0_in, rhogra_in, Av0_in, Vtnie0_in, Tlvl_in,&
      Tlsl_in, Tlvs_in, Telvs_in, Tesvs_in, Tvis_in, Eautcn_in, Eacrcn_in)
      !! Sets the dimensions and intervals for the grid.
      !! Parameters:
      !! dx1: Spatial interval (grid spacing) in the x-direction
      !! nz1: Number of points in the z-direction (first grid)
      implicit none
      real, intent(in) :: G_in, Rd_in, Rv_in, Kapa_in, T0_in, P00_in,&
         Lvl0_in, Lsl0_in, Vis0_in, rhogra_in, Av0_in, Vtnie0_in
      real, intent(in) :: Tlvl_in(:), Tlsl_in(:), Tlvs_in(:), Telvs_in(:),&
         Tesvs_in(:), Tvis_in(:), Eautcn_in(:), Eacrcn_in(:)
      G = G_in
      Rd = Rd_in
      Rv = Rv_in
      Kapa = Kapa_in
      T0 = T0_in
      P00 = P00_in
      Lvl0 = Lvl0_in
      Lsl0 = Lsl0_in
      Vis0 = Vis0_in
      rhogra = rhogra_in
      Av0 = Av0_in
      Vtnie0 = Vtnie0_in

      Tlvl = Tlvl_in
      Tlsl = Tlsl_in
      Tlvs = Tlvs_in
      Telvs = Telvs_in
      Tesvs = Tesvs_in
      Tvis = Tvis_in
      Eautcn = Eautcn_in
      Eacrcn = Eacrcn_in
   end subroutine set_constants
end module constants

