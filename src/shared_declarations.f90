module turbvar1
   !! Variables usadas en turbulencia
   use dimensions
   real :: D(3, 3, -3:nx1 + 2, -3:nx1 + 2, 3)
end module turbvar1

module cant01
   !! Cantidades auxiliares, pasos de tiempos, distancias, etc.
   implicit none
   real :: ltt  !! total simulation time
   real :: ltg  !! total save time
   real :: lte !! total statistic save time
   real :: ltb !! total backup time
   real :: ctur, cteturb, dx2, dx8, dx12, AA, ikapa, pro1, &
      pro2, pro3, pro4, cteqgot, cteqllu, cteqnie, cteqgra

   integer :: total_time, lt2, lt3
end module cant01

module turbvar
   !! Variables usadas en turbulencia
   real :: KMM, KM1, KM2, KM3
   real, dimension(3, 3) :: DD
   real, dimension(3) :: D1, D2, D3
end module turbvar

module advecs
   !! Terminos de adveccion
   use dimensions
   real, dimension(-2:nx1 + 2, -2:nx1 + 2) :: advaer1, advaer2, advgot1, advgot2, advllu1, &
      advllu2, advcri1, advcri2, advnie1, advnie2, advgra1, advgra2, advvap1, advvap2
end module advecs

module initial_z_state
   !! Non perturbed quantities
   use dimensions
   real :: temperature_z_initial(-3:nz1 + 3) !! Non perturbed temperature z initial
   real :: theta_z_initial(-3:nz1 + 3) !! Non perturbed theta z initial
   real :: Pres00(-3:nz1 + 3) !! Non perturbed pressure z initial
   real :: Presi0(-3:nz1 + 3) !! Non perturbed pressure z initial
   real :: u_z_initial(-3:nz1 + 3) !! Non perturbed u z initial
   real :: v_z_initial(-3:nz1 + 3) !! Non perturbed v z initial
   real :: cc2(-3:nz1 + 3) !! Non perturbed cc2 z initial
   real :: air_density_z_initial(-3:nz1 + 3) !! Non perturbed air density z initial
   real :: aerosol_z_initial(-3:nz1 + 3) !! Non perturbed aerosol z initial
   real :: vapor_z_initial(-3:nz1 + 3) !! Non perturbed vapor z initial

   real :: vapor_z_relative(nz1) !! Non perturbed Relative vapor z initial
   real :: aerosol_z_relative(nz1) !! Non perturbed Relative aerosol z initial
end module initial_z_state

module lmncri
   !! Posiciones en las cuales Qcri<0
   integer, dimension(2) :: lcri, mcri, ncri
end module lmncri

module lmngot
   !! Posiciones en las cuales Qgot<0
   integer, dimension(2) :: lgot, mgot, ngot
end module lmngot

module lmngra
   !! Posiciones en las cuales Qgra<0
   integer, dimension(2) :: lgra, mgra, ngra
end module lmngra

module lmnllu
   !! Posiciones en las cuales Qllu<0
   integer, dimension(2) :: lllu, mllu, nllu
end module lmnllu

module lmnnie
   !! Posiciones en las cuales Qnie<0
   integer, dimension(2) :: lnie, mnie, nnie
end module lmnnie

!!TODO: Usar para probar allocatable
module p3v3
   !! Velocidades y las presiones
   use dimensions
   real, dimension(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) :: U3, V3, W3, Pres3
end module p3v3

!!TODO: Usar para probar allocatable
module filtro01
   use dimensions
   real varia2(-5:nx1 + 5, -5:nx1 + 5, -4:nz1 + 4)
   real varx, vary, varz, fact
   integer i, j, k
end module filtro01

module sv_inhomogeneous_velocities_and_speed_pressure
   use dimensions
   real, dimension(-1:nx1 + 2, -1:nx1 + 2, -1:nz1 + 2) :: fu, fv, fw, fp
end module sv_inhomogeneous_velocities_and_speed_pressure

module microphysics_perturbation
   !! Module declaration for the microphysics_perturbation module.
   !! The microphysics_perturbation module defines variables related to cloud microphysics.
   use dimensions
   real :: vapor_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Vapor variable
   real :: vapor_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Vapor variable

   real :: drop_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Drop variable
   real :: drop_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Drop variable

   real :: aerosol_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Spray variables
   real :: aerosol_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Spray variables

   real :: rain_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Liquid cloud variables
   real :: rain_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Liquid cloud variables

   real :: crystal_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2)!! Crystal variables
   real :: crystal_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Crystal variables

   real :: snow_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Snow variables
   real :: snow_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Snow variables

   real :: hail_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Graupel (snow pellets) variables
   real :: hail_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Graupel (snow pellets) variables

   real :: Av(-3:2*nz1 + 5) !! Other related variables,
   real :: Vtnie(-3:2*nz1 + 5) !! Other related variables,
   real :: Vtgra0(-3:2*nz1 + 5) !! Other related variables
end module microphysics_perturbation

module dinamic_var_perturbation
   !! Defines perturbation-related variables for numerical simulations.
   !! These variables are used in the context of dynamic perturbations.
   use dimensions

   real :: u_perturbed_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Original velocity component
   real :: v_perturbed_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Original velocity component
   real :: w_perturbed_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Original velocity component

   real :: u_perturbed_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Perturbed velocity component
   real :: v_perturbed_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Perturbed velocity component
   real :: w_perturbed_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Perturbed velocity component

   real :: pressure_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Base perturbed pressure
   real :: pressure_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! New perturbed pressure

   real :: temperature(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Ambient temperature

   real :: theta_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Potential temperature
   real :: theta_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Potential temperature

   real :: heat_force(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2) !! Heat force
end module dinamic_var_perturbation
