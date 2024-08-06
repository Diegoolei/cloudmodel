module turbvar1
   !! Variables usadas en turbulencia
   use dimensions
   real, allocatable :: D(:,:,:,:,:)
contains
   subroutine allocate_turbvar1
      !! Allocate the turbulence variables
      allocate(D(3, 3, -3:nx1 + 2, -3:nx1 + 2, 3))
   end subroutine allocate_turbvar1

   subroutine deallocate_turbvar1
      !! Deallocate the turbulence variables
      deallocate(D)
   end subroutine deallocate_turbvar1
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
   real, allocatable, dimension(:,:) :: advaer1, advaer2, advgot1, advgot2, advllu1, &
      advllu2, advcri1, advcri2, advnie1, advnie2, advgra1, advgra2, advvap1, advvap2
contains
   subroutine allocate_advecs
      !! Allocate the advection terms
      allocate(advaer1(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advaer2(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advgot1(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advgot2(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advllu1(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advllu2(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advcri1(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advcri2(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advnie1(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advnie2(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advgra1(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advgra2(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advvap1(-2:nx1 + 2, -2:nx1 + 2))
      allocate(advvap2(-2:nx1 + 2, -2:nx1 + 2))
   end subroutine allocate_advecs

   subroutine deallocate_advecs
      !! Deallocate the advection terms
      deallocate(advaer1)
      deallocate(advaer2)
      deallocate(advgot1)
      deallocate(advgot2)
      deallocate(advllu1)
      deallocate(advllu2)
      deallocate(advcri1)
      deallocate(advcri2)
      deallocate(advnie1)
      deallocate(advnie2)
      deallocate(advgra1)
      deallocate(advgra2)
      deallocate(advvap1)
      deallocate(advvap2)
   end subroutine deallocate_advecs
end module advecs

module initial_z_state
   !! Non perturbed quantities
   use dimensions
   use constants, only: mod_nz1
   implicit none
   !!!
   !real, allocatable :: temperature_z_initial(:) !! Non perturbed temperature z initial
   !real, allocatable :: theta_z_initial(:) !! Non perturbed theta z initial
   !real, allocatable :: Pres00(:) !! Non perturbed pressure z initial
   !real, allocatable :: Presi0(:) !! Non perturbed pressure z initial
   !real, allocatable :: u_z_initial(:) !! Non perturbed u z initial
   !real, allocatable :: v_z_initial(:) !! Non perturbed v z initial
   !!!
   real :: temperature_z_initial(-3:mod_nz1 + 3) !! Non perturbed temperature z initial
   real :: theta_z_initial(-3:mod_nz1 + 3) !! Non perturbed theta z initial
   real :: Pres00(-3:mod_nz1 + 3) !! Non perturbed pressure z initial
   real :: Presi0(-3:mod_nz1 + 3) !! Non perturbed pressure z initial
   real :: u_z_initial(-3:mod_nz1 + 3) !! Non perturbed u z initial
   real :: v_z_initial(-3:mod_nz1 + 3) !! Non perturbed v z initial
   real :: air_density_z_initial(-3:mod_nz1 + 3) !! Non perturbed air density z initial
   real :: aerosol_z_initial(-3:mod_nz1 + 3) !! Non perturbed aerosol z initial
   real :: vapor_z_initial(-3:mod_nz1 + 3) !! Non perturbed vapor z initial

   real, allocatable :: cc2(:) !! Non perturbed cc2 z initial

   real, allocatable :: vapor_z_relative(:) !! Non perturbed Relative vapor z initial
   real, allocatable :: aerosol_z_relative(:) !! Non perturbed Relative aerosol z initial
contains
   subroutine set_initial_z_state(temperature_z_initial_in, u_z_initial_in,&
      v_z_initial_in, Presi0_in, air_density_z_initial_in, aerosol_z_initial_in,&
      vapor_z_initial_in)
      real, intent(in), dimension(:) :: temperature_z_initial_in, u_z_initial_in,&
         v_z_initial_in, Presi0_in, air_density_z_initial_in, aerosol_z_initial_in,&
         vapor_z_initial_in

      temperature_z_initial = temperature_z_initial_in
      Presi0 = Presi0_in
      u_z_initial = u_z_initial_in
      v_z_initial = v_z_initial_in

      air_density_z_initial = air_density_z_initial_in
      aerosol_z_initial = aerosol_z_initial_in
      vapor_z_initial = vapor_z_initial_in
   end subroutine set_initial_z_state

   subroutine allocate_initial_z_state
      !! Allocate the initial z state variables
      !!!!
      !allocate(temperature_z_initial(-3:nz1 + 3))
      !allocate(theta_z_initial(-3:nz1 + 3))
      !allocate(Pres00(-3:nz1 + 3))
      !allocate(Presi0(-3:nz1 + 3))
      !allocate(u_z_initial(-3:nz1 + 3))
      !allocate(v_z_initial(-3:nz1 + 3))
      !!!
      allocate(cc2(-3:nz1 + 3))
      !allocate(air_density_z_initial(-3:nz1 + 3))
      !allocate(aerosol_z_initial(-3:nz1 + 3))
      !allocate(vapor_z_initial(-3:nz1 + 3))

      allocate(vapor_z_relative(nz1))
      allocate(aerosol_z_relative(nz1))
   end subroutine allocate_initial_z_state

   subroutine deallocate_initial_z_state
      !! Deallocate the initial z state variables
      !!!
      !deallocate(temperature_z_initial)
      !deallocate(theta_z_initial)
      !deallocate(Pres00)
      !deallocate(Presi0)
      !deallocate(u_z_initial)
      !deallocate(v_z_initial)
      !!!
      deallocate(cc2)
      !deallocate(air_density_z_initial)
      !deallocate(aerosol_z_initial)
      !deallocate(vapor_z_initial)

      deallocate(vapor_z_relative)
      deallocate(aerosol_z_relative)
   end subroutine deallocate_initial_z_state
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
   real, allocatable, dimension(:, :, :) :: U3, V3, W3, Pres3
contains
   subroutine allocate_p3v3
      !! Allocate the p3v3 variables
      allocate(U3(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(V3(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(W3(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(Pres3(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
   end subroutine allocate_p3v3

   subroutine deallocate_p3v3
      !! Deallocate the p3v3 variables
      deallocate(U3)
      deallocate(V3)
      deallocate(W3)
      deallocate(Pres3)
   end subroutine deallocate_p3v3
end module p3v3

!!TODO: Usar para probar allocatable
module filtro01
   use dimensions
   real, allocatable :: varia2(:, :, :)
   real varx, vary, varz, fact
   integer i, j, k
contains
   subroutine allocate_filtro01
      !! Allocate the filtro01 variables
      allocate(varia2(-5:nx1 + 5, -5:nx1 + 5, -4:nz1 + 4))
   end subroutine allocate_filtro01

   subroutine deallocate_filtro01
      !! Deallocate the filtro01 variables
      deallocate(varia2)
   end subroutine deallocate_filtro01
end module filtro01

module sv_inhomogeneous_velocities_and_speed_pressure
   use dimensions
   real, allocatable, dimension(:,:,:) :: fu, fv, fw, fp
contains
   subroutine allocate_sv_inhomogeneous_velocities_and_speed_pressure()
      !! Allocate the sv_inhomogeneous_velocities_and_speed_pressure variables
      allocate(fu(-1:nx1 + 2, -1:nx1 + 2, -1:nz1 + 2))
      allocate(fv(-1:nx1 + 2, -1:nx1 + 2, -1:nz1 + 2))
      allocate(fw(-1:nx1 + 2, -1:nx1 + 2, -1:nz1 + 2))
      allocate(fp(-1:nx1 + 2, -1:nx1 + 2, -1:nz1 + 2))
   end subroutine allocate_sv_inhomogeneous_velocities_and_speed_pressure

   subroutine deallocate_sv_inhomogeneous_velocities_and_speed_pressure()
      !! Deallocate the sv_inhomogeneous_velocities_and_speed_pressure variables
      deallocate(fu)
      deallocate(fv)
      deallocate(fw)
      deallocate(fp)
   end subroutine deallocate_sv_inhomogeneous_velocities_and_speed_pressure
end module sv_inhomogeneous_velocities_and_speed_pressure

module microphysics_perturbation
   !! Module declaration for the microphysics_perturbation module.
   !! The microphysics_perturbation module defines variables related to cloud microphysics.
   use dimensions
   use constants, only: mod_nz1
   real, allocatable :: vapor_base(:, :, :) !! Vapor variable
   real, allocatable :: vapor_new(:, :, :) !! Vapor variable

   real, allocatable :: drop_base(:, :, :) !! Drop variable
   real, allocatable :: drop_new(:, :, :) !! Drop variable

   real, allocatable :: aerosol_base(:, :, :) !! Spray variables
   real, allocatable :: aerosol_new(:, :, :) !! Spray variables

   real, allocatable :: rain_base(:, :, :) !! Liquid cloud variables
   real, allocatable :: rain_new(:, :, :) !! Liquid cloud variables

   real, allocatable :: crystal_base(:, :, :)!! Crystal variables
   real, allocatable :: crystal_new(:, :, :) !! Crystal variables

   real, allocatable :: snow_base(:, :, :) !! Snow variables
   real, allocatable :: snow_new(:, :, :) !! Snow variables

   real, allocatable :: hail_base(:, :, :) !! Graupel (snow pellets) variables
   real, allocatable :: hail_new(:, :, :) !! Graupel (snow pellets) variables

   real :: Av(-3:2*mod_nz1 + 5) !! Other related variables,
   real :: Vtnie(-3:2*mod_nz1 + 5) !! Other related variables,
   real :: Vtgra0(-3:2*mod_nz1 + 5) !! Other related variables
contains
   subroutine set_microphysics_perturbation(Av_in, Vtnie_in, Vtgra0_in)
      real, intent(in), dimension(:) :: Av_in, Vtnie_in, Vtgra0_in
      Vtgra0 = Vtgra0_in
      Av = Av_in
      Vtnie = Vtnie_in
   end subroutine set_microphysics_perturbation

   subroutine allocate_microphysics_perturbation()
      !! Allocate the microphysics_perturbation variables
      allocate(vapor_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(vapor_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(drop_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(drop_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(aerosol_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(aerosol_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(rain_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(rain_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(crystal_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(crystal_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(snow_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(snow_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(hail_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(hail_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      !allocate(Av(-3:2*nz1 + 5))
      !allocate(Vtnie(-3:2*nz1 + 5))
      !allocate(Vtgra0(-3:2*nz1 + 5))
   end subroutine allocate_microphysics_perturbation

   subroutine deallocate_microphysics_perturbation()
      !! Deallocate the microphysics_perturbation variables
      deallocate(vapor_base)
      deallocate(vapor_new)

      deallocate(drop_base)
      deallocate(drop_new)

      deallocate(aerosol_base)
      deallocate(aerosol_new)

      deallocate(rain_base)
      deallocate(rain_new)

      deallocate(crystal_base)
      deallocate(crystal_new)

      deallocate(snow_base)
      deallocate(snow_new)

      deallocate(hail_base)
      deallocate(hail_new)

      !deallocate(Av)
      !deallocate(Vtnie)
      !deallocate(Vtgra0)
   end subroutine deallocate_microphysics_perturbation
end module microphysics_perturbation

module dinamic_var_perturbation
   !! Defines perturbation-related variables for numerical simulations.
   !! These variables are used in the context of dynamic perturbations.
   use dimensions

   real, allocatable :: u_perturbed_base(:, :, :) !! Original velocity component
   real, allocatable :: v_perturbed_base(:, :, :) !! Original velocity component
   real, allocatable :: w_perturbed_base(:, :, :) !! Original velocity component

   real, allocatable :: u_perturbed_new(:, :, :) !! Perturbed velocity component
   real, allocatable :: v_perturbed_new(:, :, :) !! Perturbed velocity component
   real, allocatable :: w_perturbed_new(:, :, :) !! Perturbed velocity component

   real, allocatable :: pressure_base(:, :, :) !! Base perturbed pressure
   real, allocatable :: pressure_new(:, :, :) !! New perturbed pressure

   real, allocatable :: temperature(:, :, :) !! Ambient temperature

   real, allocatable :: theta_base(:, :, :) !! Potential temperature
   real, allocatable :: theta_new(:, :, :) !! Potential temperature

   real, allocatable :: heat_force(:, :, :) !! Heat force
contains
   subroutine allocate_dinamic_var_perturbation()
      !! Allocate the dinamic_var_perturbation variables
      allocate(u_perturbed_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(v_perturbed_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(w_perturbed_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(u_perturbed_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(v_perturbed_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(w_perturbed_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(pressure_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(pressure_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(temperature(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(theta_base(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
      allocate(theta_new(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))

      allocate(heat_force(-3:nx1 + 3, -3:nx1 + 3, -2:nz1 + 2))
   end subroutine allocate_dinamic_var_perturbation

   subroutine deallocate_dinamic_var_perturbation()
      !! Deallocate the dinamic_var_perturbation variables
      deallocate(u_perturbed_base)
      deallocate(v_perturbed_base)
      deallocate(w_perturbed_base)

      deallocate(u_perturbed_new)
      deallocate(v_perturbed_new)
      deallocate(w_perturbed_new)
      deallocate(pressure_base)
      deallocate(pressure_new)
      deallocate(temperature)
      deallocate(theta_base)
      deallocate(theta_new)
      deallocate(heat_force)
   end subroutine deallocate_dinamic_var_perturbation
end module dinamic_var_perturbation


module memory_managment
   use turbvar1, only: allocate_turbvar1, deallocate_turbvar1
   use advecs, only: allocate_advecs, deallocate_advecs
   use initial_z_state, only: allocate_initial_z_state, deallocate_initial_z_state
   use p3v3, only: allocate_p3v3, deallocate_p3v3
   use filtro01, only: allocate_filtro01, deallocate_filtro01
   use sv_inhomogeneous_velocities_and_speed_pressure, only: allocate_sv_inhomogeneous_velocities_and_speed_pressure, &
      deallocate_sv_inhomogeneous_velocities_and_speed_pressure
   use microphysics_perturbation, only: allocate_microphysics_perturbation, deallocate_microphysics_perturbation
   use dinamic_var_perturbation, only: allocate_dinamic_var_perturbation, deallocate_dinamic_var_perturbation

contains
   subroutine allocate_model()
      !! Allocate the model variables
      call allocate_initial_z_state()
      call allocate_turbvar1()
      call allocate_advecs()
      call allocate_p3v3()
      call allocate_filtro01()
      call allocate_sv_inhomogeneous_velocities_and_speed_pressure()
      call allocate_microphysics_perturbation()
      call allocate_dinamic_var_perturbation()
   end subroutine allocate_model

   subroutine deallocate_model()
      !! Deallocate the model variables
      call deallocate_turbvar1()
      call deallocate_advecs()
      call deallocate_p3v3()
      call deallocate_filtro01()
      call deallocate_sv_inhomogeneous_velocities_and_speed_pressure()
      call deallocate_microphysics_perturbation()
      call deallocate_dinamic_var_perturbation()
      call deallocate_initial_z_state()
   end subroutine deallocate_model
end module memory_managment
