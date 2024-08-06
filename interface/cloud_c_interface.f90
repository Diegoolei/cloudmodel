module c_interface
   use iso_c_binding

   use config, only: init_config
   use cloud_model, only: model
   use dimensions, only: set_dimensions
   use constants, only: set_constants
   use initial_z_state, only: set_initial_z_state
   use microphysics_perturbation, only: set_microphysics_perturbation

contains
   subroutine set_dimensions_python(dx1, nz1)
      ! bind(C, name="configure_constants")
      implicit none
      integer(c_int), intent(in) :: dx1, nz1
      call set_dimensions(dx1, nz1)
   end subroutine set_dimensions_python

   subroutine set_constants_python(G_in, Rd_in, Rv_in, Kapa_in, T0_in, P00_in,&
      Lvl0_in, Lsl0_in, Vis0_in, rhogra_in, Av0_in, Vtnie0_in, Tlvl_in,&
      Tlsl_in, Tlvs_in, Telvs_in, Tesvs_in, Tvis_in, Eautcn_in, Eacrcn_in)
      implicit none
      real(c_float) :: G_in, Rd_in, Rv_in, Kapa_in, T0_in, P00_in, Lvl0_in,&
         Lsl0_in, Vis0_in, rhogra_in, Av0_in, Vtnie0_in
      real(c_float), intent(in), dimension(:) :: Tlvl_in, Tlsl_in, Tlvs_in, &
         Telvs_in, Tesvs_in, Tvis_in, Eautcn_in, Eacrcn_in

      call set_constants(G_in, Rd_in, Rv_in, Kapa_in, T0_in, P00_in, Lvl0_in,&
         Lsl0_in, Vis0_in, rhogra_in, Av0_in, Vtnie0_in, Tlvl_in, Tlsl_in,&
         Tlvs_in, Telvs_in, Tesvs_in, Tvis_in, Eautcn_in, Eacrcn_in)
   end subroutine set_constants_python

   subroutine set_initial_z_state_python(temperature_z_initial_in, u_z_initial_in,&
      v_z_initial_in, Presi0_in, air_density_z_initial_in, aerosol_z_initial_in, &
      vapor_z_initial_in)

      real(c_float), intent(in), dimension(:) :: temperature_z_initial_in, u_z_initial_in,&
         v_z_initial_in, Presi0_in, air_density_z_initial_in, aerosol_z_initial_in,&
         vapor_z_initial_in

      call set_initial_z_state(temperature_z_initial_in, u_z_initial_in,&
         v_z_initial_in, Presi0_in, air_density_z_initial_in, aerosol_z_initial_in, &
         vapor_z_initial_in)

   end subroutine set_initial_z_state_python

   subroutine set_microphysics_perturbation_python(Av_in, Vtnie_in, Vtgra0_in)
      real(c_float), intent(in), dimension(:) :: Av_in, Vtnie_in, Vtgra0_in
      call set_microphysics_perturbation(Av_in, Vtnie_in, Vtgra0_in)
   end subroutine set_microphysics_perturbation_python

   subroutine run_model_python(sim_time, save_lapse, statistic_time, backup_time, restore_backup, directory)
      ! bind(C, name="run_model")
      implicit none
      real(c_float), intent(in) :: sim_time
      real(c_float), intent(in) :: save_lapse
      real(c_float), intent(in) :: statistic_time
      real(c_float), intent(in) :: backup_time
      logical, intent(in) :: restore_backup
      character(256), intent(in) :: directory
      call init_config(sim_time=sim_time, save_lapse=save_lapse, &
         statistic_time=statistic_time, backup_time=backup_time, &
         restore_from_backup=restore_backup, directory=trim(adjustl(directory)))
      call model()
   end subroutine
end module
