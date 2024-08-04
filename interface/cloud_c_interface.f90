module c_interface
   use iso_c_binding

   use config, only: init_config
   use cloud_model, only: model
   use dimensions, only: set_dimensions
   use constants, only: set_constants

contains
   subroutine set_dimensions_python(dx1, nz1)
      ! bind(C, name="configure_constants")
      implicit none
      integer(c_int), intent(in) :: dx1, nz1
      call set_dimensions(dx1, nz1)
   end subroutine set_dimensions_python

   subroutine set_constants_python(G_in, Rd_in, Rv_in, Kapa_in, T0_in, P00_in,&
      Lvl0_in, Lsl0_in, Vis0_in, rhogra_in, Av0_in, Vtnie0_in, Tlvl_in,&
      Tlsl_in, Tlvs_in)
      implicit none
      real(c_float) :: G_in, Rd_in, Rv_in, Kapa_in, T0_in, P00_in, Lvl0_in,&
         Lsl0_in, Vis0_in, rhogra_in, Av0_in, Vtnie0_in
      real(c_float), intent(in) :: Tlvl_in(:), Tlsl_in(:), Tlvs_in(:)

      call set_constants(G_in, Rd_in, Rv_in, Kapa_in, T0_in, P00_in, Lvl0_in,&
         Lsl0_in, Vis0_in, rhogra_in, Av0_in, Vtnie0_in, Tlvl_in, Tlsl_in, Tlvs_in)
   end subroutine set_constants_python

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
