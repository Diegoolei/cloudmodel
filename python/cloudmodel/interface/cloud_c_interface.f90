module c_interface
   use iso_c_binding

   use config, only: init_config
   use cloud_model, only: model

contains

   subroutine run_model_python(sim_time, save_lapse, statistic_time, backup_time, restore_backup, directory)
      ! bind(C, name="run_model")
      IMPLICIT NONE
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
