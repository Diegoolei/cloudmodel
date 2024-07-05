module config
   character(len=:), allocatable :: output_directory
   real :: sim_time_minutes, save_lapse_minutes
   real :: statistic_time_minutes, backup_time_minutes
   logical :: restore_backup
contains
   subroutine init_config(sim_time, save_lapse, statistic_time, backup_time,&
      directory)
      !! Initializes the configuration for the simulation.
      character(len=*), intent(in) :: directory
      real, intent(in) :: sim_time !! simulation time in minutes
      real, intent(in) :: save_lapse !! save lapse in minutes
      real, intent(in) :: statistic_time !! statistic time in minutes
      real, intent(in) :: backup_time !! backup time in minutes
      sim_time_minutes = sim_time
      save_lapse_minutes = save_lapse
      statistic_time_minutes = statistic_time
      backup_time_minutes = backup_time
      restore_backup = .false.
      output_directory = directory
   end subroutine init_config
end module config
