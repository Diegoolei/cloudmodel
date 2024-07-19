program main
   use cloud_model
   use config

   implicit none
   call init_config(sim_time=45., save_lapse=3., statistic_time=3., &
                    backup_time=3., restore_from_backup=.false., &
                    directory="Data/new_code/")
   call model()
end program main

