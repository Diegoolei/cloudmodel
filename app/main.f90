program main
   use cloud_model
   use config
   use dimensions, only: set_dimensions
   use constants, only: set_constants

   implicit none
   call init_config(sim_time=45., save_lapse=3., statistic_time=3., &
                    backup_time=3., restore_from_backup=.false., &
                    directory="Data/z_initials_param/")
                    
   call set_dimensions(dx1_in=300, nz1_in=45)
   !call set_constants(G_in=9.8, Rd_in=287.04, Rv_in=461.05, Kapa_in=0.2857, &
   !                   T0_in=273.15, P00_in=101300., Lvl0_in=2.500e6, Lsl0_in=79.7, &
   !                   Vis0_in=1.718e-5, rhogra_in=500., Av0_in=1455., Vtnie0_in=.5)

   call model()
end program main

