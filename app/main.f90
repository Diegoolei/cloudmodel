program main
    use cloud_model
    use config
 
    implicit none
    call init_config(sim_time = 45., save_lapse = 3., directory= "Data/new_code/")
    call model()
 end program main
 