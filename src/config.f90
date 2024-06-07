module config
    character(len=:), allocatable :: output_directory
    real :: sim_time_minutes, save_lapse_minutes
 
 contains
    subroutine init_config()
       sim_time_minutes = 45.
       save_lapse_minutes = 3.
       output_directory ="Data/new_code/"
       call create_directory(output_directory)
    end subroutine init_config
 
    subroutine create_directory(directory)
       character(len=*), intent(in) :: directory
       call system("mkdir " // directory)
    end subroutine create_directory
 end module config
 