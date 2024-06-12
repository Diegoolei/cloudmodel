module c_interface
    use iso_c_binding

    use config, only: sim_time_minutes, save_lapse_minutes
    use model_runner, only: run_model

contains

    subroutine run_model_python(sim_time, save_time) bind(C, name="run_model")
        real(c_float), intent(in) :: sim_time
        real(c_float), intent(in) :: save_time
        call run_model()
    end subroutine
end module
