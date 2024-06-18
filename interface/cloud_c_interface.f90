module c_interface
    use iso_c_binding

    use config, only: init_config
    use cloud_model, only: model

contains

    subroutine run_model_python(sim_time, save_lapse) bind(C, name="run_model")
        real(c_float), intent(in) :: sim_time
        real(c_float), intent(in) :: save_lapse
        call init_config(sim_time=sim_time, save_lapse=save_lapse, directory= "Data/new_code/")
        call model()
    end subroutine
end module
