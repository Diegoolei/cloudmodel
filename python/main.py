from cloudmodel.cloud_read import CloudSimulation, ImageStyle


def main():
    z_initials_param = CloudSimulation(
        simulation_time_minutes=51,
        save_time_minutes=2,
        statistic_time_minutes=1 / 3,
        backup_time_minutes=2,
        restore_backup=False,
        directory="Data/demo/",
    )
    z_initials_param.run_model()
    z_initials_param.generate_cut_image("cut_image.png")

    aerosol_base_24 = z_initials_param.cloud_analytics.get_var(
        var="aerosol_base", time=24
    )  # Obtiene la variable aerosol_base del tiempo 24

    cloud = z_initials_param.cloud_analytics
    center_aerosol_base_24 = cloud._get_var_max_value_position(
        aerosol_base_24
    )  # Centro de variable por valor maximo
    var_dataframe = z_initials_param.cloud_analytics.get_var_dataframe(
        aerosol_base_24, center_aerosol_base_24, "z"
    )  # En pantalla muestra el dataframe de la var
    print(var_dataframe)

    # Generate and Save all var images in directory/img
    z_initials_param.cloud_analytics.parse_status_img(
        img_option=ImageStyle.IMAGE.value
    )


main()
