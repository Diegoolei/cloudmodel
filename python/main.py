from cloudmodel.cloud_read import CloudSimulation


def main():
    z_initials_param = CloudSimulation(
        simulation_time_minutes=45,
        save_time_minutes=3,
        statistic_time_minutes=3,
        bacup_time_minutes=3,
        restore_backup=False,
        directory="Data/z_initials_param",
    )
    z_initials_param.run_model()
    z_initials_param.cloud_analytics.parse_status_img()


main()
