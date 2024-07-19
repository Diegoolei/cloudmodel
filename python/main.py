"""main file to run the cloud model."""
from cloudmodel.cloud_read import CloudSimulation


def main():
    """Execute function to run the cloud model."""

    model = CloudSimulation(simulation_time_minutes=3,
        save_time_minutes=1,
        statistic_time_minutes=1,
        bacup_time_minutes=1,
        restore_backup=False,)
    model.run_model()

    cloud = model.run_cloud_analysis()
    cloud.parse_status_img()
    initials = model.run_initial_analysis()
    initials.parse_status_img()

    model.clean_model()
main()
