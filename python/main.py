from cloudmodel.cloud_read import FileStyle, CloudModel


def main():
    CloudModel(
        simulation_time_minutes=45,
        save_time_minutes=3,
        statistic_time_minutes=3,
        bacup_time_minutes=3,
    )
    
    iniciales = FileStyle(
        chosen_file="Inis",
        output_data_path="Data/new_code/",
        cmp_output_data_path="outputdata1/",
        img_path="img/new_code/",
        txt_path="txt/",
        cmp_txt_path="txt1/",
        vid_path="vid/",
        img_option="Contour",
        folder_handle="Delete",
    )
    nubes = FileStyle(
        chosen_file="Nube",
        output_data_path="Data/new_code/",
        cmp_output_data_path="outputdata1/",
        img_path="img/new_code/",
        txt_path="txt/",
        cmp_txt_path="txt1/",
        vid_path="vid/",
        img_option="Contour",
        folder_handle="Delete",
    )
    # nubes.check_path(nubes.output_data_path)
    iniciales.parse_status_img()
    nubes.parse_status_img()


main()
