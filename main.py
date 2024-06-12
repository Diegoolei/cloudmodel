from cloud_read import FileStyle
import numpy as np
import nubepy.nubepy as nb

nb.c_interface.run_model_python(45,3)
def main():
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
    iniciales.parse_status_img()
    nubes.parse_status_img()

main()