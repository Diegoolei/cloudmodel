from scipy.io import FortranFile
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from constants import *
import numpy as np
import os
import re  # Add the re import declaration to use regex
from subprocess import Popen, PIPE
from filecmp import cmpfiles


def main():
    compare_binaries = input("Compare binaries? (Y/N): ").upper()
    compare_text_files = input("Compare text files? (Y/N): ").upper()
    generate_cloud_animation = input("Generate cloud animation? (Y/N): ").upper()
    generate_cloud_text_files = input("Generate cloud text files? (Y/N): ").upper()
    generate_cloud_status_image = input("Generate cloud status image? (Y/N): ").upper()

    selected_file = File_style()
    if compare_binaries == "Y":
        cloud_binary_comparison(selected_file)

    if generate_cloud_animation == "Y":
        animate_variable(selected_file)

    if compare_text_files == "Y":
        cloud_text_comparison(selected_file)

    parse_bin_data(
        selected_file, generate_cloud_text_files, generate_cloud_status_image
    )


class File_style:
    def __init__(self):
        print("Supported files: ")
        print("1: *.sal")
        print("2: inis.da")
        chosen_file = int(input("Enter a number: "))
        if chosen_file == 1:
            self.var_list = nube31_var_list  # list of variables to plot
            self.var_amount = len(
                nube31_var_list
            )  # amount of variables to plot per file
            self.binary_regex = r".sal$"  # regex to match binary files
            self.data_dimension = 3  # dimension of the data
            self.var_datatype = np.float32  # datatype of the variables
            self.var_structure_size = (
                biased_nx1 * biased_nx1 * nube31_biased_nz1
            )  # size of the variables
            self.var_structure = (
                biased_nx1,
                biased_nx1,
                nube31_biased_nz1,
            )  # structure of the variables
            self.file_name = "sal"
        elif chosen_file == 2:
            self.var_list = inis_var_list
            self.var_amount = len(inis_var_list)
            self.binary_regex = r"inis.da$"
            self.data_dimension = 1
            self.var_datatype = np.float32  # datatype of the variables
            self.var_structure_size = inis_biased_nz1  # size of the variables
            self.var_structure = inis_biased_nz1  # structure of the variables
            self.file_name = "inis"
        else:
            print("Invalid option")
            exit()


def check_path(path, selected_file_name):
    if not os.path.exists(path):
        os.makedirs(path)
    elif os.path.exists(path + selected_file_name + "/"):
        print(
            path
            + selected_file_name
            + " folder already exists, files will be overwritten"
        )
        user_input = input("Press Y to continue or any other key to exit: ").upper()
        if user_input != "Y":
            exit()
        else:
            os.system(f"rm -rf {path}{selected_file_name}/*")

    else:
        os.makedirs(path + selected_file_name + "/")


def get_file_list(data_path, binary_regex):
    files = os.listdir(data_path)  # List all files in the outputdata folder
    reg = re.compile(binary_regex)  # Compile the regex to match binary files
    return sorted(filter(reg.search, files))


def cloud_binary_comparison(selected_file: File_style):
    if os.path.exists(output_data_path) and os.path.exists(cmp_output_data_path):
        file_regex = selected_file.binary_regex
        file1_list = get_file_list(output_data_path, file_regex)
        file2_list = get_file_list(cmp_output_data_path, file_regex)
        diff = True
        for file1 in file1_list:
            if file1 in file2_list:
                comparison = Popen(
                    [
                        "cmp",
                        f"outputdata/{file1}",
                        f"outputdata1/{file1}",
                        "-b",
                    ],
                    stdout=PIPE,
                )
                cmp_result = comparison.stdout.read()  # needs to be saved to a variable
                if cmp_result != b"":
                    print(f"The binaries are different in: {file1}")
                    print(cmp_result)
                    diff = False
        if diff:
            print("All the binaries are the same")
    else:
        print("The outputdata folder does not exist")


def plot_style(variable, data_dimension):
    plt.grid(False)
    if data_dimension == 3:
        plt.imshow(variable[:, :, plot_center])
    elif data_dimension == 2:
        plt.imshow(variable[:, plot_center])
    elif data_dimension == 1:
        plt.plot(variable[3:-3])  # Cleans trash data
    plt.style.use("fivethirtyeight")
    # plt.colorbar() # Generates infinite colorbars in animation
    plt.xlabel("X")
    plt.ylabel("Y")


def get_var_from_data(data, var_iterator, selected_file: File_style):
    return data[var_iterator : var_iterator + selected_file.var_structure_size].reshape(
        selected_file.var_structure, order="F"
    )


def generate_image(i, var_number, selected_file: File_style, binary_files):
    file = binary_files[i]
    f = FortranFile(f"outputdata/{file}", "r")
    data = f.read_reals(selected_file.var_datatype)
    # directory = "img/" + str(file_counter) + file.split(".")[0] + "/"
    var_iterator = var_number * selected_file.var_structure_size
    variable = get_var_from_data(data, var_iterator, selected_file)
    plt.title(f"{nube31_var_list[var_number]} {str(2 * i)}s")
    plot_style(variable, selected_file.data_dimension)


def animate_variable(selected_file: File_style):
    var_number = int(input("Enter a variable number: "))
    fargs = [var_number, selected_file]
    file_regex = selected_file.binary_regex
    binary_files = get_file_list(output_data_path, file_regex)
    fargs.append(binary_files)
    # FuncAnimation will call generate_image with the arguments in fargs
    anim = FuncAnimation(
        plt.gcf(),
        generate_image,
        interval=1000,
        frames=len(binary_files),
        fargs=fargs,
        repeat=False,
    )
    plt.tight_layout()
    # If the user wants to save the animation as mp4 before showing it
    # Will get segfault because of plt.show() implementation closing the figure

    if input("Save animation as mp4? (Y/N):").upper() == "Y":
        writervideo = FFMpegWriter(fps=2)
        check_path("vid/", nube31_var_list[var_number])
        anim.save(f"vid/{nube31_var_list[var_number]}.mp4", writer=writervideo)
    if input("Show animation? (Y/N): ").upper() == "Y":
        plt.show()


def parse_status_img(img_path, data, selected_file: File_style, file_counter):
    full_img_path = f"{img_path}{selected_file.file_name}/{str(file_counter)}/"
    check_path(full_img_path, selected_file.file_name)
    var_iterator = 0
    for structure_iterator in range(0, len(data), selected_file.var_structure_size):
        variable = get_var_from_data(data, structure_iterator, selected_file)
        plt.title(f"{str(file_counter)} {selected_file.var_list[var_iterator]}")
        plot_style(variable, selected_file.data_dimension)
        plt.savefig(
            f"img/{selected_file.file_name}/{str(file_counter)}/{selected_file.var_list[var_iterator]}.png"
        )
        plt.close()  # If not closed, images will be superimposed
        if var_iterator < selected_file.var_amount:
            var_iterator += 1
        else:
            var_iterator = 0


def show_file_diff(file, selected_file_name):
    print(f"\n\n------------------- File: {file} -------------------\n")
    with open(f"txt/{selected_file_name}/{file}", "r") as file1:
        with open(f"txt1/{selected_file_name}/{file}", "r") as file2:
            for i, l1 in enumerate(file1, start=1):
                for l2 in file2:
                    if l1 != l2:
                        print("Line ", i, ": ")
                        print("File 1: ", l1)
                        print("File 2: ", l2)
                    break


def cloud_text_comparison(selected_file: File_style):
    original_path_extists = os.path.exists(f"txt/{selected_file.file_name}/")
    compate_path_extists = os.path.exists(f"txt1/{selected_file.file_name}/")
    if original_path_extists and compate_path_extists:
        diff = get_unequal_files(selected_file)
        if diff == []:
            print("All the text files are the same")
        else:
            print(f"The text files are different in: {diff}")
            if input("Show differences between files? Y/N: ").upper() == "Y":
                for file in diff:
                    show_file_diff(file, selected_file.file_name)
    else:
        print(
            f"error: -txt is created:{original_path_extists}\n-txt1 is created:{compate_path_extists}"
        )


def get_unequal_files(selected_file: File_style):
    txt_files = os.listdir(f"txt/{selected_file.file_name}/")
    txt1_files = os.listdir(f"txt1/{selected_file.file_name}/")
    if txt_files != txt1_files:
        print("Text file names differ, cannot compare")
        exit()
    status = cmpfiles(
        f"txt/{selected_file.file_name}/",
        f"txt1/{selected_file.file_name}/",
        txt_files,
        shallow=False,
    )
    return status[1]  # return unequal files


def parse_text_files(data, selected_file_name, file_counter, txt_path, file):
    full_txt_path = f"{txt_path}{selected_file_name}"
    header = f"File: {file}\n Variable: {selected_file_name}\n File number: {str(file_counter)}\n "
    file_name = f"{full_txt_path}{str(file_counter)}.txt"
    np.savetxt(file_name, data, newline=", \n", header=header)


def parse_bin_data(
    selected_file: File_style,
    generate_cloud_text_files="N",
    generate_cloud_status_image="N",
):
    if os.path.exists(output_data_path):
        selected_file_name = f"{selected_file.file_name}/"
        file_regex = selected_file.binary_regex
        binary_files = get_file_list(output_data_path, file_regex)

        if generate_cloud_status_image == "Y":
            img_path = "img/"
            check_path(img_path, selected_file_name)

        if generate_cloud_text_files == "Y":
            txt_path = "txt/"
            check_path(txt_path, selected_file_name)
            check_path(txt_path, selected_file_name)
        for file_counter, file in enumerate(binary_files, start=1):
            with FortranFile(f"outputdata/{file}", "r") as f:
                data = f.read_reals(selected_file.var_datatype)

            if generate_cloud_text_files == "Y":
                parse_text_files(data, selected_file_name, file_counter, txt_path, file)

            if generate_cloud_status_image == "Y":
                parse_status_img(img_path, data, selected_file, file_counter)
    else:
        print("The outputdata folder does not exist")


main()
