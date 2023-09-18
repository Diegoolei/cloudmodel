from scipy.io import FortranFile
import matplotlib.pyplot as plt
from constants import *
import numpy as np
import os
import re  # Add the re import declaration to use regex
from subprocess import Popen, PIPE

output_data_path = "outputdata/"  # original data path
cmp_output_data_path = "outputdata1/"  # comparison data path
var_amount = 12  # amount of variables to plot per file
binary_regex = r".sal$"  # regex to match binary files
plot_center = 10  # center of the plot
var_datatype = np.float32  # datatype of the variables
var_structure_size = biased_nx1 * biased_nx1 * biased_nz1  # size of the variables
var_structure = (biased_nx1, biased_nx1, biased_nz1)  # structure of the variables

compare_binaries = input("Compare binaries? (Y/N): ").upper()
generate_cloud_status_image = input("Generate cloud status image? (Y/N): ").upper()


def get_file_list(data_path):
    files = os.listdir(data_path)  # List all files in the outputdata folder
    reg = re.compile(binary_regex)  # Compile the regex to match binary files
    binary_files = list(
        filter(reg.search, files)
    )  # Create iterator using filter, cast to list
    return binary_files


def cloud_binary_comparison():
    if os.path.exists(output_data_path):
        file1_list = get_file_list(output_data_path)
        file2_list = get_file_list(cmp_output_data_path)
        diff = True
        for file1 in file1_list:
            if file1 in file2_list:
                comparison = Popen(
                    ["cmp", "outputdata/" + file1, "outputdata1/" + file1], stdout=PIPE
                )
                if comparison.stdout.read() != b"":
                    print("The binaries are different in: " + file1)
                    print(
                        "d" + str(comparison.stdout.read())
                    )  # FIXME - Wont show the difference
                    diff = False
        if diff:
            print("All the binaries are the same")
    else:
        print("The outputdata folder does not exist")


if compare_binaries == "Y":
    cloud_binary_comparison()


def generate_cloud_status_img():
    if os.path.exists(output_data_path):
        if not os.path.exists("img/"):
            os.makedirs("img/")
        else:
            print("The img folder already exists, files will be overwritten")
            user_input = input("Press Y to continue or any other key to exit: ").upper()
            if user_input != "Y":
                exit()
            else:
                os.system("rm -rf img/*")

        binary_files = get_file_list(output_data_path)
        file_counter = 0
        for file in binary_files:
            f = FortranFile("outputdata/" + file, "r")
            data = f.read_reals(var_datatype)
            # directory = "img/" + str(file_counter) + file.split(".")[0] + "/"
            file_counter += 1
            if not os.path.exists("img/" + str(file_counter) + "/"):
                # TODO - Fix F77 filename strings, then replace file_counter with directory
                os.makedirs("img/" + str(file_counter) + "/")
            var_iterator = 0
            for i in range(0, len(data), var_structure_size):
                variable = data[i : i + var_structure_size].reshape(
                    (var_structure[0], var_structure[1], var_structure[2]), order="F"
                )
                plt.imshow(variable[:, :, plot_center])
                plt.savefig(
                    "img/"
                    + str(file_counter)
                    + "/"
                    + nube31_var_list[var_iterator]
                    + ".png"
                )
                plt.close()
                if var_iterator < var_amount:
                    var_iterator += 1
                else:
                    var_iterator = 0
    else:
        print("The outputdata folder does not exist")


if generate_cloud_status_image == "Y":
    generate_cloud_status_img()
