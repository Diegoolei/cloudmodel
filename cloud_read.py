from scipy.io import FortranFile
import matplotlib.pyplot as plt
from constants import *
import numpy as np
import os
import re  # Add the re import declaration to use regex

output_data_path = "outputdata/"
var_amount = 12
binary_regex = r"^nube31"
plot_center = 10
var_datatype = np.float32
var_structure_size = biased_nx1 * biased_nx1 * biased_nz1
var_structure = (biased_nx1, biased_nx1, biased_nz1)

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

    files = os.listdir(output_data_path)  # List all files in the outputdata folder
    reg = re.compile(binary_regex)  # Compile the regex to match binary files
    binary_files = list(
        filter(reg.search, files)
    )  # Create iterator using filter, cast to list
    file_counter = 0
    for file in binary_files:
        f = FortranFile("outputdata/" + file, "r")
        data = f.read_reals(var_datatype)
        directory = "img/" + str(file_counter) + file.split(".")[0] + "/"
        file_counter += 1
        if not os.path.exists("img/" + str(file_counter) + "/"):
            # TODO - Fix F77 filename strings, then replace file_counter with directory
            os.makedirs("img/" + str(file_counter) + "/")
        var_iterator = 0
        for i in range(0, len(data), var_structure_size):
            variable = data[i: i + var_structure_size].reshape(
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
