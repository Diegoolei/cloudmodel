import os
import re
from datetime import datetime
from enum import Enum
from filecmp import cmpfiles
from subprocess import PIPE, Popen

from constants import (
    biased_nx1,
    inis_biased_nz1,
    inis_var_list,
    nube31_biased_nz1,
    nube31_var_list,
    plot_center,
)

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation

import numpy as np

import pandas as pd

from scipy.io import FortranFile


def time_it(func):
    """Log the date and time of a function"""

    def wrapper(*args, **kwargs):
        date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"Function: {func.__name__}\nRun on: {date}")
        print(f'{"-"*30}')
        func(*args, **kwargs)

    return wrapper


class ImageStyle(Enum):
    """Enum to select the style of the image"""

    IMAGE = "Image"
    CONTOUR = "Contour"


class InputFileType(Enum):
    """Enum to select the input file"""

    NUBE31 = "Nube"
    INIS = "Inis"


class FolderHandle(Enum):
    """Enum to select how to respond if destiny folder exists"""

    IGNORE = "Ignore"
    DELETE = "Delete"
    CANCEL = "Cancel"


class FileStyle:
    def __init__(
        self,
        chosen_file=InputFileType.NUBE31.value,
        output_data_path="outputdata/",
        cmp_output_data_path="outputdata1/",
        img_path="img/",
        txt_path="txt/",
        cmp_txt_path="txt1/",
        vid_path="vid/",
        img_option=ImageStyle.IMAGE.value,
        folder_handle=FolderHandle.IGNORE.value,
    ):
        assert os.path.exists(output_data_path), "output_data_path does not exist"
        self.output_data_path = output_data_path
        self.cmp_output_data_path = cmp_output_data_path
        self.img_path = img_path
        self.txt_path = txt_path
        self.cmp_txt_path = cmp_txt_path
        self.vid_path = vid_path
        self.data_file = []
        self.folder_handle = folder_handle
        match chosen_file:
            case InputFileType.NUBE31.value:
                self.img_option = img_option
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
            case InputFileType.INIS.value:
                self.img_option = ImageStyle.IMAGE.value
                self.var_list = inis_var_list
                self.var_amount = len(inis_var_list)
                self.binary_regex = r"inis.da$"
                self.data_dimension = 1
                self.var_datatype = np.float32  # datatype of the variables
                self.var_structure_size = inis_biased_nz1  # var size
                self.var_structure = inis_biased_nz1  # var structure
                self.file_name = "inis"
            case _:
                raise ValueError("Invalid Input File Type")
        self.get_data()

    def get_data(self):
        selected_files = get_file_list(self.output_data_path, self.binary_regex)
        for file in selected_files:
            with FortranFile(f"{self.output_data_path}{file}", "r") as f:
                self.data_file.append(f.read_reals(self.var_datatype))

    def get_var_from_data(self, file_number, var_iterator):
        actual_file_data = self.data_file[file_number]
        return actual_file_data[
            var_iterator: var_iterator + self.var_structure_size
        ].reshape(self.var_structure, order="F")

    def list_var(self):
        for var in self.var_list:
            print(var)

    def get_var(self, var: str, time: int):
        var_index = self.var_list.index(var)
        var_iterator = var_index * self.var_structure_size
        return self.get_var_from_data(time, var_iterator)

    def show_var_dataframe(self, var_array, center: tuple, axis: str):
        df = pd.DataFrame(self.center_var(var_array, center, axis))
        print(df)

    def center_var(self, var_array, center: tuple, axis: str):
        match axis:
            case "x":
                return var_array[center[0][0], :, :]
            case "y":
                return var_array[:, center[1][0], :]
            case "z":
                return var_array[:, :, center[2][0]]
            case _:
                raise ValueError("Invalid Axis")

    def get_var_max_value_position(self, var_array):
        return np.where(var_array==np.max(var_array))

    def check_path(self, path, selected_file_name=""):
        """Checks if the path exists, if not, creates it
        If the path exists, checks if the selected_file_name folder exists,
        if not, creates it
        If no selected_file_name is given, it will only check if the path
        exists"""

        if not os.path.exists(path):
            os.makedirs(path)
        elif not os.path.exists(path + selected_file_name + "/"):
            os.makedirs(path + selected_file_name + "/")
        else:
            match self.folder_handle:
                case FolderHandle.IGNORE.value:
                    return
                case FolderHandle.DELETE.value:
                    os.system(f"rm -rf {path}{selected_file_name}/*")
                case FolderHandle.CANCEL.value:
                    raise ValueError("Folder already exists")
                case _:
                    raise ValueError("Invalid Folder Handle")

    def cloud_binary_comparison(self):
        if os.path.exists(self.cmp_output_data_path):
            file_regex = self.binary_regex
            file1_list = get_file_list(self.output_data_path, file_regex)
            file2_list = get_file_list(self.cmp_output_data_path, file_regex)
            diff = True
            for file1 in file1_list:
                if file1 in file2_list:
                    comparison = Popen(
                        [
                            "cmp",
                            f"{self.output_data_path}{file1}",
                            f"{self.cmp_output_data_path}{file1}",
                            "-b",
                        ],
                        stdout=PIPE,
                    )
                    cmp_result = (
                        comparison.stdout.read()
                    )  # needs to be saved to a variable
                    if cmp_result != b"":
                        print(f"The binaries are different in: {file1}")
                        print(cmp_result)
                        diff = False
            if diff:
                print("All the binaries are the same")
        else:
            print(f"There is no {self.cmp_output_data_path} folder to compare with")

    def live_var_animation(self, variable):
        ax = plt.axes(projection='3d')
        num_layers = 50
        num_pnt = 57
        z, x, y = np.meshgrid(np.arange(1, num_layers + 1), np.arange(num_pnt), np.arange(num_pnt), indexing='ij')

        ax.scatter3D(x, y, z, c=variable, cmap='viridis')
        plt.show()

    def plot_style(self, variable):
        match self.img_option:
            case ImageStyle.IMAGE.value:
                if self.data_dimension == 3:
                    plt.imshow(
                        self.center_var(variable, self.get_var_max_value_position(variable), "z")
                        #np.flipud(variable[:, :, plot_center])
                        )
                elif self.data_dimension == 2:
                    plt.imshow(variable[:, plot_center])
                elif self.data_dimension == 1:
                    plt.plot(variable[3:-3])  # Cleans trash data
            case ImageStyle.CONTOUR.value:
                fig, ax = plt.subplots()
                if self.data_dimension == 3:
                    cs = ax.contour(
                        #np.flipud(variable[:, :, plot_center]),
                        self.center_var(variable, self.get_var_max_value_position(variable), "z"),
                        linewidths=0.3,
                        colors="k",
                    )
                    ax.clabel(cs, fontsize=6, inline=True)
                elif self.data_dimension == 2:
                    cs = ax.contour(
                        variable[:, plot_center], linewidths=0.3, colors="k"
                    )
                    ax.clabel(cs, fontsize=6, inline=True)
            case _:
                raise ValueError("Invalid Image Style")
        plt.grid(False)
        plt.style.use("fivethirtyeight")
        # plt.colorbar() # Generates infinite colorbars in animation
        plt.xlabel("X")
        plt.ylabel("Y")

    def generate_image(self, frame, var_number):
        # directory = "img/" + str(file_counter) + file.split(".")[0] + "/"
        variable = self.get_var_from_data(frame, var_number)
        plt.title(f"{nube31_var_list[var_number]} {str(2 * frame)}s")
        self.plot_style(variable)

    def animate_variables(
        self, var_list=None, save_animation=True, show_animation=False
    ):
        if var_list is None:
            var_list = range(self.var_amount)
        self.check_path(self.vid_path)
        for var in var_list:
            self.animate_variable(var, save_animation, show_animation, check_path=False)

    def animate_variable(
        self, var_to_animate, save_animation=True, show_animation=False, check_path=True
    ):
        fargs = [var_to_animate]
        file_ammount = len(get_file_list(self.output_data_path, self.binary_regex))
        # FuncAnimation will call generate_image with the arguments in fargs
        anim = FuncAnimation(
            plt.gcf(),
            self.generate_image,
            interval=1000,
            frames=file_ammount,
            fargs=fargs,
            repeat=False,
        )
        plt.tight_layout()
        # If the user wants to save the animation as mp4 before showing it
        # Will get segfault because of plt.show() implementation closing the figure
        if save_animation:
            writervideo = FFMpegWriter(fps=2)
            if check_path:
                self.check_path(self.vid_path)
            anim.save(
                f"{self.vid_path}{nube31_var_list[var_to_animate]}.mp4",
                writer=writervideo,
            )
            plt.close()
        if show_animation:
            plt.show()

    def parse_status_img(self):
        self.check_path(f"{self.img_path}{self.file_name}")
        for file_iterator in range(len(self.data_file)):
            full_img_path = f"{self.img_path}{self.file_name}/{str(file_iterator)}/"
            self.check_path(full_img_path)
            var_iterator = 0
            for structure_iterator in range(
                0, len(self.data_file[file_iterator]), self.var_structure_size
            ):
                variable = self.get_var_from_data(file_iterator, structure_iterator)
                plt.title(f"{str(file_iterator)} {self.var_list[var_iterator]}")
                self.plot_style(variable)
                plt.savefig(
                    f"{self.img_path}{self.file_name}/{str(file_iterator)}/{self.var_list[var_iterator]}.png"
                )
                plt.close()  # If not closed, images will be superimposed
                if var_iterator < self.var_amount:
                    var_iterator += 1
                else:
                    var_iterator = 0

    def multi_var_img(self, var_1, var_2, file="inis.da"):
        self.check_path(f"{self.img_path}{self.file_name}/multivar/")
        print(len(self.data_file))
        var_1_data = self.get_var_from_data(0, var_1)
        var_2_data = self.get_var_from_data(0, var_2)
        plt.title(f"{self.var_list[var_1]} vs {self.var_list[var_2]}")
        # self.plot_style(variable, self.data_dimension)
        plt.plot(var_1_data[10:-10], var_2_data[10:-10], ".b")
        # plt.plot(var_1_data[::6], var_2_data[::6], "*")
        plt.savefig(
            f"{self.img_path}{self.file_name}/multivar/{self.var_list[var_1]}_{self.var_list[var_2]}.png"
        )
        plt.close()  # If not closed, images will be superimposed

    def show_file_diff(self, file):
        print(f"\n\n------------------- File: {file} -------------------\n")
        with open(f"{self.txt_path}{self.file_name}/{file}", "r") as file1:
            with open(f"{self.cmp_txt_path}{self.file_name}/{file}", "r") as file2:
                for i, l1 in enumerate(file1, start=1):
                    for l2 in file2:
                        if l1 != l2:
                            print("Line ", i, ": ")
                            print("File 1: ", l1)
                            print("File 2: ", l2)
                        break

    def cloud_text_comparison(self):
        original_path_extists = os.path.exists(f"{self.txt_path}{self.file_name}/")
        compate_path_extists = os.path.exists(f"{self.cmp_txt_path}{self.file_name}/")
        if original_path_extists and compate_path_extists:
            diff = self.get_unequal_files()
            if diff == []:
                print("All the text files are the same")
            else:
                print(f"The text files are different in: {diff}")
                if input("Show differences between files? Y/N: ").upper() == "Y":
                    for file in diff:
                        self.show_file_diff(file)
        else:
            print(
                f"error: -txt is created:{original_path_extists}\n-txt1 is created:{compate_path_extists}"
            )

    def get_unequal_files(self):
        txt_files = os.listdir(f"{self.txt_path}{self.file_name}/")
        txt1_files = os.listdir(f"{self.cmp_txt_path}{self.file_name}/")
        if txt_files != txt1_files:
            raise ("Text file names differ, cannot compare")
        status = cmpfiles(
            f"{self.txt_path}{self.file_name}/",
            f"{self.cmp_txt_path}{self.file_name}/",
            txt_files,
            shallow=False,
        )
        return status[1]  # return unequal files

    def parse_text_files(self):
        full_txt_path = f"{self.txt_path}{self.file_name}/"
        self.check_path(full_txt_path)
        for file_counter in range(len(self.data_file)):
            header = f"File: {file_counter}\n Variable: {self.file_name}/\n File number: {str(file_counter)}\n "
            file_name = f"{full_txt_path}{str(file_counter)}.txt"
            np.savetxt(
                file_name, self.data_file[file_counter], newline=", \n", header=header
            )


def get_file_list(data_path, binary_regex):
    files = os.listdir(data_path)  # List all files in outputdata folder
    reg = re.compile(binary_regex)  # Compile regex to match binary files
    return sorted(filter(reg.search, files))


def data_comparison(original_data: FileStyle, cmp_data: FileStyle):
    for iterator in range(len(original_data.data_file)):
        if not np.allclose(
            original_data.data_file[iterator],
            cmp_data.data_file[iterator],
            rtol=1e-05,
            atol=3e-06,
        ):
            # return False
            print(f"Data is different in data[{iterator}]: ")
            for it_var in range(len(original_data.var_list)):
                var = original_data.get_var_from_data(iterator, it_var)
                cmp_var = cmp_data.get_var_from_data(iterator, it_var)
                # Normalize var values
                if var.max() != 0:
                    norm_var = var / np.max(np.abs(var))
                    norm_cmp_var = cmp_var / np.max(np.abs(var))
                else:
                    norm_var = var
                    norm_cmp_var = cmp_var
                if not np.allclose(norm_var, norm_cmp_var, rtol=1e-05, atol=3e-06):
                    print(f"Variable {original_data.var_list[it_var]} is different ")
                    print(f"Max difference: {np.argmax(np.abs(var - cmp_var))}\n")

    return True
