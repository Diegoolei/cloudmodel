"""Module to read and compare cloud model data."""

import os
import re
from enum import Enum

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

from scipy.io import FortranFile

from .constants import (
    biased_nx1,
    inis_biased_nz1,
    inis_var_list,
    nube31_biased_nz1,
    nube31_var_list,
    plot_center,
)
from .interface import c_interface as nb


class ImageStyle(Enum):
    """Enum to select the style of the image."""

    IMAGE = "Image"
    CONTOUR = "Contour"


class InputFileType(Enum):
    """Enum to select the input file."""

    NUBE31 = "Nube"
    INIS = "Inis"


class FolderHandle(Enum):
    """Enum to select how to respond if destiny folder exists."""

    IGNORE = "Ignore"
    DELETE = "Delete"
    CANCEL = "Cancel"


class CloudSimulation:
    """
    Represents a cloud model.

    Args
    ----
        simulation_time_minutes (int): The total simulation time in minutes.
        save_time_minutes (int): The time interval in minutes at which the
            model state is saved.
        statistic_time_minutes (int): The time interval in minutes at which
            statistics are calculated.
        bacup_time_minutes (int): The time interval in minutes at which
            backups are created.
    """

    def __init__(
        self,
        simulation_time_minutes=0,
        save_time_minutes=0,
        statistic_time_minutes=0,
        bacup_time_minutes=0,
        restore_backup=False,
        directory=".temp/",
    ):
        self.simulation_time_minutes = simulation_time_minutes
        self.save_time_minutes = save_time_minutes
        self.statistic_time_minutes = statistic_time_minutes
        self.bacup_time_minutes = bacup_time_minutes
        self.restore_backup = restore_backup
        self.directory = directory
        self.initial_analytics: FileStyle = None
        self.cloud_analytics: FileStyle = None

    def run_model(self):
        """Run the cloud model."""
        check_path(FolderHandle.IGNORE.value, self.directory)
        nb.run_model_python(
            self.simulation_time_minutes,
            self.save_time_minutes,
            self.statistic_time_minutes,
            self.bacup_time_minutes,
            self.restore_backup,
            self.directory,
        )
        self.run_initial_analysis()
        self.run_cloud_analysis()

    def run_initial_analysis(self):
        """
        Run the initial analysis on the data.

        Returns
        -------
            FileStyle: An instance of the FileStyle class representing
                the analysis results.
        """
        data = FileStyle(
            chosen_file="Inis",
            output_data_path=self.directory,
            cmp_output_data_path="outputdata1/",
            img_path="img/" + self.directory,
            txt_path="txt/" + self.directory,
            cmp_txt_path="txt1/" + self.directory,
            vid_path="vid/" + self.directory,
            img_option="Contour",
            folder_handle="Delete",
        )
        self.initial_analytics = data

    def run_cloud_analysis(self):
        """
        Run the cloud analysis.

        Returns
        -------
            FileStyle: The cloud analysis object.
        """
        data = FileStyle(
            chosen_file="Nube",
            output_data_path=self.directory,
            cmp_output_data_path="outputdata1/",
            img_path="img/" + self.directory,
            txt_path="txt/" + self.directory,
            cmp_txt_path="txt1/" + self.directory,
            vid_path="vid/" + self.directory,
            img_option="Contour",
            folder_handle="Delete",
        )
        self.cloud_analytics = data

    def clean_model(self):
        """
        Clean the model by deleting the files in the specified directory.

        If the directory is set to ".temp/", this method will delete all
        the files in that directory using the `check_path` function with
        the `FolderHandle.DELETE` option.

        Note: This method assumes that the `check_path` function is defined
            elsewhere.

        Returns
        -------
            None
        """
        if self.directory == ".temp/":
            check_path(FolderHandle.DELETE.value, self.directory)


class FileStyle:
    """
    A class representing the file style for cloud data.

    Attributes
    ----------
        chosen_file (str): The chosen input file type.
        output_data_path (str): The path to the output data folder.
        cmp_output_data_path (str): The path to the comparison output
            data folder.
        img_path (str): The path to the image folder.
        txt_path (str): The path to the text folder.
        cmp_txt_path (str): The path to the comparison text folder.
        vid_path (str): The path to the video folder.
        img_option (str): The image style option.
        folder_handle (str): The folder handling option.

    Methods
    -------
        _get_data(): Get the data from the selected files.
        _get_var_from_data(file_number, var_iterator): Get a variable
            from the data.
        list_var(): List all the variables.
        get_var(var, time): Get a specific variable at a given time.
        show_var_dataframe(var_array, center, axis): Show a variable as
            a DataFrame.
        center_var(var_array, center, axis): Center a variable along a
            given axis.
        get_var_max_value_position(var_array): Get the position of the
            maximum value in a variable.
        cloud_binary_comparison(): Compare binary files in the output
            and comparison folders.
        live_var_animation(variable): Create a live variable animation.
        plot_style(variable): Plot the style of a variable.
        generate_image(frame, var_number): Generates an image for a given
            frame and variable number.
        animate_variables(var_list=None, save_animation=True,
            show_animation=False): Animate multiple variables.
        animate_variable(var_to_animate, save_animation=True,
            show_animation=False, check_path=True): Animates a variable
            from the output data.
        parse_status_img(): Animates the variables in the given var_list.
        multi_var_img(var_1, var_2): Generate a scatter plot of two variables
            from the data file.
        show_file_diff(file_name): Show the differences between two text files
            with the same name in different directories.
        cloud_text_comparison(): Compare text files in the specified paths and
            display the differences, if any.
        get_unequal_files(): Get files with identical names but unequal content
            between two directories.
        parse_text_files(): Parse the text files.
    """

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
        assert os.path.exists(
            output_data_path
        ), "output_data_path does not exist"
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
        self._get_data()

    def _get_data(self):
        """Get the data from the selected files."""
        selected_files = _get_file_list(
            self.output_data_path, self.binary_regex
        )
        for file in selected_files:
            with FortranFile(f"{self.output_data_path}{file}", "r") as f:
                self.data_file.append(f.read_reals(self.var_datatype))

    def _get_var_from_data(self, file_number, var_iterator):
        """
        Get a variable from the data.

        Args
        ----
            file_number (int): The index of the file in the data_file list.
            var_iterator (int): The iterator for the variable in the file.

        Returns
        -------
            numpy.ndarray: The variable data as a numpy array.
        """
        actual_file_data = self.data_file[file_number]
        return actual_file_data[
            var_iterator : var_iterator + self.var_structure_size
        ].reshape(self.var_structure, order="F")

    def list_var(self):
        """List all the variables."""
        for var in self.var_list:
            print(var)

    def get_var(self, var: str, time: int):
        """
        Get a specific variable at a given time.

        Args
        ----
            var (str): The name of the variable.
            time (int): The time index.

        Returns
        -------
            numpy.ndarray: The variable data as a numpy array.
        """
        var_index = self.var_list.index(var)
        var_iterator = var_index * self.var_structure_size
        return self._get_var_from_data(time, var_iterator)

    def show_var_dataframe(self, var_array, center: tuple, axis: str):
        """
        Show a variable as a DataFrame.

        Args
        ----
            var_array (numpy.ndarray): The variable data as a numpy array.
            center (tuple): The center coordinates for the variable.
            axis (str): The axis along which to center the variable.
        """
        df = pd.DataFrame(self.center_var(var_array, center, axis))
        print(df)

    def center_var(self, var_array, center: tuple, axis: str):
        """
        Center a variable along a given axis.

        Args
        ----
            var_array (numpy.ndarray): The variable data as a numpy array.
            center (tuple): The center coordinates for the variable.
            axis (str): The axis along which to center the variable.

        Returns
        -------
            numpy.ndarray: The centered variable data as a numpy array.
        """
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
        """
        Get the position of the maximum value in a variable.

        Args
        ----
            var_array (numpy.ndarray): The variable data as a numpy array.

        Returns
        -------
            tuple: The position of the maximum value in the variable.
        """
        return np.where(var_array == np.max(var_array))

    def plot_style(self, variable):
        """
        Apply a specific plot style based on the image option and data dimen.

        Parameters
        ----------
            variable (numpy.ndarray): The variable to be plotted.

        Raises
        ------
            ValueError: If the image style is invalid.

        Returns
        -------
            None
        """
        match self.img_option:
            case ImageStyle.IMAGE.value:
                if self.data_dimension == 3:
                    plt.imshow(
                        self.center_var(
                            variable,
                            self.get_var_max_value_position(variable),
                            "z",
                        )
                        # np.flipud(variable[:, :, plot_center])
                    )
                elif self.data_dimension == 2:
                    plt.imshow(variable[:, plot_center])
                elif self.data_dimension == 1:
                    plt.plot(variable[3:-3])  # Cleans trash data
            case ImageStyle.CONTOUR.value:
                fig, ax = plt.subplots()
                if self.data_dimension == 3:
                    cs = ax.contour(
                        # np.flipud(variable[:, :, plot_center]),
                        self.center_var(
                            variable,
                            self.get_var_max_value_position(variable),
                            "z",
                        ),
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

    def parse_status_img(self):
        """
        Parse the status image data and generates plots for each variable.

        This method iterates over the data files and generates plots for
        each variable in the status image. The plots are saved as PNG
        images in the specified image path.

        Returns
        -------
            None
        """
        check_path(self.folder_handle, f"{self.img_path}{self.file_name}")
        for file_iterator in range(len(self.data_file)):
            full_img_path = (
                f"{self.img_path}{self.file_name}/{str(file_iterator)}/"
            )
            check_path(self.folder_handle, full_img_path)
            var_iterator = 0
            for structure_iterator in range(
                0, len(self.data_file[file_iterator]), self.var_structure_size
            ):
                variable = self._get_var_from_data(
                    file_iterator, structure_iterator
                )
                filenumber = str(file_iterator)
                plt.title(f"{filenumber} {self.var_list[var_iterator]}")
                self.plot_style(variable)
                dirname = f"{self.img_path}{self.file_name}"
                filename = f"{self.var_list[var_iterator]}.png"
                plt.savefig(f"{dirname}/{filenumber}/{filename}")
                plt.close()  # If not closed, images will be superimposed
                if var_iterator < self.var_amount:
                    var_iterator += 1
                else:
                    var_iterator = 0

    def multi_var_img(self, var_1, var_2):
        """
        Generate a scatter plot of two variables from the data file.

        Args
        ----
            var_1 (int): Index of the first variable to plot.
            var_2 (int): Index of the second variable to plot.
            file (str, optional): Name of the data file. Defaults to "inis.da".

        Returns
        -------
            None
        """
        check_path(
            self.folder_handle, f"{self.img_path}{self.file_name}/multivar/"
        )
        print(len(self.data_file))
        var_1_data = self._get_var_from_data(0, var_1)
        var_2_data = self._get_var_from_data(0, var_2)
        plt.title(f"{self.var_list[var_1]} vs {self.var_list[var_2]}")
        # self.plot_style(variable, self.data_dimension)
        plt.plot(var_1_data[10:-10], var_2_data[10:-10], ".b")
        # plt.plot(var_1_data[::6], var_2_data[::6], "*")
        dirname = f"{self.img_path}{self.file_name}"
        filename = f"{self.var_list[var_1]}_{self.var_list[var_2]}.png"
        plt.savefig(f"{dirname}/multivar/{filename}")
        plt.close()  # If not closed, images will be superimposed


def _get_file_list(data_path, binary_regex):
    """
    Get a sorted list of files in the specified data path that match the given.

        binary regex pattern.

    Args
    ----
        data_path (str): The path to the data directory.
        binary_regex (str): The regular expression pattern to match binary
            files.

    Returns
    -------
        list: A sorted list of file names that match the binary regex pattern.
    """
    files = os.listdir(data_path)  # List all files in outputdata folder
    reg = re.compile(binary_regex)  # Compile regex to match binary files
    return sorted(filter(reg.search, files))


def check_path(folder_handle, path, selected_file_name=""):
    """
    Check if a path exists and create it if necessary.

    Args
    ----
        path (str): The path to check.
        selected_file_name (str, optional): The name of the selected file.
            Defaults to "".
    """
    if not os.path.exists(path):
        os.makedirs(path)
    elif not os.path.exists(path + selected_file_name + "/"):
        os.makedirs(path + selected_file_name + "/")
    else:
        match folder_handle:
            case FolderHandle.IGNORE.value:
                return
            case FolderHandle.DELETE.value:
                os.system(f"rm -rf {path}{selected_file_name}/*")
            case FolderHandle.CANCEL.value:
                raise ValueError("Folder already exists")
            case _:
                raise ValueError("Invalid Folder Handle")
