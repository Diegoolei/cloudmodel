# CloudModel

This repository hosts the re-engineered cloud simulation model, now implemented in Fortran 90 and wrapped for Python integration. Below are the instructions for setting up and running the model, as well as additional requirements for animation functionality and compilation with different Fortran compilers.

## Table of Contents

- [Python](#python)

  - [Installation](#instalation)
  - [Usage](#usage)
  - [Extra Functionality Requirements](#extra-functionality-requirements)

- [Fortran](#fortran-for-devs)

  - [Fortran Only](#fortran-only-compiling)
  - [F2Py](#f2py-integration)

- [Contributing](#contributing)
- [License](#license)

## Python

### Instalation

To install as a pip package:

```console
    pip install -i https://test.pypi.org/simple/ cloudmodel
```

### Usage

To run the model, first import it from the library

```python
    from cloudmodel import CloudModel
```

and call it with the desired parameters

```python
CloudModel(
    simulation_time_minutes,  #(int): The total simulation time in minutes.
    save_time_minutes,  # (int): The time interval in minutes at which the model state is saved.
    statistic_time_minutes,  # (int): The time interval in minutes at which statistics are calculated.
    bacup_time_minutes,  # (int): The time interval in minutes at which backups are created.
)
```

When the model finishes running there will be a new directory, where all the output data from the simulation will be stored.

The following section describes all the functions that can be applied to those simulation results.

### Features

```python
class FileStyle:
    """A class representing the file style for cloud data.

    Attributes:
        chosen_file (str): The chosen input file type.
        output_data_path (str): The path to the output data folder.
        cmp_output_data_path (str): The path to the comparison output data folder.
        img_path (str): The path to the image folder.
        txt_path (str): The path to the text folder.
        cmp_txt_path (str): The path to the comparison text folder.
        vid_path (str): The path to the video folder.
        img_option (str): The image style option.
        folder_handle (str): The folder handling option.

    Methods:
        _get_data(): Get the data from the selected files.
        _get_var_from_data(file_number, var_iterator): Get a specific variable from the data.
        list_var(): List all the variables.
        get_var(var, time): Get a specific variable at a given time.
        show_var_dataframe(var_array, center, axes): Show the variable data as a DataFrame.
        center_var(var_array, center, axes): Center the variable data along a specific axes.
        get_var_max_value_position(var_array): Get the position of the maximum value in the variable data.
        check_path(path, selected_file_name): Check if the path exists and create it if necessary.
        cloud_binary_comparison(): Compare the binary files in the output data and comparison output data folders.
        live_var_animation(variable): Create a live animation of a variable.
        plot_style(variable): Plot the style of a variable.
        generate_image(frame, var_number): Generate an image for a specific frame and variable.
        animate_variables(var_list, save_animation, show_animation): Animate multiple variables.
        animate_variable(var_to_animate, save_animation, show_animation, check_path): Animate a specific variable.
        parse_status_img(): Parse the status images.
        multi_var_img(var_1, var_2, file="inis.da"): Create an image with multiple variables.
        show_file_diff(file): Show the differences between two text files.
        cloud_text_comparison(): Compare text files in the specified paths and display the differences, if any.
        get_unequal_files(): Get the unequal files between two folders.
        parse_text_files(): Parse the text files.
    """
```

A run + image generation for all variables to all the simulated times would look like this:

```python
    from cloudmodel import CloudModel
    
    CloudModel(
        simulation_time_minutes=45,
        save_time_minutes=3,
        statistic_time_minutes=3,
        bacup_time_minutes=3,
    )

    
    cloud = FileStyle(
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
    cloud.parse_status_img()
    theta_base = cloud.get_var(cloud.var_list[3])
    cloud.show_var_dataframe(theta_base, 3)
```

## Extra Functionality Requirements

To display animations, follow these steps:

1. Install `ffmpeg`:

    ```console
    sudo apt install ffmpeg
    ```

2. If you're using Windows Subsystem for Linux (WSL), you can show Matplotlib plots and other GUI elements by following these [steps](https://stackoverflow.com/questions/43397162/show-matplotlib-plots-and-other-gui-in-ubuntu-wsl1-wsl2):
    - Install [Xming X Server for Windows](https://sourceforge.net/projects/xming/).
    - Run the following command in your WSL2 terminal:

      ```console
      sudo apt-get install python-tk
      ```

## Fortran for devs

This section is focused on the steps needed to re-compile the model to mantain usability with the Python package.

### F2Py Integration

To compile for F2Py, follow these steps:

1. Navigate to the `interface` directory:

    ```console
    cd interface
    ```

2. Compile using `make`:

    ```console
    make
    ```

---

### Fortran-only compiling

#### gfortran

Linux:

Install `gfortran`:

```console
sudo apt install gfortran
```

MacOS:

Install `gfortran` using Homebrew:

```console
brew install gfortran
xcode-select --install
```

---
Execution:

Run the following command with `fpm`:

```console
fpm run --profile release
```

#### nvfortran

##### Installation

1. Install the appropriate Nvidia drivers for your system.
2. Install the [Nvidia CUDA toolkit](https://developer.nvidia.com/cuda-toolkit).
3. Install the [Nvidia HPC SDK](https://developer.nvidia.com/nvidia-hpc-sdk-downloads). The installation path is usually `/opt/nvidia/hpc_sdk/Linux_x86_64/(version)/compilers/bin`. Add it to your PATH.

Execution:

Run the following command with `fpm`:

```console
fpm run --compiler "/opt/nvidia/hpc_sdk/Linux_x86_64/(version)/compilers/bin/nvfortran" --flag "-O3 -cuda"
```

## Contributing

1. Fork the repository.
2. Create a new branch: `git checkout -b feature-name`.
3. Make your changes.
4. Push your branch: `git push origin feature-name`.
5. Create a pull request.

## License

This project is licensed under the [MIT License](LICENSE).
