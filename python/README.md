# README

## Table of Contents

- [Python](#python)

-- [Installation](#instalation)
-- [Usage](#usage)

## Python

### Instalation

To install as a pip package:

```console
    pip install cloudmodel
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
    backup_time_minutes,  # (int): The time interval in minutes at which backups are created.
    directory,  # (str): Output Data Destination.
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
        _get_var_max_value_position(var_array): Get the position of the maximum value in the variable data.
        _plot_style(variable): Plot the style of a variable.
        get_var(var, time): Get a specific variable at a given time.
        get_var_dataframe(var_array, center, axes): return the variable data as a DataFrame.
        center_var(var_array, center, axes): Center the variable data along a specific axes.
        parse_status_img(): Parse the status images.
        multi_var_img(var_1, var_2, file="inis.da"): Create an image with multiple variables.
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
