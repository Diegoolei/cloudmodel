"""main file to run the cloud model."""

import os
import re

import numpy as np

import pytest

from ..cloud_read import (
    CloudSimulation,
    FileStyle,
    check_path,
    inis_var_list,
    nube31_var_list,
)


def verify_vars(model: CloudSimulation, equal: bool):
    """
    Verify the variables in the CloudSimulation model.

    Args
    ----
        model (CloudSimulation): The CloudSimulation model to verify.
        equal (bool): Flag indicating whether to check for equality or
            inequality.

    Raises
    ------
        AssertionError: If the variables do not meet the specified condition.

    Returns
    -------
        None
    """
    temp = FileStyle(
        chosen_file="Nube",
        output_data_path="cloudmodel/tests/test_data",
        img_path=f"img/{model.directory}",
        img_option="Contour",
        folder_handle="Delete",
    )
    f = True
    for var in nube31_var_list:
        for time in range(
            len(
                _get_file_list(
                    model.directory, model.cloud_analytics.binary_regex
                )
            )
        ):
            variable = model.cloud_analytics.get_var(var, time)
            temp_var = temp.get_var(var, time)

            if equal:
                assert np.allclose(
                    variable, temp_var
                ), f"Max value position of {var} is different"
            else:
                f = (
                    f and not np.allclose(variable, temp_var),
                    f"Max value position of {var} is equal",
                )
    assert f


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


def test_simulation_equal():
    """Execute function to run the cloud model."""
    model = CloudSimulation(
        simulation_time_minutes=1,
        save_time_minutes=1,
        statistic_time_minutes=1,
        bacup_time_minutes=1,
        restore_backup=False,
    )

    model.clean_model()

    model.run_model()

    verify_vars(model, True)

    model.clean_model()


def test_simulation_different():
    """Execute function to run the cloud model."""
    model = CloudSimulation(
        simulation_time_minutes=1,
        save_time_minutes=1,
        statistic_time_minutes=1,
        bacup_time_minutes=1,
        restore_backup=False,
    )
    model.clean_model()

    model.run_model()

    verify_vars(model, False)

    model.clean_model()


def test_get_data():
    """Test the get_data method."""
    model = CloudSimulation(
        simulation_time_minutes=1,
        save_time_minutes=1,
        statistic_time_minutes=1,
        bacup_time_minutes=1,
        restore_backup=False,
    )

    model.clean_model()

    model.run_model()

    model.cloud_analytics._get_data()

    model.clean_model()

    assert isinstance(model.cloud_analytics.data_file, list)
    assert len(model.cloud_analytics.data_file) > 0


def test_get_data_with_invalid_file():
    with pytest.raises(ValueError):
        FileStyle(
            chosen_file="asd",
            output_data_path="./temp",
            img_path="img",
            img_option="Contour",
            folder_handle="Delete",
        )


def test_get_var_list():
    model = CloudSimulation(
        simulation_time_minutes=1,
        save_time_minutes=1,
        statistic_time_minutes=1,
        bacup_time_minutes=1,
        restore_backup=False,
    )

    model.run_model()
    model.cloud_analytics.list_var()
    model.initial_analytics.list_var()
    model.clean_model()


def test_get_var_with_invalid_file():
    cloud = FileStyle(
        chosen_file="Nube",
        output_data_path=".temp",
        img_path="img",
        img_option="Contour",
        folder_handle="Delete",
    )

    with pytest.raises(ValueError):
        cloud.get_var("Temperature", 0)


def test_get_var_with_valid_file():
    cloud = FileStyle(
        chosen_file="Nube",
        output_data_path=".temp",
        img_path="img",
        img_option="Contour",
        folder_handle="Delete",
    )
    check_cloud_var_list(cloud)


def check_cloud_var_list(cloud: FileStyle):
    for var in nube31_var_list:
        for time in range(
            len(_get_file_list(cloud.output_data_path, cloud.binary_regex))
        ):
            variable = cloud.get_var(var, time)
            assert variable is not None
            assert isinstance(variable, np.ndarray)


def test_show_var_dataframe_with_invalid_file():
    cloud = FileStyle(
        chosen_file="Nube",
        output_data_path="cloudmodel/tests/test_data",
        img_path="img",
        img_option="Contour",
        folder_handle="Delete",
    )
    var = cloud.get_var(nube31_var_list[0], 1)
    var_max_value = cloud._get_var_max_value_position(var)
    with pytest.raises(ValueError):
        cloud.show_var_dataframe(var, var_max_value, "a")


def test_show_var_dataframe_with_valid_file():
    cloud = FileStyle(
        chosen_file="Nube",
        output_data_path="cloudmodel/tests/test_data",
        img_path="img",
        img_option="Contour",
        folder_handle="Delete",
    )
    var = cloud.get_var(nube31_var_list[0], 1)
    var_max_value = cloud._get_var_max_value_position(var)
    cloud.show_var_dataframe(var, var_max_value, "x")
    cloud.show_var_dataframe(var, var_max_value, "y")
    cloud.show_var_dataframe(var, var_max_value, "z")


def test_plot_style_invalid_type():
    cloud = FileStyle(
        chosen_file="Nube",
        output_data_path="cloudmodel/tests/test_data",
        img_path="img",
        img_option="InvalidType",
        folder_handle="Delete",
    )
    var = cloud.get_var(nube31_var_list[0], 1)
    with pytest.raises(ValueError):
        cloud._plot_style(var)


def test_plot_style_initial_contour():
    initial = FileStyle(
        chosen_file="Inis",
        output_data_path="cloudmodel/tests/test_data",
        img_path="img",
        img_option="Contour",
        folder_handle="Delete",
    )
    var = initial.get_var(inis_var_list[0], 0)
    initial._plot_style(var)


def test_plot_style_contour():
    cloud = FileStyle(
        chosen_file="Nube",
        output_data_path="cloudmodel/tests/test_data",
        img_path="img",
        img_option="Contour",
        folder_handle="Delete",
    )
    var = cloud.get_var(nube31_var_list[0], 1)
    cloud._plot_style(var)


def test_plot_style_initial_image():
    initial = FileStyle(
        chosen_file="Inis",
        output_data_path="cloudmodel/tests/test_data",
        img_path="img",
        img_option="Image",
        folder_handle="Delete",
    )
    var = initial.get_var(inis_var_list[0], 0)
    initial._plot_style(var)


def test_plot_style_cloud_image():
    cloud = FileStyle(
        chosen_file="Nube",
        output_data_path="cloudmodel/tests/test_data",
        img_path="img",
        img_option="Image",
        folder_handle="Delete",
    )
    var = cloud.get_var(nube31_var_list[0], 1)
    cloud._plot_style(var)


def test_parse_status_img():
    cloud = FileStyle(
        chosen_file="Nube",
        output_data_path="cloudmodel/tests/test_data",
        img_path=".temp/img",
        img_option="Image",
        folder_handle="Delete",
    )
    cloud.parse_status_img()
    file_list = _get_file_list(cloud.output_data_path, cloud.binary_regex)
    assert len(file_list) != 0
    check_path("Delete", ".temp/img/")


def test_multi_var_img_valid_type():
    initials = FileStyle(
        chosen_file="Inis",
        output_data_path="cloudmodel/tests/test_data",
        img_path=".temp/img",
        img_option="Image",
        folder_handle="Delete",
    )
    var1 = inis_var_list[0]
    var2 = inis_var_list[1]
    initials.multi_var_img(var1, var2)
    check_path("Delete", ".temp/img/")


def test_multi_var_invalid_type():
    cloud = FileStyle(
        chosen_file="Nube",
        output_data_path="cloudmodel/tests/test_data",
        img_path=".temp/img",
        img_option="Contour",
        folder_handle="Delete",
    )
    var1 = nube31_var_list[0]
    var2 = nube31_var_list[1]
    with pytest.raises(ValueError):
        cloud.multi_var_img(var1, var2)
    check_path("Delete", ".temp/img/")


def test_check_path_delete():
    check_path("Delete", ".temp/Delete/", "test.txt")
    check_path("Delete", ".temp/Delete/", "test.txt")
    check_path("Delete", ".temp/Delete/", "test.txt")


def test_check_path_cancel():
    check_path("Cancel", ".temp/Cancel/", "test.txt")
    check_path("Cancel", ".temp/Cancel/", "test.txt")
    with pytest.raises(ValueError):
        check_path("Cancel", ".temp/Cancel/", "test.txt")


def test_check_path_invalid():
    with pytest.raises(ValueError):
        check_path("Invalid", ".temp/img/", "test.txt")
        check_path("Invalid", ".temp/img/", "test.txt")
        check_path("Invalid", ".temp/img/", "test.txt")
