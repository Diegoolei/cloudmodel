"""main file to run the cloud model."""

import os
import re

import numpy as np

from ..cloud_read import CloudSimulation, FileStyle, nube31_var_list

import pytest

from unittest.mock import patch



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
        output_data_path="test_data",
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
    assert True


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

    assert True


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

    assert isinstance(model.cloud_analytics.data_file, list)

    model.clean_model()

    assert True