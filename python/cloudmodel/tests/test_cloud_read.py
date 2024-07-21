import pytest
from cloudmodel.cloud_read import CloudSimulation, FileStyle, get_file_list, data_comparison, check_path

def test_CloudSimulation_run_model():
    cloud_sim = CloudSimulation()
    cloud_sim.run_model()
    assert cloud_sim.initial_analytics is not None
    assert cloud_sim.cloud_analytics is not None

def test_FileStyle_get_data():
    file_style = FileStyle()
    file_style.get_data()
    assert len(file_style.data_file) > 0

def test_get_file_list():
    data_path = "outputdata/"
    binary_regex = r".sal$"
    file_list = get_file_list(data_path, binary_regex)
    assert len(file_list) > 0

def test_data_comparison():
    original_data = FileStyle()
    cmp_data = FileStyle()
    result = data_comparison(original_data, cmp_data)
    assert result == True

def test_check_path():
    folder_handle = "Ignore"
    path = "outputdata/"
    selected_file_name = ""
    check_path(folder_handle, path, selected_file_name)
    assert True  # No assertion, just checking if the function runs without errors