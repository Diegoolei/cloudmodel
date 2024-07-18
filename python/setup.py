"""Setup file for the python package."""
from setuptools import setup, find_packages

setup(name = "cloudmodel", 
      version="0.2.15",
      packages = find_packages(),
      include_package_data=True,
      package_dir={"cloudmodel": ""},
      package_data={
        '': ['*'],
    })
