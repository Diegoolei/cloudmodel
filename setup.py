"""Setup file for the python package."""
from setuptools import setup, find_packages

setup(packages = find_packages("python/*"),
      include_package_data=True,
      package_dir={"": "python"},)