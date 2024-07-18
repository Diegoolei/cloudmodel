"""Setup file for the python package."""

from setuptools import setup, find_packages

packages = find_packages()
print(packages)
setup(
    packages=packages,
    include_package_data=True,
    package_data={
        "": ["*"],
    },
)
