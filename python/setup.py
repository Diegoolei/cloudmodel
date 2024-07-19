"""Setup file for the python package."""

from setuptools import find_packages, setup

packages = find_packages()
print(packages)
setup(
    name="cloudmodel",
    version="0.2.16",
    packages=packages,
    include_package_data=True,
    package_data={
        "": ["*"],
    },
)
