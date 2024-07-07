from pathlib import Path
from setuptools import setup, find_packages

setup(
    name="cloudmodel",
    version="0.2.1",
    description="A simple, parameterized cloud model",
    long_description = Path("README.md").read_text(),
    long_description_content_type="text/markdown",
    url="https://github.com/Diegoolei/Fortran77-Cloud-Model",
    author="Diego E. Oleiarz",
    author_email= "diego.oleiarz@mi.unc.edu.ar",
    license="MIT",
    project_urls={
        "Source": "https://github.com/Diegoolei/Fortran77-Cloud-Model",
        "Documentation": "https://diegoolei.github.io/Fortran77-Cloud-Model/"
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Fortran",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6, <4",
    # install_requires=[ "numpy", "matplotlib", "pandas", "fpm", "scipy"],
    packages = find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "cloudmodel=cloudmodel.__main__:main",
        ],
    },
)