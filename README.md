# CloudModel

This repository hosts the re-engineered cloud simulation model, now implemented in Fortran 90 and wrapped for Python integration. Below are the instructions for setting up and running the model, as well as additional requirements for animation functionality and compilation with different Fortran compilers.

## Table of Contents

- [Python](#python)

-- [Installation](#instalation)
-- [Usage](#usage)
-- [Extra Functionality Requirements](#extra-functionality-requirements)

- [Fortran](#fortran-for-devs)

-- [Fortran Only](#fortran-only-compiling)
-- [F2Py](#f2py-integration)

- [Contributing](#contributing)
- [License](#license)

## Python

### Instalation

To install as a pip package:

```console
    pip install -i https://test.pypi.org/simple/ cloudmodel
```

### Usage

To run the project, use the following command:

### Features

<https://packaging.python.org/en/latest/guides/making-a-pypi-friendly-readme/>

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
