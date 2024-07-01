# Fortran77-Cloud-Model

Animation functionality requirements:

```console
sudo apt install ffmpeg
```

To show animations using Windows Subsystem for Lynux (WSL) follow this [steps](https://stackoverflow.com/questions/43397162/show-matplotlib-plots-and-other-gui-in-ubuntu-wsl1-wsl2).

- Install [Xming X Server for Windows](https://sourceforge.net/projects/xming/)
- Run in WSL2 terminal: `sudo apt-get install python-tk`

Compilation and execution fortran code requirements:

Linux:

```console
sudo apt install gfortran
```

For MacOS:

```console
brew install gfortran
xcode-select --install
```

## F2Py

### Compile for F2Py

```console
cd interface
make
```

### Execute with Python

```console
cd ..
python3 main.py
```

## nvfortran

### Installation

First of all, Install the appropriate Nvidia drivers for your system.
Install the [Nvidia CUDA toolkit.](https://developer.nvidia.com/cuda-toolkit)
Install the [Nvidia HPC SDK](https://developer.nvidia.com/nvidia-hpc-sdk-downloads) (The installation path is usually /opt/nvidia/hpc_sdk/Linux_x86_64/(version)/compilers/bin, add it to your PATH).

### Execution

fpm run --compiler "/opt/nvidia/hpc_sdk/Linux_x86_64/(version)/compilers/bin/nvfortran" --flag "-O3 -cuda"
