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

### Compile with F2Py

```console
cd interface
make
```

### Execute with Python

```console
cd ..
python3 main.py
```
