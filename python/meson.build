project('cloudmodel', 'c',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

add_languages('fortran', native: false)

meson_path = meson.current_build_dir()

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()


incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()


incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

# fpm?
fflags = '-g -fPIC -funroll-loops -fstack-arrays -Ofast -frepack-arrays -faggressive-function-elimination -fopenmp'
cflags = '-fPIC'

fpm_run = run_command(
    'fpm',
    'build',
    '--profile',
    'release',
    '--flag',
    fflags,
    '--c-flag',
    cflags
).stdout()

message('@0@'.format(fpm_run))

ls_meson_command = run_command('ls', meson_path).stdout()
message('meson_dir: @0@'.format(ls_meson_command))

# f2py?
source_path = meson.source_root() / '..' / 'interface' / 'cloud_c_interface.f90'
cloud_c_f90 = files(source_path)

cloud_source = custom_target('cloud_compiled.c',
  input: cloud_c_f90,
  output: ['cloud_compiledmodule.c', 'cloud_compiled-f2pywrappers2.f90'],
  command: [
    py,
    '-m',
    'numpy.f2py',
    '@INPUT@',
    '-m',
    'cloud_compiled',
    '--lower'
  ]
)

# Install Python sources
f_sources = ['cloudmodel/interface/__init__.py']
py.install_sources(f_sources, subdir:'cloudmodel/interface')

main_sources = ['cloudmodel/__init__.py']
py.install_sources(main_sources, subdir:'cloudmodel')

# Extension fortran module
inc_np = include_directories(incdir_numpy, incdir_f2py, '..' / '..' / 'build' / 'gfortran_AFF0773A2CF7CD9C')

py.extension_module('cloud_compiled',
  [cloud_c_f90, cloud_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: [inc_np],
  dependencies : py_dep,
  link_args: ['-L' + '..' / '..' / 'build' / 'gfortran_AFF0773A2CF7CD9C' / 'cloudsim', '-lcloudsim'],
  subdir: 'cloudmodel' / 'interface',
  install : true,
)