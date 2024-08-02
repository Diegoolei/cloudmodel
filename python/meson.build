project('cloudmodel', 'c',
  version : '0.1.0',
  license: 'MIT',
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
fpm_build_dir = 'fpm_build'

run_command(
    'fpm',
    'install',
    '--profile',
    'release',
    '--flag',
    fflags,
    '--c-flag',
    cflags,
    '--prefix',
    fpm_build_dir
)

# f2py?
source_path = meson.project_source_root() / '..' / 'interface' / 'cloud_c_interface.f90'
cloud_c_interface = files(source_path)

cloud_c_interface_source = custom_target('cloud_c_interface_compiled.c',
  input: cloud_c_interface,
  output: ['cloud_c_interface_compiledmodule.c', 'cloud_c_interface_compiled-f2pywrappers2.f90'],
  command: [
    py,
    '-m',
    'numpy.f2py',
    '@INPUT@',
    '-m',
    'cloud_c_interface_compiled',
    '--lower'
  ]
)

# Install Python sources
f_sources = ['cloudmodel/interface/__init__.py']
py.install_sources(f_sources, subdir:'cloudmodel/interface')

main_sources = ['cloudmodel/__init__.py']
py.install_sources(main_sources, subdir:'cloudmodel')

# Extension fortran module
inc_np = include_directories(incdir_numpy, incdir_f2py)
inc_ftn = include_directories( '..' / '..' / 'fpm_build' / 'include' )

py.extension_module('cloudmodel_compiled',
  [cloud_c_interface, cloud_c_interface_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: [inc_np, inc_ftn],
  dependencies : py_dep,
  link_args: ['-L' + '../../fpm_build' / 'lib', '-lcloudsim'],
  subdir: 'cloudmodel' / 'interface',
  install : true,
)