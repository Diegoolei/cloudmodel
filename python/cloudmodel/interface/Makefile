# -L path absoluto a carpeta que contiene libclodusim.a
# -I path absoluto a carpeta que contiene *.mod
# -lcloudsim linkea la lib
gfortran:
	fpm clean --all
	fpm install --profile release --prefix build/python
	f2py -m cloud_compiled \
		-L$(PWD)/../../../build/python/lib/libclodusim.a \
		-I$(PWD)/../../../build/python/include \
		-c $(PWD)/../../../interface/cloud_c_interface.f90 \
		../../build/python/lib/libclodusim.a --backend meson


nvfortran:
	#fpm clean --all
	fpm install --prefix build/python\
	   --compiler \
	   "/opt/nvidia/hpc_sdk/Linux_x86_64/2024/compilers/bin/nvfortran"\
	   --flag "-O3 -cuda"
	f2py -m cloud_compiled \
		-L$(PWD)/../build/python/lib/libclodusim.a \
		-I$(PWD)/../build/python/include/ \
		-c $(PWD)/../../../interface/cloud_c_interface.f90 \
		../../build/python/lib/libclodusim.a --backend meson
	# -lcloudsim