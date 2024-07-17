# -L path absoluto a carpeta que contiene libclodusim.a
# -I path absoluto a carpeta que contiene *.mod
# -lcloudsim linkea la lib
gfortran:
	fpm clean --all
	fpm install --profile release --prefix build/python
	f2py -m interface \
		-L$(PWD)/../../../build/python/lib/ \
		-I$(PWD)/../../../build/python/include/ \
		-c cloud_c_interface.f90 -lcloudsim


nvfortran:
	#fpm clean --all
	fpm build --compiler "/opt/nvidia/hpc_sdk/Linux_x86_64/2024/compilers/bin/nvfortran" --flag "-O3 -cuda"
	f2py -m interface \
		-L$(PWD)/../build/nvfortran_C1B52A2EA9B3E302/cloudsim/ \
		-I$(PWD)/../build/nvfortran_72B09486B8F05AD4/ \
		-c cloud_c_interface.f90 -lcloudsim
	# fpm install --profile release --prefix build/python