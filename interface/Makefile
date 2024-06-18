# -L path absoluto a carpeta que contiene libclodusim.a
# -I path absoluto a carpeta que contiene *.mod
# -lcloudsim linkea la lib
python:
	fpm clean --all
	fpm build --profile release
	f2py -m interface \
		-L$(PWD)/../build/gfortran_63E63B2712D76620/cloudsim/ \
		-I$(PWD)/../build/gfortran_E62477DA9516D2A9/ \
		-c cloud_c_interface.f90 -lcloudsim
	# fpm install --profile release --prefix build/python
