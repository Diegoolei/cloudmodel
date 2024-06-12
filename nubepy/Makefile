# -L path absoluto a carpeta que contiene libclodusim.a
# -I path absoluto a carpeta que contiene *.mod
# -lcloudsim linkea la lib
python:
	f2py -m nubepy \
		-L$(PWD)/../build/gfortran_63E63B2712D76620/cloudsim/ \
		-I$(PWD)/../build/gfortran_E62477DA9516D2A9/ \
		-c nube_c.f90 -lcloudsim
	# fpm install --profile release --prefix build/python
