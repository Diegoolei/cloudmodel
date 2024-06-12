# -L path absoluto a carpeta que contiene libclodusim.a
# -I path absoluto a carpeta que contiene *.mod
# -lcloudsim linkea la lib
python:
	f2py -m nubepy \
		-L$(PWD)/../build/gfortran_50F62D7499E64B65/cloudsim/ \
		-I$(PWD)/../build/gfortran_AC789FC96BC46CEE/ \
		-c nube_c.f90 -lcloudsim
	# fpm install --profile release --prefix build/python
