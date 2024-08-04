"""Constants for the cloud model."""


# Var list for the cloud model
# All var lists must be in the same order as the data in the binary files

nube31_var_list = [
    "original_x_velocity",
    "original_y_velocity",
    "original_z_velocity",
    "theta_base",
    "original_pressure",
    "vapor_base",
    "drop_base",
    "rain_base",
    "crystal_base",
    "snow_base",
    "hail_base",
    "aerosol_base",
]

inis_var_list = [
    "air_density_z_initial",
    "temperature_z_initial",
    "theta_z_initial",
    "Pres00",
    "vapor_z_initial",
    "cc2",
    "aerosol_z_initial",
    "u_z_initial",
    "v_z_initial",
]

# ?
n = 45 # 54

# Dimensions for the cloud model
dx1= 300
nz1 = n

# Constants for the cloud model
G = 9.8
Rd=287.04
Rv=461.05
Kapa=0.2857
T0=273.15
P00=101300.
Lvl0=2.500e6
Lsl0=79.7 #calor latente fusion a 0C en kilocalorias
Vis0 = 1.718e-5
rhogra = 500.
Av0 = 1455.
Vtnie0 = .5
dx = dx1 # 300

import numpy as np

zeta_offset = 3 #extension del dominio arriba y abajo para corregir problemas de memoria del Fortran
zeta = np.arange(0, (n+1) * dx1, dx1)
Tmin = -70 #en C
Tmax = 30 #en C

#Calulo de las velocidades horizontales
zeta_p = np.array([0, 500, 2000, 9000, 13500.])
U_p = np.array([0., 0., 4., -6., -5 ])
V_p = np.array([0., 0., 0., 3., -.5])

#Calculo de la Temperatura
T_0= 298.15
dT_p = np.array([-9e-3, -9e-3, -7e-3, -7e-3, 0, 0, (50./27)*1e-3 ])
zeta_p = np.array([0, 2000, 5500., 9000., 11000, 12000., 13500.])

rel1_p = np.array([.55, .6, .6, .35, .1, .05666667])
zeta_p = np.array([0, 500, 1500., 4000., 7000., 13500.])

nx1 = 50
biased_nx1 = nx1 + 7
nube31_biased_nz1 = nz1 + 5
inis_biased_nz1 = nz1 + 7