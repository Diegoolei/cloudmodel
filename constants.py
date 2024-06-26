total_time = 1
nx1 = 50
nz1 = 45
biased_nx1 = nx1 + 7
nube31_biased_nz1 = nz1 + 5
inis_biased_nz1 = nz1 + 7

nube31_var_list = [
    "original_x_velocity",
    "original_y_velocity",
    "original_z_velocity",
    "potential_temperature_base",
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
    "Den0",
    "Temp0",
    "Tita0",
    "Pres00",
    "Qvap0",
    "cc2",
    "aer0",
    "UU",
    "VV",
]

plot_center = 10  # center of the plot
