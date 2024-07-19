"""Constants for the cloud model."""

nx1 = 50
nz1 = 45
biased_nx1 = nx1 + 7
nube31_biased_nz1 = nz1 + 5
inis_biased_nz1 = nz1 + 7

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

plot_center = 10  # center of the plot
