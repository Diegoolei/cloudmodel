from constants import *


def MODEL31(ini):
    if ini:
        initial_conditions_setup()
    else:
        read_backup_files()

    print_var_status()
    open_file_units()

    print("----Simulation Start----")
    for tt in range(lt1 + 1):
        print("----Time Step: ", tt, "----")
        print("----Calculations without microphysics----")
        restart_vars()

        vapour_advection()

        dynamics_and_thermodynamics()

        save_dynamic_and_thermodynamic_variables()

        negative_value_correction()

        first_water_calculation()

        print("----Microphysics calculations----")

        microphysics()

        border_conditions()

        speed_and_pressure()

        borders_redefinitions_and_filters()

        negative_vapour_correction()

        print_final_status()

        save_final_status()

        print("----End of time step", tt, "from ", lt1, " ---- \n")
    print("----Simulation End----")


def initial_conditions_setup():
    condi()


def condi():
    print("*--initial_conditions_setup--")


def read_backup_files():
    print("--read_backup_files--")


def print_var_status():
    print("--print_var_status--")


def open_file_units():
    print("--open_file_units--")


def restart_vars():
    print("--restart_vars--")


def vapour_advection():
    print("--vapour_advection--")


def dynamics_and_thermodynamics():
    print("--dynamics_and_thermodynamics--")
    print("Dynamics")
    print("Thermodynamics")


def save_dynamic_and_thermodynamic_variables():
    print("--save_dynamic_and_thermodynamic_variables--")


def negative_value_correction():
    print("--negative_value_correction--")


def first_water_calculation():
    print("--first_water_calculation--")


def microphysics():
    print("--microphysics--")
    print("407 to 609")


def border_conditions():
    print("--border_conditions--")
    print("floor and ceiling")
    print("lateral")


def speed_and_pressure():
    velpre()


def velpre():
    print("*--speed_and_pressure--")


def borders_redefinitions_and_filters():
    print("--borders_redefinitions_and_filters--")
    print("Redefinition")
    print("frame floor ceiling")
    print("frame lateral")
    vapour_filter()


def vapour_filter():
    print("--vapour_filter--")


def negative_vapour_correction():
    print("--negative_vapour_correction--")


def print_final_status():
    print("--print_final_status--")


def save_final_status():
    print("--save_final_status--")
    estad03()
    posnub02()
    corrinu2()
    graba320()
    graba120()


def estad03():
    print("*interesting quantities calculation")


def posnub02():
    print("*cloud position (X, Y)")


def corrinu2():
    print("*cloud displacement")


def graba320():
    print("*3D variables saving")


def graba120():
    print("*variable saving")


MODEL31(True)
