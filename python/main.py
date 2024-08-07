from cloudmodel.cloud_read import CloudSimulation

import numpy as np
from matplotlib import pyplot as plt
import os

def generate_cut_image():
    dirout = '../cortes/'
    directory = 'Data/z_initials_param/cortes/'
    os.chdir(directory)

    ny = 45
    nx = 50
    part = ['GO', 'LL', 'CR', 'NI','GR' ]
    pos = '31y'

    vec = np.ones([ny, nx, len(part), 1])

    for i, nom in enumerate(part):
        nombre = nom + pos + "17"
        print(i, 0, nombre)
        vec[:, :, i, 0]  = np.genfromtxt(nombre)
    colores = ['r', 'tab:orange', 'c', 'b', 'g']
    plt.figure(figsize=[10,10])
    plt.ylim(0,35)
    plt.xlim(10,40)
    for i in range(len(part)):
        plt.contour(vec[: ,:, i, 0], [10, 500], colors= colores[i])

    plt.savefig(f'{dirout}PYTHON.png')
    plt.show()


def main():
    z_initials_param = CloudSimulation(
        simulation_time_minutes=45,
        save_time_minutes=3,
        statistic_time_minutes=3,
        bacup_time_minutes=3,
        restore_backup=False,
        directory="Data/z_initials_param",
    )
    z_initials_param.run_model()
    #z_initials_param.load_model()
    z_initials_param.cloud_analytics.parse_status_img()
    generate_cut_image()

main()
