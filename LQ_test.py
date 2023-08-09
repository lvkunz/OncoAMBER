import sys
sys.path.insert(0, '/PHShome/lk001/.conda/envs/amberenv/lib/python3.9/site-packages') #cluster
import amber
import numpy as np
import random
import time
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

print('Current working directory:', os.getcwd())
#print the directory of amber
print('Amber directory:', amber.__file__)
print('Amber version:', amber.__version__)

config_file = 'CONFIG_LQ'

config = amber.Config(config_file)

print('Config file')
print(config)
print(config.half_length_world)

#set seed for reproducibility
start_time = time.time()
DPI = 100

seed = config.seed
np.random.seed(seed)
random.seed(seed)


def LQ_model(dose, alpha, beta):
    dose = np.array(dose)
    return np.exp(-alpha * dose - beta * dose * dose)

def run_experiment(dose,repair_rate, radiosensitivity):
    list_of_cells = []
    for i in range(10000):
        cell = amber.TumorCell(config.radius_tumor_cells,
                                       cycle_hours=1000000.0,
                                       cycle_std=config.doubling_time_sd,
                                       intra_radiosensitivity=radiosensitivity,
                                       o2_to_vitality_factor=config.o2_to_vitality_factor,
                                       VEGF_threshold=config.VEGF_production_threshold,
                                       VEGF_rate=config.VEGF_production_per_cell
                                       )
        cell.oxygen = 1.0
        list_of_cells.append(cell)

    voxel = amber.Voxel([0,0,0], 10, 1.0, list_of_cells_in=list_of_cells, n_capillaries = 1, voxel_number = 1)

    dt = 1.0 #hours

    celldeath = amber.CellDeath(config, 'cell_death', dt,
                                            apoptosis_threshold=config.vitality_apoptosis_threshold,
                                            apoptosis_probability=config.probability_apoptosis,
                                            necrosis_threshold=config.vitality_necrosis_threshold,
                                            necrosis_probability=config.probability_necrosis,
                                            necrosis_removal_probability=config.probability_necrosis_removal,
                                            necrosis_damage_coeff=0.5,
                                            apoptosis_damage_coeff=1.0)

    cellaging = amber.CellAging(config, 'cell_aging', dt, repair_per_hour=repair_rate)

    def irradiation(voxel, dose):
            scaled_dose = dose
            if len(voxel.list_of_cells) > 0:
                print('Scaled dose: ', scaled_dose)
                for cell in voxel.list_of_cells:
                    # assume all cells get damaged the same way
                    damage = scaled_dose * cell.radiosensitivity()  # compute the damage
                    cell.damage += damage
                    cell.damage = min(cell.damage, 1.0)

    n_cells = []
    n_necrotic_cells = []


    for t in range(500):
        celldeath(voxel)
        if t == 25:
            irradiation(voxel, dose)
            #show histogram of cell damage
        cellaging(voxel)
        print('t = ', t, 'n_cells = ', len(voxel.list_of_cells))
        n_cells.append(len(voxel.list_of_cells))
        n_necrotic_cells.append(len(voxel.list_of_necrotic_cells))

    # fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=200)
    # ax.plot(n_cells, label='Cells')
    # ax.plot(n_necrotic_cells, label='Necrotic cells')
    # ax.set_xlabel('Time (hours)')
    # ax.set_ylabel('Number of cells')
    # ax.legend()
    # plt.show()

    return n_cells[-1]

repair_rates = [ 0.03 ]
doses = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3,1.4,1.5,1.6,1.7,1.8,1.9, 2.0]
# doses = [0.0, 0.5, 1.0, 1.5, 2.0]
radio_sensitivities = [0.1, 0.3, 0.5]
# doses = [0.0, 1.0]
#for each repair rate, store the number of surviving cells for each dose


for repair in repair_rates:
    surviving_cells_array = []
    for i in doses:
        print('Dose = ', i)
        surviving_cells = run_experiment(i, repair, 0.3)
        print('----------------------------------')
        surviving_cells_array.append(surviving_cells)
        np.save('LQ_tests/doses_array_repair='+str(repair)+'.npy', doses)
        np.save('LQ_tests/surviving_cells_array_repair='+str(repair)+'.npy', surviving_cells_array)

for radiosensi in radio_sensitivities
    surviving_cells_array = []
    for i in doses:
        print('Dose = ', i)
        surviving_cells = run_experiment(i, 0.03, radiosensi)
        print('----------------------------------')
        surviving_cells_array.append(surviving_cells)
        np.save('LQ_tests/doses_array_radiosensi='+str(radiosensi)+'.npy', doses)
        np.save('LQ_tests/surviving_cells_array_radiosensi='+str(radiosensi)+'.npy', surviving_cells_array)

for repair in repair_rates:

    fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=200)
    popt, pcov = curve_fit(LQ_model, doses, surviving_cells_array, p0=[0.01,0.01], maxfev=100000, bounds=(0, [1.0, 1.0])
    doses_smooth = np.linspace(0, 2, 100)
    ax.plot(doses, surviving_cells_array, 'o', label='Repair rate = ' + str(repair))
    ax.plot(doses_smooth, LQ_model(doses, *popt), label='LQ model, alpha = ' + str(popt[0]) + ', beta = ' + str(popt[1]))
    ax.set_yscale('log')
    ax.set_xlabel('Dose (Gy)')
    ax.set_ylabel('Number of surviving cells')
    ax.legend()
    plt.show()
    print('----------------------------------')

# fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=200)
# for radiosensi in radio_sensitivities:
#     surviving_cells_array = []
#     for i in doses:
#         print('Dose = ', i)
#         surviving_cells = run_experiment(i, 0.03, radiosensi)
#         print('----------------------------------')
#         surviving_cells_array.append(surviving_cells)
#     ax.plot(doses, surviving_cells_array, label='Radiosensitivity = ' + str(radiosensi))
#
# ax.set_yscale('log')
# ax.set_xlabel('Dose (Gy)')
# ax.set_ylabel('Number of surviving cells')
# ax.legend()
# plt.show()
# print('----------------------------------')




