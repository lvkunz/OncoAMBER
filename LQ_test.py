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

def LQ_wo_exp(dose, alpha, beta):
    dose = np.array(dose)
    return -alpha * dose - beta * dose * dose


def run_experiment(dose,repair_rate, radiosensitivity,n0, n_fractions):
    list_of_cells = []
    for i in range(n0):
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
        print('Irradiation')
        scaled_dose = dose
        #add small noise to dose
        for cell in voxel.list_of_cells:
            # assume all cells get damaged the same way
            damage = scaled_dose * cell.radiosensitivity()  # compute the damage
            cell.damage += (damage+ np.random.normal(0, 0.01))  # add the damage to the cell
            cell.damage = min(cell.damage, 1.0)

    n_cells = []
    n_necrotic_cells = []

    first_irradiation = 24
    time_between_fractions = 24
    fractions = [first_irradiation + i * time_between_fractions for i in range(n_fractions)]
    fractions.append(1000000)

    dose_f = dose / n_fractions

    i = 0
    for t in range(500):
        celldeath(voxel)
        if t == fractions[i]:
            irradiation(voxel, dose_f)
            i += 1
            #show histogram of cell damage
        cellaging(voxel)
        # print('t = ', t, 'n_cells = ', len(voxel.list_of_cells))
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

write = False
fit = True

repair_rates = [1, 3, 5]
doses = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]

# doses = [0.0, 0.5, 1.0, 1.5, 2.0]
radio_sensitivities = []
# doses = [0.0, 0.5, 1.0, 1.5, 2.0]
#for each repair rate, store the number of surviving cells for each dose
n0 = 5000

if write:
    for repair in repair_rates:
        surviving_cells_array = []
        for i in doses:
            print('Dose = ', i)
            surviving_cells = run_experiment(i, 0.03, 0.3, n0, n_fractions=repair)
            print('----------------------------------')
            surviving_cells_array.append(surviving_cells/n0)
            np.save('LQ_tests/doses_array_repair='+str(repair)+'.npy', doses)
            np.save('LQ_tests/surviving_cells_array_repair='+str(repair)+'.npy', surviving_cells_array)

    for radiosensi in radio_sensitivities:
        surviving_cells_array = []
        for i in doses:
            print('Dose = ', i)
            surviving_cells = run_experiment(i, 0.03, radiosensi,n0,n_fractions=1)
            print('----------------------------------')
            surviving_cells_array.append(surviving_cells/n0)
            np.save('LQ_tests/doses_array_radiosensi='+str(radiosensi)+'.npy', doses)
            np.save('LQ_tests/surviving_cells_array_radiosensi='+str(radiosensi)+'.npy', surviving_cells_array)

fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=200)
for repair in repair_rates:

    surviving_cells_array = np.load('LQ_tests/surviving_cells_array_repair='+str(repair)+'.npy')
    doses = np.load('LQ_tests/doses_array_repair='+str(repair)+'.npy')

    finite_surviving_cells_array = []
    finite_doses = []
    for i in range(len(surviving_cells_array)):
        if surviving_cells_array[i] != 0:
            finite_surviving_cells_array.append(surviving_cells_array[i])
            finite_doses.append(doses[i])

    finite_surviving_cells_array = np.log(finite_surviving_cells_array)
    popt, pcov = curve_fit(LQ_wo_exp, finite_doses, finite_surviving_cells_array, p0=[1,5], maxfev=100000, bounds=([0.0,0],[100,100]))
    alpha_beta_ratio = popt[0]/popt[1]
    print('alpha/beta ratio = ', alpha_beta_ratio)
    doses_smooth = np.linspace(0, 2, 100)
    color = next(ax._get_lines.prop_cycler)['color']
    ax.plot(finite_doses, finite_surviving_cells_array, 'x', label='Number of fractions = ' + str(repair), color=color)
    if fit: ax.plot(doses_smooth, LQ_wo_exp(doses_smooth, *popt), color = color)#, label='LQ model, alpha = ' + str(round(popt[0],3)) + ', beta = ' + str(round(popt[1],3)))

ax.set_xlabel('Dose (arb. units)')
ax.set_ylabel('Fraction of surviving cells')
#put the ticks as 10^x
ax.set_ylim(-8, 1)
ax.set_yticks([0, -1, -2, -3, -4, -5, -6, -7, -8])
ax.set_yticklabels(['$10^0$', '$10^{-1}$', '$10^{-2}$', '$10^{-3}$', '$10^{-4}$', '$10^{-5}$', '$10^{-6}$', '$10^{-7}$', '$10^{-8}$'])
ax.legend()
ax.grid()
plt.savefig('LQ_tests/repair_rate.png')
plt.show()
print('----------------------------------')

fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=200)
for radiosensi in radio_sensitivities:
    surviving_cells_array = np.load('LQ_tests/surviving_cells_array_radiosensi=' + str(radiosensi) + '.npy')
    doses = np.load('LQ_tests/doses_array_radiosensi=' + str(radiosensi) + '.npy')

    finite_surviving_cells_array = []
    finite_doses = []
    for i in range(len(surviving_cells_array)):
        if surviving_cells_array[i] != 0:
            finite_surviving_cells_array.append(surviving_cells_array[i])
            finite_doses.append(doses[i])

    finite_surviving_cells_array = np.log(finite_surviving_cells_array)
    popt, pcov = curve_fit(LQ_wo_exp, finite_doses, finite_surviving_cells_array, method = 'trf', p0=[100, 100], maxfev=100000, bounds=([0.0, 0.0], [100, 100]))
    alpha_beta_ratio = popt[0] / popt[1]
    print('alpha/beta ratio = ', alpha_beta_ratio)
    doses_smooth = np.linspace(0, 2, 100)
    color = next(ax._get_lines.prop_cycler)['color']
    ax.plot(finite_doses, finite_surviving_cells_array, 'x', color=color, label='Radiosensitivity = ' + str(radiosensi))
    if fit: ax.plot(doses_smooth, LQ_wo_exp(doses_smooth, *popt),  color=color)#, label='LQ model, alpha = ' + str(round(popt[0],3)) + ', beta = ' + str(round(popt[1],3)))

ax.set_xlabel('Dose (arb. units)')
ax.set_ylabel('Fraction of surviving cells')
ax.set_ylim(-8, 1)
ax.set_yticks([0, -1, -2, -3, -4, -5, -6, -7, -8])
ax.set_yticklabels(['$10^0$', '$10^{-1}$', '$10^{-2}$', '$10^{-3}$', '$10^{-4}$', '$10^{-5}$', '$10^{-6}$', '$10^{-7}$', '$10^{-8}$'])
ax.legend()
ax.grid()
plt.savefig('LQ_tests/radiosensitivity.png')
plt.show()
print('----------------------------------')






