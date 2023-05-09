import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob

import os
import glob

def get_repository_names(path):
    """
    Returns a list of repository names in the given path.
    """
    repo_paths = glob.glob(os.path.join(path, "*/*/.git"))
    repo_names = []
    for repo_path in repo_paths:
        repo_names.append(os.path.basename(os.path.dirname(repo_path)))
    return repo_names

def func(x, a, b, c):
    return a * (np.exp(b * x)) + c

number_cells_list = []
tumor_size_list = []
times_list = []

# dt = [10, 15, 20, 25, 30, 50, 2, 3, 4, 5]
#
# repo = 'output/CONFIG_dt_convergence_example.py_20230505_lk001_Linux/'
#
# paths = ['example.py_CONFIG_dt_convergence_dt5_iter0/DataOutput/',
#         'example.py_CONFIG_dt_convergence_dt10_iter1/DataOutput/',
#         'example.py_CONFIG_dt_convergence_dt15_iter2/DataOutput/',
#         'example.py_CONFIG_dt_convergence_dt20_iter3/DataOutput/',
#         'example.py_CONFIG_dt_convergence_dt25_iter4/DataOutput/',
#         'example.py_CONFIG_dt_convergence_dt30_iter5/DataOutput/',
#         'example.py_CONFIG_dt_convergence_dt50_iter6/DataOutput/',
#         'example.py_CONFIG_dt_convergence_dt2_iter7/DataOutput/',
#         'example.py_CONFIG_dt_convergence_dt3_iter8/DataOutput/',
#         'example.py_CONFIG_dt_convergence_dt4_iter9/DataOutput/']



CONFIG_file = 'CONFIG_dt_convergence_cycling0'
CONFIG_file = 'CONFIG_dt_convergence_no_normal_cell'

repo = 'output/'+ CONFIG_file + '_example.py_20230508_lk001_Linux'

#all repositories in repo:

dt = [2, 3, 4, 5, 7, 10, 15]

paths = []

for i, d in enumerate(dt):
    path = f'example.py_{CONFIG_file}_dt{d}_iter{i}/DataOutput/'
    full_path = f'{repo}/{path}'
    paths.append(full_path)


for i in range(len(paths)):
    for path in paths:
        number_cells = np.load(f'{path}number_tumor_cells.npy', allow_pickle=True)
        tumor_size = np.load(f'{path}tumor_size.npy', allow_pickle=True)
        times = np.load(f'{path}times.npy', allow_pickle=True)
        number_cells_list.append(number_cells)
        tumor_size_list.append(tumor_size)
        times_list.append(times)


# Fit the data to an exponential curve for each simulation and get the doubling time
doubling_times_number_cells = []
doubling_times_tumor_size = []
dpi = 300
fig, axes = plt.subplots(2, 1, figsize=(8, 10), dpi=dpi)
for i in range(len(paths)):
    print(paths[i])
    # Fit number of cells

    popt, pcov = curve_fit(func, times_list[i], number_cells_list[i], p0=(3000, 3e-3, 0))
    print(popt)
    color = axes[0].plot(times_list[i], number_cells_list[i], 'o', markersize=1, alpha=0.5, label=f'dt =' + str(dt[i]))[0].get_color()
    # axes[0].plot(times_list[i], func(times_list[i], *popt), '-', color=color, label='fit dt =' + str(dt[i]))
    doubling_time = np.log(2)/popt[1]
    print('Doubling time (Number of Cells):', doubling_time)
    doubling_times_number_cells.append(doubling_time)

    # Fit tumor size
    popt, pcov = curve_fit(func, times_list[i], tumor_size_list[i], p0=(1, 0.003, 0))
    print(popt)
    axes[1].plot(times_list[i], tumor_size_list[i], 'o', color = color, markersize = 1, alpha=0.5, label=f'dt =' + str(dt[i]))
    # axes[1].plot(times_list[i], func(times_list[i], *popt), '-', color=color, label='fit dt =' + str(dt[i]))
    doubling_time = np.log(2)/popt[1]
    doubling_times_tumor_size.append(doubling_time)

axes[0].set_title('Number of Cells Evolution')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Number of Cells')
axes[0].set_xlim(0, 200)
axes[0].set_ylim(0, 2000)
axes[0].grid(True)
axes[0].legend()

axes[1].set_title('Tumor Volume Evolution')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Tumor Volume [mm^3]')
axes[1].set_xlim(0, 200)
axes[1].set_ylim(0, 3)
axes[1].grid(True)
axes[1].legend()

plt.tight_layout()
plt.savefig('growth_evolution.png')
plt.show()

print('Doubling times (Number of Cells):', doubling_times_number_cells)
print('Doubling times (Tumor Size):', doubling_times_tumor_size)

plt.plot(dt, doubling_times_number_cells, 'bo', label='Cells doubling time')
plt.xlabel('Time step')
plt.ylabel('Doubling time [days]')
plt.title('Doubling time vs. Time step')
# plt.yscale('log')  # set y-axis to logarithmic scale
plt.legend()
plt.grid(True)
plt.savefig('doubling_time.png', dpi=300)
plt.show()

plt.plot(dt, doubling_times_tumor_size, 'ro', label='Tumor volume doubling time')
plt.xlabel('Time step')
plt.ylabel('Doubling time [days]')
plt.title('Doubling time vs. Time step')
# plt.yscale('log')  # set y-axis to logarithmic scale
plt.legend()
plt.grid(True)
plt.savefig('doubling_time_tumor_size.png', dpi=300)
plt.show()

# #plot the data and the fitted curve
# fig, axs = plt.subplots(2, 1, figsize=(14, 10))
# fig.set_dpi(300)
# # Plot number of tumor cells
# axs[0].plot(times, number_cells, 'o', linewidth=2, alpha=0.8)
# axs[0].plot(times, func(times, *popt1), 'k--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt1))
# axs[0].set_xlabel('Time (days)')
# axs[0].set_ylabel('Number of tumor cells')
# axs[0].set_title('Number of tumor cells over time')
# axs[0].set_facecolor('whitesmoke')
# axs[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
#
# # Plot tumor size
# axs[1].plot(times, tumor_size, 'o', linewidth=2, alpha=0.8)
# axs[1].plot(times, func(times, *popt2), 'k--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt2))
# axs[1].set_xlabel('Time (days)')
# axs[1].set_ylabel('Tumor volume [mm^3]')
# axs[1].set_title('Tumor volume over time')
# axs[1].set_facecolor('whitesmoke')
# axs[1].grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
#
# # Set the x-axis ticks to show 1 day, 2 days, 3 days, etc.
# xticks = [i for i in range(0, times[-1], 24)]
# xticklabels = [i for i in range(0, int(times[-1]/24))]
# axs[0].set_xticks(xticks)
# axs[0].set_xticklabels(xticklabels)
# axs[1].set_xticks(xticks)
# axs[1].set_xticklabels(xticklabels)

# Add vertical black arrows to show times of irradiation
# irradiation_times_cells = [72, 96, 120, 144, 168]  # times of irradiation in hours
# for time in irradiation_times_cells:
#     t = time  # convert hours to days
#     id = np.where(np.array(range(0, 200, 5)) == t)[0][0]  # get index of time
#     arrow = axs[0].annotate('', xy=(t, number_cells[id] + 300), xytext=(t, number_cells[id] + 301), arrowprops=dict(facecolor='black', width=1.5, headwidth=5))
#     arrow.set_zorder(-1)  # set arrow below plot line
#
# # Add vertical arrows to the tumor size plot
# irradiation_times_size = [72, 96, 120, 144, 168]  # times of irradiation in hours
# for time in irradiation_times_size:
#     t = time  # convert hours to days
#     id = np.where(np.array(range(4, 204, 4)) == t)[0][0]  # get index of time
#     arrow = axs[1].annotate('', xy=(t, tumor_size[id] + 0.5), xytext=(t, tumor_size[id] + 0.6), arrowprops=dict(facecolor='black', width=1.5, headwidth=5))
#     arrow.set_zorder(-1)  # set arrow below plot line

# legend_elements = [Line2D([0], [0], marker='>', color='black', lw=0, label='Irradiation with 6MeV photon beam')]

# # Add the custom legend to both subplots
# axs[0].legend(handles=legend_elements, loc='upper left')
# axs[1].legend(handles=legend_elements, loc='upper left')
#
# plt.tight_layout()
# plt.savefig('tumorgrowth.png')
# plt.show()