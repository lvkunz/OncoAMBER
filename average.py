import random

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os
import pandas as pd
from PIL import Image
import re

tmin = 0  # Minimum time
tmax = 5000 # Maximum time
show_fits = 1  # Show the exponential fits
fit = 'gompertz' #gompertz or exp
show_necro = 1
show_quiet_cycling = 1
show_vessels = True
local = False
param_to_plot = []

def plot_outliners(ax, x, y, y_min, y_max, color='black'):
    for i in range(len(x)):
        if y[i] < y_min[i] or y[i] > y_max[i]:
            ax.plot(x[i], y[i], '.', color=color)

repo = '20230607_lk001_Linux/CONFIG_vasculature_example.py_152741'

csv_file = ''
#all repositories in repo:
#find csv file in repo
for filename in os.listdir(repo):
    # Check if the file is a csv file
    if filename.endswith('.csv'):
        csv_file = filename

param_space = pd.read_csv(f'{repo}/{csv_file}', sep=' ', header=0)
print(param_space)

print(param_space.columns)
print(param_space.columns[1])
parameter = param_space.columns[1]

param = np.array(param_space[parameter])
number_of_iterations = len(param_space['Iteration'])

iter = []

if iter == []:
    iter = [i for i in range(number_of_iterations)]

paths = [f'{repo}/iter{i}/DataOutput/' for i in iter]
#remove paths 4
print(paths)



if local: paths = ['DataOutput/']


number_cells_list = []
necrotic_cells_list = []
cycling_cells_list = []
quiescent_cells_list = []
tumor_size_list = []
tumor_size_free_list = []
number_vessels_list = []
times_list = []

for path in paths:
    number_cells = np.load(f'{path}number_tumor_cells.npy', allow_pickle=True)
    necrotic_cells = np.load(f'{path}number_necrotic_cells.npy', allow_pickle=True)
    cycling_cells = np.load(f'{path}number_cycling_cells.npy', allow_pickle=True)
    quiescent_cells = np.load(f'{path}number_quiescent_cells.npy', allow_pickle=True)
    tumor_size = np.load(f'{path}tumor_size.npy', allow_pickle=True)
    tumor_size_free = np.load(f'{path}tumor_size_free.npy', allow_pickle=True)
    tumor_size_free = tumor_size - tumor_size_free
    number_vessels = np.load(f'{path}number_vessels.npy', allow_pickle=True)
    times = np.load(f'{path}times.npy', allow_pickle=True)

    # Find the indices of the times that are within the time range
    idx = np.where((times >= tmin) & (times<=tmax))
    # Filter the arrays to only include the data between tmin and tmax
    number_cells = number_cells[idx]
    tumor_size = tumor_size[idx]
    tumor_size_free = tumor_size_free[idx]
    necrotic_cells = necrotic_cells[idx]
    cycling_cells = cycling_cells[idx]
    quiescent_cells = quiescent_cells[idx]
    number_vessels = number_vessels[idx]
    times = times[idx]

    # Append the filtered arrays to the lists
    number_cells_list.append(number_cells)
    tumor_size_list.append(tumor_size)
    tumor_size_free_list.append(tumor_size_free)
    necrotic_cells_list.append(necrotic_cells)
    cycling_cells_list.append(cycling_cells)
    quiescent_cells_list.append(quiescent_cells)
    number_vessels_list.append(number_vessels)
    times_list.append(times)


# Fit the data to an exponential curve for each simulation and get the doubling time
doubling_times_number_cells = []
doubling_times_tumor_size = []
dpi = 300
fig, axes = plt.subplots(2, 1, figsize=(8, 10), dpi=dpi)

#for each y point, find the average and SD

max_times = []
for array in times_list:
    if len(array) > len(max_times):
        max_times = array

number_cells_average = [] ; number_cells_sd = []
necrotic_cells_average = [] ; necrotic_cells_sd = []
cycling_cells_average = [] ; cycling_cells_sd = []
quiescent_cells_average = [] ; quiescent_cells_sd = []
tumor_size_average = [] ; tumor_size_sd = []
tumor_size_free_average = [] ; tumor_size_free_sd = []
number_vessels_average = [] ; number_vessels_sd = []

for t in max_times:
    idd = np.where(max_times == t)
    number_cells = []
    necrotic_cells = []
    cycling_cells = []
    quiescent_cells = []
    tumor_size = []
    tumor_size_free = []
    number_vessels = []
    for i in range(len(paths)):
        if idd[0][0] < len(number_cells_list[i]):
            number_cells.append(number_cells_list[i][idd])
            necrotic_cells.append(necrotic_cells_list[i][idd])
            cycling_cells.append(cycling_cells_list[i][idd])
            quiescent_cells.append(quiescent_cells_list[i][idd])
            tumor_size.append(tumor_size_list[i][idd])
            tumor_size_free.append(tumor_size_free_list[i][idd])
            number_vessels.append(number_vessels_list[i][idd])
    number_cells_average.append(np.mean(number_cells))
    necrotic_cells_average.append(np.mean(necrotic_cells))
    cycling_cells_average.append(np.mean(cycling_cells))
    quiescent_cells_average.append(np.mean(quiescent_cells))
    tumor_size_average.append(np.mean(tumor_size))
    tumor_size_free_average.append(np.mean(tumor_size_free))
    number_vessels_average.append(np.mean(number_vessels))
    number_cells_sd.append(np.std(number_cells))
    necrotic_cells_sd.append(np.std(necrotic_cells))
    cycling_cells_sd.append(np.std(cycling_cells))
    quiescent_cells_sd.append(np.std(quiescent_cells))
    tumor_size_sd.append(np.std(tumor_size))
    tumor_size_free_sd.append(np.std(tumor_size_free))
    number_vessels_sd.append(np.std(number_vessels))

times_average = max_times

# Fit number of cells
if fit == 'exp':
    def func_cell(x, a, b):
        return a * (np.exp(b * x)) + number_cells_average[0] - a

    def func_volume(x, a, b):
        return a * (np.exp(b * x)) + tumor_size_average[0] - a

    p1 = (3000, 3e-3)
    p2 = (1, 0.003)

elif fit == 'gompertz':

    def func_cell(x, a, b):
        return a * np.exp(np.log(number_cells_average[0]/a)*np.exp(-b * x))

    def func_volume(x, a, b):
        return a * np.exp(np.log(tumor_size_average[0]/a)*np.exp(-b * x))

    p1 = ( 3000000, 0.01)
    p2 = ( 3000 , 0.01)

elif fit == 'exp-linear':

    def func_cell(x, a, b, c, d, e):
        return a * np.exp(b * x) + c * np.log(d * x) + e


    def func_volume(x, a, b, c, d, e):
        return a * np.exp(b * x) + c * np.log(d * x) + e


if show_fits:
    popt, pcov = curve_fit(func_cell, max_times, number_cells_average, p0=p1, maxfev=100000)
axes[0].plot(times_average, number_cells_average, '-', label = 'Number of Cells', color='blue', alpha=0.8)
axes[0].fill_between(times_average, np.array(number_cells_average) - np.array(number_cells_sd), np.array(number_cells_average) + np.array(number_cells_sd), alpha=0.2, color='blue')
# for i in range(len(paths)):
#     plot_outliners(axes[0], times_list[i], number_cells_list[i], np.array(number_cells_average) - np.array(number_cells_sd), np.array(number_cells_average) + np.array(number_cells_sd), color='black')

if show_necro:
    axes[0].plot(times_average, necrotic_cells_average, '-', label = 'Necrotic Cells', color='black', alpha=0.8)
    axes[0].fill_between(times_average, np.array(necrotic_cells_average) - np.array(necrotic_cells_sd), np.array(necrotic_cells_average) + np.array(necrotic_cells_sd), alpha=0.2, color='black')
if show_quiet_cycling:
    axes[0].plot(times_average , cycling_cells_average , '-', label = 'Cycling Cells', color='red', alpha=0.8)
    axes[0].fill_between(times_average, np.array(cycling_cells_average) - np.array(cycling_cells_sd), np.array(cycling_cells_average) + np.array(cycling_cells_sd), alpha=0.2, color='red')
    axes[0].plot(times_average , quiescent_cells_average , '-', label = 'Quiescent Cells', color='green', alpha=0.8)
    axes[0].fill_between(times_average, np.array(quiescent_cells_average) - np.array(quiescent_cells_sd), np.array(quiescent_cells_average) + np.array(quiescent_cells_sd), alpha=0.2, color='green')
if show_fits:
    axes[0].plot(times_average , func_cell(times_average, *popt), '--', color='blue')
    #diff fitting
    axes[0].plot(times_average , number_cells_average- func_cell(times_average, *popt), '--', color='orange')


    doubling_time = np.log(2) / popt[1]
    if fit == 'gompertz':
        print('Doubling time (Number of Cells):', doubling_time)
        print('Max Carrying capacity:', popt[0])
    elif fit == 'exp':
        print('Doubling time (Number of Cells):', doubling_time)

    print(popt)
    doubling_times_number_cells.append(doubling_time)

# Fit tumor size
if show_fits:
    popt, pcov = curve_fit(func_volume, times_average , tumor_size_average , p0=p2, maxfev=100000)
axes[1].plot(times_average, tumor_size_average, '-', color = 'purple', markersize = 5, alpha=0.5, label='Tumor Volume')
axes[1].fill_between(times_average, np.array(tumor_size_average) - np.array(tumor_size_sd), np.array(tumor_size_average) + np.array(tumor_size_sd), alpha=0.2, color='purple')
axes[1].plot(times_average, tumor_size_free_average, '-', color = 'black', markersize = 5, alpha=0.5, label='Necrotic Core Volume')
axes[1].fill_between(times_average, np.array(tumor_size_free_average) - np.array(tumor_size_free_sd), np.array(tumor_size_free_average) + np.array(tumor_size_free_sd), alpha=0.2, color='black')
if show_fits:
    axes[1].plot(times_list[i], func_volume(times_average, *popt), '--', color='purple')#, label='fit '+parameter+': '+str(param[i]))
    doubling_time = np.log(2)/popt[1]
    doubling_times_tumor_size.append(doubling_time)

axes[0].set_title('Number of Cells Evolution')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Number of Cells')
# axes[0].set_xlim(0, 250)
# axes[0].set_ylim(0, 5e5)
axes[0].grid(True)
axes[0].legend()

axes[1].set_title('Tumor Volume Evolution')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Tumor Volume [mm^3]')
# axes[1].set_xlim(0, 250)
# axes[1].set_ylim(0, 50)
axes[1].grid(True)
axes[1].legend()

#add a tiny text box in the corner with the repo name
plt.figtext(0.01, 0.01, repo, wrap=True, horizontalalignment='left', fontsize=6)
plt.tight_layout()
plt.savefig(repo+'/tumor_evolution_average'+str(tmax)+'.png', dpi=dpi)
plt.show()

if show_vessels:
    fig, axes = plt.subplots(1, 1, figsize=(8, 5), dpi=dpi)
    axes.plot(times_average, number_vessels_average , '-', color = 'crimson', markersize = 5, alpha=0.5, label='Number of Vessels')
    axes.fill_between(times_average, np.array(number_vessels_average) - np.array(number_vessels_sd), np.array(number_vessels_average) + np.array(number_vessels_sd), alpha=0.2, color='purple')
    axes.set_title('Number of Vessels Evolution')
    axes.set_xlabel('Time')
    axes.set_ylabel('Number of Vessels')
    axes.grid(True)
    axes.legend()
    plt.tight_layout()
    plt.savefig(repo+'/vessels_evolution_average'+str(tmax)+'.png', dpi=dpi)
    plt.show()


rates_average = []
dt = times_average[1] - times_average[0]
for t in times_average:
    idd = np.where(times_average == t)[0][0]
    if t <= 10*dt:
        rates_average.append(0)
    else:
        rate = (number_cells_average[idd] - number_cells_average[idd-10])/(10*dt)
        rates_average.append(rate)


fig, axes = plt.subplots(2, 1, figsize=(8, 5), dpi=dpi)
for i in range(len(times_average)):
    axes[0].plot(number_cells_average[i], rates_average[i], 'o', color='black', markersize=5, alpha=0.5)
axes[0].set_title('Growth Rate vs Number of Cells')
axes[0].set_xlabel('Number of Cells')
axes[0].set_ylabel('Growth Rate')
axes[0].grid(True)


axes[1].plot(times_average, rates_average, '-', color = 'crimson', markersize = 5, alpha=0.5, label='Growth Rate')
axes[1].set_title('Growth Rate Evolution')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Growth Rate')
axes[1].grid(True)
axes[1].legend()
plt.tight_layout()
plt.savefig(repo+'/growth_rate'+str(tmax)+'.png', dpi=dpi)
plt.show()

