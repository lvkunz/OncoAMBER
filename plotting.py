import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os
import pandas as pd
from PIL import Image
import re

tmin = 0  # Minimum time
tmax = 4003 # Maximum time
show_fits = 0  # Show the exponential fits
fit = 'exp' #gompertz or exp
show_necro = 0
show_quiet_cycling = 0
show_vessels = False
show_rates = True
experimental = 0
rate_choice = 'volume' #volume or number
local = 0
param_to_plot = []

repo = 'CONFIG_growth_example.py_160329_paper'
repo = repo + '/'

csv_file = ''
#all repositories in repo:
#find csv file in repo
for filename in os.listdir(repo):
    # Check if the file is a csv file
    if filename.endswith('.csv'):
        csv_file = filename
print(csv_file)
param_space = pd.read_csv(f'{repo}/{csv_file}', sep=' ', header=0)
print(param_space)

print(param_space.columns)
print(param_space.columns[1])
parameter = param_space.columns[1]

param = np.array(param_space[parameter])
iterations = np.array(param_space['Iteration'])

paths = [f'{repo}iter{i}/DataOutput/' for i in iterations]
# paths = [f'{repo}iter{i}/DataOutput/' for i in [0,1,2,3]]
# param = [param[i] for i in [0,1,2,3]]
#remove paths 4



if local:
    paths = ['DataOutput/']
    parameter = ''
    param = [1]
    repo = ''

number_cells_list = []
necrotic_cells_list = []
cycling_cells_list = []
quiescent_cells_list = []
tumor_size_list = []
tumor_size_free_list = []
number_vessels_list = []
rates_list = []
times_list = []

for path in paths:
    # number_cells = np.load(f'{path}number_tumor_cells.npy',allow_pickle=True)
    necrotic_cells = np.load(f'{path}number_necrotic_cells.npy',allow_pickle=True)
    cycling_cells = np.load(f'{path}number_cycling_cells.npy', allow_pickle=True)
    quiescent_cells = np.load(f'{path}number_quiescent_cells.npy', allow_pickle=True)
    number_cells = cycling_cells + quiescent_cells + necrotic_cells
    tumor_size = np.load(f'{path}tumor_size.npy', allow_pickle=True)
    tumor_size_free = np.load(f'{path}tumor_size_free.npy', allow_pickle=True)
    if show_vessels:
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
    if show_vessels:
        number_vessels = number_vessels[idx]
    times = times[idx]

    rates = []
    dt = times[1] - times[0]
    for t in times:
        idd = np.where(times == t)[0][0]
        if t <= 10 * dt:
            rates.append(0)
        else:
            if rate_choice == 'volume':
                rate = (tumor_size[idd] - tumor_size[idd - 10]) / (10 * dt)
            elif rate_choice == 'number':
                rate = (number_cells[idd] - number_cells[idd - 10]) / (10 * dt)
            rates.append(rate)


    # Append the filtered arrays to the lists
    number_cells_list.append(number_cells)
    tumor_size_list.append(tumor_size)
    tumor_size_free_list.append(tumor_size_free)
    necrotic_cells_list.append(necrotic_cells)
    cycling_cells_list.append(cycling_cells)
    quiescent_cells_list.append(quiescent_cells)
    if show_vessels:
        number_vessels_list.append(number_vessels)
    rates_list.append(rates)
    times_list.append(times)


# Fit the data to an exponential curve for each simulation and get the doubling time
doubling_times_number_cells = []
doubling_times_tumor_size = []
dpi = 300
#change font size
fig, axes = plt.subplots(2, 1, figsize=(8, 10), dpi=dpi)
#change font size

for i in range(len(paths)):
    print(param[i])
    if len(param_to_plot) > 0:
        if param[i] not in param_to_plot:
            continue
    print(paths[i])
    # Fit number of cells
    if fit == 'exp':
        def func_cell(x, a, b):
            return a * (np.exp(b * x)) + number_cells_list[i][0] - a

        def func_volume(x, a, b):
            return a * (np.exp(b * x)) + tumor_size_list[i][0] - a

        p1 = (3000, 3e-3)
        p2 = (1, 0.003)

    elif fit == 'gompertz':

        def func_cell(x, a, b):
            return a * np.exp(np.log(number_cells_list[i][0]/a)*np.exp(-b * x))

        def func_volume(x, a, b):
            return a * np.exp(np.log(tumor_size_list[i][0]/a)*np.exp(-b * x))

        p1 = ( 30000, 0.1)
        p2 = ( 300 , 0.1)

    if show_fits:
        popt, pcov = curve_fit(func_cell, times_list[i], number_cells_list[i], p0=p1, maxfev=100000)
    color = axes[0].plot(times_list[i], number_cells_list[i], '.', markersize=3, alpha=0.8, label=parameter+': '+str(param[i]))[0].get_color()
    if show_necro: axes[0].plot(times_list[i], necrotic_cells_list[i], 's', markersize=5, alpha=0.5, color=color)
    if show_quiet_cycling:
        axes[0].plot(times_list[i], cycling_cells_list[i], '+', markersize=3, alpha=0.5, color=color)
        axes[0].plot(times_list[i], quiescent_cells_list[i], 'D', markersize=3, alpha=0.5, color=color)
    if show_fits:
        axes[0].plot(times_list[i], func_cell(times_list[i], *popt), '-', color=color)#, label='fit '+parameter+': '+str(param[i]))



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
        popt, pcov = curve_fit(func_volume, times_list[i], tumor_size_list[i], p0=p2, maxfev=100000)
    axes[1].plot(times_list[i], tumor_size_list[i], 'o', color = color, markersize = 5, alpha=0.5, label=parameter+': '+str(param[i]))
    # axes[1].plot(times_list[i], tumor_size_list[i], 'o', color = 'black', markersize = 5, alpha=1, label='Model Values')
    axes[1].plot(times_list[i], tumor_size_free_list[i], '+', color = color, markersize = 5, alpha=0.5)
    if show_fits:
        axes[1].plot(times_list[i], func_volume(times_list[i], *popt), '-', color=color)#, label='fit '+parameter+': '+str(param[i]))
        doubling_time = np.log(2)/popt[1]
        doubling_times_tumor_size.append(doubling_time)

data = pd.read_csv('data_exp.csv', sep=',', header=0)
data2 = pd.read_csv('LLC_sc_CCSB.csv', sep=',', header=0)

data2['Time'] = data2['Time'] * 24 + 72

data['time'] = data[data['time'] <= tmax]['time'] - 95
data['ctrl'] = data[data['time'] <= tmax]['ctrl']
data['RT'] = data[data['time'] <= tmax]['RT']

# Scatter plot of 'time' vs 'ctrl'
if experimental:
    time = []
    volume = []
    sd = []
    for i in data2['Time'].unique():
        time.append(i)
        vol = []
        for j in data2[data2['Time'] == i]['Vol']:
            vol.append(j)
        volume.append(np.mean(vol))
        sd.append(np.std(vol))

    axes[1].plot(data['time'], data['ctrl'], 'x', color='red', label='Experimental Data 1 (Zou 2020)', markersize=10)
    # axes[1].plot(data['time'], data['RT'], 'x', color='blue', markersize = 10, linewidth=2)
    axes[1].errorbar(time, volume, yerr=sd, fmt='+', color='blue', label='Experimental Data 2 (Benzekry 2014)', markersize=10, linewidth=1)


axes[0].set_title('Number of Cells Evolution')
axes[0].set_xlabel('Time [h]', fontsize=16)
axes[0].set_ylabel('Number of Cells', fontsize=16)
# axes[0].set_xlim(0, 250)
# axes[0].set_ylim(0, 5e5)
axes[0].grid(True)
axes[0].legend(loc = 'best', fontsize = 12)
axes[1].set_xlabel('Time [h]', fontsize=16)
axes[1].set_ylabel('Tumor Volume [mm^3]', fontsize=16)
# axes[1].set_xlim(0, 250)
# axes[1].set_ylim(0, 50)
#change the x axis font size
# for tick in axes[1].xaxis.get_major_ticks():
#     tick.label.set_fontsize(14)
# #change the y axis font size
# for tick in axes[1].yaxis.get_major_ticks():
#     tick.label.set_fontsize(14)

axes[1].grid(True)
axes[1].legend(loc = 'best', fontsize = 12)

#add a tiny text box in the corner with the repo name
plt.figtext(0.01, 0.01, repo, wrap=True, horizontalalignment='left', fontsize=6)
plt.tight_layout()
plt.savefig(repo+'tumor_evolution_'+str(tmax)+'.png', dpi=dpi)
plt.show()

if show_vessels:
    fig, axes = plt.subplots(1, 1, figsize=(8, 5), dpi=dpi)
    for i in range(len(paths)):
        if len(param_to_plot) > 0:
            if param[i] not in param_to_plot:
                continue
        axes.plot(times_list[i], number_vessels_list[i], 'o', markersize=5, alpha=0.5, label=parameter+': '+str(param[i]))
    axes.set_title('Number of Vessels Evolution')
    axes.set_xlabel('Time')
    axes.set_ylabel('Number of Vessels')
    axes.grid(True)
    axes.legend()
    plt.tight_layout()
    plt.savefig(repo+'vessels_evolution_'+str(tmax)+'.png', dpi=dpi)
    plt.show()

if show_fits:
    if len(param_to_plot) > 0:
        param = param_to_plot
    print('Doubling times (Number of Cells):', doubling_times_number_cells)
    print('Doubling times (Tumor Size):', doubling_times_tumor_size)

    plt.plot(param, doubling_times_number_cells, 'bo', label='Cells doubling time')
    plt.xlabel(parameter)
    plt.ylabel('Doubling time [days]')
    plt.title('Doubling time vs. ' + parameter)
    # plt.yscale('log')  # set y-axis to logarithmic scale
    plt.legend()
    plt.grid(True)
    plt.savefig(repo+'doubling_time.png', dpi=300)
    plt.show()


    plt.plot(param, doubling_times_tumor_size, 'ro', label='Tumor volume doubling time')
    plt.xlabel(parameter)
    plt.ylabel('Doubling time [days]')
    plt.title('Doubling time vs. ' + parameter)
    # plt.yscale('log')  # set y-axis to logarithmic scale
    plt.legend()
    plt.grid(True)
    plt.savefig(repo+'doubling_time_tumor_size.png', dpi=300)
    plt.show()

if show_rates:
    fig, axes = plt.subplots(2, 2, figsize=(12, 12), dpi=dpi)
    for i in range(len(paths)):
        if len(param_to_plot) > 0:
            if param[i] not in param_to_plot:
                continue

        axes[0, 0].plot(number_cells_list[i], rates_list[i], 'o', markersize=1, alpha=0.5)
        axes[1,0].plot(times_list[i], rates_list[i], 'o', markersize=1, alpha=0.5, label=parameter+': '+str(param[i]))
        axes[0,1].plot(number_vessels_list[i], rates_list[i], 'o', markersize=1, alpha=0.5, label=parameter+': '+str(param[i]))
        axes[1,1].plot(cycling_cells_list[i], rates_list[i], 'o', markersize=1, alpha=0.5, label=parameter+': '+str(param[i]))

    axes[0,0].set_title('Growth Rate vs Number of Cells')
    axes[0,0].set_xlabel('Number of Cells')
    axes[0,0].set_ylabel('Growth Rate')
    axes[0,0].set_ylim(0, None)
    axes[0,0].grid(True)
    axes[1,0].set_title('Growth Rate Evolution')
    axes[1,0].set_xlabel('Time')
    axes[1,0].set_ylabel('Growth Rate')
    axes[1, 0].set_ylim(0, None)
    axes[1,0].grid(True)
    axes[1,0].legend()
    axes[0,1].set_title('Growth Rate vs Number of Vessels')
    axes[0,1].set_xlabel('Number of Vessels')
    axes[0,1].set_ylabel('Growth Rate')
    axes[0,1].set_ylim(0, None)
    axes[0,1].grid(True)
    axes[0,1].legend()
    axes[1,1].set_title('Growth Rate vs Number of Cycling Cells')
    axes[1,1].set_xlabel('Number of Cycling Cells')
    axes[1,1].set_ylabel('Growth Rate')
    axes[1,1].set_ylim(0, None)
    axes[1,1].grid(True)
    axes[1,1].legend()
    plt.tight_layout()
    plt.savefig(repo + 'growth_rate' + str(tmax) + '.png', dpi=dpi)
    plt.show()
