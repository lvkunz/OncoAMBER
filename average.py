import random

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import quad
from scipy.optimize import curve_fit
import glob
import os
import pandas as pd
from PIL import Image
import re
from scipy.integrate import odeint


tmin = 0  # Minimum time
tmax = 2700 # Maximum time
show_fits = 0  # Show the exponential fits
fit = 'gompertz' #gompertz or exp
show_necro = 1
show_quiet_cycling = 1
show_vessels = False
local = 0
irradiation = [796, 24, 5] #first, frequence, number of fractions

def plot_outliners(ax, x, y, y_min, y_max, color='black'):
    for i in range(len(x)):
        if y[i] < y_min[i] or y[i] > y_max[i]:
            ax.plot(x[i], y[i], '.', color=color)

repo = '20230814_lk001_Linux/CONFIG_vasculature_irrad_example.py_134514'

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

# iter = [0,1,2,3,4]
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
    # number_cells = np.load(f'{path}number_tumor_cells.npy', allow_pickle=True)
    necrotic_cells = np.load(f'{path}number_necrotic_cells.npy', allow_pickle=True)
    cycling_cells = np.load(f'{path}number_cycling_cells.npy', allow_pickle=True)
    quiescent_cells = np.load(f'{path}number_quiescent_cells.npy', allow_pickle=True)
    number_cells = cycling_cells + quiescent_cells + necrotic_cells
    tumor_size = np.load(f'{path}tumor_size.npy', allow_pickle=True)
    tumor_size_free = np.load(f'{path}tumor_size_free.npy', allow_pickle=True)
    # tumor_size_free = tumor_size - tumor_size_free
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
fig, axes = plt.subplots(2, 1, figsize=(16, 10), dpi=dpi)

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
axes[1].plot(times_average, tumor_size_average, '-', color = 'black', markersize = 5, alpha=0.8, label='Total Tumor Volume')
axes[1].fill_between(times_average, np.array(tumor_size_average) - np.array(tumor_size_sd), np.array(tumor_size_average) + np.array(tumor_size_sd), alpha=0.2, color='black')
axes[1].plot(times_average, tumor_size_free_average, '-', color = 'purple', markersize = 5, alpha=0.8, label='Non-Necrotic Tumor Volume')
axes[1].fill_between(times_average, np.array(tumor_size_free_average) - np.array(tumor_size_free_sd), np.array(tumor_size_free_average) + np.array(tumor_size_free_sd), alpha=0.2, color='purple')
if show_fits:
    axes[1].plot(times_list[i], func_volume(times_average, *popt), '--', color='black')#, label='fit '+parameter+': '+str(param[i]))
    doubling_time = np.log(2)/popt[1]
    doubling_times_tumor_size.append(doubling_time)



# # Set the x-axis ticks to show 1 day, 2 days, 3 days, etc.
xticks = [i for i in range(0, times_average[-1], 168)]
print(xticks)
xticklabels = [i for i in range(0, int(times_average[-1]/24)+1, 7)]
print(xticklabels)
axes[0].set_xticks(xticks)
axes[0].set_xticklabels(xticklabels)
axes[1].set_xticks(xticks)
axes[1].set_xticklabels(xticklabels)
#
# # Add vertical black arrows to show times of irradiation
irradiation_times_cells = [irradiation[0] + i*irradiation[1] for i in range(0,irradiation[2])]  # times of irradiation in hours
shift = 10000
shift2 = 3
for time in irradiation_times_cells:
    id = np.where(times_average == time)[0][0]  # get index of time
    arrow1 = axes[0].annotate('', xy=(times_average[id], number_cells_average[id] + shift), xytext=(times_average[id], number_cells_average[id] + shift+1), arrowprops=dict(facecolor='black', width=1.5, headwidth=5))
    arrow1.set_zorder(-1)  # set arrow below plot line
    arrow2 = axes[1].annotate('', xy=(times_average[id], tumor_size_average[id] + shift2), xytext=(times_average[id], tumor_size_average[id] + shift2+1), arrowprops=dict(facecolor='black', width=1.5, headwidth=5))
    arrow2.set_zorder(-1)  # set arrow below plot line

legend_elements = [Line2D([0], [0], marker='>', color='black', lw=0, label='Irradiation with 6MeV photon beam')]
#add a tiny text box in the corner with the repo name
plt.figtext(0.01, 0.01, repo, wrap=True, horizontalalignment='left', fontsize=6)

axes[0].set_title('Number of Cells Evolution', fontsize=14)
axes[0].set_xlabel('Time [Day]', fontsize=16)
axes[0].set_ylabel('Number of Cells', fontsize=16)
# axes[0].set_xlim(0, 250)
axes[0].set_ylim(0, None)
axes[0].grid(True)
axes[0].legend(fontsize=16, loc='upper left')

axes[1].set_title('Tumor Volume Evolution', fontsize=14)
axes[1].set_xlabel('Time [Day]', fontsize=16)
axes[1].set_ylabel('Tumor Volume [mm^3]', fontsize=16)
# axes[1].set_xlim(0, 250)
axes[1].set_ylim(0, None)
axes[1].grid(True)
axes[1].legend(fontsize = 16, loc = 'upper left')

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
    if t <= 3*dt:
        rates_average.append(0)
    else:
        rate = (number_cells_average[idd] - number_cells_average[idd-3])/(3*dt)
        rates_average.append(rate)

def func_gompertz(x, a, b):
    x = np.array(x)
    return x * (b - a * np.log(x))

def func_gompertz_(x,t, a, b):
    x = np.array(x)
    return x * (b - a * np.log(x))
def func_logistic(x, a, b):
    x = np.array(x)
    return a * x - b * x * x

def func_logistic_(x,t, a, b):
    x = np.array(x)
    return a * x - b * x * x

def func_holling(x, a, b, k):
    x = np.array(x)
    return (a * x)/(k + x) - b * x

def func_holling_(x, t, a, b, k):
    return (a * x)/(k + x) - b * x

param0 = np.array([1.0, 1.0, 1.0])
param0_log = np.array([1.0, 100.0])
param0_ = np.array([1.0, 1.0])
#add a random component to the initial guess
# param0 = np.array(param0) * (1 + 1.0 * np.random.random(len(param0)) - 0.1)
weights = np.zeros(len(number_cells_average))
for i in range(len(number_cells_average)):
    if number_cells_average[i] < 400000:
        weights[i] = 0.2
    else:
        weights[i] = 1

popt_g, pcov_g = curve_fit(func_gompertz, number_cells_average, rates_average, p0=param0_, maxfev=100000)
popt_l, pcov_l = curve_fit(func_logistic, number_cells_average, rates_average, p0=param0_log, maxfev=100000, sigma=weights)
popt_h, pcov_h = curve_fit(func_holling, number_cells_average, rates_average, p0=param0, maxfev=100000)
print('Fitted parameters Gompertz: ', popt_g)
print('Fitted parameters Logistic: ', popt_l)
print('Fitted parameters Holling: ', popt_h)
#N fin

N_final = popt_h[0]/popt_h[1] - popt_h[2]
print('N final: ', N_final)
#solve the ODE for the Holling model
# find at what time the tumor reaches N_final
times_final = []
for i in range(len(number_cells_average)):
    if number_cells_average[i] > N_final:
        times_final.append(times_average[i])
        break
fig, axes = plt.subplots(2, 2, figsize=(12, 12), dpi=dpi)
# popt_r = np.array([1.0, 1.0])
axes[0, 0].plot(number_cells_average, rates_average, 'o', markersize=3, color='orange', label='Model Values')
axes[0,0 ].plot(number_cells_average, func_gompertz(number_cells_average, *popt_g), linestyle = 'dotted', color='black', label='Gompertz fit')
axes[0, 0].plot(number_cells_average, func_logistic(number_cells_average, *popt_l), linestyle = 'dashed', color='black', label='Logistic fit')
axes[0, 0].plot(number_cells_average, func_holling(number_cells_average, *popt_h), linestyle = 'dashdot', color='black', label='Holling fit')
axes[0, 0].legend(fontsize = 14)
axes[1, 0].plot(times_average, rates_average, 'o', markersize=3, color='black')
axes[0, 1].plot(number_vessels_average, rates_average, 'o', markersize=3, color='crimson')
axes[1, 1].plot(cycling_cells_average, rates_average, 'o', markersize=3, color='green')

axes[0, 0].set_title('Growth Rate vs Number of Cells and fittings')
axes[0, 0].set_xlabel('Number of Cells',fontsize=14)
axes[0, 0].set_ylabel('Growth Rate',fontsize=14)
axes[0, 0].set_ylim(0, None)
axes[0, 0].grid(True)
axes[1, 0].set_title('Growth Rate Evolution')
axes[1, 0].set_xlabel('Time')
axes[1, 0].set_ylabel('Growth Rate')
axes[1, 0].set_ylim(0, None)
axes[1, 0].grid(True)
axes[0, 1].set_title('Growth Rate vs Number of Vessels')
axes[0, 1].set_xlabel('Number of Vessels')
axes[0, 1].set_ylabel('Growth Rate')
axes[0, 1].set_ylim(0, None)
axes[0, 1].grid(True)
axes[1, 1].set_title('Growth Rate vs Number of Cycling Cells')
axes[1, 1].set_xlabel('Number of Cycling Cells')
axes[1, 1].set_ylabel('Growth Rate')
axes[1, 1].set_ylim(0, None)
axes[1, 1].grid(True)
plt.tight_layout()
plt.savefig(repo + '/growth_rate' + str(tmax) + '.png', dpi=dpi)
plt.show()


#plot the integral of the growth rate to get the number of cells


# Define parameters
aH = popt_h[0]
bH = popt_h[1]
kH = popt_h[2]

aG = popt_g[0]
bG = popt_g[1]

aL = popt_l[0]
bL = popt_l[1]

# Define time points
t = times_average
# t = np.linspace(0, tmax*3, 6000)
# Set initial condition
x0 = number_cells_average[0]
# Solve the differential equation
x = odeint(func_holling_, x0, t, args=(aH, bH, kH)).flatten()
x2 = odeint(func_gompertz_, x0, t, args=(aG, bG)).flatten()
x3 = odeint(func_logistic_, x0, t, args=(aL, bL)).flatten()

fig, axes = plt.subplots(1, 1, figsize=(6, 6), dpi=dpi)
axes.plot(times_average, number_cells_average, '-', label = 'Model values', color='green', alpha=0.8)
axes.plot(t, x, linestyle = 'dashdot', label = 'Holling fit', color='black', alpha=1.0, linewidth=2)
axes.plot(t, x2, linestyle = 'dotted', label = 'Gompertz fit', color='black', alpha=1.0, linewidth=2)
axes.plot(t, x3, linestyle =  'dashed', label = 'Logistic fit', color='black', alpha=1.0, linewidth=2)
axes.fill_between(times_average, np.array(number_cells_average) - np.array(number_cells_sd), np.array(number_cells_average) + np.array(number_cells_sd), alpha=0.2, color='green')
axes.set_xlabel('Time [h]', fontsize=16)
axes.set_ylabel('Number of Cells', fontsize=16)
axes.set_title('Number of cells evolution and fittings', fontsize=16)
axes.legend(fontsize = 14)
axes.grid(True)
plt.tight_layout()
plt.savefig(repo + '/number_cells_fitted' + str(tmax) + '.png', dpi=dpi)
plt.show()




