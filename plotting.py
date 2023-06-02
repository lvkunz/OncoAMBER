import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os
import pandas as pd
from PIL import Image
import re



def create_gif(image_dir, output_path, image_suffix, image_step=1):
    images = []
    t = 0
    while True:
        image_file = os.path.join(image_dir, f"t{t}_{image_suffix}.png")
        if not os.path.exists(image_file):
            break

        print(f"Adding image {image_file}")
        img = Image.open(image_file)
        images.append(img)

        t += image_step

    if images:
        # Save the first image as the GIF background
        images[0].save(output_path, save_all=True, append_images=images[1:], loop=0, duration=image_step*10)
        print(f"GIF created successfully at {output_path}")
    else:
        print("No images found to create the GIF.")

tmin = 0  # Minimum time
tmax = 10000  # Maximum time
show_fits = False  # Show the exponential fits
fit = 'gompertz' #gompertz or exp
show_necro = True
show_quiet_cycling = False
show_vessels = True
local = False
gif = True
param_to_plot = []

repo = '20230601_lk001_Linux/CONFIG_longer_example.py_144101'

iter = [0,1,2]

for i in iter:
    print('iter', i, 'of', iter, 'is being processed')
    image_directory = repo+'/iter'+str(i)+'/Plots/CurrentPlotting'
    image_sufix1 = 'AllPlots'
    output_path1 = repo + '/'+image_sufix1+str(i)+'.gif'
    image_sufix2 = 'Vasculature'
    output_path2 = repo + '/'+image_sufix2+str(i)+'.gif'
    image_step = 20

    if gif:
        create_gif(image_directory, output_path1, image_sufix1, image_step)
        create_gif(image_directory, output_path2, image_sufix2, image_step)

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

paths = [f'{repo}/iter{i}/DataOutput/' for i in range(0, number_of_iterations)]
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
        axes[0].plot(times_list[i], func_cell(times_list[i], *popt), '-', color=color, label='fit '+parameter+': '+str(param[i]))



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
    axes[1].plot(times_list[i], tumor_size_free_list[i], '+', color = color, markersize = 5, alpha=0.5)
    if show_fits:
        axes[1].plot(times_list[i], func_volume(times_list[i], *popt), '-', color=color, label='fit '+parameter+': '+str(param[i]))
        doubling_time = np.log(2)/popt[1]
        doubling_times_tumor_size.append(doubling_time)

axes[0].set_title('Number of Cells Evolution')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Number of Cells')
# axes[0].set_xlim(0, 250)
# axes[0].set_ylim(0, 2e5)
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
plt.savefig(repo+'/tumor_evolution_'+str(tmax)+'.png', dpi=dpi)
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
    plt.savefig(repo+'/vessels_evolution_'+str(tmax)+'.png', dpi=dpi)
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
    plt.savefig(repo+'/doubling_time.png', dpi=300)
    plt.show()

    plt.plot(param, doubling_times_tumor_size, 'ro', label='Tumor volume doubling time')
    plt.xlabel(parameter)
    plt.ylabel('Doubling time [days]')
    plt.title('Doubling time vs. ' + parameter)
    # plt.yscale('log')  # set y-axis to logarithmic scale
    plt.legend()
    plt.grid(True)
    plt.savefig(repo+'/doubling_time_tumor_size.png', dpi=300)
    plt.show()






# Provide the directory containing the images and the output path for the GIF

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



