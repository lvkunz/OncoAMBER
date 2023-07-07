import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os
import pandas as pd
from PIL import Image
import amber



def create_gif(image_dir, output_path, image_suffix, image_step=1):
    images = []
    t = 0
    while True:
        image_file = os.path.join(image_dir, f"t{t}{image_suffix}.png")
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
tmax = 5000 # Maximum time
show_necro = True
show_quiet_cycling = True
show_vessels = True
local = False
generate_images = False
generate_gif = True
repo = '20230705_lk001_Linux/CONFIG_vasculature_irrad_single_example.py_162736/iter3'

image_step = 8
irradiation = [796, 24, 1]


path = f'{repo}/DataOutput/'

# number_cells = np.load(f'{path}number_tumor_cells.npy', allow_pickle=True)
necrotic_cells = np.load(f'{path}number_necrotic_cells.npy', allow_pickle=True)
cycling_cells = np.load(f'{path}number_cycling_cells.npy', allow_pickle=True)
quiescent_cells = np.load(f'{path}number_quiescent_cells.npy', allow_pickle=True)
number_cells = necrotic_cells + cycling_cells + quiescent_cells
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

#create a plot for each time step and save it as a png
if generate_images:
    for t in times:
        print(t)
        idx = np.where((times <= t))
        # Filter the arrays to only include the data between tmin and tmax
        number_cells_cut = number_cells[idx]
        tumor_size_cut = tumor_size[idx]
        tumor_size_free_cut = tumor_size_free[idx]
        necrotic_cells_cut = necrotic_cells[idx]
        cycling_cells_cut = cycling_cells[idx]
        quiescent_cells_cut = quiescent_cells[idx]
        number_vessels_cut = number_vessels[idx]
        times_cut = times[idx]

        dpi = 100
        fig, axes = plt.subplots(2, 1, figsize=(8, 10), dpi=dpi)
        axes[0].plot(times, number_cells, '.', markersize=1, alpha=0.8,color='blue')
        axes[0].plot(times_cut, number_cells_cut, '-', markersize=3, alpha=0.8, color='blue')
        axes[0].scatter(times_cut[-1], number_cells_cut[-1], s=50, alpha=0.8, color='blue', label='Total number cells')
        if show_necro:
            axes[0].plot(times, necrotic_cells, 's', markersize=1, alpha=0.8, color='black')
            axes[0].plot(times_cut, necrotic_cells_cut, '-', markersize=3, alpha=0.8, color='black')
            axes[0].scatter(times_cut[-1], necrotic_cells_cut[-1], s=50, alpha=0.8, color='black', label='Necrotic cells')
        if show_quiet_cycling:
            axes[0].plot(times, cycling_cells, '+', markersize=1, alpha=0.5, color='green')
            axes[0].plot(times_cut, cycling_cells_cut, '-', markersize=3, alpha=0.5, color='green')
            axes[0].scatter(times_cut[-1], cycling_cells_cut[-1], s=50, alpha=0.5, color='green', label='Cycling cells')
            axes[0].plot(times, quiescent_cells, 'D', markersize=1, alpha=0.5, color='red')
            axes[0].plot(times_cut, quiescent_cells_cut, '-', markersize=3, alpha=0.5, color='red')
            axes[0].scatter(times_cut[-1], quiescent_cells_cut[-1], s=50, alpha=0.5, color='red', label='Quiescent cells')

        axes[1].plot(times, tumor_size, 'o', color = 'blue', markersize = 1, alpha=0.5)
        axes[1].plot(times_cut, tumor_size_cut, '-', color = 'blue', markersize = 3, alpha=0.5)
        axes[1].scatter(times_cut[-1], tumor_size_cut[-1], s=50, alpha=0.5, color='blue', label='Tumor Volume')

        # axes[1].plot(times, tumor_size_free, '+', color = 'red', markersize = 1, alpha=0.5)
        # axes[1].plot(times_cut, tumor_size_free_cut, '-', color = 'red', markersize = 3, alpha=0.5)
        # axes[1].scatter(times_cut[-1], tumor_size_free_cut[-1], s=50, alpha=0.5, color='red')

        axes[0].set_title('Number of Cells Evolution', fontsize = 14)
        axes[0].set_xlabel('Time [h]', fontsize = 16)
        axes[0].set_ylabel('Number of Cells', fontsize = 16)
        # axes[0].set_xlim(0, 250)
        # axes[0].set_ylim(0, 2e5)
        axes[0].grid(True)
        axes[0].legend(fontsize = 14)

        axes[1].set_title('Tumor Volume Evolution', fontsize = 14)
        axes[1].set_xlabel('Time [h]', fontsize = 16)
        axes[1].set_ylabel('Tumor Volume [mm^3]')
        # axes[1].set_xlim(0, 250)
        # axes[1].set_ylim(0, 50)
        axes[1].grid(True)
        axes[1].legend(fontsize = 14)

        # # Set the x-axis ticks to show 1 day, 2 days, 3 days, etc.
        xticks = [i for i in range(0, times[-1], 168)]
        print(xticks)
        xticklabels = [i for i in range(0, int(times[-1] / 24) + 1, 7)]
        print(xticklabels)
        axes[0].set_xticks(xticks)
        axes[0].set_xticklabels(xticklabels)
        axes[1].set_xticks(xticks)
        axes[1].set_xticklabels(xticklabels)
        #
        # # Add vertical black arrows to show times of irradiation
        irradiation_times_cells = [irradiation[0] + i * irradiation[1] for i in
                                   range(0, irradiation[2])]  # times of irradiation in hours
        shift = 10000
        shift2 = 3
        for time in irradiation_times_cells:
            id = np.where(times == time)[0][0]  # get index of time
            arrow1 = axes[0].annotate('', xy=(times[id], number_cells[id] + shift),
                                      xytext=(times[id], number_cells[id] + shift + 1),
                                      arrowprops=dict(facecolor='black', width=1.5, headwidth=5))
            arrow1.set_zorder(-1)  # set arrow below plot line
            arrow2 = axes[1].annotate('', xy=(times[id], tumor_size[id] + shift2),
                                      xytext=(times[id], tumor_size[id] + shift2 + 1),
                                      arrowprops=dict(facecolor='black', width=1.5, headwidth=5))
            arrow2.set_zorder(-1)  # set arrow below plot line

        #add a tiny text box in the corner with the repo name
        plt.figtext(0.01, 0.01, repo, wrap=True, horizontalalignment='left', fontsize=6)
        plt.tight_layout()
        #create a folder for the images if it doesn't exist
        if not os.path.exists(repo+'/Plots/evolution'):
            os.makedirs(repo+'/Plots/evolution')
        #save the image
        plt.savefig(repo+'/Plots/evolution/t'+str(t)+'.png', dpi=dpi)
        plt.close()

        if show_vessels:
            fig, axes = plt.subplots(1, 1, figsize=(8, 5), dpi=dpi)
            axes.plot(times, number_vessels, 'o', markersize=1, alpha=0.8, label='vessels')
            axes.plot(times_cut, number_vessels_cut, '-', markersize=3, alpha=0.8)
            axes.scatter(times_cut[-1], number_vessels_cut[-1], s=50, alpha=0.8)
            axes.set_title('Number of Vessels Evolution', fontsize = 14)
            axes.set_xlabel('Time [h]', fontsize = 16)
            axes.set_ylabel('Number of Vessels', fontsize = 16)
            axes.grid(True)
            axes.legend(fontsize = 14)
            plt.tight_layout()
            if not os.path.exists(repo + '/Plots/vessels_evolution'):
                os.makedirs(repo + '/Plots/vessels_evolution')
            plt.savefig(repo+'/Plots/vessels_evolution/t'+str(t)+'.png', dpi=dpi)
            plt.close()

if generate_gif:
    image_directory = repo + '/Plots/CurrentPlotting'
    image_sufix1 = '_AllPlots'
    output_path1 = repo + '/' + image_sufix1 + '.gif'
    image_sufix2 = '_Vasculature'
    output_path2 = repo + '/' + image_sufix2 + '.gif'
    image_directory_evol = repo + '/Plots/evolution'
    image_directory_vessel_evol = repo + '/Plots/vessels_evolution'

    # create_gif(image_directory, output_path1, image_sufix1, image_step)
    # create_gif(image_directory, output_path2, image_sufix2, image_step)
    create_gif(image_directory_evol, repo + '/evolution.gif', '', image_step)
    create_gif(image_directory_vessel_evol, repo + '/vessels_evolution.gif', '', image_step)

#Provide the directory containing the images and the output path for the GIF

# #plot the data and the fitted curve
# fig, axs = plt.subplots(2, 1, figsize=(14, 10))
# fig.set_dpi(300)
# # Plot number of tumor cells
# axs[0].plot(times, number_cells, 'o', linewidth=2, alpha=0.8)
# # axs[0].plot(times, func(times, *popt1), 'k--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt1))
# axs[0].set_xlabel('Time (days)')
# axs[0].set_ylabel('Number of tumor cells')
# axs[0].set_title('Number of tumor cells over time')
# axs[0].set_facecolor('whitesmoke')
# axs[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
#
# # Plot tumor size
# axs[1].plot(times, tumor_size, 'o', linewidth=2, alpha=0.8)
# # axs[1].plot(times, func(times, *popt2), 'k--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt2))
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
#
# # Add vertical black arrows to show times of irradiation
# irradiation_times_cells = [72, 96, 120, 144, 168]  # times of irradiation in hours
# #read that from the config file
#
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
#
# legend_elements = [Line2D([0], [0], marker='>', color='black', lw=0, label='Irradiation with 6MeV photon beam')]
#
# # Add the custom legend to both subplots
# axs[0].legend(handles=legend_elements, loc='upper left')
# axs[1].legend(handles=legend_elements, loc='upper left')
#
# plt.tight_layout()
# plt.savefig('tumorgrowth.png')
# plt.show()



