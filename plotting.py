import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


number_cells = np.load('number_tumor_cells.npy')
tumor_size = np.load('tumor_size.npy')

fig, axs = plt.subplots(2, 1, figsize=(14, 10))
fig.set_dpi(300)
# Plot number of tumor cells
axs[0].plot(np.array(range(0, 200, 4)), number_cells, color='purple', linewidth=2, alpha=0.8)
axs[0].set_xlabel('Time (days)')
axs[0].set_ylabel('Number of tumor cells')
axs[0].set_title('Number of tumor cells over time')
axs[0].set_ylim(0, 20000)
axs[0].set_xlim(0, 200)
axs[0].set_facecolor('whitesmoke')
axs[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.5)

# Plot tumor size
axs[1].plot(np.array(range(4, 204, 4)), tumor_size, color='blue', linewidth=2, alpha=0.8)
axs[1].set_xlabel('Time (days)')
axs[1].set_ylabel('Tumor volume [mm^3]')
axs[1].set_title('Tumor volume over time')
axs[1].set_ylim(4, 37)
axs[1].set_xlim(0, 200)
axs[1].set_facecolor('whitesmoke')
axs[1].grid(True, linestyle='--', linewidth=0.5, alpha=0.5)

# Set the x-axis ticks to show 1 day, 2 days, 3 days, etc.
xticks = [0, 24, 48, 72, 96, 120, 144, 168, 192]
xticklabels = [i for i in range(9)]
axs[0].set_xticks(xticks)
axs[0].set_xticklabels(xticklabels)
axs[1].set_xticks(xticks)
axs[1].set_xticklabels(xticklabels)

# Add vertical black arrows to show times of irradiation
irradiation_times_cells = [72, 96, 120, 144, 168]  # times of irradiation in hours
for time in irradiation_times_cells:
    t = time  # convert hours to days
    id = np.where(np.array(range(0, 200, 4)) == t)[0][0]  # get index of time
    arrow = axs[0].annotate('', xy=(t, number_cells[id] + 300), xytext=(t, number_cells[id] + 301), arrowprops=dict(facecolor='black', width=1.5, headwidth=5))
    arrow.set_zorder(-1)  # set arrow below plot line

# Add vertical arrows to the tumor size plot
irradiation_times_size = [72, 96, 120, 144, 168]  # times of irradiation in hours
for time in irradiation_times_size:
    t = time  # convert hours to days
    id = np.where(np.array(range(4, 204, 4)) == t)[0][0]  # get index of time
    arrow = axs[1].annotate('', xy=(t, tumor_size[id] + 0.5), xytext=(t, tumor_size[id] + 0.6), arrowprops=dict(facecolor='black', width=1.5, headwidth=5))
    arrow.set_zorder(-1)  # set arrow below plot line

legend_elements = [Line2D([0], [0], marker='>', color='black', lw=0, label='Irradiation with 6MeV photon beam')]

# Add the custom legend to both subplots
axs[0].legend(handles=legend_elements, loc='upper left')
axs[1].legend(handles=legend_elements, loc='upper left')

plt.tight_layout()
plt.savefig('fractionation.png')
plt.show()
