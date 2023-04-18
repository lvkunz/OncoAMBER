import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit


number_cells = np.load('number_tumor_cells.npy')
tumor_size = np.load('tumor_size.npy')

plt.plot(number_cells)
plt.show()
plt.plot(tumor_size)
plt.show()

radius = (tumor_size/4)**(1/3)

#fit the data to an exponential curve
def func(x, a, b, c):
    return a* np.exp(b * x)

#fit the data
popt1, pcov1 = curve_fit(func, np.array(range(0, 200, 5)), number_cells, p0=(1, 0.01, 1))

#fit the data
popt2, pcov2 = curve_fit(func, np.array(range(0, 200, 5)), tumor_size, p0=(1, 0.01, 1))

popt3, pcov3 = curve_fit(func, np.array(range(0, 200, 5)), radius, p0=(1, 0.01, 1))


#plot the data and the fitted curve

#print the parameters
print('Parameters for number of tumor cells:')
print('a = ', popt1[0])
print('b = ', popt1[1])
#doubling time
print('Doubling time: ', np.log(2)/popt1[1])
print('c = ', popt1[2])

print('Parameters for tumor size:')
print('a = ', popt2[0])
print('b = ', popt2[1])
#doubling time
print('Doubling time: ', np.log(2)/popt2[1])
print('c = ', popt2[2])

print('Parameters for tumor radius:')
print('a = ', popt3[0])
print('b = ', popt3[1])
#doubling time
print('Doubling time: ', np.log(2)/popt3[1])
print('c = ', popt3[2])





fig, axs = plt.subplots(2, 1, figsize=(14, 10))
fig.set_dpi(300)
# Plot number of tumor cells
axs[0].plot(np.array(range(0, 200, 5)), number_cells, 'o', linewidth=2, alpha=0.8)
axs[0].plot(np.array(range(0, 200, 5)), func(np.array(range(0, 200, 5)), *popt1), 'k--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt1))
axs[0].set_xlabel('Time (days)')
axs[0].set_ylabel('Number of tumor cells')
axs[0].set_title('Number of tumor cells over time')
# axs[0].set_ylim(0, 20000)
axs[0].set_xlim(0, 200)
axs[0].set_facecolor('whitesmoke')
axs[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.5)

# Plot tumor size
axs[1].plot(np.array(range(0, 200, 5)), tumor_size, 'o', linewidth=2, alpha=0.8)
axs[1].plot(np.array(range(0, 200, 5)), func(np.array(range(0, 200, 5)), *popt2), 'k--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt2))
axs[1].set_xlabel('Time (days)')
axs[1].set_ylabel('Tumor volume [mm^3]')
axs[1].set_title('Tumor volume over time')
# axs[1].set_ylim(0, 100)
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

plt.tight_layout()
plt.savefig('tumorgrowth.png')
plt.show()
