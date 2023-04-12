import subprocess
from Terminal import *
import os
import ReadAndWrite as rw
from World import World
import matplotlib.pyplot as plt


CONFIG = rw.read_config_file('CONFIG.txt')

world = World(CONFIG['half_length_world'], CONFIG['voxel_per_side'])
world.topas_param_file(CONFIG['TOPAS_file'])
RunTopasSimulation(CONFIG['TOPAS_file'])

os.rename('TopasSimulation/MyScorer.csv', 'TopasSimulation/' + CONFIG['TOPAS_file'] + '.csv')

# read in the csv file
n, doses = rw.DoseOnWorld('TopasSimulation/' + CONFIG['TOPAS_file'] + '.csv')
world.update_dose(doses)

#plot the simulation
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
world.show_voxels_centers_dose(ax, fig)
plt.show()
