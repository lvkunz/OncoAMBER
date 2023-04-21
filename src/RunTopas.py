from src.Terminal import *
import os
import src.ReadAndWrite as rw
from src.World import World
import matplotlib.pyplot as plt
from src.config_instance import config

def RunTopasSimulation():

    world = World(config.half_length_world, config.voxel_per_side)
    world.topas_param_file(config.TOPAS_file)
    RunTopasSimulation(config.TOPAS_file)

    os.rename('TopasSimulation/MyScorer.csv', 'TopasSimulation/' + config.TOPAS_file + '.csv')

    # read in the csv file
    n, doses = rw.DoseOnWorld('TopasSimulation/' + config.TOPAS_file + '.csv')
    world.update_dose(doses)

    #plot the simulation
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    world.show_voxels_centers_dose(ax, fig)
    plt.show()
