from Voxel import Voxel
from Cell import Cell, TumorCell, HealthyCell
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from World import World
from Process import *
from Terminal import *
import os
from ReadAndWrite import *

# time the simulation
import time

start_time = time.time()

DPI = 100

show = True
Topas = True
CellDynamics = True



half_length_world = 15 # in cm
voxels_per_side = 15
world = World(half_length_world, voxels_per_side)

param_file = 'nobeam'

if Topas:
    world.topas_param_file(param_file)
    RunTopasSimulation(param_file)
    os.rename('TopasSimulation/MyScorer.csv', 'TopasSimulation/' + param_file + '.csv')

# read in the csv file
n, doses = DoseOnWorld('TopasSimulation/' + param_file + '.csv')
world.update_dose(doses)

fig3 = plt.figure()
ax3 = fig3.add_subplot(111, projection='3d')
ax3.figure.set_dpi(DPI)
# ax3.view_init(90, 0)
world.show_voxels_centers_dose(ax3, fig3)
plt.title('Dose in voxels')
plt.savefig('Plots/dose_' + param_file + '.png')
plt.show()

#########################################################################################

if CellDynamics:

    initial_number_cells = 100

    for i in range(initial_number_cells):
        # print('Adding healthy cell number: ', i)
        voxel1 = world.find_voxel(np.random.uniform(-1, 1, 3))
        voxel1.add_cell(HealthyCell(0.1, cycle_hours=40, life_expectancy=5000, color='my green'))

    # for i in range(100):
    #     # print('Adding tumor cell number: ', i)
    #     voxel1 = world.find_voxel(np.random.uniform(-1, 1, 3))
    #     voxel1.add_cell(TumorCell(0.01, cycle_hours=stable_for_le100 * 0.5, life_expectancy=10000, color='my purple'))

    # for i in centre_voxel_numbers:
    #     world.voxel_list[i].add_cell(Cell(0.003, cycle_hours=30, life_expectancy=100, color='red'))

    if show:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # set dpi
        ax.figure.set_dpi(DPI)
        # ax.view_init(90, 0)
        world.show_voxels_centers(ax, fig, colorful=True)
        plt.title('Initial cells in voxels')
        plt.savefig('Plots/initial.png')
        plt.show()

    simulation_start = time.time()

    ##########################################################################################
    end_time = 300
    dt = 20

    celldivision = CellDivision('cell_division', dt)
    cellapoptosis = CellApoptosis('cell_apoptosis', dt)
    cellaging = CellAging('cell_aging', dt)
    cellmigration = CellMigration('cell_migration', dt)
    update_cell_state = UpdateCellState('update_cell_state', dt)

    list_of_processes = [update_cell_state, cellapoptosis, cellaging, cellmigration, celldivision]

    sim = Simulator(list_of_processes, end_time, dt)
    sim.run(world, video=True)

    simulation_end = time.time()
    ##########################################################################################

    # total number of cell
    total_number_of_cells = 0
    for voxel in world.voxel_list:
        total_number_of_cells += len(voxel.list_of_cells)
    print('total number of cells: ', total_number_of_cells)

    print('ratio of cells: ', total_number_of_cells / initial_number_cells)
    print('simulation time: ', simulation_end - simulation_start, ' seconds')

    if show:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')
        ax2.figure.set_dpi(DPI)
        # view from above
        # ax2.view_init(90, 0)
        world.show_voxels_centers(ax2, fig2, colorful=True)
        plt.title('Final cells in voxels at time t = ' + str(end_time) + ' hours')
        plt.savefig('Plots/final.png')
        plt.show()

print('total time: ', time.time() - start_time, ' seconds')
