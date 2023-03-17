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
Topas = False
CellDynamics = True
Vasculature = True



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



if Vasculature:
    # world.generate_vasculature(2000)
    # world.vasculature.save('Vasculature/vasculature_current.txt')
    world.vasculature.read('Vasculature/vasculature_current.txt')
    world.compute_oxygen_map()
    figO = plt.figure()
    axO = figO.add_subplot(111, projection='3d')
    #world.show_voxels_centers_molecules(axO, figO, 'VEGF')
    world.show_voxels_centers_oxygen(axO, figO)
    world.vasculature.plot(figO, axO)
    plt.title('Oxygen in voxels')
    plt.show()

# fig3 = plt.figure()
# ax3 = fig3.add_subplot(111, projection='3d')
# ax3.figure.set_dpi(DPI)
# # ax3.view_init(90, 0)
# world.show_voxels_centers_dose(ax3, fig3)
# plt.title('Dose in voxels')
# plt.savefig('Plots/dose_' + param_file + '.png')
# plt.show()

#########################################################################################

if CellDynamics:

    initial_number_cells = int(100000)

    for i in range(initial_number_cells):
        if i % 100000 == 0: print('Adding healthy cell number: ', i)
        voxel1 = world.find_voxel(np.random.uniform(-15, 15, 3))
        voxel1.add_cell(HealthyCell(0.01, cycle_hours=30, life_expectancy=5000, color='my green'))

    # for i in range(100):
    #     print('Adding tumor cell number: ', i)
    #     voxel1 = world.find_voxel(np.random.uniform(-1, 1, 3))
    #     voxel1.add_cell(TumorCell(0.01, cycle_hours=20, life_expectancy=100000000, color='my purple'))

    # for i in centre_voxel_numbers:
    #     world.voxel_list[i].add_cell(Cell(0.003, cycle_hours=30, life_expectancy=100, color='red'))

    if show:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # set dpi
        ax.figure.set_dpi(DPI)
        # ax.view_init(90, 0)
        world.show_voxels_centers(ax, fig)
        plt.title('Initial cells in voxels')
        plt.savefig('Plots/initial.png')
        plt.show()

        figP = plt.figure()
        axP = figP.add_subplot(111, projection='3d')
        axP.figure.set_dpi(DPI)
        # axP.view_init(90, 0)
        world.show_voxels_centers_pressure(axP, figP)
        plt.title('Pressure in voxels')
        plt.show()


    simulation_start = time.time()



    ##########################################################################################

    end_time = 200 # in hours (simulation time)
    dt = 20 # in hours (time step)

    # if show:
        # matrix = world.compute_exchange_matrix(dt)
        # matrix = matrix / np.max(matrix)
        # #matrix = matrix[0:100, 0:100]
        # plt.figure()
        # plt.imshow(matrix)
        # plt.title('Initial Exchange matrix')
        # plt.colorbar()
        # plt.show()

    celldivision = CellDivision('cell_division', dt)
    cellapoptosis = CellApoptosis('cell_apoptosis', dt)
    cellaging = CellAging('cell_aging', dt)
    cellmigration = CellMigration('cell_migration', dt)
    update_cell_init_state = UpdateCellInitialState('update_cell_state', dt)
    update_cell_state = UpdateCellState('update_cell_state', dt)

    list_of_processes = [cellapoptosis, cellaging, cellmigration, celldivision, update_cell_init_state, update_cell_state]

    sim = Simulator(list_of_processes, end_time, dt)
    sim.run(world, video=False)

    simulation_end = time.time()
    ##########################################################################################

    # total number of cell
    total_number_of_cells = 0
    cells_cycling = 0
    cells_apoptotic = 0
    cells_senescent = 0
    for voxel in world.voxel_list:
        total_number_of_cells += len(voxel.list_of_cells)
        for cell in voxel.list_of_cells:
            if cell.state == 'cycling':
                cells_cycling += 1
            elif cell.state == 'apoptotic':
                cells_apoptotic += 1
            elif cell.state == 'senescent':
                cells_senescent += 1

    print('total number of cells: ', total_number_of_cells)
    print('number of cells cycling: ', cells_cycling)
    print('number of cells apoptotic: ', cells_apoptotic)
    print('number of cells senescent: ', cells_senescent)

    print('ratio of cells: ', total_number_of_cells / (initial_number_cells+1))
    print('simulation time: ', simulation_end - simulation_start, ' seconds')

    # if show:
        # matrix = world.compute_exchange_matrix(dt)
        # #matrix = matrix[0:100, 0:100]
        # plt.figure()
        # plt.imshow(matrix)
        # plt.title('Final Exchange matrix')
        # plt.colorbar()
        # plt.show()

    if show:

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')
        ax2.figure.set_dpi(DPI)
        # view from above
        # ax2.view_init(90, 0)
        world.show_voxels_centers(ax2, fig2)
        plt.title('Final cells in voxels at time t = ' + str(end_time) + ' hours')
        plt.savefig('Plots/final.png')
        plt.show()

        figP2 = plt.figure()
        axP2 = figP2.add_subplot(111, projection='3d')
        axP2.figure.set_dpi(DPI)
        # axP.view_init(90, 0)
        world.show_voxels_centers_pressure(axP2, figP2)
        plt.title('Pressure in voxels')
        plt.show()

print('total time: ', time.time() - start_time, ' seconds')
