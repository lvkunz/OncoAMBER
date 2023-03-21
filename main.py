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

#define dictionnary containing all the parameters
PARAMETERS = dict()

PARAMETERS['half_length_world'] = 15
PARAMETERS['voxel_per_side'] = 50

PARAMETERS['dt'] = 10
PARAMETERS['endtime'] = 100
PARAMETERS['TOPAS_file'] = 'nobeam'

PARAMETERS['vessel_number'] = 2000

PARAMETERS['initial_number_cells'] = 100000
PARAMETERS['initial_number_tumor_cells'] = 1000
PARAMETERS['doubling_time'] = 12
PARAMETERS['doubling_time_tumor'] = 6
PARAMETERS['radius_tumor'] = 0.013
PARAMETERS['radius_healthy'] = 0.0013

PARAMETERS['pressure_threshold_migration'] = 0

PARAMETERS['vitality_o2_threshold'] = 0.5
PARAMETERS['vitality_o2_slope'] = 30

PARAMETERS['vitality_cycling_threshold'] = 0.5
PARAMETERS['vitality_apoptosis_threshold'] = 0.1


PARAMETERS['VEGF_production_per_cell'] = 0.1

PARAMETERS['pressure_threshold_death'] = 30000
PARAMETERS['pressure_threshold_slope'] = 0.0001
PARAMETERS['Vasculature Growth Rate'] = 1.0

PARAMETERS['o2_per_volume'] = 10000

PARAMETERS['diffusion_number'] = 10

world = World(PARAMETERS['half_length_world'], PARAMETERS['voxel_per_side'])


if Topas:
    world.topas_param_file(param_file)
    RunTopasSimulation(param_file)
    os.rename('TopasSimulation/MyScorer.csv', 'TopasSimulation/' + PARAMETERS['TOPAS_file'] + '.csv')

# read in the csv file
n, doses = DoseOnWorld('TopasSimulation/' + PARAMETERS['TOPAS_file'] + '.csv')
world.update_dose(doses)

if Vasculature:
    # world.generate_vasculature(PARAMETERS['vessel_number'])
    # world.vasculature.grow_vasculature(Cube(world.half_length, [0, 0, 0]).generate_random_points(3000))
    # world.vasculature.save('Vasculature/vasculature_current.txt')

    world.read_vasculature('Vasculature/vasculature_current.txt')
    world.compute_oxygen_map()
    figO = plt.figure()
    axO = figO.add_subplot(111, projection='3d')
    #world.show_voxels_centers_molecules(axO, figO, 'VEGF')
    axO.set_xlim([-world.half_length, world.half_length])
    axO.set_ylim([-world.half_length, world.half_length])
    axO.set_zlim([-world.half_length, world.half_length])
    #world.show_voxels_centers_molecules(axO, figO, 'VEGF')
    world.vasculature.plot(figO, axO)
    plt.title('Oxygen in voxels')
    plt.show()

if Topas:
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
    for i in range(PARAMETERS['initial_number_cells']):
        if i % 100000 == 0: print('Adding healthy cell number: ', i)
        voxel1 = world.find_voxel(np.random.uniform(-world.half_length, world.half_length, 3))
        voxel1.add_cell(HealthyCell(PARAMETERS['radius_healthy'], cycle_hours=PARAMETERS['doubling_time'], life_expectancy=5000, color='my green'))

    for i in range(PARAMETERS['initial_number_tumor_cells']):
        if i % 10000: print('Adding tumor cell number: ', i)
        voxel1 = world.find_voxel(np.random.uniform(-1, 1, 3))
        voxel1.add_cell(TumorCell(PARAMETERS['radius_tumor'], cycle_hours=PARAMETERS['doubling_time_tumor'], life_expectancy=100000000, color='my purple'))

    # for i in centre_voxel_numbers:
    #     world.voxel_list[i].add_cell(Cell(0.003, cycle_hours=30, life_expectancy=100, color='red'))

    if show:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # set dpi
        ax.figure.set_dpi(DPI)
        ax.view_init(0, 0)
        world.show_voxels_centers(ax, fig, True)
        plt.title('Initial cells in voxels')
        plt.savefig('Plots/initial.png')
        plt.show()

    simulation_start = time.time()



    ##########################################################################################

    end_time = PARAMETERS['endtime'] # in hours (simulation time)
    dt = PARAMETERS['dt'] # in hours (time step)

    celldivision = CellDivision('cell_division', dt, PARAMETERS['vitality_cycling_threshold'])
    cellapoptosis = CellApoptosis('cell_apoptosis', dt, PARAMETERS['vitality_apoptosis_threshold'])
    cellaging = CellAging('cell_aging', dt)
    cellmigration = CellMigration('cell_migration', dt, pressure_threshold=PARAMETERS['pressure_threshold_migration'])
    update_cell_state = UpdateCellState('update_cell_state', dt, vitality_o2_threshold=PARAMETERS['vitality_o2_threshold'], vitality_slope=PARAMETERS['vitality_o2_slope'])
    update_molecules = UpdateVoxelMolecules('update_molecules', dt, PARAMETERS['VEGF_production_per_cell'])
    update_vessels = UpdateVasculature('update_vessels', dt, PARAMETERS['pressure_threshold_death'], PARAMETERS['pressure_threshold_slope'], PARAMETERS['Vasculature Growth Rate'], PARAMETERS['o2_per_volume'], PARAMETERS['diffusion_number'])

    list_of_processes = [update_molecules, update_vessels, update_cell_state, cellaging, cellapoptosis, celldivision, cellmigration]

    #world.update_biology_after_RT()

    sim = Simulator(list_of_processes, end_time, dt)
    sim.run(world, video=True)

    simulation_end = time.time()
    ##########################################################################################

    # total number of cell
    total_number_of_cells = 0
    cells_cycling = 0
    cells_apoptotic = 0
    cells_senescent = 0

    print(' -- Computing final number of cells')
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
        world.vasculature.plot(fig2, ax2)
        plt.title('Final cells in voxels at time t = ' + str(end_time) + ' hours')
        plt.savefig('Plots/final.png')
        plt.show()


print('total time: ', time.time() - start_time, ' seconds')
