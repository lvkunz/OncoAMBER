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
from Vessel_old import *

# time the simulation
import time

#set seed for reproducibility
np.random.seed(0)

start_time = time.time()

DPI = 100

show = True
Topas = False
CellDynamics = True
Vasculature = True
Generate = False

#define dictionnary containing all the parameters
PARAMETERS = dict()

#read from a txt file
#PARAMETERS = read_parameters('parameters.txt')

PARAMETERS['half_length_world'] = 1.5
PARAMETERS['voxel_per_side'] = 60

PARAMETERS['dt'] = 10
PARAMETERS['endtime'] = 80
PARAMETERS['TOPAS_file'] = 'nobeam'

PARAMETERS['vessel_number'] = 2000

PARAMETERS['initial_number_cells'] = 10
PARAMETERS['initial_number_tumor_cells'] = 60000
PARAMETERS['doubling_time'] = 10000000000
PARAMETERS['doubling_time_tumor'] = 6
PARAMETERS['radius_tumor'] = 0.0013
PARAMETERS['radius_healthy'] = 0.0013

PARAMETERS['pressure_threshold_migration'] = 0

PARAMETERS['spread_gaussian_o2'] = 0.1

PARAMETERS['vitality_cycling_threshold'] = 0.0001 #threshold in cell vitality for cycling [0,1]
PARAMETERS['vitality_apoptosis_threshold'] = 0.0 #threshold in cell vitality for apoptosis [0,1]
PARAMETERS['pressure_threshold_division'] = 0.68 #threshold in pressure for division [0,1]


PARAMETERS['VEGF_production_per_cell'] = 0.0 #VEGF production per TumorCell per timestep

PARAMETERS['pressure_threshold_death'] = 100000 #threshold of sigmoid where vessels start to die
PARAMETERS['pressure_threshold_slope'] = 30  #slope of the sigmoid where vessels start to die
PARAMETERS['Vasculature Growth Rate'] = 10 #how many vessels are added per 1 VEGF concentration unit

PARAMETERS['o2_per_volume'] = 1000000 #oxygen concentration in voxel per volume of vessel

PARAMETERS['diffusion_number'] = 0 #number of diffusion steps per timestep

PARAMETERS['threshold_for_VEGF_production'] = 0.3 #threshold for VEGF production

world = World(PARAMETERS['half_length_world'], PARAMETERS['voxel_per_side'])


if Topas:
    world.topas_param_file(param_file)
    RunTopasSimulation(param_file)
    os.rename('TopasSimulation/MyScorer.csv', 'TopasSimulation/' + PARAMETERS['TOPAS_file'] + '.csv')

# read in the csv file
n, doses = DoseOnWorld('TopasSimulation/' + PARAMETERS['TOPAS_file'] + '.csv')
world.update_dose(doses)

if Vasculature:
    #world.generate_vasculature(PARAMETERS['vessel_number'])
    world.vasculature = VasculatureNetwork([Vessel_old([0, 0, 0], [0, 0, 0.5], 0.1)])
    world.vasculature.grow_vasculature(Cube(0.5, [0, 0, 0]).generate_random_points(80))
    world.vasculature.grow_vasculature(Sphere(0.3, [0, 0, 0]).generate_random_points(500))
    world.vasculature.save('Vasculature/vasculature_current.txt')

    #world.read_vasculature('Vasculature/vasculature_current.txt')
    world.compute_oxygen_map(diffusion_number=PARAMETERS['diffusion_number'])
    figO = plt.figure()
    figO.set_size_inches(10, 10)
    axO = figO.add_subplot(111, projection='3d')
    axO.figure.set_dpi(100)
    #world.show_voxels_centers_molecules(axO, figO, 'VEGF')
    axO.set_xlim([-world.half_length, world.half_length])
    axO.set_ylim([-world.half_length, world.half_length])
    axO.set_zlim([-world.half_length, world.half_length])
    #world.show_voxels_centers_molecules(axO, figO, 'VEGF')
    world.vasculature.plot(figO, axO)
    plt.title('Vasculature')
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
    for i in range(world.total_number_of_voxels):
        if i %10000 == 0: print('Adding healthy cells to voxel number: ', i, ' out of ', world.total_number_of_voxels)
        for j in range(PARAMETERS['initial_number_cells']):
            world.voxel_list[i].add_cell(HealthyCell(PARAMETERS['radius_healthy'], cycle_hours=PARAMETERS['doubling_time']))


    points = Sphere(0.2, [0, 0, 0]).generate_random_points(PARAMETERS['initial_number_tumor_cells'])
    for i in range(PARAMETERS['initial_number_tumor_cells']):
        if i %10000 == 0: print('Adding tumor cells ', i, ' out of ', PARAMETERS['initial_number_tumor_cells'])
        voxel = world.find_voxel(points[i])
        voxel.add_cell(TumorCell(PARAMETERS['radius_tumor'], cycle_hours=PARAMETERS['doubling_time_tumor']))


    simulation_start = time.time()
    ##########################################################################################

    # figX = plt.figure()
    # axX = figX.add_subplot(111, projection='3d')
    # axX.figure.set_dpi(DPI)
    # world.show_tumor(axX, figX)
    # plt.title('Tumor initial')
    # plt.show()


    # axX.view_init(90, 0)
    #world.show_tumor_surface()


    end_time = PARAMETERS['endtime'] # in hours (simulation time)
    dt = PARAMETERS['dt'] # in hours (time step)

    celldivision = CellDivision('cell_division', dt, PARAMETERS['vitality_cycling_threshold'], PARAMETERS['pressure_threshold_division'])
    cellapoptosis = CellApoptosis('cell_apoptosis', dt, PARAMETERS['vitality_apoptosis_threshold'])
    cellaging = CellAging('cell_aging', dt)
    cellmigration = CellMigration('cell_migration', dt, pressure_threshold=PARAMETERS['pressure_threshold_migration'])
    update_cell_state = UpdateCellState('update_cell_state', dt, PARAMETERS['spread_gaussian_o2'])
    update_molecules = UpdateVoxelMolecules('update_molecules', dt, PARAMETERS['VEGF_production_per_cell'], PARAMETERS['threshold_for_VEGF_production'])
    update_vessels = UpdateVasculature('update_vessels', dt, PARAMETERS['pressure_threshold_death'], PARAMETERS['pressure_threshold_slope'], PARAMETERS['Vasculature Growth Rate'], PARAMETERS['o2_per_volume'], PARAMETERS['diffusion_number'])

    list_of_processes = [update_molecules, update_vessels, update_cell_state, cellaging, cellapoptosis, celldivision, cellmigration]

    #world.update_biology_after_RT()

    sim = Simulator(list_of_processes, end_time, dt)
    sim.run(world, video=True)

    simulation_end = time.time()
    ##########################################################################################

    print('simulation time: ', simulation_end - simulation_start, ' seconds')

    if show:

        figY = plt.figure()
        axY = figY.add_subplot(111, projection='3d')
        axY.figure.set_dpi(DPI)
        # axX.view_init(90, 0)
        world.show_tumor(axY, figY)
        plt.title('Tumor')
        plt.savefig('Plots/tumor_final.png')
        plt.show()

        # world.show_tumor_surface()


print('total time: ', time.time() - start_time, ' seconds')
