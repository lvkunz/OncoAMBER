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

#set seed for reproducibility
seed = random.randint(0, 1000000)
np.random.seed(seed)
print('seed: ', seed)

start_time = time.time()

DPI = 100

show = True
Topas = False
CellDynamics = True
Vasculature = True
Generate = False

#define dictionnary containing all the parameters
CONFIG = read_config_file('CONFIG.txt')

world = World(CONFIG['half_length_world'], CONFIG['voxel_per_side'])


if Topas:
    world.topas_param_file(param_file)
    RunTopasSimulation(param_file)
    os.rename('TopasSimulation/MyScorer.csv', 'TopasSimulation/' + CONFIG['TOPAS_file'] + '.csv')
    n, doses = DoseOnWorld('TopasSimulation/' + CONFIG['TOPAS_file'] + '.csv')
    world.update_dose(doses)
# read in the csv file


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
        for j in range(CONFIG['initial_number_cells']):
            world.voxel_list[i].add_cell(HealthyCell(CONFIG['radius_healthy'], cycle_hours=CONFIG['doubling_time']))

    points = Sphere(0.5, [0, 0, 0]).generate_random_points(CONFIG['initial_number_tumor_cells'])
    for i in range(CONFIG['initial_number_tumor_cells']):
        if i %10000 == 0: print('Adding tumor cells ', i, ' out of ', CONFIG['initial_number_tumor_cells'])
        voxel = world.find_voxel(points[i])
        voxel.add_cell(TumorCell(CONFIG['radius_tumor'], cycle_hours=CONFIG['doubling_time_tumor']))

    world.generate_healthy_vasculature(CONFIG['vessel_number'])
    world.vasculature.save_vessels('Vasculature/vasculature_new2')
    world.update_volume_occupied_by_vessels()
    world.update_oxygen(diffusion_number=CONFIG['diffusion_number'])

    if Vasculature:
        # Sphere = Sphere(2.0, [0, 0, 0])
        # points = Sphere.generate_random_points_on_surface(CONFIG['vessel_number'])
        # points2 = [points[i] - points[i]/20.0 for i in range(len(points))]
        # size = world.half_length

        # points_x = np.random.uniform(-size, size, PARAMETERS['vessel_number'])
        # points_y = np.random.uniform(-size, size, PARAMETERS['vessel_number'])
        # points_z = -size * np.ones(PARAMETERS['vessel_number'])
        # points = []
        # for i in range(len(points_x)):
        #     points.append([points_x[i], points_y[i], points_z[i]])
        # # create a second point towards the center
        # print(points)
        # points2 = []
        # for i in range(len(points)):
        #     # points2.append(points[i] - points[i]/20)
        #     points2.append([points[i][0], points[i][1], points[i][2] + 0.05])
        # print(points2)

        # vessels = []
        # for i in range(len(points)):
        #     vessels.append(Vessel([points[i], points2[i]], 1))

        # world.initiate_vasculature(vessels)

        # world.read_vasculature('Vasculature/vasculature_current.txt')
        # world.update_oxygen(diffusion_number=PARAMETERS['diffusion_number'])
        figO = plt.figure()
        figO.set_size_inches(10, 10)
        axO = figO.add_subplot(111, projection='3d')
        axO.figure.set_dpi(100)
        # world.show_voxels_centers_molecules(axO, figO, 'VEGF')
        axO.set_xlim([-world.half_length, world.half_length])
        axO.set_ylim([-world.half_length, world.half_length])
        axO.set_zlim([-world.half_length, world.half_length])
        # world.show_voxels_centers_molecules(axO, figO, 'VEGF')
        world.vasculature.plot(figO, axO)
        plt.title('Vasculature')
        plt.show()

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


    end_time = CONFIG['endtime'] # in hours (simulation time)
    dt = CONFIG['dt'] # in hours (time step)

    celldivision = CellDivision('cell_division', dt, CONFIG['vitality_cycling_threshold'], CONFIG['pressure_threshold_division'])
    cellapoptosis = CellApoptosis('cell_apoptosis', dt, CONFIG['vitality_apoptosis_threshold'])
    cellaging = CellAging('cell_aging', dt)
    cellmigration = CellMigration('cell_migration', dt, pressure_threshold=CONFIG['pressure_threshold_migration'])
    update_cell_state = UpdateCellState('update_cell_state', dt, CONFIG['spread_gaussian_o2'])
    update_molecules = UpdateVoxelMolecules('update_molecules', dt, CONFIG['VEGF_production_per_cell'], CONFIG['threshold_for_VEGF_production'])
    update_vessels = UpdateVasculature('update_vessels', dt, CONFIG['pressure_threshold_death'], CONFIG['o2_per_volume'], CONFIG['diffusion_number'], CONFIG['splitting_rate'], CONFIG['macro_steps'], CONFIG['micro_steps'], CONFIG['weight_direction'], CONFIG['weight_vegf'], CONFIG['pressure_threshold_death'], CONFIG['weight_pressure'], CONFIG['radius_pressure_sensitive'])


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
