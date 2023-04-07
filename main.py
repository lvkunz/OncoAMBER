from Voxel import Voxel
from Cell import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from World import World
from Process import *
from Terminal import *
import os
import ReadAndWrite as rw

# time the simulation
import time

#set seed for reproducibility
seed = random.randint(0, 1000000)
np.random.seed(seed)
print('seed: ', seed)

start_time = time.time()

DPI = 100

show = True
CellDynamics = True
Vasculature = True
Generate = False

#define dictionnary containing all the parameters
CONFIG = rw.read_config_file('CONFIG.txt')

world = World(CONFIG['half_length_world'], CONFIG['voxel_per_side'])

#########################################################################################

if CellDynamics:

    for i in range(world.total_number_of_voxels):
        if i %10000 == 0: print('Adding healthy cells to voxel number: ', i, ' out of ', world.total_number_of_voxels)
        for j in range(CONFIG['initial_number_healthy_cells']):
            world.voxel_list[i].add_cell(Cell(CONFIG['radius_healthy_cells'], cycle_hours=CONFIG['doubling_time_healthy'], type='HealthyCell'))

    points = Sphere(CONFIG['tumor_initial_radius'], [0, 0, 0]).generate_random_points(
        CONFIG['initial_number_tumor_cells'])
    for i in range(CONFIG['initial_number_tumor_cells']):
        if i % 10000 == 0: print('Adding tumor cells ', i, ' out of ', CONFIG['initial_number_tumor_cells'])
        voxel = world.find_voxel(points[i])
        voxel.add_cell(
            Cell(CONFIG['radius_tumor_cells'], cycle_hours=CONFIG['doubling_time_tumor'], type='TumorCell'))

    world.generate_healthy_vasculature(CONFIG['vessel_number'])
    world.vasculature.save_vessels('Vasculature/vasculature_new2')
    world.update_volume_occupied_by_vessels()
    print('Relative volume occupied by vessels, ratio: ', 100*(world.measure_vasculature_volume()/(world.half_length*2)**3), '%')
    print('Length of vasculature: ', 100*(world.measure_vasculature_length()/(world.half_length*2)**3), 'mm/mm^3')
    print('Area of vasculature: ', 10*(world.measure_vasculature_area()/(world.half_length*2)**3), 'mm^2/mm^3')
    world.update_oxygen(o2_per_volume=CONFIG['o2_per_volume'], diffusion_number=CONFIG['diffusion_number'])

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

    celldivision = CellDivision('cell_division', dt,
                                            cycling_threshold=CONFIG['vitality_cycling_threshold'],
                                            pressure_threshold=CONFIG['pressure_threshold_division'])

    cellapoptosis = CellApoptosis('cell_apoptosis', dt,
                                            apoptosis_threshold=CONFIG['vitality_apoptosis_threshold'],
                                            apoptosis_probability=CONFIG['probability_apoptosis'])

    cellnecrosis = CellNecrosis('cell_necrosis', dt,
                                            necrosis_threshold=CONFIG['vitality_necrosis_threshold'],
                                            necrosis_probability=CONFIG['probability_necrosis'])

    cellaging = CellAging('cell_aging', dt)

    cellmigration = CellMigration('cell_migration', dt,
                                            pressure_threshold=CONFIG['pressure_threshold_migration'])

    update_cell_state = UpdateCellOxygen('update_cell_state', dt,
                                            voxel_half_length=(CONFIG['half_length_world']/CONFIG['voxel_per_side']),
                                            effective_vessel_radius=CONFIG['effective_vessel_radius'])

    update_molecules = UpdateVoxelMolecules('update_molecules', dt,
                                            VEGF_production_per_cell=CONFIG['VEGF_production_per_cell'],
                                            threshold_for_VEGF_production=CONFIG['o2_threshold_for_VEGF_production'])

    update_vessels = UpdateVasculature('update_vessels', dt,
                                            pressure_killing_radius_threshold=CONFIG['pressure_radius_killing_threshold'],
                                            o2_per_volume=CONFIG['o2_per_volume'],
                                            diffusion_number=CONFIG['diffusion_number'],
                                            splitting_rate=CONFIG['splitting_rate_vasculature'],
                                            macro_steps=CONFIG['macro_steps'],
                                            micro_steps=CONFIG['micro_steps'],
                                            weight_direction=CONFIG['weight_direction'],
                                            weight_vegf=CONFIG['weight_vegf'],
                                            weight_pressure=CONFIG['weight_pressure'],
                                            radius_pressure_sensitive=CONFIG['radius_pressure_sensitive'])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.figure.set_dpi(DPI)
    update_cell_state.alpha_map.show_extra(fig, ax, [0.0, 0.7], [0, 10])
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.figure.set_dpi(DPI)
    update_cell_state.beta_map.show_extra(fig, ax, [0.0, 0.7], [0, 10])
    plt.show()



    list_of_processes = [update_cell_state, cellaging, cellnecrosis, cellapoptosis, update_molecules, celldivision, cellmigration, update_vessels]

    sim = Simulator(list_of_processes, end_time, dt)
    sim.run(world, video=True)

    simulation_end = time.time()
    ##########################################################################################

    print('simulation time: ', simulation_end - simulation_start, ' seconds')

    if show:

        figY = plt.figure()
        axY = figY.add_subplot(111, projection='3d0')
        axY.figure.set_dpi(DPI)
        # axX.view_init(90, 0)
        world.show_tumor(axY, figY)
        plt.title('Tumor')
        plt.savefig('Plots/tumor_final.png')
        plt.show()

        # world.show_tumor_surface()


print('total time: ', time.time() - start_time, ' seconds')
