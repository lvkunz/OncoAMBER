import sys
sys.path.insert(0, '/PHShome/lk001/.conda/envs/amberenv/lib/python3.9/site-packages') #cluster
import amber
import numpy as np
import random
import time
import os
import pyvista as pv

print('Current working directory:', os.getcwd())
#print the directory of amber
print('Amber directory:', amber.__file__)
print('Amber version:', amber.__version__)

if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    raise ValueError('No config file specified')

config = amber.Config(config_file)
print('Config file')
print(config)

def show_tumor_3D_solid(world, t):

    grid_shape = (world.number_of_voxels, world.number_of_voxels, world.number_of_voxels)
    cell_counts = np.empty(grid_shape)
    necro_counts = np.empty(grid_shape)

    # Fill the array with the number of tumor cells in each voxel
    for voxel in world.voxel_list:
        i, j, k = np.unravel_index(voxel.voxel_number, grid_shape)
        cell_counts[i, j, k] = voxel.number_of_tumor_cells()
        necro_counts[i, j, k] = voxel.number_of_necrotic_cells()


    # Create a vtkImageData object and assign the cell counts to it
    grid = pv.UniformGrid()
    grid.dimensions = grid_shape
    hl = world.half_length
    voxel_side = 2*hl/world.number_of_voxels
    grid.origin = (-hl, -hl, -hl)  # The bottom left corner of the data set
    grid.spacing = (voxel_side, voxel_side, voxel_side)  # These are the cell sizes along each axis
    grid.point_data['tumor'] = cell_counts.flatten(order="F")  # Flatten the array!


    contour_values = [1, 100, 300, 500, 800, 1000]
    max = np.max(cell_counts)
    #remove values below max from contour_values
    contour_values = [x for x in contour_values if x < max]

    # Create a Plotter object
    plotter = pv.Plotter()
    plotter.add_mesh(grid.outline_corners(), color='k')

    for i, value in enumerate(contour_values):
        opacity = 0.3 + 0.5 * i / len(contour_values)
        contour = grid.contour([value])
        plotter.add_mesh(contour, cmap='Blues', opacity= opacity, scalars='tumor')

    for vessel in world.vasculature.list_of_vessels:
        if vessel.visible:
            path = vessel.path
            if len(path) > 1:
                origin = path[0]
                end = path[-1]
                line = pv.Line(origin, end)
                plotter.add_mesh(line, color='crimson', line_width=1)

    #plotter.export_html('3D_plot_'+str(t)+'.html', 'panel')
    plotter.show(screenshot='3D_plot_'+str(t)+'.png', window_size=(1000, 1000), auto_close=True)


#set seed for reproducibility
start_time = time.time()
DPI = 100

seed = config.seed
np.random.seed(seed)
random.seed(seed)

print('python version', sys.version)
print('Config file', config_file)
print('Seed', seed)
print('#'*80)
for key, value in config.__dict__.items():
    print(key, value)
print('#'*80)


world = amber.World(config)

#########################################################################################

#add cells to the voxels (Normal stroma cells and tumor cells)
for i in range(world.total_number_of_voxels):
    if i %10000 == 0: print('Adding healthy cells to voxel number: ', i, ' out of ', world.total_number_of_voxels)
    for j in range(config.initial_number_healthy_cells):
        cell = amber.Cell(config.radius_healthy_cells, cycle_hours=config.doubling_time_healthy, cycle_std=config.doubling_time_sd, intra_radiosensitivity=config.intra_radiosensitivity, o2_to_vitality_factor=config.o2_to_vitality_factor, type='NormalCell')
        cell.time_spent_cycling = 0
        world.voxel_list[i].add_cell(cell, config.max_occupancy)

points = amber.Sphere(config.tumor_initial_radius, [0, 0, 0]).generate_random_points(config.initial_number_tumor_cells)
for i in range(config.initial_number_tumor_cells):
    if i % 10000 == 0: print('Adding tumor cells ', i, ' out of ', config.initial_number_tumor_cells)
    voxel = world.find_voxel(points[i])
    voxel.add_cell(
        amber.Cell(config.radius_tumor_cells, cycle_hours=config.doubling_time_tumor, cycle_std=config.doubling_time_sd, intra_radiosensitivity=config.intra_radiosensitivity, o2_to_vitality_factor=config.o2_to_vitality_factor, type='TumorCell'), config.max_occupancy)


#generate vasculature and print related information
world.generate_healthy_vasculature(config.vessel_number,
            splitting_rate=0.6,
            mult_macro_steps=2.0,
            micro_steps=30,
            weight_direction=1.5,
            weight_vegf=0.9,
            weight_pressure=0.0,
            )
world.update_volume_occupied_by_vessels()
print('Relative volume occupied by vessels, ratio: ', 100*(world.measure_vasculature_volume()/(world.half_length*2)**3), '%')
print('Length of vasculature: ', 100*(world.measure_vasculature_length()/(world.half_length*2)**3), 'mm/mm^3')
print('Area of vasculature: ', 10*(world.measure_vasculature_area()/(world.half_length*2)**3), 'mm^2/mm^3')
world.update_capillaries(n_capillaries_per_VVD=config.n_capillaries_per_VVD, capillary_length=config.capillary_length)

show_tumor_3D_solid(world, 0)

##########################################################################################

#prepare the simulation
end_time = config.endtime
dt = config.dt

celldivision = amber.CellDivision( config, 'cell_division', dt,
                                        cycling_threshold=config.vitality_cycling_threshold,
                                        pressure_threshold=config.max_occupancy)

celldeath = amber.CellDeath(config, 'cell_death', dt,
                                        apoptosis_threshold=config.vitality_apoptosis_threshold,
                                        apoptosis_probability=config.probability_apoptosis,
                                        necrosis_threshold=config.vitality_necrosis_threshold,
                                        necrosis_probability=config.probability_necrosis)

cellaging = amber.CellAging(config, 'cell_aging', dt)

cellmigration = amber.CellMigration(config, 'cell_migration', dt)

update_cell_state = amber.UpdateCellOxygen(config, 'update_cell_state', dt,
                                        voxel_half_length=(config.half_length_world/config.voxel_per_side))

update_molecules = amber.UpdateVoxelMolecules(config, 'update_molecules', dt,
                                        VEGF_production_per_cell=config.VEGF_production_per_cell,
                                        threshold_for_VEGF_production=config.o2_threshold_for_VEGF_production)

update_vessels = amber.UpdateVasculature(config, 'update_vessels', dt,
                                        killing_radius_threshold=config.radius_killing_threshold,
                                        killing_length_threshold=config.length_killing_threshold,
                                        n_capillaries_per_VVD=config.n_capillaries_per_VVD,
                                        capillary_length=config.capillary_length,
                                        splitting_rate=config.splitting_rate_vasculature,
                                        macro_steps=config.macro_steps,
                                        micro_steps=config.micro_steps,
                                        weight_direction=config.weight_direction,
                                        weight_vegf=config.weight_vegf,
                                        weight_pressure=config.weight_pressure,
                                        radius_pressure_sensitive=config.radius_pressure_sensitive)

list_of_processes = [update_cell_state, celldivision, celldeath, update_molecules, cellaging, cellmigration, update_vessels]

#run the simulation and time it

simulation_start = time.time()

sim = amber.Simulator(config, list_of_processes, end_time, dt)
sim.run(world, video=config.show_time_steps)

simulation_end = time.time()


##########################################################################################

print('simulation time: ', simulation_end - simulation_start, ' seconds')
print('total time: ', time.time() - start_time, ' seconds')

show_tumor_3D_solid(world,config.endtime)