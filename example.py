import sys
sys.path.insert(0, '/PHShome/lk001/.conda/envs/amberenv/lib/python3.9/site-packages') #cluster
import amber
import numpy as np
import random
import time
import os

print('Current working directory:', os.getcwd())
#print the directory of amber
print('Amber directory:', amber.__file__)
print('Amber version:', amber.__version__)

if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    raise ValueError('No config file specified')

config_dict = amber.read_config_file(config_file)
config = amber.Config(config_dict)
print('Config file')
print(config)
print(config.half_length_world)

#set seed for reproducibility
start_time = time.time()
DPI = 100

seed = config.seed
np.random.seed(seed)
random.seed(seed)

print('python version', sys.version)
print('Config file', config_file)
print('#'*80)
for key, value in config_dict.items():
    print(key, value)
print('#'*80)


world = amber.World(config)

#########################################################################################

#add cells to the voxels (Normal stroma cells and tumor cells)
for i in range(world.total_number_of_voxels):
    if i %10000 == 0: print('Adding healthy cells to voxel number: ', i, ' out of ', world.total_number_of_voxels)
    for j in range(config.initial_number_healthy_cells):
        cell = amber.Cell(config.radius_healthy_cells, cycle_hours=config.doubling_time_healthy, cycle_std=config.doubling_time_sd, radiosensitivity=config.radiosensitivity, o2_to_vitality_factor=config.o2_to_vitality_factor, type='NormalCell')
        cell.time_spent_cycling = 0
        world.voxel_list[i].add_cell(cell)

points = amber.Sphere(config.tumor_initial_radius, [0, 0, 0]).generate_random_points(config.initial_number_tumor_cells)
for i in range(config.initial_number_tumor_cells):
    if i % 10000 == 0: print('Adding tumor cells ', i, ' out of ', config.initial_number_tumor_cells)
    voxel = world.find_voxel(points[i])
    voxel.add_cell(
        amber.Cell(config.radius_tumor_cells, cycle_hours=config.doubling_time_tumor, cycle_std=config.doubling_time_sd, radiosensitivity=config.radiosensitivity, o2_to_vitality_factor=config.o2_to_vitality_factor, type='TumorCell'))


#generate vasculature and print related information
world.generate_healthy_vasculature(config.vessel_number,
            splitting_rate=0.5,
            mult_macro_steps=2.0,
            micro_steps=20,
            weight_direction=2.0,
            weight_vegf=0.9,
            weight_pressure=0.0,
            )
world.update_volume_occupied_by_vessels()
print('Relative volume occupied by vessels, ratio: ', 100*(world.measure_vasculature_volume()/(world.half_length*2)**3), '%')
print('Length of vasculature: ', 100*(world.measure_vasculature_length()/(world.half_length*2)**3), 'mm/mm^3')
print('Area of vasculature: ', 10*(world.measure_vasculature_area()/(world.half_length*2)**3), 'mm^2/mm^3')
world.update_oxygen(o2_per_volume=config.o2_per_volume, diffusion_number=config.diffusion_number)

##########################################################################################

#prepare the simulation
end_time = config.endtime
dt = config.dt

celldivision = amber.CellDivision( config, 'cell_division', dt,
                                        cycling_threshold=config.vitality_cycling_threshold,
                                        pressure_threshold=config.pressure_threshold_division)

celldeath = amber.CellDeath(config, 'cell_death', dt,
                                        apoptosis_threshold=config.vitality_apoptosis_threshold,
                                        apoptosis_probability=config.probability_apoptosis,
                                        necrosis_threshold=config.vitality_necrosis_threshold,
                                        necrosis_probability=config.probability_necrosis)

cellaging = amber.CellAging(config, 'cell_aging', dt)

cellmigration = amber.CellMigration(config, 'cell_migration', dt,
                                        pressure_threshold=config.pressure_threshold_migration)

update_cell_state = amber.UpdateCellOxygen(config, 'update_cell_state', dt,
                                        voxel_half_length=(config.half_length_world/config.voxel_per_side),
                                        effective_vessel_radius=config.effective_vessel_radius)

update_molecules = amber.UpdateVoxelMolecules(config, 'update_molecules', dt,
                                        VEGF_production_per_cell=config.VEGF_production_per_cell,
                                        threshold_for_VEGF_production=config.o2_threshold_for_VEGF_production)

update_vessels = amber.UpdateVasculature(config, 'update_vessels', dt,
                                        killing_radius_threshold=config.radius_killing_threshold,
                                        killing_length_threshold=config.length_killing_threshold,
                                        o2_per_volume=config.o2_per_volume,
                                        diffusion_number=config.diffusion_number,
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