import sys
sys.path.insert(0, '/home/users/k/kunzlo/.conda/envs/amberenv/lib/python3.8/site-packages') #cluster
import amber
import numpy as np
import random
import time
import os

print('Current working directory:', os.getcwd())
#print the directory of amber
print('Amber directory:', amber.__file__)
print('Amber version:', amber.__version__)
print('Python version:', sys.version)
print('Python path')

if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    raise ValueError('No config file specified')

config = amber.Config(config_file)

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
print('Seed', seed)

print('#'*80)

for key, value in config.__dict__.items():
    print(key, value)

print('#'*80)

if config.new_world:

    world = amber.World(config)

    #########################################################################################

    #add cells to the voxels (tumor cells)
    points = amber.Sphere(config.tumor_initial_radius, [0, 0, 0]).generate_random_points(config.initial_number_tumor_cells)
    factor = [1]
    count = 0
    max_tries = 3
    for i, point in enumerate(points):
        if count < max_tries:
            if i % 10000 == 0: print('Adding tumor cells ', i, ' out of ', config.initial_number_tumor_cells)
            voxel = world.find_voxel(point)
            doubling = factor[random.randint(0, len(factor) - 1)] * config.doubling_time_tumor
            radio_sensitivity = factor[random.randint(0, len(factor) - 1)] * config.intra_radiosensitivity
            count = 0
            was_added = False
            while not was_added and count < max_tries:
                count += 1
                was_added = voxel.add_cell(amber.TumorCell(config.radius_tumor_cells,
                                               cycle_hours=doubling,
                                               cycle_std=config.doubling_time_sd,
                                               intra_radiosensitivity=radio_sensitivity,
                                               o2_to_vitality_factor=config.o2_to_vitality_factor,
                                               VEGF_threshold = config.VEGF_production_threshold,
                                               VEGF_rate = config.VEGF_production_per_cell
                                               ), config.max_occupancy)
                if not was_added:
                    print("Failed to add a cell to the voxel, trying again, tried ", count, " times")
                    point = amber.Sphere(config.tumor_initial_radius, [0, 0, 0]).generate_random_points(1)
        else:
            print('Tried more than max times to add a cell to a voxel')
            print('Voxel is full, pressure is', voxel.pressure(), ' number of cells is', voxel.number_of_alive_cells(), ' and number of necrotic cells is', voxel.number_of_necrotic_cells())
            n = 0
            for voxel in world.voxel_list:
                n += len(voxel.list_of_cells)
            print('Total number of cells in the world is', n)
            break

    #generate vasculature and print related information
    world.generate_healthy_vasculature(config.vessel_number,
                splitting_rate=0.5,
                mult_macro_steps=1.0,
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

    # world.save(str(config.world_file) + str(config.seed) + '.pkl')

else:
    world = amber.load(str(config.world_file)+'.pkl')
    seed = world.config.seed #use the seed of the world to make sure that the same world is used
    config.seed = seed
    np.random.seed(seed)
    random.seed(seed)
    world.config = config

if config.show_3D_mesh:
    amber.show_tumor_3D_solid(world, 0)

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
                                        necrosis_probability=config.probability_necrosis,
                                        necrosis_removal_probability=config.probability_necrosis_removal,
                                        apoptosis_removal_probability=config.probability_apoptosis_removal,
                                        necrosis_damage_coeff=config.necrosis_damage_coeff,
                                        apoptosis_damage_coeff=config.apoptosis_damage_coeff)

cellaging = amber.CellAging(config, 'cell_aging', dt, repair_per_hour=config.repair_per_hour)

cellmigration = amber.CellMigration(config, 'cell_migration', dt)

update_cell_state = amber.UpdateCellOxygen(config, 'update_cell_state', dt,
                                        voxel_half_length=(config.half_length_world/config.voxel_per_side),
                                        file_prefix_alpha_beta_maps=config.file_prefix_alpha_beta_maps)

update_molecules = amber.UpdateVoxelMolecules(config, 'update_molecules', dt)

update_vessels = amber.UpdateVasculature(config, 'update_vessels', dt,
                                        killing_radius_threshold=config.radius_killing_threshold,
                                        n_capillaries_per_VVD=config.n_capillaries_per_VVD,
                                        capillary_length=config.capillary_length,
                                        splitting_rate=config.splitting_rate_vasculature,
                                        macro_steps=config.macro_steps,
                                        micro_steps=config.micro_steps,
                                        weight_direction=config.weight_direction,
                                        weight_vegf=config.weight_vegf,
                                        weight_pressure=config.weight_pressure)

cellinteraction = amber.CellInteraction(config, 'cell_interaction', dt)
#not included in the simulation

list_of_processes = [cellmigration, update_cell_state, update_molecules, cellaging, celldivision, celldeath, update_vessels]

#run the simulation and time it

simulation_start = time.time()

sim = amber.Simulator(config, list_of_processes, end_time, dt)
sim.run(world, video=config.show_time_steps)

simulation_end = time.time()


##########################################################################################

print('simulation time: ', simulation_end - simulation_start, ' seconds')
print('total time: ', time.time() - start_time, ' seconds')


# world.save('final' + str(config.world_file) + str(config.seed) + '.pkl')
if config.show_3D_mesh:
    amber.show_tumor_3D_solid(world, 0)