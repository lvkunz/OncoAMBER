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
seed = CONFIG['seed']
if seed == -1:
    seed = np.random.randint(0, 1000000)
np.random.seed(seed)
print('seed: ', seed)

start_time = time.time()

DPI = 100

#CONFIG.txt with all the parameters for the simulation
CONFIG = rw.read_config_file('CONFIG.txt')
world = World(CONFIG['half_length_world'], CONFIG['voxel_per_side'])

#########################################################################################

#add cells to the voxels (Normal stroma cells and tumor cells)
for i in range(world.total_number_of_voxels):
    if i %10000 == 0: print('Adding healthy cells to voxel number: ', i, ' out of ', world.total_number_of_voxels)
    for j in range(CONFIG['initial_number_healthy_cells']):
        world.voxel_list[i].add_cell(Cell(CONFIG['radius_healthy_cells'], cycle_hours=CONFIG['doubling_time_healthy'], type='HealthyCell'))

points = Sphere(CONFIG['tumor_initial_radius'], [0, 0, 0]).generate_random_points(CONFIG['initial_number_tumor_cells'])
for i in range(CONFIG['initial_number_tumor_cells']):
    if i % 10000 == 0: print('Adding tumor cells ', i, ' out of ', CONFIG['initial_number_tumor_cells'])
    voxel = world.find_voxel(points[i])
    voxel.add_cell(
        Cell(CONFIG['radius_tumor_cells'], cycle_hours=CONFIG['doubling_time_tumor'], type='TumorCell'))

#generate vasculature and print related information
world.generate_healthy_vasculature(CONFIG['vessel_number'])
world.update_volume_occupied_by_vessels()
print('Relative volume occupied by vessels, ratio: ', 100*(world.measure_vasculature_volume()/(world.half_length*2)**3), '%')
print('Length of vasculature: ', 100*(world.measure_vasculature_length()/(world.half_length*2)**3), 'mm/mm^3')
print('Area of vasculature: ', 10*(world.measure_vasculature_area()/(world.half_length*2)**3), 'mm^2/mm^3')
world.update_oxygen(o2_per_volume=CONFIG['o2_per_volume'], diffusion_number=CONFIG['diffusion_number'])

##########################################################################################

#prepare the simulation
end_time = CONFIG['endtime']
dt = CONFIG['dt']

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
                                        killing_radius_threshold=CONFIG['radius_killing_threshold'],
                                        killing_length_threshold=CONFIG['length_killing_threshold'],
                                        o2_per_volume=CONFIG['o2_per_volume'],
                                        diffusion_number=CONFIG['diffusion_number'],
                                        splitting_rate=CONFIG['splitting_rate_vasculature'],
                                        macro_steps=CONFIG['macro_steps'],
                                        micro_steps=CONFIG['micro_steps'],
                                        weight_direction=CONFIG['weight_direction'],
                                        weight_vegf=CONFIG['weight_vegf'],
                                        weight_pressure=CONFIG['weight_pressure'],
                                        radius_pressure_sensitive=CONFIG['radius_pressure_sensitive'])

list_of_processes = [update_cell_state, cellaging, cellnecrosis, cellapoptosis, update_molecules, celldivision, cellmigration, update_vessels]

#show alpha and beta maps to make sure there is no big discontinuities

#run the simulation and time it

simulation_start = time.time()

sim = Simulator(list_of_processes, end_time, dt)
sim.run(world, video=True)

simulation_end = time.time()
##########################################################################################

print('simulation time: ', simulation_end - simulation_start, ' seconds')
print('total time: ', time.time() - start_time, ' seconds')
