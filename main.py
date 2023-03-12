from Voxel import Voxel
from Cell import Cell, TumorCell, HealthyCell
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from World import World
from Process import *
#time the simulation
import time
start_time = time.time()

# voxel1 = Voxel(np.array([0,0,0]))
#
# for i in range(100):
#     voxel1.add_cell(Cell(0.01, color = 'red'))
# for i in range(100):
#     voxel1.add_cell(Cell(0.01))
#
# #fig,ax = voxel1.plot_vox_and_cells()
# ax.set_xlim(-0.2,0.2)
# ax.set_ylim(-0.2,0.2)
# ax.set_zlim(-0.2,0.2)
#
# plt.show()

DPI = 100

show = True

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#set dpi
ax.figure.set_dpi(DPI)
#ax.view_init(90, 0)
initial_number_cells = 10000
ratio = []
stable_for_le100 = 42.5

half_length_world = 10
voxels_per_side = 10
world = World(half_length_world, voxels_per_side)

central_voxel = world.find_voxel(np.array([0,0,0])).voxel_number

#find eight corners voxel numbers

for i in range(initial_number_cells):
    #print('Adding healthy cell number: ', i)
    voxel1 = world.find_voxel(np.random.uniform(-10,10,3))
    voxel1.add_cell(HealthyCell(0.01, cycle_hours=stable_for_le100, life_expectancy=100, color='my green'))

for i in range(100):
    print('Adding tumor cell number: ', i)
    voxel1 = world.find_voxel(np.random.uniform(-1,1,3))
    voxel1.add_cell(TumorCell(0.01, cycle_hours=stable_for_le100*0.1, life_expectancy=100, color='my purple'))

# for i in centre_voxel_numbers:
#     world.voxel_list[i].add_cell(Cell(0.003, cycle_hours=30, life_expectancy=100, color='red'))

if show:
    world.show_voxels_centers(ax,fig,colorful=True)
    plt.title('Initial cells in voxels')
    plt.savefig('Plots/initial.png')
    plt.show()


simulation_start = time.time()

##########################################################################################
end_time = 50
dt = 5

celldivision = CellDivision('cell_division', dt)
cellapoptosis = CellApoptosis('cell_apoptosis',dt)
cellaging = CellAging('cell_aging',dt)
cellmigration = CellMigration('cell_migration',dt)

list_of_processes = [celldivision, cellapoptosis, cellaging, cellmigration]


print('starting number of cells: ', len(world.voxel_list[central_voxel].list_of_cells))

sim = Simulator(list_of_processes,end_time,dt)
sim.run(world)

print('ending number of cells: ', len(world.voxel_list[central_voxel].list_of_cells))

simulation_end = time.time()

#total number of cell
total_number_of_cells = 0
for voxel in world.voxel_list:
    total_number_of_cells += len(voxel.list_of_cells)
print('total number of cells: ', total_number_of_cells)
print('ratio of cells: ', total_number_of_cells/initial_number_cells)
ratio.append(total_number_of_cells/initial_number_cells)
print('simulation time: ', simulation_end - simulation_start, ' seconds')



if show:
    print('test')
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.figure.set_dpi(DPI)
    #view from above
    #ax2.view_init(90, 0)
    world.show_voxels_centers(ax2,fig2,colorful=True)
    plt.title('Final cells in voxels at time t = ' + str(end_time) + ' hours')
    plt.savefig('Plots/final.png')
    plt.show()


print('total time: ', time.time() - start_time, ' seconds')
