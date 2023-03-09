from Voxel import Voxel
from Cell import Cell, TumorCell, HealthyCell
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from World import World
from Process import *


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

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

world = World(100, 10)

for i in range(10000):
    voxel1 = world.find_voxel(np.random.uniform(-100,100,3))
    voxel1.add_cell(Cell(0.01, cycle_hours=20, life_expectancy=48, color='red'))
    voxel1.add_cell(Cell(0.01, cycle_hours=20, life_expectancy=48, color='blue'))


world.show_voxels_centers(ax,fig)

plt.show()


##########################################################################################
end_time = 100
dt = 10

celldivision = CellDivision('cell_division', dt)
cellapoptosis = CellApoptosis('cell_apoptosis',dt)
cellaging = CellAging('cell_aging',dt)

list_of_processes = [celldivision, cellapoptosis, cellaging]


print('starting number of cells: ', len(world.voxel_list[0].list_of_cells))
#choose a few cells to print their age
for i in range(10):
    print('starting age of cells: ', world.voxel_list[0].list_of_cells[i].age)

sim = Simulator(list_of_processes,end_time,dt)
sim.run(world)

print('ending number of cells: ', len(world.voxel_list[0].list_of_cells))
for i in range(10):
    print('ending age of cells: ', world.voxel_list[0].list_of_cells[i].age)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
world.show_voxels_centers(ax2,fig2)
plt.show()
