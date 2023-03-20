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
from Vessel import *


world = World(15, 15)
vessel = Vessel([-1,0.5,0],[0.2,0.3,0.0],0.1)
points = vessel.generate_points_along_axis(10)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#dpi
fig.set_dpi(100)
limits = [-1, 1]
ax.set_xlim(limits)
ax.set_ylim(limits)
ax.set_zlim(limits)
vessel.plot(fig,ax)
ax.scatter(points[:,0], points[:,1], points[:,2], color='black')
plt.show()


vasculature = VasculatureNetwork([Vessel([0,0,0],[0,0,1],0.1)])
vasculature.grow_vasculature(Sphere(center = np.array([0,0,0]), radius = 1.0).generate_random_points(500))
# vasculature.save('Vasculature/vasculature_1000.txt')
#
# vasculature = VasculatureNetwork()
# vasculature.read('Vasculature/vasculature_1000.txt')

world.vasculature = vasculature
points = world.vasculature.list_of_vessels[0].generate_random_points_along_axis(10)
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# #dpi
# fig.set_dpi(100)
# limits = [-1, 1]
# ax.set_xlim(limits)
# ax.set_ylim(limits)
# ax.set_zlim(limits)
# vasculature.plot(fig, ax)
# ax.scatter(points[:,0], points[:,1], points[:,2], color='black')
# plt.savefig('Vasculature/vasculature.png')
# plt.show()
#
#
#
#
# # kernel_limits_x = [-1, 1]
# # kernel_limits_y = [-1, 1]
# # kernel_limits_z = [-1, 1]
# # kernel_resolution = 0.1
# # #create a matrix of zeros with the size of the kernel
# # kernel_coord_x = np.arange(kernel_limits_x[0], kernel_limits_y[1], kernel_resolution)
# # kernel_coord_y = np.arange(kernel_limits_y[0], kernel_limits_y[1], kernel_resolution)
# # kernel_coord_z = np.arange(kernel_limits_z[0], kernel_limits_z[1], kernel_resolution)
# # kernel_values = np.zeros((len(kernel_coord_x), len(kernel_coord_y), len(kernel_coord_z)))
# #
# # for i in range(len(kernel_coord_x)):
# #     for j in range(len(kernel_coord_y)):
# #         for k in range(len(kernel_coord_z)):
# #             #distance to the center of the kernel
# #             distance = np.sqrt(kernel_coord_x[i]**2 + kernel_coord_y[j]**2 + kernel_coord_z[k]**2)
# #             kernel_values[i][j][k] = 1
# #
# #
# # oxygen_kernel = Kernel(kernel_values, kernel_coord_x, kernel_coord_y, kernel_coord_z, )
# #
# # world.compute_oxygen_map_v2(oxygen_kernel,1)
# # map = world.oxygen_map
# # #heat map of the oxygen map
# # fig = plt.figure()
# # ax = fig.add_subplot(111, projection='3d')
# # #dpi
# # fig.set_dpi(100)
# # limits = [-1, 1]
# # ax.set_xlim(limits)
# # ax.set_ylim(limits)
# # ax.set_zlim(limits)
# # plt.imshow(map[:, :, 0], extent=[limits[0], limits[1], limits[0], limits[1]], cmap='hot')
# # plt.savefig('Vasculature/oxygen_map.png')
# # plt.show()
#
#
