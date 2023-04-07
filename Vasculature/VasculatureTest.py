
from Vessel import *

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import math
from BasicGeometries import *
from graphviz import Digraph
import networkx as nx
from ScalarField import *



seed = random.randint(0, 1000000)
random.seed(0)



#723384
def central_vegf_gradient(point):
    #vector that points to the center of the circle
    center = np.array([0,0,0])
    center_vector = center - point
    center_vector_norm = center_vector / np.linalg.norm(center_vector)
    return center_vector_norm

# Define random pressure function
def central_pressure(point):
    #increases towards the center
    center = np.array([0,0,0])
    center_vector = center - point
    center_vector_norm = center_vector / np.linalg.norm(center_vector)
    pressure = np.linalg.norm(center_vector_norm)
    #add a high pressure in the center
    if np.linalg.norm(point) < 1:
        pressure = 0
    return pressure



#vegf_gradient constant upwards
def upwards_vegf_gradient(point):
    if point[2] > 0:
        return np.array([0,0,1])
    else:
        return np.array([0,0,-1])
def upwards_pressure(point):
    if point[2] > 5:
        return (point[2])
    else:
        return (point[2])


#create a grid of points
points = []
values = []
n = 10
for i in range(n):
    for j in range(n):
        for k in range(n):
            i, j, k = i - n/2, j - n/2, k - n/2
            points.append([i,j,k])
            values.append(upwards_pressure([i,j,k]))

field = ScalarField3D(points, values, 1)


fig = plt.figure()
fig.set_size_inches(10, 10)
ax = fig.add_subplot(111, projection='3d')
ax.view_init(0, 0)
field.show_values(fig,ax, cmap='jet', vmin=None, vmax=None)
plt.show()

n = 30
Sphere = Sphere(20,[0,0,0])
#points = Sphere.generate_random_points_on_surface(10)
points_x = np.random.uniform(-20,20,n)
points_y = np.random.uniform(-20,20,n)
points_z = -20*np.ones(n)

points = []
for i in range(len(points_x)):
    points.append([points_x[i],points_y[i],points_z[i]])
#create a second point towards the center
print(points)
points2 = []
for i in range(len(points)):
    #points2.append(points[i] - points[i]/20)
    points2.append([points[i][0],points[i][1],points[i][2]+0.3])
print(points2)
vessels = []
for i in range(len(points)):
    vessels.append(Vessel([points[i],points2[i]],1))

network = VasculatureNetwork(vessels)

size = 20

fig = plt.figure()
fig.set_size_inches(10, 10)
fig.set_dpi(300)
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-size, size)
ax.set_ylim(-size, size)
ax.set_zlim(-size, size)
network.plot(fig, ax)
plt.show()



network.grow_and_split(splitting_rate= 0.02,
                       vegf_gradient=upwards_vegf_gradient,
                       pressure=upwards_pressure,
                       macro_steps=5,
                       micro_steps=10,
                       weight_direction=10.0,
                       weight_vegf=0.7,
                       weight_pressure=0.7,
                       pressure_threshold=0.5

                       )

network.update_vessels_radius_from_last(1.0)



# network.branching(first_vessel.id, branching_point)
#
# network.grow(
#     vegf_gradient=random_vegf_gradient,
#     pressure=random_pressure,
#     steps=100,
#     weight_direction=5.0,
#     weight_vegf=0,
#     weight_pressure=0,
#     pressure_threshold=0.5
# )

size = 20

fig2 = plt.figure()
fig2.set_size_inches(10, 10)
fig2.set_dpi(300)
ax2 = fig2.add_subplot(111, projection='3d')
#set view
#ax2.view_init(azim=0, elev=90)
ax2.set_xlim(-size, size)
ax2.set_ylim(-size, size)
ax2.set_zlim(-size, size)
network.plot(ax2)
plt.title("Vasculature Network, pressure weight = 0.7, vegf weight = 0.7")
plt.savefig("VasculatureNetwork7.png")
plt.show()
print("Seed: ", seed)

network.print_vessel_tree()

def plot_vessel_network(vessel_network):
    G = Digraph()

    # Add nodes to the graph
    for vessel in vessel_network.list_of_vessels:
        G.node(str(vessel.id), label=f"ID: {vessel.id}\nRadius: {vessel.radius:.2f}")

    # Add edges between parent and child nodes
    for vessel in vessel_network.list_of_vessels:
        for child_id in vessel.children_ids:
            G.edge(str(vessel.id), str(child_id))

    # Render and view the graph
    G.view()

def plot_vessel_network_lines(vessel_network):
    G = Digraph('VesselNetwork', format='png')
    G.attr(rankdir='TB', size='10,10')

    # Set default node attributes
    G.attr('node', shape='point', width='0', height='0')

    # Add nodes to the graph
    for vessel in vessel_network.list_of_vessels:
        G.node(str(vessel.id))

    # Add edges between parent and child nodes
    for vessel in vessel_network.list_of_vessels:
        for child_id in vessel.children_ids:
            G.edge(str(vessel.id), str(child_id), color='black', arrowsize='0')

    # Render and view the graph
    G.view()

#plot_vessel_network(network)

