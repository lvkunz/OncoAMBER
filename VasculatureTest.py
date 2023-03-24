
from Vesselv2 import *

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import math
def random_vegf_gradient(point):
    return np.array([random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)])

# Define random pressure function
def random_pressure(point):
    return random.uniform(0, 1)

first_vessel = Vessel(([1.0,-1.0,0.0],[0,0,0]), 1)

first_vessel.grow(
    vegf_gradient=random_vegf_gradient,
    pressure=random_pressure,
    steps=100,
    weight_direction=3.0,
    weight_vegf=0,
    weight_pressure=0,
    pressure_threshold=0.5
)

size = 10

fig = plt.figure()
fig.set_size_inches(10, 10)
fig.set_dpi(200)
ax = fig.add_subplot(111, projection='3d')
# ax.set_xlim(-size, size)
# ax.set_ylim(-size, size)
# ax.set_zlim(-size, size)
first_vessel.plot(ax)
plt.show()

branching_point = first_vessel.path[50]

network = VasculatureNetwork([first_vessel])

network.grow_and_split(splitting_rate= 0.07,
                       vegf_gradient=random_vegf_gradient,
                       pressure=random_pressure,
                       macro_steps=10,
                       weight_direction=3.0,
                       weight_vegf=0,
                       weight_pressure=0.0,
                       pressure_threshold=0.5
                       )

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

size = 10

fig2 = plt.figure()
fig2.set_size_inches(10, 10)
fig2.set_dpi(200)
ax2 = fig2.add_subplot(111, projection='3d')
# ax2.set_xlim(-size, size)
# ax2.set_ylim(-size, size)
# ax2.set_zlim(-size, size)
network.plot(ax2)
plt.show()