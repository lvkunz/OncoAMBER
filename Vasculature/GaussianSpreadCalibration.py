import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
from scipy.stats import beta


class Vessel:
    def __init__(self, origin, radius):
        self.origin = (origin[0], origin[1], 0)
        self.end = (origin[0], origin[1], side)
        self.radius = radius

    def closest_distance(self, point):

        # Compute the direction vector of the vessel segment
        vessel_dir = (self.end[0] - self.origin[0], self.end[1] - self.origin[1], self.end[2] - self.origin[2])

        # Compute the vector between the point and the vessel's origin
        point_vec = (point[0] - self.origin[0], point[1] - self.origin[1], point[2] - self.origin[2])

        # Compute the projection of the point vector onto the vessel direction vector
        projection = (point_vec[0] * vessel_dir[0] + point_vec[1] * vessel_dir[1] + point_vec[2] * vessel_dir[2]) / \
                     (vessel_dir[0] * vessel_dir[0] + vessel_dir[1] * vessel_dir[1] + vessel_dir[2] * vessel_dir[2])
        projection = max(0, min(1, projection))  # clamp the projection to [0, 1]

        # Compute the closest point on the vessel segment to the given point
        closest_point = (self.origin[0] + projection * vessel_dir[0],
                         self.origin[1] + projection * vessel_dir[1],
                         self.origin[2] + projection * vessel_dir[2])

        # Compute the distance between the closest point and the given point
        distance = np.sqrt((closest_point[0] - point[0]) ** 2 + (closest_point[1] - point[1]) ** 2 + (
                    closest_point[2] - point[2]) ** 2)

        return distance - self.radius

side = 6 #um/100
radius = 0.1
n_values = list(range(1, 21))

def sigmoid(x, a=1, b=0.8):
    return 1 / (1 + np.exp(-a*(x-b)))

a = -7 #um/100 give hypoxia threshold above 70 to 150um
b = 1 #um/100

# Set up the plot grid
fig, axes = plt.subplots(4, 5, figsize=(15, 12))
axes = axes.flatten()

for n_idx, n in enumerate(n_values):
    sampler = qmc.Halton(2)
    points_x = sampler.random(n)[:,0] * side
    points_y = sampler.random(n)[:,1] * side
    # points_x = np.random.uniform(0, side, n)
    # points_y = np.random.uniform(0, side, n)


    vessels = []
    for i in range(len(points_x)):
        vessels.append(Vessel([points_x[i], points_y[i]], radius))

    points = []
    for i in range(100000):
        point = [np.random.uniform(0, side), np.random.uniform(0, side), np.random.uniform(0, side)]
        distances = []
        for vessel in vessels:
            distances.append(vessel.closest_distance(point))
        points.append(min(distances))

    o2_values = []
    for point in points:
        o2_values.append(sigmoid(point, a=a, b=b))

    hist_values, bin_edges = np.histogram(o2_values, bins=100)

    # Normalize the histogram values
    #hist_values_normalized = hist_values / len(o2_values)

    # Fit a beta distribution to the data
    alpha, beta_param, _, _ = beta.fit(o2_values, floc=0, fscale=1)

    # Plot the histogram
    axes[n_idx].hist(o2_values, bins=100, density=True)

    # Plot the fitted beta distribution
    x = np.linspace(0, 1, 1000)
    y = beta.pdf(x, alpha, beta_param)
    axes[n_idx].plot(x, y, 'r-', lw=2)

    # Set title, labels, and ylim
    axes[n_idx].set_title(f"n = {n}, alpha = {alpha:.2f}, beta = {beta_param:.2f}")
    axes[n_idx].set_xlabel("O2")
    axes[n_idx].set_ylabel("Frequency")
    # axes[n_idx].set_ylim(0, 1.0)

# Adjust the layout and display the plots
plt.tight_layout()
plt.show()




