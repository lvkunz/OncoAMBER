import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
from scipy.stats import beta
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D

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

        return distance

side = 6 #um/100
radius = 0.1
n_values = list(range(1, 100))

def sigmoid(x, a=1, b=0.8):
    return 1 / (1 + np.exp(-a*(x-b)))

def model(n, A, B, C, D):
    return A * np.exp(B * n) + C * n + D


a = -7 #um/100 give hypoxia threshold above 70 to 150um
b = 1 #um/100

alpha_values = []
beta_values = []
all_n_values = []
all_side_values = []

for _ in range(1):
    pressure = 0
    print('Iteration', _+1)
    for n_idx, n in enumerate(n_values):

        alpha_data = np.genfromtxt('alpha.csv', delimiter=',', skip_header=True)
        beta_data = np.genfromtxt('beta.csv', delimiter=',', skip_header=True)
        a_row_index = np.where(alpha_data[:, 0] == side)[0][0]
        b_row_index = np.where(beta_data[:, 0] == side)[0][0]

        aA = alpha_data[a_row_index, 1]
        aB = alpha_data[a_row_index, 2]
        aC = alpha_data[a_row_index, 3]
        aD = alpha_data[a_row_index, 4]
        bA = beta_data[b_row_index, 1]
        bB = beta_data[b_row_index, 2]
        bC = beta_data[b_row_index, 3]
        bD = beta_data[b_row_index, 4]

        modelled_alpha = model(n, 5.07964613e+01, -1.76688386e-02, 8.65381485e-01, -5.06390318e+01)
        modelled_beta = model(n, 0.50852227, -0.59366832, 0.00410601, 0.37548704)

        # modelled_alpha = model(n, aA, aB, aC, aD)
        # modelled_beta = model(n, bA, bB, bC, bD)



        print('n =', n)
        print('side =', side)
        sampler = qmc.Halton(2)
        points_x = sampler.random(n)[:,0] * side
        points_y = sampler.random(n)[:,1] * side
        # points_x = np.random.uniform(0, side, n)
        # points_y = np.random.uniform(0, side, n)


        vessels = []
        for i in range(len(points_x)):
            vessels.append(Vessel([points_x[i], points_y[i]], radius))

        points = []
        for i in range(5000):
            point = [np.random.uniform(0, side), np.random.uniform(0, side), np.random.uniform(0, side)]
            distances = []
            for vessel in vessels:
                distances.append(vessel.closest_distance(point))
            points.append(min(distances))

        o2_values = []
        for point in points:
            b_ = b*(1 - pressure)
            o2_values.append(sigmoid(point, a=a, b=b_))

        # Normalize the histogram values
        #hist_values_normalized = hist_values / len(o2_values)

        # Fit a beta distribution to the data
        alpha, beta_param, _, _ = beta.fit(o2_values, floc=0, fscale=1.0)

        r = np.random.beta(alpha, beta_param, 5000)
        rbis = np.random.beta(modelled_alpha, modelled_beta, 5000)


        plt.figure()
        plt.hist(o2_values, bins=100, alpha=0.5, label='Data', color='blue')
        plt.hist(r, bins=100, alpha=0.5, label='Beta', color='red')
        plt.hist(rbis, bins=100, alpha=0.5, label='Modelled Beta', color='green')
        #plt.plot(np.linspace(0, 1, 100), beta.pdf(np.linspace(0, 1, 100), alpha, beta_param))
        plt.title('n = ' + str(n) + ', side = ' + str(side) + ', pressure = ' + str(pressure) + '\n fitted alpha = ' + str(round(alpha,4)) + ', fitted beta = ' + str(round(beta_param,4)) + ', modelled alpha ' + str(round(modelled_alpha,4)) + 'modelled beta ' + str(round(modelled_beta,4)), fontsize=8)
        plt.xlabel('O2')
        plt.ylabel('Frequency')
        plt.legend()
        plt.show()

        alpha_values = np.append(alpha_values, alpha)
        beta_values = np.append(beta_values, beta_param)
        all_side_values = np.append(all_side_values, side)
        all_n_values = np.append(all_n_values, n)

#save all the alpha and beta values
print(alpha_values)
np.save('alpha_values.npy', alpha_values)
np.save('beta_values.npy', beta_values)
np.save('all_n_values.npy', all_n_values)
np.save('all_side_values.npy', all_side_values)

#read in the alpha and beta values




