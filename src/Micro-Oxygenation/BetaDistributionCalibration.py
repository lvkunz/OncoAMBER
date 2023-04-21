import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
from scipy.stats import beta
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import seaborn as sns

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


def run_calibration(side = 6, a = -7, b = 1, max_n = 100):

    radius = 0.1

    def sigmoid(x, a=1, b=0.8):
        return 1 / (1 + np.exp(-a*(x-b)))

    def model(n, A, B, C, D):
        return A * np.exp(B * n) + C * n + D

    alpha_values = []
    beta_values = []
    all_n_values = []
    all_side_values = []
    all_pressure_values = []

    max = max_n + 1
    n_values = list(range(1, max, 1))

    for p in range(9):
        pressure = 0.1 * p
        pressure = round(pressure, 3)
        print('pressure', pressure)
        for n_idx, n in enumerate(n_values):
            print('n =', n)
            print('side =', side)
            sampler = qmc.Halton(2)
            multiple_alpha_values = []
            multiple_beta_values = []
            for _ in range(10):

                points_x = sampler.random(n)[:,0] * side #quasi random to have a distribution mimicking distance between vessels
                points_y = sampler.random(n)[:,1] * side
                # points_x = np.random.uniform(0, side, n) #pseudo random
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


                # Fit a beta distribution to the data
                alpha, beta_param, _, _ = beta.fit(o2_values, floc=0, fscale=1.0)
                multiple_alpha_values.append(alpha)
                multiple_beta_values.append(beta_param)

                # r = np.random.beta(alpha, beta_param, 5000)
                # rbis = np.random.beta(modelled_alpha, modelled_beta, 5000)

                # plt.figure()
                # plt.hist(o2_values, bins=100, alpha=0.5, label='Data', color='blue', density=True)
                # plt.hist(r, bins=100, alpha=0.5, label='Beta', color='red', density=True)
                # # plt.hist(rbis, bins=100, alpha=0.5, label='Modelled Beta', color='green', density=True)
                # # plt.plot(np.linspace(0, 1, 100), beta.pdf(np.linspace(0, 1, 100), alpha, beta_param))
                # plt.title('n = ' + str(n) + ', side = ' + str(side) + ', pressure = ' + str(pressure) + '\n fitted alpha = ' + str(round(alpha,4)) + ', fitted beta = ' + str(round(beta_param,4)))#+ ', modelled alpha ' + str(round(modelled_alpha,4)) + 'modelled beta ' + str(round(modelled_beta,4)), fontsize=8)
                # plt.xlabel('O2')
                # plt.ylabel('Frequency')
                # plt.legend()
                # plt.show()

            alpha_mean = np.mean(multiple_alpha_values)
            beta_param_mean = np.mean(multiple_beta_values)

            alpha_values = np.append(alpha_values, alpha_mean)
            beta_values = np.append(beta_values, beta_param_mean)
            all_side_values = np.append(all_side_values, side)
            all_n_values = np.append(all_n_values, n)
            all_pressure_values = np.append(all_pressure_values, pressure)

    pressure_column = np.unique(all_pressure_values)
    n_column = np.unique(all_n_values)

    alpha_values = alpha_values.reshape(len(pressure_column), len(n_column))
    beta_values = beta_values.reshape(len(pressure_column), len(n_column))

    alpha_dataframe = pd.DataFrame(alpha_values, index=pressure_column, columns=n_column)
    print(alpha_dataframe)

    beta_dataframe = pd.DataFrame(beta_values, index=pressure_column, columns=n_column)
    print(beta_dataframe)

    alpha_dataframe.to_csv('alpha_dataframe' + str(side) + '.csv')
    beta_dataframe.to_csv('beta_dataframe' + str(side) + '.csv')

    sns.heatmap(alpha_dataframe, cmap='viridis')
    plt.show()

    sns.heatmap(beta_dataframe, cmap='viridis')
    plt.show()



