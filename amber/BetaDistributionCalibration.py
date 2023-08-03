import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta
import pandas as pd

#file for running subsimulations of oxygen distribution. Not executed in the main script.
#based on Espinoza 2013. doi: 10.1118/1.4812431
class Vessel_:
    def __init__(self, origin, radius, length):
        self.origin = (origin[0], origin[1], 0)
        self.end = (origin[0], origin[1], length)
        self.radius = radius

    def closest_distance(self, point): # computes closest distance between a point and a vessel

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


def run_calibration(side = 6, a = -7, b = 1, max_n = 100): #side is the length of the cube, a and b are the parameters of the sigmoid function, max_n is the maximum number of capillaries to be simulated

    radius = 0.1 #radius of capillaries, useless. In place if needed in the future
    R0 = b       #oxygen diffusion limit for square function

    def sigmoid(x, a=1, b=0.8): #sigmoid function
        return 1 / (1 + np.exp(-a*(x-b)))

    def square(x, R = 1.0): #quadratic function
        if x > R:
            return 0
        else:
            return (x/R - 1)**2

    all_n_values = []
    all_side_values = []
    all_pressure_values = []

    max = max_n + 1
    n_values = list(range(1, max, 1))

    for p in range(8): #loop over different pressures
        pressure = 0.1 * p  # pressure is a number between 0 and 0.7 (0.7 is the max crowding for hard spheres piling)
                            # pressure reduces oxygen diffusion limit
        pressure = round(pressure, 3)
        print('pressure', pressure)
        for n_idx, n in enumerate(n_values): #loop over different number of capillaries in the cube
            print('n =', n)
            print('side =', side)
            all_o2_values = []
            for _ in range(50): #number of iterations for each n and pressure (50 is enough to have a good estimate of the mean)
                print('iteration', _)

                points_x = np.random.uniform(0, side, n) #random. Best according to Espinoza 2013. doi: 10.1118/1.4812431
                points_y = np.random.uniform(0, side, n)

                vessels = []
                for i in range(len(points_x)):
                    vessels.append(Vessel_([points_x[i], points_y[i]], radius, side))

                points = []
                for i in range(1000):
                    point = [np.random.uniform(0, side), np.random.uniform(0, side), np.random.uniform(0, side)]
                    distances = []
                    for vessel in vessels:
                        distances.append(vessel.closest_distance(point))
                    points.append(min(distances))

                # o2_values = []
                for point in points: #TODO: if function needs to be changed by sigmoid or other. Do it here.
                    # b_ = b*(1 - pressure)
                    R = R0*(1 - pressure)
                    o2 = square(point, R = R)
                    # o2_values.append(o2)
                    all_o2_values.append(o2)



            #save and show outputs from here to end
            all_o2_values = np.array(all_o2_values)
            alpha_all, beta_param_all, _, _ = beta.fit(all_o2_values)

            r_all = beta.rvs(alpha_all, beta_param_all, size = 50000)
            plt.figure()
            plt.hist(all_o2_values, bins=100, alpha=0.5, label='Data', color='green', density=True)
            plt.hist(r_all, bins=100, alpha=0.5, label='Beta_all', color='red', density=True)
            plt.title('n = ' + str(n) + ', side = ' + str(side) + ', pressure = ' + str(pressure)+ '\n fitted alpha = ' + str(round(alpha_all,4)) + ', fitted beta = ' + str(round(beta_param_all,4)))
            plt.xlabel('O2')
            plt.ylabel('Frequency')
            plt.legend()
            plt.show()

            print('alpha all', alpha_all)
            print('beta all', beta_param_all)

            alpha_values = np.append(alpha_values, alpha_all)
            beta_values = np.append(beta_values, beta_param_all)
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



