import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
from scipy.stats import beta
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D

alpha_values = np.load('alpha_values.npy')
beta_values = np.load('beta_values.npy')
all_n_values = np.load('all_n_values.npy')
all_side_values = np.load('all_side_values.npy')
all_pressure_values = np.load('all_pressure_values.npy')


#plot surface, alpha, and beta values vs. n and side

def model(n, A, B, C, D):
    return A * np.exp(B * n) + C * n + D

all_n_values = np.array(all_n_values)
all_side_values = np.array(all_side_values)
all_pressure_values = np.array(all_pressure_values)
alpha_values = np.array(alpha_values)
beta_values = np.array(beta_values)


print(all_pressure_values)
#print all existing pressure values
print(np.unique(all_pressure_values))

#select for pressure value
for i in range(16):
    pressure = 0.05*i
    pressure_array = np.array(all_pressure_values)
    tolerance = 1e-5  # Set an appropriate tolerance value based on your specific requirements

    filtered_indices = np.where(np.isclose(pressure_array, pressure, atol=tolerance))


    all_n_values = np.array(all_n_values)[filtered_indices]
    all_side_values = np.array(all_side_values)[filtered_indices]
    all_pressure_values = np.array(all_pressure_values)[filtered_indices]
    alpha_values = np.array(alpha_values)[filtered_indices]
    beta_values = np.array(beta_values)[filtered_indices]

    filtered_indices_bis = np.where(all_n_values < 30)

    all_n_values_small = np.array(all_n_values)[filtered_indices_bis]
    all_side_values_small = np.array(all_side_values)[filtered_indices_bis]
    all_pressure_values_small = np.array(all_pressure_values)[filtered_indices_bis]
    alpha_values_small = np.array(alpha_values)[filtered_indices_bis]
    beta_values_small = np.array(beta_values)[filtered_indices_bis]

    filtered_indices_bisbis = np.where(all_n_values >= 30)

    all_n_values_big = np.array(all_n_values)[filtered_indices_bisbis]
    all_side_values_big = np.array(all_side_values)[filtered_indices_bisbis]
    all_pressure_values_big = np.array(all_pressure_values)[filtered_indices_bisbis]
    alpha_values_big = np.array(alpha_values)[filtered_indices_bisbis]
    beta_values_big = np.array(beta_values)[filtered_indices_bisbis]


    print(pressure)
    print('###################################################')




    # # Fit the alpha function
    alpha_params_small, _ = curve_fit(model, all_n_values_small, alpha_values_small,[1.0, 0.03, 1.0, 0.1], maxfev=100000, bounds=[(0.0, 0.0, -3.0, -np.inf), (np.inf, np.inf, np.inf, np.inf)])
    # w
    # # Fit the beta function
    beta_params_small, _ = curve_fit(model, all_n_values_small, beta_values_small,[0.5,-0.05,-2.0,0.5], maxfev=100000, bounds=[(-np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)])

    alpha_params_big, _ = curve_fit(model, all_n_values_big, alpha_values_big,[1.0, 0.03, 1.0, 0.1], maxfev=100000, bounds=[(0.0, 0.0, -3.0, -np.inf), (np.inf, np.inf, np.inf, np.inf)])

    beta_params_big, _ = curve_fit(model, all_n_values_big, beta_values_big,[0.5,-0.05,-2.0,0.5], maxfev=100000, bounds=[(-np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)])


    #
    # # Print the estimated parameters
    print("Alpha parameters small: A, B, C, D =", alpha_params_small)
    print("Beta parameters small: E, F, G, H =", beta_params_small)
    print("Alpha parameters big: A, B, C, D =", alpha_params_big)
    print("Beta parameters big: E, F, G, H =", beta_params_big)


    #alpha_params = [4.73539119e+06, 1.44161645e-04, -6.83866915e+02, -3.00000000e+01, -4.73538663e+06]
    alpha_fitted_small = model(all_n_values, *alpha_params_small)
    beta_fitted_small = model(all_n_values, *beta_params_small)
    alpha_fitted_big = model(all_n_values, *alpha_params_big)
    beta_fitted_big = model(all_n_values, *beta_params_big)

    alpha_modeled = model(all_n_values_small, *[2.41862595e+02, 1.13938577e-02,-3.00000000e+00,-2.42338686e+02])
    beta_modeled = model(all_n_values_small, *[0.58830006, -0.05656126, 0.01705466, -0.07116887])
    #plot alpha and beta values vs. n and pressure
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(all_n_values, all_pressure_values, alpha_values, c='r', marker='o')
    # ax.set_xlabel('n')
    # ax.set_ylabel('pressure')
    # ax.set_zlabel('alpha')
    # plt.show()
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(all_n_values, all_pressure_values, beta_values, c='r', marker='o')
    # ax.set_xlabel('n')
    # ax.set_ylabel('pressure')
    # ax.set_zlabel('beta')
    # plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(all_n_values, alpha_values, c='r', marker='o', alpha = 0.5)
    ax.plot(all_n_values, alpha_fitted_small, c='b')
    ax.plot(all_n_values, alpha_fitted_big, c='g')
    # ax.set_xlim(0, 10)
    # ax.set_ylim(-1, 1)
    ax.set_xlabel('n')
    ax.set_ylabel('alpha')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(all_n_values, beta_values, c='r', marker='o', alpha = 0.5)
    ax.plot(all_n_values, beta_fitted_small, c='b')
    ax.plot(all_n_values, beta_fitted_big, c='g')
    # ax.set_xlim(0, 20)
    ax.set_xlabel('n')
    ax.set_ylabel('beta')
    plt.show()

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.scatter(all_pressure_values, alpha_values, c='r', marker='o')
    # ax.set_xlabel('pressure')
    # ax.set_ylabel('alpha')
    # plt.show()
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.scatter(all_pressure_values, beta_values, c='r', marker='o')
    # ax.set_xlabel('pressure')
    # ax.set_ylabel('beta')
    # plt.show()



