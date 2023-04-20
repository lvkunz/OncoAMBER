import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
from scipy.stats import beta
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

alpha_values = np.load('../Vasculature/alpha_values.npy')
beta_values = np.load('../Vasculature/beta_values.npy')
all_n_values = np.load('../Vasculature/all_n_values.npy')
all_side_values = np.load('../Vasculature/all_side_values.npy')
all_pressure_values = np.load('../Vasculature/all_pressure_values.npy')


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

aAt = []
aBt = []
aCt = []
aDt = []

bAt = []
bBt = []
bCt = []
bDt = []

aAs = []
aBs = []
aCs = []
aDs = []

bAs = []
bBs = []
bCs = []
bDs = []

aAb = []
aBb = []
aCb = []
aDb = []

bAb = []
bBb = []
bCb = []
bDb = []

corresponding_pressure = []

#select for pressure value
for i in range(16):
    pressure = 0.05*i
    pressure_array = np.array(all_pressure_values)
    tolerance = 1e-8 # Set an appropriate tolerance value based on your specific requirements

    filtered_indices = np.where(np.isclose(pressure_array, pressure, atol=tolerance))
    print('pressure = ', pressure, 'filtered_indices = ', filtered_indices)


    all_n_values_ = np.array(all_n_values)[filtered_indices]
    all_side_values_ = np.array(all_side_values)[filtered_indices]
    all_pressure_values_ = np.array(all_pressure_values)[filtered_indices]
    alpha_values_ = np.array(alpha_values)[filtered_indices]
    beta_values_ = np.array(beta_values)[filtered_indices]

    filtered_indices_bis = np.where(all_n_values_ < 50)

    all_n_values_small = np.array(all_n_values_)[filtered_indices_bis]
    all_side_values_small = np.array(all_side_values_)[filtered_indices_bis]
    all_pressure_values_small = np.array(all_pressure_values_)[filtered_indices_bis]
    alpha_values_small = np.array(alpha_values_)[filtered_indices_bis]
    beta_values_small = np.array(beta_values_)[filtered_indices_bis]

    filtered_indices_bisbis = np.where(all_n_values_ >= 20)

    all_n_values_big = np.array(all_n_values_)[filtered_indices_bisbis]
    all_side_values_big = np.array(all_side_values_)[filtered_indices_bisbis]
    all_pressure_values_big = np.array(all_pressure_values_)[filtered_indices_bisbis]
    alpha_values_big = np.array(alpha_values_)[filtered_indices_bisbis]
    beta_values_big = np.array(beta_values_)[filtered_indices_bisbis]

    filtered_indices_bisbisbis = np.where(all_n_values_ < 15)

    all_n_values_tiny = np.array(all_n_values_)[filtered_indices_bisbisbis]
    all_side_values_tiny = np.array(all_side_values_)[filtered_indices_bisbisbis]
    all_pressure_values_tiny = np.array(all_pressure_values_)[filtered_indices_bisbisbis]
    alpha_values_tiny = np.array(alpha_values_)[filtered_indices_bisbisbis]
    beta_values_tiny = np.array(beta_values_)[filtered_indices_bisbisbis]



    print(pressure)
    print('###################################################')




    # # Fit the alpha function
    alpha_params_small, _ = curve_fit(model, all_n_values_small, alpha_values_small,[1.0, 0.03, 1.0, 0.1], maxfev=100000)#, bounds=[(-np.inf, -np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)])
    # w
    # # Fit the beta function
    beta_params_small, _ = curve_fit(model, all_n_values_small, beta_values_small,[0.5,-0.05,-2.0,0.5], maxfev=100000)#, bounds=[(-np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)])

    alpha_params_big, _ = curve_fit(model, all_n_values_big, alpha_values_big,[1.0, 0.03, 1.0, 0.1], maxfev=100000)#, bounds=[(-np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)])

    beta_params_big, _ = curve_fit(model, all_n_values_big, beta_values_big,[0.5,-0.05,-2.0,0.5], maxfev=100000)#, bounds=[(-np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)])

    alpha_params_tiny, _ = curve_fit(model, all_n_values_tiny, alpha_values_tiny,[1.0, 0.03, 1.0, 0.1], maxfev=100000)#, bounds=[(-np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)])

    beta_params_tiny, _ = curve_fit(model, all_n_values_tiny, beta_values_tiny,[0.5,-0.05,-2.0,0.5], maxfev=100000)#, bounds=[(-np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)])

    #
    # # Print the estimated parameters
    print("Alpha parameters small: A, B, C, D =", alpha_params_small)
    print("Beta parameters small: E, F, G, H =", beta_params_small)
    print("Alpha parameters big: A, B, C, D =", alpha_params_big)
    print("Beta parameters big: E, F, G, H =", beta_params_big)

    aAs.append(alpha_params_small[0])
    aBs.append(alpha_params_small[1])
    aCs.append(alpha_params_small[2])
    aDs.append(alpha_params_small[3])

    bAs.append(beta_params_small[0])
    bBs.append(beta_params_small[1])
    bCs.append(beta_params_small[2])
    bDs.append(beta_params_small[3])

    aAb.append(alpha_params_big[0])
    aBb.append(alpha_params_big[1])
    aCb.append(alpha_params_big[2])
    aDb.append(alpha_params_big[3])

    bAb.append(beta_params_big[0])
    bBb.append(beta_params_big[1])
    bCb.append(beta_params_big[2])
    bDb.append(beta_params_big[3])

    aAt.append(alpha_params_tiny[0])
    aBt.append(alpha_params_tiny[1])
    aCt.append(alpha_params_tiny[2])
    aDt.append(alpha_params_tiny[3])

    bAt.append(beta_params_tiny[0])
    bBt.append(beta_params_tiny[1])
    bCt.append(beta_params_tiny[2])
    bDt.append(beta_params_tiny[3])


    corresponding_pressure.append(pressure)


    #alpha_params = [4.73539119e+06, 1.44161645e-04, -6.83866915e+02, -3.00000000e+01, -4.73538663e+06]
    alpha_fitted_small = model(all_n_values_small, *alpha_params_small)
    beta_fitted_small = model(all_n_values_small, *beta_params_small)
    alpha_fitted_big = model(all_n_values_big, *alpha_params_big)
    beta_fitted_big = model(all_n_values_big, *beta_params_big)
    alpha_fitted_tiny = model(all_n_values_tiny, *alpha_params_tiny)
    beta_fitted_tiny = model(all_n_values_tiny, *beta_params_tiny)


    # alpha_modeled = model(all_n_values_small, *[2.41862595e+02, 1.13938577e-02,-3.00000000e+00,-2.42338686e+02])
    # beta_modeled = model(all_n_values_small, *[0.58830006, -0.05656126, 0.01705466, -0.07116887])
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
    ax.scatter(all_n_values_, alpha_values_, c='r', marker='o', alpha=0.5)
    ax.plot(all_n_values_big, alpha_fitted_big, c='g')
    ax.plot(all_n_values_small, alpha_fitted_small, c='b')
    ax.plot(all_n_values_tiny, alpha_fitted_tiny, c='y')
    ax.set_xlabel('n')
    ax.set_ylabel('alpha, pressure = ' + str(pressure))
    plt.show()



    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(all_n_values_, alpha_values_, c='r', marker='o', alpha = 0.5)
    ax.plot(all_n_values_big, alpha_fitted_big, c='g')
    ax.plot(all_n_values_small, alpha_fitted_small, c='b')
    ax.plot(all_n_values_tiny, alpha_fitted_tiny, c='y')
    ax.set_xlim(0, 20)
    ax.set_ylim(0.0, 1.5)
    ax.set_xlabel('n')
    ax.set_ylabel('alpha, pressure = ' + str(pressure))
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(all_n_values_, alpha_values_, c='r', marker='o', alpha=0.5)
    ax.plot(all_n_values_big, alpha_fitted_big, c='g')
    ax.plot(all_n_values_small, alpha_fitted_small, c='b')
    ax.plot(all_n_values_tiny, alpha_fitted_tiny, c='y')
    ax.set_xlim(10, 30)
    ax.set_ylim(0.0, 10)
    ax.set_xlabel('n')
    ax.set_ylabel('alpha, pressure = ' + str(pressure))
    plt.show()

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.scatter(all_n_values_, beta_values_, c='r', marker='o', alpha = 0.5)
    # ax.plot(all_n_values_big, beta_fitted_big, c='g')
    # ax.plot(all_n_values_small, beta_fitted_small, c='b')
    # ax.plot(all_n_values_tiny, beta_fitted_tiny, c='y')
    # # ax.set_xlim(0, 10)
    # # ax.set_ylim(0.0, 1)
    # ax.set_xlabel('n')
    # ax.set_ylabel('beta, pressure = ' + str(pressure))
    # plt.show()

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

data_frame_small = pd.DataFrame({'pressure': corresponding_pressure, 'aA': aAs, 'aB': aBs, 'aC': aCs, 'aD': aDs, 'bA': bAs, 'bB': bBs, 'bC': bCs, 'bD': bDs})
data_frame_big = pd.DataFrame({'pressure': corresponding_pressure, 'aA': aAb, 'aB': aBb, 'aC': aCb, 'aD': aDb, 'bA': bAb, 'bB': bBb, 'bC': bCb, 'bD': bDb})
data_frame_tiny = pd.DataFrame({'pressure': corresponding_pressure, 'aA': aAt, 'aB': aBt, 'aC': aCt, 'aD': aDt, 'bA': bAt, 'bB': bBt, 'bC': bCt, 'bD': bDt})

data_frame_small.to_csv('small.csv')
data_frame_big.to_csv('big.csv')
data_frame_tiny.to_csv('tiny.csv')


# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title('alphaSmall')
# ax.scatter(corresponding_pressure, aAs, c='r', marker='o', label = 'aA')
# ax.scatter(corresponding_pressure, aBs, c='g', marker='o', label = 'aB')
# ax.scatter(corresponding_pressure, aCs, c='b', marker='o', label = 'aC')
# ax.scatter(corresponding_pressure, aDs, c='y', marker='o', label = 'aD')
# plt.legend()
# ax.set_xlabel('pressure')
# plt.show()
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title('betaSmall')
# ax.scatter(corresponding_pressure, bAs, c='r', marker='o', label = 'bA')
# ax.scatter(corresponding_pressure, bBs, c='g', marker='o', label = 'bB')
# ax.scatter(corresponding_pressure, bCs, c='b', marker='o', label = 'bC')
# ax.scatter(corresponding_pressure, bDs, c='y', marker='o', label = 'bD')
# plt.legend()
# ax.set_xlabel('pressure')
# plt.show()
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title('alphaBig')
# ax.scatter(corresponding_pressure, aAb, c='r', marker='o', label = 'aA')
# ax.scatter(corresponding_pressure, aBb, c='g', marker='o', label = 'aB')
# ax.scatter(corresponding_pressure, aCb, c='b', marker='o', label = 'aC')
# ax.scatter(corresponding_pressure, aDb, c='y', marker='o', label = 'aD')
# plt.legend()
# ax.set_xlabel('pressure')
# plt.show()
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title('betaBig')
# ax.scatter(corresponding_pressure, bAb, c='r', marker='o', label = 'bA')
# ax.scatter(corresponding_pressure, bBb, c='g', marker='o', label = 'bB')
# ax.scatter(corresponding_pressure, bCb, c='b', marker='o', label = 'bC')
# ax.scatter(corresponding_pressure, bDb, c='y', marker='o', label = 'bD')
# plt.legend()
# ax.set_xlabel('pressure')
# plt.show()





