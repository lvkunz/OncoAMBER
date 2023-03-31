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



#plot surface, alpha, and beta values vs. n and side

def model_alpha(n, A, B, C, D):
    return A * np.exp(B * n) + C * n + D

def model_beta(n, A, B, C, D):
    return A * np.exp(-B * n) + C * n + D
all_n_values = np.array(all_n_values)
all_side_values = np.array(all_side_values)
alpha_values = np.array(alpha_values)
beta_values = np.array(beta_values)

# all_n_values = all_n_values[:300]
# all_side_values = all_side_values[:300]
# alpha_values = alpha_values[:300]
# beta_values = beta_values[:300]

#compute average values for each n
new_alpha_values = []
new_beta_values = []

max = 50
for i in range(1, max+1):
    new_alpha_values.append(np.mean(alpha_values[all_n_values == i]))
    new_beta_values.append(np.mean(beta_values[all_n_values == i]))

alpha_values = np.array(new_alpha_values)
beta_values = np.array(new_beta_values)
all_n_values = np.array(range(1, max+1))



# plt.figure(figsize=(8, 6))
# ax = plt.axes(projection='3d')
# ax.scatter(all_n_values, alpha_values)
# ax.set_xlabel('n')
# ax.set_title('Alpha values depending on n and side')
# plt.show()
#
# plt.figure(figsize=(8, 6))
# ax = plt.axes(projection='3d')
# ax.scatter(all_n_values, beta_values)
# ax.set_xlabel('n')
# ax.set_title('Beta values depending on n and side')
# plt.show()


# # Fit the alpha function
alpha_params, _ = curve_fit(model_alpha, all_n_values, alpha_values,[1.0, 0.03, -5, 0.1], maxfev=100000)
# w
# # Fit the beta function
beta_params, _ = curve_fit(model_beta, all_n_values, beta_values,[0.5,0.05,-2.0,0.5], maxfev=100000)
#
# # Print the estimated parameters
print("Alpha parameters: A, B, C, D =", alpha_params)
print("Beta parameters: E, F, G, H =", beta_params)


#alpha_params = [4.73539119e+06, 1.44161645e-04, -6.83866915e+02, -3.00000000e+01, -4.73538663e+06]
alpha_fitted = model_alpha(all_n_values, *alpha_params)
beta_fitted = model_beta(all_n_values, *beta_params)

# Plot the original alpha and beta values and the fits
plt.figure(figsize=(8, 6))
plt.scatter(all_n_values, beta_values, label='Beta (data)', color='C1')
plt.plot(all_n_values, beta_fitted, label='Beta (fit)', linestyle='--', color='C1')
plt.xlabel('n')
plt.ylabel('Values')
plt.title('Beta values and fits depending on n')
plt.legend()
plt.savefig('beta_values.png')
plt.show()

plt.figure(figsize=(8, 6))
plt.scatter(all_n_values, alpha_values, label='Alpha (data)', color='C0')
plt.plot(all_n_values, alpha_fitted, label='Alpha (fit)', linestyle='--', color='C0')
plt.xlabel('n')
plt.ylabel('Values')
plt.title('Alpha values and fits depending on n')
plt.legend()
plt.savefig('alpha_values.png')
plt.show()
