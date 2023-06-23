import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

data2 = pd.read_csv('LLC_sc_CCSB.csv', sep=',', header=0)

data2['Time'] = data2['Time'] * 24

# Scatter plot of 'time' vs 'ctrl'
def func_exp (x, a, b):
    x = np.array(x)
    return a * np.exp(b * x) - a


time = []
volume = []
sd = []
for i in data2['Time'].unique():
    time.append(i)
    vol = []
    for j in data2[data2['Time'] == i]['Vol']:
        vol.append(j)
    volume.append(np.mean(vol))
    sd.append(np.std(vol))
order = np.argsort(time)
time = np.array(time)[order]
volume = np.array(volume)[order]

popt, pcov = curve_fit(func_exp, time, volume, p0=(1, 1e-6))

print(popt)
#doubling time
print(np.log(2)/popt[1])

fig, axes = plt.subplots(1, 1, figsize=(10, 5))
axes.scatter(time, volume, color='black', label='data')
axes.plot(time, func_exp(time, *popt), color='red', label='fit')
axes.set_xlabel('Time (days)')
axes.set_ylabel('Volume (mm^3)')
axes.legend()
plt.show()
