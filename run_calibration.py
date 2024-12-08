import sys
import amber
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, '/PHShome/lk001/.conda/envs/amberenv/lib/python3.9/site-packages') #cluster
#amber.run_calibration(side=8, a=-7, b=1.0, max_n=50)

o2_big = amber.run_calibration_single(side=6, a=-7, b=1.0, n=27, p = 0.3)

o2_1 = amber.run_calibration_single(side=4, a=-7, b=1.0, n=12, p = 0.3)
o2_2 = amber.run_calibration_single(side=4, a=-7, b=1.0, n=12, p = 0.3)
o2_3 = amber.run_calibration_single(side=4, a=-7, b=1.0, n=12, p = 0.3)
o2_4 = amber.run_calibration_single(side=4, a=-7, b=1.0, n=12, p = 0.3)
o2_5 = amber.run_calibration_single(side=4, a=-7, b=1.0, n=12, p = 0.3)
o2_6 = amber.run_calibration_single(side=4, a=-7, b=1.0, n=12, p = 0.3)
o2_7 = amber.run_calibration_single(side=4, a=-7, b=1.0, n=12, p = 0.3)
o2_8 = amber.run_calibration_single(side=4, a=-7, b=1.0, n=12, p = 0.3)

o2_small = [o2_1, o2_2, o2_3, o2_4, o2_5, o2_6, o2_7, o2_8]

o2_small_flat = [item for sublist in o2_small for item in sublist]


#histogram of o2 levels


plt.hist(o2_big, bins=20, alpha=0.5, label='big', density=True)
plt.hist(o2_small_flat, bins=20, alpha=0.5, label='small', density=True)
plt.legend(loc='upper right')
plt.show()
