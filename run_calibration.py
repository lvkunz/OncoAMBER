import sys
import amber

sys.path.insert(0, '/PHShome/lk001/.conda/envs/amberenv/lib/python3.9/site-packages') #cluster
amber.run_calibration(side=6, a=-7, b=1.5, max_n=200)

