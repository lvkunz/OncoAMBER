import sys
import amber

sys.path.insert(0, '/PHShome/lk001/.conda/envs/amberenv/lib/python3.9/site-packages') #cluster
amber.run_calibration(side=2, a=-7, b=1.0, max_n=10)
