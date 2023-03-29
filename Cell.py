import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import ReadAndWrite as rw

CONFIG = rw.read_config_file('CONFIG.txt')


class Cell (object):
    def __init__(self, radius, cycle_hours = 10, type = 'HealthyCell'):
        self.radius = radius
        self.necrotic = False
        self.doubling_time = cycle_hours
        self.volume = 4/3 * np.pi * self.radius**3
        self.oxygen = 0.0
        self.type = type

    def duplicate(self):
        return Cell(self.radius, self.doubling_time, self.type)

    def vitality(self):
        factor = CONFIG['o2_to_vitality_factor']
        vitality = self.oxygen * factor
        vitality = min(vitality, 1.0)
        return vitality #needs to be normalized between 0 and 1




