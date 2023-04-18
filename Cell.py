import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import ReadAndWrite as rw
from config_instance import config

class Cell (object):
    def __init__(self, radius, cycle_hours = 10, type = 'NormalCell'):
        self.radius = radius
        self.necrotic = False
        self.doubling_time = cycle_hours
        self.volume = 4/3 * np.pi * self.radius**3
        self.oxygen = 1.0
        if type != 'NormalCell' and type != 'TumorCell':
            print(type + ' is not a valid cell type')
            raise ValueError('Cell type must be either NormalCell or TumorCell')
            return
        self.type = type
        self.time_before_death = None
    def duplicate(self):
        return Cell(self.radius, self.doubling_time, self.type)

    def vitality(self):
        factor = config.o2_to_vitality_factor
        vitality = self.oxygen * factor
        vitality = min(vitality, 1.0)
        return vitality #needs to be normalized between 0 and 1

    def radiosensitivity(self):
        return config.radiosensitivity*self.oxygen




