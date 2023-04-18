import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import ReadAndWrite as rw

CONFIG = rw.read_config_file('CONFIG.txt')
seed = CONFIG['seed']
if seed == -1:
    seed = np.random.randint(0, 1000000)
np.random.seed(seed)
print('seed: ', seed)

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
        factor = CONFIG['o2_to_vitality_factor']
        vitality = self.oxygen * factor
        vitality = min(vitality, 1.0)
        return vitality #needs to be normalized between 0 and 1

    def radiosensitivity(self):
        return CONFIG['radiosensitivity']*self.oxygen




