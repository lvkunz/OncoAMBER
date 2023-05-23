import numpy as np
import scipy.stats as stats
from scipy.stats import gamma
class Cell (object):
    def __init__(self, radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, type = 'NormalCell'):
        self.radius = radius
        self.necrotic = False
        self.usual_cycle_length = cycle_hours
        self.cycle_length_std = cycle_std
        self.gamma_shape = (self.usual_cycle_length / self.cycle_length_std) ** 2
        self.gamma_scale = self.cycle_length_std ** 2 / self.usual_cycle_length
        self.doubling_time = self.random_doubling_time() #new cells have a random doubling time
        self.volume = 4/3 * np.pi * self.radius**3
        self.oxygen = 1.0
        self.intra_radiosensitivity = intra_radiosensitivity
        self.o2_to_vitality_factor = o2_to_vitality_factor
        if type != 'NormalCell' and type != 'TumorCell':
            print(type + ' is not a valid cell type')
            raise ValueError('Cell type must be either NormalCell or TumorCell')
        self.type = type
        self.time_before_death = None
        self.time_spent_cycling = self.random_time_spent_cycling() #new cells have already been cycling for a random amount of time
    def duplicate(self): #returns a new cell with the same properties
        cell = Cell(self.radius, self.usual_cycle_length, self.cycle_length_std, self.intra_radiosensitivity, self.o2_to_vitality_factor, self.type)
        cell.time_spent_cycling = 0 #the new cell has not been cycling yet
        return cell

    def vitality(self):
        factor = self.o2_to_vitality_factor
        vitality = self.oxygen * factor
        vitality = min(vitality, 1.0)
        return vitality #needs to be normalized between 0 and 1

    def radiosensitivity(self):
        print('#'*80)
        print(self.intra_radiosensitivity, self.oxygen)
        value = self.intra_radiosensitivity * self.oxygen
        if self.necrotic:
            value = 1
        return value

    def random_doubling_time(self):
        return np.random.gamma(shape = self.gamma_shape, scale = self.gamma_scale)

    def random_time_spent_cycling(self):
        longest_possible_time = self.doubling_time
        return np.random.uniform(0, longest_possible_time)
