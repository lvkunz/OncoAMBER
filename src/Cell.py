import numpy as np
from src.config_instance import config

class Cell (object):
    def __init__(self, radius, cycle_hours = 10, type = 'NormalCell'):
        self.radius = radius
        self.necrotic = False
        self.usual_cycle_length = cycle_hours
        self.doubling_time = self.random_doubling_time(self.usual_cycle_length, config.doubling_time_sd) #new cells have a random doubling time
        self.volume = 4/3 * np.pi * self.radius**3
        self.oxygen = 1.0
        if type != 'NormalCell' and type != 'TumorCell':
            print(type + ' is not a valid cell type')
            raise ValueError('Cell type must be either NormalCell or TumorCell')
        self.type = type
        self.time_before_death = None
        self.time_spent_cycling = self.random_doubling_time(self.usual_cycle_length, config.doubling_time_sd) #new cells have already been cycling for a random amount of time
    def duplicate(self): #returns a new cell with the same properties
        cell = Cell(self.radius, self.usual_cycle_length, self.type)
        cell.time_spent_cycling = 0 #the new cell has not been cycling yet
        return cell

    def vitality(self):
        factor = config.o2_to_vitality_factor
        vitality = self.oxygen * factor
        vitality = min(vitality, 1.0)
        return vitality #needs to be normalized between 0 and 1

    def radiosensitivity(self):
        return config.radiosensitivity*self.oxygen

    def random_doubling_time(self, mu = None, sigma = config.doubling_time_sd):
        if mu is None:
            mu = self.usual_cycle_length
        #return a sample of a gamma distribution with mean cycle_hours
        shape = (mu / sigma) ** 2
        scale = sigma ** 2 / mu
        return np.random.gamma(shape = shape, scale = scale)

