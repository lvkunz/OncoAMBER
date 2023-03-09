import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Cell (object):
    def __init__(self, radius, cycle_hours = 10, life_expectancy = 1000, color = 'blue'):
        self.radius = radius
        self.state = 'cycling' # 'cycling', 'senescent', 'dead'
        self.doubling_time = cycle_hours
        self.volume = 4/3 * np.pi * self.radius**3
        self.color = color
        self.age = 0
        self.life_expectancy = life_expectancy

class TumorCell(Cell):
    def __init__(self, radius, color = 'red'):
        Cell.__init__(self, radius, color)

class HealthyCell (Cell):
    def __init__(self, radius, color = 'blue'):
        Cell.__init__(self, radius, color)