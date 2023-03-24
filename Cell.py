import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Cell (object):
    def __init__(self, radius, cycle_hours = 10):
        self.radius = radius
        self.necrotic = False
        self.doubling_time = cycle_hours
        self.volume = 4/3 * np.pi * self.radius**3
        self.oxygen = 0.0

    def duplicate(self):
        return Cell(self.radius, self.doubling_time)

    def vitality(self):
        factor = 1.0
        vitality = self.oxygen * factor
        vitality = min(vitality, 1.0)
        return vitality #needs to be normalized between 0 and 1

class TumorCell(Cell):
    def __init__(self, radius, cycle_hours = 5, life_expectancy = 1000, color = 'my purple'):
        Cell.__init__(self, radius, cycle_hours)

class HealthyCell (Cell):
    def __init__(self, radius, cycle_hours=10, life_expectancy=1000, color='my green'):
        Cell.__init__(self, radius, cycle_hours)

