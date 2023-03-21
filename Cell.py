import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Cell (object):
    def __init__(self, radius, cycle_hours = 10, life_expectancy = 1000, color = 'blue'):
        self.radius = radius
        self.vitality = 1.0 #below certain threshold cell enter senescence and eventually appoptosis
        self.necrosis_score = 0 #above certain threshold cell enter necrosis
        self.necrotic = False
        self.doubling_time = cycle_hours
        self.volume = 4/3 * np.pi * self.radius**3
        self.color = color
        self.age = 0
        self.life_expectancy = life_expectancy
        self.diffusion_coefficient = 1e-3/self.radius   # probably not the right value

    def duplicate(self):
        return Cell(self.radius, self.doubling_time, self.life_expectancy, self.color)

class TumorCell(Cell):
    def __init__(self, radius, cycle_hours = 5, life_expectancy = 1000, color = 'my purple'):
        Cell.__init__(self, radius, cycle_hours, life_expectancy, color)

class HealthyCell (Cell):
    def __init__(self, radius, cycle_hours=10, life_expectancy=1000, color='my green'):
        Cell.__init__(self, radius, cycle_hours, life_expectancy, color)

