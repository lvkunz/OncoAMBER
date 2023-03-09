from World import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Cell import *
from Voxel import *

class Simulator:
    def __init__(self, list_of_process : list, finish_time, dt):
        self.list_of_process = list_of_process
        self.finish_time = finish_time
        self.dt = dt
        self.time = 0

    def run(self, world: World):
        print('Running simulation for {} hours'.format(self.finish_time))
        exchange_matrix = world.compute_exchange_matrix()
        while self.time < self.finish_time:
            print('Time: {} hours'.format(self.time))
            for voxel in world.voxel_list:
                for process in self.list_of_process:
                    process(voxel)
            self.time = self.time + self.dt


class Process:
    def __init__(self, name, dt):
        self.name = name
        self.dt = dt
    def __call__(self, voxel):
        pass

class CellDivision(Process):
    def __init__(self, name, dt):
        super().__init__('CellDivision', dt)
    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            if cell.state == 'cycling':
                probability = 1 - np.exp(-self.dt/cell.doubling_time) #this is wrong
                if np.random.random() < probability:
                    voxel.add_cell(cell)
                    cell.age = 0 #the cell is still old but the new one is aged 0
class CellApoptosis(Process):
    def __init__(self, name, dt):
        super().__init__('CellDeath', dt)

    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            probability = 1 - np.exp(-cell.age/cell.life_expectancy) #this is wrong
            #print(probability)
            if np.random.random() < probability:
                voxel.remove_cell(cell)
class CellAging(Process):
    def __init__(self, name, dt):
        super().__init__('CellAging', dt)
    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            cell.age = cell.age + 1

class CellMigration(Process):
    def __init__(self, name, dt):
        super().__init__('CellMigration', dt)
    def __call__(self, voxel, exchange_matrix):
        for cell in voxel.list_of_cells:
            list_of_neighbors = find_neighbors(voxel)
            scrambled_list_of_neighbors = np.random.shuffle(list_of_neighbors)
            for neighbor in scrambled_list_of_neighbors:
                probability = exchange_matrix[voxel.voxel_number, neighbor.voxel_number]
                if np.random.random() < probability:
                    neighbor.add_cell(cell)
                    voxel.remove_cell(cell)
                    break