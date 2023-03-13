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
        while self.time <= self.finish_time:
            print('Time: {} hours'.format(self.time))
            process_local = []
            process_global = []
            for process in self.list_of_process:
                if process.is_global:
                    process_global.append(process)
                else:
                    process_local.append(process)
            #loop in random order
            copy_voxel_list = world.voxel_list.copy()
            np.random.shuffle(copy_voxel_list)

            for voxel in copy_voxel_list:
                for process in process_local:
                        process(voxel)
            for process in process_global:
                process(world)
            self.time = self.time + self.dt


class Process:
    def __init__(self, name, dt):
        self.name = name
        self.dt = dt
        self.is_global = False
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
            effective_age = cell.age + self.dt/2
            probability = 1 - np.exp(-effective_age / cell.life_expectancy)
            #print(probability)
            if np.random.random() < probability:
                voxel.remove_cell(cell)
class CellAging(Process):
    def __init__(self, name, dt):
        super().__init__('CellAging', dt)
    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            cell.age = cell.age + self.dt

class CellMigration(Process):
    def __init__(self, name, dt):
        super().__init__('CellMigration', dt)
        self.is_global = True
    def __call__(self, world : World):
        #print('voxel', voxel.voxel_number)
        exchange_matrix = world.compute_exchange_matrix(self.dt)
        for voxel in world.voxel_list:
            for cell in voxel.list_of_cells:
                list_of_neighbors = world.find_neighbors(voxel)
                #print('voisins', list_of_neighbors[0].voxel_number)
                np.random.shuffle(list_of_neighbors)
                for neighbor in list_of_neighbors:
                    probability = exchange_matrix[voxel.voxel_number, neighbor.voxel_number]
                    if np.random.random() < probability:
                        neighbor.add_cell(cell)
                        voxel.remove_cell(cell)
                        break