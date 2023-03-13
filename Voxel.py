from Cell import Cell, TumorCell, HealthyCell
from BasicPlots import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random


class Voxel(object):
        def __init__(self, position = np.array([0,0,0]), half_length = 0.1, list_of_cells_in = [], oxygen = 0, voxel_number = 0):
                self.position = position
                self.half_length = half_length
                self.list_of_cells = list_of_cells_in
                self.oxygen = oxygen
                self.volume = 8*half_length**3
                self.free_space = self.volume - sum([cell.volume for cell in self.list_of_cells])
                self.voxel_number = voxel_number
                self.dose = 0
        def random_points_in_voxel(self, n):
                points = np.random.uniform(-self.half_length, self.half_length, (n,3))
                points = points + self.position
                return points
        def add_cell(self, cell):
                #print('adding cell to voxel number : ', self.voxel_number)
                self.list_of_cells = np.append(cell,self.list_of_cells)
                self.free_space = self.free_space - cell.volume
                if self.free_space < 0:
                        raise ValueError('Voxel is full')

        def remove_cell(self, cell):
                self.list_of_cells = np.delete(self.list_of_cells, np.where(self.list_of_cells == cell))
                self.free_space = self.free_space + cell.volume

        def plot_vox(self, ax, fig, color='black'):
                plot_cube(ax, fig, self.position, self.half_length, color,)
                return fig, ax

        def plot_cells(self):
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                points = self.random_points_in_voxel(len(self.list_of_cells))
                color = sum([cell.color for cell in self.list_of_cells])
                color = color/len(self.list_of_cells)
                ax.scatter(points[:,0], points[:,1], points[:,2], color = color, alpha=0.5)

        def plot_vox_and_cells(self):
                fig, ax = self.plot_vox()
                points = self.random_points_in_voxel(len(self.list_of_cells))
                i = 0
                for cell in self.list_of_cells:
                        ax.scatter(points[i, 0], points[i, 1], points[i, 2], color=cell.color, alpha=0.5)
                        i = i + 1
                return fig, ax

        def update_cells_afterRT(self):
               RADIOSENSITIVITY = 0.1 #1% of cells killed per Gy
               expected_deaths = self.dose*RADIOSENSITIVITY
               expected_senescent = expected_deaths * 10
               # Use a Poisson distribution to model the number of cell deaths
               num_deaths = np.random.poisson(expected_deaths)
               num_deaths = min(num_deaths, len(self.list_of_cells))
               cells_to_remove = np.random.choice(self.list_of_cells, size=num_deaths, replace=False)
               for cell in cells_to_remove:
                       cell.state = 'dead' # kill the cell



