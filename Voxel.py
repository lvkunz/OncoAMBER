from Cell import Cell, TumorCell, HealthyCell
from BasicPlots import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
def sigmoid(L, x, x0, k):
    return L/(1 + np.exp(-k*(x-x0)))

class Voxel(object):
        def __init__(self, position = np.array([0,0,0]), half_length = 0.1, list_of_cells_in = [], oxygen = 0, voxel_number = 0):
                self.position = position
                self.half_length = half_length
                self.list_of_cells = list_of_cells_in
                self.oxygen = oxygen
                self.volume = 8*half_length**3
                self.occupied_volume = sum([cell.volume for cell in self.list_of_cells])
                self.voxel_number = voxel_number
                self.dose = 0
                self.number_cells = len(self.list_of_cells)
                self.molecular_factors = {'EGF': 0, 'FGF': 0, 'HGF': 0, 'IGF': 0, 'TGF': 0, 'VEGF': 0, 'WNT': 0}
                self.molecular_factors['VEGF'] = 0.1
                # if np.linalg.norm(self.position) < 7.0:
                #         self.molecular_factors['VEGF'] = 1.0
                self.list_of_vessels_ids = []
        def pressure(self): #units of pressure are pascals
                # print('x', len(self.list_of_cells))
                # print('y', self.occupied_volume)
                # B = 5.7e-4 #system of hard spheres interacting via hard core repulsion at 37 degrees C
                # kT = 4.1e-3 #4.1e-21 #kT at 37 degrees C
                # if len(self.list_of_cells) == 0:
                #         return 0
                # packing_density = self.occupied_volume/self.volume
                # print('packing density', packing_density)
                # number_density = packing_density
                # print('number density', number_density)
                # ratio = number_density/(1-number_density)
                # print('ratio', ratio)
                # ratio2 = ratio**2
                # pressure = kT * ratio - B * ratio2
                # print('pressure', pressure)
                # return pressure
                NkT = 1e5 #use somethig sensitive to make sense of this
                if len(self.list_of_cells) == 0:
                        return 0
                packing_density = self.occupied_volume/self.volume
                ratio = (1+packing_density+packing_density**2-packing_density**3)/(1-packing_density**3) #using Carrahan-Sterling equation
                pressure = (NkT*ratio)/self.volume
                return pressure

        def random_points_in_voxel(self, n):
                points = np.random.uniform(-self.half_length, self.half_length, (n,3))
                points = points + self.position
                return points
        def add_cell(self, cell):
                #print('adding cell to voxel number : ', self.voxel_number)
                self.list_of_cells = np.append(cell,self.list_of_cells)
                self.occupied_volume = self.occupied_volume + cell.volume
                self.number_cells = self.number_cells + 1
                if self.occupied_volume > self.volume:
                        raise ValueError('Voxel is full')

        def remove_cell(self, cell):
                self.occupied_volume = self.occupied_volume - cell.volume
                id = np.where(self.list_of_cells == cell)
                self.list_of_cells = np.delete(self.list_of_cells, id)
                self.number_cells = self.number_cells - 1
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

        def update_cells_for_oxygen_state(self):
                #print('updating cells for oxygen state')
                #define sigmoid function depending on oxygen (logistic function)
                death = lambda pO2: sigmoid(1, pO2, 0.15, 30)
                senescence = lambda pO2: sigmoid(1, pO2, 0.4, 10)
                for cell in self.list_of_cells:
                        sample = np.random.random()
                        if sample > death(self.oxygen):
                                cell.state = 'dead'
                        elif sample > senescence(self.oxygen):
                                cell.state = 'senescent'
                        else:
                                cell.state = 'cycling'

        def update_molecules(self, dt):
                #print('updating molecules')
                count = 0
                for cell in self.list_of_cells:
                        if cell.state == 'cycling':
                                count = count + 1
                self.molecular_factors['VEGF'] = min(1, count*0.001)










