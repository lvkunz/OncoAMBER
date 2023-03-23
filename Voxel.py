from Cell import Cell, TumorCell, HealthyCell
from BasicPlots import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
def sigmoid(L, x, x0, k):
    return L/(1 + np.exp(-k*(x-x0)))

class Voxel(object):
        def __init__(self, position = np.array([0,0,0]), half_length = 0.1, list_of_cells_in=None, oxygen = 0, voxel_number = 0):
                if list_of_cells_in is None:
                        list_of_cells_in = []
                self.position = position
                self.half_length = half_length
                self.list_of_cells = list_of_cells_in
                self.oxygen = oxygen
                self.volume = 8*half_length**3
                self.voxel_number = voxel_number
                self.dose = 0
                self.molecular_factors = {'EGF': 0, 'FGF': 0, 'HGF': 0, 'IGF': 0, 'TGF': 0, 'VEGF': 0, 'WNT': 0}
                # if np.linalg.norm(self.position) < 5:
                #         self.molecular_factors['VEGF'] = 1.0
                self.list_of_vessels_ids = []
                self.vessel_volume = 0
                self.viscosity = 10

        def number_of_cells(self):
                return len(self.list_of_cells)

        def occupied_volume(self):
                volume = 0.0
                for cell in self.list_of_cells:
                        volume = volume + cell.volume
                return volume

        def pressure(self):
                packing_density = ( self.occupied_volume() /self.volume)
                # NkT = 1 #use somethig sensitive to make sense of this
                # if len(self.list_of_cells) == 0:
                #         return 0
                # ratio = (1+packing_density+packing_density**2-packing_density**3)/((1-packing_density)**3) #using Carrahan-Sterling equation
                # ratio = ratio - 1
                # pressure = (NkT*ratio)/self.volume
                return packing_density

        def random_points_in_voxel(self, n):
                points = np.random.uniform(-self.half_length, self.half_length, (n,3))
                points = points + self.position
                return points
        def add_cell(self, cell):
                if self.pressure() > 0.998:
                        print('Voxel is full, pressure is', self.pressure())
                        return False
                else:
                        self.list_of_cells = np.append(cell, self.list_of_cells)
                        return True

        def remove_cell(self, cell):
                id = np.where(self.list_of_cells == cell)
                self.list_of_cells = np.delete(self.list_of_cells, id)
                return True


        def update_cells_afterRT(self, radio_sensitivity):
               expected_deaths = self.dose*radio_sensitivity
               # Use a Poisson distribution to model the number of cell deaths
               num_deaths = np.random.poisson(expected_deaths)
               num_deaths = min(num_deaths, len(self.list_of_cells))
               cells_to_remove = np.random.choice(self.list_of_cells, size=num_deaths, replace=False)
               for cell in cells_to_remove:
                       cell.state = 'dead' # kill the cell

        def update_cells_for_oxygen_state(self, vitality_threshold, vitality_slope):
                #print('updating cells for oxygen state')
                #define sigmoid function depending on oxygen (logistic function)
                vitality = lambda pO2: sigmoid(1, pO2, vitality_threshold, vitality_slope)
                for cell in self.list_of_cells:
                        cell.vitality = vitality(self.oxygen)

        def update_occupied_volume(self):
                #print('updating occupied volume')
                volume = 0
                for cell in self.list_of_cells:
                        volume = volume + cell.volume
                self.occupied_volume = volume
                return volume


        def update_molecules(self, dt, VEGF_production_per_cell):
                #print('updating molecules')
                count = 0
                for cell in self.list_of_cells:
                        if isinstance(cell, TumorCell):
                                count = count + 1
                self.molecular_factors['VEGF'] = min(1, count*VEGF_production_per_cell)










