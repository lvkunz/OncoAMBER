import numpy as np
from amber.config_instance import config

def sigmoid(L, x, x0, k):
    return L/(1 + np.exp(-k*(x-x0)))

class Voxel(object): #extra parameters are max_occupancy, viscosity
        def __init__(self, position = np.array([0,0,0]), half_length = 0.1, list_of_cells_in=None, oxygen = 0, voxel_number = 0):
                if list_of_cells_in is None:
                        list_of_cells_in = []
                self.position = position
                self.half_length = half_length
                self.list_of_cells = list_of_cells_in
                self.list_of_necrotic_cells = []
                self.oxygen = oxygen
                self.volume = 8*half_length**3
                self.voxel_number = voxel_number
                self.dose = 0
                self.molecular_factors = {'EGF': 0, 'FGF': 0, 'HGF': 0, 'IGF': 0, 'TGF': 0, 'VEGF': 0, 'WNT': 0}
                # if np.linalg.norm(self.position) < 5:
                #         self.molecular_factors['VEGF'] = 1.0
                self.viscosity = config.viscosity
                self.vessel_volume = 0

        def number_of_tumor_cells(self):
                number = 0
                for cell in self.list_of_cells:
                        if cell.type == 'TumorCell':
                                number = number + 1
                return number

        def number_of_necrotic_cells(self):
                return len(self.list_of_necrotic_cells)
        def number_of_alive_cells(self):
                return len(self.list_of_cells)
        def occupied_volume(self):
                volume = 0.0
                for cell in self.list_of_cells:
                        volume = volume + cell.volume
                for cell in self.list_of_necrotic_cells:
                        volume = volume + cell.volume
                return volume

        def pressure(self):
                packing_density = (self.occupied_volume() /self.volume)
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
                max_occupancy = config.max_occupancy #hard spheres is 0.64, + consider a little bit of compression
                if self.pressure() > max_occupancy:
                        #print('Voxel is full, pressure is', self.pressure(), ' number of cells is', self.number_of_alive_cells(), ' and number of necrotic cells is', self.number_of_necrotic_cells())
                        return False
                else:
                        self.list_of_cells = np.append(cell, self.list_of_cells)
                        return True

        def remove_cell(self, cell):
                id = np.where(self.list_of_cells == cell)
                self.list_of_cells = np.delete(self.list_of_cells, id)
                return True

        def remove_necrotic_cell(self, cell):
                id = np.where(self.list_of_necrotic_cells == cell)
                self.list_of_necrotic_cells = np.delete(self.list_of_necrotic_cells, id)
                return True

        def cell_becomes_necrotic(self, cell):
                self.list_of_cells = np.delete(self.list_of_cells, np.where(self.list_of_cells == cell))
                cell.necrotic = True
                self.list_of_necrotic_cells = np.append(cell, self.list_of_necrotic_cells)
                return True

        def oxygen_histogram(self, ax, fig):
                oxygen = []
                for cell in self.list_of_cells:
                        oxygen = np.append(oxygen, cell.oxygen)
                ax.hist(oxygen, bins = 50, color = 'blue', alpha = 0.5, range = (0,1))
                ax.set_xlim(0, 1)
                ax.set_ylabel('Number of cells')
                ax.set_title('Oxygen histogram')
                return ax, fig

        def vitality_histogram(self, ax, fig):
                vitality = []
                for cell in self.list_of_cells:
                        vitality.append(cell.vitality())
                ax.hist(vitality, bins=50, color='orange', alpha=0.5, range=(0, 1))
                ax.vlines(config.vitality_cycling_threshold, 0, ax.get_ylim()[1], colors='darkgreen',
                          linestyles='dashed')
                ax.vlines(config.vitality_apoptosis_threshold, 0, ax.get_ylim()[1], colors='darkred',
                          linestyles='dashed')
                ax.vlines(config.vitality_necrosis_threshold, 0, ax.get_ylim()[1], colors='black',
                          linestyles='dashed')
                ax.vlines(config.o2_threshold_for_VEGF_production, 0, ax.get_ylim()[1], colors='purple', linestyles='dashed')
                ax.set_xlim(0, 1)
                ax.set_xlabel(
                        f'Vitality, Oxygen in voxel was = {self.oxygen} and number of cells = {self.number_of_alive_cells()}')
                ax.set_ylabel('Number of cells')
                ax.set_title('Vitality histogram')
                return ax, fig

        def cycling_time_and_age_histogram(self, ax, fig):
                cycling_time = []
                age = []
                for cell in self.list_of_cells:
                        cycling_time.append(cell.doubling_time)
                        age.append(cell.time_spent_cycling)
                ax.hist(cycling_time, bins=30, color='green', alpha=0.5, label='Cycling time')
                ax.hist(age, bins=30, color='red', alpha=0.5, label='time_spent_cycling')
                ax.set_xlabel('Cycling time')
                ax.set_ylabel('Number of cells')
                ax.set_title('Cycling time histogram')
                ax.legend()
                return ax, fig










