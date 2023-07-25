import numpy as np
import scipy.sparse as sparse


def sigmoid(L, x, x0, k):
    return L/(1 + np.exp(-k*(x-x0)))

class Voxel(object): #extra parameters are max_occupancy, viscosity
        def __init__(self, position, half_length, viscosity, list_of_cells_in=None, n_capillaries = 0, voxel_number = 0):
                if list_of_cells_in is None:
                        list_of_cells_in = np.array([])
                self.position = position
                self.half_length = half_length
                self.list_of_cells = list_of_cells_in
                self.list_of_necrotic_cells = np.array([])
                self.n_capillaries = n_capillaries
                self.volume = 8*half_length**3
                self.voxel_number = voxel_number
                self.dose = 0
                self.molecular_factors = {'VEGF': 0}
                self.viscosity = viscosity
                self.vessel_volume = 0
                self.vessel_length = 0
                self.bifurcation_density = 0
                self.pH = 7.4

        def number_of_tumor_cells(self):
                number = 0
                for cell in self.list_of_cells:
                        if cell.type == 'TumorCell':
                                number = number + 1
                return number

        def number_of_necrotic_cells(self):
                number = 0
                for cell in self.list_of_necrotic_cells:
                        if cell.type == 'TumorCell':
                                number = number + 1
                return number
        def number_of_alive_cells(self):
                return len(self.list_of_cells)
        def occupied_volume(self):
                volume = 0.0
                for cell in self.list_of_cells:
                        volume = volume + cell.volume
                for cell in self.list_of_necrotic_cells:
                        volume = volume + cell.volume
                return volume

        def vessel_volume_density(self):
                side = 2*self.half_length
                capillary_volume = side * np.pi * 0.002 ** 2
                vessel_volume_density = ((self.vessel_volume + self.n_capillaries * capillary_volume) / self.volume)*100
                return vessel_volume_density
        def vessel_length_density(self):
                side = 2*self.half_length
                vessel_length_density = (self.vessel_length + self.n_capillaries * side) / self.volume
                return vessel_length_density

        def pressure(self):
                packing_density = (self.occupied_volume() /self.volume)
                return packing_density

        def random_points_in_voxel(self, n):
                points = np.random.uniform(-self.half_length, self.half_length, (n,3))
                points = points + self.position
                return points
        def add_cell(self, cell, max_occupancy):
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
                        oxygen = np.append(oxygen, cell.capillaries)
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
                ax.set_xlim(0, 1)
                ax.set_xlabel(
                        f'Vitality, Number of capillaries in voxel was = {self.n_capillaries} and number of cells = {self.number_of_alive_cells()}')
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

        def compute_cell_interaction_matrix(self, dt):
                cell_interaction_matrix = sparse.lil_matrix((len(self.list_of_cells), len(self.list_of_cells)), dtype=np.float32)
                for i in range(len(self.list_of_cells)):
                        for j in range(len(self.list_of_cells)):
                                if i != j:
                                        cell_interaction_matrix[i,j] = self.list_of_cells[i].probability_of_interaction(self.list_of_cells[j], dt)
                cell_interaction_matrix = cell_interaction_matrix.tocsr()
                return cell_interaction_matrix

        def average_cell_damage(self):
                if len(self.list_of_cells) == 0:
                        return -1
                damage = 0
                for cell in self.list_of_cells:
                        damage = damage + cell.damage
                return damage / len(self.list_of_cells)









