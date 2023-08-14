import numpy as np
import scipy.sparse as sparse


def sigmoid(L, x, x0, k): #sigmoid function
    return L/(1 + np.exp(-k*(x-x0)))

class Voxel(object): #extra parameters are max_occupancy, viscosity
        def __init__(self, position, half_length, viscosity, list_of_cells_in=None, n_capillaries = 0, voxel_number = 0):
                if list_of_cells_in is None:
                        list_of_cells_in = np.array([])
                self.position = position #position is a 3D vector
                self.half_length = half_length #half_length is a scalar
                self.list_of_cells = list_of_cells_in #list of cells in the voxel
                self.list_of_dead_cells = [] #list of necrotic cells in the voxel
                self.n_capillaries = n_capillaries #number of capillaries in the voxel
                self.volume = 8*half_length**3 #volume of the voxel
                self.voxel_number = voxel_number #id of the voxel
                self.dose = 0 #dose in the voxel from last irradiation
                self.molecular_factors = {'VEGF': 0} #molecular factors in the voxel. Only VEGF is implemented
                self.viscosity = viscosity #viscosity of the voxel
                self.vessel_volume = 0 #volume of the vessels in the voxel
                self.vessel_length = 0 #length of the vessels in the voxel
                self.bifurcation_density = 0 #bifurcation density in the voxel
                self.pH = 7.4 #pH in the voxel

        def number_of_tumor_cells(self): #returns the number of tumor cells in the voxel
                number = 0
                for cell in self.list_of_cells:
                        if cell.type == 'TumorCell':
                                number = number + 1
                return number

        def number_of_dead_cells(self): #returns the number of necrotic cells in the voxel
                return len(self.list_of_dead_cells)

        def number_of_necrotic_cells(self): #returns the number of apoptotic cells in the voxel
                number = 0
                for cell in self.list_of_dead_cells:
                        if cell.necrotic == True:
                                number = number + 1
                return number

        def number_of_apoptotic_cells(self): #returns the number of apoptotic cells in the voxel
                number = 0
                for cell in self.list_of_dead_cells:
                        if cell.necrotic == False:
                                number = number + 1
                return number

        def number_of_alive_cells(self): #returns the number of alive cells in the voxel
                return len(self.list_of_cells)
        def occupied_volume(self): #returns the volume occupied by the cells in the voxel
                volume = 0.0
                for cell in self.list_of_cells:
                        volume = volume + cell.volume
                for cell in self.list_of_dead_cells:
                        volume = volume + cell.volume
                return volume

        def vessel_volume_density(self): #returns the vessel volume density in the voxel
                side = 2*self.half_length
                capillary_volume = side * np.pi * 0.002 ** 2
                vessel_volume_density = ((self.vessel_volume + self.n_capillaries * capillary_volume) / self.volume)*100
                return vessel_volume_density
        def vessel_length_density(self): #returns the vessel length density in the voxel
                side = 2*self.half_length
                vessel_length_density = (self.vessel_length + self.n_capillaries * side) / self.volume
                return vessel_length_density

        def pressure(self): #returns the pressure in the voxel
                packing_density = (self.occupied_volume() /self.volume)
                return packing_density

        def random_points_in_voxel(self, n): #returns n random points in the voxel
                points = np.random.uniform(-self.half_length, self.half_length, (n,3))
                points = points + self.position
                return points
        def add_cell(self, cell, max_occupancy): #try to add a cell to the voxel
                if self.pressure() > max_occupancy:
                        #print('Voxel is full, pressure is', self.pressure(), ' number of cells is', self.number_of_alive_cells(), ' and number of necrotic cells is', self.number_of_necrotic_cells())
                        return False
                else:
                        self.list_of_cells = np.append(cell, self.list_of_cells)
                        return True

        def remove_cell(self, cell): #remove a cell from the voxel
                id = np.where(self.list_of_cells == cell)
                self.list_of_cells = np.delete(self.list_of_cells, id)
                return True

        def remove_dead_cell(self, cell):
                indices = np.nonzero(self.list_of_dead_cells == cell)[0]
                if indices.size == 0:
                        return False  # Cell not found in the list
                self.list_of_dead_cells = np.delete(self.list_of_dead_cells, indices)
                return True

        def cell_becomes_necrotic(self, cell): #remove a cell from the voxel and add it to the list of necrotic cells
                self.list_of_cells = np.delete(self.list_of_cells, np.where(self.list_of_cells == cell))
                cell.necrotic = True
                self.list_of_dead_cells = np.append(cell, self.list_of_dead_cells)
                return True

        def cell_becomes_apoptotic(self, cell):
                self.list_of_cells = np.delete(self.list_of_cells, np.where(self.list_of_cells == cell))
                self.list_of_dead_cells = np.append(cell, self.list_of_dead_cells)
                return True

        def oxygen_histogram(self, ax, fig): #plot the oxygen histogram of the voxel
                oxygen = []
                for cell in self.list_of_cells:
                        oxygen = np.append(oxygen, cell.capillaries)
                ax.hist(oxygen, bins = 50, color = 'blue', alpha = 0.5, range = (0,1))
                ax.set_xlim(0, 1)
                ax.set_ylabel('Number of cells')
                ax.set_title('Oxygen histogram')
                return ax, fig

        def vitality_histogram(self, ax, fig): #plot the vitality histogram of the voxel
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

        def cycling_time_and_age_histogram(self, ax, fig): #plot the cycling time and age histogram of the voxel

                if self.number_of_tumor_cells() == 0:
                        return ax, fig

                from scipy.stats import gamma
                gamma_scale = self.list_of_cells[0].gamma_scale
                gamma_shape = self.list_of_cells[0].gamma_shape
                cycling_time = []
                age = []
                for cell in self.list_of_cells:
                        cycling_time.append(cell.doubling_time)
                        age.append(cell.time_spent_cycling)
                ax.hist(cycling_time, bins=20, color='green', alpha=0.5, label='Time before next division', density=True)
                x = np.linspace(0, 70, 100)
                ax.hist(age, bins=20, color='red', alpha=0.5, label='Time spent cycling', density=True)
                ax.plot(x, gamma.pdf(x, gamma_shape, scale=gamma_scale), linestyle = 'dashed', color = 'black', lw=4, alpha=1.0, label='Gamma distribution')
                ax.set_xlabel('Time [h]')
                ax.set_ylabel('Frequency')
                ax.set_title('Cycling time histogram')
                ax.legend()
                return ax, fig

        def compute_cell_interaction_matrix(self, dt): #compute the cell interaction matrix in the voxel
                cell_interaction_matrix = sparse.lil_matrix((len(self.list_of_cells), len(self.list_of_cells)), dtype=np.float32)
                for i in range(len(self.list_of_cells)):
                        for j in range(len(self.list_of_cells)):
                                if i != j:
                                        cell_interaction_matrix[i,j] = self.list_of_cells[i].probability_of_interaction(self.list_of_cells[j], dt)
                cell_interaction_matrix = cell_interaction_matrix.tocsr()
                return cell_interaction_matrix

        def average_cell_damage(self): #returns the average damage of the cells in the voxel
                if len(self.list_of_cells) == 0:
                        return -1
                damage = 0
                for cell in self.list_of_cells:
                        damage = damage + cell.damage
                return damage / len(self.list_of_cells)









