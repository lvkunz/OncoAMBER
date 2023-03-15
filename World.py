from Voxel import Voxel
from Cell import Cell, TumorCell, HealthyCell
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BasicPlots import *
from scipy.sparse import csr_matrix
import sys
from Vessel import *

#np.set_printoptions(threshold=sys.maxsize)


class World:
    def __init__(self, half_length, number_of_voxels : int = 20):
        self.half_length = half_length
        self.voxel_list = []
        self.total_number_of_voxels = number_of_voxels**3
        for i in range(number_of_voxels):
            for j in range(number_of_voxels):
                for k in range(number_of_voxels):
                    position = np.array([i*2*half_length/number_of_voxels - half_length, j*2*half_length/number_of_voxels - half_length, k*2*half_length/number_of_voxels - half_length])
                    self.voxel_list.append(Voxel(position, half_length/number_of_voxels, voxel_number = i*number_of_voxels**2 + j*number_of_voxels + k))
        self.number_of_voxels = number_of_voxels
        self.vasculature = VasculatureNetwork()

    def generate_vasculature(self, num_vessels):
        self.vasculature.generate_vasculature(num_vessels)
        return

    def compute_oxygen_map(self):
        for voxel in self.voxel_list:
            distance = self.vasculature.closest_distance(voxel.position)

    def find_voxel_number(self, position):
        #doesn't handle the case where the position is outside the world
        num_voxels = self.number_of_voxels
        i = int((position[0] + self.half_length) * num_voxels / (2 * self.half_length))
        j = int((position[1] + self.half_length) * num_voxels / (2 * self.half_length))
        k = int((position[2] + self.half_length) * num_voxels / (2 * self.half_length))
        n = i * num_voxels ** 2 + j * num_voxels + k
        if n >= self.total_number_of_voxels or n < 0:
            raise ValueError('position ' +str(position) + ' is outside the world (voxel number = ' + str(n) + ')' )
        return n
    def find_voxel(self, position):
        #find the voxel that contains the position given
        #position is a numpy array
        #returns the voxel object
        voxel_number = self.find_voxel_number(position)
        return self.voxel_list[voxel_number]
    def find_neighbors(self, voxel):
        #finds the neighbors of a voxel
        voxel_number = voxel.voxel_number
        num_voxels = self.number_of_voxels
        i = voxel_number // (num_voxels ** 2)
        j = (voxel_number // num_voxels) % num_voxels
        k = voxel_number % num_voxels

        neighbors = []
        if i > 0:
            neighbors.append((i - 1) * num_voxels ** 2 + j * num_voxels + k)
        if i < num_voxels - 1:
            neighbors.append((i + 1) * num_voxels ** 2 + j * num_voxels + k)
        if j > 0:
            neighbors.append(i * num_voxels ** 2 + (j - 1) * num_voxels + k)
        if j < num_voxels - 1:
            neighbors.append(i * num_voxels ** 2 + (j + 1) * num_voxels + k)
        if k > 0:
            neighbors.append(i * num_voxels ** 2 + j * num_voxels + (k - 1))
        if k < num_voxels - 1:
            neighbors.append(i * num_voxels ** 2 + j * num_voxels + (k + 1))

        neighbors = [n for n in neighbors if 0 <= n < num_voxels ** 3]
        neighbors_voxels = [self.voxel_list[n] for n in neighbors]
        return neighbors_voxels
    def show_voxels(self, ax, fig):
        #plots all the voxels in the world
        for voxel in self.voxel_list:
            voxel.plot_vox(ax,fig)
        return fig, ax
    def show_voxels_centers(self, ax, fig, colorful = False):
        #plots all the voxel centers in the world
        for voxel in self.voxel_list:
            # print(voxel)
            # print(voxel.position)
            # print(len(voxel.list_of_cells))
            average_color = (0, 0, 0)
            if colorful:
                for cell in voxel.list_of_cells:
                    #transform the string color into rgb
                    rgb = to_rgb(cell.color)
                    average_color = (average_color[0] + rgb[0], average_color[1] + rgb[1], average_color[2] + rgb[2])
                if len(voxel.list_of_cells) > 0: #if the voxel is empty we'd have an error
                    average_color = (average_color[0]/len(voxel.list_of_cells), average_color[1]/len(voxel.list_of_cells), average_color[2]/len(voxel.list_of_cells))
                    average_color = (int(average_color[0]), int(average_color[1]), int(average_color[2]))
            color = rgb_to_hex(average_color)
            ax.scatter(voxel.position[0], voxel.position[1], voxel.position[2], s= len(voxel.list_of_cells) , color=color, alpha = 0.3)
        return fig, ax

    def show_voxels_centers_dose(self, ax, fig):
        doses = []
        positions = []
        # collect doses and positions for all voxels
        for voxel in self.voxel_list:
            doses.append(voxel.dose)
            positions.append(voxel.position)
        # plot all voxels with their respective doses
        #alpha between 0.1 and 0.5 depending on dose
        ax.scatter(
            [p[0] for p in positions],
            [p[1] for p in positions],
            [p[2] for p in positions],
            c=doses, cmap='BuPu', alpha=0.3, vmin=min(doses), vmax=max(doses)
        )
        # add colorbar
        fig.colorbar(ax.collections[0])
        return fig, ax
    def compute_exchange_matrix(self, dt):
        #computes the matrix that describes the exchange of cells between voxels
        #returns the matrix
        #the matrix is a list of lists
        #the first index is the voxel number
        exchange_matrix = csr_matrix((self.total_number_of_voxels, self.total_number_of_voxels)).toarray()
        #exchange_matrix = np.zeros((self.total_number_of_voxels, self.total_number_of_voxels))
        d = 2*self.voxel_list[0].half_length
        volume = d**3
        for voxel in self.voxel_list:
            i = voxel.voxel_number
            N = len(voxel.list_of_cells)
            f1 = voxel.free_space/voxel.volume
            neighbors_voxel = self.find_neighbors(voxel)
            for neighbor in neighbors_voxel:
                f2 = neighbor.free_space/neighbor.volume
                j = neighbor.voxel_number
                exchange_matrix[i][j] = 1 - np.exp((-2*(d**2)*N*f1*dt)/(volume*(f1+f2))) #this is the probability of a cell moving from voxel i to voxel j
        #print(exchange_matrix)
        return exchange_matrix

    def topas_param_file(self, name : str = 'Params'):
        print('creating parameter file for topas')
        #returns a txt file that can be used as a parameter file for topas
        name = name + '.txt'
        #the file is saved in the TopasSimulation folder
        file = open('TopasSimulation/' + name, 'w')
        file.write('IncludeFile = BasicParameters.txt \n' +
                                                'd:Ge/MyBox/HLX      = ' + str(self.half_length*2) +' cm \n' +
                                                'd:Ge/MyBox/HLY      = ' + str(self.half_length*2) +' cm \n' +
                                                'd:Ge/MyBox/HLZ      = ' + str(self.half_length*2) +' cm \n' +
                                                'd:Ge/MyBox/TransX   = 0. m \n' +
                                                'd:Ge/MyBox/TransY   = 0. m \n' +
                                                'd:Ge/MyBox/TransZ   = 0. m \n' +
                                                'd:Ge/MyBox/RotX     = 0. deg \n' +
                                                'd:Ge/MyBox/RotY     = 0. deg \n' +
                                                'd:Ge/MyBox/RotZ     = 0. deg \n' +
                                                'i:Ge/MyBox/ZBins = ' + str(self.number_of_voxels) +' \n' +
                                                'i:Ge/MyBox/XBins = ' + str(self.number_of_voxels) +' \n' +
                                                'i:Ge/MyBox/YBins = ' + str(self.number_of_voxels))
        file.close()
        print('File saved as ' + name)
        return

    def update_dose(self, doses):
        #updates the dose in each voxel
        if len(doses) != self.total_number_of_voxels:
            print('Error: the number of doses is not equal to the number of voxels, Probably due to a unmatching Topas simulation')
            return
        for voxel in self.voxel_list:
            voxel.dose = doses[voxel.voxel_number]
        return
    #function to plot cells
