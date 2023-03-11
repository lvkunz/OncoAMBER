from Voxel import Voxel
from Cell import Cell, TumorCell, HealthyCell
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BasicPlots import *
from scipy.sparse import csr_matrix
import sys

#np.set_printoptions(threshold=sys.maxsize)


class World:
    def __init__(self, half_length, number_of_voxels : int = 100):

        #raise error if half length is not a multiple of number of voxels
        if half_length % number_of_voxels != 0:
            raise ValueError('half_length must be a multiple of number_of_voxels')

        self.half_length = half_length
        self.voxel_list = []
        self.total_number_of_voxels = number_of_voxels**3
        for i in range(number_of_voxels):
            for j in range(number_of_voxels):
                for k in range(number_of_voxels):
                    self.voxel_list.append(Voxel(np.array([i,j,k]), half_length/number_of_voxels, voxel_number = i*number_of_voxels**2 + j*number_of_voxels + k))
        #self.vascular_network = VascularNetwork()
        self.number_of_voxels = number_of_voxels
    def find_voxel(self, position):
        #find the voxel that contains the position given
        #position is a numpy array
        #returns the voxel object
        voxel_number = int((position[0] + self.half_length)/(2*self.half_length)*self.total_number_of_voxels)
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
    def show_voxels_centers(self, ax, fig):
        #plots all the voxel centers in the world
        for voxel in self.voxel_list:
            # print(voxel)
            # print(voxel.position)
            # print(len(voxel.list_of_cells))
            #choose the color of the main type of cell in the voxel
            ax.scatter(voxel.position[0], voxel.position[1], voxel.position[2], s= len(voxel.list_of_cells) , color='green', alpha=0.5)
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
            f1 = voxel.free_space
            neighbors_voxel = self.find_neighbors(voxel)
            for neighbor in neighbors_voxel:
                f2 = neighbor.free_space
                j = neighbor.voxel_number
                exchange_matrix[i][j] = 1 - np.exp((-2*(d**2)*N*f1*dt)/(volume*(f1+f2))) #this is the probability of a cell moving from voxel i to voxel j
        print(exchange_matrix)
        return exchange_matrix
