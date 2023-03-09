from Voxel import Voxel
from Cell import Cell, TumorCell, HealthyCell
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BasicPlots import *


class World:
    def __init__(self, half_length, number_of_voxels = 100):

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

    def find_voxel(self, position):
        #find the voxel that contains the position given
        #position is a numpy array
        #returns the voxel object
        voxel_number = int((position[0] + self.half_length)/(2*self.half_length)*self.total_number_of_voxels)
        return self.voxel_list[voxel_number]

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
        exchange_matrix = np.zeros((self.total_number_of_voxels, self.total_number_of_voxels))
        d = 2*self.voxel_list[0].half_length
        volume = d**3
        for voxel in self.voxel_list:
            i = voxel.voxel_number
            N = len(voxel.list_of_cells)
            f1 = voxel.free_space
            neigbors_voxel = find_neigbors(voxel)
            for neigbor in neigbors_voxel:
                f2 = neigbor.free_space
                j = neigbor.voxel_number
                    exchange_matrix[i][j] = 1 - np.exp((-2*(d**2)*N*f1*dt)/(volume*(f1+f2))) #this is the probability of a cell moving from voxel i to voxel j
        print(exchange_matrix)
        return exchange_matrix
    def find_neighbors(voxel):
        #finds the neigbors of a voxel
        #returns a list of voxel objects
        #voxel is a voxel object
        neigbors_voxel = []
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    if i == 0 and j == 0 and k == 0:
                        continue
                    else:
                        neigbors_voxel.append(self.find_voxel(voxel.position + np.array([i,j,k])))
        return neigbors_voxel