from Voxel import *
from Cell import Cell, TumorCell, HealthyCell
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BasicPlots import *
from scipy.sparse import csr_matrix
import sys
from Vessel import *
from BasicGeometries import *
#np.set_printoptions(threshold=sys.maxsize)
from ScalarField import *
from matplotlib.colors import Normalize


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
        self.vasculature : VasculatureNetwork

    def generate_vasculature(self, num_vessels):
        points = self.random_points_for_voxels_concentration(num_vessels+1, 'VEGF')
        self.vasculature = VasculatureNetwork()
        self.vasculature.add_vessel(points[0], points[1], 0.01)
        self.vasculature.grow_vasculature(points[2:])
        self.update_list_of_vessels_in_voxels()
        self.update_volume_occupied_by_vessels()
        return self.vasculature

    def vasculature_growth(self, num_vessels):
        print('-- Growing ' + str(num_vessels) + ' new vessels')
        points = self.random_points_for_voxels_concentration(num_vessels, 'VEGF')
        print('-- Adding new vessels')
        self.vasculature.grow_vasculature(points)
        self.update_list_of_vessels_in_voxels()
        self.update_volume_occupied_by_vessels()
        return

    def read_vasculature(self, path):
        self.vasculature = VasculatureNetwork()
        self.vasculature.read(path)
        self.update_list_of_vessels_in_voxels()
        self.update_volume_occupied_by_vessels()
        return

    def random_points_for_voxels_concentration(self, num_points, molecule : str):
        if num_points == 0:
            return []
        if num_points % 1000 == 0:
            print('-- Generating ' + str(num_points) + ' random points for ' + molecule + ' concentration, if no progress is shown, VEGF concentration might be too low')
        points = []
        copy = self.voxel_list.copy()
        while len(points) < num_points:
            copy = [voxel for voxel in copy if voxel.molecular_factors[molecule] > 0]
            if len(copy) == 0:
                print('No more voxels with ' + molecule + ' concentration')
                return []
            copy = np.array(copy)
            np.random.shuffle(copy)
            for voxel in copy:
                mol = voxel.molecular_factors[molecule]
                if np.random.random() < mol:
                    points.append(voxel.random_points_in_voxel(1)[0])
        points = points[0:num_points]
        points = np.array(points)
        np.random.shuffle(points)
        return points


    def update_volume_occupied_by_vessels(self, resolution = 10):
        scorer = np.zeros((self.total_number_of_voxels))
        print('-- Computing volume occupied by vessels in each voxel')
        for vessel in self.vasculature.list_of_vessels:
            points = vessel.generate_points_along_axis(resolution)
            for point in points:
                voxel_n = self.find_voxel_number(point)
                scorer[voxel_n] += np.pi*(vessel.length/resolution)*vessel.radius**2
        print('-- Finishing up')
        for voxel in self.voxel_list:
            voxel.vessel_volume = scorer[voxel.voxel_number]
        return

    def update_biology_after_RT(self):
        print('-- Updating biology after RT')
        for voxel in self.voxel_list:
            voxel.update_cells_afterRT()
        self.vessels_killed_by_dose()
        self.compute_oxygen_map()
    def vessels_killed_by_pressure(self, pressure_killing_threshold, pressure_killing_slope):
        death = lambda pressure: sigmoid(1, pressure, pressure_killing_threshold, pressure_killing_slope)
        count = 0
        for voxel in self.voxel_list:
            for vessel in voxel.list_of_vessels_ids:
                if np.random.random() < death(voxel.pressure()):
                    self.vasculature.remove_vessel(vessel)
                    voxel.list_of_vessels_ids.remove(vessel)
                    count += 1
        return count
    def vessels_killed_by_dose(self):
        death = lambda dose: sigmoid(1, dose, 6e-8, 200000000)
        count = 0
        for voxel in self.voxel_list:
            for vessel in voxel.list_of_vessels_ids:
                if np.random.random() < death(voxel.dose):
                    self.vasculature.remove_vessel(vessel)
                    voxel.list_of_vessels_ids.remove(id)
                    count += 1
        return count

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

    def update_list_of_vessels_in_voxels(self):
        for vessel in self.vasculature.list_of_vessels:
            current_list = self.find_voxel(vessel.center).list_of_vessels_ids
            if vessel.id not in current_list:
                current_list.append(vessel.id)
        return
    def compute_exchange_matrix(self, dt, pressure_threshold):
        hL = self.voxel_list[0].half_length #half length of the voxel
        #computes the matrix that describes the exchange of cells between voxels
        #returns the matrix
        #the matrix is a list of lists
        #the first index is the voxel number
        exchange_matrix = csr_matrix((self.total_number_of_voxels, self.total_number_of_voxels)).toarray()
        for voxel in self.voxel_list:
            neighbors_voxel = self.find_neighbors(voxel)
            for neighbor in neighbors_voxel:
                DeltaP = voxel.pressure() - neighbor.pressure()
                if DeltaP <= pressure_threshold: #if the pressure in the neighbor is higher than in the voxel, no exchange,
                    exchange_matrix[voxel.voxel_number][neighbor.voxel_number] = 0
                else:
                    t_res = 3 * hL ** 2 / (DeltaP * voxel.viscosity)
                    ratio = t_res / dt
                    exchange_matrix[voxel.voxel_number][neighbor.voxel_number] = ratio * np.exp(-ratio)#this is the probability of a cell moving from voxel i to voxel j
        return exchange_matrix

    def topas_param_file(self, name : str):
        print('-- Creating parameter file for topas simulation, file :', name)
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
        print('-- File saved as ' + name)
        return

    def update_dose(self, doses):
        #updates the dose in each voxel
        if len(doses) != self.total_number_of_voxels:
            print('Error: the number of doses is not equal to the number of voxels, Probably due to a unmatching Topas simulation')
            return
        for voxel in self.voxel_list:
            voxel.dose = doses[voxel.voxel_number]
        return

    def compute_oxygen_map(self, o2_per_volume=10, diffusion_number=5):
        print('-- Computing oxygen map')
        for voxel in self.voxel_list:
            voxel.oxygen = voxel.vessel_volume * o2_per_volume
        for i in range(diffusion_number):
            new_oxygen_map = np.zeros(self.total_number_of_voxels)
            print('--- o2 map computing', i, 'out of', diffusion_number)
            for voxel in self.voxel_list:
                sum = voxel.oxygen
                for neighbor in self.find_neighbors(voxel):
                    sum += neighbor.oxygen
                new_oxygen_map[voxel.voxel_number] = sum / (1 + len(self.find_neighbors(voxel)))
            for voxel in self.voxel_list:
                voxel.oxygen = new_oxygen_map[voxel.voxel_number]
        return
    #function to plot cells

    def show_voxels_centers(self, ax, fig, slice = False):
        print('-- Plotting Cell Population')

        #central z slice of world first voxel is at z = 0
        if slice:
            middle_slice = self.number_of_voxels // 2
            first_voxel = middle_slice * self.number_of_voxels**2
            last_voxel = (middle_slice + 1) * self.number_of_voxels**2 -1
            list = self.voxel_list[first_voxel:last_voxel]
            size = 250/self.number_of_voxels
        else:
            list = self.voxel_list
            size = 1
        number = []
        positions = []
        # collect doses and positions for all voxels
        for voxel in list:
            number.append(voxel.number_cells)
            positions.append(voxel.position)

        ax.scatter(
        [p[0] for p in positions],
        [p[1] for p in positions],
        [p[2] for p in positions],
        c=number, cmap='YlGn', alpha= 0.5, vmin=min(number), vmax=max(number), s= size)
        # add colorbar
        fig.colorbar(ax.collections[0])
        return fig, ax

    def show_voxels_centers_celltype(self, ax, fig):
        for voxel in self.voxel_list:
            # print(voxel)
            # print(voxel.position)
            # print(len(voxel.list_of_cells))
            average_color = (0, 0, 0)
            for cell in voxel.list_of_cells:
                # transform the string color into rgb
                rgb = to_rgb(cell.color)
                average_color = (average_color[0] + rgb[0], average_color[1] + rgb[1], average_color[2] + rgb[2])
            if len(voxel.list_of_cells) > 0:  # if the voxel is empty we'd have an error
                average_color = (
                average_color[0] / len(voxel.list_of_cells), average_color[1] / len(voxel.list_of_cells),
                average_color[2] / len(voxel.list_of_cells))
                average_color = (int(average_color[0]), int(average_color[1]), int(average_color[2]))
            color = rgb_to_hex(average_color)
            ax.scatter(voxel.position[0], voxel.position[1], voxel.position[2], s=len(voxel.list_of_cells), color=color,
                       alpha=0.3)
        return fig, ax

    def show_voxels_centers_dose(self, ax, fig):
        print('-- Plotting Dose')
        doses = []
        positions = []
        # collect doses and positions for all voxels
        for voxel in self.voxel_list:
            doses.append(voxel.dose)
            positions.append(voxel.position)
        ax.scatter(
            [p[0] for p in positions],
            [p[1] for p in positions],
            [p[2] for p in positions],
            c=doses, cmap='BuPu', alpha=0.5, vmin=min(doses), vmax=max(doses)
        )
        # add colorbar
        fig.colorbar(ax.collections[0])
        return fig, ax

    def show_voxels_centers_oxygen(self, ax, fig, slice = False):
        print('-- Plotting Oxygen')
        if slice:
            middle_slice = self.number_of_voxels // 2
            first_voxel = middle_slice * self.number_of_voxels**2
            last_voxel = (middle_slice + 1) * self.number_of_voxels**2 -1
            list = self.voxel_list[first_voxel:last_voxel]
            size = 250/self.number_of_voxels
        else:
            list = self.voxel_list
            size = 1
        oxygen = []
        positions = []
        # collect doses and positions for all voxels
        for voxel in list:
            oxygen.append(voxel.oxygen)
            positions.append(voxel.position)
        ax.scatter(
            [p[0] for p in positions],
            [p[1] for p in positions],
            [p[2] for p in positions],
            c=oxygen, cmap='Blues', alpha=0.5, vmin=min(oxygen), vmax=max(oxygen), s= size
        )
        # add colorbar
        fig.colorbar(ax.collections[0])
        return fig, ax
    def show_voxels_centers_pressure(self, ax, fig, slice = False):
        print('-- Plotting Pressure')
        if slice:
            middle_slice = self.number_of_voxels // 2
            first_voxel = middle_slice * self.number_of_voxels ** 2
            last_voxel = (middle_slice + 1) * self.number_of_voxels ** 2 - 1
            list = self.voxel_list[first_voxel:last_voxel]
            size = 250 / self.number_of_voxels
        else:
            list = self.voxel_list
            size = 1
        pressure = []
        positions = []
        for voxel in list:
            pressure.append(voxel.pressure())
            positions.append(voxel.position)
        ax.scatter(
            [p[0] for p in positions],
            [p[1] for p in positions],
            [p[2] for p in positions],
            c=pressure, cmap='RdPu', alpha=0.5, vmin=min(pressure), vmax=max(pressure), s=size
        )
        # add colorbar
        fig.colorbar(ax.collections[0])
        return fig, ax

    def show_voxels_centers_molecules(self, ax, fig, molecule : str, slice = False):
        print('-- Plotting Molecules')
        if slice:
            middle_slice = self.number_of_voxels // 2
            first_voxel = middle_slice * self.number_of_voxels ** 2
            last_voxel = (middle_slice + 1) * self.number_of_voxels ** 2 - 1
            list = self.voxel_list[first_voxel:last_voxel]
            size = 250 / self.number_of_voxels
        else:
            list = self.voxel_list
            size = 1
        molecules = []
        positions = []
        # collect doses and positions for all voxels
        for voxel in list:
            molecules.append(voxel.molecular_factors[molecule])
            positions.append(voxel.position)
        ax.scatter(
            [p[0] for p in positions],
            [p[1] for p in positions],
            [p[2] for p in positions],
            c=molecules, cmap='Oranges', alpha=0.5, vmin=min(molecules), vmax=max(molecules), s=size
        )
        # add colorbar
        fig.colorbar(ax.collections[0])
        return fig, ax