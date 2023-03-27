from Voxel import *
from Cell import Cell, TumorCell, HealthyCell
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BasicPlots import *
from scipy.sparse import csr_matrix
import sys
from Vesselv2 import *
from BasicGeometries import *
#np.set_printoptions(threshold=sys.maxsize)
from ScalarField import *
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
from scipy import ndimage

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

    def initiate_vasculature(self, list_of_mother_vessels):
        self.vasculature = VasculatureNetwork(list_of_mother_vessels)
        return

    def vasculature_growth(self, dt, splitting_rate, macro_steps=1, micro_steps=10, weight_direction=0.5, weight_vegf=0.5, pressure_threshold=0.5, weight_pressure=0.5, radius_pressure_sensitive = False):
        print('Vasculature growth')
        pressure_map = self.pressure_map()
        print(1)
        def pressure(point): return pressure_map.evaluate(point)
        print(2)
        vegf = self.vegf_map()
        print(3)
        def vegf_gradient(point): return vegf.gradient(point)

        self.vasculature.grow_and_split(dt, splitting_rate, vegf_gradient, pressure, macro_steps, micro_steps, weight_direction, weight_vegf, pressure_threshold, weight_pressure)
        self.vasculature.update_vessels_radius(1.0,radius_pressure_sensitive, pressure)
        return

    def read_vasculature(self, path):
        pass

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


    def update_volume_occupied_by_vessels(self):
        scorer = np.zeros((self.total_number_of_voxels))
        print('-- Computing volume occupied by vessels in each voxel')
        for vessel in self.vasculature.list_of_vessels:
            for point in vessel.path:
                voxel_n = self.find_voxel_number(point)
                if voxel_n == -1:
                    continue
                scorer[voxel_n] += vessel.volume_per_point()
        print('-- Finishing up -- CHECK IF NECESSARY')
        for voxel in self.voxel_list:
            voxel.vessel_volume = scorer[voxel.voxel_number]
        return

    def update_biology_after_RT(self):
        print('-- Updating biology after RT')
        for voxel in self.voxel_list:
            voxel.update_cells_afterRT()
        #self.vessels_killed_by_dose()
        #self.update_oxygen()
    def vessels_killed_by_pressure(self, radius_threshold = 1e-3):
        for vessel in self.vasculature.list_of_vessels:
            if vessel.radius < radius_threshold:
                self.vasculature.kill_vessel(vessel.id)
        return
    def find_voxel_number(self, position):
        #doesn't handle the case where the position is outside the world
        num_voxels = self.number_of_voxels
        i = int((position[0] + self.half_length) * num_voxels / (2 * self.half_length))
        j = int((position[1] + self.half_length) * num_voxels / (2 * self.half_length))
        k = int((position[2] + self.half_length) * num_voxels / (2 * self.half_length))
        n = i * num_voxels ** 2 + j * num_voxels + k
        if n >= self.total_number_of_voxels or n < 0:
            print('position ' +str(position) + ' is outside the world (voxel number = ' + str(n) + ')' )
            return -1
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

    def compute_exchange_matrix(self, dt, pressure_threshold):
        V = self.voxel_list[0].volume
        total_voxels = self.total_number_of_voxels

        # Extract pressure and viscosity values for all voxels
        pressures = np.array([voxel.pressure() for voxel in self.voxel_list])
        viscosities = np.array([voxel.viscosity for voxel in self.voxel_list])

        # Initialize migration matrix with zeros
        migration_matrix = np.zeros((total_voxels, total_voxels))

        for i in range(total_voxels):
            voxel_i = self.voxel_list[i]
            voxel_pressure = pressures[i]
            viscosity = viscosities[i]

            # Find neighbors of the current voxel
            neighbors = self.find_neighbors(voxel_i)

            for neighbor in neighbors:
                j = neighbor.voxel_number

                pressure_diff = voxel_pressure - pressures[j]
                if pressure_diff > 0:
                    t_res = (V / pressure_diff) * viscosity
                    n_events = dt / t_res
                    migration_matrix[i, j] = n_events

        return migration_matrix


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

    def update_oxygen(self, o2_per_volume=10, diffusion_number=5):
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

    def oxygen_map(self):
        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.oxygen)
            positions.append(voxel.position)

        step = self.half_length*2/self.number_of_voxels
        o2_map = ScalarField(positions, values, step)
        return o2_map

    def pressure_map(self):
        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.pressure())
            positions.append(voxel.position)

        step = self.half_length*2/self.number_of_voxels
        print(step)
        pressure_map = ScalarField(positions, values, step)
        print('ok')
        return pressure_map

    def vegf_map(self):
        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.molecular_factors['VEGF'])
            positions.append(voxel.position)

        step = (self.half_length*2/self.number_of_voxels)*4
        vegf_map = ScalarField(positions, values, step)
        return vegf_map

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
            number.append(voxel.number_of_cells())
            positions.append(voxel.position)

        ax.scatter(
        [p[0] for p in positions],
        [p[1] for p in positions],
        [p[2] for p in positions],
        c=number, cmap='YlGn', alpha= 0.5, vmin=0, vmax=2000, s= size)
        # add colorbar
        fig.colorbar(ax.collections[0])
        return fig, ax

    def show_tumor(self, ax, fig, slice = False):
        print('-- Plotting Cell Population')

        #central z slice of world first voxel is at z = 0
        if slice:
            middle_slice = self.number_of_voxels // 2
            first_voxel = middle_slice * self.number_of_voxels**2
            last_voxel = (middle_slice + 1) * self.number_of_voxels**2 -1
            list = self.voxel_list[first_voxel:last_voxel]
        else:
            list = self.voxel_list
            size = 1
        size = 500 / self.number_of_voxels
        number = []
        positions = []
        # collect doses and positions for all voxels
        for voxel in list:
            count = voxel.number_of_tumor_cells()
            if count > 0:
                number.append(count)
                positions.append(voxel.position)

        ax.scatter(
        [p[0] for p in positions],
        [p[1] for p in positions],
        [p[2] for p in positions],
        c=number, cmap='viridis', alpha= 0.8, vmin=30, vmax=500, s= size*3)
        # add colorbar
        fig.colorbar(ax.collections[0])
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

    def show_tumor_surface(self):
        print('-- Plotting Tumor Surface')

        threshold = 20

        voxel_data = np.zeros((self.number_of_voxels, self.number_of_voxels, self.number_of_voxels))
        for i in range(self.number_of_voxels):
            for j in range(self.number_of_voxels):
                for k in range(self.number_of_voxels):
                    voxel = self.voxel_list[i * self.number_of_voxels ** 2 + j * self.number_of_voxels + k]
                    if voxel.number_of_tumor_cells() > threshold:
                        voxel_data[i, j, k] = 1
        print('ok')
        # Label connected regions of the tumor cells
        labels, num_features = ndimage.label(voxel_data)

        # Convert the labeled volume data to a StructuredGrid object
        grid = pv.wrap(labels)

        # Extract the surface using the contour method
        mesh = grid.contour([0.5])

        # Create a PyVista plotter
        plotter = pv.Plotter()

        # Add the mesh to the plotter
        plotter.add_mesh(mesh, cmap="viridis")

        # Set plot title
        plotter.add_title("Tumor Surface")
        #add axis labels
        plotter.add_axes()

        # Show the plot
        #save
        plotter.show(screenshot='tumor_surface.png')

        return

