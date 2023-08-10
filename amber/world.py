import numpy as np
from amber.voxel import *
from amber.vessel import *
from amber.ScalarField import *
from amber.BasicGeometries import *
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.sparse as sparse
import pickle
import os
from time import time as current_time


class World: # class that contains the voxels and the vasculature
    def __init__(self, config):
        print('Initializing world')
        print(config.voxel_per_side, config.half_length_world)
        self.half_length = config.half_length_world #half length of the world
        self.voxel_list = [] #list of voxels in the world
        self.number_of_voxels = config.voxel_per_side #number of voxels per side
        self.total_number_of_voxels = self.number_of_voxels ** 3 #total number of voxels
        self.config = config #configuration class
        voxel_length = 2 * self.half_length / self.number_of_voxels #length of the voxels

        for i in range(self.number_of_voxels): #creates the voxels
            for j in range(self.number_of_voxels):
                for k in range(self.number_of_voxels):
                    position = np.array([
                        i * voxel_length - self.half_length + voxel_length / 2,
                        j * voxel_length - self.half_length + voxel_length / 2,
                        k * voxel_length - self.half_length + voxel_length / 2
                    ])
                    self.voxel_list.append(Voxel(position, self.half_length / self.number_of_voxels, viscosity= self.config.viscosity, voxel_number=i * self.number_of_voxels ** 2 + j * self.number_of_voxels + k))
        self.vasculature = VasculatureNetwork(self.config) #creates the vasculature
        self.o_diameters = []; self.o_length_values = []; self.o_bifurcation_values = []; self.o_VSL_values = [] #lists to store the values of the vasculature at beginning of the simulation
        self.center_of_mass = np.array([0.0, 0.0, 0.0]) #center of mass of the tumor

    def initiate_vasculature(self, list_of_mother_vessels): #initiates the vasculature
        self.vasculature = VasculatureNetwork(self.config, list_of_mother_vessels) #creates the vasculature
        for vessel in self.vasculature.list_of_vessels: #for each vessel in the vasculature
            vessel.must_be_updated = True #the vessel must be updated
        return

    def vasculature_growth(self, dt, splitting_rate, macro_steps=1, micro_steps=10, weight_direction=0.5, weight_vegf=0.5, weight_pressure=0.5, radius_pressure_sensitive = False):
        print('Vasculature growth')
        pressure_map = self.pressure_map(step_gradient=5) #creates the pressure map
        def pressure(point): return pressure_map.evaluate(point) #defines the pressure function
        vegf = self.vegf_map(step_gradient= self.config.vegf_map_step_gradient) #creates the vegf map

        def vegf_gradient(point): return vegf.gradient(point) #defines the vegf gradient function

        #grows the vasculature and updates the radius of the vessels
        self.vasculature.grow_and_split(dt = dt,
                                        splitting_rate = splitting_rate,
                                        macro_steps = macro_steps,
                                        micro_steps = micro_steps,
                                        weight_direction = weight_direction,
                                        weight_vegf = weight_vegf,
                                        weight_pressure = weight_pressure,
                                        pressure = pressure,
                                        vegf_gradient = vegf_gradient)
        self.vasculature.update_vessels_radius_from_last(self.config.radius_root_vessels, radius_pressure_sensitive, pressure)
        return

    def generate_healthy_vasculature(self, initial_vessel_number, splitting_rate =0.3, mult_macro_steps=1, micro_steps=8, weight_direction=3.0, weight_vegf=1.0, weight_pressure=0.0, extra_step = True):
        #generates the healthy vasculature at the beginning of the simulation
        initial_vessel_number = int(initial_vessel_number * 4 * self.half_length ** 2) #number of vessels depending on area
        points_z = np.random.uniform(-self.half_length, self.half_length, initial_vessel_number) #z coordinates of the vessels
        points_y = np.random.uniform(-self.half_length, self.half_length, initial_vessel_number) #y coordinates of the vessels
        points_x = -self.half_length * np.ones(initial_vessel_number) #x coordinates of the vessels
        points = [] #list of points
        points_bis = [] #list of points of the other face
        for i in range(len(points_x)):
            points.append([points_x[i], points_y[i], points_z[i]])
            points_bis.append([-points_x[i], points_y[i], points_z[i]])
        points = np.array(points)
        points_bis = np.array(points_bis)
        points2 = []
        points2_bis = []
        for i in range(len(points)):
            points2.append([points[i][0] + 0.01, points[i][1], points[i][2]])
            points2_bis.append([points_bis[i][0] - 0.01, points_bis[i][1], points_bis[i][2]])

        points2 = np.array(points2)
        points2_bis = np.array(points2_bis)
        print('Initial vessels: ', points[0], points2[0])
        print('Initial vessels other face: ', points_bis[0], points2_bis[0])
        list_of_vessels = []
        list_of_vessels_bis = []
        for j in range(len(points)): #creates the vessels on both opposite side of the world
            list_of_vessels.append(Vessel([points[j], points2[j]], radius = 0.01, step_size=self.config.vessel_step_size, parent_id=None, children_ids=None, in_growth=True, intra_radiosensitivity= self.config.vessel_radiosensitivity))
            list_of_vessels_bis.append(Vessel([points_bis[j], points2_bis[j]], radius = 0.01, step_size=self.config.vessel_step_size, parent_id=None, children_ids=None, in_growth=True, intra_radiosensitivity= self.config.vessel_radiosensitivity))

        self.initiate_vasculature(list_of_vessels)
        vasculature_bis = VasculatureNetwork(self.config, list_of_vessels_bis)

        def pressure(point): #defines the pressure function for the healthy vasculature growth
            return (self.half_length - abs(point[0]))*0.06
        def pressure_bis(point): #defines the pressure function for the healthy vasculature growth on the other face
            return (-self.half_length + abs(point[0]))*0.06
        def vegf_gradient(point): #defines the vegf gradient function for the healthy vasculature growth
            return np.array([1,0,0])

        def vegf_gradient_bis(point): #defines the vegf gradient function for the healthy vasculature growth on the other face
            return np.array([-1,0,0])

        #growth of the healthy vasculature on both opposite faces
        vasculature_bis.grow_and_split(
            dt=1,
            splitting_rate=splitting_rate,
            vegf_gradient= vegf_gradient_bis,
            pressure= pressure_bis,
            macro_steps=int(0.9*self.half_length*mult_macro_steps),
            micro_steps=micro_steps,
            weight_direction=weight_direction,
            weight_vegf=weight_vegf,
            weight_pressure=weight_pressure
        )

        self.vasculature.grow_and_split(
            dt=1,
            splitting_rate=splitting_rate,
            vegf_gradient= vegf_gradient,
            pressure= pressure,
            macro_steps=int(0.9*self.half_length*mult_macro_steps),
            micro_steps=micro_steps,
            weight_direction=weight_direction,
            weight_vegf=weight_vegf,
            weight_pressure=weight_pressure
        )

        #merge the two vasculatures
        self.vasculature.add_multiple_vessels(vasculature_bis.list_of_vessels)

        #update the radius of the vessels

        for vessel in self.vasculature.list_of_vessels:
            vessel.must_be_updated = True
            vessel.in_growth = False
            vessel.maturity = 1.0
            vessel.visible = self.config.visible_original_vessels

        self.vasculature.update_vessels_radius_from_last(self.config.radius_root_vessels, True, pressure)

        return

    def update_volume_occupied_by_vessels(self): #updates the volume occupied by the vessels in each voxel
        volume_score = np.zeros(self.total_number_of_voxels) #volume occupied by the vessels in each voxel
        bifurcation_score = np.zeros(self.total_number_of_voxels) #number of bifurcations in each voxel
        length_score = np.zeros(self.total_number_of_voxels) #length of the vessels in each voxel
        print('-- Computing volume and length occupied by vessels in each voxel')
        for vessel in self.vasculature.list_of_vessels: #for each vessel
            if len(vessel.path) > 0: #if the vessel has a path
                voxel_n_ = self.find_voxel_number(vessel.path[0]) #find the voxel number of the first point of the vessel
                bifurcation_score[voxel_n_] += 1 #add 1 to the bifurcation score of the voxel
            for i in range(1, len(vessel.path)): #for each point of the vessel
                start_point = vessel.path[i - 1] #start point of the line segment
                end_point = vessel.path[i] #end point of the line segment
                start_voxel = self.find_voxel_number(start_point) #find the voxel number of the start point
                if start_voxel == -1: #if the voxel number is -1, the point is outside the domain
                    continue
                line_segment = end_point - start_point #line segment between the two points
                line_segment_length = np.linalg.norm(line_segment) #length of the line segment
                length_score[start_voxel] += line_segment_length #add the length of the line segment to the length score of the voxel
                volume_score[start_voxel] += vessel.volume_per_point() #add the volume of the line segment to the volume score of the voxel
        print('-- Finishing up')
        for voxel in self.voxel_list: #for each voxel
            voxel.vessel_volume = volume_score[voxel.voxel_number] #update the volume occupied by the vessels in the voxel
            voxel.bifurcation_density = bifurcation_score[voxel.voxel_number] #update the bifurcation density in the voxel
            voxel.vessel_length = length_score[voxel.voxel_number] #update the length of the vessels in the voxel
        return

    def vessels_killed(self, radius_threshold): #kills the vessels with a radius below a threshold
        for vessel in self.vasculature.list_of_vessels:
            if vessel.radius < radius_threshold:
                self.vasculature.kill_vessel(vessel.id)
                print('vessel ' + str(vessel.id) + ' killed by radius')
        return

    def find_voxel_number(self, position): #finds the voxel number of a position
        # doesn't handle the case where the position is outside the world
        num_voxels = self.number_of_voxels
        voxel_length = 2 * self.half_length / num_voxels
        i = int(round((position[0] + self.half_length - voxel_length / 2) / voxel_length))
        j = int(round((position[1] + self.half_length - voxel_length / 2) / voxel_length))
        k = int(round((position[2] + self.half_length - voxel_length / 2) / voxel_length))
        n = i * num_voxels ** 2 + j * num_voxels + k
        if n >= self.total_number_of_voxels or n < 0:
            # print('position ' +str(position) + ' is outside the world (voxel number = ' + str(n) + ')' )
            return -1
        return n

    def find_voxel(self, position): #finds the voxel of a position
        voxel_number = self.find_voxel_number(position)
        return self.voxel_list[voxel_number]

    def find_moor_neighbors(self, voxel): #finds the Moore's neighbors of a voxel (26 neighbors)
        # finds the Moore's neighbors of a voxel
        voxel_number = voxel.voxel_number
        num_voxels = self.number_of_voxels
        i = voxel_number // (num_voxels ** 2)
        j = (voxel_number // num_voxels) % num_voxels
        k = voxel_number % num_voxels

        neighbors = []
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                for dk in [-1, 0, 1]:
                    if di == 0 and dj == 0 and dk == 0:
                        continue
                    if (0 <= i + di < num_voxels) and (0 <= j + dj < num_voxels) and (0 <= k + dk < num_voxels):
                        neighbors.append((i + di) * num_voxels ** 2 + (j + dj) * num_voxels + (k + dk))

        neighbors_voxels = [self.voxel_list[n] for n in neighbors]
        # print('Moore neighbors found, number of neighbors = ' + str(len(neighbors_voxels)) + ' (should be 26)')
        return neighbors_voxels


    def find_neighbors(self, voxel): #finds the direct neighbors of a voxel (6 neighbors)
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


    #this function needs to be changed into something that is consistant for different voxel sizes and different dt
    #Use PhysiCell to find the dependencies?
    def compute_exchange_matrix(self, dt): #computes the exchange matrix
        print('-- Computing exchange matrix')
        start = current_time()
        V = self.voxel_list[0].volume #volume of a voxel
        side = self.voxel_list[0].half_length * 2 #length of a voxel
        total_voxels = self.total_number_of_voxels #total number of voxels
        # Extract pressure and viscosity values for all voxels
        pressures = np.array([voxel.pressure() for voxel in self.voxel_list]) #pressure of each voxel
        viscosities = np.array([voxel.viscosity for voxel in self.voxel_list]) #viscosity of each voxel
        # Initialize migration matrix with zeros
        migration_matrix = sparse.lil_matrix((total_voxels, total_voxels), dtype=np.float32) #migration matrix

        #find center of mass of the tumor
        center_of_mass = np.array([0.0,0.0,0.0]) #center of mass
        total_cells = 0
        for voxel in self.voxel_list:
            cells = voxel.number_of_tumor_cells()
            center_of_mass += voxel.position * cells
            total_cells += voxel.number_of_tumor_cells()
        if total_cells < 1:
            center_of_mass = np.array([0.0,0.0,0.0])
        else:
            center_of_mass /= total_cells
        print('center of mass real', center_of_mass) #center of mass
        #we don't let the tumor move too far from the center
        if np.linalg.norm(center_of_mass) > self.half_length/3:
            center_of_mass = center_of_mass / np.linalg.norm(center_of_mass) * self.half_length/3 #center of mass
        print('center of mass moved', center_of_mass) #new center of mass
        self.center_of_mass = center_of_mass

        for i in range(total_voxels): #for each voxel
            voxel_i = self.voxel_list[i]
            voxel_pressure = pressures[i]
            viscosity = viscosities[i]

            # Find neighbors of the current voxel
            neighbors = self.find_moor_neighbors(voxel_i)
            for neighbor in neighbors:
                j = neighbor.voxel_number
                pressure_diff = voxel_pressure - pressures[j] #pressure difference between the current voxel and its neighbor
                distance = np.linalg.norm(voxel_i.position - neighbor.position) #distance between the current voxel and its neighbor
                coeff = (side/distance) * dt
                if pressure_diff > 0: #if the pressure of the current voxel is higher than the pressure of the neighbor
                    t_res = (V / pressure_diff) * viscosity
                    if self.config.verbose:
                        print('V, pressure diff, viscosity ', V, pressure_diff, viscosity)
                        print('t_res = ', t_res)
                    n_events = coeff / t_res
                    migration_matrix[i, j] = n_events

            # pressure pushing cells to move towards the center of the tumor\

            vector_to_center = voxel_i.position - center_of_mass #vector from the center of the tumor to the current voxel
            distance = np.linalg.norm(vector_to_center) #distance between the center of the tumor and the current voxel
            if distance > 0:
                neighbor_towards_center = self.find_voxel_number(voxel_i.position - (side/distance) * vector_to_center)
                migration_matrix[i, neighbor_towards_center] += self.config.pressure_coefficient_central_migration * dt * (distance**2)

        # Convert the lil_matrix to a csr_matrix for faster arithmetic operations
        migration_matrix = migration_matrix.tocsr()
        #show the matrix
        # plt.figure(figsize=(20, 20))
        # #only show the first 100 voxels
        # matrix = migration_matrix.toarray()
        # matrix = matrix[:500, :500]
        # plt.imshow(matrix, cmap='inferno')
        # plt.colorbar()
        # plt.title('Migration Matrix')
        # plt.xlabel('Destination Voxel')
        # plt.ylabel('Source Voxel')
        # plt.tight_layout()
        # plt.show()
        end = current_time()
        print('Time to compute exchange matrix : ', end - start)
        return migration_matrix

    def topas_param_file(self, name : str):
        print('-- Creating parameter file for topas simulation, file :', name)

        #returns a txt file that can be used as a parameter file for topas. Gets most of the parameters from TOPAS file specified in the config file
        #the tumor geometry is added to the file
        repo_path = self.config.working_directory
        print('repo_path = ', repo_path)
        file = open(repo_path + '/' + name + '_geom.txt', 'w')
        print('file = ', file)
        file.write('IncludeFile = '+ repo_path + '/' + name + '.txt \n' +
                                                'd:Ge/Tumor/HLX      = ' + str(self.half_length*2) +' mm \n' +
                                                'd:Ge/Tumor/HLY      = ' + str(self.half_length*2) +' mm \n' +
                                                'd:Ge/Tumor/HLZ      = ' + str(self.half_length*2) +' mm \n' +
                                                'd:Ge/Tumor/TransX   = 0. m \n' +
                                                'd:Ge/Tumor/TransY   = 0. m \n' +
                                                'd:Ge/Tumor/TransZ   = 0. m \n' +
                                                'd:Ge/Tumor/RotX     = 0. deg \n' +
                                                'd:Ge/Tumor/RotY     = 0. deg \n' +
                                                'd:Ge/Tumor/RotZ     = 0. deg \n' +
                                                'i:Ge/Tumor/ZBins = ' + str(self.number_of_voxels) +' \n' +
                                                'i:Ge/Tumor/XBins = ' + str(self.number_of_voxels) +' \n' +
                                                'i:Ge/Tumor/YBins = ' + str(self.number_of_voxels))
        file.close()
        print('-- File saved as ' + name + '_geom')
        return name + '_geom'

    def update_dose(self, doses): #updates the dose in each voxel
        #updates the dose in each voxel
        if len(doses) != self.total_number_of_voxels:
            raise ValueError('Error: the number of doses is not equal to the number of voxels, Probably due to a unmatching Topas simulation')
        for voxel in self.voxel_list:
            voxel.dose = doses[voxel.voxel_number]
        return

    def update_capillaries(self, n_capillaries_per_VVD=1, capillary_length=5): #updates the number of capillaries in each voxel
        print('-- Computing capillaries map')
        side = self.voxel_list[0].half_length * 2
        for voxel in self.voxel_list:
            voxel.n_capillaries = int((voxel.vessel_volume / side ** 3) * n_capillaries_per_VVD) #vessel volume density in voxel * number of capillaries per VVD
            # voxel.bifurcation_density += voxel.capillaries
        diffusion_number = int(capillary_length / side) #number of diffusion steps to reach the capillary length
        for i in range(diffusion_number): # "diffusion" of the capillaries
            new_oxygen_map = np.zeros(self.total_number_of_voxels)
            print('--- o2 map computing', i, 'out of', diffusion_number)
            for voxel in self.voxel_list: #for each voxel, compute the average number of capillaries in neighbors
                sum = voxel.n_capillaries  # the voxel itself is taken into account
                list_neighbors = self.find_moor_neighbors(voxel) #list of the 6 neighbors of the voxel
                for neighbor in list_neighbors:
                    sum += neighbor.n_capillaries
                new_oxygen_map[voxel.voxel_number] = sum / (1 + len(list_neighbors))
            for voxel in self.voxel_list:
                voxel.n_capillaries = int(new_oxygen_map[voxel.voxel_number])
        return

    def oxygen_map(self):
        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.n_capillaries)
            positions.append(voxel.position)

        step = self.half_length*2/self.number_of_voxels
        o2_map = ScalarField3D(positions, values, step) #creates a scalar field with the oxygen values
        return o2_map

    def pressure_map(self, step_gradient = 3): #creates a scalar field with the pressure values
        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.pressure())
            positions.append(voxel.position)

        step = step_gradient
        pressure_map = ScalarField3D(positions, values, step)
        return pressure_map

    def vegf_map(self, step_gradient = 3): #creates a scalar field with the vegf values

        diffusion_number = 3  # number of diffusion steps to reach the capillary length
        for i in range(diffusion_number):  # "diffusion" of the capillaries
            new_vegf_map = np.zeros(self.total_number_of_voxels)
            print('--- vegf diffusion computing', i, 'out of', diffusion_number)
            for voxel in self.voxel_list:  # for each voxel, compute the average number of capillaries in neighbors
                sum = voxel.molecular_factors['VEGF'] # the voxel itself is taken into account
                list_neighbors = self.find_moor_neighbors(voxel)  # list of the 6 neighbors of the voxel
                for neighbor in list_neighbors:
                    sum += neighbor.molecular_factors['VEGF']
                new_vegf_map[voxel.voxel_number] = sum / (1 + len(list_neighbors))
            for voxel in self.voxel_list:
                voxel.molecular_factors['VEGF'] = new_vegf_map[voxel.voxel_number]


        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.molecular_factors['VEGF'])
            positions.append(voxel.position)

        step = step_gradient
        vegf_map = ScalarField3D(positions, values, step)
        return vegf_map

    def dose_map(self): #creates a scalar field with the dose values
        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.dose)
            positions.append(voxel.position)

        step = 1
        dose_map = ScalarField3D(positions, values, step)
        return dose_map

    def vegf_total(self):
        total = 0
        for voxel in self.voxel_list:
            total += voxel.molecular_factors['VEGF']
        return total

    def show_tumor_slice(self, ax, fig, voxel_attribute, factor=None, cmap='viridis', vmin=None, vmax=None, norm=None,
                         levels=None, refinement_level=0, extend = 'both', slice = None, round_n = None): #plots a slice of the tumor with the voxel_attribute as color
        print('-- Plotting Tumor Slice')

        if slice == None:
            slice = self.config.slice


        middle_slice = self.number_of_voxels // 2
        first_voxel = middle_slice * self.number_of_voxels ** 2
        last_voxel = (middle_slice + 1) * self.number_of_voxels ** 2
        voxel_list_z = self.voxel_list[first_voxel:last_voxel]

        middle_slice_x = self.number_of_voxels // 2
        voxel_list_x = []
        for i in range(self.number_of_voxels):
            for j in range(self.number_of_voxels):
                index = middle_slice_x + j * self.number_of_voxels + i * self.number_of_voxels ** 2
                voxel_list_x.append(self.voxel_list[index])

        middle_slice_y = self.number_of_voxels // 2
        voxel_list_y = []
        for i in range(self.number_of_voxels):
            for j in range(self.number_of_voxels):
                index = j + middle_slice_y * self.number_of_voxels + i * self.number_of_voxels ** 2
                voxel_list_y.append(self.voxel_list[index])

        if slice == 'x':
            voxel_list = voxel_list_x
        elif slice == 'z':
            voxel_list = voxel_list_z
        elif slice == 'y':
            voxel_list = voxel_list_y
        else:
            raise ValueError('Slice must be x, y or z')

        values = []
        positions = []
        for voxel in voxel_list:
            if factor is not None and isinstance(getattr(voxel, voxel_attribute), dict):
                value = getattr(voxel, voxel_attribute)[factor]
            else:
                value = getattr(voxel, voxel_attribute)() if callable(getattr(voxel, voxel_attribute)) else getattr(
                    voxel, voxel_attribute)
            values.append(value)
            positions.append(voxel.position)

        if slice == 'x':
            x = np.array([p[0] for p in positions])
            y = np.array([p[1] for p in positions])
            z = np.array(values)
        elif slice == 'z':
            x = np.array([p[1] for p in positions])
            y = np.array([p[2] for p in positions])
            z = np.array(values)
        elif slice == 'y':
            x = np.array([p[0] for p in positions])
            y = np.array([p[2] for p in positions])
            z = np.array(values)
        # Create a triangulation of the data
        triangulation = mtri.Triangulation(x, y)

        # Refine the triangulation to increase smoothness
        refiner = mtri.UniformTriRefiner(triangulation)
        triangulation, z = refiner.refine_field(z, subdiv=refinement_level)

        if levels is None:
            if vmin is None: vmin = np.min(z)
            if vmax is None: vmax = np.max(z)
            ticks = np.linspace(vmin, vmax, 11)
            if round_n is not None: ticks = np.round(ticks, round_n)
            contour = ax.tricontourf(triangulation, z, cmap=cmap, alpha=0.7, vmin=vmin, vmax=vmax, norm=norm, extend = extend)
        else:
            ticks = np.linspace(levels[0], levels[-1], 11)
            if round_n is not None: ticks = np.round(ticks, round_n)
            contour = ax.tricontourf(triangulation, z, cmap=cmap, alpha=0.7, levels=levels, norm=norm, extend = extend)

        # #add a line for scale bar in the bottom right corner
        # hf = self.half_length
        # bar_length = self.half_length/self.number_of_voxels * 5
        # text = int(bar_length * 1000)
        #
        # ax.plot([-hf, -hf], [hf-bar_length, hf], color='black', linewidth=1)
        # ax.text(hf - bar_length/2, hf - bar_length/2, str(text) + ' $\mu$m', fontsize=10)


        #add a point at the center of mass
        if slice == 'x':
            ax.scatter(self.center_of_mass[0], self.center_of_mass[1], color='black', s=10)
        elif slice == 'z':
            ax.scatter(self.center_of_mass[1], self.center_of_mass[2], color='black', s=10)
        elif slice == 'y':
            ax.scatter(self.center_of_mass[0], self.center_of_mass[2], color='black', s=10)


        fig.colorbar(contour, ax=ax, shrink=0.5, ticks = ticks)

        return fig, ax

    def show_tumor_3D(self, ax, fig, voxel_attribute, factor=None, cmap='viridis', vmin=None, vmax=None): #plots the tumor in 3D with the voxel_attribute as color
        print('-- Plotting Cell Population 3D')

        #central z slice of world first voxel is at z = 0
        list = self.voxel_list
        size = 300
        values = []
        positions = []
        # collect doses and positions for all voxels
        for voxel in list:
            if factor is not None and isinstance(getattr(voxel, voxel_attribute), dict):
                value = getattr(voxel, voxel_attribute)[factor]
            else:
                value = getattr(voxel, voxel_attribute)() if callable(getattr(voxel, voxel_attribute)) else getattr(
                    voxel, voxel_attribute)
            if value > 0:
                values.append(value)
                positions.append(voxel.position)

        ax.scatter(
        [p[0] for p in positions],
        [p[1] for p in positions],
        [p[2] for p in positions],
        c=values, cmap=cmap, alpha= 0.7, vmin=vmin, vmax=vmax, s=size, marker='o', edgecolors= 'none')
        # add colorbar
        fig.colorbar(ax.collections[0], ax=ax, shrink=0.5)
        fig.tight_layout()
        return fig, ax

    # def show_tumor_surface(self):
    #     print('-- Plotting Tumor Surface')
    #
    #     threshold = 20
    #     voxel_data = np.zeros((self.number_of_voxels, self.number_of_voxels, self.number_of_voxels))
    #     for i in range(self.number_of_voxels):
    #         for j in range(self.number_of_voxels):
    #             for k in range(self.number_of_voxels):
    #                 voxel = self.voxel_list[i * self.number_of_voxels ** 2 + j * self.number_of_voxels + k]
    #                 if voxel.number_of_tumor_cells() > threshold:
    #                     voxel_data[i, j, k] = 1
    #     # Label connected regions of the tumor cells
    #     labels, num_features = ndimage.label(voxel_data)
    #     grid = pv.wrap(labels)
    #     mesh = grid.contour([0.5])
    #     plotter = pv.Plotter()
    #     plotter.add_mesh(mesh, cmap="viridis")
    #     plotter.add_title("Tumor Surface")
    #     plotter.add_axes()
    #     plotter.show(screenshot='tumor_surface.png')
    #
    #     return

    def is_inside(self, point, cube_half_length = None): #checks if a point is inside the tumor
        if cube_half_length == None:
            cube_half_length = self.half_length
        if point[0] < cube_half_length and point[0] > -cube_half_length and point[1] < cube_half_length and point[1] > -cube_half_length and point[2] < cube_half_length and point[2] > -cube_half_length:
            return True
        else:
            return False

    def measure_vasculature_length(self): #measures the length of the vasculature inside the tumor
        length = 0
        for vessel in self.vasculature.list_of_vessels:
            for point in vessel.path:
                if self.is_inside(point): #add condition for some specific regions
                    length += vessel.step_size
        return length
    def measure_vasculature_area(self): #measures the area of the vasculature inside the tumor
        area = 0
        for vessel in self.vasculature.list_of_vessels:
            for point in vessel.path:
                if self.is_inside(point):
                    area += 2*np.pi*vessel.radius*vessel.step_size
        return area
    def measure_vasculature_volume(self): #measures the volume of the vasculature inside the tumor
        volume = 0
        for vessel in self.vasculature.list_of_vessels:
            for point in vessel.path:
                if self.is_inside(point):
                    volume += vessel.volume_per_point()
        return volume

    def measure_tumor_volume(self): #measures the volume of the tumor
        volume_necrotic_free = 0
        volume_total = 0
        for voxel in self.voxel_list:
            if voxel.number_of_tumor_cells() > 50: #threshold to be considered tumor
                volume_necrotic_free += voxel.volume
            if voxel.number_of_tumor_cells() > 50 or voxel.number_of_necrotic_cells() > 50:
                volume_total += voxel.volume
        return volume_total, volume_necrotic_free

    def show_angiogenesis_metrics(self, t, real=True): #plots the angiogenesis metrics
        # Extract the voxel values for each parameter
        volume = self.voxel_list[0].volume
        side = self.voxel_list[0].half_length * 2
        length_values = []
        volume_values = []
        voxel_list_iter = []

        if self.config.visible_original_vessels:
            voxel_list_iter = self.voxel_list
        else:
            for voxel in self.voxel_list:
                if voxel.number_of_tumor_cells() > 0:
                    voxel_list_iter.append(voxel)

        for voxel in voxel_list_iter:
            if real:
                vessel_volume_density = (voxel.vessel_volume) * 100 / volume
                vessel_length_density = (voxel.vessel_length) / volume
            else:
                vessel_volume_density = voxel.vessel_volume_density()
                vessel_length_density = voxel.vessel_length_density()
            volume_values.append(vessel_volume_density)
            length_values.append(vessel_length_density)

        VSL_values = self.vasculature.compute_VSL()
        diameters_values = self.vasculature.compute_diameters() * 1000

        # save values
        if t == 0:
            self.o_diameters = diameters_values
            self.o_length_values = length_values
            self.o_bifurcation_values = volume_values
            self.o_VSL_values = VSL_values

        # Compute statistics for current distribution
        diameter_mean = np.mean(diameters_values)
        diameter_median = np.median(diameters_values)
        length_mean = np.mean(length_values)
        length_median = np.median(length_values)
        volume_mean = np.mean(volume_values)
        volume_median = np.median(volume_values)
        VSL_mean = np.mean(VSL_values)
        VSL_median = np.median(VSL_values)

        # Plot histograms for each parameter
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))

        # Plot histogram for original distribution
        axes[0, 0].hist(self.o_diameters, bins=20, color='gray', alpha=0.5, label='Original',density = True)
        axes[0, 1].hist(self.o_length_values, bins=20, color='gray', alpha=0.5, label='Original',density = True)
        axes[1, 0].hist(self.o_bifurcation_values, bins=20, color='gray', alpha=0.5, label='Original',density = True)
        axes[1, 1].hist(self.o_VSL_values, bins=20, color='gray', alpha=0.5, label='Original',density = True)

        # Plot histogram for current distribution
        axes[0, 0].hist(diameters_values, bins=20, color='blue', alpha=0.5, label='Current',density = True)
        axes[0, 1].hist(length_values, bins=20, color='green', alpha=0.5, label='Current',density = True)
        axes[1, 0].hist(volume_values, bins=20, color='red', alpha=0.5, label='Current',density = True)
        axes[1, 1].hist(VSL_values, bins=20, color='purple', alpha=0.5, label='Current',density = True)

        # Add titles and legend to the histograms
        axes[0, 0].set_title(f'Vessel diameters [um]\nMean: {diameter_mean:.2f}, Median: {diameter_median:.2f}')
        axes[0, 0].legend()
        axes[0, 1].set_title(f'Vessel length density [mm/mm^3]\nMean: {length_mean:.2f}, Median: {length_median:.2f}')
        axes[1, 0].legend()
        axes[1, 0].set_title(f'Vessel Volume density [%] \nMean: {volume_mean:.2f}, Median: {volume_median:.2f}')
        axes[0, 1].legend()
        axes[1, 1].set_title(f'Vessel segment length [mm] \nMean: {VSL_mean:.2f}, Median: {VSL_median:.2f}')
        axes[1, 1].legend()

        fig.suptitle('Angiogenesis Metrics at time ' + str(t) + ' hours, number of vessels = ' + str(len(self.vasculature.list_of_vessels)))
        fig.tight_layout()
        # Show the plot
        plt.savefig('Plots/CurrentPlotting/t' + str(t) +'_Angiogenesis.png')
        if self.config.running_on_cluster:
            plt.close()
        else:
            plt.show()

        #print a table with the values of the metrics for the current distribution, mean, median, std, min, max
        print('Angiogenesis Metrics')

        length_values = np.array(length_values)
        volume_values = np.array(volume_values)
        VSL_values = np.array(VSL_values)
        diameters_values = np.array(diameters_values)

        if diameters_values.any():
            diameter_std = np.std(diameters_values)
            diameter_min = np.min(diameters_values)
            diameter_max = np.max(diameters_values)
            print('Diameters')
            print(
                f'Mean: {diameter_mean:.2f}, Median: {diameter_median:.2f}, Std: {diameter_std:.2f}, Min: {diameter_min:.2f}, Max: {diameter_max:.2f}')
        else:
            print('Diameters: No data available.')

        if length_values.any():
            length_std = np.std(length_values)
            length_min = np.min(length_values)
            length_max = np.max(length_values)
            print('Length density')
            print(
                f'Mean: {length_mean:.2f}, Median: {length_median:.2f}, Std: {length_std:.2f}, Min: {length_min:.2f}, Max: {length_max:.2f}')
        else:
            print('Length density: No data available.')

        if volume_values.any():
            volume_std = np.std(volume_values)
            volume_min = np.min(volume_values)
            volume_max = np.max(volume_values)
            print('Volume density')
            print(
                f'Mean: {volume_mean:.2f}, Median: {volume_median:.2f}, Std: {volume_std:.2f}, Min: {volume_min:.2f}, Max: {volume_max:.2f}')
        else:
            print('Volume density: No data available.')

        if VSL_values.any():
            VSL_std = np.std(VSL_values)
            VSL_min = np.min(VSL_values)
            VSL_max = np.max(VSL_values)
            print('Vessel segment length')
            print(
                f'Mean: {VSL_mean:.2f}, Median: {VSL_median:.2f}, Std: {VSL_std:.2f}, Min: {VSL_min:.2f}, Max: {VSL_max:.2f}')
        else:
            print('Vessel segment length: No data available.')

    def save(self, path): #save the world object to a file
        print('Saving world object to file: ' + path)
        # Save the tumor object to a file
        with open(path, 'wb') as f:
            pickle.dump(self, f)

def load(path): #load the world object from a file
    # Load the tumor object from a file
    with open(path, 'rb') as f:
        world = pickle.load(f)
    return world



