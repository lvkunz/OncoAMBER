from amber.voxel import *
from amber.vessel import *
from amber.ScalarField import *
from amber.BasicGeometries import *
#np.set_printoptions(threshold=sys.maxsize)
# from scipy.stats import qmc
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import scipy.sparse as sparse
import os


class World:
    def __init__(self, config):
        print('Initializing world')
        print(config.voxel_per_side, config.half_length_world)
        self.half_length = config.half_length_world
        self.voxel_list = []
        self.number_of_voxels = config.voxel_per_side
        self.total_number_of_voxels = self.number_of_voxels ** 3
        self.config = config
        voxel_length = 2 * self.half_length / self.number_of_voxels

        for i in range(self.number_of_voxels):
            for j in range(self.number_of_voxels):
                for k in range(self.number_of_voxels):
                    position = np.array([
                        i * voxel_length - self.half_length + voxel_length / 2,
                        j * voxel_length - self.half_length + voxel_length / 2,
                        k * voxel_length - self.half_length + voxel_length / 2
                    ])
                    self.voxel_list.append(Voxel(position, self.half_length / self.number_of_voxels, viscosity= self.config.viscosity, voxel_number=i * self.number_of_voxels ** 2 + j * self.number_of_voxels + k))
        self.vasculature = VasculatureNetwork(self.config)
        self.o_diameters = []
        self.o_length_values = []
        self.o_bifurcation_values = []
        self.o_VSL_values = []

    def initiate_vasculature(self, list_of_mother_vessels):
        self.vasculature = VasculatureNetwork(self.config, list_of_mother_vessels)
        return

    def vasculature_growth(self, dt, splitting_rate, macro_steps=1, micro_steps=10, weight_direction=0.5, weight_vegf=0.5, weight_pressure=0.5, radius_pressure_sensitive = False):
        print('Vasculature growth')
        pressure_map = self.pressure_map(step_voxel=5)
        def pressure(point): return pressure_map.evaluate(point)
        vegf = self.vegf_map(step_voxel= self.config.vegf_map_step_voxel)

        def vegf_gradient(point): return vegf.gradient(point)

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
        initial_vessel_number = int(initial_vessel_number * 4 * self.half_length ** 2)
        # sampler = qmc.Halton(2, seed=config.seed)
        # points_z = (sampler.random(initial_vessel_number)[:,0] - 0.5) * self.half_length*2
        # points_y = (sampler.random(initial_vessel_number)[:,1] - 0.5) * self.half_length*2
        points_z = np.random.uniform(-self.half_length, self.half_length, initial_vessel_number)
        points_y = np.random.uniform(-self.half_length, self.half_length, initial_vessel_number)
        points_x = -self.half_length * np.ones(initial_vessel_number)
        points = []
        for i in range(len(points_x)):
            points.append([points_x[i], points_y[i], points_z[i]])
        points = np.array(points)
        points2 = []
        for i in range(len(points)):
            if points[i][0] == self.half_length:
                points2.append([points[i][0] - 0.01, points[i][1], points[i][2]])
            else:
                points2.append([points[i][0] + 0.01, points[i][1], points[i][2]])
        points2 = np.array(points2)
        print('Initial vessels: ', points[0], points2[0])

        list_of_vessels = []
        for j in range(len(points)):
            list_of_vessels.append(Vessel([points[j], points2[j]], 0.5, step_size=self.config.vessel_step_size))
        self.initiate_vasculature(list_of_vessels)
        def pressure(point):
            return (self.half_length - abs(point[0]))*0.06
        def vegf_gradient(point):
            return np.array([1,0,0])
            # if point[0] > 0:
            #     return np.array([-point[0],0,0])
            # else:
            #     return np.array([point[0],0,0])

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

        self.vasculature.update_vessels_radius_from_last(self.config.radius_root_vessels, False, pressure)
        for vessel in self.vasculature.list_of_vessels:
            vessel.in_growth = False
            vessel.visible = self.config.visible_original_vessels
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

    def update_volume_occupied_by_vessels(self):
        volume_score = np.zeros(self.total_number_of_voxels)
        bifurcation_score = np.zeros(self.total_number_of_voxels)
        length_score = np.zeros(self.total_number_of_voxels)
        print('-- Computing volume and length occupied by vessels in each voxel')
        for vessel in self.vasculature.list_of_vessels:
            if len(vessel.path) > 0:
                voxel_n_ = self.find_voxel_number(vessel.path[0])
                bifurcation_score[voxel_n_] += 1
            for i in range(1, len(vessel.path)):
                start_point = vessel.path[i - 1]
                end_point = vessel.path[i]
                start_voxel = self.find_voxel_number(start_point)
                if start_voxel == -1:
                    continue
                line_segment = end_point - start_point
                line_segment_length = np.linalg.norm(line_segment)
                length_score[start_voxel] += line_segment_length
                volume_score[start_voxel] += vessel.volume_per_point()
        print('-- Finishing up')
        for voxel in self.voxel_list:
            voxel.vessel_volume = volume_score[voxel.voxel_number]
            voxel.bifurcation_density = bifurcation_score[voxel.voxel_number]
            voxel.vessel_length = length_score[voxel.voxel_number]
        return


    def update_biology_after_RT(self):
        print('-- Updating biology after RT')
        for voxel in self.voxel_list:
            voxel.update_cells_afterRT()
        #self.vessels_killed_by_dose()
        #self.update_oxygen()
    def vessels_killed(self, radius_threshold = 1e-3, length_threshold = 0):
        for vessel in self.vasculature.list_of_vessels:
            if vessel.radius < radius_threshold:
                self.vasculature.kill_vessel(vessel.id)
                print('vessel ' + str(vessel.id) + ' killed by radius')
            #
            # elif vessel.length() < length_threshold:
            #     self.vasculature.kill_vessel(vessel.id)
            #     print('vessel ' + str(vessel.id) + ' killed by length')
        return

    def find_voxel_number(self, position):
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

    def find_voxel(self, position):
        voxel_number = self.find_voxel_number(position)
        return self.voxel_list[voxel_number]

    def find_moor_neighbors(self, voxel):
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
        return neighbors_voxels


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

    def compute_exchange_matrix(self, dt):
        V = self.voxel_list[0].volume
        side = self.voxel_list[0].half_length * 2
        total_voxels = self.total_number_of_voxels
        # Extract pressure and viscosity values for all voxels
        pressures = np.array([voxel.pressure() for voxel in self.voxel_list])
        viscosities = np.array([voxel.viscosity for voxel in self.voxel_list])
        # Initialize migration matrix with zeros
        migration_matrix = sparse.lil_matrix((total_voxels, total_voxels), dtype=np.float32)

        for i in range(total_voxels):
            voxel_i = self.voxel_list[i]
            voxel_pressure = pressures[i]
            viscosity = viscosities[i]

            # Find neighbors of the current voxel
            neighbors = self.find_moor_neighbors(voxel_i)
            for neighbor in neighbors:
                j = neighbor.voxel_number
                pressure_diff = voxel_pressure - pressures[j]
                distance = np.linalg.norm(voxel_i.position - neighbor.position)
                coeff = (side/distance) * dt
                if pressure_diff > 0:
                    t_res = (V / pressure_diff) * viscosity
                    if self.config.verbose:
                        print('V, pressure diff, viscosity ', V, pressure_diff, viscosity)
                        print('t_res = ', t_res)
                    n_events = coeff / t_res
                    migration_matrix[i, j] = n_events

            # pressure pushing cells to move towards the center of the tumor
            voxel_distance = np.linalg.norm(voxel_i.position)
            neighbor_towards_center = self.find_voxel_number(voxel_i.position * (1 - side/voxel_distance))
            migration_matrix[i, neighbor_towards_center] += self.config.pressure_coefficient_central_migration * dt * (voxel_distance**2)
        # Convert the lil_matrix to a csr_matrix for faster arithmetic operations
        migration_matrix = migration_matrix.tocsr()
        #show the matrix
        # plt.spy(migration_matrix)
        # plt.show()
        return migration_matrix

    def topas_param_file(self, name : str):
        print('-- Creating parameter file for topas simulation, file :', name)

        #returns a txt file that can be used as a parameter file for topas
        #the file is saved in the TopasSimulation folder
        repo_path = self.config.working_directory
        print('repo_path = ', repo_path)
        file = open(repo_path + '/' + name + '_geom.txt', 'w')
        print('file = ', file)
        file.write('IncludeFile = '+ repo_path + '/' + name + '.txt \n' +
                                                'd:Ge/MyBox/HLX      = ' + str(self.half_length*2) +' mm \n' +
                                                'd:Ge/MyBox/HLY      = ' + str(self.half_length*2) +' mm \n' +
                                                'd:Ge/MyBox/HLZ      = ' + str(self.half_length*2) +' mm \n' +
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
        print('-- File saved as ' + name + '_geom')
        return name + '_geom'

    def update_dose(self, doses):
        #updates the dose in each voxel
        if len(doses) != self.total_number_of_voxels:
            raise ValueError('Error: the number of doses is not equal to the number of voxels, Probably due to a unmatching Topas simulation')
        for voxel in self.voxel_list:
            voxel.dose = doses[voxel.voxel_number]
        return

    def update_oxygen(self, n_capillaries_per_VVD=1, capillary_length=5):
        print('-- Computing oxygen map')
        side = self.voxel_list[0].half_length * 2
        for voxel in self.voxel_list:
            voxel.oxygen = int((voxel.vessel_volume / side**3) * n_capillaries_per_VVD)
            # voxel.bifurcation_density += voxel.oxygen
        diffusion_number = int(capillary_length / side)
        for i in range(diffusion_number):
            new_oxygen_map = np.zeros(self.total_number_of_voxels)
            print('--- o2 map computing', i, 'out of', diffusion_number)
            for voxel in self.voxel_list:
                sum = voxel.oxygen
                list_neighbors = self.find_moor_neighbors(voxel)
                for neighbor in list_neighbors:
                    sum += neighbor.oxygen
                new_oxygen_map[voxel.voxel_number] = sum / (1 + len(list_neighbors))
            for voxel in self.voxel_list:
                voxel.oxygen = int(new_oxygen_map[voxel.voxel_number])
        return

    def oxygen_map(self):
        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.oxygen)
            positions.append(voxel.position)

        step = self.half_length*2/self.number_of_voxels
        o2_map = ScalarField3D(positions, values, step)
        return o2_map

    def pressure_map(self, step_voxel = 3):
        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.pressure())
            positions.append(voxel.position)

        step = (self.half_length*2/self.number_of_voxels)*step_voxel
        pressure_map = ScalarField3D(positions, values, step)
        return pressure_map

    def vegf_map(self, step_voxel = 3):
        values = []
        positions = []
        for voxel in self.voxel_list:
            values.append(voxel.molecular_factors['VEGF'])
            positions.append(voxel.position)

        step = (self.half_length*2/self.number_of_voxels)*step_voxel
        vegf_map = ScalarField3D(positions, values, step)
        return vegf_map

    def show_tumor_slice(self, ax, fig, voxel_attribute, factor=None, cmap='viridis', vmin=None, vmax=None, norm=None,
                         levels=None, refinement_level=0, extend = 'both'):
        print('-- Plotting Tumor Slice')

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

        if self.config.slice == 'x':
            voxel_list = voxel_list_x
        else:
            voxel_list = voxel_list_z
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

        if self.config.slice == 'x':
            x = np.array([p[0] for p in positions])
            y = np.array([p[1] for p in positions])
            z = np.array(values)
        else:
            x = np.array([p[1] for p in positions])
            y = np.array([p[2] for p in positions])
            z = np.array(values)
        # Create a triangulation of the data
        triangulation = mtri.Triangulation(x, y)

        # Refine the triangulation to increase smoothness
        refiner = mtri.UniformTriRefiner(triangulation)
        triangulation, z = refiner.refine_field(z, subdiv=refinement_level)

        if levels is None:
            contour = ax.tricontourf(triangulation, z, cmap=cmap, alpha=0.7, vmin=vmin, vmax=vmax, norm=norm, extend = extend)
        else:
            contour = ax.tricontourf(triangulation, z, cmap=cmap, alpha=0.7, levels=levels, norm=norm, extend = extend)

        fig.colorbar(contour, ax=ax, shrink=0.5)

        return fig, ax

    def show_tumor_3D(self, ax, fig, voxel_attribute, factor=None, cmap='viridis', vmin=None, vmax=None):
        print('-- Plotting Cell Population 3D')

        #central z slice of world first voxel is at z = 0
        list = self.voxel_list
        size = 750 / self.number_of_voxels
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
        c=values, cmap=cmap, alpha= 0.7, vmin=vmin, vmax=vmax, s=size, marker='h', edgecolors= 'none')
        # add colorbar
        fig.colorbar(ax.collections[0], ax=ax, shrink=0.5)
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

    def is_inside(self, point, cube_half_length = None):
        if cube_half_length == None:
            cube_half_length = self.half_length
        if point[0] < cube_half_length and point[0] > -cube_half_length and point[1] < cube_half_length and point[1] > -cube_half_length and point[2] < cube_half_length and point[2] > -cube_half_length:
            return True
        else:
            return False

    def measure_vasculature_length(self):
        length = 0
        for vessel in self.vasculature.list_of_vessels:
            for point in vessel.path:
                if self.is_inside(point): #add condition for some specific regions
                    length += vessel.step_size
        return length
    def measure_vasculature_area(self):
        area = 0
        for vessel in self.vasculature.list_of_vessels:
            for point in vessel.path:
                if self.is_inside(point):
                    area += 2*np.pi*vessel.radius*vessel.step_size
        return area
    def measure_vasculature_volume(self):
        volume = 0
        for vessel in self.vasculature.list_of_vessels:
            for point in vessel.path:
                if self.is_inside(point):
                    volume += vessel.volume_per_point()
        return volume

    def measure_tumor_volume(self):
        volume_necrotic_free = 0
        volume_total = 0
        for voxel in self.voxel_list:
            if voxel.number_of_tumor_cells() > 0:
                volume_necrotic_free += voxel.volume
            if voxel.number_of_tumor_cells() > 0 or voxel.number_of_necrotic_cells() > 0:
                volume_total += voxel.volume
        return volume_total, volume_necrotic_free

    def show_angiogenesis_metrics(self, t, real=True):
        # Extract the voxel values for each parameter
        volume = self.voxel_list[0].volume
        side = self.voxel_list[0].half_length * 2
        length_values = []
        volume_values = []
        for voxel in self.voxel_list:
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
        plt.show()

        #print a table with the values of the metrics for the current distribution, mean, median, std, min, max
        print('Angiogenesis Metrics')
        print('Diameters')
        print(f'Mean: {diameter_mean:.2f}, Median: {diameter_median:.2f}, Std: {np.std(diameters_values):.2f}, Min: {np.min(diameters_values):.2f}, Max: {np.max(volume_values):.2f}')
        print('Length density')
        print(f'Mean: {length_mean:.2f}, Median: {length_median:.2f}, Std: {np.std(length_values):.2f}, Min: {np.min(length_values):.2f}, Max: {np.max(length_values):.2f}')
        print('Volume density')
        print(f'Mean: {volume_mean:.2f}, Median: {volume_median:.2f}, Std: {np.std(volume_values):.2f}, Min: {np.min(volume_values):.2f}, Max: {np.max(volume_values):.2f}')
        print('Vessel segment length')
        print(f'Mean: {VSL_mean:.2f}, Median: {VSL_median:.2f}, Std: {np.std(VSL_values):.2f}, Min: {np.min(VSL_values):.2f}, Max: {np.max(VSL_values):.2f}')




