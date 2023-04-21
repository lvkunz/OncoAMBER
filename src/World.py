from src.Voxel import *
from src.Vessel import *
from src.ScalarField import *
from src.BasicGeometries import *
#np.set_printoptions(threshold=sys.maxsize)
from scipy.stats import qmc
import matplotlib.tri as mtri
import scipy.sparse as sparse
from src.config_instance import config



class World:
    def __init__(self, half_length, number_of_voxels : int = 20):
        self.half_length = half_length
        self.voxel_list = []
        self.total_number_of_voxels = number_of_voxels ** 3
        voxel_length = 2 * half_length / number_of_voxels

        for i in range(number_of_voxels):
            for j in range(number_of_voxels):
                for k in range(number_of_voxels):
                    position = np.array([
                        i * voxel_length - half_length + voxel_length / 2,
                        j * voxel_length - half_length + voxel_length / 2,
                        k * voxel_length - half_length + voxel_length / 2
                    ])
                    self.voxel_list.append(Voxel(position, half_length / number_of_voxels, voxel_number=i * number_of_voxels ** 2 + j * number_of_voxels + k))
        self.number_of_voxels = number_of_voxels
        self.vasculature = VasculatureNetwork()

    def initiate_vasculature(self, list_of_mother_vessels):
        self.vasculature = VasculatureNetwork(list_of_mother_vessels)
        return

    def vasculature_growth(self, dt, splitting_rate, macro_steps=1, micro_steps=10, weight_direction=0.5, weight_vegf=0.5, weight_pressure=0.5, radius_pressure_sensitive = False):
        print('Vasculature growth')
        pressure_map = self.pressure_map(step_voxel=5)
        def pressure(point): return pressure_map.evaluate(point)
        vegf = self.vegf_map(step_voxel= config.vegf_map_step_voxel)

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
        self.vasculature.update_vessels_radius_from_last(config.radius_root_vessels, radius_pressure_sensitive, pressure)
        return

    def generate_healthy_vasculature(self, initial_vessel_number):
        initial_vessel_number = int(initial_vessel_number * 4 * self.half_length ** 2)
        sampler = qmc.Halton(2, seed=config.seed)
        points_z = (sampler.random(initial_vessel_number)[:,0] - 0.5) * self.half_length*2
        points_y = (sampler.random(initial_vessel_number)[:,1] - 0.5) * self.half_length*2
        points_x = np.append(self.half_length * np.ones(initial_vessel_number//2), -self.half_length * np.ones(initial_vessel_number//2))
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
            list_of_vessels.append(Vessel([points[j], points2[j]], 0.5))
        self.initiate_vasculature(list_of_vessels)
        def pressure(point):
            return (self.half_length - abs(point[0]))
        def vegf_gradient(point):
            if point[0] > 0:
                return np.array([-point[0],0,0])*0.1
            else:
                return np.array([-point[0],0,0])*0.1

        self.vasculature.grow_and_split(
            dt=1,
            splitting_rate=0.3,
            vegf_gradient= vegf_gradient,
            pressure= pressure,
            macro_steps=int(8.5*self.half_length),
            micro_steps=8,
            weight_direction=3.0,
            weight_vegf=1.0,
            weight_pressure=0.0
        )
        self.vasculature.grow_and_split(
            dt=1,
            splitting_rate=0.3,
            vegf_gradient=vegf_gradient,
            pressure=pressure,
            macro_steps=1,
            micro_steps=8,
            weight_direction=5.0,
            weight_vegf=0.1,
            weight_pressure=0.0
        )


        self.vasculature.update_vessels_radius_from_last(config.radius_root_vessels, False, pressure)
        for vessel in self.vasculature.list_of_vessels:
            vessel.in_growth = False
            vessel.visible = config.visible_original_vessels
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
        migration_matrix = sparse.lil_matrix((total_voxels, total_voxels), dtype=np.float32)

        for i in range(total_voxels):
            voxel_i = self.voxel_list[i]
            voxel_pressure = pressures[i]
            viscosity = viscosities[i]

            # Find neighbors of the current voxel
            neighbors = self.find_neighbors(voxel_i)

            for neighbor in neighbors:
                j = neighbor.voxel_number

                pressure_diff = voxel_pressure - pressures[j]
                if pressure_diff > pressure_threshold:
                    t_res = (V / pressure_diff) * viscosity
                    if config.verbose: print('t_res = ', t_res)
                    n_events = dt / t_res
                    migration_matrix[i, j] = n_events

        # Convert the lil_matrix to a csr_matrix for faster arithmetic operations
        migration_matrix = migration_matrix.tocsr()

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

    def update_oxygen(self, o2_per_volume=1, diffusion_number=5):
        print('-- Computing oxygen map')
        for voxel in self.voxel_list:
            voxel.oxygen = int((voxel.vessel_volume * o2_per_volume) / (np.pi * config.effective_vessel_radius**2 * voxel.half_length * 2))
        for i in range(diffusion_number):
            new_oxygen_map = np.zeros(self.total_number_of_voxels)
            print('--- o2 map computing', i, 'out of', diffusion_number)
            for voxel in self.voxel_list:
                sum = voxel.oxygen
                for neighbor in self.find_neighbors(voxel):
                    sum += neighbor.oxygen
                new_oxygen_map[voxel.voxel_number] = sum / (1 + len(self.find_neighbors(voxel)))
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
                         levels=None, refinement_level=0):
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

        if config.slice == 'x':
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

        if config.slice == 'x':
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
            contour = ax.tricontourf(triangulation, z, cmap=cmap, alpha=0.7, vmin=vmin, vmax=vmax, norm=norm)
        else:
            contour = ax.tricontourf(triangulation, z, cmap=cmap, alpha=0.7, levels=levels, norm=norm)

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
        volume = 0
        for voxel in self.voxel_list:
            if voxel.number_of_tumor_cells() > 3:
                volume += voxel.volume
        return volume




