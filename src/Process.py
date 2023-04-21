from src.World import *
from src.Voxel import *
from src.ScalarField import *
import src.Terminal as term
import src.ReadAndWrite as rw
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from src.config_instance import config
import matplotlib.pyplot as plt
import os

class Simulator: #this class is used to run the whole simulation

    def __init__(self, list_of_process : list, finish_time, dt):
        self.list_of_process = list_of_process
        self.finish_time = finish_time
        self.dt = dt
        self.time = 0
    def show_cell_and_tumor_volume(self, number_tumor_cells, number_necrotic_cells, tumor_size, times):
        # plot number of cells evolution
        plt.plot(times, number_tumor_cells, 'purple')
        plt.plot(times, number_necrotic_cells, 'black')
        plt.title('Number of cells evolution')
        plt.xlabel('Time')
        plt.ylabel('Number of cells')
        plt.grid(True)
        plt.savefig('Plots/Number_cells_evolution.png')
        plt.show()

        # plot tumor size evolution
        fig = plt.figure()
        plt.plot(times, tumor_size, 'red')
        plt.title('Tumor volume evolution')
        plt.xlabel('Time')
        plt.ylabel('Tumor volume [mm^3]')
        plt.grid(True)
        plt.savefig('Plots/Tumor_size_evolution.png')
        plt.show()
    def show(self, world: World, t = 0): #this function is used to show the world at a certain time

        DPI = 100
        size = world.half_length

        #plot vasculature
        if config.show_tumor_and_vessels_3D:
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(25, 25), dpi=150, subplot_kw={'projection': '3d'})
            fig.suptitle('Visualization at time t = ' + str(t) + ' hours', fontsize=16)
            axes.view_init(0, 90)
            size = size/2
            axes.set_xlim(-size, size)
            axes.set_ylim(-size, size)
            axes.set_zlim(-size, size)
            world.show_tumor_3D(axes, fig, 'number_of_tumor_cells', cmap='viridis', vmin=0, vmax=1000)
            world.vasculature.plot(fig, axes)
            axes.set_title('Vasculature')
            plt.savefig('Plots/CurrentPlotting/t' + str(t) + '_Vasculature.png')
            plt.show()
            size = size*2

        if config.show_slices:
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(25, 20), dpi=DPI)
            fig.suptitle('Visualization at time t = ' + str(t) + ' hours', fontsize=16)

            axes[0, 0].set_xlim(-size, size)
            axes[0, 0].set_ylim(-size, size)
            world.show_tumor_slice(axes[0, 0], fig, 'number_of_tumor_cells', levels= np.linspace(1, 1001, 11))
            axes[0,0].grid(True)
            axes[0,0].set_facecolor('whitesmoke')
            axes[0, 0].set_title('Cells in voxels')

            norm = TwoSlopeNorm(vmin=0, vcenter=15, vmax=100)

            axes[0, 1].set_xlim(-size, size)
            axes[0, 1].set_ylim(-size, size)
            world.show_tumor_slice(axes[0, 1], fig, 'oxygen', cmap = 'RdBu', norm = norm, levels= np.linspace(0, 110, 16))
            axes[0, 1].grid(True)
            axes[0, 1].set_facecolor('whitesmoke')
            axes[0, 1].set_title('Oxygen in voxels')

            axes[1, 0].set_xlim(-size, size)
            axes[1, 0].set_ylim(-size, size)
            world.show_tumor_slice(axes[1, 0], fig, 'molecular_factors', factor='VEGF', levels= np.linspace(0.001, 1.0, 11), cmap='Oranges')
            axes[1, 0].grid(True)
            axes[1, 0].set_facecolor('whitesmoke')
            axes[1, 0].set_title('VEGF in voxels')

            axes[1, 1].set_xlim(-size, size)
            axes[1, 1].set_ylim(-size, size)
            world.show_tumor_slice(axes[1, 1], fig, 'number_of_necrotic_cells', levels= np.linspace(1, 1001, 11))
            axes[1, 1].grid(True)
            axes[1, 1].set_facecolor('whitesmoke')
            axes[1, 1].set_title('Necrosis in voxels')

            plt.tight_layout()
            plt.savefig('Plots/CurrentPlotting/t' + str(t) + '_AllPlots.png')
            plt.show()

        if config.show_o2_vitality_histograms:
            voxels_positions = [[0,0,0],[0.06, 0.06, 0.0], [0.12,0.12,0.0]]
            fig, axes = plt.subplots(nrows=2, ncols=len(voxels_positions), figsize=(20, 10), dpi=100)
            fig.suptitle('Visualization at time t = ' + str(t) + ' hours', fontsize=16)
            for i in range(len(voxels_positions)):
                #show histograms for the three voxels
                axes[0, i].set_title('Voxel ' + str(i))
                axes[0, i].set_xlabel('Oxygen')
                axes[0, i].set_ylabel('Number of cells')
                voxel = world.find_voxel(voxels_positions[i])
                voxel.cycling_time_and_age_histogram(axes[0, i], fig)
                axes[1, i].set_xlabel('Vitality')
                axes[1, i].set_ylabel('Number of cells')
                voxel.vitality_histogram(axes[1, i], fig)
            plt.show()

    def run(self, world: World, video=False):
        print(f'Running simulation for {self.finish_time} hours with dt={self.dt}')
        process_local = [process for process in self.list_of_process if not process.is_global]
        process_global = [process for process in self.list_of_process if process.is_global]

        irradiations_times = [config.first_irradiation_time + i * config.time_between_fractions for i in
                              range(config.number_fractions)]
        applied_fractions = 0

        number_tumor_cells = []; number_necrotic_cells = []; tumor_size = []; times = []

        while self.time < self.finish_time:
            print(f'\033[1;31;47mTime: {self.time} hours / {self.finish_time} hours\033[0m')

            if applied_fractions < config.number_fractions and self.time >= irradiations_times[applied_fractions]:
                irrad = Irradiation('irrad', self.dt, config.topas_file, config.first_irradiation_time,
                                    config.irradiation_intensity, world)
                irrad(world)
                applied_fractions += 1

            for voxel in world.voxel_list:
                for process in process_local:
                    process(voxel)
            for process in process_global:
                print('Currently running global process:', process.name)
                process(world)

            if video:
                self.show(world, self.time)

            number_tumor_cells.append(sum([voxel.number_of_tumor_cells() for voxel in world.voxel_list]))
            number_necrotic_cells.append(sum([voxel.number_of_necrotic_cells() for voxel in world.voxel_list]))
            tumor_size_ = world.measure_tumor_volume()
            tumor_size.append(tumor_size_ * 1000)
            times.append(self.time)

            np.save('DataOutput/number_tumor_cells.npy', number_tumor_cells)
            np.save('DataOutput/number_necrotic_cells.npy', number_necrotic_cells)
            np.save('DataOutput/tumor_size.npy', tumor_size)
            np.save('DataOutput/times.npy', times)

            if config.show_cell_and_tumor_volume:
                self.show_cell_and_tumor_volume(number_tumor_cells, number_necrotic_cells, tumor_size, times)

            self.time += self.dt

        print('Simulation finished')
        if config.show_final:
            self.show(world, self.time)

        return

class Process: #abstract class, represents all the processes that can happen in the simulation
    def __init__(self, name, dt):
        self.name = name
        self.dt = dt
        self.is_global = False
    def __call__(self, voxel):
        pass


class CellDivision(Process): #cell division process, cells divide in a voxel if they have enough vitality
    def __init__(self, name, dt, cycling_threshold, pressure_threshold = np.inf):
        super().__init__(name, dt)
        self.dt = dt
        self.cycling_threshold = cycling_threshold
        self.pressure_threshold = pressure_threshold

    def __call__(self, voxel):
        if len(voxel.list_of_cells) > 0:
            for cell in voxel.list_of_cells:
                if cell.time_spent_cycling >= cell.doubling_time:
                    time_diff = cell.time_spent_cycling - cell.doubling_time
                    leftover_time = self.dt - time_diff
                    new_cell = cell.duplicate() #create a new cell (start cycling at 0)
                    new_cell.time_spent_cycling = leftover_time #reset the time spent cycling
                    cell.doubling_time = cell.random_doubling_time() #sample a new doubling time for the old cell
                    cell.time_spent_cycling = leftover_time #reset the time spent cycling
                    voxel.add_cell(new_cell) #add the new cell to the voxel
        return

class CellDeath(Process): #cell necrosis process, cells die in a voxel if they have too low vitality
    def __init__(self, name, dt, necrosis_threshold, necrosis_probability, apoptosis_threshold, apoptosis_probability):
        super().__init__('CellNecrosis', dt)
        self.necrosis_threshold = necrosis_threshold
        self.necrosis_probability = necrosis_probability
        self.apoptosis_threshold = apoptosis_threshold
        self.apoptosis_probability = apoptosis_probability
        if self.necrosis_threshold > self.apoptosis_threshold:
            raise ValueError('necrosis threshold must be smaller or equal to apoptosis threshold. you can set apoptosis probability to 0 if you want to avoid apoptosis.')
    def necrosis_curve(self, x):
        r = self.necrosis_probability - x / self.necrosis_threshold
        if r < 0: r = 0
        return r

    def apoptosis_curve(self, x):
        r_ = (self.apoptosis_probability / self.apoptosis_threshold) * x
        if r_ < 0: r_ = 0
        return r_

    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            vitality = cell.vitality()
            if vitality < self.apoptosis_threshold:
                sample = np.random.uniform(0, 1)
                n = self.necrosis_curve(vitality)
                a = self.apoptosis_curve(vitality) + n
                if sample < n:
                    #necrosis
                    voxel.cell_becomes_necrotic(cell)
                elif sample < a:
                    #apoptosis
                    voxel.remove_cell(cell)


class CellAging(Process): #cell aging process, cells age in a voxel
    def __init__(self, name, dt):
        super().__init__('CellAging', dt)
    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            if cell.time_before_death is not None:
                cell.time_before_death -= self.dt
                if cell.time_before_death < 0:
                    voxel.remove_cell(cell)
            if cell.type == 'TumorCell' and cell.vitality() > config.vitality_cycling_threshold:
                cell.time_spent_cycling += self.dt

        for n_cell in voxel.list_of_necrotic_cells:
            if n_cell.time_before_death is not None:
                n_cell.time_before_death -= self.dt
                if n_cell.time_before_death < 0:
                    voxel.remove_necrotic_cell(n_cell)
        pass


class CellMigration(Process): #cell migration process, cells migrate in the world
    def __init__(self, name, dt, pressure_threshold):
        super().__init__('CellMigration', dt)
        self.is_global = True #run on the whole world, after the other processes
        self.pressure_threshold = pressure_threshold

    def __call__(self, world: World):
        exchange_matrix = world.compute_exchange_matrix(self.dt, pressure_threshold=self.pressure_threshold)
        for voxel in world.voxel_list:
            voxel_num = voxel.voxel_number
            if voxel_num % 10000 == 0: print('voxel number = ', voxel_num)
            list_of_neighbors = world.find_neighbors(voxel)
            np.random.shuffle(list_of_neighbors) #shuffle the list to avoid bias
            for neighbor in list_of_neighbors:
                n_events = exchange_matrix[voxel_num, neighbor.voxel_number] #number of expected events in the time step
                n_moving_cells = np.random.poisson(n_events)
                n_moving_cells = min(n_moving_cells, int(round(len(voxel.list_of_cells)*0.5)))
                list_of_moving_cells = np.random.choice(voxel.list_of_cells, n_moving_cells, replace=False)
                for cell in list_of_moving_cells:
                    if neighbor.add_cell(cell):
                        voxel.remove_cell(cell)

class UpdateCellOxygen(Process):
    def __init__(self, name, dt, voxel_half_length, effective_vessel_radius):
        super().__init__('UpdateState', dt)
        self.voxel_side = int(voxel_half_length*200) #um/100
        self.effective_vessel_radius = effective_vessel_radius

        alpha_file_name = 'Micro-Oxygenation/save_alpha_dataframe' + str(self.voxel_side) + '.csv'
        beta_file_name = 'Micro-Oxygenation/save_beta_dataframe' + str(self.voxel_side) + '.csv'
        if not os.path.isfile(alpha_file_name) or not os.path.isfile(beta_file_name):
            print('voxel side', self.voxel_side)
            raise ValueError('alpha/beta file not found! It might be in the wrong directory or information for chosen voxel size is not stored. Check "BetaDistributionCalibration.py" to generate the file for the chosen voxel size.')

        alpha_dataframe = pd.read_csv(alpha_file_name, index_col=0)
        beta_dataframe = pd.read_csv(beta_file_name, index_col=0)
        # Read the saved dataframes from the CSV files
        pressure_column = alpha_dataframe.index.values
        n_column = alpha_dataframe.columns.values.astype(float)

        # Create a 2D grid of points (pressure, n)
        points = []
        values_alpha = []
        values_beta = []

        for p in pressure_column:
            for n in n_column:
                alpha_value = alpha_dataframe.at[p, str(n)]
                beta_value = beta_dataframe.at[p, str(n)]
                points.append([p, n])
                values_alpha.append(alpha_value)
                values_beta.append(beta_value)

        self.alpha_map = ScalarField2D(points, values_alpha, bounds_error=False, fill_value= None)
        self.beta_map = ScalarField2D(points, values_beta, bounds_error=False, fill_value= None)

        if config.show_alpha_beta_maps:
            # Plot the alpha and beta maps
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.figure.set_dpi(100)
            self.alpha_map.show_extra(fig, ax, [min(pressure_column), max(pressure_column)], [min(n_column), max(n_column)])
            plt.show()

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.figure.set_dpi(100)
            self.beta_map.show_extra(fig, ax, [min(pressure_column), max(pressure_column)], [min(n_column), max(n_column)])
            plt.show()

    def __call__(self, voxel: Voxel):

        n_vessels = voxel.oxygen
        n_cells = voxel.number_of_alive_cells()
        pressure = voxel.pressure()

        if n_vessels == 0:
            o2_values = np.zeros(n_cells)

        elif n_vessels > 100:
            o2_values = np.ones(n_cells)

        else:
            alpha_ = self.alpha_map.evaluate((pressure, n_vessels))
            beta_ = self.beta_map.evaluate((pressure, n_vessels))

            if alpha_ < 0 or beta_ < 0:
                print('pressure', pressure, 'n_vessels', n_vessels)
                print('alpha_', alpha_, 'beta_', beta_)

            o2_values = np.random.beta(alpha_, beta_, size=n_cells)

        for i in range(n_cells):
            voxel.list_of_cells[i].oxygen = o2_values[i]
class UpdateVoxelMolecules(Process): #update the molecules in the voxel (VEGF), other not implemented yet
    def __init__(self, name, dt, VEGF_production_per_cell, threshold_for_VEGF_production):
        super().__init__('UpdateMolecules', dt)
        self.VEGF_production_per_cell = VEGF_production_per_cell
        self.threshold_for_VEGF_production = threshold_for_VEGF_production
    def __call__(self, voxel: Voxel):
        VEGF = 0
        for cell in voxel.list_of_cells:
            if cell.vitality() < self.threshold_for_VEGF_production and cell.type == 'TumorCell':
                VEGF = VEGF + self.VEGF_production_per_cell*(1-cell.vitality())
        VEGF = min(VEGF, 1.0)
        voxel.molecular_factors['VEGF'] = VEGF
        return
class UpdateVasculature(Process): #update the vasculature
    def __init__(self, name, dt, killing_radius_threshold, killing_length_threshold, o2_per_volume, diffusion_number, splitting_rate, macro_steps, micro_steps, weight_direction, weight_vegf, weight_pressure, radius_pressure_sensitive):
        super().__init__('UpdateVasculature', dt)
        self.is_global = True
        self.killing_radius_threshold = killing_radius_threshold
        self.killing_length_threshold = killing_length_threshold
        self.o2_per_volume = o2_per_volume
        self.diffusion_number = diffusion_number
        self.dt = dt
        self.splitting_rate = splitting_rate
        self.macro_steps = macro_steps
        self.micro_steps = micro_steps
        self.weight_direction = weight_direction
        self.weight_vegf = weight_vegf
        self.weight_pressure = weight_pressure
        self.radius_pressure_sensitive = radius_pressure_sensitive


    def __call__(self, world: World):
        #print in separate thread
        n_killed = world.vessels_killed(self.killing_radius_threshold, self.killing_length_threshold)
        print('Killed vessels: ', n_killed)
        print('Growing vessels')
        total_VEGF = 0
        for voxel in world.voxel_list:
            total_VEGF += voxel.molecular_factors['VEGF']
        print('Total VEGF: ', total_VEGF)
        vessels = world.vasculature.list_of_vessels
        n_new_vessels = int(config.new_vessels_per_hour * self.dt)
        chosen_vessels = random.sample(vessels, n_new_vessels)

        for vessel in chosen_vessels:
            if len(vessel.path) > 2:
                point = vessel.choose_random_point()
                world.vasculature.branching(vessel.id, point)

        # for i in range(n_new_vessels):
        #     #order by vessels radius, so that the largest vessels are more likely to branch
        #     if i >= len(vessels):
        #         i = i - len(vessels)
        #     vessel = vessels[i]
        #     #choose random point on vessel
        #     if len(vessel.path) > 2:
        #         point = vessel.choose_random_point()
        #         #branch from this point
        #         world.vasculature.branching(vessel.id, point)

        world.vasculature_growth(self.dt, self.splitting_rate, self.macro_steps, self.micro_steps, self.weight_direction, self.weight_vegf, self.weight_pressure, self.radius_pressure_sensitive)
        world.update_volume_occupied_by_vessels()
        # world.vasculature.print_vessel_tree()
        world.update_oxygen(o2_per_volume = self.o2_per_volume, diffusion_number=self.diffusion_number)

class Irradiation(Process): #irradiation
    def __init__(self, name, dt, topas_file, irradiation_time, irradiation_intensity, world: World):
        super().__init__('Irradiation', dt)
        self.irradiation_time = irradiation_time
        self.irradiation_intensity = irradiation_intensity

        #check if the file exists
        if not os.path.isfile('TopasSimulation/' + topas_file + '.csv'):
            world.topas_param_file(topas_file)
            term.RunTopasSimulation(topas_file)
            # rename the file

            os.rename('TopasSimulation/MyScorer.csv', 'TopasSimulation/' + topas_file + '.csv')
        _, self.doses = rw.DoseOnWorld('TopasSimulation/' + topas_file + '.csv')
        world.update_dose(self.doses)

        # plot the simulation
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        world.show_voxels_centers_dose(ax, fig, slice=True)
        plt.show()
    def __call__(self, world: World):
        for voxel in world.voxel_list:
            count = 0
            probability = self.doses[voxel.voxel_number]*self.irradiation_intensity*config.death_rate_radiation
            for cell in voxel.list_of_cells:
                if random.random() < probability:
                    if cell.time_before_death is None:
                        cell.time_before_death = random.lognormvariate(1, 1)
            for n_cell in voxel.list_of_necrotic_cells:
                if random.random() < probability:
                    n_cell.time_before_death = random.lognormvariate(1, 1)
        return
