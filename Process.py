from World import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Cell import *
from Voxel import *
from scipy.stats import beta
import ReadAndWrite as rw
from BasicGeometries import *
import pandas as pd
from ScalarField import *
import os

CONFIG = rw.read_config_file('CONFIG.txt')


class Simulator:
    def __init__(self, list_of_process : list, finish_time, dt):
        self.list_of_process = list_of_process
        self.finish_time = finish_time
        self.dt = dt
        self.time = 0


    def show(self, world: World, t = 0, slice = False):

        Vasculature_show = True

        if slice:
            angle = 0
            angle2 = 0
        else:
            angle = 30
            angle2 = 60

        DPI = 100
        size = world.half_length

        #plot vasculature
        if Vasculature_show:
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(25, 25), dpi=DPI, subplot_kw={'projection': '3d'})
            fig.suptitle('Visualization at time t = ' + str(t) + ' hours', fontsize=16)
            axes.view_init(30, 60)
            axes.set_xlim(-size/2, size/2)
            axes.set_ylim(-size/2, size/2)
            axes.set_zlim(-size/2, size/2)
            world.show_tumor(axes, fig)
            world.vasculature.plot(fig, axes)
            axes.set_title('Vasculature')
            plt.savefig('Plots/Video/t' + str(t) + '_Vasculature.png')
            plt.show()


        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(20, 20), dpi=DPI, subplot_kw={'projection': '3d'})
        fig.suptitle('Visualization at time t = ' + str(t) + ' hours', fontsize=16)

        axes[0, 0].view_init(angle, angle2)
        axes[0, 0].set_xlim(-size, size)
        axes[0, 0].set_ylim(-size, size)
        axes[0, 0].set_zlim(-size, size)
        world.show_tumor(axes[0, 0], fig, slice=slice)
        #world.vasculature.plot(fig, axes[1, 1])
        axes[0, 0].set_title('Cells in voxels')

        axes[0, 1].view_init(angle, angle2)
        axes[0, 1].set_xlim(-size, size)
        axes[0, 1].set_ylim(-size, size)
        axes[0, 1].set_zlim(-size, size)
        world.show_voxels_centers_oxygen(axes[0, 1], fig, slice=slice)
        axes[0, 1].set_title('Oxygen in voxels')

        axes[1, 0].view_init(angle, angle2)
        axes[1, 0].set_xlim(-size, size)
        axes[1, 0].set_ylim(-size, size)
        axes[1, 0].set_zlim(-size, size)
        world.show_voxels_centers_molecules(axes[1, 0], fig, slice=slice, molecule='VEGF')
        #show gradient of VEGF
        #world.show_voxels_centers_molecules(axes[1, 0], fig, slice=slice, molecule='VEGF', gradient=True)
        axes[1, 0].set_title('VEGF in voxels')

        axes[1, 1].view_init(angle, angle2)
        axes[1, 1].set_xlim(-size, size)
        axes[1, 1].set_ylim(-size, size)
        axes[1, 1].set_zlim(-size, size)
        world.show_necrosis(axes[1, 1], fig, slice=slice)
        #world.vasculature.plot(fig, axes[1, 1])
        axes[1, 1].set_title('Necrosis in voxels')

        plt.savefig('Plots/Video/t' + str(t) + '_AllPlots.png')
        plt.show()

        voxels_positions = [[0,0,0],[0.06, 0.06, 0.0], [0.12,0.12,0.0]]
        fig, axes = plt.subplots(nrows=2, ncols=len(voxels_positions), figsize=(20, 10), dpi=100)
        fig.suptitle('Visualization at time t = ' + str(t) + ' hours', fontsize=16)
        for i in range(len(voxels_positions)):
            #show histograms for the three voxels
            axes[0, i].set_title('Voxel ' + str(i))
            axes[0, i].set_xlabel('Oxygen')
            axes[0, i].set_ylabel('Number of cells')
            voxel = world.find_voxel(voxels_positions[i])
            voxel.oxygen_histogram(axes[0, i], fig)
            axes[1, i].set_xlabel('Vitality')
            axes[1, i].set_ylabel('Number of cells')
            voxel.vitality_histogram(axes[1, i], fig)
        plt.show()




    def run(self, world: World, video = False):
        slice = True
        print('Running simulation for {} hours'.format(self.finish_time), ' with dt = ', self.dt)
        number_cells = []
        process_local = [process for process in self.list_of_process if (not process.is_global)]
        process_global = [process for process in self.list_of_process if process.is_global]

        self.show(world, self.time, slice = slice)

        while self.time < self.finish_time:

            count_cells = 0
            print('\033[1;31;47mTime: {} hours'.format(self.time) + ' / ' + str(self.finish_time) + ' hours\033[0m')
            #loop in random order
            copy_voxel_list = world.voxel_list.copy()
            np.random.shuffle(copy_voxel_list)
            for voxel in world.voxel_list:
                count_cells += voxel.number_of_tumor_cells()
                for process in process_local:
                    process(voxel)
            for process in process_global:
                print('Currently running global process : ', process.name)
                process(world)
            self.time = self.time + self.dt
            if video: self.show(world, self.time, slice = slice)
            number_cells.append(count_cells)
            print('Vessel volume in center voxel: ', world.find_voxel([0, 0, 0]).vessel_volume)

        print('Simulation finished')
        self.show(world, self.time, slice = slice)
        fig_final = plt.figure()
        #plot number of cells evolution
        plt.plot(np.linspace(0, self.finish_time, len(number_cells)), number_cells, 'purple')
        plt.title('Number of cells evolution')
        plt.xlabel('Time (hours)')
        plt.ylabel('Number of cells')
        plt.savefig('Plots/Number_cells_evolution.png')
        plt.show()
        return


class Process:
    def __init__(self, name, dt):
        self.name = name
        self.dt = dt
        self.is_global = False
    def __call__(self, voxel):
        pass


class CellDivision(Process):
    def __init__(self, name, dt, cycling_threshold, pressure_threshold = np.inf):
        super().__init__(name, dt)
        self.cycling_threshold = cycling_threshold
        self.pressure_threshold = pressure_threshold

    def __call__(self, voxel):
        #print('CellDivision')
        if voxel.pressure() < self.pressure_threshold:
            for cell in voxel.list_of_cells:
                if cell.vitality() > self.cycling_threshold:
                    # Calculate the expected number of cell divisions in the time step
                    expected_divisions = self.dt / cell.doubling_time
                    # Use a Poisson distribution to model the number of cell divisions
                    num_divisions = np.random.poisson(expected_divisions)
                    for i in range(num_divisions):
                        new_cell = cell.duplicate()
                        voxel.add_cell(new_cell)
        else:
            if CONFIG['verbose']: print('pressure = ', voxel.pressure(), ' > ', self.pressure_threshold, ' so no cell division')
        return
class CellApoptosis(Process):
    def __init__(self, name, dt, apoptosis_threshold, apoptosis_probability):
        super().__init__('CellDeath', dt)
        self.apoptosis_threshold = apoptosis_threshold
        self.apoptosis_probability = apoptosis_probability

    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            if cell.vitality() < self.apoptosis_threshold and np.random.rand() < self.apoptosis_probability:
                print('Cell apoptosis')
                voxel.remove_cell(cell)
class CellNecrosis(Process):
    def __init__(self, name, dt, necrosis_threshold, necrosis_probability):
        super().__init__('CellNecrosis', dt)
        self.necrosis_threshold = necrosis_threshold
        self.necrosis_probability = necrosis_probability

    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            if cell.vitality() < self.necrosis_threshold:
                if cell.vitality() < self.necrosis_threshold and np.random.rand() < self.necrosis_probability:
                    if CONFIG['verbose']: print('Cell necrosis in voxel ', voxel.voxel_number)
                    voxel.cell_becomes_necrotic(cell)

class CellAging(Process):
    def __init__(self, name, dt):
        super().__init__('CellAging', dt)
    def __call__(self, voxel):
        #print('CellAging')
        pass


class CellMigration(Process):
    def __init__(self, name, dt, pressure_threshold):
        super().__init__('CellMigration', dt)
        self.is_global = True
        self.pressure_threshold = pressure_threshold

    def __call__(self, world: World):
        exchange_matrix = world.compute_exchange_matrix(self.dt, pressure_threshold=self.pressure_threshold)
        for voxel in world.voxel_list:
            voxel_num = voxel.voxel_number
            if voxel_num % 10000 == 0: print('voxel number = ', voxel_num)
            list_of_neighbors = world.find_neighbors(voxel)
            np.random.shuffle(list_of_neighbors)
            for neighbor in list_of_neighbors:
                n_events = exchange_matrix[voxel_num, neighbor.voxel_number]
                n_moving_cells = np.random.poisson(n_events)
                n_moving_cells = min(n_moving_cells, len(voxel.list_of_cells))
                if n_moving_cells > 0:
                    list_of_moving_cells = np.random.choice(voxel.list_of_cells, n_moving_cells, replace=False)
                    for cell in list_of_moving_cells:
                        if not cell.necrotic:
                            if neighbor.add_cell(cell):
                                voxel.remove_cell(cell)

class UpdateCellOxygen(Process):
    def __init__(self, name, dt, voxel_half_length, effective_vessel_radius):
        super().__init__('UpdateState', dt)
        self.voxel_side = int(voxel_half_length*200) #um/100

        #read alpha and beta from files
        alpha = np.genfromtxt('alpha.csv', delimiter=',', skip_header=True)
        beta = np.genfromtxt('beta.csv', delimiter=',', skip_header=True)
        a_row_index = np.where(alpha[:, 0] == self.voxel_side)[0][0]
        b_row_index = np.where(beta[:, 0] == self.voxel_side)[0][0]
        #add check that the voxel size is stored in alpha/beta file
        self.aA = alpha[a_row_index, 1]
        self.aB = alpha[a_row_index, 2]
        self.aC = alpha[a_row_index, 3]
        self.aD = alpha[a_row_index, 4]
        self.bA = beta[b_row_index, 1]
        self.bB = beta[b_row_index, 2]
        self.bC = beta[b_row_index, 3]
        self.bD = beta[b_row_index, 4]
        self.effective_vessel_radius = effective_vessel_radius

        alpha_file_name = 'alpha_dataframe' + str(self.voxel_side) + '.csv'
        beta_file_name = 'beta_dataframe' + str(self.voxel_side) + '.csv'
        if not os.path.isfile(alpha_file_name) or not os.path.isfile(beta_file_name):
            raise ValueError('alpha/beta file not found! It might be in the wrong directory or information for chosen voxel size is not stored. Check "BetaDistributionCalibration.py" to generate the file for the chosen voxel size.')

        alpha_dataframe = pd.read_csv(alpha_file_name, index_col=0)
        beta_dataframe = pd.read_csv(beta_file_name, index_col=0)
        print('alpha_dataframe', alpha_dataframe)
        # Read the saved dataframes from the CSV files
        pressure_column = alpha_dataframe.index.values
        n_column = alpha_dataframe.columns.values.astype(float)
        print('pressure_column', pressure_column)
        print('n_column', n_column)

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

    def __call__(self, voxel: Voxel):

        n_vessels = voxel.oxygen
        n_cells = voxel.number_of_alive_cells()
        pressure = voxel.pressure()


        if n_vessels == 0:
            o2_values = np.zeros(n_cells)

        else:
            alpha_ = self.alpha_map.evaluate((pressure, n_vessels))
            beta_ = self.beta_map.evaluate((pressure, n_vessels))

            if alpha_ < 0 or beta_ < 0:
                print('pressure', pressure, 'n_vessels', n_vessels)
                print('alpha_', alpha_, 'beta_', beta_)

            o2_values = np.random.beta(alpha_, beta_, size=n_cells)

        for i in range(n_cells):
            voxel.list_of_cells[i].oxygen = o2_values[i]
class UpdateVoxelMolecules(Process):
    def __init__(self, name, dt, VEGF_production_per_cell, threshold_for_VEGF_production):
        super().__init__('UpdateMolecules', dt)
        self.VEGF_production_per_cell = VEGF_production_per_cell
        self.threshold_for_VEGF_production = threshold_for_VEGF_production
    def __call__(self, voxel: Voxel):
        VEGF = 0
        for cell in voxel.list_of_cells:
            if cell.type == 'TumorCell' and cell.vitality() < self.threshold_for_VEGF_production:
                VEGF = VEGF + self.VEGF_production_per_cell*(1-cell.vitality())
        VEGF = min(VEGF, 1.0)
        voxel.molecular_factors['VEGF'] = VEGF
        return
class UpdateVasculature(Process):
    def __init__(self, name, dt, pressure_killing_radius_threshold, o2_per_volume, diffusion_number, splitting_rate, macro_steps, micro_steps, weight_direction, weight_vegf, weight_pressure, radius_pressure_sensitive):
        super().__init__('UpdateVasculature', dt)
        self.is_global = True
        self.pressure_killing_radius_threshold = pressure_killing_radius_threshold
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
        n_killed = world.vessels_killed_by_pressure(self.pressure_killing_radius_threshold)
        print('Killed vessels: ', n_killed)
        print('Growing vessels')
        total_VEGF = 0
        for voxel in world.voxel_list:
            total_VEGF += voxel.molecular_factors['VEGF']
        print('Total VEGF: ', total_VEGF)
        for i in range(int(total_VEGF)):
            #choose random vessel
            vessel = np.random.choice(world.vasculature.list_of_vessels)
            #choose random point on vessel
            point = vessel.choose_random_point()
            #branch from this point
            world.vasculature.branching(vessel.id, point)

        world.vasculature_growth(self.dt, self.splitting_rate, self.macro_steps, self.micro_steps, self.weight_direction, self.weight_vegf, self.weight_pressure, self.radius_pressure_sensitive)
        world.update_volume_occupied_by_vessels()
        world.update_oxygen(o2_per_volume = self.o2_per_volume, diffusion_number=self.diffusion_number)