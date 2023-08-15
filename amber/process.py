from amber.world import *
from amber.voxel import *
from amber.ScalarField import *
import amber.terminal as term
import amber.ReadAndWrite as rw
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
import matplotlib.pyplot as plt
import os
from time import time as current_time

class Simulator: #this class is used to run the whole simulation

    def __init__(self, config, list_of_process : list, finish_time, dt):
        self.list_of_process = list_of_process
        self.list_of_process_names = [process.name for process in list_of_process]
        self.finish_time = finish_time
        self.dt = dt
        self.time = 0
        self.config = config

        if not os.path.exists('DataOutput/'):
            os.makedirs('DataOutput/')

        if not os.path.exists('Plots/'):
            os.makedirs('Plots/')

    def show_center_of_mass(self, center_of_mass, times): #3D plot of the center of mass
        #3D plot of the center of mass
        #set dpi to 300 for high quality
        size = self.config.half_length_world
        center_of_mass = np.array(center_of_mass)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(center_of_mass[:, 0], center_of_mass[:, 1], center_of_mass[:, 2], 'black')
        ax.scatter(center_of_mass[0, 0], center_of_mass[0, 1], center_of_mass[0, 2], color ='green', label='start')
        ax.scatter(center_of_mass[-1, 0], center_of_mass[-1, 1], center_of_mass[-1, 2], color ='red', label='end')
        ax.set_xlim(-size, size)
        ax.set_ylim(-size, size)
        ax.set_zlim(-size, size)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Center of mass path')
        ax.legend()
        plt.savefig('Plots/Center_of_mass.png', dpi=100)
        if self.config.running_on_cluster: #if running on cluster, save plot to file and do not show
            plt.close()
        else:
            plt.show()

    def show_cell_and_tumor_volume(self, number_tumor_cells, number_necrotic_cells, number_quiescent_cells, number_cycling_cells, tumor_size, tumor_size_free, number_vessels, times): #plot the number of cells and tumor volume
        # plot number of cells evolution
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 12))
        #change font size
        plt.rcParams.update({'font.size': 14})
        # Plot number of cells evolution
        ax1.plot(times, number_tumor_cells, 'blue', label='All cells')
        ax1.plot(times, number_cycling_cells, 'red', label='Cycling cells')
        ax1.plot(times, number_quiescent_cells, 'green', label='Quiescent cells')
        ax1.plot(times, number_necrotic_cells, 'black', label='Necrotic cells')
        ax1.set_title('Number of cells evolution')
        ax1.set_xlabel('Time [h]', fontsize=14)
        ax1.set_ylabel('Number of cells', fontsize=14)
        ax1.grid(True)
        ax1.legend(fontsize=14)

        # Plot tumor size evolution
        ax2.plot(times, tumor_size, 'red')
        ax2.plot(times, tumor_size_free, 'blue')
        ax2.set_title('Tumor volume evolution', fontsize=14)
        ax2.set_xlabel('Time [h]', fontsize=14)
        ax2.set_ylabel('Tumor volume [mm^3]', fontsize=14)
        ax2.grid(True)

        # Plot number of vessels evolution
        ax3.plot(times, number_vessels, 'black')
        ax3.set_title('Number of vessels evolution', fontsize=14)
        ax3.set_xlabel('Time [h]', fontsize=14)
        ax3.set_ylabel('Number of vessels', fontsize=14)
        ax3.grid(True)

        # Adjust the spacing between subplots
        plt.tight_layout()

        # Save the figure
        fig.savefig('Plots/Combined_plots_tumor_evolution.png', dpi=100)
        if self.config.running_on_cluster:
            plt.close()
        else:
            plt.show()
    def show(self, world: World, t = 0): #this function is used to show the world at a certain time
        print('Showing world at time : ', t)
        start = current_time()

        if not os.path.exists('Plots/CurrentPlotting/'):
            os.makedirs('Plots/CurrentPlotting/')

        size = world.half_length

        if self.config.show_angiogenesis_metrics: #if angiogenesis metrics are to be shown, show them
            print('Showing angiogenesis metrics')
            world.show_angiogenesis_metrics(t, self.config.true_vasculature)

            #
            # fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(40, 20))
            # axes[0].set_xlim(-size, size)
            # axes[0].set_ylim(-size, size)
            # world.show_tumor_slice(axes[0], fig, 'vessel_length_density', levels= np.linspace(0, 200, 20), refinement_level=3, cmap='jet')
            # axes[0].grid(True)
            # axes[0].set_facecolor('whitesmoke')
            # axes[0].set_title('Vessel length density [mm/mm^3]')
            #
            # axes[1].set_xlim(-size, size)
            # axes[1].set_ylim(-size, size)
            # world.show_tumor_slice(axes[1], fig, 'vessel_volume_density', refinement_level=3, cmap='jet')
            # axes[1].grid(True)
            # axes[1].set_facecolor('whitesmoke')
            # axes[1].set_title('Vessel volume density [%]')
            #
            # plt.tight_layout()
            # plt.savefig('Plots/CurrentPlotting/t' + str(t) + '_VesselMetricsMaps.png', dpi=100)
            # if self.config.running_on_cluster:
            #     plt.close()
            # else:
            #     plt.show()


        #plot vasculature
        if self.config.show_tumor_and_vessels_3D:
            print('Showing tumor and vessels 3D')
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': '3d'})
            plt.rcParams.update({'font.size': 10})
            plt.title('Visualization at time t = ' + str(t) + ' hours', fontsize=16)
            axes.set_xlim(-size, size)
            axes.set_ylim(-size, size)
            axes.set_zlim(-size, size)
            #remove z axis ticks
            if self.config.slice == 'x':
                axes.set_zticks([])
            #change text size
            # axes.view_init(90, -90)
            # world.show_tumor_3D(axes, fig, 'number_of_tumor_cells', cmap='viridis', vmin=0, vmax=1000)
            world.vasculature.plot(fig, axes)
            plt.tight_layout()
            plt.savefig('Plots/CurrentPlotting/t' + str(t) + '_Vasculature.png', dpi=300)
            if self.config.running_on_cluster:
                plt.close()
            else:
                plt.show()

        if self.config.show_slices:
            font = 22
            print('Showing slices')
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 12))
            fig.suptitle('t = ' + str(t) + 'h', fontsize=16)

            axes[0, 0].set_xlim(-size, size)
            axes[0, 0].set_ylim(-size, size)
            world.show_tumor_slice(axes[0, 0], fig, 'number_of_tumor_cells', levels= np.linspace(1, 1001, 11), cmap='viridis', extend = 'neither')
            axes[0,0].grid(True)
            axes[0,0].set_facecolor('whitesmoke')
            axes[0, 0].set_title('Number of Cells', fontsize=font)

            norm = TwoSlopeNorm(vmin=0, vcenter=10, vmax=110)

            axes[0, 1].set_xlim(-size, size)
            axes[0, 1].set_ylim(-size, size)
            world.show_tumor_slice(axes[0, 1], fig, 'n_capillaries', cmap = 'RdBu', norm = norm, levels= np.linspace(0, 100, 11), round_n = 0)
            axes[0, 1].grid(True)
            axes[0, 1].set_facecolor('whitesmoke')
            axes[0, 1].set_title('Number of Capillaries', fontsize=font)

            axes[1, 0].set_xlim(-size, size)
            axes[1, 0].set_ylim(-size, size)
            world.show_tumor_slice(axes[1, 0], fig, 'molecular_factors', factor='VEGF', levels= np.linspace(0.001, 1.0, 11), cmap='Oranges', round_n = 2)
            axes[1, 0].grid(True)
            axes[1, 0].set_facecolor('whitesmoke')
            axes[1, 0].set_title('VEGF concentration', fontsize=font)

            axes[1, 1].set_xlim(-size, size)
            axes[1, 1].set_ylim(-size, size)
            world.show_tumor_slice(axes[1, 1], fig, 'number_of_necrotic_cells', levels= np.linspace(1, 1001, 11), cmap='viridis', extend = 'neither')
            axes[1, 1].grid(True)
            axes[1, 1].set_facecolor('whitesmoke')
            axes[1, 1].set_title('Number of Necrotic Cells', fontsize=font)

            plt.tight_layout()
            plt.savefig('Plots/CurrentPlotting/t' + str(t) + '_AllPlots.png', dpi=100)
            if self.config.running_on_cluster:
                plt.close()
            else:
                plt.show()

        if self.config.show_o2_vitality_histograms:
            print('Showing histograms')
            voxel_side = 2*self.config.half_length_world / self.config.voxel_per_side

            voxels_positions = [[0,0,0], [3*voxel_side, 3*voxel_side, 3*voxel_side], [5*voxel_side, 5*voxel_side, 5*voxel_side]]
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
            if self.config.running_on_cluster:
                plt.close()
            else:
                plt.show()

        if self.config.show_cell_damage:
            print('Showing cell damage')
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), dpi=100)
            fig.suptitle('Visualization at time t = ' + str(t) + ' hours', fontsize=16)
            world.show_tumor_slice(axes, fig, 'average_cell_damage', levels= np.linspace(0, 1, 10), cmap='cool', extend = 'neither')
            if self.config.running_on_cluster:
                plt.close()
            else:
                plt.show()

        end = current_time()
        print('Time elapsed for showing graphs: ' + str(end - start) + ' seconds')
    def run(self, world: World, video=False): #run the simulation! (the main function)

        print(f'Running simulation for {self.finish_time} hours with dt={self.dt}')
        process_local = [process for process in self.list_of_process if not process.is_global] #list of local processes
        process_global = [process for process in self.list_of_process if process.is_global] #list of global processes

        #create an array of execution times for each process (local and global)
        execution_times = [[] for _ in range(len(self.list_of_process))] #list of execution times for each process (local and global)
        sums = []

        # if not self.config.skip_weekends:
        #     irradiations_times = [self.config.first_irradiation_time + i * self.config.time_between_fractions for i in
        #                           range(self.config.number_fractions)] #list of irradiation times
        #read the fractionation schedule from a file
        irradiations_times = []
        with open(str(self.config.fractionation_schedule), 'r') as f:
            for line in f:
                irradiations_times.append(float(line))

        print('Irradiation times: ' + str(irradiations_times))
        can_irradiate = False
        if not self.config.irradiation_cell_based:
            can_irradiate = True
            for i in range(len(irradiations_times)):
                irradiations_times[i] = irradiations_times[i] + self.config.first_irradiation_time

        number_of_fractions = len(irradiations_times) #number of fractions
        applied_fractions = 0 #number of applied fractions

        number_cycling_cells = []; number_quiescent_cells = []; number_necrotic_cells = [];
        tumor_size = []; tumor_size_free = []; times = []; number_tumor_cells = []; number_vessels = [];
        center_of_mass = []

        while self.time < self.finish_time: #while the simulation is running
            print(f'\033[1;31;47mTime: {self.time} hours / {self.finish_time} hours\033[0m')
            #show the simulation
            if video:
                self.show(world, self.time) #show the simulation

            #do irradiation if needed
            if can_irradiate and applied_fractions < number_of_fractions and self.time >= irradiations_times[applied_fractions]:
                irrad = Irradiation(self.config, 'irrad', self.dt, self.config.topas_file,
                                    self.config.irradiation_intensity, world)
                irrad(world)
                applied_fractions += 1

            print('Currently running local processes:')
            start = current_time()
            for voxel in world.voxel_list:
                for process in process_local:
                    process(voxel) #run local processes on each voxel
            end = current_time()
            print('Time for locals processes:', end - start)
            for process in process_global: #run global processes on the whole world
                print('Currently running global process:', process.name)
                start = current_time()
                process(world)
                end = current_time()
                print('Time for process', process.name, ':', end - start)

            cycling_cells = 0
            quiescent_cells = 0
            necrotic_cells = 0
            start = current_time()
            for voxel in world.voxel_list: #count the number of cells in each state
                necrotic_cells += voxel.number_of_necrotic_cells()
                for cell in voxel.list_of_cells:
                    if cell.type == 'TumorCell':
                        if cell.vitality() > self.config.vitality_cycling_threshold:
                            cycling_cells += 1
                        else:
                            quiescent_cells += 1
            end = current_time()
            print('Time for counting cells:', end - start)

            number_cycling_cells.append(cycling_cells)
            number_quiescent_cells.append(quiescent_cells)
            number_necrotic_cells.append(necrotic_cells)
            number_tumor_cells.append(cycling_cells + quiescent_cells + necrotic_cells)
            tumor_size_, tumor_size_free_ = world.measure_tumor_volume()
            tumor_size.append(tumor_size_)
            tumor_size_free.append(tumor_size_free_)
            number_vessels.append(len(world.vasculature.list_of_vessels))
            center_of_mass.append(world.center_of_mass)
            times.append(self.time)

            np.save('DataOutput/number_tumor_cells.npy', number_tumor_cells)
            np.save('DataOutput/number_necrotic_cells.npy', number_necrotic_cells)
            np.save('DataOutput/number_cycling_cells.npy', number_cycling_cells)
            np.save('DataOutput/number_quiescent_cells.npy', number_quiescent_cells)
            np.save('DataOutput/tumor_size.npy', tumor_size)
            np.save('DataOutput/tumor_size_free.npy', tumor_size_free)
            np.save('DataOutput/number_vessels.npy', number_vessels)
            np.save('DataOutput/center_of_mass.npy', center_of_mass)
            np.save('DataOutput/times.npy', times)

            if self.time % self.config.save_world_every == 0 and self.time != 0:
                world.save('t'+str(self.time)+str(self.config.world_file) + str(self.config.seed) + '.pkl')

            if self.config.show_cell_and_tumor_volume:
                self.show_cell_and_tumor_volume(number_tumor_cells, number_necrotic_cells, number_quiescent_cells, number_cycling_cells, tumor_size, tumor_size_free, number_vessels, times)


            if self.config.show_center_of_mass:
                self.show_center_of_mass(center_of_mass, times)

            sum = 0
            for i, process in enumerate(self.list_of_process):
                sum += process.time_spent_doing_process
                execution_times[i].append(process.time_spent_doing_process)
                process.time_spent_doing_process = 0
            sums.append(sum)

            if not can_irradiate and self.config.irradiation_cell_based and number_tumor_cells[-1] >= self.config.critical_n_cells:
                can_irradiate = True
                for i in range(len(irradiations_times)):
                    irradiations_times[i] = irradiations_times[i] + self.time

            #print('Time spent doing processes:', execution_times)
            plt.figure()
            for i, process in enumerate(self.list_of_process):
                plt.plot(times, execution_times[i], label = process.name)
            plt.plot(times, sums, label = 'Total')
            plt.legend()
            plt.xlabel('Time (hours)')
            plt.ylabel('Time spent doing process (s)')
            plt.yscale('log')
            plt.grid()
            plt.savefig('Plots/execution_times.png')
            plt.close()

            self.time += self.dt

        print('Simulation finished')

        if self.config.show_final:
            self.show(world, self.time)


        return

class Process: #abstract class, represents all the processes that can happen in the simulation

    time_spent_doing_process = 0
    def __init__(self, config, name, dt):
        self.name = name
        self.dt = dt
        self.is_global = False
        self.config = config

    def __call__(self, voxel):
        pass

    @classmethod #decorator to measure the time spent doing a process
    def timeit(cls, process_func):
        def wrapped_process(self, voxel):
            start_time = current_time()
            result = process_func(self, voxel)
            end_time = current_time()
            elapsed_time = end_time - start_time
            self.time_spent_doing_process += elapsed_time
            return result
        return wrapped_process


class CellDivision(Process): #cell division process, cells divide in a voxel if they have enough vitality
    def __init__(self, config, name, dt, cycling_threshold, pressure_threshold = np.inf):
        super().__init__(config, 'CellDivision', dt)
        self.dt = dt
        self.cycling_threshold = cycling_threshold
        self.pressure_threshold = pressure_threshold

    @Process.timeit
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
                    voxel.add_cell(new_cell, self.config.max_occupancy) #add the new cell to the voxel
        return

class CellDeath(Process): #cell necrosis process, cells die in a voxel if they have too low vitality
    def __init__(self, config, name, dt, necrosis_threshold, necrosis_probability, apoptosis_threshold, apoptosis_probability,necrosis_removal_probability, apoptosis_removal_probability, necrosis_damage_coeff, apoptosis_damage_coeff):
        super().__init__(config, 'CellNecrosis', dt)
        self.necrosis_threshold = necrosis_threshold
        self.necrosis_probability = necrosis_probability
        self.apoptosis_threshold = apoptosis_threshold
        self.apoptosis_probability = apoptosis_probability
        self.necrosis_damage_coeff = necrosis_damage_coeff
        self.apoptosis_damage_coeff = apoptosis_damage_coeff
        if self.necrosis_threshold > self.apoptosis_threshold:
            raise ValueError('necrosis threshold must be smaller or equal to apoptosis threshold. you can set apoptosis probability to 0 if you want to avoid apoptosis.')
        self.necrosis_removal_probability = necrosis_removal_probability
        self.apoptosis_removal_probability = apoptosis_removal_probability

    @Process.timeit
    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            sample = np.random.uniform(0, 1)
            #probability of necrosis and apoptosis. use math to get sampling every hour
            p_necro = (1 - ((1-cell.necrosis_probability(self.necrosis_probability,self.necrosis_threshold, self.necrosis_damage_coeff))**self.dt))
            p_apopt = (1 - ((1-cell.apoptosis_probability(self.apoptosis_probability, self.apoptosis_threshold, self.apoptosis_damage_coeff))**self.dt))
            if self.config.verbose: print('probability necro:', p_necro, 'probability apopto:', p_apopt)
            p_tot = p_necro + p_apopt
            if p_tot > 1:
                p_necro = p_necro/p_tot; p_apopt = p_apopt/p_tot

            if sample < p_necro:
                #necrosis
                voxel.cell_becomes_necrotic(cell)
            elif sample < p_necro + p_apopt:
                #apoptosis
                voxel.cell_becomes_apoptotic(cell)

        for dead in voxel.list_of_dead_cells: #remove dead cells with a certain probability
            if dead.necrotic: p = self.necrosis_removal_probability
            else: p = self.apoptosis_removal_probability
            proba = (1 - ((1-p)**self.dt))
            if random.random() < proba:
                voxel.remove_dead_cell(dead)

class CellAging(Process): #cell aging process, cells age in a voxel
    def __init__(self, config, name, dt, repair_per_hour):
        super().__init__(config,'CellAging', dt)
        self.repair_per_hour = repair_per_hour

    @Process.timeit
    def __call__(self, voxel):
        for cell in voxel.list_of_cells:

            if cell.is_cycling(self.config.vitality_cycling_threshold, 0.0): #cell aging towards duplication
                cell.time_spent_cycling += self.dt

            if cell.damage > 0: #cell being repaired
                cell.damage_repair(self.dt, self.repair_per_hour)

        pass

class CellInteraction(Process):

    def __init__(self, config, name, dt):
        super().__init__(config, 'CellInteraction', dt)
        self.dt = dt

    @Process.timeit
    def __call__(self, voxel): #cell interaction process, cells interact in a voxel
        matrix = voxel.compute_cell_interaction_matrix(self.dt)
        number_cells = len(voxel.list_of_cells)
        #find non zero elements of the sparse matrix
        idd = matrix.nonzero()
        for i in range(len(idd[0])):
            if np.random.random() < matrix[idd[0][i], idd[1][i]]:
                cell1 = voxel.list_of_cells[idd[0][i]]
                cell2 = voxel.list_of_cells[idd[1][i]]
                cell1.interact(cell2, self.dt) #do the interaction between the two cells

class CellMigration(Process): #cell migration process, cells migrate in the world
    def __init__(self, config, name, dt):
        super().__init__(config, 'CellMigration', dt)
        self.is_global = True #run on the whole world, after the other processes
    @Process.timeit
    def __call__(self, world: World):
        exchange_matrix = world.compute_exchange_matrix(self.dt) #compute the exchange matrix for the time step
        for voxel in world.voxel_list:
            voxel_num = voxel.voxel_number
            if voxel_num % 10000 == 0: print('voxel number = ', voxel_num)
            list_of_neighbors = world.find_moor_neighbors(voxel)
            np.random.shuffle(list_of_neighbors) #shuffle the list to avoid bias
            for neighbor in list_of_neighbors:
                n_events = exchange_matrix[voxel_num, neighbor.voxel_number] #number of expected events in the time step
                n_moving_cells = np.random.poisson(n_events)
                n_moving_cells = min(n_moving_cells, int(round(len(voxel.list_of_cells))))
                list_of_moving_cells = np.random.choice(voxel.list_of_cells, n_moving_cells, replace=False) #choose the cells to move randomly
                for cell in list_of_moving_cells: #move the cells
                    if neighbor.add_cell(cell, self.config.max_occupancy):
                        voxel.remove_cell(cell)


class UpdateCellOxygen(Process):
    def __init__(self, config, name, dt, voxel_half_length, file_prefix_alpha_beta_maps):
        super().__init__(config, 'UpdateState', dt)
        self.voxel_side = int(voxel_half_length*20) #um/100

        #read alpha and beta maps from csv files
        amber_dir = os.path.abspath(os.path.dirname(__file__))
        alpha_file_name = str(file_prefix_alpha_beta_maps) + '_alpha_dataframe'+str(self.voxel_side)+'.csv'
        beta_file_name = str(file_prefix_alpha_beta_maps) + '_beta_dataframe'+str(self.voxel_side)+'.csv'
        alpha_file_name = os.path.join(amber_dir, alpha_file_name)
        beta_file_name = os.path.join(amber_dir, beta_file_name)

        if not os.path.isfile(alpha_file_name) or not os.path.isfile(beta_file_name):
            print('voxel side', self.voxel_side)
            print('alpha file name is ', alpha_file_name)
            print('beta file name is ', beta_file_name)
            raise ValueError('alpha/beta file not found! It might be in the wrong directory or information for chosen voxel size is not stored. Check "BetaDistributionCalibration.py" to generate the file for the chosen voxel size.')

        alpha_dataframe = pd.read_csv(alpha_file_name, index_col=0)
        beta_dataframe = pd.read_csv(beta_file_name, index_col=0)

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

        if self.config.show_alpha_beta_maps:
            # Plot the alpha and beta maps
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(111, projection='3d')
            ax.view_init(azim=45, elev=30)
            self.alpha_map.show_extra(fig, ax, [min(pressure_column), max(pressure_column)], [min(n_column), max(n_column)])
            ax.axes.set_xlabel('Cell density')
            ax.axes.set_ylabel('Number of capillaries')
            ax.title.set_text('Alpha map')
            plt.savefig('alpha_map.png', dpi=300)
            if self.config.running_on_cluster:
                plt.close()
            else:
                plt.show()

            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(111, projection='3d')
            ax.view_init(azim=135, elev=30)
            self.beta_map.show_extra(fig, ax, [min(pressure_column), max(pressure_column)], [min(n_column), max(n_column)])
            ax.axes.set_xlabel('Cell density')
            ax.axes.set_ylabel('Number of capillaries')
            ax.title.set_text('Beta map')
            plt.savefig('beta_map.png', dpi=300)
            if self.config.running_on_cluster:
                plt.close()
            else:
                plt.show()
    @Process.timeit
    def __call__(self, voxel: Voxel):

        n_vessels = voxel.n_capillaries
        n_cells = voxel.number_of_alive_cells()
        pressure = voxel.pressure()
        if n_vessels == 0:
            o2_values = np.zeros(n_cells) #if there are no vessels, all cells have 0 oxygen

        elif n_vessels > 100:
            o2_values = np.ones(n_cells) #if there are too many vessels, all cells have 1 oxygen

        else:
            alpha_ = self.alpha_map.evaluate((pressure, n_vessels))
            beta_ = self.beta_map.evaluate((pressure, n_vessels))

            if alpha_ < 0 or beta_ < 0:
                print('pressure', pressure, 'n_vessels', n_vessels)
                print('alpha_', alpha_, 'beta_', beta_)

            o2_values = np.random.beta(alpha_, beta_, size=n_cells) #sample from beta distribution
        for i in range(n_cells):
            voxel.list_of_cells[i].oxygen = o2_values[i]
class UpdateVoxelMolecules(Process): #update the molecules in the voxel (VEGF), other not implemented yet
    def __init__(self, config, name, dt):
        super().__init__(config, 'UpdateMolecules', dt)

    def update_VEGF(self, voxel: Voxel): #update the VEGF concentration in the voxel
        VEGF = 0.0
        for cell in voxel.list_of_cells: #sum the VEGF secreted by each cell. Oxygen and and damage play a role in the secretion
            VEGF += cell.VEGF_secretion(self.config.metabolic_damage_threshold)
        VEGF = min(VEGF, 1.0)
        voxel.molecular_factors['VEGF'] = VEGF
        return
    def update_fiber_density(self, voxel: Voxel): #update the fiber density in the voxel
        fiber_density = 0.0
        for cell in voxel.list_of_cells:
            fiber_density += cell.fiber_secretion
        fiber_density = min(fiber_density, 1.0)
        voxel.molecular_factors['fiber_density'] = fiber_density
        return

    @Process.timeit
    def __call__(self, voxel: Voxel):
        self.update_VEGF(voxel) #update the VEGF concentration
        # self.update_fiber_density(voxel)
        return
class UpdateVasculature(Process): #update the vasculature
    def __init__(self, config, name, dt, killing_radius_threshold, n_capillaries_per_VVD, capillary_length, splitting_rate, macro_steps, micro_steps, weight_direction, weight_vegf, weight_pressure):
        super().__init__(config, 'UpdateVasculature', dt)
        self.is_global = True
        self.killing_radius_threshold = killing_radius_threshold
        self.n_capillaries_per_VVD = n_capillaries_per_VVD
        self.capillary_length = capillary_length
        self.dt = dt
        self.splitting_rate = splitting_rate
        self.macro_steps = macro_steps
        self.micro_steps = micro_steps
        self.weight_direction = weight_direction
        self.weight_vegf = weight_vegf
        self.weight_pressure = weight_pressure

    @Process.timeit
    def __call__(self, world: World):
        #print in separate thread
        n_killed = world.vessels_killed(self.killing_radius_threshold) #kill vessels that have a radius smaller than the threshold



        print('Killed vessels: ', n_killed)
        print('Growing vessels')

        vessels = world.vasculature.list_of_vessels
        volume_world = 8*world.half_length**3
        n_new_vessels = int(self.config.new_vessels_per_hour * self.dt * volume_world)
        n_new_vessels = min(n_new_vessels, len(vessels)) #the number of new vessels cannot be larger than the number of existing vessels
        vegf = world.vegf_map(step_gradient= self.config.vegf_map_step_gradient) #compute the gradient of the VEGF map

        # figure = plt.figure()
        # ax = figure.add_subplot(111)
        # vegf.show_values(figure,ax, 'viridis', 0.0, 1.0)
        # plt.show()

        def vegf_gradient(point): return vegf.gradient(point) #define the gradient of the VEGF map

        for _ in range(n_new_vessels): #create a tEC on some vessels randomly
            random_vessel = random.choice(vessels)
            if len(random_vessel.path) > 2:
                point = random_vessel.choose_random_point(self.config.seed)
                if vegf.evaluate(point) > self.config.vegf_scalar_threshold: #if the VEGF concentration is too low, tEC does not start growing
                    if np.linalg.norm(vegf_gradient(point)) > self.config.vegf_gradient_threshold: #if the gradient of the VEGF map is large enough, tEC starts growing
                        world.vasculature.branching(random_vessel.id, point)

        #grow the vessels and update the volume occupied by the vessels
        world.vasculature_growth(self.dt, self.splitting_rate, self.macro_steps, self.micro_steps, self.weight_direction, self.weight_vegf, self.weight_pressure)
        world.update_volume_occupied_by_vessels()
        #update the capillary map
        world.update_capillaries(n_capillaries_per_VVD= self.n_capillaries_per_VVD, capillary_length = self.capillary_length)

class Irradiation(Process): #irradiation
    def __init__(self, config, name, dt, topas_file, irradiation_intensity, world: World):
        super().__init__(config, 'Irradiation', dt)
        self.irradiation_intensity = irradiation_intensity

        #check if the file exists
        if not os.path.isfile(topas_file + '.csv'):
            #if it does not exist, run the simulation
            print('Running Topas simulation')
            print('Topas file: ', topas_file)
            file_with_geom = world.topas_param_file(topas_file)
            print('Topas param file: ', file_with_geom)
            term.RunTopasSimulation(file_with_geom, cluster=self.config.running_on_cluster)
            os.rename('MyScorer.csv', topas_file + '.csv')

        #read the dose from the file
        _, read_doses = rw.DoseOnWorld(topas_file + '.csv')

        #transform numpy aray into a list
        read_doses = read_doses.tolist()

        self.doses = np.zeros(len(read_doses)) #create an array of zeros

        #added to compensate for some irregularities in the Topas simulation, TODO: fix this in the Topas simulation and remove this
        for i in range(len(read_doses)):
            self.doses[i] = 0.5*(read_doses[i] + read_doses[len(read_doses) - i - 1])

        world.update_dose(self.doses) #update the dose on the world

        # plot the simulation
        fig, ax = plt.subplots(1, 3, figsize=(18, 5))
        world.show_tumor_slice(ax[0], fig, 'dose', cmap='RdPu',refinement_level=2, slice='x', round_n=2, vmin=0.2, vmax=0.7)
        world.show_tumor_slice(ax[1], fig, 'dose', cmap='RdPu',refinement_level=2, slice='y', round_n=2, vmin=0.2, vmax=0.7)
        world.show_tumor_slice(ax[2], fig, 'dose', cmap='RdPu',refinement_level=2, slice='z', round_n=2, vmin=0.2, vmax=0.7)
        fig.suptitle('Dose (arb. units)', fontsize=20)
        plt.tight_layout()
        plt.savefig('dose.png', dpi=300)
        if self.config.running_on_cluster:
            plt.close()
        else:
            plt.show()
    @Process.timeit
    def __call__(self, world: World):
        for voxel in world.voxel_list:
            scaled_dose = self.doses[voxel.voxel_number]*self.irradiation_intensity
            if len(voxel.list_of_cells) > 0:
                print('Scaled dose: ', scaled_dose)
                for cell in voxel.list_of_cells:
                    #assume all cells get damaged the same way
                    damage = scaled_dose * cell.radiosensitivity() #compute the damage
                    cell.damage += damage
                    cell.damage = min(cell.damage, 1.0)

        for vessel in world.vasculature.list_of_vessels:
            path = vessel.path
            total_dose = 0

            for point in path: #compute the mean dose on the vessel
                current_voxel = world.find_voxel_number(point)
                total_dose += world.voxel_list[current_voxel].dose

            if len(path) > 0:
                mean_dose = total_dose / len(path)
                vessel.must_be_updated = True
            else:
                mean_dose = 0

            damage_vessel = mean_dose * self.irradiation_intensity * vessel.radiosensitivity()
            vessel.maturity -= damage_vessel
            if vessel.maturity < 0:
                vessel.maturity = 0

        return
