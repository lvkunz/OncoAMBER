from World import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Cell import *
from Voxel import *

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

        DPI = 200
        size = world.half_length

        fig = plt.figure()
        # set fig size
        fig.set_size_inches(10, 10)
        fig.set_dpi(DPI)
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(angle, angle2)
        ax.set_xlim(-size, size)
        ax.set_ylim(-size, size)
        ax.set_zlim(-size, size)
        world.show_tumor(ax, fig, slice = slice)
        if Vasculature_show: world.vasculature.plot(fig, ax)
        plt.title('Cells in voxels at time t = ' + str(t) + ' hours')
        plt.savefig('Plots/Video/t' + str(t) + '.png')
        plt.show()



        # fig = plt.figure()
        # #set fig size
        # fig.set_size_inches(10, 10)
        # fig.set_dpi(DPI)
        # ax = fig.add_subplot(111, projection='3d')
        # ax.view_init(angle, angle2)
        # world.show_voxels_centers(ax, fig, slice = slice)
        # if Vasculature_show: world.vasculature.plot(fig, ax)
        # plt.title('Cells in voxels at time t = ' + str(t) + ' hours')
        # plt.savefig('Plots/Video/t'+ str(t) + '.png')
        # plt.show()

        # figA = plt.figure()
        # figA.set_size_inches(10, 10)
        # figA.set_dpi(DPI)
        # axA = figA.add_subplot(111, projection='3d')
        # axA.view_init(angle, angle2)
        # #world.show_voxels_centers_molecules(axA, figA, 'VEGF', slice = slice)
        # print('de')
        # vegf_map = world.vegf_map()
        # print('de2')
        # #vegf_map.show_gradient(figA, axA, length_scale=0.2)
        # print('de3')
        # if Vasculature_show: world.vasculature.plot(figA, axA)
        # plt.title('VEGF in voxels at time t = ' + str(t) + ' hours')
        # plt.savefig('Plots/Video/t'+ str(t) + '_VEGF.png')
        # plt.show()
        #
        # figB = plt.figure()
        # figB.set_size_inches(10, 10)
        # figB.set_dpi(DPI)
        # axB = figB.add_subplot(111, projection='3d')
        # axB.view_init(angle, angle2)
        # world.show_voxels_centers_oxygen(axB, figB, slice = slice)
        # if Vasculature_show: world.vasculature.plot(figB, axB)
        # plt.title('Oxygen in voxels at time t = ' + str(t) + ' hours')
        # plt.savefig('Plots/Video/t'+ str(t) + '_Oxygen.png')
        # plt.show()
        #
        # figC = plt.figure()
        # figC.set_size_inches(10, 10)
        # figC.set_dpi(DPI)
        # axC = figC.add_subplot(111, projection='3d')
        # axC.view_init(angle, angle2)
        # world.show_voxels_centers_pressure(axC, figC, slice = slice)
        # if Vasculature_show: world.vasculature.plot(figC, axC)
        # plt.title('Pressure in voxels at time t = ' + str(t) + ' hours')
        # plt.savefig('Plots/Video/t'+ str(t) + '_Pressure.png')
        # plt.show()

        # figO = plt.figure()
        # axO = figO.add_subplot(111, projection='3d')
        # axO.set_xlim([-world.half_length, world.half_length])
        # axO.set_ylim([-world.half_length, world.half_length])
        # axO.set_zlim([-world.half_length, world.half_length])
        # world.vasculature.plot(figO, axO)
        # plt.title('Vasculature at time t = ' + str(t) + ' hours')
        # plt.show()

    def run(self, world: World, video = False):
        print('Running simulation for {} hours'.format(self.finish_time), ' with dt = ', self.dt)
        number_cells = []
        process_local = [process for process in self.list_of_process if (not process.is_global)]
        process_global = [process for process in self.list_of_process if process.is_global]

        self.show(world, self.time, slice = False)

        while self.time < self.finish_time:
            count_cells = 0
            print('Time: {} hours'.format(self.time) + ' / ' + str(self.finish_time) + ' hours')
            for process in process_global:
                print('Currently running global process : ', process.name)
                process(world)
            #loop in random order
            # copy_voxel_list = world.voxel_list.copy()
            # np.random.shuffle(copy_voxel_list)
            for voxel in world.voxel_list:
                count_cells += voxel.number_of_tumor_cells()
                for process in process_local:
                    process(voxel)
            self.time = self.time + self.dt
            if video: self.show(world, self.time, slice = False)
            number_cells.append(count_cells)
        print('Simulation finished')
        self.show(world, self.time, slice = True)
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
            print('pressure = ', voxel.pressure(), ' > ', self.pressure_threshold, ' so no cell division')
        return
class CellApoptosis(Process):
    def __init__(self, name, dt, apoptosis_threshold):
        super().__init__('CellDeath', dt)
        self.apoptosis_threshold = apoptosis_threshold

    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            if cell.vitality() < self.apoptosis_threshold:
                voxel.remove_cell(cell)

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
                        if neighbor.add_cell(cell):
                            voxel.remove_cell(cell)


class UpdateCellState(Process):
    def __init__(self, name, dt, spread_gaussian_o2):
        super().__init__('UpdateState', dt)
        self.spread_gaussian_o2 = spread_gaussian_o2

    def __call__(self, voxel: Voxel):
        voxel.update_cells_oxygen_state(self.spread_gaussian_o2)
class UpdateVoxelMolecules(Process):
    def __init__(self, name, dt, VEGF_production_per_cell, threshold_for_VEGF_production):
        super().__init__('UpdateMolecules', dt)
        self.VEGF_production_per_cell = VEGF_production_per_cell
        self.threshold_for_VEGF_production = threshold_for_VEGF_production
    def __call__(self, voxel: Voxel):
        voxel.update_vegf(self.dt, self.VEGF_production_per_cell, self.threshold_for_VEGF_production)
class UpdateVasculature(Process):
    def __init__(self, name, dt, pressure_killing_threshold, o2_per_volume, diffusion_number, splitting_rate, macro_steps, micro_steps, weight_direction, weight_vegf, pressure_threshold, weight_pressure, radius_pressure_sensitive):
        super().__init__('UpdateVasculature', dt)
        self.is_global = True
        self.pressure_killing_radius_threshold = pressure_killing_threshold
        self.o2_per_volume = o2_per_volume
        self.diffusion_number = diffusion_number
        self.dt = dt
        self.splitting_rate = splitting_rate
        self.macro_steps = macro_steps
        self.micro_steps = micro_steps
        self.weight_direction = weight_direction
        self.weight_vegf = weight_vegf
        self.pressure_threshold = pressure_threshold
        self.weight_pressure = weight_pressure
        self.radius_pressure_sensitive = radius_pressure_sensitive


    def __call__(self, world: World):
        print('Killing vessels')
        n_killed = world.vessels_killed_by_pressure()
        print('Killed vessels: ', n_killed)
        print('Growing vessels')
        total_VEGF = 0
        for voxel in world.voxel_list:
            total_VEGF += voxel.molecular_factors['VEGF']
        world.vasculature_growth(self.dt, self.splitting_rate, self.macro_steps, self.micro_steps, self.weight_direction, self.weight_vegf, self.pressure_threshold, self.weight_pressure, self.radius_pressure_sensitive)
        world.update_volume_occupied_by_vessels()
        world.update_oxygen(o2_per_volume = self.o2_per_volume, diffusion_number=self.diffusion_number)