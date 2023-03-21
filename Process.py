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

        Vasculature_show = False


        if slice:
            angle = 0
            angle2 = 0

        else:
            angle = 30
            angle2 = 60

        DPI = 100

        fig = plt.figure()
        #set fig size
        fig.set_size_inches(10, 10)
        fig.set_dpi(DPI)
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(angle, angle2)
        world.show_voxels_centers(ax, fig, slice = slice)
        if Vasculature_show: world.vasculature.plot(fig, ax)
        plt.title('Cells in voxels at time t = ' + str(t) + ' hours')
        plt.savefig('Plots/Video/t'+ str(t) + '.png')
        plt.show()

        figA = plt.figure()
        figA.set_size_inches(10, 10)
        figA.set_dpi(DPI)
        axA = figA.add_subplot(111, projection='3d')
        axA.view_init(angle, angle2)
        world.show_voxels_centers_molecules(axA, figA, 'VEGF', slice = slice)
        if Vasculature_show: world.vasculature.plot(figA, axA)
        plt.title('VEGF in voxels at time t = ' + str(t) + ' hours')
        plt.savefig('Plots/Video/t'+ str(t) + '_VEGF.png')
        plt.show()

        figB = plt.figure()
        figB.set_size_inches(10, 10)
        figB.set_dpi(DPI)
        axB = figB.add_subplot(111, projection='3d')
        axB.view_init(angle, angle2)
        world.show_voxels_centers_oxygen(axB, figB, slice = slice)
        if Vasculature_show: world.vasculature.plot(figB, axB)
        plt.title('Oxygen in voxels at time t = ' + str(t) + ' hours')
        plt.savefig('Plots/Video/t'+ str(t) + '_Oxygen.png')
        plt.show()

        figC = plt.figure()
        figC.set_size_inches(10, 10)
        figC.set_dpi(DPI)
        axC = figC.add_subplot(111, projection='3d')
        axC.view_init(angle, angle2)
        world.show_voxels_centers_pressure(axC, figC, slice = slice)
        if Vasculature_show: world.vasculature.plot(figC, axC)
        plt.title('Pressure in voxels at time t = ' + str(t) + ' hours')
        plt.savefig('Plots/Video/t'+ str(t) + '_Pressure.png')
        plt.show()

    def run(self, world: World, video = False):
        print('Running simulation for {} hours'.format(self.finish_time), ' with dt = ', self.dt)
        process_local = [process for process in self.list_of_process if (not process.is_global)]
        process_global = [process for process in self.list_of_process if process.is_global]

        while self.time <= self.finish_time:
            print('Time: {} hours'.format(self.time) + ' / ' + str(self.finish_time) + ' hours')

            #loop in random order
            copy_voxel_list = world.voxel_list.copy()
            np.random.shuffle(copy_voxel_list)
            for process in process_global:
                print('Currently running global process : ', process.name)
                process(world)
            for voxel in copy_voxel_list:
                for process in process_local:
                    process(voxel)
            if video: self.show(world, self.time, slice = True)
            self.time = self.time + self.dt
        return


class Process:
    def __init__(self, name, dt):
        self.name = name
        self.dt = dt
        self.is_global = False
    def __call__(self, voxel):
        pass


class CellDivision(Process):
    def __init__(self, name, dt, cycling_threshold):
        super().__init__(name, dt)
        self.cycling_threshold = cycling_threshold

    def __call__(self, voxel):
        #print('CellDivision')
        for cell in voxel.list_of_cells:
            if cell.vitality > self.cycling_threshold:
                # Calculate the expected number of cell divisions in the time step
                expected_divisions = self.dt / cell.doubling_time
                # Use a Poisson distribution to model the number of cell divisions
                num_divisions = np.random.poisson(expected_divisions)
                for i in range(num_divisions):
                    new_cell = cell.duplicate()
                    voxel.add_cell(new_cell)
        return
class CellApoptosis(Process):
    def __init__(self, name, dt, apoptosis_threshold):
        super().__init__('CellDeath', dt)
        self.apoptosis_threshold = apoptosis_threshold

    def __call__(self, voxel):
        for cell in voxel.list_of_cells:
            if cell.vitality < self.apoptosis_threshold :
                voxel.remove_cell(cell)

class CellAging(Process):
    def __init__(self, name, dt):
        super().__init__('CellAging', dt)
    def __call__(self, voxel):
        #print('CellAging')
        for cell in voxel.list_of_cells:
            cell.age = cell.age + self.dt

class CellMigration(Process):
    def __init__(self, name, dt, pressure_threshold):
        super().__init__('CellMigration', dt)
        self.is_global = True
        self.pressure_threshold = pressure_threshold
    def __call__(self, world : World): #maybe make it cleaner
        #print('CellMigration')
        #print('voxel', voxel.voxel_number)
        exchange_matrix = world.compute_exchange_matrix(self.dt, pressure_threshold=self.pressure_threshold)
        for voxel in world.voxel_list:
            list_of_neighbors = world.find_neighbors(voxel)
            for cell in voxel.list_of_cells:
                np.random.shuffle(list_of_neighbors)
                for neighbor in list_of_neighbors:
                    probability = exchange_matrix[voxel.voxel_number, neighbor.voxel_number]
                    if np.random.random() < probability:
                        neighbor.add_cell(cell)
                        voxel.remove_cell(cell)
                        break


class UpdateCellState(Process):
    def __init__(self, name, dt, vitality_o2_threshold, vitality_slope):
        super().__init__('UpdateState', dt)
        self.vitality_o2_threshold = vitality_o2_threshold
        self.vitality_slope = vitality_slope

    def __call__(self, voxel: Voxel):
        #print('UpdateState')
        voxel.update_cells_for_oxygen_state(vitality_slope=self.vitality_slope, vitality_threshold=self.vitality_o2_threshold)
class UpdateVoxelMolecules(Process):
    def __init__(self, name, dt, VEGF_production_per_cell):
        super().__init__('UpdateMolecules', dt)
        self.VEGF_production_per_cell = VEGF_production_per_cell
    def __call__(self, voxel: Voxel):
        #print('UpdateMolecules')
        voxel.update_molecules(self.dt, self.VEGF_production_per_cell)
class UpdateVasculature(Process):
    def __init__(self, name, dt, pressure_killing_threshold, pressure_killing_slope, vasculature_growth_factor, o2_per_volume, diffusion_number):
        super().__init__('UpdateVasculature', dt)
        self.is_global = True
        self.pressure_killing_threshold = pressure_killing_threshold
        self.pressure_killing_slope = pressure_killing_slope
        self.vasculature_growth_factor = vasculature_growth_factor
        self.o2_per_volume = o2_per_volume
        self.diffusion_number = diffusion_number
    def __call__(self, world: World):
        #print('UpdateVasculature')
        print('Killing vessels')
        n = world.vessels_killed_by_pressure(pressure_killing_threshold=self.pressure_killing_threshold, pressure_killing_slope=self.pressure_killing_slope)
        print('Killed vessels: ', n)
        print('Growing vessels')
        world.vasculature_growth(n * self.vasculature_growth_factor)
        print('Computing oxygen map')
        world.compute_oxygen_map(o2_per_volume = self.o2_per_volume, diffusion_number=self.diffusion_number)