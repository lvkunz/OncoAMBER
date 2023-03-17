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


    def show(self, world: World, t = 0):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(90, 0)
        world.show_voxels_centers(ax, fig)
        world.vasculature.plot(fig, ax)
        plt.title('Cells in voxels at time t = ' + str(t) + ' hours')
        plt.savefig('Plots/Video/t'+ str(t) + '.png')
        plt.show()

        figA = plt.figure()
        axA = figA.add_subplot(111, projection='3d')
        axA.view_init(90, 0)
        world.show_voxels_centers_molecules(axA, figA, 'VEGF')
        world.vasculature.plot(figA, axA)
        plt.title('VEGF in voxels at time t = ' + str(t) + ' hours')
        plt.savefig('Plots/Video/t'+ str(t) + '_VEGF.png')
        plt.show()

        figB = plt.figure()
        axB = figB.add_subplot(111, projection='3d')
        axB.view_init(90, 0)
        world.show_voxels_centers_oxygen(axB, figB)
        world.vasculature.plot(figB, axB)
        plt.title('Oxygen in voxels at time t = ' + str(t) + ' hours')
        plt.savefig('Plots/Video/t'+ str(t) + '_Oxygen.png')
        plt.show()

        figC = plt.figure()
        axC = figC.add_subplot(111, projection='3d')
        axC.view_init(90, 0)
        world.show_voxels_centers_pressure(axC, figC)
        world.vasculature.plot(figC, axC)
        plt.title('Pressure in voxels at time t = ' + str(t) + ' hours')
        plt.savefig('Plots/Video/t'+ str(t) + '_Pressure.png')
        plt.show()

    def run(self, world: World, video = False):
        print('Running simulation for {} hours'.format(self.finish_time), ' with dt = ', self.dt)
        process_local = [process for process in self.list_of_process if (not process.is_global)]
        process_global = [process for process in self.list_of_process if process.is_global]

        if video: self.show(world, 0)

        while self.time <= self.finish_time:
            print('Time: {} hours'.format(self.time) + ' / ' + str(self.finish_time) + ' hours')

            #loop in random order
            copy_voxel_list = world.voxel_list.copy()
            np.random.shuffle(copy_voxel_list)

            for voxel in copy_voxel_list:
                for process in process_local:
                    process(voxel)
            for process in process_global:
                process(world)
            if video: self.show(world, self.time)
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
    def __init__(self, name, dt):
        super().__init__(name, dt)

    def __call__(self, voxel):
        #print('CellDivision')
        for cell in voxel.list_of_cells:
            if cell.state == 'cycling':
                # Calculate the expected number of cell divisions in the time step
                expected_divisions = self.dt / cell.doubling_time
                # Use a Poisson distribution to model the number of cell divisions
                num_divisions = np.random.poisson(expected_divisions)
                for i in range(num_divisions):
                    new_cell = cell.duplicate()
                    voxel.add_cell(new_cell)
        return
class CellApoptosis(Process):
    def __init__(self, name, dt):
        super().__init__('CellDeath', dt)

    def __call__(self, voxel):
        #print('CellDeath')
        # for cell in voxel.list_of_cells:
        #     effective_age = cell.age + self.dt/2
        #     probability = 1 - np.exp(-effective_age / cell.life_expectancy)
        #     #print(probability)
        #     if np.random.random() < probability:
        #         voxel.remove_cell(cell)

        # for cell in voxel.list_of_cells:
        #     # Calculate the expected number of cell deaths in the time step
        #     expected_deaths = (self.dt / cell.life_expectancy) * len(voxel.list_of_cells)
        #     # Use a Poisson distribution to model the number of cell deaths
        #     num_deaths = np.random.poisson(expected_deaths)
        #     num_deaths = min(num_deaths, len(voxel.list_of_cells))
        #     #print('num_deaths ', num_deaths)
        #     # Remove cells from the voxel at random
        #     cells_to_remove = np.random.choice(voxel.list_of_cells, size=num_deaths, replace=False)
        #     for cell in cells_to_remove:
        #         voxel.remove_cell(cell)

        for cell in voxel.list_of_cells:
            if cell.state == 'dead':
                voxel.remove_cell(cell)

class CellAging(Process):
    def __init__(self, name, dt):
        super().__init__('CellAging', dt)
    def __call__(self, voxel):
        #print('CellAging')
        for cell in voxel.list_of_cells:
            cell.age = cell.age + self.dt

class CellMigration(Process):
    def __init__(self, name, dt):
        super().__init__('CellMigration', dt)
        self.is_global = True
    def __call__(self, world : World): #maybe make it cleaner
        #print('CellMigration')
        #print('voxel', voxel.voxel_number)
        exchange_matrix = world.compute_exchange_matrix(self.dt)
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
    def __init__(self, name, dt):
        super().__init__('UpdateState', dt)
    def __call__(self, voxel: Voxel):
        #print('UpdateState')
        voxel.update_cells_for_oxygen_state()

class UpdateVoxelMolecules(Process):
    def __init__(self, name, dt):
        super().__init__('UpdateMolecules', dt)
    def __call__(self, voxel: Voxel):
        #print('UpdateMolecules')
        voxel.update_molecules(self.dt)
class UpdateVasculature(Process):
    def __init__(self, name, dt):
        super().__init__('UpdateVasculature', dt)
        self.is_global = True
    def __call__(self, world: World):
        #print('UpdateVasculature')
        n = world.vessels_killed_by_pressure()
        world.vasculature_growth(n)
        world.compute_oxygen_map()