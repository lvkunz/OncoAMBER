import random

import numpy as np

from amber.BasicGeometries import *
import sys

sys.setrecursionlimit(1500)

class Vessel:
    def __init__(self, path, radius, step_size, parent_id=None, children_ids=None, in_growth=True, intra_radiosensitivity=0.00001):
        if children_ids is None:
            children_ids = []
        for point in path:
            if not isinstance(point, np.ndarray) or len(point.shape) != 1 or point.shape[0] != 3:
                raise ValueError("Each point in the 'path' list must be a 3D array with shape (3,)")
        self.path = np.array(path)
        self.radius = radius
        self.id = id(self)
        self.parent_id = parent_id
        self.children_ids = children_ids
        self.step_size = step_size
        self.in_growth = in_growth
        self.visible = True
        self.must_be_updated = False
        self.intra_radiosensitivity = intra_radiosensitivity

    def __iter__(self):
        return self

    def length(self):
        return np.sum(np.linalg.norm(np.diff(self.path, axis=0), axis=1))

    def grow(self, limit_half_length, lower_pressure_threshold_step, higher_pressure_threshold_step, step_stop_threshold, vegf_gradient, pressure, steps, weight_direction, weight_vegf, weight_pressure):
        length = 0
        for i in range(steps):
            step_size = self.step(limit_half_length, lower_pressure_threshold_step, higher_pressure_threshold_step, vegf_gradient, pressure, weight_direction, weight_vegf, weight_pressure)
            if step_size < step_stop_threshold:
                self.in_growth = False
                break
            length += step_size
        return length
        # if step_size < config.growth_stop_threshold'] * self.step_size:
        #     self.in_growth = False

    def step(self, half_length_world, lower_pressure_threshold_step, higher_pressure_threshold_step, vegf_gradient, pressure, weight_direction=0.5, weight_vegf=0.5, weight_pressure=0.5):
        if not isinstance(self.path, np.ndarray) or len(self.path.shape) != 2 or self.path.shape[1] != 3:
            raise ValueError("The 'self.path' array must be a 3D array with shape (n, 3)")
        last_point = np.array(self.path[-1])
        first_point = np.array(self.path[0])

        if last_point[0] < -half_length_world or last_point[0] > half_length_world or \
                last_point[1] < -half_length_world or last_point[1] > half_length_world or \
                last_point[2] < -half_length_world or last_point[2] > half_length_world:
            return 0

        direction = last_point - first_point

        # Get VEGF gradient at the last point
        vegf_grad = np.array(vegf_gradient(last_point))
        local_pressure = pressure(last_point)

        # Normalize direction and VEGF gradient
        direction_norm = direction / np.linalg.norm(direction) if np.linalg.norm(direction) != 0 else direction
        vegf_grad_norm = vegf_grad / np.linalg.norm(vegf_grad) if np.linalg.norm(vegf_grad) != 0 else vegf_grad

        # Calculate the weighted direction based on vessel direction and VEGF gradient
        weighted_dir = weight_direction * direction_norm + weight_vegf * vegf_grad_norm

        pressure_vec = -weighted_dir
        weighted_dir += weight_pressure * pressure_vec * local_pressure

        # Add random noise
        noise = np.array([random.random()*2 - 1, random.random()*2-1, random.random()*2 - 1])
        weighted_dir += noise

        # Normalize the weighted direction
        weighted_dir /= np.linalg.norm(weighted_dir)

        # Calculate the new point and add it
        step_size = self.step_size

        if weight_pressure > 0 and local_pressure > lower_pressure_threshold_step:
            step_size = step_size * (1 - (local_pressure - lower_pressure_threshold_step) / higher_pressure_threshold_step)

        new_point = last_point + step_size * weighted_dir
        self.path = np.append(self.path, [new_point], axis=0)
        return step_size

    def volume_per_point(self):
        return np.pi * self.radius ** 2 * self.step_size

    def plot(self,fig, ax, color='crimson'):
        if self.visible:
            if self.in_growth:
                ax.plot(self.path[:, 0], self.path[:, 1], self.path[:, 2], color=color, alpha=0.7, linewidth= self.radius*100)
            else:
                ax.plot(self.path[:, 0], self.path[:, 1], self.path[:, 2], color='mediumblue', alpha=0.1, linewidth= self.radius*100)
        return fig, ax

    def choose_random_point(self, seed):
        # choose random point on the path, not the first or last point
        if len(self.path) < 3:
            print(self.path)
            raise ValueError("The vessel path is too short to choose a random point")
        random_index = random.randint(1, len(self.path) - 2)
        return self.path[random_index]

    def mean_pressure(self, pressure):
        if len(self.path) < 2:
            return 0
        else:
            return np.mean([pressure(point) for point in self.path])

    def max_pressure(self, pressure):
        return np.max([pressure(point) for point in self.path])

    def max_dose(self, dose_map):
        max_dose = 0
        for point in self.path:
            if dose_map.evaluate(point) > max_dose:
                max_dose = dose_map.evaluate(point)
        return max_dose
    def radiosensitivity(self):
        radius = self.radius
        radiosensitivity = self.intra_radiosensitivity/radius
        return radiosensitivity
class VasculatureNetwork:
    def __init__(self, config, list_of_vessels=None):
        if list_of_vessels is None:
            list_of_vessels = []
        self.list_of_vessels = list_of_vessels
        self.config = config

    def branching(self, vessel_id, branching_point):
        if branching_point == [] or branching_point is None:
            print("No branching point found, branching point is", branching_point)
            return
        mother_vessel = self.get_vessel(vessel_id)
        if mother_vessel is None:
            raise ValueError("Vessel with id {} does not exist".format(vessel_id))
        # split the path in two parts
        path = mother_vessel.path
        if len(path) < 3:
            raise ValueError("Vessel with id {} has a path with less than 3 points".format(vessel_id))
        visible = mother_vessel.visible
        split_index = np.where((path == branching_point))[0][0]
        path_begin, path_end = np.split(path, [split_index])
        #remove the first element of path_end to keep the same point distribution in space
        path_end = np.delete(path_end, 0, 0)
        mother_vessel.path = path_begin  # reduce the mother vessel path. It stops growing
        mother_vessel.in_growth = False

        # print(mother_vessel)
        # print('radio', mother_vessel.intra_radiosensitivity)

        # create two new vessels
        vessel_end = Vessel(path_end, mother_vessel.radius, self.config.vessel_step_size, parent_id= mother_vessel.id, children_ids=mother_vessel.children_ids, in_growth= mother_vessel.in_growth, intra_radiosensitivity= mother_vessel.intra_radiosensitivity)  # the radius has to be updated later
        for child_id in vessel_end.children_ids:
            child = self.get_vessel(child_id)
            child.parent_id = vessel_end.id
        #choose random point around branching point
        #random_point = branching_point + np.array([random.uniform(-mother_vessel.step, mother_vessel.step), random.uniform(-mother_vessel.step, mother_vessel.step), random.uniform(-mother_vessel.step, mother_vessel.step)])
        # the radius has to be updated later
        mother_vessel.visible = visible
        vessel_end.visible = visible

        vessel_new = Vessel([branching_point], self.config.radius_root_vessels, self.config.vessel_step_size, parent_id= mother_vessel.id, children_ids=None, in_growth=True, intra_radiosensitivity=mother_vessel.intra_radiosensitivity)  # the radius has to be updated later
        vessel_new.visible = True
        mother_vessel.children_ids = [vessel_end.id, vessel_new.id]
        self.list_of_vessels.append(vessel_end)
        self.list_of_vessels.append(vessel_new)

        #find the root vessel:
        root_vessel = mother_vessel
        while root_vessel.parent_id is not None:
            root_vessel = self.get_vessel(root_vessel.parent_id)

        root_vessel.must_be_updated = True
        # update the mother vessel


    def add_vessel(self, vessel: Vessel):
        self.list_of_vessels.append(vessel)

    def add_multiple_vessels(self, vessels):
        for vessel in vessels:
            self.add_vessel(vessel)

    def get_vessel(self, vessel_id):
        for vessel in self.list_of_vessels:
            if vessel.id == vessel_id:
                return vessel
        return None

    def kill_vessel(self, vessel_id):
        vessel = self.get_vessel(vessel_id)
        if vessel is None:
            raise ValueError("Vessel with id {} does not exist".format(vessel_id))
        # remove the children
        for child_id in list(vessel.children_ids):
            self.kill_vessel(child_id)
        # remove from the parent list of id
        parent_vessel = self.get_vessel(vessel.parent_id)
        if parent_vessel is not None:
            if vessel_id not in parent_vessel.children_ids:
                print("Vessel {} is not in the list of children of vessel {}".format(vessel_id, parent_vessel.id))
            parent_vessel.children_ids.remove(vessel_id)
        # remove from the list of vessels
        #root vessel to be updated
        root_vessel = vessel
        while root_vessel.parent_id is not None:
            root_vessel = self.get_vessel(root_vessel.parent_id)
        root_vessel.must_be_updated = True
        self.list_of_vessels.remove(vessel)

    def update_vessels_radius_from_last(self, final_radius, pressure_sensitive, pressure):
        print("Updating vessels radius from last")
        def update_radius_recursive(vessel_id):
            vessel = self.get_vessel(vessel_id)
            if not vessel.children_ids:
                vessel.radius = final_radius
            else:
                for child_id in vessel.children_ids:
                    update_radius_recursive(child_id)

                r_cubed_sum = sum([self.get_vessel(child_id).radius**3 for child_id in vessel.children_ids])
                vessel.radius = r_cubed_sum**(1/3)

        # Find all root vessels with parent_id=None
        root_vessels = [v for v in self.list_of_vessels if v.parent_id is None]
        root_vessels = [v for v in root_vessels if v.must_be_updated is True]
        print("Number of root vessels to update: ", len(root_vessels))

        # Call update_radius_recursive for each root vessel
        for root_vessel in root_vessels:
            update_radius_recursive(root_vessel.id)
            root_vessel.must_be_updated = False

        if pressure_sensitive:
            print("Updating for pressure")
            pow = self.config.radius_decrease_exponent
            coeff = self.config.max_occupancy ** (-pow)
            for vessel in self.list_of_vessels:
                vessel.radius = vessel.radius * (1 - coeff * (vessel.mean_pressure(pressure)**pow))

    def volume_occupied(self):
        points = []
        volume = []
        for vessel in self.list_of_vessels:
            for point in vessel.path:
                points.append(point)
                volume.append(vessel.volume_per_point())
        return points, volume

    def grow(self, vegf_gradient, pressure, steps, weight_direction, weight_vegf, weight_pressure):
        for vessel in self.list_of_vessels:
            if vessel.in_growth:
                vessel.grow(self.config.half_length_world, self.config.lower_pressure_threshold_step, self.config.higher_pressure_threshold_step, self.config.step_stop_threshold, vegf_gradient, pressure, steps, weight_direction, weight_vegf, weight_pressure)
    def plot(self, fig, ax, color='crimson'):
        for vessel in self.list_of_vessels:
            vessel.plot(fig, ax, color)
        return fig, ax


    def grow_and_split(self, dt, splitting_rate, vegf_gradient, pressure, macro_steps, micro_steps, weight_direction, weight_vegf, weight_pressure):
        micro_steps = micro_steps

        macro_steps_ = int(macro_steps * dt)

        if macro_steps_ == 0:
            ValueError("Macro steps must be at least 1")
            if macro_steps == 0:
                ValueError("Angiogenesis deactivation")
                macro_steps_ = 0
            else:
                macro_steps_ = 1

        for i in range(macro_steps_):
            print("Macro step {}".format(i))
            j = 0
            for vessel in self.list_of_vessels:
                splitting_rate_ = splitting_rate
                if j % 1000 == 0: print('current number of vessels {}'.format(len(self.list_of_vessels)))
                if vessel.in_growth:
                    # total_path_length = vessel.step_size * micro_steps
                    new_path_length = vessel.grow(self.config.half_length_world, self.config.lower_pressure_threshold_step, self.config.higher_pressure_threshold_step, self.config.growth_step_stop_threshold, vegf_gradient, pressure, micro_steps, weight_direction, weight_vegf, weight_pressure)
                    if len(vessel.path) > 3:
                        splitting_rate_ = vessel.length() * splitting_rate_
                        if random.uniform(0, 1) < splitting_rate_:
                            branching_point = vessel.choose_random_point(self.config.seed)
                            self.branching(vessel.id, branching_point)
                j += 1
    def print_vessel_tree_recursive(self, vessels, children_ids, indent):
        for child_id in children_ids:
            child_vessel = next((v for v in vessels if v.id == child_id), None)
            if child_vessel is not None:
                print(' ' * indent, f"ID: {child_vessel.id}  Radius: {child_vessel.radius:.5f}")
                if child_vessel.children_ids:
                    self.print_vessel_tree_recursive(vessels, child_vessel.children_ids, indent + 2)


    def print_vessel_tree(self, indent=0):
        vessels = self.list_of_vessels
        # get the root vessels (i.e. vessels with no parent)
        root_vessels = [v for v in vessels if v.parent_id is None]

        # recursively print the vessel tree
        for root_vessel in root_vessels:
            print(' ' * indent, f"ID: {root_vessel.id}  Radius: {root_vessel.radius:.5f}")
            if root_vessel.children_ids:
                self.print_vessel_tree_recursive(vessels, root_vessel.children_ids, indent + 2)

    def compute_VSL(self):
        list_VSL = np.array([])
        for vessel in self.list_of_vessels:
            if vessel.visible:
                list_VSL = np.append(list_VSL, vessel.length())
        return list_VSL

    def compute_diameters(self):
        list_diameters = np.array([])
        for vessel in self.list_of_vessels:
            if vessel.visible:
                list_diameters = np.append(list_diameters, vessel.radius*2)
        return list_diameters
