import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import math
from BasicPlots import *
from BasicGeometries import *
import networkx as nx


class Vessel:
    def __init__(self, path, radius, parent_id=None):
        self.path = np.array(path)
        self.radius = radius
        self.id = id(self)
        self.parent_id = parent_id
        self.children_ids = []
        self.step_size = 0.05
        self.in_growth = True

    def __iter__(self):
        return self

    def grow(self, vegf_gradient, pressure, steps=1, weight_direction=0.5, weight_vegf=0.5, pressure_threshold=0.5,
                weight_pressure=0.5):
        for i in range(steps):
            self.step(vegf_gradient, pressure, weight_direction, weight_vegf, pressure_threshold, weight_pressure)

    def step(self, vegf_gradient, pressure, weight_direction=0.5, weight_vegf=0.5, pressure_threshold=0.5,
             weight_pressure=0.5):
            last_point = np.array(self.path[-1])
            prev_point = np.array(self.path[-2]) if len(self.path) > 1 else last_point
            direction = last_point - prev_point

            # Get VEGF gradient at the last point
            vegf_grad = np.array(vegf_gradient(last_point))

            local_pressure = pressure(last_point)

            # Normalize direction and VEGF gradient
            direction_norm = direction / np.linalg.norm(direction) if np.linalg.norm(direction) != 0 else direction
            vegf_grad_norm = vegf_grad / np.linalg.norm(vegf_grad) if np.linalg.norm(vegf_grad) != 0 else vegf_grad

            # Calculate the weighted direction based on vessel direction and VEGF gradient
            weighted_dir = weight_direction * direction_norm + weight_vegf * vegf_grad_norm

            if local_pressure > pressure_threshold:
                pressure_vec = -weighted_dir
                weighted_dir += weight_pressure * pressure_vec

                # Add random noise
            noise = np.array([random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)])
            weighted_dir += noise

            # Normalize the weighted direction
            weighted_dir /= np.linalg.norm(weighted_dir)

            # Calculate the new point and add it
            new_point = last_point + self.step_size * weighted_dir
            self.path = np.append(self.path, [new_point], axis=0)
    def volume_per_point(self):
        return np.pi * self.radius ** 2 * self.step_size

    def plot(self, ax, color='blue'):
        ax.plot(self.path[:, 0], self.path[:, 1], self.path[:, 2], color=color, alpha=0.5, linewidth= self.radius)

    def choose_random_point(self):
        #choose random point on the path, not the first or last point
        if len(self.path) < 3:
            return
        return random.choice(self.path[1:-1])
class VasculatureNetwork:
    def __init__(self, list_of_vessels=None):
        if list_of_vessels is None:
            list_of_vessels = []
        self.list_of_vessels = list_of_vessels

    def branching(self, vessel_id, branching_point):
        if branching_point == [] or branching_point is None:
            print("No branching point found")
            return
        mother_vessel = self.get_vessel(vessel_id)
        if mother_vessel is None:
            raise ValueError("Vessel with id {} does not exist".format(vessel_id))
        # split the path in two parts
        path = mother_vessel.path
        split_index = np.where((path == branching_point))[0][0]
        path1, path2 = np.split(path, [split_index])
        #remove the first element of path2
        path2 = np.delete(path2, 0, 0)
        mother_vessel.path = path1  # reduce the mother vessel path. It stops growing
        mother_vessel.in_growth = False
        # create two new vessels
        vessel1 = Vessel(path2, mother_vessel.radius, mother_vessel.id)
        #choose random point around branching point
        #random_point = branching_point + np.array([random.uniform(-mother_vessel.step, mother_vessel.step), random.uniform(-mother_vessel.step, mother_vessel.step), random.uniform(-mother_vessel.step, mother_vessel.step)])
        # the radius has to be updated later
        vessel2 = Vessel([branching_point], mother_vessel.radius, mother_vessel.id)  # the radius has to be updated later
        self.list_of_vessels.append(vessel1)
        self.list_of_vessels.append(vessel2)
        # update the mother vessel
        mother_vessel.children_ids.append(vessel1.id)
        mother_vessel.children_ids.append(vessel2.id)

    def add_vessel(self, vessel: Vessel):
        self.list_of_vessels.append(vessel)

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
        for child_id in vessel.children_ids:
            self.kill_vessel(child_id)
        # remove from the parent list of id
        parent_vessel = self.get_vessel(vessel.parent_id)
        if parent_vessel is not None:
            parent_vessel.children_ids.remove(vessel_id)
        # remove from the list of vessels
        self.list_of_vessels.remove(vessel)

    def update_vessels_radius(self, final_radius):
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

        # Call update_radius_recursive for each root vessel
        for root_vessel in root_vessels:
            update_radius_recursive(root_vessel.id)

    def volume_occupied(self):
        points = []
        volume = []
        for vessel in self.list_of_vessels:
            for point in vessel.path:
                points.append(point)
                volume.append(vessel.volume_per_point())
        return points, volume

    def grow(self, vegf_gradient, pressure, steps=1, weight_direction=0.5, weight_vegf=0.5, pressure_threshold=0.5,
             weight_pressure=0.5):
        for vessel in self.list_of_vessels:
            if vessel.in_growth:
                vessel.grow(vegf_gradient, pressure, steps, weight_direction, weight_vegf, pressure_threshold, weight_pressure)

    def plot(self, ax, color='blue'):
        for vessel in self.list_of_vessels:
            vessel.plot(ax, color)


    def grow_and_split(self,splitting_rate, vegf_gradient, pressure, macro_steps=1, micro_steps=10, weight_direction=0.5, weight_vegf=0.5, pressure_threshold=0.5, weight_pressure=0.5):
        dt = 10
        micro_steps = micro_steps*dt
        for i in range(macro_steps):
            print("Macro step {}".format(i))
            for vessel in self.list_of_vessels:
                print('current number of vessels {}'.format(len(self.list_of_vessels)))
                if vessel.in_growth:
                    vessel.grow(vegf_gradient, pressure, micro_steps, weight_direction, weight_vegf, pressure_threshold, weight_pressure)
                    if vessel.path.shape[0] > 1:
                        if random.uniform(0, 1) < splitting_rate*dt:
                            branching_point = vessel.choose_random_point()
                            self.branching(vessel.id, branching_point)

    def print_vessel_tree_recursive(self, vessels, children_ids, indent):
        for child_id in children_ids:
            child_vessel = next((v for v in vessels if v.id == child_id), None)
            if child_vessel is not None:
                print(' ' * indent, f"ID: {child_vessel.id}  Radius: {child_vessel.radius:.2f}")
                if child_vessel.children_ids:
                    self.print_vessel_tree_recursive(vessels, child_vessel.children_ids, indent + 2)


    def print_vessel_tree(self, indent=0):
        vessels = self.list_of_vessels
        # get the root vessels (i.e. vessels with no parent)
        root_vessels = [v for v in vessels if v.parent_id is None]

        # recursively print the vessel tree
        for root_vessel in root_vessels:
            print(' ' * indent, f"ID: {root_vessel.id}  Radius: {root_vessel.radius:.2f}")
            if root_vessel.children_ids:
                self.print_vessel_tree_recursive(vessels, root_vessel.children_ids, indent + 2)



