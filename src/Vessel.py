import random
from src.BasicGeometries import *
import sys
from src.config_instance import config

sys.setrecursionlimit(1500)
rng = np.random.default_rng(config.seed)

class Vessel:
    def __init__(self, path, radius, parent_id=None, children_ids=None, in_growth=True):
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
        self.step_size = config.vessel_step_size
        self.in_growth = in_growth
        self.healthy = False
        self.visible = True

    def to_dict(self):
        return {
            "path": self.path.tolist(),
            "radius": self.radius,
            "id": self.id,
            "parent_id": self.parent_id,
            "children_ids": self.children_ids,
            "step_size": self.step_size,
            "in_growth": self.in_growth,
        }
    def __iter__(self):
        return self

    def length(self):
        return np.sum(np.linalg.norm(np.diff(self.path, axis=0), axis=1))

    def grow(self, vegf_gradient, pressure, steps=1, weight_direction=0.5, weight_vegf=0.5, weight_pressure=0.5):
        length = 0
        for i in range(steps):
            step_size = self.step(vegf_gradient, pressure, weight_direction, weight_vegf, weight_pressure)
            if step_size < config.growth_step_stop_threshold:
                self.in_growth = False
                break
            length += step_size
        return length
        # if step_size < config.growth_stop_threshold'] * self.step_size:
        #     self.in_growth = False

    def step(self, vegf_gradient, pressure, weight_direction=0.5, weight_vegf=0.5, weight_pressure=0.5):
        if not isinstance(self.path, np.ndarray) or len(self.path.shape) != 2 or self.path.shape[1] != 3:
            raise ValueError("The 'self.path' array must be a 3D array with shape (n, 3)")
        last_point = np.array(self.path[-1])
        prev_point = np.array(self.path[-2]) if len(self.path) > 1 else last_point

        if last_point[0] < -config.half_length_world or last_point[0] > config.half_length_world or \
                last_point[1] < -config.half_length_world or last_point[1] > config.half_length_world or \
                last_point[2] < -config.half_length_world or last_point[2] > config.half_length_world:
            return 0

        direction = last_point - prev_point

        # Get VEGF gradient at the last point
        vegf_grad = np.array(vegf_gradient(last_point))
        local_pressure = pressure(last_point)

        # Normalize direction and VEGF gradient
        direction_norm = direction / np.linalg.norm(direction) if np.linalg.norm(direction) != 0 else direction
        vegf_grad_norm = vegf_grad / np.linalg.norm(vegf_grad) if np.linalg.norm(vegf_grad) != 0 else vegf_grad
        vegf_grad_norm_scalar = np.linalg.norm(vegf_grad_norm)

        # Calculate the weighted direction based on vessel direction and VEGF gradient
        weighted_dir = weight_direction * direction_norm + weight_vegf * vegf_grad_norm

        pressure_vec = -weighted_dir
        weighted_dir += weight_pressure * pressure_vec * local_pressure

        # Add random noise
        noise = np.array([random.random()*2 - 1, random.random()*2-1, random.random()*2 - 1])
        weighted_dir += noise

        # Normalize the weighted direction
        weighted_dir /= np.linalg.norm(weighted_dir)

        if config.verbose:
            print('vegf_grad_norm_scalar: ', vegf_grad_norm_scalar)
            print('local_pressure: ', local_pressure)

        # Calculate the new point and add it
        step_size = self.step_size
        if weight_vegf > 0 and vegf_grad_norm_scalar < config.reference_vegf_gradient:
            step_size = step_size * (vegf_grad_norm_scalar / config.reference_vegf_gradient)
        if weight_pressure > 0 and local_pressure > config.reference_pressure:
            step_size = step_size * (config.reference_pressure / local_pressure)

        new_point = last_point + step_size * weighted_dir
        self.path = np.append(self.path, [new_point], axis=0)
        return step_size

    def volume_per_point(self):
        return np.pi * self.radius ** 2 * self.step_size

    def plot(self,fig, ax, color='crimson'):
        if self.visible:
            if self.in_growth:
                ax.plot(self.path[:, 0], self.path[:, 1], self.path[:, 2], color=color, alpha=0.7, linewidth= self.radius*300)
            else:
                ax.plot(self.path[:, 0], self.path[:, 1], self.path[:, 2], color='mediumblue', alpha=0.1, linewidth= self.radius*300)
        return fig, ax
    def choose_random_point(self):
        #choose random point on the path, not the first or last point
        if len(self.path) < 3:
            print(self.path)
            raise ValueError("The vessel path is too short to choose a random point")
        return rng.choice(self.path[1:-1])
    def mean_pressure(self, pressure):
        if len(self.path) < 2:
            return 0
        else:
            return np.mean([pressure(point) for point in self.path])

    def max_pressure(self, pressure):
        return np.max([pressure(point) for point in self.path])
class VasculatureNetwork:
    def __init__(self, list_of_vessels=None):
        if list_of_vessels is None:
            list_of_vessels = []
        self.list_of_vessels = list_of_vessels

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

        split_index = np.where((path == branching_point))[0][0]
        path_begin, path_end = np.split(path, [split_index])
        #remove the first element of path_end to keep the same point distribution in space
        path_end = np.delete(path_end, 0, 0)
        mother_vessel.path = path_begin  # reduce the mother vessel path. It stops growing
        mother_vessel.in_growth = False

        # create two new vessels
        vessel_end = Vessel(path_end, mother_vessel.radius, parent_id= mother_vessel.id, children_ids=mother_vessel.children_ids, in_growth= mother_vessel.in_growth)
        for child_id in vessel_end.children_ids:
            child = self.get_vessel(child_id)
            child.parent_id = vessel_end.id
        #choose random point around branching point
        #random_point = branching_point + np.array([random.uniform(-mother_vessel.step, mother_vessel.step), random.uniform(-mother_vessel.step, mother_vessel.step), random.uniform(-mother_vessel.step, mother_vessel.step)])
        # the radius has to be updated later

        vessel_new = Vessel([branching_point], config.radius_root_vessels, parent_id= mother_vessel.id)  # the radius has to be updated later
        mother_vessel.children_ids = [vessel_end.id, vessel_new.id]
        self.list_of_vessels.append(vessel_end)
        self.list_of_vessels.append(vessel_new)
        # update the mother vessel


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
        for child_id in list(vessel.children_ids):
            self.kill_vessel(child_id)
        # remove from the parent list of id
        parent_vessel = self.get_vessel(vessel.parent_id)
        if parent_vessel is not None:
            if vessel_id not in parent_vessel.children_ids:
                print("Vessel {} is not in the list of children of vessel {}".format(vessel_id, parent_vessel.id))
            parent_vessel.children_ids.remove(vessel_id)
        # remove from the list of vessels
        self.list_of_vessels.remove(vessel)

    #THIS FUNCTION HAS NOT BEEN TESTED
    def update_vessels_radius_from_last(self, final_radius, pressure_sensitive=False, pressure=None):
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

        # Call update_radius_recursive for each root vessel
        for root_vessel in root_vessels:
            update_radius_recursive(root_vessel.id)

        if pressure_sensitive:
            for vessel in self.list_of_vessels:
                vessel.radius = vessel.radius / ((1 + (vessel.mean_pressure(pressure)))**config.radius_decrease_exponent)

    def update_vessels_radius_from_root(self, root_radius, pressure_sensitive=False, pressure=None):
        print("Updating vessels radius from root")

        def update_radius_recursive(vessel_id):
            vessel = self.get_vessel(vessel_id)

            if vessel.children_ids:
                # Determine the scaling factor for the children based on Murray's law
                r_cubed_sum = sum([self.get_vessel(child_id).radius ** 3 for child_id in vessel.children_ids])
                scaling_factor = (vessel.radius ** 3 / r_cubed_sum) ** (1 / 3)

                for child_id in vessel.children_ids:
                    child_vessel = self.get_vessel(child_id)
                    child_vessel.radius *= scaling_factor
                    update_radius_recursive(child_id)

        # Find all root vessels with parent_id=None
        root_vessels = [v for v in self.list_of_vessels if v.parent_id is None]

        # Call update_radius_recursive for each root vessel
        for root_vessel in root_vessels:
            root_vessel.radius = root_radius
            update_radius_recursive(root_vessel.id)

        if pressure_sensitive:
            for vessel in self.list_of_vessels:
                vessel.radius = vessel.radius / (
                            (1 + (vessel.mean_pressure(pressure))) ** config.radius_decrease_exponent)

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

    def plot(self, fig, ax, color='crimson'):
        for vessel in self.list_of_vessels:
            vessel.plot(fig, ax, color)
        return fig, ax


    def grow_and_split(self, dt, splitting_rate, vegf_gradient, pressure, macro_steps=1, micro_steps=10, weight_direction=0.5, weight_vegf=0.5, weight_pressure=0.5):
        micro_steps = micro_steps
        for i in range(macro_steps * dt):
            print("Macro step {}".format(i))
            j = 0
            for vessel in self.list_of_vessels:
                splitting_rate_ = splitting_rate
                if j % 1000 == 0: print('current number of vessels {}'.format(len(self.list_of_vessels)))
                if vessel.in_growth:
                    total_path_length = vessel.step_size * micro_steps
                    new_path_length = vessel.grow(vegf_gradient, pressure, micro_steps, weight_direction, weight_vegf, weight_pressure)
                    if len(vessel.path) > 3:
                        splitting_rate_ = (new_path_length / total_path_length) * splitting_rate_
                        if random.uniform(0, 1) < splitting_rate_ * dt:
                            branching_point = vessel.choose_random_point()
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
