import random
from amber.BasicGeometries import *
import numpy as np
import sys

sys.setrecursionlimit(1500)

class Vessel: #class for vessels
    def __init__(self, path, radius, step_size, parent_id=None, children_ids=None, in_growth=True, intra_radiosensitivity=0.00001):
        if children_ids is None:
            children_ids = []
        for point in path:
            if not isinstance(point, np.ndarray) or len(point.shape) != 1 or point.shape[0] != 3: #check that the path is a list of 3D vectors
                raise ValueError("Each point in the 'path' list must be a 3D array with shape (3,)")

        self.path = np.array(path) #path is a list of 3D vectors giving the whole path of the vessel
        self.radius = radius #radius is the radius of the vessel
        self.id = id(self) #vessels id
        self.parent_id = parent_id #id of the vessel that gave rise to this one
        self.children_ids = children_ids #list of ids of vessels that arose from this one (max of 2)
        self.step_size = step_size #step size when growing the vessel
        self.in_growth = in_growth #boolean to indicate if the vessel is still growing
        self.visible = True #boolean to indicate if the vessel is visible when plotting
        self.must_be_updated = True #boolean to indicate if the vessel radius must be updated
        self.maturity = 1.0 # 1.0 is fully mature, changes vessels radius and increases slowly over time
        self.intra_radiosensitivity = intra_radiosensitivity #intra_radiosensitivity is the radiosensitivity of the vessel
    def __iter__(self): #iterator for the vessel
        return self

    def length(self): #returns the length of the vessel
        return np.sum(np.linalg.norm(np.diff(self.path, axis=0), axis=1))

    def grow(self, limit_half_length, lower_pressure_threshold_step, higher_pressure_threshold_step, step_stop_threshold, #grows the vessel by calling the step function enough times
             vegf_gradient, pressure, steps, weight_direction, weight_vegf, weight_pressure):
        length = 0
        for i in range(steps):
            step_size = self.step(limit_half_length, lower_pressure_threshold_step, higher_pressure_threshold_step, vegf_gradient, pressure, weight_direction, weight_vegf, weight_pressure)
            if step_size < step_stop_threshold: #if the last step was too small, then stop growing the vessel
                self.in_growth = False
                break
            length += step_size
        return length


    def step(self, half_length_world, lower_pressure_threshold_step, higher_pressure_threshold_step, #grows the vessel by one step
             vegf_gradient, pressure, weight_direction=0.5, weight_vegf=0.5, weight_pressure=0.5):

        if not isinstance(self.path, np.ndarray) or len(self.path.shape) != 2 or self.path.shape[1] != 3:
            raise ValueError("The 'self.path' array must be a 3D array with shape (n, 3)")
        last_point = np.array(self.path[-1]) #last point of the vessel path
        first_point = np.array(self.path[0]) #first point of the vessel path

        #check that the last point is within the world limits, if not, stop growing the vessel.
        if last_point[0] < -half_length_world or last_point[0] > half_length_world or \
                last_point[1] < -half_length_world or last_point[1] > half_length_world or \
                last_point[2] < -half_length_world or last_point[2] > half_length_world:
            return 0
        direction = last_point - first_point

        # Get VEGF gradient and pressure at the last point
        vegf_grad = np.array(vegf_gradient(last_point))
        local_pressure = pressure(last_point)

        # Normalize direction and VEGF gradient
        direction_norm = direction / np.linalg.norm(direction) if np.linalg.norm(direction) != 0 else direction
        vegf_grad_norm = vegf_grad / np.linalg.norm(vegf_grad) if np.linalg.norm(vegf_grad) != 0 else vegf_grad

        # Calculate the weighted direction based on vessel direction and VEGF gradient
        weighted_dir = weight_direction * direction_norm + weight_vegf * vegf_grad_norm

        pressure_vec = -weighted_dir #pressure vector is the opposite of the weighted direction. antagonist effect of pressure on vessel growth
        weighted_dir += weight_pressure * pressure_vec * local_pressure #add the pressure vector weighted by the local pressure

        # Add random noise
        noise = np.array([random.random()*2 - 1, random.random()*2-1, random.random()*2 - 1])
        weighted_dir += noise

        # Normalize the weighted direction
        weighted_dir /= np.linalg.norm(weighted_dir)

        # Calculate the new point and add it
        step_size = self.step_size

        if weight_pressure > 0 and local_pressure > lower_pressure_threshold_step: #if the local pressure is higher than the lower threshold, then decrease the step size
            step_size = step_size * (1 - (local_pressure - lower_pressure_threshold_step) / higher_pressure_threshold_step)

        new_point = last_point + step_size * weighted_dir #calculate the new point
        self.path = np.append(self.path, [new_point], axis=0)
        return step_size

    def volume_per_point(self): #returns the volume of the vessel associated to each point
        return np.pi * self.radius ** 2 * self.step_size

    def plot(self,fig, ax, color='crimson'): #plots the vessel
        if self.visible:
            if self.in_growth:
                ax.plot(self.path[:, 0], self.path[:, 1], self.path[:, 2], color='green', alpha=0.9, linewidth= self.radius*15)
            elif self.must_be_updated:
                ax.plot(self.path[:, 0], self.path[:, 1], self.path[:, 2], color='crimson', alpha=0.9, linewidth= self.radius*15)
            else :
                ax.plot(self.path[:, 0], self.path[:, 1], self.path[:, 2], color='blue', alpha=0.9, linewidth= self.radius*15)
        return fig, ax

    def choose_random_point(self, seed):
        # choose random point on the path, not the first or last point
        if len(self.path) < 3:
            print(self.path)
            raise ValueError("The vessel path is too short to choose a random point")
        random_index = random.randint(1, len(self.path) - 2)
        return self.path[random_index]

    def mean_pressure(self, pressure): #returns the mean pressure of the vessel
        if len(self.path) < 2:
            return 0
        else:
            return np.mean([pressure(point) for point in self.path])

    def radiosensitivity(self): #returns the radiosensitivity of the vessel
        radiosensitivity = self.intra_radiosensitivity
        return radiosensitivity
class VasculatureNetwork: #class that contains the list of vessels
    def __init__(self, config, list_of_vessels=None):
        if list_of_vessels is None:
            list_of_vessels = []
        self.list_of_vessels = list_of_vessels
        self.config = config

    def branching(self, vessel_id, branching_point): #branches a vessel at a given point
        if branching_point == [] or branching_point is None:
            print("No branching point found, branching point called was", branching_point)
            return
        mother_vessel = self.get_vessel(vessel_id) #get the vessel that will be branched
        if mother_vessel is None:
            raise ValueError("Vessel with id {} does not exist".format(vessel_id))

        path = mother_vessel.path #get the path of the mother vessel
        if len(path) < 3:
            raise ValueError("Vessel with id {} has a path with less than 3 points".format(vessel_id)) #the path must have at least 3 points
        visible = mother_vessel.visible
        split_index = np.where((path == branching_point))[0][0] #get the index of the branching point
        path_begin, path_end = np.split(path, [split_index]) #split the path in two parts
        #remove the first element of path_end to keep the same point distribution in space
        path_end = np.delete(path_end, 0, 0)
        mother_vessel.path = path_begin  # reduce the mother vessel path. It stops growing
        mother_vessel.in_growth = False

        # create two new vessels. One will be the previous end of the mother vessel and the other is the new vessel that will grow
        vessel_end = Vessel(path_end, mother_vessel.radius, self.config.vessel_step_size, parent_id= mother_vessel.id, children_ids=mother_vessel.children_ids, in_growth= mother_vessel.in_growth, intra_radiosensitivity= mother_vessel.intra_radiosensitivity)  # the radius has to be updated later
        for child_id in vessel_end.children_ids: #update the children of the vessel end
            child = self.get_vessel(child_id) #get the child
            child.parent_id = vessel_end.id #update the parent id of the child

        # the radius has to be updated later
        mother_vessel.visible = visible #update the visibility of the mother vessel
        vessel_end.visible = visible #update the visibility of the vessel end

        vessel_new = Vessel([branching_point], self.config.radius_root_vessels, self.config.vessel_step_size, parent_id= mother_vessel.id, children_ids=None, in_growth=True, intra_radiosensitivity=mother_vessel.intra_radiosensitivity)  # the radius has to be updated later
        vessel_new.visible = True #the new vessel is always visible
        vessel_new.maturity = self.config.maturity_new_vessels #the new vessel has maturity set up by user

        mother_vessel.children_ids = [vessel_end.id, vessel_new.id] #update the children of the mother vessel
        self.list_of_vessels.append(vessel_end) #add the vessel end to the list of vessels
        self.list_of_vessels.append(vessel_new) #add the new vessel to the list of vessels

        #find the root vessel:
        root_vessel = mother_vessel #start with the mother vessel
        while root_vessel.parent_id is not None: #while the vessel has a parent it is not the root vessel
            root_vessel = self.get_vessel(root_vessel.parent_id)

        root_vessel.must_be_updated = True #the root vessel tree needs to be updated as we added a vessel

    def add_vessel(self, vessel: Vessel): #adds a vessel to the list of vessels
        self.list_of_vessels.append(vessel)

    def add_multiple_vessels(self, vessels): #adds multiple vessels to the list of vessels
        for vessel in vessels:
            self.add_vessel(vessel)

    def get_vessel(self, vessel_id): #returns the vessel with the given id
        for vessel in self.list_of_vessels:
            if vessel.id == vessel_id:
                return vessel
        return None

    def kill_vessel(self, vessel_id): #kills a vessel and its children
        vessel = self.get_vessel(vessel_id)
        if vessel is None:
            raise ValueError("Vessel with id {} does not exist".format(vessel_id))
        # remove the children
        for child_id in list(vessel.children_ids):
            self.kill_vessel(child_id) #recursive call
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
        root_vessel.must_be_updated = True #the root vessel tree needs to be updated as we removed a vessel
        self.list_of_vessels.remove(vessel) #remove the vessel from the list of vessels

    def update_vessels_radius_from_last(self, final_radius, pressure_sensitive, pressure):
        print("Updating vessels radius from last") #update the radius of the vessels from the last to the first
        def update_radius_recursive(vessel_id): #recursive function to update the radius of the vessels
            vessel = self.get_vessel(vessel_id) #get the vessel
            if not vessel.children_ids: #if the vessel has no children, it is a terminal vessel and its radius is the final radius set up by the user
                vessel.radius = final_radius
            else: #if the vessel has children, its radius is given by Murray's law
                for child_id in vessel.children_ids:
                    update_radius_recursive(child_id)

                r_cubed_sum = sum([self.get_vessel(child_id).radius**3 for child_id in vessel.children_ids])
                vessel.radius = r_cubed_sum**(1/3)

        # Find all root vessels with parent_id=None
        root_vessels = [v for v in self.list_of_vessels if v.parent_id is None] #get the root vessels
        # root_vessels = [v for v in root_vessels if v.must_be_updated is True] #get the root vessels that need to be updated TODO: bugs with this, but it would improve efficiency
        print("Number of root vessels to update: ", len(root_vessels))

        # Call update_radius_recursive for each root vessel
        for root_vessel in root_vessels: #update the radius of each root vessel
            update_radius_recursive(root_vessel.id) #recursive call
            root_vessel.must_be_updated = True

        print("Updating for pressure") #update the radius of the vessels from the first to the last
        pow = self.config.radius_decrease_exponent #the exponent of the radius decrease
        coeff = self.config.max_occupancy ** (-pow) #the coefficient of the radius decrease
        for vessel in self.list_of_vessels: #update the radius of each vessel
            vessel.radius = vessel.radius * (1 - coeff * (vessel.mean_pressure(pressure)**pow)) #the radius is decreased by a factor that depends on the pressure
            vessel.radius = vessel.radius * vessel.maturity #the radius is decreased by a factor that depends on the maturity


    def volume_occupied(self): #returns the total volume occupied by the vessels
        points = []
        volume = []
        for vessel in self.list_of_vessels:
            for point in vessel.path:
                points.append(point)
                volume.append(vessel.volume_per_point())
        return points, volume

    def grow(self, vegf_gradient, pressure, steps, weight_direction, weight_vegf, weight_pressure): #grows all the vessels
        for vessel in self.list_of_vessels:
            if vessel.in_growth:
                vessel.grow(self.config.half_length_world, self.config.lower_pressure_threshold_step, self.config.higher_pressure_threshold_step, self.config.step_stop_threshold, vegf_gradient, pressure, steps, weight_direction, weight_vegf, weight_pressure)
    def plot(self, fig, ax, color='crimson'): #plots all the vessels
        for vessel in self.list_of_vessels:
            vessel.plot(fig, ax, color)
        return fig, ax


    def grow_and_split(self, dt, splitting_rate, vegf_gradient, pressure, macro_steps, micro_steps, weight_direction, weight_vegf, weight_pressure): #grows and splits all the vessels
        micro_steps = micro_steps #number of steps for each vessel growth cycle
        macro_steps_ = int(macro_steps * dt) #number of vessel growth cycles (has a probability of branching)

        if macro_steps_ == 0: #if user set 0 then deativate angiogenesis. otherwise if user just put a very small number, set it to 1
            ValueError("Macro steps must be at least 1")
            if macro_steps == 0:
                ValueError("Angiogenesis deactivation")
                macro_steps_ = 0
            else:
                macro_steps_ = 1

        for i in range(macro_steps_): #for each vessel growth cycle
            print("Macro step {}".format(i))
            j = 0
            for vessel in self.list_of_vessels: #for each vessel
                splitting_rate_ = splitting_rate
                if j % 1000 == 0: print('current number of vessels {}'.format(len(self.list_of_vessels)))
                if vessel.in_growth: #if the vessel is in growth
                    # total_path_length = vessel.step_size * micro_steps
                    new_path_length = vessel.grow(self.config.half_length_world, self.config.lower_pressure_threshold_step, self.config.higher_pressure_threshold_step, self.config.growth_step_stop_threshold, vegf_gradient, pressure, micro_steps, weight_direction, weight_vegf, weight_pressure)
                    if len(vessel.path) > 3: #if the vessel has at least 3 points
                        splitting_rate_ = vessel.length() * splitting_rate_ #the splitting rate is proportional to the length of the vessel
                        if random.uniform(0, 1) < splitting_rate_: #if the splitting rate is higher than a random number
                            branching_point = vessel.choose_random_point(self.config.seed) #choose a random point
                            self.branching(vessel.id, branching_point) #branch the vessel
                j += 1
    def print_vessel_tree_recursive(self, vessels, children_ids, indent): #used to print the tree of vessels for debugging
        for child_id in children_ids:
            child_vessel = next((v for v in vessels if v.id == child_id), None)
            if child_vessel is not None:
                print(' ' * indent, f"ID: {child_vessel.id}  Radius: {child_vessel.radius:.5f}")
                if child_vessel.children_ids:
                    self.print_vessel_tree_recursive(vessels, child_vessel.children_ids, indent + 2)


    def print_vessel_tree(self, indent=0): #prints the tree of vessels for debugging
        vessels = self.list_of_vessels
        # get the root vessels (i.e. vessels with no parent)
        root_vessels = [v for v in vessels if v.parent_id is None]

        # recursively print the vessel tree
        for root_vessel in root_vessels:
            print(' ' * indent, f"ID: {root_vessel.id}  Radius: {root_vessel.radius:.5f}")
            if root_vessel.children_ids:
                self.print_vessel_tree_recursive(vessels, root_vessel.children_ids, indent + 2)

    def compute_VSL(self): #computes all the Vascular Segment Length
        list_VSL = np.array([])
        for vessel in self.list_of_vessels:
            if vessel.visible:
                list_VSL = np.append(list_VSL, vessel.length())
        return list_VSL

    def compute_diameters(self): #computes all the diameters of the vessels
        list_diameters = np.array([])
        for vessel in self.list_of_vessels:
            if vessel.visible:
                list_diameters = np.append(list_diameters, vessel.radius*2)
        return list_diameters
