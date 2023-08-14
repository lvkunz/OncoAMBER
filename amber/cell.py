import numpy as np
import scipy.stats as stats
from scipy.stats import gamma

#file for defining the cell class and subclasses
class Cell:
    def __init__(self, radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold, VEGF_rate, type = 'NormalCell'):
        self.radius = radius #radius in mm of the cell
        self.necrotic = False #whether the cell is necrotic or not
        self.usual_cycle_length = cycle_hours #average time in hours for the cell to divide
        self.cycle_length_std = cycle_std #standard deviation of the cycle length
        self.gamma_shape = (self.usual_cycle_length / self.cycle_length_std) ** 2 #shape parameter for the gamma distribution
        self.gamma_scale = self.cycle_length_std ** 2 / self.usual_cycle_length #scale parameter for the gamma distribution
        self.doubling_time = self.random_doubling_time() #new cells have a random doubling time
        self.volume = 4/3 * np.pi * self.radius**3 #volume of the cell in mm^3
        self.oxygen = np.random.uniform(0, 1) #random oxygen level between 0 and 1 when cell is created
        self.intra_radiosensitivity = intra_radiosensitivity  # radiosensitivity of the cell
        self.o2_to_vitality_factor = o2_to_vitality_factor #factor for how much oxygen affects the vitality of the cell
        self.VEGF_threshold = VEGF_threshold #threshold for vitality of the cell to secrete VEGF
        self.VEGF_rate = VEGF_rate #rate at which the cell secretes VEGF
        self.type = type #type of cell
        self.time_before_death = None #time in hours before the cell dies
        self.time_spent_cycling = self.random_time_spent_cycling() #new cells have already been cycling for a random amount of time
        self.pH = 7.4 #pH of the cell
        self.damage = 0 #amount of damage of the cell
    def duplicate(self): #returns a copy of the cell
        cell_class = type(self)  # Get the class of the current instance dynamically
        cell = cell_class(radius = self.radius,
                          cycle_hours = self.usual_cycle_length,
                          cycle_std = self.cycle_length_std,
                          intra_radiosensitivity = self.intra_radiosensitivity,
                          o2_to_vitality_factor = self.o2_to_vitality_factor,
                          VEGF_threshold = self.VEGF_threshold,
                          VEGF_rate= self.VEGF_rate,
                          )
        cell.time_spent_cycling = 0  # The new cell has not been cycling yet
        return cell

    def metabolic_rate(self, threshold): #returns the metabolic rate of the cell. This is a function of the damage
        if self.damage <= threshold:
            return 1.0 - self.damage
        else:
            return 0.0

    def vitality(self): #returns the vitality of the cell. This is a function of the oxygen
        factor = self.o2_to_vitality_factor
        vitality = self.oxygen * factor
        vitality = min(vitality, 1.0)
        return vitality #needs to be normalized between 0 and 1

    def radiosensitivity(self): #returns the radiosensitivity of the cell. This is a function of the pH and the oxygen
        coeff_pH = 7.4 / self.pH
        value = self.intra_radiosensitivity * (0.5 + 0.5 * self.oxygen) * coeff_pH

        if self.necrotic == True: #if the cell is necrotic, it is not very radiosensitive. Might have to be modified
            value = 0.1
        return value

    def necrosis_probability(self, proba0, threshold, coeff_damage): #returns the probability of necrosis of the cell. This is a function of the vitality and the damage
        p = 0

        r = proba0
        if self.vitality() < threshold:
            p += r

        if self.damage > 0:
            p += coeff_damage * self.damage
        p = min(p, 1.0)
        return p

    def apoptosis_probability(self, proba0, threshold, coeff_damage): #returns the probability of apoptosis of the cell. This is a function of the vitality and the damage
        p = 0
        if self.vitality() < threshold:
            p += proba0
        if self.damage > 0:
            p += coeff_damage * self.damage
        p = min(p, 1.0)
        return p

    def damage_repair(self, dt, repair_per_hour): #repairs the damage of the cell. This is a function of the damage and the repair rate
        def repair_amount(damage):
            if damage == 0:
                return 0
            else:
                return max(repair_per_hour * (1 - damage), 0.0001)

        if dt < 1:
            raise ValueError("dt must be greater than 1")
        
        for i in range(int(dt)):
            self.damage -= repair_amount(self.damage)

        if self.damage < 0:
            self.damage = 0

    def random_doubling_time(self): #returns a random time of division for the cell
        return np.random.gamma(shape = self.gamma_shape, scale = self.gamma_scale)

    def random_time_spent_cycling(self): #returns a random time spent cycling for the cell
        longest_possible_time = self.doubling_time
        return np.random.uniform(0, longest_possible_time)

    def VEGF_secretion(self, damage_threshold): #returns the VEGF secretion rate of the cell. This is a function of the vitality and the metabolic rate
        if self.vitality() < self.VEGF_threshold:
            return self.VEGF_rate * (1 - self.vitality()) * self.metabolic_rate(damage_threshold)
        else:
            return 0

    def is_cycling(self, vitality_threshold, damage_threshold): #returns True if the cell is cycling, False otherwise. This is a function of the vitality and the damage
        vit = self.vitality() > vitality_threshold
        dmg = self.damage > damage_threshold
        return vit and not dmg


    def probability_of_interaction(self, cell, dt): #returns the probability of interaction with another cell. Placeholder function for now
        return 0.0

    def interact(self, cell, dt): #interacts with another cell. Placeholder function for now
        pass

    def fiber_secretion(self, dt): #returns the fiber secretion rate of the cell. Placeholder function for now
        return 0.0


#subclass for different types of cells. Can be used to add new cell types and methods for each cell
class TumorCell(Cell):
    def __init__(self, radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold, VEGF_rate):
        super().__init__(radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold, VEGF_rate,  type =  'TumorCell')

class ImmuneCell(Cell):
    def __init__(self, radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor):
        super().__init__(radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold=0.0, VEGF_rate=0.0, type = 'ImmuneCell')
    #immune cells don't cycle, so they don't need a doubling time or time spent cycling
    def probability_of_interaction(self, cell, dt): #returns the probability of interaction with another cell.
        return 0.001 * dt

    def interact(self, cell, dt): #interacts with another cell. If the other cell is a tumor cell, it kills it
        if cell.type == 'TumorCell':
            cell.time_before_death = 0

class Fibroblast(Cell):

    def __init__(self, radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, fiber_secretion_rate):
        super().__init__(radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold = 0.0, VEGF_rate= 0.0, type = 'Fibroblast')
        self.fiber_secretion_rate = fiber_secretion_rate
    def fiber_secretion(self, dt): #returns the fiber secretion rate of the cell.
        return self.fiber_secretion_rate * dt

