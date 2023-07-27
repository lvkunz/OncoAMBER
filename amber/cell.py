import numpy as np
import scipy.stats as stats
from scipy.stats import gamma
class Cell:
    def __init__(self, radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold, VEGF_rate, type = 'NormalCell'):
        self.radius = radius
        self.necrotic = False
        self.usual_cycle_length = cycle_hours
        self.cycle_length_std = cycle_std
        self.gamma_shape = (self.usual_cycle_length / self.cycle_length_std) ** 2
        self.gamma_scale = self.cycle_length_std ** 2 / self.usual_cycle_length
        self.doubling_time = self.random_doubling_time() #new cells have a random doubling time
        self.volume = 4/3 * np.pi * self.radius**3
        self.oxygen = np.random.uniform(0, 1)
        self.intra_radiosensitivity = intra_radiosensitivity
        self.o2_to_vitality_factor = o2_to_vitality_factor
        self.VEGF_threshold = VEGF_threshold
        self.VEGF_rate = VEGF_rate
        if type != 'NormalCell' and type != 'TumorCell':
            print(type + ' is not a valid cell type')
            raise ValueError('Cell type must be either NormalCell or TumorCell')
        self.type = type
        self.time_before_death = None
        self.time_spent_cycling = self.random_time_spent_cycling() #new cells have already been cycling for a random amount of time
        self.pH = 7.4
        self.damage = 0
    def duplicate(self):
        cell_class = type(self)  # Get the class of the current instance dynamically
        cell = cell_class(self.radius,
                          self.usual_cycle_length,
                          self.cycle_length_std,
                          self.intra_radiosensitivity,
                          self.o2_to_vitality_factor,
                          self.VEGF_threshold,
                          self.VEGF_rate,
                          )
        cell.time_spent_cycling = 0  # The new cell has not been cycling yet
        return cell

    def metabolic_rate(self):
        if self.damage < 0.5:
            return 1.0
        elif self.damage < 0.8:
            return 0.5
        else:
            return 0.0

    def vitality(self):
        factor = self.o2_to_vitality_factor
        vitality = self.oxygen * factor
        vitality = min(vitality, 1.0)
        return vitality #needs to be normalized between 0 and 1

    def radiosensitivity(self):
        coeff_pH = 7.4 / self.pH
        value = self.intra_radiosensitivity * (0.5 + 0.5 * self.oxygen) * coeff_pH
        # if self.necrotic:
        #     value = self.intra_radiosensitivity * 0.1
        return value

    def necrosis_probability(self, proba0, threshold, coeff_damage):
        p = 0
        r = proba0 - self.vitality() / threshold
        if r < 0: r = 0
        p += r
        if self.damage > 0:
            p += coeff_damage * self.damage
        return p

    def apoptosis_probability(self, proba0, threshold, coeff_damage):
        p = 0
        if self.vitality() < threshold:
            p += proba0
        if self.damage > 0:
            p += coeff_damage * self.damage
        return p

    def damage_repair(self, dt, repair_per_hour):
        def repair_amount(damage):
            if damage == 0:
                return 0
            else:
                return max(repair_per_hour * (1 - damage), 0.0001)

        for i in range(int(dt)):
            self.damage -= repair_amount(self.damage)
        if self.damage < 0:
            self.damage = 0

    def random_doubling_time(self):
        return np.random.gamma(shape = self.gamma_shape, scale = self.gamma_scale)

    def random_time_spent_cycling(self):
        longest_possible_time = self.doubling_time
        return np.random.uniform(0, longest_possible_time)

    def VEGF_secretion(self):
        if self.vitality() < self.VEGF_threshold:
            return self.VEGF_rate * (1 - self.vitality()) * self.metabolic_rate()
        else:
            return 0

    def is_cycling(self, vitality_threshold, damage_threshold):
        vit = self.vitality() > vitality_threshold
        dmg = self.damage > damage_threshold
        return vit and not dmg


    def probability_of_interaction(self, cell, dt):
        return 0.0

    def interact(self, cell, dt):
        pass

    def fiber_secretion(self, dt):
        return 0.0

class TumorCell(Cell):
    def __init__(self, radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold, VEGF_rate):
        super().__init__(radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold, VEGF_rate,  type =  'TumorCell')

class ImmuneCell(Cell):
    def __init__(self, radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor):
        super().__init__(radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold=0.0, VEGF_rate=0.0, type = 'ImmuneCell')
    #immune cells don't cycle, so they don't need a doubling time or time spent cycling
    def probability_of_interaction(self, cell, dt):
        return 0.001 * dt

    def interact(self, cell, dt):
        if cell.type == 'TumorCell':
            cell.time_before_death = 0

class Fibroblast(Cell):

    def __init__(self, radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, fiber_secretion_rate):
        super().__init__(radius, cycle_hours, cycle_std, intra_radiosensitivity, o2_to_vitality_factor, VEGF_threshold = 0.0, VEGF_rate= 0.0, type = 'Fibroblast')
        self.fiber_secretion_rate = fiber_secretion_rate
    def fiber_secretion(self, dt):
        return self.fiber_secretion_rate * dt

