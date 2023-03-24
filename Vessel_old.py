import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import math
from BasicPlots import *
from BasicGeometries import *

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    if any(v):  # if not all zeros then
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    else:
        return np.eye(3)  # cross of all zeros only occurs on identical directions


class Vessel_old:
    def __init__(self, origin, end, radius, id=0):
        self.origin = np.array(origin)
        self.end = np.array(end)
        self.radius = radius
        self.direction_vector = self.end - self.origin
        self.length = np.linalg.norm(self.direction_vector)
        self.direction_vector = self.direction_vector / self.length
        self.rotation_matrix = np.array(rotation_matrix_from_vectors(self.direction_vector, [0.0, 0.0, 1.0]))
        self.id = id
        self.center = (self.origin + self.end) / 2  # center of the vessel

    def __iter__(self):
        return self

    def print(self):
        print('Vessel ::' + 'origin: ' + str(self.origin) + ' end: ' + str(self.end) + ' radius: ' + str(self.radius))

    def generate_random_points(self, n):
        # Generate random points in a bounding box that encloses the cylinder
        x = np.random.uniform(-self.radius, self.radius, n)
        y = np.random.uniform(-self.radius, self.radius, n)
        z = np.random.uniform(0, self.length, n)
        points = np.array([x, y, z]).T
        # Rotate the points to align with the cylinder
        points = np.dot(points, self.rotation_matrix)
        print(points)
        # Translate the points to the correct position
        return points

    def generate_points_along_axis(self, n):
        # generate random points along the axis of the vessel
        z = np.linspace(0, self.length, n)
        x = np.zeros(n)
        y = np.zeros(n)
        points = np.array([x, y, z]).T
        # Rotate the points to align with the cylinder
        points = np.dot(points, self.rotation_matrix)
        # Translate the points to the correct position
        points = points + self.origin
        return points

    def closest_point(self, p):
        """
        Get the point on the vessel that is closest to point p.
        """
        v = p - self.origin
        projection = np.dot(v, self.direction_vector)
        projection = np.clip(projection, 0, self.length)
        closest_point = self.origin + projection * self.direction_vector
        return closest_point

    def distance(self, p):
        """
        Get the distance from point p to the vessel.
        """
        closest_point = self.closest_point(p)
        return np.linalg.norm(p - closest_point)

    def plot(self, fig, ax):
        x = [self.origin[0], self.end[0]]
        y = [self.origin[1], self.end[1]]
        z = [self.origin[2], self.end[2]]
        ax.plot(x, y, z, color='orangered', linewidth=1, alpha=0.3)
        return fig, ax


class VasculatureNetwork:
    def __init__(self, list_of_vessels=None):
        if list_of_vessels is None:
            list_of_vessels = []
        self.list_of_vessels = list_of_vessels
        self.next_vessel_number = len(list_of_vessels)

    def closest_point(self, p):
        """
        Get the point on the vasculature that is closest to point p.
        """
        closest_point = None
        closest_distance = None
        for vessel in self.list_of_vessels:
            current_point = vessel.closest_point(p)
            current_distance = vessel.distance(p)
            if closest_point is None or current_distance < closest_distance:
                closest_point = current_point
                closest_distance = current_distance
        return closest_point

    def closest_distance(self, p):
        closest_distance = None
        for vessel in self.list_of_vessels:
            current_distance = vessel.distance(p)
            if closest_distance is None or current_distance < closest_distance:
                closest_distance = current_distance
        return closest_distance

    def grow_vasculature(self, random_points=None):
        if random_points is None:
            random_points = []
        if len(self.list_of_vessels) == 0:
            raise ValueError('You need at least one initial vessel to grow a vasculature network')
        for i in range(0, len(random_points)):
            if i % 100 == 0: print('Growing Vasculature, current number of vessels added: ', i)
            end = random_points[i]
            self.list_of_vessels.append(Vessel_old(self.closest_point(end), end, 0.1, id=self.next_vessel_number))
            self.next_vessel_number = self.next_vessel_number + 1
        return self.list_of_vessels

    def add_vessel(self, origin, end, radius):
        self.list_of_vessels.append(Vessel_old(origin, end, radius, id=self.next_vessel_number))
        self.next_vessel_number += 1

    def remove_vessel(self, id):
        for vessel in self.list_of_vessels:
            if vessel.id == id:
                self.list_of_vessels.remove(vessel)
                break

    def find_vessel(self, id):
        for vessel in self.list_of_vessels:
            if vessel.id == id:
                return vessel
        return None

    def plot(self, fig, ax):
        print('-- Plotting Vasculature')
        for vessel in self.list_of_vessels:
            x = [vessel.origin[0], vessel.end[0]]
            y = [vessel.origin[1], vessel.end[1]]
            z = [vessel.origin[2], vessel.end[2]]
            ax.plot(x, y, z, color='crimson', linewidth=1.5, alpha=0.7)

    def save(self, filename):
        with open(filename, 'w') as f:
            for vessel in self.list_of_vessels:
                f.write(str(vessel.origin[0]) + ' ' + str(vessel.origin[1]) + ' ' + str(vessel.origin[2]) + ' ' \
                        + str(vessel.end[0]) + ' ' + str(vessel.end[1]) + ' ' + str(vessel.end[2]) + ' ' + str(
                    vessel.radius) + '\n')
        print('-- Vasculature saved to file: ', filename)

    def read(self, filename):
        print('-- Reading vasculature from file: ', filename)
        self.list_of_vessels = []
        self.next_vessel_number = 0
        with open(filename, 'r') as f:
            for line in f:
                line = line.split()
                self.list_of_vessels.append(Vessel_old([float(line[0]), float(line[1]), float(line[2])],
                                                       [float(line[3]), float(line[4]), float(line[5])], float(line[6]),
                                                       id=self.next_vessel_number))
                self.next_vessel_number += 1
        print('-- Vasculature read from file: ', filename)
