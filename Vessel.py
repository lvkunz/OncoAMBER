import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import math
from BasicPlots import *
from BasicGeometries import *

def rotation_matrix_from_vectors(u, v):
    # Find the rotation matrix that aligns vec1 to vec2
    # https://stackoverflow.com/questions/43507479/how-to-calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    c = np.dot(u, v)
    s = np.sqrt(1 - c ** 2)
    k = np.cross(u, v)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.identity(3) + s * K + (1 - c) * np.dot(K, K)
    return R
class Vessel:
    def __init__(self, origin, end, radius, id = 0):


        self.origin = np.array(origin)
        self.end = np.array(end)
        self.radius = radius
        self.direction_vector = self.end - self.origin
        self.length = np.linalg.norm(self.direction_vector)
        self.direction_vector = self.direction_vector / self.length
        self.rotation_matrix = np.array(rotation_matrix_from_vectors( u = [0, 0, 1], v = self.direction_vector))
        self.id = id
    def __iter__(self):
        return self

    def print(self):
        print('Vessel ::' + 'origin: '+ str(self.origin) + ' end: '+ str(self.end) + ' radius: '+  str(self.radius))

    def generate_random_points(self, n):
        # Generate random points in a bounding box that encloses the cylinder
        x = np.random.uniform(-self.radius, self.radius, n)
        y = np.random.uniform(-self.radius,self.radius, n)
        z = np.random.uniform(0, self.length, n)
        points = np.array([x, y, z]).T
        # Rotate the points to align with the cylinder
        points = np.dot(points, self.rotation_matrix)
        print(points)
        # Translate the points to the correct position
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

class VasculatureNetwork:
    def __init__(self, bounds : Shape = Sphere(center = np.array([0,0,0]), radius = 3.0), list_of_vessels = None):
        self.bounds = bounds
        self.next_vessel_number = 0
        if list_of_vessels is None:
            self.list_of_vessels = self.generate_vasculature(20)
            self.next_vessel_number = len(self.list_of_vessels)
        else:
            self.list_of_vessels = list_of_vessels
            self.next_vessel_number = len(self.list_of_vessels)

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
    def generate_vasculature(self, n, bounds : Shape = Sphere(center = np.array([0,0,0]), radius = 3.0)):
        self.bounds = bounds
        points = self.bounds.generate_random_points(n)
        starting = points[0]
        ending = points[1]
        self.list_of_vessels = [Vessel(starting, ending, 0.1, id = 0)]
        for i in range(2, n):
            if i % 100 == 0: print('Generating Vasculature, current number of vessels: ', i)
            end = points[i]
            self.list_of_vessels.append(Vessel(self.closest_point(end), end, 0.1, id = i-1))
        self.next_vessel_number = len(self.list_of_vessels)
        return self.list_of_vessels
    def add_vessel(self, origin, end, radius):
        self.list_of_vessels.append(Vessel(origin, end, radius, id = self.next_vessel_number))
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
        print('Plotting Vasculature')
        for vessel in self.list_of_vessels:
            x = [vessel.origin[0], vessel.end[0]]
            y = [vessel.origin[1], vessel.end[1]]
            z = [vessel.origin[2], vessel.end[2]]
            ax.plot(x, y, z, color='red', linewidth=1, alpha=0.5)

    def save(self, filename):
        with open(filename, 'w') as f:
            for vessel in self.list_of_vessels:
                f.write(str(vessel.origin[0]) + ' ' + str(vessel.origin[1]) + ' ' + str(vessel.origin[2]) + ' ' \
                        + str(vessel.end[0]) + ' ' + str(vessel.end[1]) + ' ' + str(vessel.end[2]) + ' ' + str(vessel.radius) + '\n')
    def read(self, filename):
        self.list_of_vessels = []
        with open(filename, 'r') as f:
            for line in f:
                line = line.split()
                self.list_of_vessels.append(Vessel([float(line[0]), float(line[1]), float(line[2])], [float(line[3]), float(line[4]), float(line[5])], float(line[6])))

