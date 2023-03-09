import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import math
from BasicPlots import *
class Vessel:
    def __init__(self, origin, end, radius, mesh_density=300):


        self.origin = np.array(origin)
        self.end = np.array(end)
        self.radius = radius
        self.direction_vector = self.end - self.origin
        self.length = np.linalg.norm(self.direction_vector)
        self.direction_vector = self.direction_vector / self.length
        self.rotation_matrix = np.array(rotation_matrix_from_vectors( u = [0, 0, 1], v = self.direction_vector))

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