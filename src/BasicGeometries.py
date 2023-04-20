import numpy as np

class Shape(object):
    def __init__(self, name: str = "Shape"):
        self.name = name
    def generate_random_points(self, n):
        raise NotImplementedError("generate_random_points not implemented")

class Sphere(Shape):
    def __init__(self, radius: float = 1.0, center = np.array([0,0,0]) , name: str = "Sphere"):
        super().__init__(name)
        self.radius = radius
        self.center = center
        self.volume = 4/3 * np.pi * self.radius**3
    def generate_random_points(self, n):
        # generate n random points within the tumor
        # the points are generated in a cube of side 2*radius
        # then the points outside the tumor are discarded
        # initialize the points array as a list of 3D points
        points = np.array([[0, 0, 0]])
        while points.shape[0] <= n+2:
            points_new = np.random.uniform(-self.radius, self.radius, (n, 3))
            points_new = points_new + self.center
            points_new = points_new[np.linalg.norm(points_new - self.center, axis=1) < self.radius]
            points = np.concatenate((points, points_new))
        # make sure the number of points is n
        points = points[1:]
        points = points[:n]
        return points

    def generate_random_points_on_surface(sphere, n):
        points = []
        for i in range(n):
            # generate a random point on a unit sphere
            point = np.random.normal(size=3)
            point /= np.linalg.norm(point)

            # scale the point by the sphere's radius and translate it to the sphere's center
            point *= sphere.radius
            point += sphere.center

            points.append(point)

        return points

class Cube(Shape):
    def __init__(self, half_side: float = 1.0, center = np.array([0,0,0]) , name: str = "Cube"):
        super().__init__(name)
        self.half_side = half_side
        self.side = 2*half_side
        self.center = center
        self.volume = self.side**3
    def generate_random_points(self, n):
        # generate n random points within the cube
        # the points are generated in a cube of side side
        points = np.random.uniform(-self.half_side, self.half_side, (n, 3))
        return points