import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

class ScalarField:
    def __init__(self, values, x_coords, y_coords, z_coords):
        self.values = values
        self.x_coords = x_coords
        self.y_coords = y_coords
        self.z_coords = z_coords
        self.interpolator = RegularGridInterpolator((x_coords, y_coords, z_coords), values, bounds_error=False, fill_value=0)

    def __call__(self, x, y, z):
        return self.interpolator((x, y, z))

    def plot(self, z_stack=0, range = [0,1]):
        # plot the scalar field
        plt.figure()
        plt.title('Scalar field, z = ' + str(self.z_coords[z_stack]) + ' mm')
        plt.imshow(self.values[:,:,z_stack], vmin=range[0], vmax=range[1], extent=[self.x_coords[0], self.x_coords[-1], self.y_coords[0], self.y_coords[-1]])
        plt.colorbar()
        plt.show()


class Kernel(ScalarField):
    def __init__(self, values, x_coords, y_coords, z_coords):
        super().__init__(values, x_coords, y_coords, z_coords)
        self.values = values
        self.x_coords = x_coords
        self.y_coords = y_coords
        self.z_coords = z_coords
        self.interpolator = RegularGridInterpolator((x_coords, y_coords, z_coords), values, bounds_error=False, fill_value=0)

    def __call__(self, x, y, z, center_kernel = np.array([0,0,0]), scale = 1 ): #pas sur de cette methode
        x = x - center_kernel[0]
        y = y - center_kernel[1]
        z = z - center_kernel[2]
        x = x * scale
        y = y * scale
        z = z * scale
        return self.interpolator((x, y, z))

    def plot(self, z_stack=0, range = [0,1]):
        # plot the scalar field
        plt.figure()
        plt.title('Kernel, z = ' + str(self.z_coords[z_stack]) + ' mm')
        plt.imshow(self.values[:,:,z_stack], vmin=range[0], vmax=range[1], extent=[self.x_coords[0], self.x_coords[-1], self.y_coords[0], self.y_coords[-1]])
        plt.colorbar()
        plt.show()
