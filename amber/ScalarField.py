import numpy as np
from scipy.interpolate import RegularGridInterpolator


class ScalarField3D: #this class is used to interpolate a scalar field in 3D from a set of points. Used for VEGF, pressure, etc.
    def __init__(self, points, values, gradient_step_size=1e-5, bounds_error=False, fill_value=None):
        self.points = np.array(points)
        self.values = np.array(values)

        x, y, z = np.unique(self.points[:, 0]), np.unique(self.points[:, 1]), np.unique(self.points[:, 2])

        self.values = self.values.reshape(len(x), len(y), len(z))
        self.interpolator = RegularGridInterpolator((x, y, z), self.values, method='linear', bounds_error=bounds_error,
                                                    fill_value=fill_value)
        self.gradient_step_size = gradient_step_size

    def evaluate(self, point): #return value at any point
        return self.interpolator((point[0], point[1], point[2]))

    def gradient(self, point): #return gradient at any point
        return np.array([(self.evaluate(point + d) - self.evaluate(point - d)) / (2 * self.gradient_step_size) for d in
                         np.identity(len(point)) * self.gradient_step_size]).T

    def show_values(self, fig, ax, cmap='jet', vmin=None, vmax=None): #show values as a scatter plot
        if vmin is None:
            vmin = self.values.min()
        if vmax is None:
            vmax = self.values.max()
        ax.scatter(self.points[:, 0], self.points[:, 1], self.points[:, 2], c=self.values, cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.5)
        return fig, ax

    def show_gradient(self, fig, ax, length_scale=1, color='black'): #show gradient as a quiver plot
        for point in self.points:
            print(point)
            ax.quiver(*point, *self.gradient(point), color='black', length=length_scale, normalize=True)
        return fig, ax

    def show_all(self, fig, ax, length_scale=1, cmap='jet', vmin=None, vmax=None): #show both values and gradient
        fig, ax = self.show_values(fig, ax, cmap=cmap, vmin=vmin, vmax=vmax)
        fig, ax = self.show_gradient(fig, ax, length_scale=length_scale)
        return fig, ax

class ScalarField2D: #this class is used to interpolate a scalar field in 2D from a set of points. Used for alpha, beta maps
    def __init__(self, points, values, gradient_step_size=1e-5, bounds_error=False, fill_value=None):
        self.points = np.array(points)
        self.values = np.array(values)

        x, y = np.unique(self.points[:, 0]), np.unique(self.points[:, 1])

        self.values = self.values.reshape(len(x), len(y))
        self.interpolator = RegularGridInterpolator((x, y), self.values, method='linear', bounds_error=bounds_error,fill_value=fill_value)
        self.gradient_step_size = gradient_step_size
        self.bounds_error = bounds_error
    def evaluate(self, point):
        if self.bounds_error:
            if point[0] < min(self.points[:, 0]) or point[0] > max(self.points[:, 0]) or point[1] < min(self.points[:, 1]) or point[1] > max(self.points[:, 1]):
                print('Point', point, 'is out of bounds')

        return self.interpolator((point[0], point[1]))

    def show(self, fig, ax, cmap='jet', vmin=None, vmax=None):
        ax.scatter(self.points[:, 0], self.points[:, 1], c=self.values, cmap=cmap, s=15, vmin=vmin, vmax=vmax, alpha=0.5)
        return fig, ax

    def show_extra(self, fig, ax, range_x, range_y): #show values as a surface plot
        x, y = np.meshgrid(np.linspace(range_x[0], range_x[1], 100), np.linspace(range_y[0], range_y[1], 100))
        z = np.zeros_like(x)
        for i in range(len(x)):
            for j in range(len(x[0])):
                z[i, j] = self.evaluate((x[i, j], y[i, j]))
        ax.plot_surface(x, y, z, cmap='gray', alpha=0.9, vmin=0, vmax=15)
        return fig, ax
