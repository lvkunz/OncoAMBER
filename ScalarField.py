import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt



class ScalarField:
    def __init__(self, points, values, gradient_step_size=1e-5):
        self.points = np.array(points)
        self.values = np.array(values)
        x = []
        y = []
        z = []
        for i in range(len(self.points)):
            x.append(self.points[i][0])
            y.append(self.points[i][1])
            z.append(self.points[i][2])
        x = np.unique(x)
        y = np.unique(y)
        z = np.unique(z)
        self.values = self.values.reshape(len(x), len(y), len(z))
        self.interpolator = RegularGridInterpolator((x,y,z), self.values, method='linear', bounds_error=False, fill_value=None)
        self.gradient_step_size = gradient_step_size

    def evaluate(self, point):
        return self.interpolator((point[0], point[1], point[2]))

    def gradient(self, point):
        return np.array([(self.evaluate(point + d) - self.evaluate(point - d)) / (2 * self.gradient_step_size) for d in
                         np.identity(len(point)) * self.gradient_step_size]).T

    def show_values(self, fig, ax, cmap='jet', vmin=None, vmax=None):
        if vmin is None:
            vmin = self.values.min()
        if vmax is None:
            vmax = self.values.max()
        ax.scatter(self.points[:, 0], self.points[:, 1], self.points[:, 2], c=self.values, cmap=cmap, s=15, vmin=vmin, vmax=vmax, alpha=0.5)
        return fig, ax

    def show_gradient(self, fig, ax, length_scale=1, color='black'):
        for point in self.points:
            print(point)
            ax.quiver(*point, *self.gradient(point), color='black', length=length_scale, normalize=True)
        return fig, ax

    def show_all(self, fig, ax, length_scale=1, cmap='jet', vmin=None, vmax=None):
        fig, ax = self.show_values(fig, ax, cmap=cmap, vmin=vmin, vmax=vmax)
        fig, ax = self.show_gradient(fig, ax, length_scale=length_scale)
        return fig, ax
