import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf


class ScalarField:
    def __init__(self, points, values, gradient_step_size=1e-5):
        self.points = np.array(points)
        self.values = np.array(values)
        self.rbf = Rbf(*self.points.T, self.values)
        self.gradient_step_size = gradient_step_size

    def evaluate(self, point):
        return self.rbf(*point)

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
            length = np.linalg.norm(self.gradient(point))
            ax.quiver(*point, *self.gradient(point), color='black', length=length_scale * length, normalize=True)
        return fig, ax

    def show_all(self, fig, ax, length_scale=1, cmap='jet', vmin=None, vmax=None):
        fig, ax = self.show_values(fig, ax, cmap=cmap, vmin=vmin, vmax=vmax)
        fig, ax = self.show_gradient(fig, ax, length_scale=length_scale)
        return fig, ax
