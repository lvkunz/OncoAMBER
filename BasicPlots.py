import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import ReadAndWrite as rw
from matplotlib.colors import LinearSegmentedColormap


CONFIG = rw.read_config_file('CONFIG.txt')
seed = CONFIG['seed']
if seed == -1:
    seed = np.random.randint(0, 1000000)
np.random.seed(seed)
print('seed: ', seed)

def plot_sphere(ax, fig, radius, center, color='blue'):
    # show the cell
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
    y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
    ax.plot_surface(x, y, z, color=color, alpha=0.5)
    return fig, ax


def plot_cube(ax, fig, center, half_length, color='black'):
    # Create a figure and add a 3D axes
    # Get the coordinates of the corners of the cube
    x = center[0] + half_length*np.array([-1, 1, 1, -1, -1, 1, 1, -1])
    y = center[1] + half_length*np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    z = center[2] + half_length*np.array([-1, -1, -1, -1, 1, 1, 1, 1])

    # Plot the cube
    lw = 0.1
    ax.plot([x[0], x[4]], [y[0], y[4]], [z[0], z[4]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[1], x[5]], [y[1], y[5]], [z[1], z[5]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[2], x[6]], [y[2], y[6]], [z[2], z[6]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[3], x[7]], [y[3], y[7]], [z[3], z[7]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[0], x[1]], [y[0], y[1]], [z[0], z[1]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[1], x[2]], [y[1], y[2]], [z[1], z[2]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[2], x[3]], [y[2], y[3]], [z[2], z[3]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[3], x[0]], [y[3], y[0]], [z[3], z[0]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[4], x[5]], [y[4], y[5]], [z[4], z[5]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[5], x[6]], [y[5], y[6]], [z[5], z[6]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[6], x[7]], [y[6], y[7]], [z[6], z[7]], color, linewidth=lw, alpha=0.5)
    ax.plot([x[7], x[4]], [y[7], y[4]], [z[7], z[4]], color, linewidth=lw, alpha=0.5)
    return fig, ax

def to_rgb(color_name: str):
    color_map = {
        "black": (0, 0, 0),
        "red": (255, 0, 0),
        "green": (0,255,0),
        "blue": (0, 0, 255),
        "yellow": (255, 255, 0),
        "cyan": (0, 255, 255),
        "magenta": (255, 0, 255),
        "my green": (19, 136, 8),
        "my purple": (191, 64, 191)

        # Add more colors and their corresponding RGB values as needed
    }
    if color_name in color_map:
        return color_map[color_name]
    else:
        return color_map["black"]

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

def create_custom_RdBu(separation):
    # Define the RdBu colormap
    RdBu = plt.cm.get_cmap('RdBu')

    # Create the custom colormap segments
    cdict = {
        'red': [],
        'green': [],
        'blue': []
    }

    # Adjust the separation point
    for i in range(256):
        r, g, b, _ = RdBu(i)
        pos = i / 255.0
        if pos < separation:
            new_pos = pos / separation
        else:
            new_pos = (pos - separation) / (1 - separation)

        cdict['red'].append((new_pos, r, r))
        cdict['green'].append((new_pos, g, g))
        cdict['blue'].append((new_pos, b, b))

    # Sort the colormap segments by position
    for key in cdict:
        cdict[key].sort(key=lambda x: x[0])

    # Create and return the custom colormap
    return LinearSegmentedColormap('Custom_RdBu', cdict)
