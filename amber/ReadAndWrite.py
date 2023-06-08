import numpy as np
import pickle
def DoseOnWorld(file_path : str):
    # Get the current dose on the world
    # file_path is the path to the file containing the dose
    nn = np.array([])
    doses = np.array([])
    arr = np.loadtxt(file_path, skiprows=1, delimiter=",")
    num_voxels = len(np.unique(arr[:, 0]))
    for i in range(len(arr)):
        doses = np.append(doses, arr[i][3])
        n = arr[i][0] * num_voxels ** 2 + arr[i][1] * num_voxels + arr[i][2]
        nn = np.append(nn, n)
    if np.max(doses) != 0:
        doses = doses / np.max(doses)
    return nn, doses

def show_tumor_3D_solid(world, t):

    import pyvista as pv

    grid_shape = (world.number_of_voxels, world.number_of_voxels, world.number_of_voxels)
    cell_counts = np.empty(grid_shape)
    necro_counts = np.empty(grid_shape)

    # Fill the array with the number of tumor cells in each voxel
    for voxel in world.voxel_list:
        i, j, k = np.unravel_index(voxel.voxel_number, grid_shape)
        cell_counts[i, j, k] = voxel.number_of_tumor_cells()
        necro_counts[i, j, k] = voxel.number_of_necrotic_cells()


    # Create a vtkImageData object and assign the cell counts to it
    grid = pv.UniformGrid()
    grid2 = pv.UniformGrid()
    grid.dimensions = grid_shape
    grid2.dimensions = grid_shape
    hl = world.half_length
    voxel_side = 2*hl/world.number_of_voxels
    grid.origin = (-hl, -hl, -hl)  # The bottom left corner of the data set
    grid2.origin = (-hl, -hl, -hl)  # The bottom left corner of the data set
    grid.spacing = (voxel_side, voxel_side, voxel_side)  # These are the cell sizes along each axis
    grid2.spacing = (voxel_side, voxel_side, voxel_side)  # These are the cell sizes along each axis
    grid.point_data['tumor'] = cell_counts.flatten(order="F")  # Flatten the array!
    grid2.point_data['necro'] = necro_counts.flatten(order="F")  # Flatten the array!

    contour_values = [1,950]
    opacities = [0.2, 0.7]
    max = np.max(cell_counts)
    #remove values below max from contour_values
    contour_values = [x for x in contour_values if x < max]

    # Create a Plotter object
    plotter = pv.Plotter(shape=(1, 2))
    plotter.subplot(0, 0)
    plotter.add_mesh(grid.outline_corners(), color='k')
    contour_necro = grid2.contour([950])
    plotter.add_mesh(contour_necro, cmap='Reds', opacity=0.9, scalars='necro', clim=[0, 1000])


    for i, value in enumerate(contour_values):
        contour = grid.contour([value])
        plotter.add_mesh(contour, cmap='Greens', opacity= opacities[i], scalars='tumor')

    plotter.subplot(0, 1)
    plotter.add_mesh(grid.contour([1]), cmap='Greens', opacity= opacities[0], scalars='tumor')
    #cmap limits

    plotter.add_mesh_slice(grid, normal='x', cmap='Greens', scalars='tumor', clim=[0, 1000])
    plotter.add_mesh_slice(grid2, normal='y', cmap='Reds', scalars='necro', clim=[0, 1000])
    # for vessel in world.vasculature.list_of_vessels:
    #     if vessel.visible:
    #         path = vessel.path
    #         if len(path) > 1:
    #             origin = path[0]
    #             end = path[-1]
    #             line = pv.Line(origin, end)
    #             plotter.add_mesh(line, color='crimson', line_width=1)

    # plotter.export_html('3D_plot_' + str(t) + '.html', 'panel')
    plotter.show(window_size=(4000, 2000), screenshot='3D_plot_' + str(t) + '.png')


