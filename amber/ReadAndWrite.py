import numpy as np

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




