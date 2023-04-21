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

def read_config_file(file_name):
    config_dict = {}
    file_name = file_name + '.txt'

    with open(file_name, 'r') as f:
        for line in f:
            # Remove comments and strip leading/trailing whitespaces
            line = line.split('#')[0].strip()

            # Skip blank lines
            if not line:
                continue

            key, value = line.split(': ')

            # Check if the value is a number (integer or float)
            if value.replace('.', '', 1).isdigit() or value.lstrip('-').replace('.', '', 1).isdigit():
                if '.' in value:
                    config_dict[key] = float(value)
                else:
                    config_dict[key] = int(value)
            # Check if the value is a boolean
            elif value.lower() == 'true' or value.lower() == 'false':
                config_dict[key] = bool(value.lower() == 'true')
            # Otherwise, treat the value as a string
            else:
                config_dict[key] = value

    return config_dict




