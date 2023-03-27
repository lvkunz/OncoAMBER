import subprocess


def RunTopasSimulation(file: str):
    print('Running Topas Simulation')
    file = file + '.txt'
    run = '/Applications/topas/bin/topas ' + file

    #open current folder
    sequence_of_commands = ['export TOPAS_G4_DATA_DIR=/Applications/G4Data',
                            'export QT_QPA_PLATFORM_PLUGIN_PATH=/Applications/topas/Frameworks',
                            'cd /Users/louiskunz/PycharmProjects/TumorVoxelBased-Git/TopasSimulation',
                            run]
    command = sequence_of_commands[0] + ' && ' + sequence_of_commands[1] + ' && ' + sequence_of_commands[2] + ' && ' + sequence_of_commands[3]
    print(command)
    subprocess.call(command, shell=True)

    print('Topas Simulation Done')


def read_config_file(file_name):
    config_dict = {}

    with open(file_name, 'r') as f:
        for line in f:
            key, value = line.strip().split(': ')

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