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


