import subprocess
import os

def RunTopasSimulation(file: str, cluster = True):
    print('Running Topas Simulation')
    file = file + '.txt'
    working_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    working_dir = working_dir + '/'
    if cluster: run = 'topas ' + working_dir + file
    else: run = '/Applications/topas/bin/topas '+ working_dir + file

    #open current folder
    sequence_of_commands = ['export TOPAS_G4_DATA_DIR=/Applications/G4Data',
                            'export QT_QPA_PLATFORM_PLUGIN_PATH=/Applications/topas/Frameworks',
                            'cd ' + working_dir,
                            run]
    command = sequence_of_commands[0] + ' && ' + sequence_of_commands[1] + ' && ' + sequence_of_commands[2] + ' && ' + sequence_of_commands[3]
    print(command)
    subprocess.call(command, shell=True)

    print('Topas Simulation Done')


