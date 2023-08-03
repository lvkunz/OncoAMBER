import subprocess
import os

def RunTopasSimulation(file: str, cluster = True): #runs a topas simulation using the given file on MACOS. TODO: check if path are correct for your environment
    print('Running Topas Simulation')
    file = file + '.txt'
    working_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    working_dir = working_dir + '/'
    if cluster: run = 'topas ' + working_dir + file
    else: run = '/Applications/topas/bin/topas ' + working_dir + file

    #open current folder
    sequence_of_commands = ['export TOPAS_G4_DATA_DIR=/Applications/G4Data',
                            'export QT_QPA_PLATFORM_PLUGIN_PATH=/Applications/topas/Frameworks', #this is to set up the environment for topas
                            'cd ' + working_dir, #this is to change the directory to the working directory
                            run] #this is to run the topas simulation
    command = sequence_of_commands[0]
    for i in range(1, len(sequence_of_commands)):
        command = command + ' && ' + sequence_of_commands[i]
    print(command)
    subprocess.call(command, shell=True)

    print('Topas Simulation Done')


