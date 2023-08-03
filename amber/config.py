import numpy as np
import os

class Config: #this is a class that reads the config file and stores the parameters as attributes of the class
    def __init__(self, file_name):

        config_dict = self.read_config_file(file_name)

        self.__dict__ = {**self.__dict__, **config_dict}

        ## RULES CAN BE ADDED HERE TO CHECK THE PARAMETERS, FOR EXAMPLE:
        if self.seed == -1: #if the seed is -1, generate a random seed
            self.seed = np.random.randint(0, 1000000)

        ##ADDITIONAL ATTRIBUTES NOT SPECIFIED IN CONFIG FILE CAN BE ADDED HERE, FOR EXAMPLE:
        self.working_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        print('Working directory: ' + self.working_directory)

    def read_config_file(self, file_name):
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