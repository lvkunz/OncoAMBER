import numpy as np
import os

class Config:
    def __init__(self, config_dict):
        self.__dict__ = {**self.__dict__, **config_dict}

        if self.seed == -1:
            self.seed = np.random.randint(0, 1000000)

        self.working_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        print('Working directory: ' + self.working_directory)