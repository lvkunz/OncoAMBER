import numpy as np

class Config:
    def __init__(self, config_dict):
        self.__dict__ = {**self.__dict__, **config_dict}

        if self.seed == -1:
            self.seed = np.random.randint(0, 1000000)
