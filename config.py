import numpy as np

class Config:
    def __init__(self, config_dict):

        for key, value in config_dict.items():
            setattr(self, key, value)

        if self.seed == -1:
            self.seed = np.random.randint(0, 1000000)
    @classmethod
    def from_dict(cls, config_dict):
        return cls(config_dict)
