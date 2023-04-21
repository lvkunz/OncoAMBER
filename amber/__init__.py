from .config import *
from .config_instance import *
from .ReadAndWrite import *

#get path to amber folder
import os
amber_dir = os.path.abspath(os.path.dirname(__file__))

config_path = os.path.join(amber_dir, 'CONFIG_default')
config_dict = read_config_file(config_path)
config = Config.from_dict(config_dict)
set_config_instance(config)

from .BasicGeometries import *
from .Cell import *
from .Voxel import *
from .Vessel import *
from .World import *
from .Process import *
from .Terminal import *
from .RunTopas import *
from .ScalarField import *
from .BetaDistributionCalibration import *