import src.ReadAndWrite as rw

config = rw.read_config_file('src/CONFIG_default')

def set_config_instance(config_instance):
    global config
    config = config_instance
