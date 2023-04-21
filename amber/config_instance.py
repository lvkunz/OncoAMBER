import amber.ReadAndWrite as rw

config = None

def set_config_instance(config_instance):
    global config
    config = config_instance
