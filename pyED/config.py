from yaml import safe_load

def load_params(container):
    with open('params.yaml', 'r') as f:
        for key, value in safe_load(f).items():
            setattr(container, key, value)

class params:
    def reload():
        load_params(params)
            
try:
    load_params(params)
except FileNotFoundError:
    print('WARNING: module loaded without params.yaml file.')
