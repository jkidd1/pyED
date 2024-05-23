from pyED.filo import cd, writer
from itertools import product
from importlib import import_module
import os
import numpy as np

def make_inputs(grid, calc_dir = os.getcwd() + '/calc', template_dir = os.getcwd() + '/template', write_job_list = True):
    '''
    Generate params.yaml automatically from a template file.
    Useful when running high-throughput calculations.
    `grid` is a `dict` of (str, list) pairs with the following format:
    {'@L': [L1, L2, ...]}
    See examples/Chain for an example.
    '''
    labels = list(grid.keys())
    values = [grid[label] for label in labels]
    calc_configs = list(product(*values))
    folder_pattern = '/'.join([f'{label[1:]}={label}' for label in labels])
    job_list = []
    if not os.path.exists(calc_dir):
        os.makedirs(calc_dir)
    with cd(calc_dir):
        for params in calc_configs:
            folder = folder_pattern
            for i, label in enumerate(labels):
                if type(params[i]) == float:
                    folder = folder.replace(label, f'{params[i]:0.6f}')
                else:
                    folder = folder.replace(label, str(params[i]))
            job_list.append(folder)
            if not os.path.exists(folder):
                os.makedirs(folder)
            with cd(folder):
                os.system(f'cp {template_dir}/* .')
                for i, label in enumerate(labels):
                    if type(params[i]) == float:
                        writer('params.yaml').replace(label, f'{params[i]:0.6f}')
                    else:
                        writer('params.yaml').replace(label, str(params[i]))
    if write_job_list:
        with open('job.list', 'a') as f:
            f.write('\n'.join(job_list) + '\n')
        job_folder = 'slurmIO'
        if not os.path.exists(job_folder):
            os.makedirs(job_folder)
        os.system(f'mv job.list {job_folder}')
    return list(zip(calc_configs, job_list))

def diag(result_path = None):
    '''
    Diagonalize a model on the lattice.
    '''
    from pyED.config import params
    try:
        # custom model type
        Model = getattr(import_module(params.model_type), params.model_type)
    except ModuleNotFoundError:
        # built-in model type
        Model = getattr(import_module('.model', package = 'pyED'), params.model_type)
    if not result_path:
        print('WARNING: only energies will be saved!')
        print('To save eigenstates, use -s ... [cli mode] or diag(...) [script mode]')
        m = Model()
        m.solve()
    else:
        m = Model()
        m.solve()
        m.save_result(result_path)

    with open('energy', 'w') as f:
        f.write('\n'.join(np.sort(m.E).astype(str)))
    print('Energy saved to file.')
    
    return m

