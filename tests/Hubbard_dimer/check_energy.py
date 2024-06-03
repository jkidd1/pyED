from pyED.run import make_inputs, diag
from pyED.filo import cd
import numpy as np
from sys import stderr

def E_ref(U):
    # arXiv:0807.4878
    E_plus  = U/2 + np.sqrt((U/2)**2 + 4)
    E_minus = U/2 - np.sqrt((U/2)**2 + 4)
    return np.array(sorted([0, 0, 0, U, E_plus, E_minus]))

U_vals = np.linspace(0, 20, 100)
grid = {'@U': U_vals}
jobs = make_inputs(grid, write_job_list = False)

error = np.zeros((len(jobs), 6))
for i, (cfg, folder) in enumerate(jobs):
    U = cfg[0]
    ref = E_ref(U)
    with cd(f'calc/{folder}'):
        E_vals = np.sort(diag().E)
    for j, E in enumerate(E_vals):
        error[i, j] = E.real - ref[j]

with open('max_error', 'w') as f:
    f.write(f'{np.max(error)}')

if np.max(error) < 1e-10:
    print('Hubbard dimer test passed.', file = stderr)

else:
    print('Hubbard dimer test failed, check installation.', file = stderr)
