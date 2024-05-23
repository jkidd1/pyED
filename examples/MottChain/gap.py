from pyED.config import params as p
from pyED.filo import cd, reader
from importlib import reload
import pandas as pd
import matplotlib.pyplot as plt

df = {'charge': pd.DataFrame(), 'spin': pd.DataFrame()}
folders = reader('job_files/job.list').contents

# Read energy
E0 = dict()
for folder in folders:
    with cd(f'calc/{folder}'):
        p.reload()
        U, L, N, Sz = p.U, p.LX, p.N, p.Sz
        E0[(U, L, N, Sz)] = float(reader('energy').contents[0])

# Compute charge gap and spin gap
U_vals = sorted(list(set(key[0] for key in E0.keys())))
L_vals = sorted(list(set(key[1] for key in E0.keys())))
for L in L_vals:
    df['charge'][f'L={L}'] = [E0[(U, L, L + 1, 1)] + E0[(U, L, L - 1, 1)] - 2 * E0[(U, L, L, 0)] for U in U_vals]
    df['spin'][f'L={L}'] = [E0[(U, L, L, 2)] - E0[(U, L, L, 0)] for U in U_vals]

# Plot results
def plot_gap(df, ylabel, title, fname):
    fig, ax = plt.subplots()
    for L in L_vals:
        ax.plot(U_vals, df[f'L={L}'], marker = 'o', label = f'L = {L}')
    ax.legend()
    ax.set_ylabel(ylabel, fontsize = 16, labelpad = 15, rotation = 0)
    ax.set_xlabel('U/t', fontsize = 16)
    ax.set_title(title, fontsize = 16)
    fig.tight_layout()
    fig.savefig(fname, dpi = 600)

plot_gap(df['charge'], r'$\Delta_c$', 'Charge gap (OBC)', 'figures/charge_gap.png')
plot_gap(df['spin'], r'$\Delta_s$', 'Spin gap (OBC)', 'figures/spin_gap.png')
