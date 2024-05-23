from dataclasses import dataclass
from .lattice import Lattice
from .basis import OperatorString
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh
import numpy as np
from .config import params
import pickle

@dataclass
class ModelBase(Lattice):
    
    def __post_init__(self):
        super().__post_init__()
        print(f'Hilbert space dimension = {self.dim}')
        self.H = csr_matrix((self.dim, self.dim), dtype = complex)

    def add_coupling(self, coupling: float, op: OperatorString):
        row_idx = list()
        col_idx = list()
        H_data  = list()
        for col in range(self.dim):
            phase_fac, f_n = self.get_fermion_matrix_element(col, op)
            if phase_fac != 0:
                col_idx.append(col)
                row_idx.append(f_n)
                H_data.append(coupling * phase_fac)
        self.H += coo_matrix((H_data, (row_idx, col_idx)), shape = (self.dim, self.dim), dtype = complex).tocsr()

    def add_constant(self, constant: float):
        row_idx = list(range(self.dim))
        col_idx = list(range(self.dim))
        H_data  = [constant for _ in range(self.dim)]
        self.H += coo_matrix((H_data, (row_idx, col_idx)), shape = (self.dim, self.dim), dtype = complex).tocsr()
        
    def solve(self):
        print('Solving...')
        if params.solver == 'full':
            #self.E, self.psi = eigh(self.H.toarray())
            import numpy as np
            self.E, self.psi = np.linalg.eig(self.H.toarray())
        elif params.solver == 'Lanczos':
            # Get only the first `k` eigenstates
            num_eigs = params.num_eigs if params.num_eigs else self.dim // 2
            self.E, self.psi = eigsh(self.H, k = num_eigs, which = 'SA')
        else:
            raise ValueError(f'{params.solver} is not a valid solver!')

    def save_result(self, fname = 'result.p'):
        with open(fname, 'wb') as f:
            pickle.dump((self.E, self.psi), f)
        print(f'Result saved to {fname}')

@dataclass
class nnHopping(ModelBase):

    def __post_init__(self):
        super().__post_init__()
        self.add_Ht()

    def add_Ht(self):
        print('Building hopping matrix...')
        for w_i, w_j in self.nn_pairs:
            for sigma in ['up', 'down']:
                CdC = OperatorString([+1, -1], [w_i, w_j], [sigma, sigma])
                self.add_coupling(-params.t, CdC)

@dataclass
class nnHubbard(nnHopping):
    
    def __post_init__(self):
        super().__post_init__()
        self.add_HU()

    def add_HU(self):
        print('Building Hubbard matrix...')
        for w_i in self.sites:
            nu_nd = OperatorString([+1, -1, +1, -1], [w_i, w_i, w_i, w_i], ['up', 'up', 'down', 'down'])
            self.add_coupling(params.U, nu_nd)

@dataclass
class nnHubbardPH(nnHubbard):
    
    def __post_init__(self):
        super().__post_init__()
        self.add_PH()

    def add_PH(self):
        print('Adding particle-hole correction...')
        mu = -params.U / 2
        for w_i in self.sites:
            for sigma in ['up', 'down']:
                n = OperatorString([+1, -1], [w_i, w_i], [sigma, sigma])
                self.add_coupling(mu, n)
        L = params.LX * params.LY * params.N_sub
        self.add_constant(params.U * L / 4)
