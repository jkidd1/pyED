from dataclasses import dataclass
import numpy as np
from .config import params
from typing import Optional
from itertools import combinations
from copy import deepcopy

@dataclass
class WannierState:
    '''
    cellY = layer index, i.e., unit cell along y-direction (0, ..., LY - 1)
    cellX = unit cell along x-direction (0, ..., LX - 1)
    sublattice = 0, ..., N_sub - 1
    '''
    cellY: int
    cellX: int
    sublattice: int
        
    def __post_init__(self):
        self.get_position()
        self.get_index()
    
    def get_position(self):
        '''
        Find the position of the site in the occupation array. 
        '''
        row = self.cellY
        col = params.N_sub * self.cellX + self.sublattice
        self.position = (row, col)
    
    def get_index(self):
        '''
        Find the (spatial) index of the site.
        '''
        self.index = params.N_sub * (params.LX * self.cellY + self.cellX) + self.sublattice


@dataclass
class ProductState:
    '''
    Each spin state is an `LY` by `N_sub * LX` array of occupation numbers (0 or 1).
    --> We follow the convention of `WannierState.get_position`.
    Each basis state is normally ordered.
    --> We follow the convention of `WannierState.get_index`.
    '''
    up: np.ndarray
    down: np.ndarray
    null: Optional[bool] = False
    phase: Optional[float] = 0
    
    def __post_init__(self):
        self.assign_string()
        self.get_occupied_orbitals()
    
    def assign_string(self):
        self.up_string = ''.join(self.up.flatten().astype(str))
        self.down_string = ''.join(self.down.flatten().astype(str))
    
    def get_occupied_orbitals(self):
        self.occ_up = np.array([i for i, ni in enumerate(self.up_string) if ni == '1'])
        self.occ_down = np.array([i for i, ni in enumerate(self.down_string) if ni == '1'])
    
    def c_(self, Delta_nj: int, w_j: WannierState, sigma: str):
        '''
        Apply raising (`Delta_nj` = +1) or lowering (`Delta_nj` = -1) operator to the state (in-place).
        Include phase (+1 or -1) based on normal ordering.
        Becomes a null state for c_j\ket{n_j = 0} or c_j^\dagger\ket{n_j = 1}.
        '''
        # Record current occupied states
        occ_states = {'up': self.occ_up, 'down': self.occ_down}
        # Update occupation numbers
        getattr(self, sigma)[w_j.position] += Delta_nj
        self.__post_init__()
        # Update phase
        if sigma == 'down':
            # First, anti-commute with all $c^\dagger_{i\uparrow}$
            self.phase += len(occ_states['up']) * np.pi
        # Anti-commute with $c^\dagger_{i\sigma}$ until normal ordering is reached
        self.phase += sum(occ_states[sigma] < w_j.index) * np.pi
        # If necessary, flag state as null
        if np.any(getattr(self, sigma) < 0) or np.any(getattr(self, sigma) > 1):
            self.null = True

@dataclass
class OperatorString:
    Deltas: list
    sites:  list
    spins:  list
    
    def __iter__(self):
        return iter(zip(self.Deltas, self.sites, self.spins))

@dataclass
class Basis:
        
    def __post_init__(self):
        self.get_sites()
        self.Lin = list()
        self.generate_Lin_table()
    
    def get_sites(self):
        '''
        The list of all spatial sites (i.e., `WannierStates`) on the lattice.
        '''
        self.sites = list()
        for cellY in range(params.LY):
            for cellX in range(params.LX):
                for sublattice in range(params.N_sub):
                    self.sites.append(WannierState(cellY, cellX, sublattice))  
        self.num_sites = len(self.sites)
        
    def generate_Lin_table(self):
        '''
        Write a list of basis states, ordered by their binary representation.
        Depends on whether conserved quantities (N and/or Sz) are considered.
        '''
        L = params.LX * params.LY * params.N_sub

        cons_N  = (params.N  != None)
        cons_Sz = (params.Sz != None)
        
        # There are 4 cases total.

        # 1. Both N and Sz conserved
        if cons_N and cons_Sz:
            self.update_Lin_table(params.N, params.Sz)

        # 2. N conserved, Sz not conserved
        if cons_N and (not cons_Sz):
            Sz_min  = abs(params.N - L) - L
            Sz_max  = L - abs(params.N - L)
            Sz_vals = range(Sz_min, Sz_max + 1, 2)
            for Sz in Sz_vals:
                self.update_Lin_table(params.N, Sz)

        # 3. N not conserved, Sz conserved
        if (not cons_N) and cons_Sz:
            N_min  = abs(params.Sz)
            N_max  = 2*L - abs(params.Sz)
            N_vals = range(N_min, N_max + 1, 2)
            for N in N_vals:
                self.update_Lin_table(N, params.Sz)

        # 4. Neither N nor Sz conserved
        if (not cons_N) and (not cons_Sz):
            N_vals = range(2*L + 1)
            for N in N_vals:
                Sz_min  = abs(N - L) - L
                Sz_max  = L - abs(N - L)
                Sz_vals = range(Sz_min, Sz_max + 1, 2)
                for Sz in Sz_vals:
                    self.update_Lin_table(N, Sz)

        # Optionally write basis file
        if params.verbose == True:
            with open('basis', 'w') as f:
                f.write('\n'.join(self.Lin))
        
        # Hilbert space dimension
        self.dim = len(self.Lin)
    

    def update_Lin_table(self, N, Sz):
        if (N + Sz) % 2 != 0 or (N - Sz) % 2 != 0:
            err_msg = 'Number of fermions in each spin sector must be an integer!'
            raise ValueError(err_msg)
        N_up   = (N + Sz) // 2
        N_down = (N - Sz) // 2
        self.Lin += self.get_binary_strings(N_up, N_down)


    def get_binary_strings(self, N_up, N_down):
        '''
        Mathematical function to generate all possible spin configurations on the lattice.
        '''
        bin_strings = list()
        for x in combinations([1<<i for i in range(self.num_sites)], N_up):
            for y in combinations([1<<i for i in range(self.num_sites)], N_down):
                bin_strings.append(f'{sum(x):0{self.num_sites}b}' + f'{sum(y):0{self.num_sites}b}')
        return sorted(bin_strings, key = lambda b: int(b, 2))
    
    def get_label(self, state: ProductState):
        '''
        Look up basis state on the Lin table and assign a unique label.
        '''
        full_string = state.up_string + state.down_string
        return self.Lin.index(full_string)

    def get_state(self, label: int):
        '''
        Inverse of `get_label`, i.e., find a `ProductState` from its basis label.
        '''
        full_string = self.Lin[label]
        states = list()
        for string in (full_string[:self.num_sites], full_string[self.num_sites:]):
            state = np.array(list(string), dtype = int).reshape(params.LY, params.N_sub * params.LX)
            states.append(state)
        return ProductState(*states)

    def get_fermion_matrix_element(self, label: int, op: OperatorString):
        ket_state = self.get_state(label)
        for idx_set in op:
            ket_state.c_(*idx_set)
        if ket_state.null:
            return 0, 0
        else:
            f_n = self.get_label(ket_state)
            phi_n = ket_state.phase
            phase_fac = np.exp(1j * phi_n)
            return phase_fac, f_n

