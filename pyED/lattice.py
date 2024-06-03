from dataclasses import dataclass
from .basis import Basis, WannierState
from importlib import import_module
from .config import params
import sys
import numpy as np

@dataclass
class Lattice(Basis):
    
    def __post_init__(self):
        super().__post_init__()
        bonds = ['nearest_neighbors', 'next_nearest_neighbors']
        if hasattr(params, 'bonds'):
            bonds += params.bonds
        for neighbor_type in bonds:
            getter = f'get_{neighbor_type}'
            # custom lattice
            try:
                try:
                    neighbor_func = getattr(import_module(params.lattice_type), getter)
                except AttributeError:
                    neighbor_func = lambda w_j: []
            # built-in lattice
            except ModuleNotFoundError:
                try:
                    built_in = getattr(import_module('.lattice', package = 'pyED'), params.lattice_type)()
                    try:
                        neighbor_func = getattr(built_in, getter)
                    except AttributeError:
                        neighbor_func = lambda w_j: []
                except AttributeError:
                    print(f'ERROR: \'{params.lattice_type}\' is not a built-in lattice type.\nFor custom lattices, the script name must match lattice_type in params.yaml.')
                    sys.exit()
            setattr(self, getter, neighbor_func)

        if hasattr(self, 'get_nearest_neighbors'):
            self.get_nn_pairs()

    def get_nn_pairs(self):
        '''
        Write a list of all nearest neighbor pairs <ij> on the lattice.
        '''
        self.nn_pairs = list()
        for w_j in self.sites:
            for w_i in self.get_nn(w_j):
                self.nn_pairs.append((w_i, w_j))
    
    def get_nn(self, w_j: WannierState):
        nn = self.get_nearest_neighbors(w_j)
        self.apply_BC(nn)
        return nn

    def get_nnn(self, w_j: WannierState):
        nnn, nu = self.get_next_nearest_neighbors(w_j)
        self.apply_BC(nnn)
        return nnn, nu

    def apply_BC(self, neighbor_list):
        for direction in ('X', 'Y'):
            side_length = getattr(params, f'L{direction}')
            # OBC --> remove bonds outside supercell
            if getattr(params, f'open_{direction}'):
                for n in neighbor_list.copy():
                    cell = getattr(n, f'cell{direction}')
                    if cell < 0 or cell > side_length - 1:
                        neighbor_list.remove(n)
            # PBC --> identify the cells modulo L
            else:
                for n in neighbor_list:
                    cell = getattr(n, f'cell{direction}') % getattr(params, f'L{direction}')
                    setattr(n, f'cell{direction}', cell)
                    n.__post_init__()
        
@dataclass
class Chain:
    def get_nearest_neighbors(self, w_j: WannierState):
        return [WannierState(w_j.cellY, w_j.cellX + dX, 0) for dX in (-1, +1)]

@dataclass
class ZigzagHoneycomb:
    def get_nearest_neighbors(self, w_j: WannierState):
        nn = list()
        Y, X, subl = w_j.cellY, w_j.cellX, w_j.sublattice
        nn.append(WannierState(Y, X, (subl + 1) % 2))
        nn.append(WannierState(Y + (-1)**subl, X, (subl + 1) % 2))
        nn.append(WannierState(Y, X - (-1)**subl, (subl + 1) % 2))
        return nn
    def get_next_nearest_neighbors(self, w_j: WannierState):
        nnn = list()
        Y, X, subl = w_j.cellY, w_j.cellX, w_j.sublattice
        nnn.append(WannierState(Y + 1, X, subl)) # nu = -1
        nnn.append(WannierState(Y - 1, X, subl)) # nu = +1
        nnn.append(WannierState(Y, X + 1, subl)) # nu = -1
        nnn.append(WannierState(Y, X - 1, subl)) # nu = +1
        nnn.append(WannierState(Y + 1, X + 1, subl)) # nu = +1
        nnn.append(WannierState(Y - 1, X - 1, subl)) # nu = -1
        nu = np.array([-1, +1, -1, +1, +1, -1])
        if subl == 1:
            nu *= -1
        return nnn, nu

@dataclass
class Hexagon:
    def get_nearest_neighbors(self, w_j: WannierState):
        nn = list()
        Y, X, subl = w_j.cellY, w_j.cellX, w_j.sublattice
        for ds in (-1, +1):
            nn.append( WannierState( Y, X, (subl + ds)%6 ) )
        return nn
    def get_next_nearest_neighbors(self, w_j: WannierState):
        nnn = list()
        Y, X, subl = w_j.cellY, w_j.cellX, w_j.sublattice
        for ds in (-2, +2):
            nnn.append( WannierState( Y, X, (subl + ds)%6 ) )
        nu = [+1, -1]
        return nnn, nu

@dataclass
class Square:
    def get_nearest_neighbors(self, w_j: WannierState):
        nn = list()
        for dR in (-1, +1):
            nn.append(WannierState(w_j.cellY + dR, w_j.cellX, 0))
            nn.append(WannierState(w_j.cellY, w_j.cellX + dR, 0))
        return nn

