from .config import params
import numpy as np
from .basis import Basis, OperatorString
import pickle

class SpectralProps(Basis):

    def __init__(self, solved_model = None, save_file = None):
        super().__post_init__()
        if solved_model:
            self.E, self.Psi = solved_model.E, solved_model.psi
        elif save_file:
            with open(save_file, 'rb') as f:
                self.E, self.Psi = pickle.load(f)
        else:
            raise FileNotFoundError('Model could not be found!')
        # Ground state energy
        self.E0 = sorted(self.E)[0].real

    def Z_tilde(self, beta: float):
        return np.sum( np.exp( - beta * (self.E.real - self.E0) ) )


class GroundStateProps(Basis):

    def __init__(self, solved_model = None, save_file = None):
        super().__post_init__()
        if solved_model:
            self.E, self.Psi = solved_model.E, solved_model.psi
        elif save_file:
            with open(save_file, 'rb') as f:
                self.E, self.Psi = pickle.load(f)
        else:
            raise FileNotFoundError('Model could not be found!')
        self.get_energy()
        self.get_ground_state()

    def get_energy(self, dtol = 1e-4):
        self.E0 = sorted(self.E)[0]
        self.GS_degeneracy = np.sum(np.abs(self.E - self.E0) < dtol)
        if self.GS_degeneracy > 1:
            print(f'Warning: ground state is degenerate, D = {self.GS_degeneracy}')

    def get_ground_state(self, state_idx = None):
        GS_idx = np.argmin(self.E)
        self.Psi0 = self.Psi[:, GS_idx]
        if state_idx != None:
            print(f'Energy of state {state_idx}: {self.E[state_idx]}')
            self.Psi0 = self.Psi[:, state_idx]

    def corr_func(self, op: OperatorString, beta = None):
        '''
        Compute the n-pt correlation function
        < c^\\pm_{i_1\\sigma_1} ... c^\pm_{i_n\\sigma_n} >
        '''
        corr = 0
        for n, alpha_n in enumerate(self.Psi0):
            phase_fac, f_n = self.get_fermion_matrix_element(n, op)
            if phase_fac != 0:
                alpha_f = self.Psi0[f_n]
                corr += np.conj(alpha_f) * alpha_n * phase_fac
        return corr

    def get_spectral_func(self, q: float, sigma: str, omega: float, eta: float):
        '''
        The spectral function A_q(\omega).
        - q = momentum slice (only in 1D for now --> 0D point)
        - sigma = spin ("up" or "down")
        - omega = frequency
        - eta = finite broadening (keep small)
        '''
        if not hasattr(self, 'c_0k'):
            self.get_fermion_matrix_elements()
        retarded_green = 0
        #for k in np.where(self.E > self.E0):
        for k in range(len(self.E)):
            dE_0k = self.E0   - self.E[k]
            dE_k0 = self.E[k] - self.E0
            # Discrete Fourier transform
            c_0k = 0
            c_k0 = 0
            for i, w_i in enumerate(self.sites):
                phase_factor = np.exp(1j * q * w_i.cellX)
                c_0k += phase_factor * self.c_0k[(i, sigma)][k]
                c_k0 += phase_factor * self.c_k0[(i, sigma)][k]
            # Normalize
            L = len(self.sites)
            c_0k *= 1/np.sqrt(L)
            c_k0 *= 1/np.sqrt(L)
            term1 = c_0k * np.conj(c_0k) / (omega + dE_0k + 1j * eta)
            term2 = c_k0 * np.conj(c_k0) / (omega + dE_k0 + 1j * eta)
            retarded_green += term1 - term2
        return -2 * retarded_green.imag

    def get_fermion_matrix_elements(self):
        c_k0 = dict()
        c_0k = dict()
        for i, w_i in enumerate(self.sites):
            for sigma in ['up', 'down']:
                print(f'Computing matrix element for site {i}, spin {sigma}')
                c_k0[(i, sigma)] = list()
                c_0k[(i, sigma)] = list()
                # First compute k0
                for k in range(len(self.E)):
                    Psi_k = self.Psi[:, k]
                    element = 0
                    for n, alpha_n in enumerate(self.Psi0):
                        ket_state = self.get_state(n)
                        ket_state.c_(-1, w_i, sigma)
                        if not ket_state.null:
                            f_n = self.get_label(ket_state)
                            phi_n = ket_state.phase
                            alpha_f = Psi_k[f_n]
                            element += np.conj(alpha_f) * alpha_n * np.exp(1j * phi_n)
                    c_k0[(i, sigma)].append(element)
                # Then compute 0k
                for k in range(len(self.E)):
                    Psi_k = self.Psi[:, k]
                    element = 0
                    for n, alpha_n in enumerate(Psi_k):
                        ket_state = self.get_state(n)
                        ket_state.c_(-1, w_i, sigma)
                        if not ket_state.null:
                            f_n = self.get_label(ket_state)
                            phi_n = ket_state.phase
                            alpha_f = self.Psi0[f_n]
                            element += np.conj(alpha_f) * alpha_n * np.exp(1j * phi_n)
                    c_0k[(i, sigma)].append(element)
        self.c_k0 = c_k0
        self.c_0k = c_0k


class FiniteTempProps(Basis):

    def __init__(self, solved_model = None, save_file = None):
        super().__post_init__()
        if solved_model:
            self.E, self.Psi = solved_model.E, solved_model.psi
        elif save_file:
            with open(save_file, 'rb') as f:
                self.E, self.Psi = pickle.load(f)
        else:
            raise FileNotFoundError('Model could not be found!')
        self.E0 = sorted(self.E)[0].real
        #print(self.E0)


#    def partition_func(self, beta: float):
#        return np.sum( np.exp( - beta * self.E ) )

    def Z_tilde(self, beta: float):
        return np.sum( np.exp( - beta * (self.E.real - self.E0) ) )

    def total_energy(self, beta: float):
        Z_tilde = self.Z_tilde(beta)
        return np.sum(self.E.real * np.exp( -beta * (self.E.real - self.E0))) / Z_tilde

    def corr_func(self, op: OperatorString, beta: float):
        '''
        Compute the n-pt correlation function
        < c^\\pm_{i_1\\sigma_1} ... c^\pm_{i_n\\sigma_n} >
        '''
        Z_tilde = self.Z_tilde(beta)
        #print(f'Z_tilde: {Z_tilde}')
        trace = 0
        for k, Psi_k in enumerate(self.Psi.T):
            boltz_fac = np.exp( - beta * (self.E[k].real - self.E0) )
            #print(f'BF: {boltz_fac}')
            for n, alpha_n in enumerate(Psi_k):
                phase_fac, f_n = self.get_fermion_matrix_element(n, op)
                if phase_fac != 0:
                    alpha_f = Psi_k[f_n]
                    trace += boltz_fac * np.conj(alpha_f) * alpha_n * phase_fac
        #print(f'Trace: {trace}')
        return trace/Z_tilde

