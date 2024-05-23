<h1 align="center">Python Exact Diagonalization (pyED)</h1>
  
- [Overview](#python-exact-diagonalization-pyed)
- [Installation](#how-to-install)
- [Pre-defined systems](#pre-defined-systems)
  - [Lattices](#lattices)
  - [Models](#models)
- [Basic usage](#how-to-use)
  - [Calculation parameters](#calculation-parameters)
  - [Minimal run](#minimal-run)
  - [Job array](#job-array)
  - [Custom lattice](#custom-lattice)
  - [Custom model](#custom-model)


# Overview
_pyED_ is a simple toolkit for numerically diagonalizing an interacting spinful fermion model on a small 2D lattice. The model always has the form
### $$\hat H = \hat H_0 + \hat H_{\text{int}},$$
where
### $$\hat H_0 = \sum \ T(i\alpha,\ j\beta)\  c^\dagger({i\alpha}) \ c({j\beta}),$$
### $$\hat H_{\text{int}} = \sum \ V(i\alpha,\ j\beta,\ k\gamma, l\sigma)\  c^\dagger({i\alpha}) \ c^\dagger({j\beta}) \ c({k\gamma}) \ c({l\sigma}), $$
Latin letters label spatial Wannier orbitals (spanning the lattice), and Greek letters label spin ($\uparrow$ or $\downarrow$). 
### Features
- [ ] Compute observables at arbitrary temperature
- [X] Built-in model/lattice types
- [X] Support for custom models/lattices
- [X] Intuitive workflows via command line interface (CLI)

# How to install
Navigate to your preferred installation folder. Then, run the following:
```bash
git clone https://github.com/jkidd1/pyED.git
cd pyED
bash install.bash
```
To check the installation, you can run `bash tests/run_tests.bash`

# Pre-defined systems
### Lattices
- `Chain`: an effective 1D system that couples only along the $x$-direction. One sublattice.
- `Square`: like `Chain`, except that coupling is allowed in the $y$-direction. One sublattice.
- `ZigzagHoneycomb`: Zigzag edge termination for the hexagonal graphene-type structure, i.e., rhombic unit cell with two sublattices.

### Models
- `nnHopping`: free hopping chain ($H_{\text{int}} = 0$) with nearest neighbor coupling $t$:
### $$\hat H_t = -t\sum_{\langle ij \rangle} \sum_{\sigma} \ c^\dagger({i\sigma}) \ c({j\sigma}) $$
- `nnHubbard`: adds Hubbard interaction $U$ to `nnHopping`, where
### $$\hat H_U = U \sum_i \ c^\dagger({i\uparrow}) \ c({i\uparrow}) \ c^\dagger({i\downarrow}) \ c({i\downarrow})$$

# How to use
### Calculation parameters
A basic run requires only one input file, `params.yaml`, containing the following parameters.
|Category| Parameter |Description |
|-| ----------- | ----------- |
|**Run**|  str: _**solver**_ | The method used to diagonalize the Hamiltonian. Options are _**full**_ or _**Lanczos**_. <ul><li>_**full**_: convert sparse matrix to a NumPy array and call `np.linalg.eig`. This has poor scaling, so it is not recommended unless the total number of sites $L\leq 6$.</li><li> _**Lanczos**_: use the sparse matrix algorithm from `scipy.sparse.linalg.eigsh`.</li></ul> |
|| int: _**num_eigs**_ | In the _**Lanczos**_ solver, only the _**num_eigs**_ lowest-energy eigenstates will be found. If set to `null`, then the number of eigenstates is automatically set to half the dimension of the Hilbert space. |
|| float: _**beta**_ | Inverse temperature. For ground state properties, this can be ignored. |
|**Lattice**| str: _**lattice_type**_ | Name of either (1) a built-in type or (2) the Python script containing the neighbor function(s). _Do not include '.py' at the end!_ |
| |  int: _**LX**_  |  Number of unit cells in the $x$-direction. |
| |  int: _**LY**_  |  Number of unit cells in the $y$-direction. |
| |  int: _**N_sub**_ | Number of sublattices per unit cell. |
| **Model** | str: _**model_type**_ | Name of either (1) a built-in type or (2) the Python script containing the couplings. _Do not include '.py' at the end!_ |
| | float: _[couplings]_ | The values for all coupling parameters present in the model, e.g., _**t**_ and _**U**_ for `nnHubbard`. |
| **Conserved quantities** |  int: _**N**_  | Number of fermions, i.e., $N_{\uparrow} + N_{\downarrow}$. For half-filling, set equal to $L_X \cdot L_Y \cdot N_{\text{sub}}$.  |
| | int: _**Sz**_ | Net spin, i.e., $N_{\uparrow} - N_{\downarrow}$. Options are the following: <ul><li>If the model conserves $S_z$, set to one of $\\{-N, \ -N + 2, \  \dots,\  N-2,\  N\\}.$</li><li>If the model does not conserve $S_z$, set to `null`. _The Hamiltonian will include matrix elements that couple different sectors of_ $S_z$. _Be mindful of system size scaling!_ </li></ul>|
| **Boundary conditions** |  bool: _**open_X**_ | If `True`, remove all bonds between cell $0$ and cell $L_X-1$. |
|  |  bool: _**open_Y**_ | If `True`, remove all bonds between cell $0$ and cell $L_Y-1$. |

### Minimal run
To diagonalize the model and (optionally) save the eigenstates to a file, run 
```bash 
pyED -s result.p
```
This is equivalent to the Python script
```python
from pyED.run import diag
solved_model = diag('result.p')
```
### Job array
To instead perform a high-throughput calculation with different sets of parameters, save `params.yaml` in a folder called `template`. Tag each variable parameter using `@`, i.e.,
```yaml
# Variables
var1: @var1
var2: @var2
```
Then, specify the parameter grid using `prepED`, e.g.,
```bash
prepED  -var1  1 2 3  -var2  4.5 6.7
```
This is equivalent to the Python script
```python
from pyED.run import make_inputs
param_grid = {'@var1': [1, 2, 3], '@var2': [4.5, 6.7]}
make_inputs(param_grid)
```
The folder names will be written to the file `job.list`. If you are using SLURM, you can run the calculation in parallel via job arrays. The job script should have the following basic structure (where `N` is the number of folders):
```bash
#!/bin/bash
#SBATCH --array=1-N
#[other SLURM settings]
subfolder=$(awk "NR==$SLURM_ARRAY_TASK_ID" job.list)
cd calc/$subfolder
pyED -s result.p
```

### Custom lattice
In the code, lattice sites are instantiated using the `WannierState` class of `pyED.basis`. This class takes three positional arguments [^1]: 
1. `cellY`, the $y$-coordinate of the unit cell,
2. `cellX`, the $x$-coordinate of the unit cell, and
3. `sublattice`, the sublattice index of the site.

[^1]: Be mindful of the convention: the _y_-coordinate comes before the _x_-coordinate.

To customize the lattice, you must specify all bonds relevant to the model being studied. Typically, this might only include nearest neighbor (n.n.) pairs. In that case, the lattice is nothing but a prescription for obtaining all n.n. of a given `WannierState`. In the working directory [^2], create the following Python script (the name should match `params.lattice_type`):
[^2]: If using a `prepED` workflow, place the script in the `template` folder.
```python
from pyED.basis import WannierState
def get_nearest_neighbors(w_i: WannierState):
    nn = list()
    # Define the lattice here.
    return nn
```
For example, the built-in `Square` lattice is written as follows:
```python
for dR in (-1, +1):
    nn.append(WannierState(w_i.cellY + dR, w_i.cellX, 0))
    nn.append(WannierState(w_i.cellY, w_i.cellX + dR, 0))
```

### Custom model
Similar to custom lattices, custom models require a separate Python script. Create a script with the same name as `params.model_type` [^3]. Start by importing the following:
[^3]: If the name conflicts with a built-in model, the custom script takes priority.
```python
from dataclasses import dataclass
from pyED.config import params
from pyED.basis import OperatorString
from pyED.model import Parent
```
In the last line, `Parent` is a placeholder for one of the following:
- If you are defining the model from scratch, import `ModelBase`
- If your model can be written as one of the built-in Hamiltonians plus additional terms, instead import the inital model as a base (see [pre-defined models](#models)).
  
Having imported the correct `Parent`, define a new model class that inherits it, with the following basic structure:
```python
@dataclass
class MyModel(Parent):
    def __post_init__(self):
        self.add_term()
    def add_term():
        # Define the model here.
```
For example, suppose we wanted to add a staggered potential $\Delta$ to the 1D Hubbard model, i.e.,
### $$\hat H = \hat H_t + \hat H_U + \hat H_\Delta$$
where
### $$\hat H_\Delta = \Delta \sum_i \sum_\sigma \ (-1)^i \ c^\dagger({i\sigma}) \ c({i\sigma})  $$

The resulting model is sometimes called the ionic Hubbard model, or IHM. In `params`, we set `N_sub` equal to 2. We then define the model like so:
```python
@dataclass
class IHM(nnHubbard):
    def __post_init__(self):
        self.add_term()
    def add_term():
        for w_i in self.sites:
            for sigma in ['up', 'down']:
                sub = w_i.sublattice
                E_i = params.Delta * (-1)**sub
                n_i = OperatorString([+1, 1], [w_i, w_i], [sigma, sigma])
                self.add_coupling(E_i, n_i)
```
A full pyED calculation can be found in `examples/IHM`.

# Contact
For questions about how to use the software or to request that a certain feature be added, please post an issue to GitHub or send me an email: [jkidd@tulane.edu](mailto:jkidd@tulane.edu)
