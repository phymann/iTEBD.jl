# itebd.jl

Classical two-dimensional Ising model in a longitudinal field, given by the following Hamiltonian,
$$
H = -J \sum_{\langle i,j \rangle} \sigma_i^z \sigma_j^z - h \sum_i \sigma_i^z
$$

is solved numerically using the iTEBD algorithm proposed in [PRB 78, 155117 (2008)](https://link.aps.org/doi/10.1103/PhysRevB.78.155117).

## Usage

Start from `runiTEBD.jl`