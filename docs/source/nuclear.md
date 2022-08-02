# Nuclear Structure Input

Users must provide nuclear one-body density matrix elements of the form:
\begin{equation}
    \rho_{K,T}^{\Psi}(a,b) = \langle \Psi| [\hat{c}_a^\dagger \hat{c}_b]_{K,T}|\Psi \rangle ,
\end{equation}
where $\Psi$ is the nuclear-target wave function and $\hat{c}^\dagger$,
$\hat{c}$ are the one-body creation, destruction operators. The matrix elements
must be stored in a file in a standard format produced by shell-model codes like
BIGSTICK.

## List of nuclear targets

Table of nuclear data we include with the program. Each corresponds to a (.dres)
density matrix file. The source indicates the nuclear Hamiltonian that was used
to generate the wave function data:

| Nuclei | Isotopes | Source |
| -- | -- | ------ |
| He | 4 | [@PhysRevC.68.041001;@shirokov2016n3lo]|
| C  | 12 |  [@cohen1965effective; @PhysRevC.68.041001; @shirokov2016n3lo]|
| F | 19 | [@PhysRevC.74.034315]|
| Si | 28, 29 | [@PhysRevC.74.034315]|
| Ar | 40 |  [@PhysRevC.86.051301]|
| Ge | 70, 72, 73, 74, 76 | [@PhysRevC.80.064323] |
| I | 127 |  [@GCN5082], used in [@PhysRevLett.100.052503;@PhysRevC.82.064304]|
| Xe | 128, 129, 130, 131, 132, 134, 136 | [@GCN5082], used in [@PhysRevLett.100.052503; @PhysRevC.82.064304] |



## Filling core orbitals for phenomenological interactions
Since standard one-body density matrices in phenomenological model spaces
contain only matrix elements for orbitals in the valence space, it is necessary
to infer the matrix elements for the core orbitals. Our code does this by
default, but the user can disable this option using the `fillnuclearcore`
control word.

For phenomenological interactions one typically has a `frozen' core of nucleons
which do not participate in the two-body forces of the Hamiltonian. In such
cases the single-particle space listed in the .dres file consists only of the
valence orbitals and the one-body density matrices are only specified for the
valence orbitals.

DMFormFactor reads the valence space orbitals from the .dres file and infers the
number of core nucleons by subtracting the number of valence protons and
neutrons from the number of nucleons in the target nucleus. The core orbitals
are assumed to be one of the standard shell model orbital sets associated with
possible cores: He-4, O-16, Ca-40, Ni-56, Sn-100.

The one-body density matrix elements for the core orbitals are then determined
from the (full) occupation of the core orbitals. In proton-neutron format:

\begin{equation}
\rho_{J,x=p,n}^{\Psi}(a,b)_{(core)} = \delta_{a,b}[j_a][J],
\end{equation}

where $[y] \equiv \sqrt{2y+1}$ and $j_a$ is the angular momentum of $a$-orbit.
$J$ is the total spin of the nuclear target state $\Psi$. And in isospin format
for a target state with total isospin $T$:
\begin{align}
    \rho_{J,\tau=0}^{\Psi}(a,b)_{(core)} &= \delta_{a,b} [1/2][j_a][J][T],\\
    \rho_{J,\tau=1}^{\Psi}(a,b)_{(core)} &= 0.0.
\end{align}

## Nuclear density matrix format
We adopt the output format from the BIGSTICK shell-model code. The output
one-body densities are written to a file with extension `.dres`. We provide a
full specification of this plain-text-file format in the `docs` directory. Here,
we show the form of the file and explain its contents.  
```
  State      E        Ex         J       T
    1   -330.17116   0.00000     1.500  11.500
  Single particle state quantum numbers
ORBIT      N     L   2 x J
     1     0     2     3
     2     0     2     5
     3     1     0     1
 Initial state #    1 E = -330.17117 2xJ, 2xT =    3  23
 Final state   #    1 E = -330.17117 2xJ, 2xT =    3  23
 Jt =   0, proton      neutron
    1    1   1.55844   5.40558
```
The file is comprised of three sections: 

1. Many-body state information
2. Single-particle state quantum numbers
3. Density matrix element blocks

Only the ground state is needed for inelastic WIMP-nucleus scattering
calculations.  The single-particle state quantum numbers specify the quantum
numbers for the simple-harmonic oscillator states involved in the one-body
operators. 

Finally, the one-body density matrix elements are listed in nested blocks with
three layers: 

* (i) the initial and final state specification (corresponding to the many-body states listed in section (1) of the file), 
* (ii) the angular momentum carried by the one-body density matrix operator, labeled {\tt Jt} here, 
* (iii) the single-particle state labels `a`, `b` in columns 1 and 2
  (corresponding to the single-particle state labels listed in section (2) of
  the file) and the proton and neutron (isospin-0 and isospin-1) density matrix
  elements in columns 3 and 4. 

Both (i) and (ii) must be specified along with columns 1 and 2 of (iii) in order
to fully determine a matrix element $\rho^{f,i}_K(a,b)$, where $K=J_t$. Note
that the values of $K$ are restricted by conservation of angular momentum; both
between the many-body states labeled $i$ and $f$, and the single-particle states
labeled $a$ and $b$.

