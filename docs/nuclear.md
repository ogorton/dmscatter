# Nuclear structure input

Users must provide nuclear one-body density matrix elements of the form:
$$
    \rho_{K,T}^{\Psi}(a,b) = \langle \Psi| [\hat{c}_a^\dagger \hat{c}_b]_{K,T}|\Psi \rangle ,
$$
where $\Psi$ is the nuclear-target wave function and $\hat{c}^\dagger$,
$\hat{c}$ are the one-body creation, destruction operators. The matrix elements
must be stored in a file in a standard format produced by shell-model codes like
BIGSTICK.

Table of nuclear data we include with the program. Each corresponds to a (.dres)
density matrix file. The source indicates the nuclear Hamiltonian that was used
to generate the wave function data:

|Nuclei | Isotopes | Source |
| ----- | ------------ | ------ |
|Si     | 28, 29   | [@PhysRevC.74.034315] |
|Xe     | 128, 129, 130, 131, 132, 134, 136 | |
|Ar     | 40       | |
|C      | 12       |  [@cohen1965effective] |
|He     | 4        | |


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

$$
\rho_{J,x=p,n}^{\Psi}(a,b)_{(core)} = \delta_{a,b}[j_a][J],
$$

where $[y] \equiv \sqrt{2y+1}$ and $j_a$ is the angular momentum of $a$-orbit.
$J$ is the total spin of the nuclear target state $\Psi$. And in isospin format
for a target state with total isospin $T$:
\begin{align*}
    \rho_{J,\tau=0}^{\Psi}(a,b)_{(core)} &= \delta_{a,b} [1/2][j_a][J][T],\\
    \rho_{J,\tau=1}^{\Psi}(a,b)_{(core)} &= 0.0.
\end{align*}

