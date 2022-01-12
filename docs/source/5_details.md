# Details of computation

We present the equations necessary to reproduce the code. For a more complete
description of the theory, see
\href{https://link.aps.org/doi/10.1103/PhysRevC.89.065501}{Phys. Rev. C 89.065501.}

## Differential event rate
$$\label{ER}
	\frac{dR}{dE_r}(E_r)
	 = N_T n_\chi \int_{v_{min}}^{v_{escape}} \frac{d\sigma}{dE_r}(v,E_r) \tilde{f}(\vec{v})vd^3v,
$$
where $E_r$ is the recoil energy of the WIMP-nucleus scattering event, $N_T$ is
the number of target nuclei, $n_\chi = \rho_\chi/m_\chi$ is the local dark
matter number density, $\sigma$ is the WIMP-nucleus cross section.  The dark
matter velocity distribution in the lab frame, $\tilde{f}(\vec{v})$, is obtained
by boosting the Galactic-frame distribution $f(\vec{v})$: $\tilde{f}(\vec{v}) =
f(\vec{v} + \vec{v}_{earth})$, where $\vec{v}_{earth}$ is the velocity of the
earth in the galactic rest frame. The simplest model is a three-dimensional
Maxwell distribution:
$$
	f(\vec{v}) \propto e^{-\vec{v}^2/v_0^2},
$$
where $v_0$ is some scaling factor (typically taken to be around $220\ km/s$).

In order to evaluate the integral over the dark matter distribution, we make the
conversion to spherical coordinates. We need to evaluate an integral of the
form:
$$
I = \int_{v_{min}}^{v_{max}} d^3v \frac{f(\vec{v} + \vec{v}_{earth})}{v} = \int_{v_{min}}^{v_{max}} d^3v \frac{1}{v} e^{-(\vec{v}+\vec{v}_{earth})^2/v_0^2}
$$
Noting that $(\vec{v}+\vec{v}_{earth})^2 = \vec{v}^2 + \vec{v}^2_{earth} + 2vv_{earth}\cos(\theta)$,
with $\|\vec{v}\|\equiv v$ and $\theta$ defining
the angle between the two vectors, it's convenient to make the substitution
$d^3v = d\phi d(\cos \theta) v^2 dv$:

\begin{align}
	I &=  \int_0^{2\pi} d\phi \int_{v_{min}}^{v_{max}} dv \int_{-1}^1 d(\cos \theta) e^{-2vv_{earth}\cos\theta/v_0^2} v^2 \frac{1}{v} e^{-(\vec{v}^2+\vec{v}^2_{earth})/v_0^2}\\
	&= 2\pi \int_{v_{min}}^{v_{max}} dv v e^{-(\vec{v}^2+\vec{v}^2_{earth})/v_0^2} \left(-\frac{v_0^2}{2vv_{earth}} e^{-2vv_{earth}\cos\theta/v_0^2}\right)_{-1}^1\\
	&= \frac{\pi v_0^2}{v_{earth} }\int_{v_{min}}^{v_{max}} dv e^{-(\vec{v}^2+\vec{v}^2_{earth})/v_0^2}
		\left(- e^{-2vv_{earth}/v_0^2} + e^{+2vv_{earth}/v_0^2}\right)\\
	&= \frac{\pi v_0^2}{v_{earth} }\int_{v_{min}}^{v_{max}} dv
		\left(- e^{(v+v_{earth})^2/v_0^2} + e^{(v-v_{earth})^2/v_0^2}\right)\\
	&= \frac{\pi v_0^2}{v_{earth} }\int_{v_{min}}^{v_{max}} dv
		\left( g(v-v_{earth}) - g(v+v_{earth}) \right)
\end{align}

where in the last equality, we have defined a one-dimensional Gaussian form
$$
g(v) \propto e^{-v^2/v_0^2}.
$$

The final expression for $I$ can be trivially generalized to other spherically
symmetric velocity-dependent forms of the differential cross section. What's
important is the reduction of the velocity-boosted $d^3v$ integral to a radial
integral which can be carried out with one-dimensional quadrature:
\begin{align}
\int_{v_{min}}^{v_{max}} d^3v \sigma(v) e^{-(\vec{v}+\vec{v}_{earth})^2/v_0^2} \\
	= \frac{\pi v_0^2}{v_{earth} }\int_{v_{min}}^{v_{max}} dv \sigma(v) v^2\left( g(v-v_{earth}) - g(v+v_{earth}) \right).
\end{align}
The Fortran code uses this equation to evaluate the event rate integral with
quadrature. Analytic solutions exist in the form of error functions; we use
quadrature since it makes easy to later modify the velocity distribution (as
long as it remains spherically symmetric). For example, adding a velocity
cut-off is as easy as changing the limit on the quadrature, with no need to
write a whole new subroutine.

## Differential cross section

The differential cross section for the target nucleus can be expressed in terms
of either the nuclear recoil energy $E_R$, or the momentum transfer $q$:
$$
\frac{d\sigma(v,E_R)}{dE_R} = 2m_T \frac{d\sigma(v)(v,\vec{q}^2)}{d\vec{q}^2} = 2m_T\frac{1}{4\pi v^2}T(v,q),
$$
Where $v$ is the velocity of the dark matter particles in the lab-frame, $q$
is the momentum transfer of the scattering event, $m_T$ is the mass of the
target nucleus, and $T(v,q)$ is the transition or scattering probability. We
can see here that the differential cross section has an explicit $1/v^2$
dependence, independent of any velocity dependence of $T(v,q)$.


## Transition probability
The scattering probability is
$$
T(v,q) = \frac{1}{2j_\chi+1}\frac{1}{2j_T+1}\sum_{spins}|\mathcal{M}(v,q)|^2
$$
where $j_\chi$ is the spin of the WIMP, $j_T$ is the spin angular momentum of
the target nucleus, and $\mathcal{M}$ Galilean invariant amplitude, which is
defined by
\begin{align}
	T(v,q) = \frac{4\pi}{2j_T+1}\frac{1}{(4m_\chi)^2}
		\sum_{x=p,n}\sum_{x'=p,n}^1\sum_{i=1}^8 R_i^{xx'}(v^2,q^2)
		W_i^{xx'}(q)
\end{align}
where $m_\chi$ is the mass of the dark matter particle and $x$ is an index used
to sum over isospin couplings. The coefficients $R_i^{x,x'}$ are dark matter
particle response functions, to be define in another section. The operators
$W_i^{xx'}(q)$ are nuclear response functions, which are sums over matrix
elements of nuclear operators constructed from Bessel spherical harmonics and
vector spherical harmonics.

## Dark matter response functions
There are 8 dark matter response functions which group 15 operator coefficients
$c_i^x$ according the pair of nuclear response functions which they multiply.

As a shorthand, $cl(j) \equiv 4j(j+1)/3$, and $v^{\perp 2}\equiv v^2 - (q/2\mu_t)^2$.
\begin{align}
R_{M}^{xx'}(v,q) &= \frac{1}{4}cl(j_\chi) [ v^{\perp 2}
        (c_5^{x}c_5^{x'}q^2 + c_8^{x}c_8^{x'}) + c_{11}^{x}c_{11}^{x'}q^2 ]\\&
        + (c_1^{x} + c_2^{x}v^{\perp 2} ) (c_1^{x'}
        + c_2^{x'}v^{\perp 2} ) \\
R_{\Sigma''}^{xx'}(v,q) &= \frac{1}{16}cl(j_\chi) [c_6^{x}c_6^{x'}q^4
    + (c_{13}^{x}c_{13}^{x'}q^2 + c_{12}^{x} c_{12}^{x'} ) v^{\perp 2} + 2c_4^xc_6^{x'}q^2 + c_4^xc_4^{x'}]
    + \frac{1}{4}c_{10}^xc_{10}^{x '}q^2\\
R_{\Sigma'}^{xx'}(v,q) &= \frac{1}{32} cl(j_\chi) [ 2c_{9}^{x}c_{9}^{x'}q^2
        + ( c_{15}^{x}c_{15}^{x'}q^4 + c_{14}^{x}c_{14}^{x'}q^2 \\&
        - 2c_{12}^{x}c_{15}^{x'} q^2 + c_{12}^{x}c_{12}^{x'}) v^{\perp 2}
        + 2c_{4}^{x}c_{4}^{x'} ]
        +\frac{1}{8}(c_{3}^{x}c_3^{x'}q^2 + c_{7}^{x}c_{7}^{x'})v^{\perp 2}\\
R_{\Phi''}^{xx'}(v,q) &= \frac{q^2}{16m_N^2}cl(j_\chi) (c_{12}^x - c_{15}^{x}q^2
        )(c_{12}^{x '}-c_{15}^{x '}q^2 )
    + \frac{q^4}{4m_N^2}c_3^x c_3^{x'} \\
R_{\tilde{\Phi}'}^{xx'}(v,q) &= \frac{q^2}{16m_N^2}cl(j_\chi)(
        c_{13}^xc_{13}^{x'}q^2 + c_{12}^x c_{12}^{x'})\\
        %was erroneously c13c12q^2\\
R_{\Delta}^{xx'}(v,q) &= \frac{q^2}{4m_N^2}cl(j_\chi) (c_{5}^{x}c_{5}^{x'}q^2
        + c_{8}^{x}c_{8}^{x'})
        + 2\frac{q^2}{m_N^2}c_{2}^{x}c_{2}^{x'}v^{\perp 2}\\
R_{\Delta \Sigma'}^{xx'}(v,q) &= \frac{q^2}{4m_N}cl(j_\chi)
        (c_{4}^{x}c_{5}^{x'} - c_{8}^{x}c_{9}^{x'})
        - \frac{q^2}{m_N} c_{2}^{x}c_{3}^{x'} v^{\perp 2}\\
R_{\Phi''M}^{xx'}(v,q) &= \frac{q^2}{4m_N}cl(j_\chi)c_{11}^{x}
        (c_{12}^{x'} - c_{15}^{x'} q^2)
        + \frac{q^2}{m_N}c_{3}^{x'}  (c_{1}^{x} + c_{2}^{x} v^{\perp 2})\\
\end{align}

## Cross terms
Previous work has focused on setting limits on a single operator coupling at a
time. But of course, multiple couplings may exist simultaneously, and in fact,
some nuclear response functions are only activated with specific pairs of EFT
coefficients.

To create a minimal list of inputs to validate all possible nonzero
couplings, we need to test each coefficient on its own ($i=1,...,15$), and
also test the following 9 unique combinations: (1,2), (1,3), (2,3), (4, 5), (5,6), (8,9),
(11,12), (11,15), (12,15).

Table of EFT coefficient interactions. Shows which coefficients
    multiply each coefficient in addition to itself.

| Coefficient | Couples to |
| ---------- | ---------- |
| 1           |   2, 3   |
| 2           |   1, 3   |
| 3           |   1, 2   |
| 4           |   5, 6   |
| 5           |   4      |
| 6           |   4      |
| 7           |          |
| 8           |   9      |
| 9           |   8      |
| 10          |          |
| 11          |   12, 15 |
| 12          |   11, 15 |
| 13          |          |
| 14          |          |
| 15          |   11, 12 |

## Operators
The code uses the EFT coefficients in explicit proton-neutron couplings, i.e.
the interaction is defined by:
$$
    \mathcal{H} = \sum_{x=p,n}\sum_{i=1,15} c^x_i \mathcal{O} _i
$$
and the 15 momentum-dependent operators are:

\begin{align}
    \mathcal{O} _1 &= 1_\chi 1_N\\
    \mathcal{O} _2 &= (v^\perp)^2\\
    \mathcal{O} _3 &= i\vec{S}_N \cdot \left(\frac{\vec{q}}{m_N}\times
        \vec{v}^\perp\right)\\
    \mathcal{O} _4 &= \vec{S}_\chi \cdot \vec{S}_N\\
    \mathcal{O} _5 &= i\vec{S}_\chi \cdot \left(\frac{\vec{q}}{m_N}\times
        \vec{v}^\perp\right)\\
    \mathcal{O} _6 &= \left(\vec{S}_\chi \cdot \frac{\vec{q}}{m_N} \right)
        \left(\vec{S}_N \cdot \frac{\vec{q}}{m_N} \right) \\
    \mathcal{O} _7 &= \vec{S}_N\cdot \vec{v}^\perp \\
    \mathcal{O} _8 &= \vec{S}_\chi\cdot \vec{v}^\perp \\
    \mathcal{O} _9 &= i\vec{S}_\chi \cdot \left(\vec{S}_N \times
        \frac{\vec{q}}{m_N}\right)\\
    \mathcal{O} _{10} &= i\vec{S}_N \cdot \frac{\vec{q}}{m_N}\\
    \mathcal{O} _{11} &= i\vec{S}_\chi \cdot \frac{\vec{q}}{m_N}\\
    \mathcal{O} _{12} &= \vec{S}_\chi \cdot \left( \vec{S}_N\times
        \vec{v}^\perp\right)\\
    \mathcal{O} _{13} &= i\left( \vec{S}_\chi \cdot \vec{v}^\perp \right)
        \left(\vec{S}_N\cdot \frac{\vec{q}}{m_N}\right )\\
    \mathcal{O} _{14} &= i\left( \vec{S}_\chi \cdot \frac{\vec{q}}{m_N} \right)
        \left(\vec{S}_N\cdot \vec{v}^\perp \right )\\
    \mathcal{O} _{15} &= -\left(\vec{S}_\chi \cdot \frac{\vec{q}}{m_N} \right )
        \left( \left( \vec{S}_N\times \vec{v}^\perp\right)\cdot
        \frac{\vec{q}}{m_N} \right)
\end{align}

## Nuclear response functions

The EFT physics has been grouped into eight WIMP response functions
$R_i^{x,x'}$, and eight nuclear response functions $W_i^{x,x'}$.
The first six nuclear response functions have
the following form:
$$
    W_{X}^{x,x'} = \sum_{J}\bra{\Psi} X^{x}_J \ket{\Psi}\bra{\Psi} X^{x'}_J \ket{\Psi},
$$
with $X$ selecting one of the six electroweak operators,
$$
    X_J=M_J, \Delta_J, \Sigma_J', \Sigma_J'', \tilde{\Phi}_J', \Phi_J'',
$$
and $\Psi$ being the nuclear wave function for the ground state of the target
nucleus.  The sum over operators spins $J$ is restricted to even or odd values
of $J$, depending on restrictions from conservation of parity and charge
conjugation parity (CP) symmetry.

Two additional response functions add interference-terms:
$$
W_{M\Phi''}^{x,x'} =
\sum_{J}\bra{\Psi} M_{J}^{x} \ket{\Psi}\bra{\Psi} \Phi_{J}^{''x'} \ket{\Psi},
$$
$$
W_{\Delta\Sigma'}^{x,x'} = 
\sum_{J}\bra{\Psi} \Sigma_{J}^{'x} \ket{\Psi}\bra{\Psi} \Delta_{J}^{x'} \ket{\Psi}.
$$
The indices $i$ in equation (\ref{eq:T2}) correspond to these operators as:
$i\to X$ for $i=1,..,6$, and $i=7 \to M\Phi''$, $i=8\to \Delta\Sigma'$.

DMFortFactor can print the nuclear form factors to a file over a range of either
transfer momenta or recoil energy.

## Nuclear (electroweak) operators
There are six parity-and-CP-conserving nuclear operators, $M_J, \Delta_J,
\Sigma_J', \Sigma_J'', \tilde{\Phi}_J', \Phi_J''$, describing the electro-weak
coupling of the WIMPs to the nucleon degrees of freedom.  These are constructed
from Bessel spherical and vector harmonics [@DONNELLY1979103]:
\begin{align}
    M_{JM}(q\vec{x})\equiv j_J(qx)Y_{JM}(\Omega_x)
\end{align}
\begin{align}
    \vec{M}_{JML}(q\vec{x}) \equiv j_L(qx) \vec{Y}_{JLM}(\Omega_x),
\end{align}
where, using unit vectors $\vec{e}_{\lambda = -1, 0, +1}$,
\begin{align}
    Y_{JLM}(\Omega_x) = \sum_{m\lambda} \bra{Lm1\lambda}\ket{(L1)JM_J} Y_{Lm}(\Omega_x)\vec{e}_\lambda.
\end{align}
The six multipole operators are defined as:
\begin{align}
    \label{oplist}
M_{JM}\ \ &\\
\Delta_{JM} \equiv& \vec{M}_{JJM}\cdot \frac{1}{q}\vec{\nabla}\\
\Sigma'_{JM} \equiv& -i \left \{\frac{1}{q}\vec{\nabla}\times \vec{M}_{JJM}  \right\}\cdot \vec{\sigma}\\
\Sigma''_{JM} \equiv& \left \{ \frac{1}{q}\vec{\nabla}M_{JM} \right \}\cdot \vec{\sigma}\\
\tilde{\Phi}'_{JM} \equiv& \left( \frac{1}{q} \vec{\nabla} \times \vec{M}_{JJM}\right)\cdot \left(\vec{\sigma}\times \frac{1}{q}\vec{\nabla} \right) + \frac{1}{2}\vec{M}_{JJM}\cdot \vec{\sigma}\\
\Phi''_{JM}\equiv& i\left(\frac{1}{q}\vec{\nabla}M_{JM} \right)\cdot \left(\vec{\sigma}\times \frac{1}{q}\vec{\nabla} \right)
\end{align}

The matrix elements of these operators can be calculated for standard wave
functions from second-quantized shell model calculations:
$$
    \bra{\Psi_f} X_J \ket{\Psi_i} = \Tr(X_J \rho^{f,i}_J )
$$
$$
 = \sum_{a,b} \bra{a} |X_J| \ket{b} \rho^{fi}_J(ab),
$$
where single-particle orbital labels $a$ imply shell model quantum number $n_a,
l_a, j_a$, and the double-bar $||$ indicates reduced matrix
elements [@edmonds1996angular]. For elastic collisions, only the ground
state is involved, i.e. $\Psi_f=\Psi_i=\Psi_{g.s.}$.

We assume a harmonic oscillator single-particle basis, with the important
convention that the radial nodal quantum number $n_a$ starts at 0, that is, we
label the orbitals as $0s, 0p, 1s0d$, etc..,
and _not_ starting with $1s, 1p,$ etc.
Then, the one-body matrix elements for operators $\bra{a} |X^{(f)}_J| \ket{b}$,
built from spherical Bessel functions and vector spherical harmonics,  have
closed-form expressions in terms of confluent hypergeometric
functions [@DONNELLY1979103].

The nuclear structure input is in the form of one-body density matrices between
many-body eigenstates,
$$
\rho^{fi}_J(ab) = \frac{1}{\sqrt{2J+1} }\langle \Psi_f || [ \hat{c}^\dagger_a \otimes \tilde{c}_b ]_J
|| \Psi_i \rangle, \label{eqn:denmat}
$$
where $\hat{c}^\dagger_a$ is the fermion creation operator (with good angular
momentum quantum numbers), $\tilde{c}_b$ is the
time-reversed [@edmonds1996angular] fermion destruction operator.  Here the
matrix element is reduced in angular momentum but not isospin, and so are in
proton-neutron format. These density matrices are the product of a many-body
code, in our case BIGSTICK [@BIGSTICK1,@BIGSTICK2], although one could use
one-body density matrices, appropriately formatted, from any many-body code.


## Electroweak matrix elements
To compute the matrix elements of the electroweak operators in a harmonic
oscillator basis, we use the derivations from [@DONNELLY1979103]. Namely,
equations (1a) - (1f) and (3a) - (3d), which express the necessary geometric
matrix elements in terms of matrix elements of the spherical Bessel functions.
Here, we write out the remaining explicit formulas for obtaining matrix elements
of the Bessel functions $j_L(y)$ in a harmonic oscillator basis in terms of the
confluent hypergeometric function:
$$
    _1F_1(a,b,z) = \sum_{n=0}^\infty \frac{a^{(n)}z^n}{b^{(n)}n!},
$$
which makes use of the rising factorial function:
$$
    m^{(n)} = \frac{(m+n-1)!}{(m-1)!}.
$$

The first additional relation is:
\begin{align*}
\bra{n'l'j'} j_L(y) \ket{nlj} = \frac{2^L}{(2L+1)!!} y^{L/2} e^{-y}
    \sqrt{(n'-1)!(n-1)!}
    \sqrt{\Gamma(n'+l'+1/2)\Gamma(n+l+1/2)}\\
    \times
    \sum_{m=0}^{n-1}\sum_{m'=0}^{n'-1}
    \frac{(-1)^{m+m'}}{m!m'!(n-m-1)!(n'-m'-1)!} \\
    \times
    \frac{\Gamma[(l+l'+L+2m+2m'+3)/2]}{\Gamma(l+m+3/2)\Gamma(l'+m'+3/2)}
    \ _1F_1[(L-l'-l-2m'-2m)/2; L+3/2; y],
\end{align*}
which is computed in DMFortFactor by the function `BesselElement`.

The two additional relations are needed. As computed by `BesselElementMinus`:
\begin{align*}
\bra{n'l'j'} j_L(y) (\frac{d}{dy}-\frac{l}{y}) \ket{nlj}
    = \frac{2^(L-1)}{(2L+1)!!} y^{(L-1)/2} e^{-y}
    \sqrt{(n'-1)!(n-1)!}
    \sqrt{\Gamma(n'+l'+1/2)\Gamma(n+l+1/2)} \\
    \times
    \sum_{m=0}^{n-1}\sum_{m'=0}^{n'-1}
    \frac{(-1)^{m+m'}}{m!m'!(n-m-1)!(n'-m'-1)!}
    \frac{\Gamma[(l+l'+L+2m+2m'+2)/2]}{\Gamma(l+m+3/2)\Gamma(l'+m'+3/2)} \\
    \times
    \Big\{ -\frac{1}{2}(l+l'+L+2m+2m'+2)\ _1F_1[(L-l'-l-2m'-2m-1)/2; L+3/2; y]\\
    + 2m\ _1F_1[(L-l'-l-2m'-2m+1)/2; L+3/2; y] \Big\}.
\end{align*}
As computed by `BesselElementPlus`:
\begin{align*}
\bra{n'l'j'} j_L(y) (\frac{d}{dy}+\frac{l}{y}) \ket{nlj}
    = \frac{2^(L-1)}{(2L+1)!!} y^{(L-1)/2} e^{-y}
    \sqrt{(n'-1)!(n-1)!}
    \sqrt{\Gamma(n'+l'+1/2)\Gamma(n+l+1/2)} \\
    \times
    \sum_{m=0}^{n-1}\sum_{m'=0}^{n'-1}
    \frac{(-1)^{m+m'}}{m!m'!(n-m-1)!(n'-m'-1)!}
    \frac{\Gamma[(l+l'+L+2m+2m'+2)/2]}{\Gamma(l+m+3/2)\Gamma(l'+m'+3/2)} \\
    \times
    \Big\{ -\frac{1}{2}(l+l'+L+2m+2m'+2)\ _1F_1[(L-l'-l-2m'-2m-1)/2; L+3/2; y]\\
    + (2l+2m+1)\ _1F_1[(L-l'-l-2m'-2m+1)/2; L+3/2; y] \Big\}.
\end{align*}


## Wigner vector coupling functions
We implement a standard set of functions and subroutines for computing the
vector-coupling 3-j, 6-j, and 9-j symbols using the Racah alebraic expressions
[@edmonds1996angular].

One method we use to improve  compute time  is to cache Wigner 3-$j$ and 6-$j$
symbols~\cite{edmonds1996angular} (used to evaluate electro-weak matrix
elements) in memory at the start of run-time. As a side effect, our tests show
that this adds a constant compute time to any given calculation of roughly 0.3
seconds in serial execution and uses roughly 39 MB of memory (for the default
table size). As a point of comparison, the $^{131}$Xe example with all-nonzero
EFT coefficients has a run-time of 30 seconds in parallel execution. If we
disable the table caching, the run-time is roughly 150 seconds, 5 times longer.
The size of the table stored in memory can be controlled via the control file
with the keywords `sj2tablemin` and `sj2tablemax`.

For the 3-j symbol, we use the relation to the Clebsh-Gordon vector-coupling
coefficients:
\begin{align*}
    \begin{pmatrix}
        j_1 & j_2 & J\\
        m_1 & m_1 & M
    \end{pmatrix}
    = (-1)^{j_1-j_2-M}(2J+1)^{-1/2}\\
    (j_1j_2m_1m_2 | j_1 j_2; J, -M).
\end{align*}
The vector coupling coefficients are computed as:
\begin{align*}
    (j_1j_2 & m_1m_2 | j_1 j_2; J, M) = \delta(m_1+m_1,m) (2J+1)^{1/2}\Delta(j_1j_2J)\\
    & \times[(j_1+m_1)(j_1-m_1)(j_2+m_2)(j_2-m_2)(J+M)(J-M)]^{1/2}\sum_z (-1)^z \frac{1}{f(z)},
\end{align*}
where
\begin{align*}
    f(z) &= z!(j_1+j_2-J-z)!(j_1-m_2-z)!\\
    & \times(j_2+m_2-z)!(J-j_2+m_1+z)!(J-m_1-m_2+z)!,
\end{align*}
and
\begin{align*}
    \Delta(abc) = \left[\frac{(a+b-c)!(a-b+c)!(-a+b+c)!}{(a+b+c+1)!} \right]^{1/2}.
\end{align*}
The sum over $z$ is over all integers such that the factorials are well-defined
(non-negative-integer arguments).

Similarly, for the 6-j symbols:
\begin{align*}
    \begin{Bmatrix}
        j_1 & j_2 & j_3\\
        m_1 & m_1 & m_3
    \end{Bmatrix}
    &= \Delta(j_1j_2j_3)\Delta(j_1m_2m_3)\Delta(m_1j_2m_3)\\
    &\times \Delta(m_1m_2j_3) \sum_z (-1)^z\frac{(z+1)!}{g(z)},
\end{align*}
with
\begin{align*}
    g(z) &= (\alpha - z)!(\beta-z)!(\gamma-z)!\\
    &\times (z-\delta)!(z-\epsilon)!(z-\zeta)!(z-\eta)!
\end{align*}
\begin{align*}
    \alpha &= j_1+j_1+m_1+m_2 & \beta  &= j_2+j_3+m_2+m_3\\
    \gamma &= j_3+j_1+m_3+m_1 \\
    \delta &= j_1+j_2+j_3 & \epsilon &= j_1+m_2+m_3 \\
    \zeta &= m_1+j_2+m_3 & \eta &= m_1+m_2+j_3.
\end{align*}

For the 9-j symbol, we use the relation to the 6-j symbol:
\begin{align*}
        \begin{Bmatrix}
        j_1 & j_2 & j_3\\
        j_4 & j_5 & j_6\\
        j_7 & j_8 & j_9
    \end{Bmatrix}
    &= \sum_k (-1)^{2k} (2k+1) \\
        &\times \begin{Bmatrix}
        j_1 & j_4 & j_7\\
        j_8 & j_9 & z
        \end{Bmatrix}
        \begin{Bmatrix}
        j_2 & j_5 & j_8\\
        j_4 & z & j_6
        \end{Bmatrix}
        \begin{Bmatrix}
        j_3 & j_6 & j_9\\
        z & j_1 & j_2
        \end{Bmatrix}.        
\end{align*}
The 6-j symbols used to calculate the 9-j symbol are first taken from any
tabulated values. Otherwise, they are computed as previously described.
