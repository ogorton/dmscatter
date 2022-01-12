# Control file keywords

| Keyword | Symbol | Meaning | Units | Default |
| --------- | --- | --------------------------- | ---- | ----- |
| dmdens  | $\rho_\chi$ | Local dark matter density. | GeV/cm$^3$ | 0.3 |
| dmspin  | $j_\chi$    | Instrinsic spin of WIMP particles. | $\hbar$ | $\frac{1}{2}$ |
| fillnuclearcore | | Logical flag (enter 0 for False, 1 for True) to fill the inert-core single-particle orbitals in the nuclear level densities. Phenomenological shell model calculations typically provide only the density matrices for the active valence-space orbitals, leaving it to the user to infer the core-orbital densities. This option automatically assigns these empty matrix elements assuming a totally filled core. | | 1 (true) |
| gaussorder | | Order of the Gauss-Legendre quadrature to use when using Type 2 quadrature. (See quadtype.) An n-th order routine will perform n function evaluations.  Naturally, a higher order will result in higher precision, but longer compute time. | | 12 |
| hofrequency | $\hbar \omega$ | Set the harmonic oscillator length by specifying the harmonic oscillator frequency. (b = 6.43/sqrt($\hbar\omega$)). If using an \textit{ab initio} interaction, $\hbar \omega$ should be set to match the value used in the interaction. | MeV | See hoparameter. |
| hoparameter | $b$ | Harmonic oscillator length. Determines the scale of the nuclear wavefunction interaction. | fm | See eqn. (\ref{bho}). |
|  maxwellv0 | $v_0$ | Maxwell-Boltzman velocity distribution scaling factor. | km/s | 220.0 |
|  mnucleon | $m_N$ | Mass of a nucleon. It's assumed that $m_p\approx m_n$. |
  GeV | 0.938272 |
|  ntscale | $N_t$ | Effective number of target nuclei scaling factor. The differential event rate is multiplied by this constant in units of kilogram-days. For example, if the detector had a total effective exposure of 2500 kg days, one might enter 2500 for this value. | kg days | 1.0 |
|  quadrelerr |  | Desired relative error for the adaptive numerical quadrature routine (quadtype 1).  | | $10^{-6}$ |
|  quadtype | | Option for type of numerical quadrature. (Type 1 = adaptive 8th order Gauss-Legendre quadrature.  Type 2 = static n-th order Gauss-Legendre quadrature.) || 1 (type 1) |
|  sj2tablemax | | Maximum value of $2\times J$ used when caching Wigner 3-J and 6-J functions into memory. | | 12 |
|  sj2tablemin | | Minimum value of $2\times J$ used when caching Wigner 3-J and 6-J functions into memory. | | -2 |
|  useenergyfile | | Logical flag (enter 0 for False, 1 for True) to read energy grid used for calculation from a user-provided file intead of specifying a range. | | 0 (false) |
|  usemomentum | | Logical flag (enter 0 for False, 1 for True) to use momentum transfer intead of recoil energy as the independent variable. | |0 (false) |
|  vearth | $v_{earth}$ | Speed of the earth in the galactic frame. | km/s | 232.0 |
|  vescape | $v_{escape}$ | Galactic escape velocity. Particles moving faster than this speed will escape the galaxy, thus setting an upper limit on the WIMP velocity distribution. | km/s | 12 $\times\ v_{scale}$ |
|  weakmscale | $m_v$ | Weak interaction mass scale. User defined EFT coefficients are divided by $m_v^2$. | GeV | 246.2 |
|  wimpmass | $m_\chi$ | WIMP particle mass. | GeV | 50.0 |

# References

::: {#refs}
:::
