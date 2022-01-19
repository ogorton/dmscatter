---
title: User Manual for DMFortFactor
subtitle: A Fast Fortran Code for WIMP-Nucleus Form Factors
author:
- Oliver C. Gorton
- Changfeng Jiao
- Calvin W. Johnson
geometry:
- margin=1in
toc: true
numbersections: true
header-includes: |
    \usepackage{amsmath}
    \usepackage{physics}
    \def\be{\begin{align}}
    \def\ee{\end{align}}
bibliography: source/refs.bib
csl: source/american-physics-society.csl
linkcolor: blue
urlcolor: blue
link-citations: yes
---

# Introduction
We present here  a fast modern Fortran code, DMFortFactor, for computing
WIMP-nucleus scattering event rates using a previously studied theoretical
framework[@PhysRevC.89.065501; @Fitzpatrick_2013], now with advanced algorithmic
and numerical implementation, including the ability to take advantage of
multi-core CPUs.  Furthermore, we enhance accessibility by including Python
wrappers with example scripts.

The code is available on the public 
[GitHub repository](https://github.com/ogorton/dmfortfactor).

This program is principally concerned with computing the dark matter-nucleus
differential event rate as a function of the nuclear recoil energy $E_R$:
$$
\frac{dR_D}{dE_R}
	= N_T\frac{\rho_\chi}{m_\chi}\int_{v_{min}}^{v_{escape}}
	\frac{2m_T}{4\pi v^2}\frac{1}{2j_\chi+1}\frac{1}{2j_T+1}
	\sum_{spins}|\mathcal{M}(v,q)|^2  \tilde{f}(\vec{v})vd^3v
$$
This quantity has units of events/GeV and is implicitly multiplied by
an effective exposure of 1 Kilogram-Day of target nuclei. This is done by
taking $N_t = 1\ kilogram\cdot day / m_T$, where $m_T$ is the mass of the target
nucleus in GeV. Recoil energies $E_R$ are given in keV.

## SuperQuickstart guide

- Navigate to the `build` directory from wherever you have stored `dmfortfactor/`
- Run the command: `make openmp`
- Navigate to the directory `runs` (e.g. `cd ../runs/`)
- Run the command: `python3 ../examples/exampleXe.py`

This should generate Figure 1.

![Example output graph.](source/exampleXe.pdf){width=80%}
