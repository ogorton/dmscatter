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

The key product of the code is the differential event rate for WIMP-nucleus scattering events in number of events per MeV. This is obtained by integrating the differential WIMP-nucleus cross section over the velocity distribution of the WIMP-halo in the galactic frame:
\begin{equation}\label{ER}
	\frac{dR}{dE_r}(E_r)
	 = N_T n_\chi \int \frac{d\sigma}{dE_r}(v,E_r)\ \tilde{f}(\vec{v})\ v\ d^3v,
\end{equation}
where $E_r$ is the recoil energy of the WIMP-nucleus scattering event, $N_T$ is the number of target nuclei, $n_\chi = \rho_\chi/m_\chi$ is the local dark matter number density, $\sigma$ is the WIMP-nucleus cross section.  The dark matter velocity distribution in the lab frame,
$\tilde{f}(\vec{v})$, is obtained by boosting the Galactic-frame distribution $f(\vec{v})$: $\tilde{f}(\vec{v}) = f(\vec{v} + \vec{v}_{earth})$, where $\vec{v}_{earth}$ is the velocity of the earth in the galactic rest frame. 

## SuperQuickstart guide

- Navigate to the `build` directory from wherever you have stored `dmfortfactor/`
- Run the command: `make openmp`
- Navigate to the directory `runs` (e.g. `cd ../runs/`)
- Run the command: `python3 ../examples/exampleXe.py`

This should generate Figure 1.

![Example output graph.](source/exampleXe.pdf){width=80%}
