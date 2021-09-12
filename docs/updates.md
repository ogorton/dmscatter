# A Fortran Program for Experimental WIMP Analysis
*Oliver Gorton, Changfeng Jiao, and Calvin Johnson*

### Version 1.7 update (May 5, 2021)
* Improved speed by re-writing Gaussian-Legendre routine in modern fortran with
  fixed (non-adaptive) number of evaluation points. The new implementation also
  supports parallelization, which the previously referenced library did not.
### Version 1.6 update (Apr. 19, 2021)
* Improved speed by caching Wigner coefficients in memory
* Added options to compute transition probabilities and differential cross
  sections. (For fixed recoil energy, for a range of velocities.)
* More options in the Makefile
* Previous versions claimed compatibility with isospin-formalism density
  matrices. This turns out not to be the case. An appropriate error trap has
been added.
### Version 1.5 update (Feb. 22, 2021)
* Renamed executable to dmf90factor.x (previously darkmatter.x)
* A ".sps" file is no longer required; the code now deduces this information
  using data provided by the (still required) ".dres" file.
* Manual now has quick-start guides for the Fortran and Python interfaces

### Version 1.4 update (Jan. 14, 2021)
New features:
* Compute integrated event rate spectra (total events) using adaptive
  integration routine
* Computing an event rate spectra will also report the total integrated event
  rate
* EFT coefficients can now be provided as either proton/neutron couplings, or as
  scaler/vector isospin couplings

Bugfixes:
* Updated definition of proton-neutron to isospin transformation to be
  consistent with Mathematica script definition (script, not paper)
### Version 1.3 update (Jan. 13, 2021)
New features:
* ~~Now supports nuclear density matrix files in either isospin or proton-neutron
  formalism~~ _see Version 1.6 notes._
* Inputs and outputs now carry specified units

Bugfixes:
* Fixed bug involving illegal sqrt() evaluations
* Updated numerical quadrature routine to library (instead of 'homebrew')
* Fixed error in denisty-matrix core-filler

### Version 1.2 update (Nov. 24, 2020)
* Data reorganized to support future extension to multiple target species.

* Added python script which easily compares event rate spectra for different dark
  matter masses. To run the example, cd to sample/ and run:

        python ../python/masscompare.py

  Script is general and can be used to compare runs for any variable which can 
  be modified in the control file.

### Version 1.1 update (Nov. 20, 2020)
* Now supports computing event rate spectra (event rate versus recoil energy in kev)
* Energy range is entered either (a) as a linear grid by specifying Emin, Emax, 
  Estep, or (b) from a file specifying energies
* Now takes advantage of multi-core systems using openMP


### Developer contacts
* Calvin Johnson (PI) cjohnson@sdsu.edu
* Oliver Gorton (Grad. student) ogorton@sdsu.edu
