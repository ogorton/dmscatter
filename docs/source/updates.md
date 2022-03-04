# Update log

__Version 0.11 update (Feb. 28, 2022)__

* More robust module dependencies
* Legacy density matrices with conversion script

__Version 0.10 update (Jan. 12, 2022)__

* Validation test scripts
* Improved organization of code
* Comprehensive documentation

__Version 0.7 update (May 5, 2021)__

* Improved speed by re-writing Gaussian-Legendre routine in modern fortran with
  fixed (non-adaptive) number of evaluation points. The new implementation also
  supports parallelization, which the previously referenced library did not.

__Version 0.6 update (Apr. 19, 2021)__

* Improved speed by caching Wigner coefficients in memory
* Added options to compute transition probabilities and differential cross
  sections. (For fixed recoil energy, for a range of velocities.)
* More options in the Makefile
* Previous versions claimed compatibility with isospin-formalism density
  matrices. This turns out not to be the case. An appropriate error trap has
been added.

__Version 0.5 update (Feb. 22, 2021)__

* Renamed executable to dmf90factor.x (previously darkmatter.x)
* A ".sps" file is no longer required; the code now deduces this information
  using data provided by the (still required) ".dres" file.
* Manual now has quick-start guides for the Fortran and Python interfaces

__Version 0.4 update (Jan. 14, 2021)__

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

__Version 0.3 update (Jan. 13, 2021)__

New features:

* ~~Now supports nuclear density matrix files in either isospin or proton-neutron formalism~~ _see Version 0.6 notes._
* Inputs and outputs now carry specified units

Bugfixes:

* Fixed bug involving illegal sqrt() evaluations
* Updated numerical quadrature routine to library (instead of 'homebrew')
* Fixed error in denisty-matrix core-filler

__Version 0.2 update (Nov. 24, 2020)__

* Data reorganized to support future extension to multiple target species.

* Added python script which easily compares event rate spectra for different dark
  matter masses. To run the example, cd to sample/ and run:

        python ../python/masscompare.py

  Script is general and can be used to compare runs for any variable which can 
  be modified in the control file.

__Version 0.1 update (Nov. 20, 2020)__

* Now supports computing event rate spectra (event rate versus recoil energy in kev)
* Energy range is entered either (a) as a linear grid by specifying Emin, Emax, 
  Estep, or (b) from a file specifying energies
* Now takes advantage of multi-core systems using openMP

