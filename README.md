# A Fortran Program for Experimental WIMP Analysis
*Oliver Gorton, Changfeng Jiao, and Calvin Johnson*

Computes the differential event rate per recoil energy for WIMP-nucleon 
scattering events. This code based on the Mathematica package described 
in [this arXiv link](https://arxiv.org/abs/1308.6288).

Source code is located in src directory. To compile the program, you must have
the gfortran compiler. From the /src/ directory, run:

    make darkmatter.x
    
This places the darkmatter.x executable in the src directory. Move it to your bin
directory if desired.


An example is provided in the "sample" directory. To run this example, move to
the /sample/ directory and run:

    ../src/darkmatter.x < input.si28

This computes the event rate spectra for a Si28 target. 

Two example python script for running multiple jobs with different settings can
be found in the /python/ directory. For example, from within the /sample/
directory, try running:

    python3 ../python/masscompare.py

More extensive documentation can be found in the manual document.

### To-do:
* Update the manual

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
* Now supports nuclear density matrix files in either isospin or proton-neutron
  formalism
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
