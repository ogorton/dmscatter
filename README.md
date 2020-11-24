# A Fortran Program for Experimental WIMP Analysis
*Oliver Gorton, Changfeng Jiao, and Calvin Johnson*

Computes the differential event rate for WIMP-nucleon scattering events. This code 
is based on the Mathematica package described in [this arXiv link](https://arxiv.org/abs/1308.6288).

Source code is located in src directory. To compile run from the src directory:

    make darkmatter.x
    
This places the darkmatter.x executable in the src directory. Move it to your bin
directory if desired.


An example is provided in the "sample" directory. cd to this directory
and run:

    ../src/darkmatter.x < input.si28
    
to run the sample input file.

More extensive documentation can be found in the manual document.

### Version 1.2 update (Nov. 24, 2020)
* Data reorganized to support future extension to multiple target species.
* Add python script which easily compares event rate spectra for different dark
  matter masses. Script is general and can be used to compare runs for any
  variable which can be modified in the control file.

### Version 1.1 update (Nov. 20, 2020)
* Now supports computing event rate spectra (event rate versus recoil energy in kev). 
* Energy range is entered either (a) as a linear grid by specifying Emin, Emax, Estep, or (b) from a file specifying energies
* Now takes advantage of multi-core systems using openMP.


### Developer contacts
* Calvin Johnson (PI) cjohnson@sdsu.edu
* Oliver Gorton (Grad. student) ogorton@sdsu.edu
