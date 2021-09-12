# A Fortran Program for Experimental WIMP Analysis
*Oliver Gorton, Changfeng Jiao, and Calvin Johnson*

Computes the differential event rate per recoil energy for WIMP-nucleon 
scattering events. This code based on the Mathematica package described 
in [this arXiv link](https://arxiv.org/abs/1308.6288).

Source code is located in src directory. To compile the program, you must have
the gfortran compiler. From the /src/ directory, run:

    make dmfortfactor
    
This places the dmfortfactor.x executable in the src directory. Move it to your bin
directory if desired.


An example is provided in the "sample/si" directory. To run this example, move to
the /sample/si directory and run:

    ../../src/dmfortfactor.x < input.si28

This computes the event rate spectra for a Si28 target. 

Two example python scripts for running multiple jobs with different settings can
be found in the /python/ directory. For example, from within the /sample/xe/
directory, try running:

    python3 ../../python/masscompare.py

More extensive documentation can be found in the manual document.

### Developer contacts
* Calvin Johnson (PI) cjohnson@sdsu.edu
* Oliver Gorton (Grad. student) ogorton@sdsu.edu
