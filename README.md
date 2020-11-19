# A Fortran Program for Experimental WIMP Analysis

Computes the differential event rate for WIMP-nucleon scattering events.

Source code is located in src directory. To compile run from the src directory:

    make darkmatter.x
    
This places the darkmatter.x executable in the src directory. Move it to your bin
directory if desired.


An example is provided in the "sample" directory. cd to this directory
and run:

    ../src/darkmatter.x < input.si28
    
to run the sample input file.

More extensive documentation can be found in the manual document.


## Developer contacts
* Calvin Johnson (PI) cjohnson@sdsu.edu
* Oliver Gorton (Grad. student) ogorton@sdsu.edu
