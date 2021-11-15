# A Fortran Program for Experimental WIMP Analysis
*Oliver Gorton, Changfeng Jiao, and Calvin Johnson*

Computes the differential event rate per recoil energy for WIMP-nucleon
scattering events. This code based on the Mathematica package described in [this
arXiv link](https://arxiv.org/abs/1308.6288).

Source code is located in `src` directory. To compile the program, you must have
the gfortran compiler and Cmake installed. From the `build` directory, run:

    cmake ../src
    cmake --build .

This places the dmfortfactor.x executable in the `build` directory. Move it to
your bin directory if desired.

To run one of the example Python scripts, try moving to `runs` and running:

    python3 ../examples/exampleXe.py

More extensive documentation can be found in the manual document.

### Developer contacts
* Calvin Johnson (PI) cjohnson@sdsu.edu
* Oliver Gorton (Grad. student) ogorton@sdsu.edu
