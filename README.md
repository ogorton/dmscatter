# A Fortran Program for Experimental WIMP Analysis
*Oliver Gorton, Changfeng Jiao, and Calvin Johnson*

**Dev Contacts**

* Oliver Gorton (PhD student) ogorton@sdsu.edu
* Calvin Johnson (PI) cjohnson@sdsu.edu

Computes the differential event rate per recoil energy for WIMP-nucleon
scattering events. This code based on the Mathematica package described in [this
arXiv link](https://arxiv.org/abs/1308.6288).

Source code is located in `src` directory. To compile the program, you must have
the gfortran compiler and make installed. From the `src` directory, run:

    make dmfortfactor

to compile a serial version of the code. To compile an OpenMP-parallel version
of the code, instead run:

    make openmp

This places the dmfortfactor executable in the `src` directory. Move it to
your path's bin directory if desired. Note that if you want to switch between
and serial or parallel version, you must run `make clean' to clear out the old
object files.

To run one of the example Python scripts, try moving to `runs` and running:

    python3 ../examples/exampleXe.py

More extensive documentation can be found in the manual document.

## Validation plots
We include a script for generating validation tests against data generated with
DMFormFactor, the Mathematica package. After compiling the code, navigate to
`test` and run

    make test

After a few minutes a file called `validation.pdf` should be generated with a
number of plots comparing the output of DMFortFactor against DMFormFactor.
