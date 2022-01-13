# DMFortFactor: A Fast Fortran Program for WIMP-Nucleus Form Factors
*Oliver Gorton, Changfeng Jiao, and Calvin Johnson*

Computes the differential event rate per recoil energy for WIMP-nucleon
scattering events. This code is based on the Mathematica package,
[DMFormFactor](https://www.ocf.berkeley.edu/~nanand/software/dmformfactor/)
described in [this arXiv link](https://arxiv.org/abs/1308.6288).

Can be used as a standard Fortran program via the command line, or via our
convenient Python interfaces, which offer a minimal interface, optionally 
extensible to the full capabilities of the underlying program.

**Dev Contacts**

* Oliver Gorton (PhD student) ogorton@sdsu.edu
* Calvin Johnson (PI) cjohnson@sdsu.edu

## Compile
Source code is located in `src` directory. To compile the program, you must have
the `gfortran` compiler and `make` installed. From the `src` directory, run:

    make dmfortfactor

to compile a serial version of the code. To compile an OpenMP-parallel version
of the code, instead run:

    make openmp

Both options create the executable `dmfortfactor` in the `src` directory. (Many
of the example Python scripts expect to find it there.) Note that if you want to
switch between and serial or parallel version, you must run `make clean`
in-between compiles.

## Run with Python
The main Python interface to the code can be imported and called:
```Python
import sys
sys.path.append("../python")
import dmfortfactor as dm
```
If necessary, replace "../python" with the path to the `dmfortfactor/python` diretory on your
system.

The main function in this module, used to compute the differential event-rate
spectra, can be called like this:
```Python
Recoilenergykev, Eventrate = dm.EventrateSpectra(
            Z = 54,
            N = 77,
            dres = "../dres/xe131gcn",
            cn = [0.00048, 0,0,0,0,0,0,0,0,0,0,0,0,0,0] )
```
To run one of the example Python scripts, try moving to `runs` and running:

    python3 ../examples/exampleXe.py

More extensive documentation can be found in the manual document.

## Validation plots
We include a script for generating validation tests against data generated with
DMFormFactor (the Mathematica package). After compiling the code, navigate to
`test` and run

    make test

After a few minutes, and if you have [pandoc](https://pandoc.org/index.html) 
installed, a file called `validation.pdf` should be generated with a number of
plots comparing the output of DMFortFactor against DMFormFactor. If you don't
have pandoc, you can still view the plots individually in the test directory. 

## Directory

| Directory | Description |
| --------- | ----------- |
| docs      | Manual and other documentation | 
| dres      | Nuclear structure information: the density matrix files (.dres) |
| examples  | Example Python scripts using our Python wrapper |
| python    | Python modules containing the wrapper for DMFortFactor |
| runs      | A good place to run the code and example scripts. (Some of the examples expect to be run here. |
| src       | Fortran source code and Makefile. Compile the code here. |
| test      | Test script and validation data. |
