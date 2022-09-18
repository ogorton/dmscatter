# dmscatter: A Fast Fortran Program for WIMP-Nucleus Scattering
*Oliver Gorton, Changfeng Jiao, and Calvin Johnson*

Computes the differential event rate per recoil energy for WIMP-nucleon
scattering events. This code is based on the Mathematica package,
[DMFormFactor](https://www.ocf.berkeley.edu/~nanand/software/dmformfactor/)
described in [this arXiv link](https://arxiv.org/abs/1308.6288).

Can be used as a standard Fortran program via the command line, or via our
convenient Python interfaces, which offer a minimal interface, optionally 
extensible to the full capabilities of the underlying program.

Read the documentation:
* [PDF](docs/dmscatter_User_Manual.pdf)

**Dev Contacts**

* Oliver Gorton, PhD student (email: ogorton@sdsu.edu)
* Calvin Johnson, PI (email: cjohnson@sdsu.edu)

## Compile
To compile the program, you must have the `gfortran` compiler and `make`
installed. to compile an OpenMP-parallel version of the code, navigate to the 
`build` directory and run:

    make openmp

To compile a serial version of the code, instead run:

    make dmscatter

Both options create the executable `dmscatter` in the `bin` directory. (Many
of the example Python scripts expect to find it there.) Note that if you want to
switch between and serial or parallel version, you must run `make clean`
in-between compiles.

## Run with Python
The main Python interface to the code can be imported and called:
```Python
import sys
sys.path.append("../python")
import dmscatter as dm
```
If necessary, replace "../python" with the path to the `dmscatter/python` diretory on your
system.

To run one of the example Python scripts, try moving to `examples` and running:

    python3 exampleXe.py

More extensive documentation can be found in the manual document.

### EventrateSpectra
The main function in this module, used to compute the differential event-rate
spectra. The interface looks like this:
```Python
EventrateSpectra(Z, N, dres=None, target=None, epmin=1, epmax=1000, epstep=1,
        controlwords={}, cp=None, cn=None, cs=None, cv=None,
        exec_path='dmscatter', name=None)
```
Here is an example call to this function:
```Python
import sys
sys.path.append("../python")
import dmscatter as dm
Recoilenergykev, Eventrate = dm.EventrateSpectra(
            Z = 54,
            N = 77,
            target = "../targets/Xe/xe131gcn",
            cn = [0.00048, 0,0,0,0,0,0,0,0,0,0,0,0,0,0] )
```
A more detailed example is given in `examples/exampleXe.py`.

### NucFormFactor
The other availale function is for computing nuclear form factors. The interface
is
```Python
NucFormFactor(Z, N, dres=None, target=None, epmin=1, epmax=1000, epstep=1,
    controlwords={}, exec_path='dmscatter', name=".nucFFspectra")
```
Here is an example call to this function:
```Python
Wfunc = dm.NucFormFactor(
        Z = 54,
        N = 77,
        dres = "../targets/Xe/xe131gcn",
        controlwords = cwords,
        epmin = 0.001,
        epmax = 10.0,
        epstep = 0.001,
        exec_path='../bin/dmscatter')
print(Wfunc(q=0.001))
```
A more detailed example is given in `examples/exampleXeFF.py`.

## Validation plots
We include a script for generating validation tests against data generated with
DMFormFactor (the Mathematica package). After compiling the code, navigate to
`test` and run

    make test

After a few minutes, and if you have [pandoc](https://pandoc.org/index.html) 
installed, a file called `validation.pdf` should be generated with a number of
plots comparing the output of dmscatter against DMFormFactor. If you don't
have pandoc, you can still view the plots individually in the test directory. 

## Directory

| Directory | Description |
| --------- | ----------- |
| bin       | This directory is created once you compile the code from the build directory. |
| make      | Compile the code with the included Makefile here. |
| docs      | Manual and other documentation | 
| data      | Nuclear structure information: the density matrix files (.dres) |
| examples  | Example Python scripts using our Python wrapper |
| python    | Python modules containing the wrapper for dmscatter |
| src       | Fortran source code and Makefile. Compile the code here. |
| test      | Test script and validation data. |
