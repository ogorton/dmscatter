# Usage Guide

DMFortFactor can be used interactively from the command line,
where the user is guided by prompts for a small number of datafiles and
parameters.  Naturally, this interactive process can be expedited by piping a
pre-written input file into the command line interface (CLI).

The command line interface (CLI) to the code prompts the user for
the type of calculation they wish to perform from a menu of options, then the
target nucleus, given by the number of protons $Z$ and neutrons $N$, two input
files, and finally other CLI inputs depending on the compute-option chosen.
These final inputs are typically three numbers specifying the range of recoil
energies or momentum for which to compute the output. After running, the code
prints the results to a plain text file in tabulated format.

The two files are:

1. Control file
2. Nuclear structure input

The control file specifies the EFT interaction and any optional settings
desired.

The nuclear structure inputs needed are one-body density matrices,  defined
below. We supply a library of density matrices for many of the common expected
targets, as listed [here](#nuclear-structure-input). The density matrix files
are written in plain ASCII, using the format output by the nuclear
configuration-interaction code BIGSTICK [@BIGSTICK1; @BIGSTICK2]. The only
assumption is that the single-particle basis states are harmonic oscillator 
states; the user must supply the harmonic oscillator single particle basis
frequency $\Omega$, typically given in MeV as $\hbar \Omega$, or the related
length parameter $b= \sqrt{\hbar/M\Omega}$, where $M$ is the nucleon mass. 
[^spb]

[^spb]: If density matrices are generated in some other single-particle basis, such as
    those from a Woods-Saxon potential or a Hartree-Fock calculation, that basis
    must be expanded into harmonic oscillator states.  By using harmonic oscillator
    basis states one can efficiently compute the matrix elements.  One can use
    either phenomenological or _ab initio_ model spaces and interactions; as
    an example, we provide density matrices for $^{12}$C, both from the
    phenomenological Cohen-Kurath shell model interaction [@cohen1965effective],
    and from an interaction derived from chiral effective field
    theory [@PhysRevC.68.041001].

While the DMFortFactor executable can be run by itself, we provide Python APIs
for integrating the Fortran program into Python work flows; see section .

We also provide a generic application programming interface (API) for the Python
language.  This API essentially offers a prescribed and easy-to-use way to run
DMFortFactor in a programmatic way.  Any sufficiently experienced linux user
could probably write a script to produce any possible set of inputs to our code.
But our API removes the need by making it easy for anyone who can use a Python
function to write their own advanced scripts allowing them to perform parameter
studies and comparisons of different inputs to the theory.

## Compiling with make

There are multiple directories in the project. All of the Fortran code which
needs to be compiled is found in the `src` directory.  

The easiest way to get started is simply to navigate to the `build` directory and run
```
make dmfortfactor
```
This will compile `DMFortFactor` using `gfortran`.  If you want to use
a different compiler, you must edit the following line in the Makefile:
```
FC = gfortran
```
changing `gfortran` to your compiler of choice.

If you want a OpenMP parallelized version of the code, you can compile with:
```
make openmp
```

Both of these options will compile the source code and leave the executable,
called `dmfortfactor` in the `bin` directory. (I.e.
`dmfortfactor/bin/dmfortfactor`.)
Note that if you change from a serial executable to a parallel executable (or
vice versa) you should run:
```
make clean
```

## Required files
There are two files required for any calculation:

1. Control file (.control)
2. Nuclear density matrix file (.dres)

Additionally, if the user enables the option `usenergyfile`, then a file
containing the input energies or momentum will also be required.

### Control file
Each EFT parameter is written on its own line in [mycontrolfile].control, with
four values: the keyword "coefnonrel", the operator number (integer 1..16), the
coupling type ("p"=proton, "n"=neutron, "s"=scalar, "v"=vecctor), and the
coefficient value. For example,
```
coefnonrel    1    s     3.1
```
would set $c_1^{\tau=0} = 3.1$. We take the isospin convention:
$$
		c^0 = c^p + c^n
$$
$$
		c^1 = c^p - c^n
$$
Thus, the previous example is equivalent to:
```
coefnonrel    1    p     1.55
coefnonrel    1    n     1.55
```

The control file also serves a more general but optional function: to set any
parameter in the program to a custom value.  
Simply add an entry to the control file with two values: the first
should be the keyword for the parameter and the
second should be the value to set that parameter to. For example, to set the
velocity of the earth in the galactic frame to $240\ km/s$, you should add the line:

```
vearth  240.0
```

As an example, here is the complete control file used to calculate the event
rate for the $c^n_1$ coupling to $^{131}$Xe:

```
# Coefficient matrix (non-relativistic)
# Ommitted values are assumed to be 0.0.
# c_i^t
# i = 1,...,16
# t: p=proton n=neutron s=scaler v=vector
coefnonrel  1  n  0.00048
wimpmass 150.0
vearth 232.0
maxwellv0 220.0
dmdens 0.3
usemomentum 0
useenergyfile 0
ntscale  2500.0
printdensities 0
#vescape 550.
```
Uncommenting the last line would set the escape velocity to 550 km/s.
A complete list of keywords is given [here](#control-file-keywords).

### Nuclear density matrix file (.dres)
We adopt the output format from the { BIGSTICK} shell-model code. The output
one-body densities are written to a file with extension { .dres}. We provide a
full specification of this plain-text-file format in the { docs} directory.
Here, we show the form of the file and explain its contents.

```
  State      E        Ex         J       T
    1   -330.17116   0.00000     1.500  11.500
  Single particle state quantum numbers
ORBIT      N     L   2 x J
     1     0     2     3
     2     0     2     5
     3     1     0     1
 Initial state #    1 E = -330.17117 2xJ, 2xT =    3  23
 Final state   #    1 E = -330.17117 2xJ, 2xT =    3  23
 Jt =   0, proton      neutron
    1    1   1.55844   5.40558
```

The file is comprised of three sections:

1. many-body state information
2. single-particle state quantum numbers
3. density matrix element blocks  

Only the ground state is needed for inelastic WIMP-nucleus scattering
calculations. The single-particle state quantum numbers specify the quantum
numbers for the simple-harmonic oscillator states involved in the one-body
operators.

Finally, the one-body density matrix elements are listed in nested blocks with
three layers: i. the initial and final state specification (corresponding to the
many-body states listed in section (1) of the file),
ii. the angular momentum carried by the one-body density matrix operator,
labeled { Jt} here, and iii. the single-particle state labels { a}, { b} in
columns 1 and 2 (corresponding to the single-particle state labels listed in
section (2) of the file) and the proton and neutron (isospin-0 and isospin-1)
density matrix elements in columns 3 and 4.

Both (i) and (ii) must be specified along with columns 1 and 2 of (iii) in order
to fully determine a matrix element $\rho^{f,i}_K(a,b)$, where $K=J_t$. Note
that the values of $K$ are restricted by conservation of angular momentum; both
between the many-body states labeled $i$ and $f$, and the single-particle states
labeled $a$ and $b$.

## Command-line interface
The program will prompt the user for the minimum necessary inputs to run a
calculation with default parameter values, including the name of a control file
which contains the EFT coefficients, and optionally, additional changes
to the calculation parameters.

After selecting the option {[er]} to compute an event rate spectra, there
are six further lines of input. These will be explained by an example:

```
 Enter the target proton number
54
 Enter the target neutron number
77
 Enter name of control file (.control):
xe131
...
  Enter name of one-body density file (.dres)
xe131gcn
...
 What is the range of recoil energies in kev?
 Enter starting energy, stoping energy, step size:
0.0001 250. 1.0
```

The first two entries are self-evident: we specifiy the number of protons and
neutrons in the target nucleus. In this case, 54 and 74, respectively, for
$^{131}$Xe.

Third is the name of the control file containing the EFT coefficients and other,
optional, settings. The `.control' file extension should be omitted. This
contents of this file will be explained in more detail later.

Fourth is the file containing the nuclear wave functions in the form of one-body
density matrices. Only the ground-state need be included. The '.dres' file
extension is omitted.

Fifth and finally are three numbers specifying the range of recoil energies
$E_R$ that the differential scattering rate should be computed for.

The event rate spectra will be written to a file, and as a side effect of the
calculation, the total event rate for the energy range requested will be
estimated by numerical integration. Note that the accuracy of this result will
depend on the choice of the step size.

## Compute options explained
There are a handful of options available from the main menu of the code:

* [er] Differential event rate, for a range of recoil energies or transfer momenta
* [cs] Differential cross section at a fixed recoil energy over a range of speeds
* [tp] Transition/scattering probability at a fixed recoil energy over a range of speeds
* [te] Total integrated events (without producing spectra file)
* [wd] Nuclear response functions at a given value of $y=(qb/2)^2$
* [ws] Nuclear response function spectra (for a range of $q$ or $E$)

The string in square brackets [x] is the compute-option.
Once a compute-option is chosen,  subsequent CLI inputs are the same up until
the density matrix file (.dres) has been read-in. Then, the inputs depend on the
compute-option chosen.  

- [cs] Differential cross section per recoil energy. Four additional inputs:
    1. E-recoil (keV)
	2. v-start (km/s)
	3. v-stop (km/s)
	4. v-step (km/s)
- [tp] Scattering probability. Same as [2].
- [te] Total scattering events per detector (does not produce spectra data).  This option uses adaptive quadrature to perform the integral of the event rate spectra with the fewest number of evaluations to reach the desired relative error.  This will be much faster than the result from options [1]. Three additional inputs:
	1. E-start (keV)
	2. E-stop (keV)
	3. Desired relative error (decimal value)
- [wd] Nuclear response function test. This compute-option allows the user to evaluate the nuclear response functions $W_i^{x,x'}(y)$ for a provided value of $y$. All combinations of $x$ and $x'$ will be printed for both isospin and proton-neutron couplings.  Two additional inputs are required:
	1. Function number (1 - 8)
	2. Value of $y=(qb/2)^2$ (dimensionless)
- [ws] Nuclear response function spectra. Enter a range of recoil energy or momentum values to evaluate the nuclear response funcions on. Tabulates the data to a file - one momentum per line or energy per line. Momentum/energy is written to the first column. The inputs are:
    1. E-start (keV)
    2. E-stop (keV)
    3. E-step (keV)

The following 32 columns store the (8,2,2)-dimensional response functions $W_i^{x,x'}$:
```
q_1  W_1^00 W_2^00 ... W_8^00 W_1^10 ... W_8^11
q_2  W_1^00 W_2^00 ... W_8^00 W_1^10 ... W_8^11
...
q_m  W_1^00 W_2^00 ... W_8^00 W_1^10 ... W_8^11
```

For the event-rate spectra [er] and for the nuclear response function
spectra [ws], the range of values is either over recoil energy $E_r$ (kev) (the
default) or over the transfer momentum $q$ (Gev/c). To use $q$ instead of $E_r$,
use the control word { usemomentum} set to 1.

