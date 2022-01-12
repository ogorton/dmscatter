# Python interface

We provide a Python interface (a wrapper) for the Fortran code and a number of
example scripts demonstrating its use. The wrapper comes with two Python
functions `EventrateSpectra`, and `NucFormFactor` which can be imported
from `dmfortfactor.py` in the Python directory. 

Each function has three required arguments: 

1. Z the number of protons in the target nucleus, 
2. N the number of neutrons in the target nucleus, and 
3. The filename for the one-body density matrix file describing the nuclear structure. If no other arguments are provided, default values will be used for all of the remaining necessary parameters, including zero interaction strength. 


## Event rate spectra
To calculate an event rate with a nonzero interaction, the user should also
provide one or more of the optional EFT coupling coefficient arrays: `cpvec,
cnvec, csvec, cvvec`. These set the couplings to protons, neutrons, isoscalar,
and isovector, respectively. The $0^{th}$ index sets the first operator
coefficient: `cpvec[0]`$= c_1^p$, etc.  Finally, the user can also pass a
dictionary of valid control keywords and values to the function in order to set
any of the control words defined in the manual.

To compute the event-rate spectra for $^{131}$Xe with a WIMP mass of 50 GeV and
a $c_3^v=0.0048$ coupling, one might call: 
```Python
import dmfortfactor as dm
control_dict = {"wimpmass" : 50.0}
cvvec = np.zeros(15)
cvvec[2] = 0.0048 
Erkev, ER = dm.EventrateSpectra(
            Z = 54,
            N = 77,
            dres = "../dres/xe131gcn",
            controlwords = control_dict,
            cvvec = cvvec,
            exec_path = "../build/dmfortfactor")
```
This will return the differential event rate spectra for recoil energies from 1
keV to 1 MeV in 1 keV steps. 

The file `xe131gcn.dres` must be accessible at the relative or absolute path
name specified (in this case `../dres/`), and contain a valid one-body density
matrix for $^{131}$Xe. Similarly, the DMFortFactor executable (`dmfortfactor`) 
should be accessible from the user's default path - or else the
path to the executable should be specified, as in the above example (`
exec_path = "../build/dmfortfactor"`).

## Nuclear response functions for WimPyDD

There are published codes which compute the WIMP-nucleus event rate spectra,
etc., but which rely on nuclear form factors from an external source. Once such
code is WimPyDD [@jeong2021wimpydd]. We have provided an additional option in
DMFortFactor which computes these nuclear form factors from the target
density-matrix and exports the results to a data file. 

The Python wrapper-function `NucFormFactor` runs the DMFortFactor option
to export the nuclear response functions to file, and additionally creates and
returns an interpolation function $W(q)$ which can be called. In the following
code listing, the nuclear response function for $^{131}$Xe is generated for
transfer momentum from 0.001 to 10.0 GeV/c.  
```Python
import dmfortfactor as dm
cwords = {
        "usemomentum": 1} # epmin/max/step sets momentum instead of energy
Wfunc = dm.NucFormFactor(
        Z = 54,
        N = 77,
        dres = "../dres/xe131gcn",
        controlwords = cwords,
        epmin = 0.001,
        epmax = 10.0,
        epstep = 0.001)
Wfunc(0.001)
```
The final line returns an (8,2,2)-shaped array with the evaluate nuclear
response functions at $q=0.001$ GeV/c. Note that, had we not set the keyword
`usemomentum` to 1, the function input values would have been specified in
terms of recoil energy (the default) instead of transfer momentum. 
