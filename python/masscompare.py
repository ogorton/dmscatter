import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import runCustomControl

exec_name = "../src/darkmatter.x"

# This script works by replacing MASS_KEYWORD with the desired value of the
# WIMP mass, then running the code with a control file "si28.control".

control_template = "si28.control" # Filename ends in .template, this is the prefix.
inputfile = 'input.si28'
workdir = os.getcwd()
label = "si28"
# Where to find the resulting event rate spectra:
resultfile = "eventrate_spectra.dat"

plt.figure()

for wimp_mass in (5., 50., 500.):

    param_dict = { "MASS_KEYWORD" : wimp_mass }
    runCustomControl.runCustomControl(exec_name, inputfile, control_template, param_dict,
        workdir, label)

    RecoilE, EventRate = np.loadtxt(resultfile, unpack=True, skiprows=1)
    plt.plot(RecoilE, EventRate, label="$m_\chi=$%2.2f"%wimp_mass)

plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("%s.WIMPmassCompare.pdf"%label)
plt.show()
