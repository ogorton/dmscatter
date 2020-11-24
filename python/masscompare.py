import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import runner

exec_name = "../src/darkmatter.x"

# The following is the name of a control file which has the keyword pair:
#    wimpmass MASS_KEWORD
# This script works by replacing MASS_KEYWORD with the desired value of the
# WIMP mass, then running the code with a control file "si28.control".
template = "si28.control.template" 
pnames = ["MASS_KEYWORD"]
workdir = os.getcwd()

# This string determines the input file to be used (the command-line inputs,
# not the control file) for each run. The script will therefore run
#     <exec_name> < input.<label>
# 
label = "si28"

# Where to find the resulting event rate spectra:
resultfile = "eventrate_spectra.dat"

plt.figure()
for wimp_mass in (5., 50., 500.):
    pvalues = [wimp_mass]
    runner.runCustom(exec_name, template, pnames, pvalues, workdir, label)
    RecoilE, EventRate = np.loadtxt(resultfile, unpack=True, skiprows=1)
    plt.plot(RecoilE, EventRate, label="$m_\chi=$%2.2f"%wimp_mass)
plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
