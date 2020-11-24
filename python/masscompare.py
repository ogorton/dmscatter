import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import runner

exec_name = "../src/darkmatter.x"
template = "si28.control.template"
pnames = ["MASS_KEYWORD"]
workdir = os.getcwd()
label = "si28"
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
