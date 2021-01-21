import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import runCustomInput

exec_name = "../src/darkmatter.x"
input_template = "input.xeAAA" # Filename ends in .template, this is the prefix.
workdir = os.getcwd()
plt.figure()

Z = 54
weights = [.26401, .04071, .21232, .26909, .10436, .08857]
for i,isotope in enumerate([129, 130, 131, 132, 134, 136]):
    neutrons = isotope - Z
    label = "xe"+str(isotope)
    param_dict = { "NEUTRONS" : neutrons, "AAA" : isotope }

    RecoilE, EventRate = runCustomInput.runCustomInput(exec_name, 
            input_template, param_dict, workdir, label)

    if (i==0):
        weightedsum = EventRate * weights[i]
    else:
        weightedsum += EventRate * weights[i]
    plt.plot(RecoilE, EventRate, label="$^{%s}$Xe"%isotope)

plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate/MeV")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("isotopecompare.pdf")
plt.figure()
plt.title("Differential event rate spectra by natural abundances m$_{chi}$=500GeV")
plt.plot(RecoilE, weightedsum)
plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate/MeV")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("weightedspectra.pdf")
plt.show()
