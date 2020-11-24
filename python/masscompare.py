import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

exec_name = "../src/darkmatter.x"

def runCustom(controltemplate, pnames, pvalues, workdir, label):

    path = os.getcwd()
    os.chdir(workdir)

    # Replace pnames in the control file with pvalues
    findandreplace = "sed '"
    for i,keyword in enumerate(pnames):
        findandreplace += "s/%s/%s/g;"%(keyword,pvalues[i])
    findandreplace += "' %s > %s.control"%(controltemplate,label)

    returned = subprocess.call(findandreplace, shell=True)

    #Run with given settings
    command = "%s < input.%s > output.%s"%(exec_name,label,label)
    print(command)
    returned = subprocess.call(command,shell=True)
    os.chdir(path)

    return


template = "si28.control.template"
pnames = ["MASS_KEYWORD"]
pvalues = [[]]
workdir = os.getcwd()
label = "si28"
resultfile = "eventrate_spectra.dat"

plt.figure()
for wimp_mass in (5., 50., 500.):
    pvalues = [wimp_mass]
    runCustom(template, pnames, pvalues, workdir, label)
    RecoilE, EventRate = np.loadtxt(resultfile, unpack=True, skiprows=1)
    plt.plot(RecoilE, EventRate, label="$m_\chi=$%2.2f"%wimp_mass)
plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
