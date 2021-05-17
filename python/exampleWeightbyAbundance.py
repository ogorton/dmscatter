import os
import matplotlib.pyplot as plt
import dmfortfactor as dm

exec_name = "dmfortfactor.x"
input_template = "input.xeAAA" # Filename ends in .template, this is the prefix.
workdir = os.getcwd()
plt.figure()

Z = 54
isotopes = [128, 129, 130, 131, 132, 134, 136]
weights = [.01910, .26401, .04071, .21232, .26909, .10436, .08857]
weightedsum = 0.0

for i,isotope in enumerate(isotopes):
    N = isotope - Z
    label = "xe"+str(isotope)
    param_dict = { "NEUTRONS" : N, "AAA" : isotope }

    RecoilE, EventRate = dm.runCustomInput(exec_name, 
            input_template, param_dict, workdir, label)

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
