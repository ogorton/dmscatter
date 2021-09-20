import sys
sys.path.append("/Users/oliver/projects/darkmatter/python")
import numpy as np
import matplotlib.pyplot as plt
import dmfortfactor as dm

plt.figure()
Z = 54
isotopes = [128, 129, 130, 131, 132, 134, 136]
weights = [.01910, .26401, .04071, .21232, .26909, .10436, .08857]
weightedsum = 0.0

for i,isotope in enumerate(isotopes):
    N = isotope - Z
    controls = {
            "wimpmass" :500.0,
            "vescape" : 550.0}
    cp = np.zeros(15)
    cp[0] = 1
    RecoilE, EventRate = dm.EventrateSpectra(
            Z,
            N,
            dres = "xe%igcn"%isotope,
            controlwords = controls,
            cpvec = cp )

    weightedsum += EventRate * weights[i]

    label = "xe"+str(isotope)
    plt.plot(RecoilE, EventRate, label="$^{%s}$Xe"%isotope)

plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate/MeV")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("isotopecompare.pdf")
plt.figure()
plt.title(r"Differential event rate spectra by natural abundances m$_{chi}$=50GeV")
plt.plot(RecoilE, weightedsum)
plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate/MeV")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("weightedspectra.pdf")
plt.show()
