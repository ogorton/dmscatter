import sys
sys.path.append("../python")
import dmscatter as dm
import numpy as np
import matplotlib.pyplot as plt

cwords = {
         "wimpmass" : 150.0,
         "ntscale" : 2500.0
         }
cv = np.zeros(15)
cs = np.zeros(15)
cp = np.zeros(15)
cn = np.zeros(15)

cn[:] = 4.8e-4
cp[:] = 4.8e-4

Erkev, ER = dm.EventrateSpectra(
        Z = 54,
        N = 77,
        dres = "../targets/Xe/xe131gcn",
        controlwords = cwords,
        epmin = 1,
        epmax = 1000.0,
        epstep = 1.0,
        cs = cs,
        cv = cv,
        cp = cp,
        cn = cn,
        exec_path='../bin/dmscatter')
plt.plot(Erkev,ER, label="Xe131, GCN", linestyle=":")

Erkev, ER = dm.EventrateSpectra(
        Z = 54,
        N = 77,
        dres = "../targets/Xe/xe131jj55",
        controlwords = cwords,
        epmin = 1,
        epmax = 1000.0,
        epstep = 1.0,
        cs = cs,
        cv = cv,
        cp = cp,
        cn = cn,
        exec_path='../bin/dmscatter')
plt.plot(Erkev,ER, label="Xe131, jj55", linestyle="-")

Erkev, ER = np.loadtxt("test.xe131.dat", unpack=True)
plt.plot(Erkev,ER, label="Xe131, DMFormFactor", linestyle='--')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$E_{recoil}$ (keV)")
plt.ylabel(r"$dR/dE_r$ (1/MeV)")
plt.legend()
plt.savefig("exampleXe.pdf")
