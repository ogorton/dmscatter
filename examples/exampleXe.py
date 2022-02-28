import sys
sys.path.append("../python")
import dmfortfactor as dm
import numpy as np
import matplotlib.pyplot as plt

cwords = {"wimpmass" : 150.0}
cv = np.zeros(15)
cs = np.zeros(15)
cp = np.zeros(15)
cn = np.zeros(15)

cn[:] = 0.0048
cp[:] = 0.0048

Erkev, ER = dm.EventrateSpectra(
        Z = 54,
        N = 77,
        dres = "../dres/Xe/xe131gcn",
        controlwords = cwords,
        epmin = 1,
        epmax = 1000.0,
        epstep = 1.0,
        cs = cs,
        cv = cv,
        cp = cp,
        cn = cn,
        exec_path='../bin/dmfortfactor')
plt.plot(Erkev,ER, label="Xe131, GCN", linestyle=":")

Erkev, ER = dm.EventrateSpectra(
        Z = 54,
        N = 77,
        dres = "../dres/Xe/xe131jj55",
        controlwords = cwords,
        epmin = 1,
        epmax = 1000.0,
        epstep = 1.0,
        cs = cs,
        cv = cv,
        cp = cp,
        cn = cn,
        exec_path='../bin/dmfortfactor')
plt.plot(Erkev,ER, label="Xe131, jj55", linestyle="-")
Erkev, ER = dm.EventrateSpectra(
        Z = 54,
        N = 77,
        dres = "../dres/Legacy/sdgXe131",
        controlwords = cwords,
        epmin = 1,
        epmax = 1000.0,
        epstep = 1.0,
        cs = cs,
        cv = cv,
        cp = cp,
        cn = cn,
        exec_path='../bin/dmfortfactor')
plt.plot(Erkev,ER, label="Xe131, legacy", linestyle='--')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$E_{recoil}$ (keV)")
plt.ylabel(r"$dR/dE_r$ (1/MeV)")
plt.legend()
plt.savefig("exampleXe.pdf")
