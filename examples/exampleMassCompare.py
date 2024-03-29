import sys
sys.path.append("../python")
import matplotlib.pyplot as plt
import dmscatter as dm
import numpy as np

label = "xe131"

plt.figure()
cv = np.zeros(15)
cs = np.zeros(15)
cp = np.zeros(15)
cn = np.zeros(15)
cn[0] = 0.0048

for wimp_mass in (50.,150.0, 500.,5000.):

    cwords = { "wimpmass" : wimp_mass }

    RecoilE, EventRate = dm.EventrateSpectra(
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

    plt.plot(RecoilE, EventRate, label="$m_\chi=$%2.2f"%wimp_mass)

plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("%s.WIMPmassCompare.pdf"%label)
plt.show()
