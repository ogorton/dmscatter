import sys
sys.path.append('/Users/oliver/projects/darkmatter/python')
import matplotlib.pyplot as plt
import dmfortfactor as dm
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
        dres = "../dres/xe131gcn",
        controlwords = cwords,
        ermin = 1,
        ermax = 1000.0,
        erstep = 1.0,
        csvec = cs,
        cvvec = cv,
        cpvec = cp,
        cnvec = cn)

    plt.plot(RecoilE, EventRate, label="$m_\chi=$%2.2f"%wimp_mass)

plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("%s.WIMPmassCompare.pdf"%label)
plt.show()
