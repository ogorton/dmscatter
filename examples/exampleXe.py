import sys
sys.path.append('/Users/oliver/projects/darkmatter/python')
import dmfortfactor as dm
import numpy as np

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
        dres = "xe131gcn",
        controlwords = cwords,
        ermin = 1,
        ermax = 1000.0,
        erstep = 1.0,
        csvec = cs,
        cvvec = cv,
        cpvec = cp,
        cnvec = cn)
