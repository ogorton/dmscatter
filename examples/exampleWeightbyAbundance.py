import sys
sys.path.append("../python")
import numpy as np
import matplotlib.pyplot as plt
import dmfortfactor as dm

plt.figure()
Z = 54
isotopes = [128, 129, 130, 131, 132, 134, 136]
weights = [.01910, .26401, .04071, .21232, .26909, .10436, .08857]
weightedsum = 0.0
mchi = 1000.
vesc = 550.
paths = ("../targets/Xe/xe","../targets/Xe/xe","../targets/Legacy/sdgXe")
ints =  ("jj55","gcn","")

plt.title(r"Differential event rate spectra by natural abundances"+
        " m$_{\chi}$=%sGeV"%mchi+"\n"+r"$v_{esc} = %s$km/s, $O_8$"%vesc+
        "$(c_s=1)$")
for j in range(len(ints)):

    weightedsum = 0.0

    for i,isotope in enumerate(isotopes):
        N = isotope - Z
        controls = {
                "wimpmass" :mchi,
                "vescape" : vesc}
        cs = np.zeros(15)
        cs[0] = 1
        dres = paths[j]+"%s"%isotope+ints[j]
        RecoilE, EventRate = dm.EventrateSpectra(
                Z,
                N,
                dres = dres,
                controlwords = controls,
                cs = cs,
                exec_path='../bin/dmfortfactor')
    
        weightedsum += EventRate * weights[i]
    plt.plot(RecoilE, weightedsum, label=paths[j]+ints[j])
    
plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate/MeV")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("weightedspectra.pdf")
