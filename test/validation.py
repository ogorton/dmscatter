import sys
sys.path.append('../python')
import dmscatter as dm

from scipy.interpolate import interp1d
import shutil
import numpy as np
import matplotlib.pyplot as plt
import dmscatter as dm

with open("operators.txt") as f: operatorsymbols = f.readlines()

# These are the pairs of linearly independent operator coefficients needed to
# completely test all operator terms.
pairs = [
        (1,1), (1,2), (1,3), 
        (2,2), (2,3),
        (3,3),
        (4,4), (4,5),
        (5,5), (5,6),
        (6,6), 
        (7,7),
        (8,8), (8,9),
        (9,9),
        (10,10),
        (11,11), (11,12), (11,15),
        (12,12), (12,15),
        (13,13),
        (14,14),
        (15,15)]
vesc = 550.0
cwords = {
        "wimpmass" : 150.0,
        "ntscale"  : 2500.0,
        "vearth" : 232.0,
        "maxwellv0" : 220.0,
        "vescape" : vesc,
        "quadrelerr" : 1e-6,
        "printdensities": 1
        }

log = True
typ = "lin"
if log: typ="log"

dreslbl = ("5 decimal .dres","7 decimal .dres")

for coupleto in ("n","p"):

    for opi, pair in enumerate(pairs):

        plt.figure()
        o1, o2 = pair[:]
        e_cut, er_cut = np.loadtxt("data/dmformfactor.vesc550.si29.%s.%i.%i.dat"%(
                coupleto, o1, o2), unpack=True)

        e_mb, er_mb = np.loadtxt("data/dmformfactor.si29.%s.%i.%i.dat"%(
                coupleto, o1, o2), unpack=True)

        for idres, dres in enumerate(("si29usd-iso","si29nb-iso")):

            print("%s-coupling, operator %i-%i, density %s"%(coupleto,o1,o2,dres))
            cv = np.zeros(15)
            cs = np.zeros(15)
            cp = np.zeros(15)
            cn = np.zeros(15)
            if coupleto=='n':
                cn[o1-1] = 1.0
                cn[o2-1] = 1.0
            else:
                cp[o1-1] = 1.0
                cp[o2-1] = 1.0

            energy, eventrate = dm.EventrateSpectra(
                Z = 14,
                N = 15,
                dres = "%s"%dres,
                controlwords = cwords,
                epmin = 1,
                epmax = 1000.0,
                epstep = 1.0,
                exec_path='../bin/dmscatter',
                cs = cs,
                cv = cv,
                cp = cp,
                cn = cn)

            # Plot the error w.r.t. DMFormFactor V6
            y = er_cut
            x = e_cut
            R = interp1d(energy, eventrate)
            Rx = R(x)
            perr = np.where(y!=0,abs((y - Rx)/y),0)
            ms=["--", ":", "-."]
            plt.plot(x, perr, label = "%s"%dreslbl[idres],
                    linestyle=ms[idres],marker="x")

            # Save dmscatter results
            shutil.copy("eventrate_spectra.dat.old",
                    "dmscatter.si29.%s.%i.%i.dat"%(coupleto,o1,o2))

        if o1==o2:
            plt.title("$^{29}$Si + 150 GeV WIMP [%s] %s"%(coupleto,
                operatorsymbols[o1-1]))
        else:
            plt.title("$^{29}$Si + 150 GeV WIMP [%s] $c_{%i}\cdot c_{%i}$"%(coupleto,
                o1,o2))
        plt.xscale("log")
        if log: plt.yscale("log")
        plt.xlabel("Recoil energy (keV)", fontsize=12)
        plt.ylabel("Ratio comparison of differential event rate", fontsize=12)
        plt.ylim(1e-6,1e1)

        if (all(eventrate==0)): plt.yscale("linear")
        plt.legend()
        plt.savefig("perr.%s.usd.%s.%i.%i.pdf" %(typ, coupleto, o1, o2))
        plt.close()            

        # Plot both sources
        plt.figure()
        if o1==o2:
            plt.title("$^{29}$Si + 150 GeV WIMP [%s] %s"%(coupleto,
                operatorsymbols[o1-1]))
        else:
            plt.title("$^{29}$Si + 150 GeV WIMP [%s] $c_{%i}\cdot c_{%i}$"%(coupleto,
                o1,o2))        
        plt.xscale("log")
        if log: plt.yscale("log")
        plt.xlabel("Recoil energy (keV)", fontsize=12)
        plt.ylabel("Event rate (1/MeV)", fontsize=12)
        if (all(eventrate==0)): plt.yscale("linear")  

        plt.plot(energy, eventrate, label="dmscatter, MB vesc %skm/s"%vesc,lw=2)
        emax500 = 265.2987
        plt.axvline(x=emax500, label="$E_{R,max}$ @ 550 km/s", linestyle=":",
                color='red')
        plt.plot(e_cut, er_cut,ms=1,
                linestyle="-.", label="DMFormFactor, MB vesc 550km/s",lw=2)
        plt.plot(e_mb, er_mb,ms=1,
                linestyle="--", label="DMFormFactor, MB")
        plt.legend(loc=3, fontsize=12)
        plt.savefig("%s.usd.%s.%i.%i.pdf"%(typ,coupleto,o1,o2))
        plt.close()
print("Validation script complete.")        
