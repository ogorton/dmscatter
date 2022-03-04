import sys
sys.path.append("../python")

import numpy as np
import matplotlib.pyplot as plt
import dmfortfactor as dm

datafiles = [
        '../data/C/c12Nmax8chi20hw.data',
        '../data/C/c12Nmax6DAEhw22.5.data',
        '../data/C/c12Nmax6chi20hw.data',
        '../data/C/c12Nmax8DAEhw22.5.data', 
        '../data/C/c12ck.data']
hofrequencies = [20.0, 22.5, 20.0, 22.5, "na"]
hofrequencies = len(datafiles) * ["na"]
operators = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

with open("operators.txt") as f: operatorsymbols = f.readlines()

for opi,operator in enumerate(operators):
    print("Operator-%i"%operator)
    first=True
    plt.figure(1)
    plt.clf()
    plt.figure(2)
    plt.clf()
    for i, datafile in enumerate(datafiles):

        rlabel = 'run_%s_%s'%(operator,hofrequencies[i])
        
        control_dict = {
                "hofrequency" : hofrequencies[i],
                "wimpmass" : 500.0,
                "vearth" : 232.0,
                "maxwellv0" : 220.0,
                "vescape" : 550.0,
                "dmdens" : 0.3,
                "ntscale" : 2500.0
            }

        cn = np.zeros(15)
        cn[operator-1] = 1.0

        E, R = dm.EventrateSpectra(
                Z = 6,
                N = 6,
                dres = datafile[:-5],
                controlwords = control_dict,
                cn = cn,
                exec_path='../bin/dmfortfactor')


        if (first):
            first=False
            R0 = R
        plt.figure(1)
        plt.plot(E,R,label=datafile[3:-5])
        plt.figure(2)
        plt.plot(E,abs(R - R0)/R0,label=datafile[3:-5])
    plt.figure(1) 
    plt.title("C-12, 500 GeV WIMP, Op: %s, neutron 4.8E-4"%operator) 
    plt.xlabel('$E_{recoil}$ (keV)')
    plt.ylabel('Events/MeV')
    plt.legend(loc=1)
    plt.xscale('log')
    if (not all(R==0)):
        plt.yscale('log')
        zeros=False
    else:
        zeros=True        
    plt.savefig('c12-o%s.pdf'%operator)
    if (zeros): continue
    plt.figure(2)
    plt.title("%s C-12, 500 GeV WIMP, $c_{%i}^n=$ 4.8E-4"%(operators[opi],operator))
    plt.xlabel('$E_{recoil}$ (keV)')
    plt.ylabel('Relative error w.r.t. c12Nmax8chi20hw')
    plt.legend(loc=2)
    plt.xscale('log')
    if (not all(R==0)):plt.yscale('log')
    plt.savefig('c12-o%s-relerr.pdf'%operator)
