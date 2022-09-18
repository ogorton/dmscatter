import sys
sys.path.append("../python")

import numpy as np
import matplotlib.pyplot as plt
import dmscatter as dm

dresfiles = [
        '../targets/He/he4n3lo_nmax16hw26']
hofrequencies = [26.0]
operators = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

with open("operators.txt") as f: operatorsymbols = f.readlines()

for opi,operator in enumerate(operators):
    print("Operator-%i"%operator)
    plt.figure(1)
    plt.clf()
    for i, dresfile in enumerate(dresfiles):

        rlabel = 'run_%s_%s'%(operator,hofrequencies[i])
        control_dict = {
                "hofrequency" : hofrequencies[i],
                "wimpmass" : 50.0,
                "vearth" : 232.0,
                "maxwellv0" : 220.0,
                "vescape" : 550.0,
                "dmdens" : 0.3,
                "ntscale" : 2500.0
            }

        cn = np.zeros(15)
        cn[operator-1] = 1.0

        E, R = dm.EventrateSpectra(
                Z = 2,
                N = 2,
                dres = dresfile,
                controlwords = control_dict,
                cp = cn,
                exec_path='../bin/dmscatter')

        plt.figure(1)
        plt.plot(E,R,label=dresfile[3:])
    plt.figure(1) 
    plt.title("He-4, 50 GeV WIMP, Op: %s, neutron 4.8E-4"%operator) 
    plt.xlabel('$E_{recoil}$ (keV)')
    plt.ylabel('Events/MeV')
    plt.legend(loc=1)
    plt.xscale('log')
    if (not all(R==0)):
        plt.yscale('log')
        zeros=False
    else:
        zeros=True        
    plt.savefig('he4-o%s.pdf'%operator)
