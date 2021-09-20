import sys
sys.path.append("/Users/oliver/projects/darkmatter/python")

import numpy as np
import matplotlib.pyplot as plt
import dmfortfactor as dm

exec_name = "dmfortfactor.x"
input_template = "c12.input.template"
dresfile = 'c12Nmax8chi20hw.dres'

hofrequencies = np.arange(15., 25., 1)

control_template = "c12.control.template"
operators = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

for operator in operators:
    print("Operator-%i"%operator)
    first=True
    plt.figure(1)
    plt.clf()
    plt.figure(2)
    plt.clf()
    for i, hofrequency in enumerate(hofrequencies):

        rlabel = 'run_%s_%s'%(operator,hofrequency)
        control_dict = {"OPERATOR" : operator,
            "HOFREQUENCY" : hofrequency}
        input_dict = {"DRESFILE" : dresfile[:-5],
                "LABEL": rlabel} # remove .dres extension
        
        E, R = dm.runTemplates(exec_name, input_template, control_template, input_dict,
            control_dict, workdir='./',
            label=rlabel,
            resultfile='eventrate_spectra.dat')
        if (first):
            first=False
            R0 = R
        plt.figure(1)
        plt.plot(E,R,label="hw = %s MeV"%hofrequency)
        plt.figure(2)
        plt.plot(E,abs(R - R0)/R0,label="hw = %s MeV"%hofrequency)
    plt.figure(1) 
    plt.title("C-12, 500 GeV WIMP, Op: %s, neutron 4.8E-4"%operator) 
    plt.xlabel('$E_{recoil}$ (keV)')
    plt.ylabel('Events/MeV')
    plt.legend(loc=3)
    plt.xscale('log')
    if (not all(R==0)):
        plt.yscale('log')
        zeros=False
    else:
        zeros=True        
    plt.savefig('c12-hw-o%s.pdf'%operator)
    if (zeros): continue
    plt.figure(2)
    plt.title("C-12, 500 GeV WIMP, Op: %s, neutron 4.8E-4"%operator)
    plt.xlabel('$E_{recoil}$ (keV)')
    plt.ylabel('Relative error w.r.t. hw=%s'%hofrequencies[0])
    plt.legend(loc=2)
    plt.xscale('log')
    if (not all(R==0)):plt.yscale('log')
    plt.savefig('c12-hw-o%s-relerr.pdf'%operator)
