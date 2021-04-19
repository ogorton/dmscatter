import numpy as np
import matplotlib.pyplot as plt
import dmf90factor as dm

exec_name = "dmf90factor.x"
input_template = "xe131.input.template"
dresfiles = [
    'xe131gcn.dres',
    'xe131jj55.dres',
    'pnism.xe131.dres']
control_template = "xe131.control.template"
operators = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
with open("operators.txt") as f: operatorsymbols = f.readlines()

for opi, operator in enumerate(operators):
    print("Operator-%i"%operator)
    control_dict = {"OPERATOR" : operator}
    first=True
    plt.figure(1)
    plt.clf()
    plt.figure(2)
    plt.clf()
    for dresfile in dresfiles:
        input_dict = {"DRESFILE" : dresfile[:-5]} # remove .dres extension
        
        E, R = dm.runTemplates(exec_name, input_template, control_template, input_dict,
            control_dict, workdir='./', label='runCustom',
            resultfile='eventrate_spectra.dat')
        if (first):
            first=False
            R0 = R
        plt.figure(1)
        plt.plot(E,R,label=dresfile[:-5])
        plt.figure(2)
        plt.plot(E,abs(R - R0)/R0,label=dresfile[:-5])
    plt.figure(1) 
    plt.title("%s Xe-131, 500 GeV WIMP, $c^n_%s=$4.8E-4"%(operatorsymbols[opi],
        operator)) 
    plt.xlabel('$E_{recoil}$ (keV)')
    plt.ylabel('Events/MeV')
    plt.legend(loc=1)
    plt.xscale('log')
    if (not all(R==0)):
        plt.yscale('log')
        zeros=False
    else:
        zeros=True        
    plt.savefig('xe131-o%s.pdf'%operator)
    if (zeros): continue
    plt.figure(2)
    plt.title("Op: %s, Xe-131, 150 GeV WIMP, neutron 4.8E-4"%operator)
    plt.xlabel('$E_{recoil}$ (keV)')
    plt.ylabel('Relative error w.r.t. gcn')
    plt.legend(loc=2)
    plt.xscale('log')
    if (not all(R==0)):plt.yscale('log')
    plt.savefig('xe131-o%s-relerr.pdf'%operator)
