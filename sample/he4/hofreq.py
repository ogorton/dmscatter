import numpy as np
import matplotlib.pyplot as plt
import dmf90factor as dm

exec_name = "dmf90factor.x"
input_template = "he4.input.template"

bs = np.arange(.1,60,.5) # range of h.o. parameters

control_template = "he4.control.template"
operators = [1]#, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
op = [ "$1_\chi 1_N$"]

for o, operator in enumerate(operators):
    print("Operator-%i"%operator)
    first=True
    plt.figure(1)
    plt.clf()
    plt.figure(2)
    plt.clf()
    ERvsHO = []
    #for i, hofrequency in enumerate(hofrequencies):
    for i, b in enumerate(bs):

        # Setup
        hofrequency = (6.43/b)**2
        rlabel = 'run_%s_%s'%(operator,b)
        control_dict = {"OPERATOR" : operator,
            "HOFREQUENCY" : hofrequency}
        input_dict = {"LABEL": rlabel} # remove .dres extension
        
        # Execute
        E, R = dm.runTemplates(exec_name, input_template, control_template, input_dict,
            control_dict, workdir='./',
            label=rlabel,
            resultfile='eventrate_spectra.dat')

        ERvsHO.append(R)

    plt.figure(1) 
    plt.scatter(bs, ERvsHO)
    plt.title("4He filled 0s1/2, 500 GeV WIMP, O=%s, neutron 4.8E-4"%op[o]) 
    plt.xlabel('$b$ (fm)')
    plt.ylabel('Events/MeV')
    plt.savefig('he4_hodependence_o%s.pdf'%operator)
