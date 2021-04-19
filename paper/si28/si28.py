import shutil
import numpy as np
import matplotlib.pyplot as plt
import dmf90factor as dm

exec_name = "dmf90factor.x"
input_template = "si28.input"
control_template = "si28.control.template"
with open("operators.txt") as f: operatorsymbols = f.readlines()

operators = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
#with open("operators.txt") as f: operatorsymbols = f.readlines()

for opi, operator in enumerate(operators):
    print("Operator-%i"%operator)
    plt.figure(opi)
    plt.title("$^{28}$Si + 150 GeV WIMP %s"%operatorsymbols[operator-1])
    control_dict = {"OPERATOR" : operator}
    input_dict = {}
    energy, eventrate = dm.runTemplates(exec_name, input_template,
            control_template, input_dict, control_dict)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Recoil energy (keV)", fontsize=12)
    plt.ylabel("Event rate (1/MeV)", fontsize=12)
    if (all(eventrate==0)): plt.yscale("linear")
    plt.plot(energy, eventrate, label="DMf90Factor")
    energy, eventrate = np.loadtxt("dmformfactor.si28.%i.dat"%(operator),
            unpack=True)
    plt.plot(energy, eventrate,
            linestyle=":", label="DMFormFactor")
    plt.legend(loc=3, fontsize=12)
    plt.savefig("compare.%i.pdf"%(operator))

    shutil.copy("eventrate_spectra.dat",
            "dmf90factor.si28.On-%i.dat"%operator)
