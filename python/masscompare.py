import os
import matplotlib.pyplot as plt
import runCustomControl

exec_name = "dmfortfactor.x"

# This script works by replacing MASS_KEYWORD with the desired value of the
# WIMP mass, then running the code with a control file "<>.control".

control_template = "xe131.control" # Filename ends in .template, this is the prefix.
inputfile = 'input.xe131'
workdir = os.getcwd()
label = "xe131"

plt.figure()

for wimp_mass in (50.,150.0, 500.,5000.):

    param_dict = { "MASS_KEYWORD" : wimp_mass }
    RecoilE, EventRate = runCustomControl.runCustomControl(exec_name, 
            inputfile, control_template, param_dict, workdir, label)

    plt.plot(RecoilE, EventRate, label="$m_\chi=$%2.2f"%wimp_mass)

plt.xlabel("Recoil energy [kev]")
plt.ylabel("Event rate")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("%s.WIMPmassCompare.pdf"%label)
plt.show()
