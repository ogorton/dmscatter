import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import dmfortfactor as dm

###############################################################################
# This script generates a heat map of differential event rate spectra for a 
# range of WIMP masses.
###############################################################################

masses = np.arange(1.0, 300., 1) # GeV

# Options to do compute the data and plot the data independently
compute =  True
plot = True

# I will store the data we compute in files ending in this
filename_base = "GeV-rate.dat"

if compute:

    # Path to dmfortfactor executable (mine is in ~/bin/)
    exec_name = "dmfortfactor.x"
    
    # Template input file. In this case, it's just an input file since nothing is
    # modified.
    input_template = "xe131.input"
    
    workdir = "./"
    label="runCustom" # prefix for temporary files generatd by dm API
    
    # Template control file. This one does get modified and will be renamed to drop
    # the .template suffix.
    control_template = "xe131.control.template"
    
    # I'm doing an event rate calculation, so this is the filename this script
    # should look for.
    resultfile = "eventrate_spectra.dat"
    
    for WIMPMASS in masses:

        print("m = %s"%WIMPMASS)
    
        
        input_dict = {}
        control_dict = {"WIMPMASS" : WIMPMASS}
    
        E, R = dm.runTemplates(exec_name, 
            input_template, control_template,
            input_dict, control_dict, 
            workdir, label, resultfile)

        f = open(str(WIMPMASS)+filename_base, "w+")
        f.write("# Event rate (events/GeV) for %s GeV WIMPs."%WIMPMASS)
        f.write("# Recoil energy (keV)                # event rate")
        for ii, energy in enumerate(E):
            f.write("%10.5f    %s\n"%(energy, R[ii]))
        f.close()

if plot:

    # Just a text time with tex commands for convenient plot labels
    with open("operators.txt") as f: operatorsymbols = f.readlines()

    for ii,WIMPMASS in enumerate(masses):
        E, R = np.loadtxt(str(WIMPMASS)+filename_base, unpack = True)
        z = np.where(R==0, 1e-100, R) #prevent divide by zero in log10
        z = np.log10(z)
        if (ii==0): Z = z
        if (ii>0) : Z = np.vstack((Z,z))

    plt.figure()
    (X,Y) = np.meshgrid((E),masses)
    plt.pcolormesh(X,Y,Z,shading='nearest',cmap=cm.gnuplot2, vmin = -12)
    plt.colorbar(extend="both", label=r" log$_{10}$ Event rate (1/GeV)"
            )
    plt.ylabel("WIMP mass (GeV)")
    plt.xlabel(r"Recoil energy (keV)")
    plt.xscale('log')
    plt.savefig("wimp_mass_heatmap.png", dpi=600)
