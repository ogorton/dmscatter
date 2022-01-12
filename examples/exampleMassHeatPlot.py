import sys
sys.path.append("../python")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import dmfortfactor as dm
import timeit

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

    g = open(filename_base+".csv", "w+")
    g.write("# WIMP_mass, Recoil_energy, Event_rate\n")

    tzero = timeit.default_timer()
    for WIMPMASS in masses:

        print("m = %s"%WIMPMASS)
        # Control words with a variable wimp mass
        control_dict = {
                "wimpmass" : WIMPMASS,
                "vearth" : 232.0,
                "ntscale" : 2500.0,
                "maxwellv0" : 220.0
        }
        tstart = timeit.default_timer()

        # EFT coupling: isotor
        cv = np.zeros(15)
        cv[0] = 0.00048

        # Evaluate the model
        E, R = dm.EventrateSpectra(
            Z = 54, 
            N = 77, 
            dres = "../dres/xe131gcn", 
            epmin = 1.0, 
            epmax = 1000.0, 
            epstep = 1.0,
            controlwords = control_dict,
            cv = cv,
            exec_path='../src/dmfortfactor')

        exectime = timeit.default_timer() - tstart
        tcycle = timeit.default_timer()
        
        print("DMFortFactor exec time: %s s"%(exectime))
        print("Python exec time:       %s s"%(tcycle - tzero - exectime))
        tzero = tcycle 

        f = open(str(WIMPMASS)+filename_base, "w+")
        f.write("# Event rate (events/GeV) for %s GeV WIMPs."%WIMPMASS)
        f.write("# Recoil energy (keV)                # event rate")
        for ii, energy in enumerate(E):
            f.write("%10.5f    %s\n"%(energy, R[ii]))
        f.close()

        for ii, energy in enumerate(E):
            g.write("%s, %10.5f, %s\n"%(WIMPMASS,energy,R[ii]))
    g.close()

if plot:

    # Just a text file with tex commands for convenient plot labels
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
