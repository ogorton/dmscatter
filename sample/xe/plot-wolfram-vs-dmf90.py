import numpy as np
import matplotlib.pyplot as plt

fs=12

wolframresult = "xe131.spectra.wolfram.20210407.dat"
dmf90result = "xe131.spectra.dmf90.20210407B.dat"

werecoil, weventrate = np.loadtxt(wolframresult, comments="#", unpack=True)
ferecoil, feventrate = np.loadtxt(dmf90result, comments="#", unpack=True)

plt.figure()
plt.title("150 GeV WIMPs on $^{131}$Xe with $c_1^n$=0.00048")
plt.xlabel("Recoil energy (keV)", fontsize=fs)
plt.ylabel("Event rate (1/MeV)", fontsize=fs)
plt.xscale("log")
plt.yscale("log")
plt.xlim(1,1000)

plt.plot(werecoil, weventrate, linestyle=":", 
    label="DMFormFactor\n(Wolfram Language Code)")
plt.plot(ferecoil, feventrate, linestyle="--", label="DMF90Factor\n(This work)")

plt.legend(fontsize=fs)
plt.savefig("xe131.spectra.pdf")


