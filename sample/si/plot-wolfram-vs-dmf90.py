import numpy as np
import matplotlib.pyplot as plt

wolframresult = "test.si29.dat"
dmf90result = "eventrate_spectra.dat"

werecoil, weventrate = np.loadtxt(wolframresult, comments="#", unpack=True)
ferecoil, feventrate = np.loadtxt(dmf90result, comments="#", unpack=True)

plt.figure()
plt.title("150 GeV WIMPs on $^{29}$Si with $c_1^n$=1.0")
plt.xlabel("Recoil energy (keV)")
plt.ylabel("Event rate (1/MeV)")
plt.xscale("log")
plt.yscale("log")
plt.xlim(1,1000)

plt.plot(werecoil, weventrate, linestyle=":", 
    label="DMFormFactor\n(Wolfram Language Code)")
plt.plot(ferecoil, feventrate, linestyle="--", label="DMF90Factor\n(This work)")

plt.plot(werecoil, abs(feventrate-weventrate), label='Abs. error')

plt.legend()
plt.show()


