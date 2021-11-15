import sys
sys.path.append('/Users/oliver/projects/darkmatter/python')
import dmfortfactor as dm
import numpy as np
import matplotlib.pyplot as plt
import random

cwords = {
        "wimpmass" : 150.0,
        "usemomentum": 1}

Wfunc = dm.NucFormFactor(
        Z = 54,
        N = 77,
        dres = "../dres/xe131gcn",
        controlwords = cwords,
        epmin = 0.001,
        epmax = 10.0,
        epstep = 0.001)

q = 0.001
print("q = %10.5f"%q)
print("W_i^{tau,tau_prime}(q) = ")
print(Wfunc(q))
