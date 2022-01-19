#This code is designed to repeat the results of section 7.5.1 of
#Larsen's PhD thesis
import sys
sys.path.append("../python")
import dmfortfactor as dm
import numpy as np
import matplotlib.pyplot as plt
import random
#inistantiate couplings 
cp = np.zeros(15)
cn = np.zeros(15)
c1_p=1
ratio=0.6
c1_n=ratio*c1_p
#only M and PHI''M interct
#M is the only non-zero term with all others set to 0
#keep spin jx=1/2
#--incase I want others to be non-zero----
jx=0.5#default
cl_jx=4*jx*(jx+1)/3
m_n=0.938272#default
#---------
cwords = {
        "wimpmass" : 50.0,
        "usemomentum" : 0 
        }#"usemomentum": 1<- want default 0 as use recoil energy.

Wfunc = dm.NucFormFactor(
        Z = 54,
        N = 77,
        dres = "../dres/xe131gcn",
        controlwords = cwords,
        epmin = 1.0, 
        epmax = 1000.0, 
        epstep = 1.0,
        exec_path='../bin/dmfortfactor')

test_ratios=np.arange(0,100,1)/100
#print("assuming kev")
#print("0,0"+str(Wfunc(30)[0,0,0]))
#print("0,1"+str(Wfunc(30)[0,0,1]))
#print("1,0"+str(Wfunc(30)[0,1,0]))
#print("1,1"+str(Wfunc(30)[0,1,1]))
print("assuming Kev")
print("0,0 "+str(Wfunc(30)[0,0,0]))
print("0,1 "+str(Wfunc(30)[0,0,1]))
print("1,0 "+str(Wfunc(30)[0,1,0]))
print("1,1 "+str(Wfunc(30)[0,1,1]))



