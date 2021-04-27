import numpy as np
import matplotlib.pyplot as plt
import dmfortfactor as dm

with open("operators.txt") as f: operatorsymbols = f.readlines()

exec_name = "dmfortfactor.x"
input_template = "xe131-c7c8.input.template"
control_template = "xe131-c7c8.control.template"


angle = np.arange(0,np.pi/2.0,np.pi/2.0/9)
angle = np.linspace(0, np.pi/2, 20)
#angle = [0, 30, 60, 85, 90]
angle = np.array(angle) * 2*np.pi / 360.

angle = np.linspace(0, np.pi/2, 20)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

for i, a in enumerate(angle):

    r = 0.00048
    c7 = np.sin(a)
    c8 = np.cos(a)

    control_dict = {"c7" : c7, "c8" : c8}
    input_dict = {}

    Er, dRdE = dm.runTemplates(exec_name, input_template, control_template, input_dict,
        control_dict, workdir='./')

    theta = np.ones(len(Er)) * a

    #plt.plot(Er,dRdE,label=r"$\theta=%1.2f^{\circ}$"%(a*180/np.pi))
    print("min %s max %s"%(min(dRdE),max(dRdE)))
    ax.plot(theta, np.log(Er), np.log(dRdE),color='black')

#plt.title(r"$c_7 \mathcal{O}_7 + c_8 \mathcal{O}_8$")
#plt.text(10.,10**(-6),"%s %s with $c_7^2+c_8^2=1$"%(operatorsymbols[7-1],operatorsymbols[8-1]))

ax.set_xlabel(r"$\theta$")
ax.set_ylabel(r"$log(E_r)$ (keV)")
ax.set_zlabel(r"$log(dR_D/dE_r)$ (1/GeV)")

plt.show()
plt.savefig("xe131c7c8.pdf")
