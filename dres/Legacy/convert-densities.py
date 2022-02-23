import numpy as np
'''
    This script converts the pre-processed nuclear density matrix elements
    stored in the file legacy.data into the .dres format used by the BIGSTICK
    shell model code.

    Legacy density matrix format:
        ME, 2Jt, Tt, Nbra, 2Jbra, Nket, 2Jket
    Example for Si29:
        (* Jt = 0, Tt = 0 *)
        sdSi29[1]={sdSi2900x0101=4.,         0,0,0,1,0,1};
        sdSi29[2]={sdSi2900x1111=4.,         0,0,1,1,1,1};
        sdSi29[3]={sdSi2900x1313=5.65685,    0,0,1,3,1,3};
        sdSi29[4]={sdSi2900x2121=1.6791266,  0,0,2,1,2,1};
        sdSi29[5]={sdSi2900x2323=0.8589358,  0,0,2,3,2,3};
        sdSi29[6]={sdSi2900x2525=5.8344745,  0,0,2,5,2,5};
        (* Jt = 1, Tt = 1 *)
        sdSi29[7]={sdSi2901x2121=0.5065078,  0,1,2,1,2,1};
        sdSi29[8]={sdSi2901x2323=0.1242974,  0,1,2,3,2,3};
        sdSi29[9]={sdSi2901x2525=0.1833635,  0,1,2,5,2,5};
        (* Jt = 1, Tt = 0 *)
        sdSi29[10]={sdSi2910x2121=0.5720814, 1,0,2,1,2,1};
        sdSi29[11]={sdSi2910x2321=0.0249074, 1,0,2,3,2,1};
        sdSi29[12]={sdSi2910x2123=-0.0249074,1,0,2,1,2,3};
        sdSi29[13]={sdSi2910x2323=-0.0803281,1,0,2,3,2,3};
        sdSi29[14]={sdSi2910x2523=0.1966218, 1,0,2,5,2,3};
        sdSi29[obj.n]={sdSi2910x2325=-0.1966218,1,0,2,3,2,5};
        sdSi29[16]={sdSi2910x2525=0.1obj.n2554, 1,0,2,5,2,5};
        (* Jt = 1, Tt = 1 *)
        sdSi29[17]={sdSi2911x2121=0.5528034, 1,1,2,1,2,1};
        sdSi29[18]={sdSi2911x2321=0.0142132, 1,1,2,3,2,1};
        sdSi29[19]={sdSi2911x2123=-0.0142132,1,1,2,1,2,3};
        sdSi29[20]={sdSi2911x2323=-0.0875390,1,1,2,3,2,3};
        sdSi29[21]={sdSi2911x2523=0.1856207, 1,1,2,5,2,3};
        sdSi29[22]={sdSi2911x2325=-0.1856207,1,1,2,3,2,5};
        sdSi29[23]={sdSi2911x2525=0.1110496, 1,1,2,5,2,5};    
    Equivalent .dres format
         Initial state #    1 E = -144.35397 2xJ, 2xT =    1   1
         Final state   #    1 E = -144.35397 2xJ, 2xT =    1   1
         Jt =   0, Tt = 0        1
            1    1   0.8589358   0.1242974
            2    2   5.8344745   0.1833635
            3    3   1.6791266   0.5065078
         Jt =   1, Tt = 0        1
            1    1  -0.0803281  -0.0875390
            1    2  -0.1966218  -0.1856207
            1    3   0.0249074   0.0142132
            2    1   0.1966218   0.1856207
            2    2   0.1obj.n2554   0.1110496
            3    1  -0.0249074  -0.0142132
            3    3   0.5720814   0.5528034    
'''
elements = (
        "sdF",
        "sdSi28",
        "sdSi29",
        "sdSi30",
        "sdgI",
        "pfGe70",
        "pfGe72",
        "pfGe73",
        "pfGe74",
        "pfGe76",
        "sdNa",
        "sdgXe128",
        "sdgXe129",
        "sdgXe130",
        "sdgXe131",
        "sdgXe132",
        "sdgXe134",
        "sdgXe136")
spaces = {
        "sdF"      :"sd",
        "sdSi28"   :"sd",
        "sdSi29"   :"sd",
        "sdSi30"   :"sd",
        "sdgI"     :"sdg",
        "pfGe70"   :"pf",
        "pfGe72"   :"pf",
        "pfGe73"   :"pf",
        "pfGe74"   :"pf",
        "pfGe76"   :"pf",
        "sdNa"     :"sd",
        "sdgXe128" :"sdg",
        "sdgXe129" :"sdg",
        "sdgXe130" :"sdg",
        "sdgXe131" :"sdg",
        "sdgXe132" :"sdg",
        "sdgXe134" :"sdg",
        "sdgXe136" :"sdg"}        
particles = {
        "sdF"      : (9,10),
        "sdSi28"   : (14,14),
        "sdSi29"   : (14,15),
        "sdSi30"   : (14,16),
        "sdgI"     : (53,74),
        "pfGe70"   : (32,38),
        "pfGe72"   : (32,40),
        "pfGe73"   : (32,41),
        "pfGe74"   : (32,42),
        "pfGe76"   : (32,44),
        "sdNa"     : (11,12),
        "sdgXe128" : (54,74),
        "sdgXe129" : (54,75),
        "sdgXe130" : (54,76),
        "sdgXe131" : (54,77),
        "sdgXe132" : (54,78),
        "sdgXe134" : (54,80),
        "sdgXe136" : (54,82)}        
spins = {
        "sdF"      : 0.5,
        "sdSi28"   : 0.0,
        "sdSi29"   : 0.5,
        "sdSi30"   : 0.0,
        "sdgI"     : 2.5,
        "pfGe70"   : 0.0,
        "pfGe72"   : 0.0,
        "pfGe73"   : 4.5,
        "pfGe74"   : 0.0,
        "pfGe76"   : 0.0,
        "sdNa"     : 1.5,
        "sdgXe128" : 0.0,
        "sdgXe129" : 0.5,
        "sdgXe130" : 0.0,
        "sdgXe131" : 1.5,
        "sdgXe132" : 0.0,
        "sdgXe134" : 0.0,
        "sdgXe136" : 0.0}
class dres:
    def __init__(self,name):
        self.name = name
        self.space = spaces[self.name]
        self.maxjt = -1
        self.minjt = 1000
        self.spos = []
        self.n = 100
        self.nspo = 100
        self.declaredlength = ""
        self.matrixby2j = np.zeros((self.n,self.n,self.n,self.n))
        self.matrixbyspo = np.zeros((self.n,self.n,self.nspo,self.nspo))        
    def setme(self,datastring):
        jt, tt, lf, jf, li, ji = [int(i) for i in datastring[4:10]]
        self.maxjt = max(self.maxjt,jt)
        self.minjt = min(self.minjt,jt)
        me = float(datastring[3])
        self.matrixby2j[jt, tt, jf, ji] = me
        a = self.idspo(lf, jf)
        b = self.idspo(li, ji)
        self.matrixbyspo[jt,tt,a,b] = me
    def addspo(self,spo):
        self.spos.append(spo)
    def idspo(self,l,jx2):
        ID = -1
        for a, spo in enumerate(self.spos):
            if spo.matches(l,jx2): 
                ID = a
        return ID

class spo:
    def __init__(self,l,jx2,):
        self.l = l
        self.jx2 = jx2
        self.j = jx2 / 2.0
        self.weight = 1
    def matches(self,l,jx2):
        if (self.l==l and self.jx2==jx2) :
            return True
        else: 
            return False
def parsespo(datastring,fi):
    jt, tt, lf, jf, li, ji = [int(i) for i in datastring[4:10]]
    if fi=='i': 
        return li, ji
    elif fi=='f':
        return lf, jf
    else:
        return -1
def parsetransition(datastring):
    jt, tt, lf, jf, li, ji = [int(i) for i in datastring[4:10]]
    return jt, tt

alldres = {}
for iso in elements:
    obj = dres(iso)
    alldres[iso] = obj

with open("legacy.data") as f:
    lines = f.readlines()
    for line in lines:
        words = line.split(",")
        iso = words[0]
        if iso not in elements: 
            continue
        key = words[1]
        if key=='length':
            alldres[iso].declaredlength = words[2]
            continue        
        for fi in ('f', 'i'):
            l, jx2 = parsespo(words,fi)
            spoid = alldres[iso].idspo(l,jx2)
            if spoid == -1:
                alldres[iso].addspo(spo(l,jx2))
        alldres[iso].setme(words)

nme = 0
for iso in elements:
    print("\n")
    obj = alldres[iso]
    nspos = len(obj.spos)
    fdres= open("%s.dres"%obj.name,"w+")
    nme = 0    
    print("Converting dres data for %s..."%obj.name)
    fdres.write("Legacy dres file for %s\n"%obj.name)
    Z, N = particles[obj.name]
    T = (N - Z)/2.0
    gsj = int(2*spins[obj.name])
    gst = int(2*T)    
    fdres.write("  %3i  %3i\n"%(Z, N))
#    fdres.write("Converted from legacy format by automated script.\n")
#    fdres.write("The following data is not provided in the legacy files:\n")
#    fdres.write("  Nuclear state Energy, 2xJ, 2xT\n")
#    fdres.write("  Single particle state N, L quantum numbers\n")
    fdres.write("\n")
    fdres.write("  State      E        Ex         J       T\n")
    fdres.write("    1      0.0      0.00000     %s      %s\n"%(spins[obj.name],T))
    fdres.write("\n\n\n")

    fdres.write("  Single particle state quantum numbers\n")
    fdres.write("ORBIT      N     L   2 x J\n")
    for a in range(0,nspos):
        spoa = obj.spos[a]
        fdres.write("%6i  %4i  %4i  %4i\n"%(a+1, 0, spoa.l, spoa.jx2)) 

    fdres.write("\n\n\n")
    fdres.write(" Initial state #    1 E = 0.0 2xJ, 2xT =   %2i  %2i\n"%(gsj,gst))
    fdres.write(" Final state   #    1 E = 0.0 2xJ, 2xT =   %2i  %2i\n"%(gsj,gst))
    for jt in range(obj.minjt,obj.maxjt+1):
        fdres.write(" Jt =%4i, Tt = 0        1\n"%jt)
        for a in range(0,nspos):
            for b in range(0,nspos):
                x = obj.matrixbyspo[jt,0,a,b]
                y = obj.matrixbyspo[jt,1,a,b]
                if (x==y and y==0): 
                    continue
                nme += np.count_nonzero((x,y))
                fdres.write("%5i  %5i  %10.8f  %10.8f\n"%(a+1,b+1,x,y))
    fdres.close()
    print("Number of declared matrix elements: %s"%obj.declaredlength)
    print("Number of found matrix elements:    %s"%nme)









