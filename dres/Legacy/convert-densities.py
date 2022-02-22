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
        jt, tt, nf, jf, ni, ji = [int(i) for i in datastring[4:10]]
        self.maxjt = max(self.maxjt,jt)
        self.minjt = min(self.minjt,jt)
        me = float(datastring[3])
        self.matrixby2j[jt, tt, jf, ji] = me
        a = self.idspo(nf, jf)
        b = self.idspo(ni, ji)
        self.matrixbyspo[jt,tt,a,b] = me
    def addspo(self,spo):
        self.spos.append(spo)
    def idspo(self,nr,jx2):
        ID = -1
        for a, spo in enumerate(self.spos):
            if spo.matches(nr,jx2): 
                ID = a
        return ID

class spo:
    def __init__(self,nr,jx2,):
        self.nr = nr
        self.jx2 = jx2
        self.j = jx2 / 2.0
        self.weight = 1
        self.l = -1
    def matches(self,nr,jx2):
        if (self.nr==nr and self.jx2==jx2) :
            return True
        else: 
            return False
def parsespo(datastring,fi):
    jt, tt, nf, jf, ni, ji = [int(i) for i in datastring[4:10]]
    if fi=='i': 
        return ni, ji
    elif fi=='f':
        return nf, jf
    else:
        return -1
def parsetransition(datastring):
    jt, tt, nf, jf, ni, ji = [int(i) for i in datastring[4:10]]
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
            nr, jx2 = parsespo(words,fi)
            spoid = alldres[iso].idspo(nr,jx2)
            if spoid == -1:
                alldres[iso].addspo(spo(nr,jx2))
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
    fdres.write("Converted from legacy format by automated script.\n")
    fdres.write("The following data is not provided in the legacy files:\n")
    fdres.write("  Nuclear state Energy, 2xJ, 2xT\n")
    fdres.write("  Single particle state N, L quantum numbers\n")
    fdres.write("\n\n\n")
    fdres.write("  Single particle state quantum numbers\n")
    fdres.write("ORBIT      N     L   2 x J\n")
    for a in range(0,nspos):
        spoa = obj.spos[a]
        fdres.write("%6i  %4i  %4i  %4i\n"%(a+1, spoa.nr, spoa.l, spoa.jx2)) 

    fdres.write("\n\n\n")
    fdres.write(" Initial state #    1 E = 0.0 2xJ, 2xT =   -1  -1\n")
    fdres.write(" Final state   #    1 E = 0.0 2xJ, 2xT =   -1  -1\n")
    for jt in range(obj.minjt,obj.maxjt+1):
        fdres.write(" Jt = %4i, Tt = 0        1\n"%jt)
        for a in range(0,nspos):
            for b in range(0,nspos):
                x = obj.matrixbyspo[jt,0,a,b]
                y = obj.matrixbyspo[jt,1,a,b]
                if (x==y and y==0): 
                    continue
                nme += np.count_nonzero((x,y))
                fdres.write("%4i  %4i  %10.8f  %10.8f\n"%(a+1,b+1,x,y))
    fdres.close()
    print("Number of declared matrix elements: %s"%obj.declaredlength)
    print("Number of found matrix elements:    %s"%nme)









