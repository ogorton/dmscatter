import os
import subprocess
import numpy as np

def EventrateSpectra(Z, N, dres, epmin=1, epmax=1000, epstep=1, 
        controlwords={}, cp=None, cn=None, cs=None, cv=None,
        exec_path='dmfortfactor', name=".eventratespectra"):

    '''
    Calls the dmfortfactor eventrate spectra function, which computes the
    differential WIMP-nucleus scattering event rate as a function of nuclear
    recoil energy.

    Required arguments:
        Z 
            Number of protons in the target nucleus
        N   
            Number of neutron in the target nucleus
        dres
            Filename ending in .dres (not including .dres) for the file
            containing the reduced one-body density matrices for the target
            nucleus wave function

    Optional arguments:
        epmin
            Recoil energy minimum (keV). Or, transfer momentum (gev/c) if 
            usemomentum control word set to 1.
        epmax
            Recoil energy maximum (keV). Or, transfer momentum (gev/c) if
            usemomentum control word set to 1.
        epstep
            Recoil energy step size (keV). Or, transfer momentum (gev/c) if 
            usemomentum control word set to 1. Spectra will be produced for 
            recoil energies from epmin to epmax in steps of epstep.    
        controlwords
            Dictionary of control words. Keys must be valid dmfortfactor
            control keywords, and values must be numbers.
        cp
            Length-15 array of nonrelativistic proton- coupling coefficients
        cn
            Length-15 array of nonrelativistic neutron- coupling coefficients
        cs
            Length-15 array of nonrelativistic isoscalar- coupling coefficients
        cv
            Length-15 array of nonrelativistic isotor- coupling coefficients
        exec_path
            Path to the executable for dmfortfactor
        name
            Name string assigned to temporary files

    Returns:
        RecoilE
            Array of recoil energies (keV)
        EventRate
            Array of differential event rates (events/GeV)
    '''

    inputfile = writeinput('er', name, Z, N, dres, epmin, epmax, epstep)
    controlfile = writecontrol(name, controlwords, cp, cn, cs, cv)

    RecoilE, EventRate = runTemplates(
        exec_path, 
        inputfile,
        controlfile, 
        label=name
        )

    return RecoilE, EventRate

def NucFormFactor(Z, N, dres, epmin=1, epmax=1000, epstep=1,
    controlwords={}, exec_path='dmfortfactor', name=".nucFFspectra"):

    from scipy.interpolate import interp1d
    
    inputfile = writeinput('ws', name, Z, N, dres, epmin, epmax, epstep)
    controlfile = writecontrol(name, controlwords)


    columns  = runTemplates(
        exec_path,
        inputfile,
        controlfile,
        label=name,
        resultfile='nucresponse_spectra.dat'
        )

    q = columns[0]
    W = columns[1:]

    function_lst = []
    for operator in range(0,8):
        for tau in range(0,2):
            for tau_prime in range(0,2):
                windx = operator + tau * 8 + tau_prime * 8 * 2
                function_lst.append(interp1d(q, W[windx,:]))

    def Wfunc(qq):

        result = np.zeros((8,2,2))
        for operator in range(0,8):
            for tau in range(0,2):
                for tau_prime in range(0,2):
                    windx = operator + tau * 8 + tau_prime * 8 * 2
                    f = function_lst[windx]
                    result[operator, tau, tau_prime] = f(qq)
        return result

    return Wfunc
    
##
## Helper funtions
##
def writeinput(option, name, Z, N, dres, epmin, epmax, epstep):
    # Create input file
    CSspectra_inputfilename = name + ".input"
    CSspectra_inputfile = open(CSspectra_inputfilename, "w+")
    CSspectra_inputfile.write("%s\n"%option)
    CSspectra_inputfile.write("%i\n"%Z)
    CSspectra_inputfile.write("%i\n"%N)
    CSspectra_inputfile.write("%s\n"%name)
    CSspectra_inputfile.write("%s\n"%dres)
    CSspectra_inputfile.write("%s %s %s\n"%(epmin, epmax, epstep))
    CSspectra_inputfile.close()
    return CSspectra_inputfilename

def writecontrol(name, controlword_dict, cp=None, cn=None, cs=None,
        cv=None):

    # Create control file
    filename = name + ".control"

    f = open(filename, "w+")

    nonzero = False
    allcouplings = {'p':cp, 'n':cn, 's':cs, 'v':cv}
    for coupling in allcouplings:
        c = allcouplings[coupling]
        if c is not None:
            if len(c) != 15:
                print("Error: '%s' coupling tor is length-%i. Should be 15."%coupling)
                exit()
            nonzeroOps = np.flatnonzero(c)
            for operator in nonzeroOps:
                nonzero = True
                f.write("coefnonrel  %2i  %s  %20.10f\n"%(operator+1, coupling,
                    c[operator]))
    if not nonzero: print("Warning: there were no nonzero EFT couplings!")

    for key in controlword_dict:
        f.write("%s    %s\n"%(key, controlword_dict[key]))
    f.close()
    return filename

def runTemplates(exec_name, input_template, control_template, input_dict={},
    control_dict={}, workdir='./', label='runCustom',
    resultfile='eventrate_spectra.dat'):

    """
     Author: Oliver Gorton, 2021

     This is a function for running an executable program with a custom input
     file and a custom control file. An input file is one which is used like
        exec_name < input_file
     while a control file is one which is read by exec_name during runtime.

     The input_template and control_template and modified by replacing keywords
     in a dictionary with their respective values.

     exec_name: (string) containing the path to and name of the
         executable you want to run
     input_template: (string) this is an inputfile template used by the 
         executable and is modified by this function
     control_template: (string) the name of a file ending in '.control'. which
         this function will edit.
     input_dict: (dictionary) of keywords and the values to which those keywords
         will be changed in the input file. The keywords are strings which should 
         appear in the input_template file. 
     control_dict: (dictionary) of keywords and the values to which those keywords
         will be changed in the control file. The keywords are strings which 
         should appear in the control_template file.         
     workdir: (string) name of the working directory where you want to call
         the executable and send the output.
     label: (string) name for this job, to label the input and output files 
         generated by this function.

    """    
    path = os.getcwd()
    os.chdir(workdir)

    if (not os.path.exists(input_template)):
        print("File %s does not exist."%input_template)
        exit()
    if (not os.path.exists("%s"%control_template)):
        print("File %s does not exist."%control_template)
        exit()

    # Control file
    if control_dict == {}:
        pass
    else:
        findandreplace = "sed '"
        for keyword in control_dict.keys():
            findandreplace += "s/%s/%s/g;"%(keyword,control_dict[keyword])
        findandreplace += "' %s > %s.control"%(control_template, label)
        returned = subprocess.call(findandreplace, shell=True)    

    # Input file
    if input_dict == {}:
        custominput = input_template
    else:
        findandreplace = "sed '"
        custominput = "%s.input"%label
        for keyword in input_dict.keys():
            findandreplace += "s/%s/%s/g;"%(keyword,input_dict[keyword])
        findandreplace += "' %s > %s"%(input_template,custominput)
        returned = subprocess.call(findandreplace, shell=True)    

    # Execute
    command = "%s < %s > %s.output"%(exec_name,custominput,label)
    print(command)
    returned = subprocess.call(command,shell=True)
    if(returned != 0): print("Return code: %s"%returned)    

    # Collect output
    columns = np.loadtxt(resultfile, unpack=True)
    os.chdir(path)
    return columns
