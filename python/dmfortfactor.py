import os
import subprocess
import numpy as np

def EventrateSpectra(Z, N, dres, ermin=1, ermax=1000, erstep=1, 
        controlwords={}, cpvec=None, cnvec=None, csvec=None, cvvec=None,
        exec_path='dmfortfactor.x', name="CSspectra"):

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
            ermin
                Recoil energy minimum (keV)
            ermax
                Recoil energy maximum (keV)
            erstep
                Recoil energy step size (keV). Spectra will be produced for recoil
                energies from ermin to ermax in steps of erstep.
            controlwords
                Dictionary of control words. Keys must be valid dmfortfactor
                control keywords, and values must be numbers.
            cpvec
                Length-15 array of nonrelativistic proton- coupling coefficients
            cnvec
                Length-15 array of nonrelativistic neutron- coupling coefficients
            csvec
                Length-15 array of nonrelativistic isoscalar- coupling coefficients
            cvvec
                Length-15 array of nonrelativistic isovector- coupling coefficients
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

    inputfile = writeinput(name, Z, N, dres, ermin, ermax, erstep)
    controlfile = writecontrol(name, controlwords, cpvec, cnvec, csvec, cvvec)

    RecoilE, EventRate = runTemplates(
        exec_path, 
        inputfile,
        controlfile, 
        )

    return RecoilE, EventRate

def writeinput(name, Z, N, dres, ermin, ermax, erstep):
    # Create input file
    CSspectra_inputfilename = name + ".input"
    CSspectra_inputfile = open(CSspectra_inputfilename, "w+")
    CSspectra_inputfile.write("1\n")
    CSspectra_inputfile.write("%i\n"%Z)
    CSspectra_inputfile.write("%i\n"%N)
    CSspectra_inputfile.write("CSspectra\n")
    CSspectra_inputfile.write("%s\n"%dres)
    CSspectra_inputfile.write("%s %s %s\n"%(ermin, ermax, erstep))
    CSspectra_inputfile.close()
    return CSspectra_inputfilename

def writecontrol(name, controlword_dict, cpvec, cnvec, csvec, cvvec):

    # Create control file
    filename = name + ".control"

    f = open(filename, "w+")

    nonzero = False
    allcouplings = {'p':cpvec, 'n':cnvec, 's':csvec, 'v':cvvec}
    for coupling in allcouplings:
        cvec = allcouplings[coupling]
        if cvec is not None:
            if len(cvec) != 15:
                print("Error: '%s' coupling vector is length-%i. Should be 15."%coupling)
                exit()
            nonzeroOps = np.flatnonzero(cvec)
            for operator in nonzeroOps:
                nonzero = True
                f.write("coefnonrel  %2i  %s  %20.10f\n"%(operator+1, coupling,
                    cvec[operator]))
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
    findandreplace = "sed '"
    for keyword in control_dict.keys():
        findandreplace += "s/%s/%s/g;"%(keyword,control_dict[keyword])
    findandreplace += "' %s > %s.control"%(control_template, label)
    returned = subprocess.call(findandreplace, shell=True)    

    # Input file
    findandreplace = "sed '"
    for keyword in input_dict.keys():
        findandreplace += "s/%s/%s/g;"%(keyword,input_dict[keyword])
    findandreplace += "' %s > %s.input"%(input_template,label)
    returned = subprocess.call(findandreplace, shell=True)    

    # Execute
    command = "%s < %s.input > %s.output"%(exec_name,label,label)
    print(command)
    returned = subprocess.call(command,shell=True)
    if(returned != 0): print("Return code: %s"%returned)    

    # Collect output
    RecoilE, EventRate = np.loadtxt(resultfile, unpack=True, skiprows=1)
    os.chdir(path)
    return RecoilE, EventRate    
