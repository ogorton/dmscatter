import os
import subprocess

def runCustom(exec_name, controltemplate, pnames, pvalues, workdir, label):

    path = os.getcwd()
    os.chdir(workdir)

    # Replace pnames in the control file with pvalues
    findandreplace = "sed '"
    for i,keyword in enumerate(pnames):
        findandreplace += "s/%s/%s/g;"%(keyword,pvalues[i])
    findandreplace += "' %s > %s.control"%(controltemplate,label)

    returned = subprocess.call(findandreplace, shell=True)

    #Run with given settings
    command = "%s < input.%s > output.%s"%(exec_name,label,label)
    print(command)
    returned = subprocess.call(command,shell=True)
    os.chdir(path)

    return
