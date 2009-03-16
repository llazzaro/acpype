import sys, os, time

from subprocess import Popen

pythonCommand = 'python2.5'
numCpu = 20

def runConversionJobs(chemCompVarFiles,scriptName):

    _timeStamp = time.strftime("%Y_%m_%d_%H_%M_%S")

    startCode = 'start'

    currentChemCompVarFile = None
    currentProcesses = {startCode: None}
    currentJobOut = {}
    currentIndex = -1
    endChemCompVarFile = chemCompVarFiles[-1]

    outputHandle = sys.__stdout__

    while (currentProcesses):

        if startCode in currentProcesses.keys():
            del(currentProcesses[startCode])

        if len(currentProcesses.keys()) < numCpu:

            tempIndex = currentIndex + 1
            for _i in range(currentIndex,currentIndex + numCpu - len(currentProcesses.keys())):
                # Don't start a job if it's at the end!
                if currentChemCompVarFile != endChemCompVarFile:

                    chemCompVarFile = chemCompVarFiles[tempIndex]

                    #
                    # TODO: ensure that stdout and stdin OK for this job!! Might want to reroute!!
                    #

                    currentOutFile = chemCompVarFile.replace('.mol2','.out')
                    jobOut = open(currentOutFile,'w')

                    process = Popen(['nice', '-19', pythonCommand, scriptName, '-i', chemCompVarFile, '-d'], stdout = jobOut, stderr = jobOut)

                    currentJobOut[chemCompVarFile] = jobOut
                    currentProcesses[chemCompVarFile] = process
                    currentChemCompVarFile = chemCompVarFile

                    outputHandle.write("\n*** Job %s started ***\n\n" % chemCompVarFile)

                    tempIndex += 1

                    # Allow jobs to start one by one... nicer for output
                    time.sleep(1)

            currentIndex = tempIndex - 1

        time.sleep(3)

        #
        # Check if finished
        #

        for chemCompVarFile in currentProcesses.keys():

            # Finished...
            if currentProcesses[chemCompVarFile].poll() != None:
                del(currentProcesses[chemCompVarFile])
                currentJobOut[chemCompVarFile].close()
                del(currentJobOut[chemCompVarFile])
                outputHandle.write("\n *** Job %s finished ***\n" % chemCompVarFile)
                if not currentProcesses:
                    currentProcesses = {startCode: None}

        outputHandle.flush()

if __name__ == '__main__':

    # Only run this on 'other'!
    chemCompVarFiles = []

    curDir = os.getcwd()

    ccpCodes = os.listdir('other')[:5]
    ccpCodes.sort()

    for ccpCode in ccpCodes:
        # HACK to restart after powercut
        #if ccpCode < 'NA':
        #  continue

        ccvNames = os.listdir(os.path.join('other',ccpCode))

        for ccvName in  ccvNames:
            if ccvName[-5:] == '.mol2' and not ccvName.count("bcc_gaff"):
                chemCompVarFiles.append(os.path.join(curDir,'other',ccpCode,ccvName))

    print chemCompVarFiles
#    runConversionJobs(chemCompVarFiles,'../acpypi.py')
