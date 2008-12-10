#!/usr/bin/env python

import os#, sys

results = {}

molDirs = os.listdir('.')

totalPdbOkCount = 0
totalIdealOkCount = 0
totalIdealMol2Count = 0
totalPdbMol2Count = 0

emptyDir = []
failedPdb = []
failedIdeal = []
failedBoth = []
missPdbMol2 = []
missIdealMol2 = []
multIdealMol2 = []
multPdbMol2 = []

ET1 = []
ET2 = []
ET3 = []
ET4 = []
ET6 = []
ET7 = []
ET8 = []
ET9 = set()

WT6 = []
WT7 = set()

def analyseFile(mol, structure, file):
    '''
        warnType 1 = 'no 'babel' executable, no PDB file as input can be used!'
        warnType 2 = "In ..., residue name will be 'R instead of ... elsewhere"
        warnType 3 = 'no charge value given, trying to guess one...'
        warnType 4 = 'residue label ... is not all UPPERCASE'
        warnType 5 = 'this may raise problem with some applications like CNS'
        warnType 6 = "Couldn't determine all parameters"
        warnType 7 = "The unperturbed charge of the unit ..."

        errorType 1 = 'guessCharge failed'
        errorType 2 = 'Atoms with same coordinates in'
        errorType 3 = 'Atoms TOO close'
        errorType 4 = 'Atoms TOO alone'
        errorType 5 = "Use '-f' option if you want to proceed anyway. Aborting"
        errorType 6 = 'Antechamber failed'
        errorType 7 = 'Parmchk failed'
        errorType 8 = 'Tleap failed'
        errorType 9 = 'syntax error'
    '''

    countWarn = 0
    warnTypes = ''
    errorTypes = ''
    countError = 0
    tmpFile = open(file, 'r')
    content = tmpFile .readlines()
    for line in content:
        if 'WARNING: ' in line:
            countWarn += 1
            if "no 'babel' executable, no PDB" in line:
                warnTypes += '_1'
            elif "residue name will be 'R" in line:
                warnTypes += '_2'
            elif "no charge value given" in line:
                warnTypes += '_3'
            elif "not all UPPERCASE" in line:
                warnTypes += '_4'
            elif "applications like CNS" in line:
                warnTypes += '_5'
            elif "Couldn't determine all parameters" in line:
                warnTypes += '_6'
                WT6.append('%s_%s'% (mol, structure))
            elif "The unperturbed charge of the" in line:
                warnTypes += '_7'
                charge = round(float(line[44:54]))
                WT7.add('%s: %s'% (mol, charge))
            else:
                print file, line
        if 'ERROR: ' in line:
            countError += 1
            if 'guessCharge failed' in line:
                errorTypes += '_1'
                ET1.append('%s_%s'% (mol, structure))
            elif 'Atoms with same coordinates in' in line:
                errorTypes += '_2'
                ET2.append('%s_%s'% (mol, structure))
            elif 'Atoms TOO close' in line:
                errorTypes += '_3'
                ET3.append('%s_%s'% (mol, structure))
            elif 'Atoms TOO alone' in line:
                errorTypes += '_4'
                ET4.append('%s_%s'% (mol, structure))
            elif "Use '-f' option if you want to proceed anyway. Aborting" in line:
                errorTypes += '_5'
            elif 'Antechamber failed' in line:
                errorTypes += '_6'
                ET6.append('%s_%s'% (mol, structure))
            elif 'Parmchk failed' in line:
                errorTypes += '_7'
                ET7.append('%s_%s'% (mol, structure))
            elif 'Tleap failed' in line:
                errorTypes += '_8'
                ET8.append('%s_%s'% (mol, structure))
            elif 'syntax error' in line:
                errorTypes += '_9'
                ET9.add('%s_%s'% (mol, structure))
            else:
                print file, line
    return "%i W, %i E, WT %s, ET %s" % (countWarn, countError, warnTypes, errorTypes)

#for molDir in molDirs:
#    count = 0
#    files = os.listdir(molDir)
#    for file in files:
#        if '.pdb.mol2' in file:
#            count += 1
#        if count > 1: print file
#
#sys.exit()

for molDir in molDirs[:]:
    files = os.listdir(molDir)
    results[molDir] = []
    pdb = False
    ideal = False
    pdbMol2 = False
    idealMol2 = False
    countPdbMol2 = 1
    countIdealMol2 = 1
    if not files:
        emptyDir.append(molDir)
    for file in files:
        if 'pdb_NEW.pdb' in file:
            pdb = True
            results[molDir].append('pdb_OK')
            totalPdbOkCount += 1
        elif 'ideal_NEW.pdb' in file:
            ideal = True
            results[molDir].append('ideal_OK')
            totalIdealOkCount += 1
        elif 'ideal.out' in file:
            out = analyseFile(molDir, 'ideal', os.path.join(molDir, file))
            results[molDir].append('ideal: '+out)
        elif 'pdb.out' in file:
            out = analyseFile(molDir, 'pdb', os.path.join(molDir, file))
            results[molDir].append('pdb: '+out)
        elif '.ideal.mol2' in file:
            if idealMol2:
                countIdealMol2 +=1
            else:
                idealMol2 = True
                totalIdealMol2Count += 1
        elif '.pdb.mol2' in file:
            if pdbMol2:
                countPdbMol2 +=1
            else:
                pdbMol2 = True
                totalPdbMol2Count += 1
    if countIdealMol2 > 1: multIdealMol2.append(molDir)
    if countPdbMol2 > 1: multPdbMol2.append(molDir)
    if not pdbMol2: missPdbMol2.append(molDir)
    if not idealMol2: missIdealMol2.append(molDir)
    if not pdb: failedPdb.append(molDir)
    if not ideal: failedIdeal.append(molDir)
    if not pdb and not ideal: failedBoth.append(molDir)

numDirs = len(molDirs)
#invPerNumDirs = 1.0/numDirs * 100.0
perPdb = totalPdbOkCount / (totalPdbMol2Count * 0.01)
perIdeal= totalIdealOkCount / (totalIdealMol2Count * 0.01)

print "Total molecules: %i" % numDirs
print "\nPDB OK: %i of %i (%3.2f%%)" % (totalPdbOkCount,totalPdbMol2Count, perPdb)
print "\nIdeal OK: %i of %i (%3.2f%%)" % (totalIdealOkCount, totalIdealMol2Count, perIdeal)
print '\nEmpty Dirs: %i\t %s' % (len(emptyDir), str(emptyDir))
a, b, c = set(emptyDir), set(missPdbMol2), set(missIdealMol2)
diff1 = list(b.difference(a))
print '\nMissing pdb.mol2 (- EmpytDirs): %i\t %s' % (len(diff1), str(diff1))
#diff2 = list(b.difference(a))
#print '\nMissing ideal.mol2 (- EmpytDirs): %i\t %s' % (len(diff2), str(diff2))
diff3 = list(c.difference(b))
print '\nMissing ideal.mol2 (- pdb.mol2): %i\t %s' % (len(diff3), str(diff3))

print '\nMulti pdb.mol2: %i\t %s' % (len(multPdbMol2), str(multPdbMol2))
print '\nMulti ideal.mol2: %i\t %s' % (len(multIdealMol2), str(multIdealMol2))

print "\nFailed Both: %i\t %s" % (len(failedBoth), str(failedBoth))
d, e, f = set(failedBoth), set(failedPdb), set(failedIdeal)
diff4 = list(e.difference(d))
print "\nFailed PDB (- Both): %i\t %s" % (len(diff4), str(diff4))
diff5 = list(f.difference(d))
print "\nFailed Ideal (- Both): %i\t %s" % (len(diff5), str(diff5))

print "\nguessCharge failed: %i\t %s" % (len(ET1), str(ET1))

print '\nSERIOUS ERRORS'
print "\nAtoms same coords: %i\t %s" % (len(ET2), str(ET2))
print "\nAtoms TOO close: %i\t %s" % (len(ET3), str(ET3))
print "\nAtoms TOO alone: %i\t %s" % (len(ET4), str(ET4))
print '\nAntechamber failed: %i\t %s' % (len(ET6), str(ET6))
diff6 = set(ET8).difference(set(ET6))
print '\nTleap failed (- AC): %i\t %s' % (len(diff6), str(diff6))
print '\nSintax error: %i\t %s' % (len(ET9), str(list(ET9)))

print '\nIT SHOULD BE OK, BUT ONE MAY WANT TO DOUBLE CHECK THE CHARGES'
print '\nCharge not Zero: %i\t %s' % (len(WT7), str(list(WT7)))

print '\nSERIOUS WARNINGS (acpypi fineshed but results are likely unreliable)'
print '\nNot all parameters: %i\t %s' % (len(WT6), str(WT6))

#print '\nMissing pdb.mol2: %i\t %s' % (len(missPdbMol2), str(missPdbMol2))
#print '\nMissing ideal.mol2: %i\t %s' % (len(missIdealMol2), str(missIdealMol2))

#print results
