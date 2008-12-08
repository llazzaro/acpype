#!/usr/bin/env python

import os

results = {}

molDirs = os.listdir('.')

totalPdbOkCount = 0
totalIdealOkCount = 0
emptyDir = []
failedPdb = []
failedIdeal = []
failedBoth = []

ET1 = []
ET2 = []
ET3 = []
ET4 = []
ET6 = []
ET7 = []
ET8 = []

def analyseFile(mol, structure, file):
    '''
        warnType 1 = 'no 'babel' executable, no PDB file as input can be used!'
        warnType 2 = "In ..., residue name will be 'R instead of ... elsewhere"
        warnType 3 = 'no charge value given, trying to guess one...'
        warnType 4 = 'residue label ... is not all UPPERCASE'
        warnType 5 = 'this may raise problem with some applications like CNS'

        errorType 1 = 'guessCharge failed'
        errorType 2 = 'Atoms with same coordinates in'
        errorType 3 = 'Atoms TOO close'
        errorType 4 = 'Atoms TOO alone'
        errorType 5 = "Use '-f' option if you want to proceed anyway. Aborting"
        errorType 6 = 'Antechamber failed'
        errorType 7 = 'Parmchk failed'
        errorType 8 = 'Tleap failed'
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

            else:
                print file, line
    return "%i W, %i E, WT %s, ET %s" % (countWarn, countError, warnTypes, errorTypes)

for molDir in molDirs[:1000]:
    files = os.listdir(molDir)
    results[molDir] = []
    pdb = False
    ideal = False
    if not file:
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
    if not pdb: failedPdb.append(molDir)
    if not ideal: failedIdeal.append(molDir)
    if not pdb and not ideal: failedBoth.append(molDir)

numDirs = len(molDirs)
invPerNumDirs = 1/numDirs * 100.0
perPdb = totalPdbOkCount * invPerNumDirs
perIdeal= totalIdealOkCount * invPerNumDirs

print "Total molecules: %i\nPDB OK: %i (%3.2f%%)\nIdeal OK: %i (%3.2f%%)" % (numDirs, totalPdbOkCount, perPdb, totalIdealOkCount, perIdeal)
print "Failed Both: %i\t %s" % (len(failedBoth), str(failedBoth))
print "Failed PDB: %i\t\t %s" % (len(failedPdb), str(failedPdb))
print "Failed Ideal: %i\t %s" % (len(failedIdeal), str(failedIdeal))
print "guessCharge failed: %i\t %s" % (len(ET1), str(ET1))
print "Atoms same coords: %i\t %s" % (len(ET2), str(ET2))
print "Atoms TOO close: %i\t %s" % (len(ET3), str(ET3))
print "Atoms TOO alone: %i\t %s" % (len(ET4), str(ET4))
print 'Antechamber failed: %i\t %s' % (len(ET6), str(ET6))
#print 'Parmchk failed: %i\t %s' % (len(ET7), str(ET7))
print 'Tleap failed: %i\t %s' % (len(ET8), str(ET8))

#print results