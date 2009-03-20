#!/usr/bin/env python

import os, sys

results = {}

ccpCodes = os.listdir('other')
ccpCodes.sort()

DirsPassed = []
DirsFailed = []

goodResults = ['2 W, 0 E, WT _1_3, ET ', '3 W, 0 E, WT _1_2_3, ET ',
               '6 W, 0 E, WT _1_3_4_5_4_5, ET ', '4 W, 0 E, WT _1_2_3_7, ET ',
               '7 W, 0 E, WT _1_3_7_4_5_4_5, ET ', '5 W, 0 E, WT _3_4_5_4_5, ET ',
               '6 W, 0 E, WT _2_3_4_5_4_5, ET ', '1 W, 0 E, WT _3, ET ',
               '2 W, 0 E, WT _2_3, ET ']

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
WT8 = set()

def analyseFile(mol, structure, file):
    '''
        warnType 0
        warnType 1 = 'no 'babel' executable, no PDB file as input can be used!'
        warnType 2 = "In ..., residue name will be 'R instead of ... elsewhere"
        warnType 3 = 'no charge value given, trying to guess one...'
        warnType 4 = 'residue label ... is not all UPPERCASE'
        warnType 5 = 'this may raise problem with some applications like CNS'
        warnType 6 = "Couldn't determine all parameters"
        warnType 7 = "The unperturbed charge of the unit ..."
        warnType 8 = 'There is a bond of ... angstroms between'

        errorType 0
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
                WT7.add('%s: %s' % (mol, charge))
            elif "There is a bond of" in line:
                warnTypes += '_8'
                _dist = round(float(line[27:37]))
                #WT8.add('%s_%s: %s' % (mol, structure, _dist))
                WT8.add('%s_%s' % (mol, structure))
            else:
                print "UNKNOWN WARN:", file, line
                warnTypes += '_0'
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
                print "UNKNOWN ERROR:", file, line
                errorTypes += '_0'
    out = "%i W, %i E, WT %s, ET %s" % (countWarn, countError, warnTypes, errorTypes)
    if out not in goodResults:
        print "%s *%s*" % (mol, out)

    return out

#for molDir in molDirs:
#    count = 0
#    files = os.listdir(molDir)
#    for file in files:
#        if '.pdb.mol2' in file:
#            count += 1
#        if count > 1: print file
#
#sys.exit()
os.chdir('other')

# Only run this on 'other'!
if len(sys.argv) > 1:
    args = sys.argv[1:]
    if args[0].startswith('n='):
        num = int(args[0][2:])
        ccpCodes = ccpCodes[:num]
    else:
        args.sort()
        ccpCodes = args


for molDir in ccpCodes:
#    files = os.listdir(molDir)
    files = []
    for dirpath, dirnames, filenames in os.walk(molDir):
        files += filenames
    results[molDir] = []
    pdb = False
    ideal = False
    pdbMol2 = False
    idealMol2 = False
    countPdbMol2 = 1
    countIdealMol2 = 1
    out1 = ''
    out2 = ''
    pass1 = False
    pass2 = False
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
            out1 = analyseFile(molDir, 'ideal', os.path.join(molDir, file))
            results[molDir].append('ideal: '+out1)
        elif 'pdb.out' in file:
            out2 = analyseFile(molDir, 'pdb', os.path.join(molDir, file))
            results[molDir].append('pdb: '+out2)
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
    if countIdealMol2 > 2: multIdealMol2.append(molDir)
    if countPdbMol2 > 2: multPdbMol2.append(molDir)
    if not pdbMol2: missPdbMol2.append(molDir)
    if not idealMol2: missIdealMol2.append(molDir)
    if not pdb: failedPdb.append(molDir)
    if not ideal: failedIdeal.append(molDir)
    if not pdb and not ideal: failedBoth.append(molDir)

    if out1 in goodResults and out2 in goodResults:
        DirsPassed.append(molDir)


numDirs = len(ccpCodes)
#invPerNumDirs = 1.0/numDirs * 100.0
perPdb = totalPdbOkCount / (totalPdbMol2Count * 0.01)
perIdeal= totalIdealOkCount / (totalIdealMol2Count * 0.01)

print "\n--------------------------------\nTotal molecules (Dirs): %i" % numDirs
print "\nPDB OK: %i of %i (%3.2f%%)" % (totalPdbOkCount,totalPdbMol2Count, perPdb)
print "\nIdeal OK: %i of %i (%3.2f%%)" % (totalIdealOkCount, totalIdealMol2Count, perIdeal)
emptyDir.sort()
print '\nComplete Empty Dirs (NO mol2 input files for either PDB or IDEAL):\n%i\t %s' % (len(emptyDir), str(emptyDir))

a, b, c = set(emptyDir), set(missPdbMol2), set(missIdealMol2)

diff1 = list(b.difference(a))
diff1.sort()
print '\nDirs missing pdb.mol2 input files (besides Empty Dirs above):\n%i\t %s' % (len(diff1), str(diff1))

diff3 = list(c.difference(a))
diff3.sort()
print '\nDirs missing ideal.mol2 input files (besides Empty Dirs above):\n%i\t %s' % (len(diff3), str(diff3))

intersec1 = list(set(diff3).intersection(set(diff1)))
intersec1.sort()
print '\nIntersection between Dirs (missing pdb.mol2 x missing ideal.mol2): %i\t %s' % (len(intersec1), str(intersec1))

multPdbMol2.sort()
print '\nMore than one entry pdb.mol2 per Dir: %i\t %s' % (len(multPdbMol2), str(multPdbMol2))
multIdealMol2.sort()
print '\nMore than one entry ideal.mol2 per Dir: %i\t %s' % (len(multIdealMol2), str(multIdealMol2))

d, e, f = set(failedBoth), set(failedPdb), set(failedIdeal)

diff2 = list(d.difference(a))
diff2.sort()
print "\nDirs failed for Both pdb.mol2 and ideal.mol2 input files (besides Empty Dirs above):\n%i\t %s" % (len(diff2), str(diff2))

diff4 = list(e.difference(d))
diff4.sort()
print "\nDirs failed for pdb.mol2 besides Both above:\n%i\t %s" % (len(diff4), str(diff4))

diff5 = list(f.difference(d))
diff5.sort()
print "\nDirs failed for ideal.mol2 besides Both above:\n%i\t %s" % (len(diff5), str(diff5))

intersec2 = list(set(diff5).intersection(set(diff4)))
intersec2.sort()
print '\nIntersection between Dirs (failed pdb.mol2 x failed ideal.mol2): %i\t %s' % (len(intersec2), str(intersec2))

ET1_pdb = []
ET1_ideal = []
for i in ET1:
    if 'ideal' in i:
        ET1_ideal.append(i.replace('_ideal', ''))
    else:
        ET1_pdb.append(i.replace('_pdb', ''))
Both_ET1 = list(set(ET1_pdb).intersection(set(ET1_ideal)))
Both_ET1.sort()
print "\nguessCharge failed (but ACPYPI may have fineshed) for Both PDB and Ideal:\n%i\t %s" % (len(Both_ET1), str(Both_ET1))

diff7 = list(set(ET1_pdb).difference(set(Both_ET1)))
diff7.sort()
print "\nguessCharge failed (but ACPYPI may have fineshed) for PDB (besides Both above):\n%i\t %s" % (len(diff7), str(diff7))

diff8 = list(set(ET1_ideal).difference(set(Both_ET1)))
diff8.sort()
print "\nguessCharge failed (but ACPYPI may have fineshed) for Ideal (besides Both above):\n%i\t %s" % (len(diff8), str(diff8))

ET2_pdb = []
ET2_ideal = []
for i in ET2:
    if 'ideal' in i:
        ET2_ideal.append(i.replace('_ideal', ''))
    else:
        ET2_pdb.append(i.replace('_pdb', ''))
Both_ET2 = list(set(ET2_pdb).intersection(set(ET2_ideal)))
Both_ET2.sort()

print "\nSERIOUS ERRORS: Atoms same coords for Both PDB and Ideal:\n%i\t %s" % (len(Both_ET2), str(Both_ET2))

diff11 = list(set(ET2_pdb).difference(set(Both_ET2)))
diff11.sort()
print "\nSERIOUS ERRORS: Atoms same coords for PDB (besides Both above):\n%i\t %s" % (len(diff11), str(diff11))

diff12 = list(set(ET2_ideal).difference(set(Both_ET2)))
diff12.sort()
print "\nSERIOUS ERRORS: Atoms same coords for Ideal (besides Both above):\n%i\t %s" % (len(diff12), str(diff12))

ET3_pdb = []
ET3_ideal = []
for i in ET3:
    if 'ideal' in i:
        ET3_ideal.append(i.replace('_ideal', ''))
    else:
        ET3_pdb.append(i.replace('_pdb', ''))
Both_ET3 = list(set(ET3_pdb).intersection(set(ET3_ideal)))
Both_ET3.sort()

print "\nSERIOUS ERRORS: Atoms TOO close:\n%i\t %s" % (len(Both_ET3), str(Both_ET3))

diff13 = list(set(ET3_pdb).difference(set(Both_ET3)))
diff13.sort()
print "\nSERIOUS ERRORS: Atoms TOO close for PDB (besides Both above):\n%i\t %s" % (len(diff13), str(diff13))

diff14 = list(set(ET3_ideal).difference(set(Both_ET3)))
diff14.sort()
print "\nSERIOUS ERRORS: Atoms TOO close for Ideal (besides Both above):\n%i\t %s" % (len(diff14), str(diff14))

ET4_pdb = []
ET4_ideal = []
for i in ET4:
    if 'ideal' in i:
        ET4_ideal.append(i.replace('_ideal', ''))
    else:
        ET4_pdb.append(i.replace('_pdb', ''))
Both_ET4 = list(set(ET4_pdb).intersection(set(ET4_ideal)))
Both_ET4.sort()

print "\nSERIOUS ERRORS: Atoms TOO alone:\n%i\t %s" % (len(Both_ET4), str(Both_ET4))

diff15 = list(set(ET4_pdb).difference(set(Both_ET4)))
diff15.sort()
print "\nSERIOUS ERRORS: Atoms TOO alone for PDB (besides Both above):\n%i\t %s" % (len(diff15), str(diff15))

diff16 = list(set(ET4_ideal).difference(set(Both_ET4)))
diff16.sort()
print "\nSERIOUS ERRORS: Atoms TOO alone for Ideal (besides Both above):\n%i\t %s" % (len(diff16), str(diff16))

WT8_pdb = []
WT8_ideal = []
for i in WT8:
    if 'ideal' in i:
        WT8_ideal.append(i.replace('_ideal', ''))
    else:
        WT8_pdb.append(i.replace('_pdb', ''))
Both_WT8 = list(set(WT8_pdb).intersection(set(WT8_ideal)))
Both_WT8.sort()

print "\nSERIOUS ERRORS: ACPYPI ran but bond TOO long for Both PDB and Ideal:\n%i\t %s" % (len(Both_WT8), str(Both_WT8))

diff20 = list(set(WT8_pdb).difference(set(Both_WT8)))
diff20.sort()
print "\nSERIOUS ERRORS: ACPYPI ran but bond TOO long for PDB:\n%i\t %s" % (len(diff20), str(diff20))

diff21 = list(set(WT8_ideal).difference(set(Both_WT8)))
diff21.sort()
print "\nSERIOUS ERRORS: ACPYPI ran but bond TOO long for Ideal:\n%i\t %s" % (len(diff21), str(diff21))

ET6_pdb = []
ET6_ideal = []
for i in ET6:
    if 'ideal' in i:
        ET6_ideal.append(i.replace('_ideal', ''))
    else:
        ET6_pdb.append(i.replace('_pdb', ''))
Both_ET6 = list(set(ET6_pdb).intersection(set(ET6_ideal)))
Both_ET6.sort()

print "\nSERIOUS ERRORS: Antechamber failed:\n%i\t %s" % (len(Both_ET6), str(Both_ET6))

diff17 = list(set(ET6_pdb).difference(set(Both_ET6)))
diff17.sort()
print "\nSERIOUS ERRORS: Antechamber failed for PDB (besides Both above):\n%i\t %s" % (len(diff17), str(diff17))

diff18 = list(set(ET6_ideal).difference(set(Both_ET6)))
diff18.sort()
print "\nSERIOUS ERRORS: Antechamber failed for Ideal (besides Both above):\n%i\t %s" % (len(diff18), str(diff18))

diff6 = list(set(ET8).difference(set(ET6)))
diff6.sort()
print '\nSERIOUS ERRORS: Tleap failed: %i\t %s' % (len(diff6), str(diff6))
lET9 = list(ET9)
lET9.sort()
print '\nSERIOUS ERRORS: Sintax error: %i\t %s' % (len(lET9), str(lET9))

WT6_pdb = set()
WT6_ideal = set()
for i in WT6:
    if '_ideal' in i:
        WT6_ideal.add(i.replace('_ideal',''))
    else:
        WT6_pdb.add(i.replace('_pdb', ''))
Both_WT6 = list(WT6_ideal.intersection(WT6_pdb))
diff9 = list(WT6_pdb.difference(Both_WT6))
diff10 = list(WT6_ideal.difference(Both_WT6))
Both_WT6 = list(Both_WT6)
Both_WT6.sort()
diff9.sort()
diff10.sort()
print '\nSERIOUS WARNINGS (acpypi fineshed but results are likely unreliable)'
print '\nNot all parameters calculated for Both PDB and Ideal:\n%i\t %s' % (len(Both_WT6), str(Both_WT6))
print '\nNot all parameters calculated for PDB (besides Both above):\n%i\t %s' % (len(diff9), str(diff9))
print '\nNot all parameters calculated for Ideal (besides Both above):\n%i\t %s' % (len(diff10), str(diff10))

print '\nIT SHOULD BE OK, BUT ONE MAY WANT TO DOUBLE CHECK THE CHARGES'
lWT7 = list(WT7)
lWT7.sort()
print '\nCharge not Zero: %i\t %s' % (len(lWT7), str(lWT7))

DirsPassed.sort()
print '\nMol (Dirs) OK:\n%i\t %s' % (len(DirsPassed), str(DirsPassed))
