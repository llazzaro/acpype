#!/usr/bin/env python

import os, sys

global id
id = 1

results = {}

ccpCodes = os.listdir('other')
ccpCodes.sort()

DirsPassed = []
DirsFailed = []

goodResults = [
                '0 E, 1 W, ET , WT _0','0 E, 2 W, ET , WT _0_1',
                '0 E, 2 W, ET , WT _0_2','0 E, 3 W, ET , WT _0_1_2',
                '0 E, 2 W, ET , WT _0_3','0 E, 3 W, ET , WT _0_1_3',
                '0 E, 3 W, ET , WT _0_2_3','0 E, 4 W, ET , WT _0_1_2_3'
               ]

guessFailedOnly = [
                   '1 E, 1 W, ET _0, WT _0', '1 E, 2 W, ET _0, WT _0_1',
                   '1 E, 2 W, ET _0, WT _0_2', '1 E, 3 W, ET _0, WT _0_1_2'
                   ]

closeContact = [
                '0 E, 2 W, ET , WT _0_7', '0 E, 3 W, ET , WT _0_1_7',
                '0 E, 3 W, ET , WT _0_2_7', '0 E, 3 W, ET , WT _0_3_7',
                '0 E, 4 W, ET , WT _0_1_2_7', '0 E, 4 W, ET , WT _0_1_3_7',
                '0 E, 4 W, ET , WT _0_2_3_7', '0 E, 5 W, ET , WT _0_1_2_3_7'
                ]

bondIrreg = [
             '0 E, 2 W, ET , WT _0_5', '0 E, 3 W, ET , WT _0_1_5',
             '0 E, 3 W, ET , WT _0_2_5', '0 E, 4 W, ET , WT _0_1_2_5'
            ]

# BOND IRREGULAR / ATOMS IN CLOSE CONTACT
biAndAicc = [
             '0 E, 3 W, ET , WT _0_5_7', '0 E, 4 W, ET , WT _0_1_5_7',
             '0 E, 4 W, ET , WT _0_2_5_7', '0 E, 5 W, ET , WT _0_1_2_5_7',
             '0 E, 6 W, ET , WT _0_1_2_3_5_7'
             ]

missingPar = [
              '0 E, 3 W, ET , WT _0_2_4'
              ]

error_warn_messages = \
'''
    warnType 0 = 'no charge value given, trying to guess one...'
    warnType 1 = "In ..., residue name will be 'R instead of ... elsewhere"
    warnType 2 = 'residue label ... is not all UPPERCASE'

    # mild warning (need to be sure the charge guessed is correct)
    warnType 3 = "The unperturbed charge of the unit ... is not zero"

    # serious warnings (topology not reliable):
    warnType 4 = "Couldn't determine all parameters"
    warnType 5 = 'There is a bond of ... angstroms between'
    warnType 6 = 'atom type may be wrong'
    warnType 7 = 'Close contact of ... angstroms between ...'

    warnType 8 = 'no 'babel' executable, no PDB file as input can be used!'
    warnType 9 = 'UNKNOWN WARN'

    errorType 0 = 'guessCharge failed' # (not so serious if only err    or because mopac worked with charge Zero)
    errorType 1 = 'Atoms with same coordinates in'
    errorType 2 = 'Atoms TOO close'
    errorType 3 = 'Atoms TOO alone'
    errorType 4 = 'Antechamber failed'
    errorType 5 = 'Parmchk failed'
    errorType 6 = 'Tleap failed'
    errorType 7 = "No such file or directory: 'tmp'" # can be bondtyes wrong or wrong frozen atom type
    errorType 8 = 'syntax error'
    errorType 9 = 'UNKNOWN ERROR'
'''

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

ET0 = []
ET1 = []
ET2 = []
ET3 = []
ET4 = []
ET5 = []
ET6 = []
ET7 = set()
ET8 = set()

WT3 = set()
WT4 = []
WT5 = set()
WT6 = set()
WT7 = set()

mapResults = {}

def analyseFile(mol, structure, file):
    warnTypes = ''
    errorTypes = ''
    tmpFile = open(file, 'r')
    content = tmpFile .readlines()
    for line in content:
        if 'WARNING: ' in line:
            if "no charge value given" in line:
                warnTypes += '_0'
            elif "residue name will be 'R" in line:
                warnTypes += '_1'
            elif "not all UPPERCASE" in line:
                warnTypes += '_2'
            elif "applications like CNS" in line:
                pass
            elif "The unperturbed charge of the" in line:
                warnTypes += '_3'
                charge = round(float(line[44:54]))
                WT3.add('%s: %s' % (mol, charge))
            elif "Couldn't determine all parameters" in line:
                warnTypes += '_4'
                WT4.append('%s_%s'% (mol, structure))
            elif "There is a bond of" in line:
                warnTypes += '_5'
                _dist = round(float(line[27:37]))
                WT5.add('%s_%s' % (mol, structure))
            elif ' atom type of ' in line:
                warnTypes += '_6'
                WT6.add('%s_%s'% (mol, structure))
            elif "no 'babel' executable, no PDB" in line:
                warnTypes += '_8'
            else:
                print "UNKNOWN WARN:", file, line
                warnTypes += '_9'
        if 'Warning: Close contact of ' in line:
            warnTypes += '_7'
            WT7.add('%s_%s'% (mol, structure))

        if 'ERROR: ' in line:
            if 'guessCharge failed' in line:
                errorTypes += '_0'
                ET0.append('%s_%s'% (mol, structure))
            elif 'Atoms with same coordinates in' in line:
                errorTypes += '_1'
                ET1.append('%s_%s'% (mol, structure))
            elif 'Atoms TOO close' in line:
                errorTypes += '_2'
                ET2.append('%s_%s'% (mol, structure))
            elif 'Atoms TOO alone' in line:
                errorTypes += '_3'
                ET3.append('%s_%s'% (mol, structure))
            elif 'Antechamber failed' in line:
                errorTypes += '_4'
                ET4.append('%s_%s'% (mol, structure))
            elif 'Parmchk failed' in line:
                errorTypes += '_5'
                ET5.append('%s_%s'% (mol, structure))
            elif 'Tleap failed' in line:
                errorTypes += '_6'
                ET6.append('%s_%s'% (mol, structure))
            elif 'syntax error' in line:
                errorTypes += '_8'
                ET8.add('%s_%s'% (mol, structure))
            elif "Use '-f' option if you want to proceed anyway. Aborting" in line:
                pass
            else:
                print "UNKNOWN ERROR:", file, line
                errorTypes += '_9'
        if "No such file or directory: 'tmp'" in line:
            errorTypes += '_7'
            ET7.add('%s_%s'% (mol, structure))
    out = parseSummurisedLine(warnTypes, errorTypes)
    if mapResults.has_key(out):
        if mapResults[out].has_key(mol):
            mapResults[out][mol].append(structure)
        else:
            mapResults[out][mol] = [structure]
    else:
        mapResults[out] = {mol:[structure]}
    return out

def parseSummurisedLine(warnTypes, errorTypes):
    wt = list(set(warnTypes.split('_')))
    if wt != ['']:
        wt = wt[1:]
        wt.sort(cmp=lambda x,y: int(x)-int(y))
    warnTypes = ''
    for i in wt:
        if i: warnTypes += '_'+i
    countWarn = warnTypes.count('_')
    et = list(set(errorTypes.split('_')))
    et.sort()
    errorTypes = ''
    for j in et:
        if j: errorTypes += '_'+j
    countError = errorTypes.count('_')
    return "%i E, %i W, ET %s, WT %s" % (countError, countWarn, errorTypes, warnTypes)

def parseChargeList(WT, dList):
    listOk = []
    for wt in WT:
        if wt.split(':')[0] in dList:
            listOk.append(wt)
    listOk.sort()
    return listOk

def myComp(vx,vy):
    ix = int(vx.split()[0])
    iy = int(vy.split()[0])
    tx = vx[-1:]
    ty = vy[-1:]
    if tx.isdigit(): x = int(tx)
    else: x = 0
    if ty.isdigit(): y = int(ty)
    else: y = 0

    if ix>iy:
            return 1
    elif ix==iy:
        if x>y:
            return 1
        elif x==y:
            return 0
        else:
            return -1
    else:
        if x>y:
            return 1
        elif x==y:
            return 0
        else:
            return -1

def sortList(d,p,i,typeMess):
    for mol, l in mapResults[typeMess].items():
        if len(l) == 2:
            d.append(mol)
        elif len(l) == 1:
            if l[0] == 'pdb': p.append(mol)
            elif l[0] == 'ideal': i.append(mol)
        else:
            print "problem with", typeMess, mol, l
    return d,p,i

def printResults(dList, pList, iList, subHead, header=None):
    global id
    dList.sort()
    pList.sort()
    iList.sort()
    id1 = id + 1
    id2 = id + 2
    print 80*'-'
    if header:
        print "\n*** For results [%i], [%i], [%i], %s:" % (id,id1,id2,header)
    print '\n[%i] %s for both PDB and Ideal:\n%i\t %s' % (id, subHead, len(dList), str(dList))
    print '\n[%i] %s with PDB ONLY, besides [%i]:\n%i\t %s' % (id1, subHead, id, len(pList), str(pList))
    print '\n[%i] %s with IDEAL ONLY, besides [%i]:\n%i\t %s' % (id2, subHead, id, len(iList), str(iList))
    pTotal = len(dList)+len(pList)
    iTotal = len(dList)+len(iList)
    total = pTotal+iTotal
    per = total / (allTotal * 0.01)
    perPdb = pTotal / (allTotal * 0.01)
    perIdeal = iTotal / (allTotal * 0.01)
    print "Total: %i of %i (%3.2f%%)" % (total, allTotal, per)
    print "PDB Total: %i of %i (%3.2f%%)" % (pTotal, allTotal, perPdb)
    print "IDEAL Total: %i of %i (%3.2f%%)" % (iTotal, allTotal, perIdeal)
    print 80*'-'
    id += 3
    return total

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

keys = mapResults.keys()
keys.sort()
keys.sort(cmp=myComp)
contador = 0
dGood, pdbGood, idealGood = [],[],[]
dGuessFail, pdbGuessFail, idealGuessFail = [],[],[]
dCloseContact, pdbCloseContact, idealCloseContact = [],[],[]
dBondIrreg, pdbBondIrreg, idealBondIrreg = [],[],[]
dBiAndAicc, pdbBiAndAicc, idealBiAndAicc = [],[],[]
dMissPar, pdbMissPar, idealMissPar = [],[],[]
for typeMess in keys:
    size = len(typeMess)
    txt = str(mapResults[typeMess])
    Nboth = txt.count("['pdb', 'ideal']") + txt.count("['ideal', 'pdb']")
    Npdb = txt.count("['pdb']")
    Nideal = txt.count("['ideal']")
    contador += 2*Nboth + Npdb + Nideal
    if typeMess in goodResults:
        dGood, pdbGood, idealGood = sortList(dGood, pdbGood, idealGood, typeMess)
    elif typeMess in guessFailedOnly:
        dGuessFail, pdbGuessFail, idealGuessFail = sortList(dGuessFail, pdbGuessFail, idealGuessFail, typeMess)
    elif typeMess in closeContact:
        dCloseContact, pdbCloseContact, idealCloseContact = sortList(dCloseContact, pdbCloseContact, idealCloseContact, typeMess)
    elif typeMess in bondIrreg:
        dBondIrreg, pdbBondIrreg, idealBondIrreg = sortList(dBondIrreg, pdbBondIrreg, idealBondIrreg, typeMess)
    elif typeMess in biAndAicc:
        dBiAndAicc, pdbBiAndAicc, idealBiAndAicc = sortList(dBiAndAicc, pdbBiAndAicc, idealBiAndAicc, typeMess)
    elif typeMess in missingPar:
        dMissPar, pdbMissPar, idealMissPar = sortList(dMissPar, pdbMissPar, idealMissPar, typeMess)
    else:
        print '*%s*%s%i %i %i' % (typeMess,(40-size)*' ', 2*Nboth, Npdb, Nideal)

allTotal = len(ccpCodes) *2

emptyDir.sort()
print '\n[A] Totally Empty Dirs (NO mol2 input files for either PDB or IDEAL):\n%i\t %s' % (len(emptyDir), str(emptyDir))
a, b, c = set(emptyDir), set(missPdbMol2), set(missIdealMol2)
diff1 = list(b.difference(a))
diff1.sort()
print '\n[B] Dirs missing pdb.mol2 input files (besides [A]):\n%i\t %s' % (len(diff1), str(diff1))
diff3 = list(c.difference(a))
diff3.sort()
print '\n[C] Dirs missing ideal.mol2 input files (besides [A]):\n%i\t %s' % (len(diff3), str(diff3))
missPdbTotal = len(missPdbMol2)
missIdealTotal = len(missIdealMol2)
missTotal = missPdbTotal + missIdealTotal
perMiss = missTotal / (allTotal * 0.01)
perPdbMiss = missPdbTotal / (allTotal * 0.01)
perIdealMiss = missIdealTotal / (allTotal * 0.01)
print "Missing Total: %i of %i (%3.2f%%)" % (missTotal, allTotal, perMiss)
print "Missing PDB Total: %i of %i (%3.2f%%)" % (missPdbTotal, allTotal, perPdbMiss)
print "Missing IDEAL Total: %i of %i (%3.2f%%)" % (missIdealTotal, allTotal, perIdealMiss)


header = "mols clean, no erros or warnings"
subHead = "Mols clean"
goodTotal = printResults(dGood, pdbGood, idealGood, subHead, header)

header = "only guessCharge failed, but running mopac with charge = 0 finished fine"
subHead = "Mols only guessCharge failed"
guessTotal = printResults(dGuessFail, pdbGuessFail, idealGuessFail, subHead, header)

header = "atoms in close contact"
subHead = "Mols have atoms in close contact"
closeTotal = printResults(dCloseContact, pdbCloseContact, idealCloseContact, subHead, header)

header = "irregular bonds "
subHead = "Mols have irregular bonds"
bondTotal = printResults(dBondIrreg, pdbBondIrreg, idealBondIrreg, subHead, header)

header = "irregular bonds and atoms in close contact"
subHead = "Mols have irregular bonds and atoms in close contact"
biAndAiccTotal = printResults(dBiAndAicc, pdbBiAndAicc, idealBiAndAicc, subHead, header)

header = "couldn't determine all parameters"
subHead = "Mols have missing parameters"
missParTotal = printResults(dMissPar, pdbMissPar, idealMissPar, subHead, header)

print "\nmissTotal + goodTotal + guessTotal + closeTotal + bondTotal + biAndAiccTotal + missParTotal =" ,\
     missTotal + goodTotal + guessTotal + closeTotal + bondTotal + biAndAiccTotal + missParTotal

#print "\n*** For results [J], [K] and [L], guessCharge failed, but running mopac with charge = 0 finished fine:"
#dBondProb.sort()
#print '\n[M] Mols bond problems for both PDB and Ideal:\n%i\t %s' % (len(dBondProb), str(dBondProb))
#pdbBondProb.sort()
#print '\n[N] Mols bond problems with PDB ONLY, besides [M]:\n%i\t %s' % (len(pdbBondProb), str(pdbBondProb))
#idealBondProb.sort()
#print '\n[O] Mols bond problems with IDEAL ONLY, besides [M]:\n%i\t %s' % (len(idealBondProb), str(idealBondProb))


#lChargeOk = parseChargeList(WT6, dGood)
#print '\n[G] Mols with charge not ZERO that are OK: %i\t %s' % (len(lChargeOk), str(lChargeOk))
#lPdbChargeOk = parseChargeList(WT6, pdbGood)
#print '\n[H] Mols with charge not ZERO that are OK for PDB ONLY, besides [G]: %i\t %s' % (len(lPdbChargeOk), str(lPdbChargeOk))
#lIdealChargeOk = parseChargeList(WT6, idealGood)
#print '\n[I] Mols with charge not ZERO that are OK for IDEAL ONLY, besides [G]: %i\t %s' % (len(lIdealChargeOk), str(lIdealChargeOk))

numDirs = len(ccpCodes)
totalPdbGood = len(pdbGood) + len(dGood)
totalIdealGood = len(idealGood) + len(dGood)
perPdb = totalPdbGood / (totalPdbMol2Count * 0.01)
perIdeal= totalIdealGood / (totalIdealMol2Count * 0.01)

#print "\n--------------------------------\nTotal molecules (Dirs): %i" % numDirs
#print "\nPDB believed to be OK: %i of %i (%3.2f%%)" % (totalPdbGood,totalPdbMol2Count, perPdb)
#print "\nIdeal believed to be OK: %i of %i (%3.2f%%)" % (totalIdealGood, totalIdealMol2Count, perIdeal)
#
#perPdb = totalPdbOkCount / (totalPdbMol2Count * 0.01)
#perIdeal= totalIdealOkCount / (totalIdealMol2Count * 0.01)
#print "\nPDB for which ACPYPI generated results: %i of %i (%3.2f%%)" % (totalPdbOkCount,totalPdbMol2Count, perPdb)
#print "\nIDEAL for which ACPYPI generated results: %i of %i (%3.2f%%)" % (totalIdealOkCount, totalIdealMol2Count, perIdeal)

print contador

