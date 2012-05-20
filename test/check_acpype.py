#!/usr/bin/env python
import sys #@UnusedImport
import os
from acpype import ACTopol
from acpype import _getoutput
import shutil
from glob import glob
#from ccpnmr.format.general.Conversion import FormatConversion

usePymol = True
useGmx45 = False
ffType = 'amber' # gaff
cType = 'bcc' #'gas'
doOpls = False
debug = False

water = ''
if useGmx45:
    water = ' -water none'

print('usePymol: %s, useGmx45: %s, ffType: %s, cType: %s, doOpls: %s' % (usePymol, useGmx45, ffType, cType, doOpls))

tmpDir = '/tmp/testAcpype'

delList = ['topol.top', 'posre.itp']

#create Dummy PDB
tpdb = '/tmp/tmp.pdb'
dummyLine = 'ATOM      1  N   ALA A   1      -1.188  -0.094   0.463  1.00  0.00           N\n'
open(tpdb, 'w').writelines(dummyLine)
tempObj = ACTopol(tpdb, chargeVal = 0, verbose = False)

# Binaries for ambertools
acExe = tempObj.acExe
tleapExe = tempObj.tleapExe
sanderExe = os.path.join(os.path.dirname(acExe), 'sander')
ambpdbExe = os.path.join(os.path.dirname(acExe), 'ambpdb')

exePymol = '/sw/bin/pymol'

# Binaries for gromacs
gpath = '/sw/'
if useGmx45:
    gpath = '/usr/local/'
gmxTopDir = gpath + "share"
pdb2gmx = gpath + '/bin/pdb2gmx'
grompp = gpath + '/bin/grompp'
mdrun = gpath + '/bin/mdrun'
g_energy = '/usr/local/bin/g_energy' #'/sw/bin/g_energy'
pdb2gmx45 = '/usr/local/bin/pdb2gmx'

amberff = "leaprc.ff99SB"

genPdbTemplate = \
'''
source %(amberff)s
mol = sequence {N%(res1)s %(res)s C%(res1)s}
savepdb mol %(aai3)s.pdb
saveamberparm mol prmtop inpcrd
quit
'''

spe_mdp = \
"""
# Create SPE.mdp file #single point energy
cat << EOF >| SPE.mdp
define = -DFLEXIBLE
integrator               = md
nsteps                   = 0
dt                       = 0.001
constraints              = none
emtol                    = 10.0
emstep                   = 0.01
nstcomm                  = 1
ns_type                  = simple
nstlist                  = 0
rlist                    = 0
rcoulomb                 = 0
rvdw                     = 0
Tcoupl                   = no
Pcoupl                   = no
gen_vel                  = no
nstxout                  = 1
pbc                      = no
nstlog = 1
nstenergy = 1
nstvout = 1
nstfout = 1
nstxtcout = 1
comm_mode = ANGULAR
continuation = yes
EOF
"""

aa_dict = {'A':'ala', 'C':'cys', 'D':'asp', 'E':'glu', 'F':'phe', 'G':'gly',
           'H':'hie', 'J':'hip', 'O':'hid', 'I':'ile', 'K':'lys', 'L':'leu',
           'M':'met', 'N':'asn', 'P':'pro', 'Q':'gln', 'R':'arg', 'S':'ser',
           'T':'thr', 'V':'val', 'W':'trp', 'Y':'tyr'}

hisDictConv = {'HHH':['HIE', 'HISB'], 'JJJ':['HIP', 'HISH'], 'OOO':['HID', 'HISA']}

if useGmx45:
    hisDictConv = {'HHH':['HIE', 'HISE'], 'JJJ':['HIP', 'HISH'], 'OOO':['HID', 'HISD']}

pymolScript = 'aa_dict = %s' % str(aa_dict) + \
'''
def build3(object_name, sequence, first_residue = "1"):
    if len(sequence):
        codeNT, code, codeCT = sequence
        if code in ['R']:
            codeNT = codeCT = 'G'
        cmd.fragment('nt_'+aa_dict[codeNT],object_name)
        cmd.alter(object_name,'resi="%s"'%first_residue)
        cmd.edit(object_name+" and name C")
        for code in sequence[1:-1]:
            editor.attach_amino_acid("pk1",aa_dict[code])
        editor.attach_amino_acid("pk1",'ct_'+aa_dict[codeCT])
        cmd.edit()

seqi=aa_dict.keys()

count=0

cmd.delete('all')

for aai in seqi:
    aai3=3*aai
    build3(aai3,aai3)
    cmd.alter(aai3,"chain='%s'"%aai3)
    cmd.translate([15*count,0,0],aai3)
    cmd.sculpt_activate(aai3)
    for n in range(100):
        cmd.sculpt_iterate(aai3)
    count=count+1
    aafile = aai3+'.pdb'
    cmd.set('pdb_conect_all')
    cmd.save(aafile,'all')
    cmd.delete('all')
'''

def pdebug(msg):
    if debug:
        print('DEBUG: %s' % msg)

def perror(msg):
    if debug:
        print('ERROR: %s' % msg)

def build_residues():
    '''Usually failing thanks to Pymol'''
    import __main__
    __main__.pymol_argv = [ 'pymol', '-Gicq'] # Quiet and no GUI

    from pymol import finish_launching
    from pymol import editor
    from pymol import cmd

    finish_launching()

    def build3(object_name, sequence, first_residue = "1"):
        if len(sequence):
            codeNT, code, codeCT = sequence
            if code in ['J', 'R']:
                codeNT = codeCT = 'G'
            cmd.fragment('nt_' + aa_dict[codeNT], object_name)
            cmd.alter(object_name, 'resi="%s"' % first_residue)
            cmd.edit(object_name + " and name C")
            for code in sequence[1:-1]:
                editor.attach_amino_acid("pk1", aa_dict[code])
            editor.attach_amino_acid("pk1", 'ct_' + aa_dict[codeCT])
            cmd.edit()

    seqi = aa_dict.keys()

    count = 0

    cmd.delete('all')

    for aai in seqi:
        aai3 = 3 * aai
        print(aai)
        build3(aai3, aai3)
        cmd.alter(aai3, "chain='%s'" % aai3)
        cmd.translate([15 * count, 0, 0], aai3)
        cmd.sculpt_activate(aai3)
        for dummy in range(100):
            cmd.sculpt_iterate(aai3)
        count = count + 1
        aafile = aai3 + '.pdb'
        cmd.set('pdb_conect_all')
        cmd.save(aafile, 'all')
        cmd.delete('all')
    #cmd.quit()
    return "PDB Residues generated"

def hisConv(fName, his, opt):
    ''' opt = 0 : pymol pdb to opls
        opt = 1 : apdb (from opdb) to amber pdb
        opt = 3 : put HIS
    '''
    hisA, hisO = hisDictConv.get(his)
    if opt == 0:
        cmd = "sed s/%s\ /%s/g %s > tMpFile; mv tMpFile %s" % (hisA, hisO, fName, fName)
    elif opt == 1:
        cmd = "sed s/%s/%s/g %s > tMpFile; mv tMpFile %s" % ('HIS', hisA, fName, fName)
    elif opt == 2:
        cmd = "sed s/%s\ /%s/g %s > tMpFile; mv tMpFile %s" % (hisA, 'HIS', fName, fName)
    pdebug(cmd)
    os.system(cmd)

def createOplsPdb(fpdb):
    '''
    pdb -> opdb (opls)
    fix atom names for OPLS in GMX
    '''
    fLines = file(fpdb).readlines()
    fName = 'o' + fpdb
    fopls = open(fName, 'w')
    for line in fLines:
        if 'ATOM  ' in line:
            if 'AAA' not in fpdb:
                line = line.replace(' 3HB ', ' 1HB ')
            line = line.replace(' 3HG ', ' 1HG ')
            line = line.replace('  GLY', '1HA  GLY')
            line = line.replace(' HA  GLY', '2HA  GLY')
            line = line.replace(' 3HD ', ' 1HD ')
            line = line.replace('3HE  LYS', '1HE  LYS')
            line = line.replace('3HG1 ILE', '1HG1 ILE')
            line = line.replace('3H   PRO', '1H   PRO')
        fopls.write(line)
    fopls.close()
    return fName

def createOplsPdb2(fpdb):
    '''
    pdb -> opdb (opls)
    fix atom names for OPLS in GMX
    special case tleap
    '''
    fLines = file(fpdb).readlines()
    fName = 'o' + fpdb
    fopls = open(fName, 'w')
    for line in fLines:
        if 'ATOM  ' in line:
            if 'AAA' not in fpdb:
                line = line.replace(' HB3 ', ' HB1 ')
            line = line.replace(' HG3 ', ' HG1 ')
            line = line.replace(' HA3 GLY', ' HA1 GLY')
            #line = line.replace(' HA  GLY', '2HA  GLY')
            line = line.replace(' HD3 ', ' HD1 ')
            line = line.replace('HE3 LYS', 'HE1 LYS')
            line = line.replace('HG13 ILE', 'HG11 ILE')
            line = line.replace('H3  PRO', 'H1  PRO')
        fopls.write(line)
    fopls.close()
    return fName

def createAmberPdb(fpdb, opt):
    '''pdb -> apdb (amber)
       opt = 0 : pymol pdb to apdb
       opt = 1 : ogpdb to apdb
    '''
    fLines = file(fpdb).readlines()
    fName = 'a' + fpdb
    famb = open(fName, 'w')
    for line in fLines:
        if 'ATOM  ' in line:
            if not useGmx45:
                line = line.replace('CYS', 'CYN')
                line = line.replace('LYS', 'LYP')
            if opt == 0:
                if 'AAA' not in fpdb:
                    if doOpls:
                        line = line.replace(' 3HB ', ' 1HB ')
                    else:
                        line = line.replace(' HB3 ', ' HB1 ')
                line = line.replace('CD1 ILE', 'CD  ILE')
                if doOpls:
                    line = line.replace(' 3HG ', ' 1HG ')
                    line = line.replace('  ', ' 1HA ')
                    line = line.replace(' HA  GLY', '2HA  GLY')
                    line = line.replace('3HG1 ILE', '1HG1 ILE')
                    line = line.replace('1HD1 ILE', ' HD1 ILE')
                    line = line.replace('2HD1 ILE', ' HD2 ILE')
                    line = line.replace('3HD1 ILE', ' HD3 ILE')
                    line = line.replace(' 3HD ', ' 1HD ')
                    line = line.replace('3HE  LYP', '1HE  LYP')
                    line = line.replace('3H   PRO', '1H   PRO')
                else:
                    line = line.replace(' HG3 ', ' HG1 ')
                    line = line.replace(' HA3 GLY', ' HA1 GLY')
                    line = line.replace('HG13 ILE', 'HG11 ILE')
                    line = line.replace('HG13 ILE', 'HG11 ILE')
                if line[22:26] == '   3':
                    line = line.replace('  O  ', '  OC1').replace('  OXT', '  OC2')
            elif opt == 1:
                if line[22:26] == '   3':
                    line = line.replace(' O1 ', ' OC1').replace(' O2 ', ' OC2')
            if line[22:26] == '   1':
                res = line[17:20]
                line = line.replace('%s ' % res, 'N%s' % res)
            if line[22:26] == '   3':
                res = line[17:20]
                line = line.replace('%s ' % res, 'C%s' % res)
        famb.write(line) ### pay attention here!
    famb.close()
    return fName

def createAmberPdb3(fpdb):
    '''Add N and C XXX --> NXXX, CXXX'''
    fLines = file(fpdb).readlines()
    fName = fpdb.replace('_', 'a')
    famb = open(fName, 'w')
    for line in fLines:
        if 'ATOM  ' in line:
            if not useGmx45:
                line = line.replace('CYS', 'CYN')
                line = line.replace('LYS', 'LYP')
            if line[22:26] == '   1':
                res = line[17:20]
                line = line.replace('%s ' % res, 'N%s' % res)
            if line[22:26] == '   3':
                res = line[17:20]
                line = line.replace('%s ' % res, 'C%s' % res)
        famb.write(line) ### pay attention here!
    famb.close()
    return fName

def createAmberPdb2(fpdb, opt):
    '''
    use formatConverter
    Add N and C XXX --> NXXX, CXXX and OC1 & OC2
    '''
    projectName = 'testImport'
    if os.path.exists(projectName):
        shutil.rmtree(projectName)
    #fc = FormatConversion(identifier = projectName, useGui = False)
    #models = fc.importFile('coordinates', 'pseudoPdb', fpdb)
    _fpdb = "_%s" % fpdb
    #fc.exportFile('coordinates', 'pseudoPdb', _fpdb, addKeywords = {'structures': models, 'forceNamingSystemName': 'PDB'})

    fLines = file(_fpdb).readlines()
    fName = 'a' + fpdb
    famb = open(fName, 'w')
    for line in fLines:
        if 'ATOM  ' in line:
            line = line.replace('CYS', 'CYN')
            line = line.replace('LYS', 'LYP')
            if opt == 0:
                if line[22:26] == '   3':
                    line = line.replace('  O   ', '  OC1 ').replace('  OXT ', '  OC2 ')
            elif opt == 1:
                if line[22:26] == '   3':
                    line = line.replace(' O1 ', ' OC1').replace(' O2 ', ' OC2')
            if line[22:26] == '   1':
                res = line[17:20]
                line = line.replace('%s ' % res, 'N%s' % res)
            if line[22:26] == '   3':
                res = line[17:20]
                line = line.replace('%s ' % res, 'C%s' % res)
        famb.write(line) ### pay attention here!
    famb.close()
    return fName

def createOldPdb(fpdb):
    '''using pdb2gmx 4.5 -ignh to get the atom names and in file _fpdb
    '''
    os.system('rm -fr \#*')
    fName = '_' + fpdb
    cmd = "%s -f %s -o %s -ff amber99sb -water none -ignh >& /dev/null" % (pdb2gmx45, fpdb, fName)
    out = os.system(cmd)
    atList = []
    if out == 0:
        for item in delList:
            os.remove(item)
        _fLines = file(fName).readlines()
        for line in _fLines:
            if 'ATOM  ' in line:
                atList.append(line[11:17])
    else:
        perror(out)
    fLines = file(fpdb).readlines()
    nLines = []
    n = 0
    for line in fLines:
        if 'ATOM  ' in line:
            nLines.append(line[:11] + atList[n] + line[17:])
            n += 1
    file(fName, 'w').writelines(nLines)
    os.system('rm -fr \#*')
    return fName

def createOldPdb2(fpdb):
    '''using my own dict for e.g. 2HB = HB2 -> HB1'''
    defHB1 = [' HB2', '2HB ', ' HB1']
    defHB2 = [' HB3', '3HB ', ' HB2']
    defHG1 = [' HG2', '2HG ', ' HG1']
    defHG2 = [' HG3', '3HG ', ' HG2']
    defHD1 = [' HD2', '2HD ', ' HD1']
    defHD2 = [' HD3', '3HD ', ' HD2']
    def1HG2 = ['HG21', '1HG2']
    def2HG2 = ['HG22', '2HG2']
    def3HG2 = ['HG23', '3HG2']

    dictPdb2GmxAtomNames = {'ALA':(), 'VAL':()
            , 'CYS':(defHB1, defHB2)
            , 'ASP':(defHB1, defHB2)
            , 'GLU':(defHB1, defHB2, defHG1, defHG2)
            , 'PHE':(defHB1, defHB2)
            , 'GLY':([' HA2 ', ' HA  ', ' HA1 '], [' HA3', '3HA ', ' HA2'])
            , 'HIE':(defHB1, defHB2)
            , 'ILE':(['CD1', 'CD '], ['HG12', '2HG1', '1HG1'], ['HG13', '3HG1', '2HG1'], def1HG2, def2HG2, def3HG2, ['HD11', '1HD1', ' HD1'], ['HD12', '2HD1', ' HD2'], ['HD13', '3HD1', ' HD3'])
            , 'HIP':(defHB1, defHB2)
            , 'HID':(defHB1, defHB2)
            , 'LYS':(defHB1, defHB2, defHG1, defHG2, defHD1, defHD2, [' HE2', '2HE ', ' HE1'], [' HE3', '3HE ', ' HE2'])
            , 'LEU':(defHB1, defHB2)
            , 'MET':(defHB1, defHB2, defHG1, defHG2)
            , 'ASN':(defHB1, defHB2)
            , 'PRO':(defHB1, defHB2, defHG1, defHG2, defHD1, defHD2, [' H3 ', '3H  ', ' H1 '])
            , 'GLN':(defHB1, defHB2, defHG1, defHG2, ['HE21', '1HE2'], ['HE22', '2HE2'])
            , 'ARG':(defHB1, defHB2, defHG1, defHG2, defHD1, defHD2)
            , 'SER':(defHB1, defHB2)
            , 'THR':(def1HG2, def2HG2, def3HG2)
            , 'TRP':(defHB1, defHB2)
            , 'TYR':(defHB1, defHB2)
            }
    fName = '_' + fpdb
    fLines = file(fpdb).readlines()
    aLines = []
    for line in fLines:
        if 'ATOM  ' in line:
            aLines.append(line)
    nRes = int(aLines[-1][22:26])
    nLines = []
    for line in aLines:
        res = line[17:20]
        nResCur = int(line[22:26])
        for item in dictPdb2GmxAtomNames.get(res):
            #line = line.replace('%s  %s' % (item[0], res), '%s  %s' % (item[2], res)).replace('%s  %s' % (item[1], res), '%s  %s' % (item[2], res))
            atAim = item[-1]
            for at in item[:-1]:
                #print(at, atAim)
                #print(line)
                line = line.replace('%s' % at, '%s' % atAim)
                #print(line)
        if nResCur == nRes:
            line = line.replace('O   %s' % res, 'OC1 %s' % res)
            line = line.replace('OXT %s' % res, 'OC2 %s' % res)
        nLines.append(line)
    file(fName, 'w').writelines(nLines)
    return fName

def parseTopFile(lista):
    ''' lista = top file in list format
        return a (dict,dict) with fields'''
    defDict = {}
    parDict = {}
    for line in lista:
        line = line.split(';')[0]
        line = line.strip()
        if line.startswith('#def'):
            vals = line.split()
            defDict[vals[1]] = map(eval, vals[2:])
        elif line and not line.startswith('#'):
            if line.startswith('[ '):
                flag = line.split()[1]
                if not parDict.has_key(flag):
                    parDict[flag] = []
            else:
                u = []
                t = line.split()
                for i in t:
                    try: v = eval(i)
                    except: v = i
                    u.append(v)
                parDict[flag].append(u)
    return parDict, defDict

def nbDict(lista):
    tdict = {}
    for line in lista:
        line = line.split(';')[0]
        line = line.strip()
        if line and not line.startswith('#') and not line.startswith('['):
            name, type = line.split()[0:2]
            tdict[name] = type
    return tdict

def checkTopOplsVsAmber(res):

    agRes = parseTopFile(file('ag%s.top' % res).readlines())
    ogRes = parseTopFile(file('og%s.top' % res).readlines())
#    mRes = parseTopFile(file('g%s.acpype/g%s_GMX.itp' % (res,res)).readlines())

#    flags = ['pairs', 'bonds', 'angles', 'dihedrals'], ['dihedraltypes', 'angletypes', 'bondtypes']

    agPairs = agRes[0]['pairs']
    ogPairs = ogRes[0]['pairs']
#    mPairs = mRes[0]['pairs']

    if agPairs != ogPairs:
        print("!!!!!!", res)

def roundAllFloats(lista, l):
    """Round to 3 decimals"""
    nlista = []
    for ii in lista:
        tt = ii[:l + 1]
        for jj in ii[l + 1:]:
            if jj > 100.0:
                jj = round(jj, -1)
            nn = round(jj, 3)
            tt.append(nn)
        nlista.append(tt)
    return nlista

def checkTopAcpype(res, ff):
    '''Compare acpype gmx itp against amber or opls pdb2gmx results'''
    global amber2oplsDict, ac2opls
    os.chdir(tmpDir)
    def addParam(l, item):
        ''' l : index for bonded types
            item : set of atomtypes'''
        dict = {2:'bondtypes', 3:'angletypes', 4:'dihedraltypes'}
        dType = {} # dict
        for type in ffBon[0][dict[l]]:
            i = type[:l + 1]
            j = type[l + 1:]
            dType[str(i)] = j # dict {[atomtypes,funct] : parameters}
        entries = []
        lNum = item[l] # funct
        ent = [ffDictAtom[x] for x in item[:l]] # convert atomtypes ids to atomtypes names
        rent = ent[:]
        rent.reverse()
        entries.append(ent + [lNum])
        entries.append(rent + [lNum])
        if l == 4:
            if len(item) == 6:
                par = ffBon[1][item[5]]
                return par
            tent = ent[:]
            rtent = rent[:]
            if lNum in [3, 9]: # dih proper
                tent[0] = 'X'
                entries.append(tent + [lNum])
                tent.reverse()
                entries.append(tent + [lNum])
                tent[0] = 'X'
                entries.append(tent + [lNum])
                tent.reverse()
                entries.append(tent + [lNum])
                rtent[0] = 'X'
                entries.append(rtent + [lNum])
                rtent.reverse()
                entries.append(rtent + [lNum])
            if lNum in [1, 4]: # dih improp
                tent[0] = 'X'
                entries.append(tent + [lNum])
                tent[1] = 'X'
                entries.append(tent + [lNum])
                rtent[0] = 'X'
                entries.append(rtent + [lNum])
                rtent[1] = 'X'
                entries.append(rtent + [lNum])
        found = False
        for e in entries:
            try:
                par = dType[str(e)]
                found = True
                break
            except: pass
        if not found:
            print("%s %s %s not found in %s Bon" % (dict[l], ent, item, ffType))
        #print(item, e, par)
        return par

    compareParameters = False
    if ff == 'a' : compareParameters = True

    agRes = parseTopFile(file('ag%s.top' % (res)).readlines())
    if doOpls:
        ogRes = parseTopFile(file('og%s.top' % (res)).readlines())
    #print(ogRes)
    #ffgRes = parseTopFile(file('%sg%s.top' % (ff,res)).readlines())
    #print 111, os.getcwd()
    acRes = parseTopFile(file('%sg%s.acpype/%sg%s_GMX.itp' % (ff, res, ff, res)).readlines())
    if ff == 'a':
        ffNb = aNb
        ffBon = aBon
        ffgRes = agRes
    if doOpls:
        if ff == 'o':
            ffNb = oNb
            ffBon = oBon
            ffgRes = ogRes

    atError = False
    if ff == 'a':
        print("    ==> Comparing atomtypes AC x AMBER")

    ffDictAtom = {} # dict link res atom numbers to amber or opls atom types
    for item in ffgRes[0]['atoms']:
        i, j = item[:2] # id, at
        ambat = ffNb[j]
        if useGmx45:
            ambat = j
        ffDictAtom[i] = ambat
        # compare atom types AC x Amber
        if ff == 'a':
            acat = acRes[0]['atoms'][i - 1][1]
            if ambat != acat:
                print("    %i %s AC: %s   x   AMB: %s" % (i, item[4], acat, ambat))
                atError = True

    if ff == 'a' and not atError:
        print("        atomtypes OK")

    acDictAtom = {}
    for item in acRes[0]['atoms']:
        i, j = item[:2]
        acDictAtom[i] = j # dict for atom id -> atom type from acpype itp file

    for k in acDictAtom.keys():
        if amber2oplsDict.has_key(acDictAtom[k]):
            amber2oplsDict[acDictAtom[k]].add(ffDictAtom[k])
        else:
            amber2oplsDict[acDictAtom[k]] = set([ffDictAtom[k]])

    # to build a dict acpype atom types (amber or gaff) to [[opls_code],[opls_mass]]
    #ac2opls
    if doOpls:
        if ff == 'o':
            for item in acRes[0]['atoms']:
                id, at = item[:2]
                it = ogRes[0]['atoms'][id - 1]
                ocode, omass = it[1], it[-1]
                #print item, it, ocode, omass
                if ac2opls.has_key(at):
                    ac2opls[at][0].add(ocode)
                    ac2opls[at][1].add(omass)
                else:
                    ac2opls[at] = [set([ocode]), set([omass])]

    if doOpls:
        # to build dict a2oD: amber atom type -> opls atom type
        idx = 0
        for aAT in agRes[0]['atoms']:
            oAT = ogRes[0]['atoms'][idx]
            aAT, aPdb = aAT[1], aAT[4]
            oAT, oPdb = oAT[1], oAT[4]
            if aPdb == oPdb:
                if a2oD.has_key(aAT):
                    a2oD[aAT].add(oAT)
                else:
                    a2oD[aAT] = set([oAT])
            else:
                #print agRes[0]['atoms'][idx], ogRes[0]['atoms'][idx]
                pass
            idx += 1

    flags = [('pairs', 2), ('bonds', 2), ('angles', 3), ('dihedrals', 4)] #, ['dihedraltypes', 'angletypes', 'bondtypes']

    for flag, l in flags:
        print("    ==> Comparing %s" % flag)
        if flag != flags[0][0] and compareParameters: # not 'pairs'
            agres = []
            tAgRes = ffgRes[0][flag] # e.g. dic 'bonds' from gmx top file
            for item in tAgRes:
                if flag == flags[1][0]: # 'bonds'
                    par = addParam(l, item)
                elif flag == flags[2][0]: # 'angles'
                    par = addParam(l, item)
                elif flag == flags[3][0]: # 'dihedrals', e.g. item = [2, 1, 5, 6, 9]
                    par = addParam(l, item)
                    if len(item) == 6:
                        item.pop()
                agres.append(item + par)
            if compareParameters: acres = acRes[0][flag]
        else:
            if flag == flags[3][0]:
                agres = [x[:l + 1] for x in ffgRes[0][flag]]
            else: agres = ffgRes[0][flag]
            if compareParameters: acres = acRes[0][flag]
            else: acres = [x[:l + 1] for x in acRes[0][flag]]

        act = acres[:]
        agt = agres[:]

        if compareParameters:
            if flag != flags[0][0]:
                # round parameters
                act = roundAllFloats(act, l)
                agt = roundAllFloats(agt, l)
                act2 = act[:]
                agt2 = agt[:]

        if not compareParameters or flag == flags[0][0]:
            act2 = act[:]
            agt2 = agt[:]

        for ac in act:
            if str(ac) in str(agt):
                act2.remove(ac)
                agt2.remove(ac)
            else:
                t = ac[:-1]
                t.reverse()
                acr = t + [ac[-1]]
                if str(acr) in str(agt):
                    act2.remove(ac)
                    agt2.remove(acr)

        act3 = act2[:]
        agt3 = agt2[:]

        if act2 and agt2:
            # specially for dih since it may need to resort indexes
            agl = {}
            for ag in agt3:
                t = ag[:-1]
                t.sort()
                ags = t + [ag[-1]]
                agl[str(ags)] = ag
            for ac in act2:
                t = ac[:-1]
                t.sort()
                acs = t + [ac[-1]]
                if str(acs) in str(agl.keys()):
                    act3.remove(ac)
                    agt3.remove(agl[str(acs)])

        if act3:
            act3.sort()
            print("    ac: ", act3)
        if agt3:
            agt3.sort()
            print("    %sg: " % ff, agt3)

        if not act3 and not agt3:
            print("        %s OK" % flag)

def parseOut(out):
    lines = out.splitlines()
    #for line in lines:
    count = 0
    while count < len(lines):
        line = lines[count]
        if 'WARNING' in line.upper():
            if 'will be determined based' in line:
                pass
            elif 'Problems reading a PDB file' in line:
                pass
            elif 'Open Babel Warning' in line:
                pass
            elif 'no charge value given' in line:
                pass
            elif 'audit log messages' in line:
                pass
            elif 'all CONECT' in line:
                pass
            elif 'with zero occupancy' in line:
                pass
            elif 'check:  Warnings:' in line:
                pass
            else:
                print(line)
        elif 'Total charge' in line:
            print("    ", line)
        elif 'ERROR' in line.upper():
            if 'Fatal error' in line:
                print(line, lines[count + 1])
            else:
                print(line)
        count += 1

def fixRes4Acpype(fpdb):
    code = fpdb[2]
    fLines = file(fpdb).readlines()
    famb = open(fpdb, 'w')
    for line in fLines:
        if 'ATOM  ' in line:
            line = line[:17] + aa_dict[code].upper() + line[20:]
        famb.write(line)
    famb.close()

def build_residues_tleap():
    """Build residues tripeptides with tleap and minimise with sander"""
    mdin = '''Minimization\n&cntrl\nimin=1, maxcyc=200, ntmin=2, ntb=0, igb=0,cut=999,/\n'''
    open('mdin', 'w').writelines(mdin)
    seqi = aa_dict.keys()
    for aai in seqi:
        #if aai != 'H': continue
        aai3 = 3 * aai
        res = aa_dict.get(aai).upper()
        if res in ['ARG']:
            res1 = 'GLY'
        else:
            res1 = res
        leapDict = {'amberff' : amberff, 'res' : res, 'aai3' : aai3, 'res1':res1}
        tleapin = genPdbTemplate % leapDict
        open('tleap.in', 'w').writelines(tleapin)
        cmd = "%s -f tleap.in" % (tleapExe)
        _getoutput(cmd)
        #cmd = "%s -O; %s < restrt > %s.pdb; mv mdinfo %s.mdinfo" % (sanderExe, ambpdbExe, aai3, aai3) # -i mdin -o mdout -p prmtop -c inpcrd" % (sanderExe)
        cmd = "%s -O; %s < restrt > %s.pdb" % (sanderExe, ambpdbExe, aai3)
        _getoutput(cmd)
    _getoutput('rm -f mdout mdinfo mdin restrt tleap.in prmtop inpcrd leap.log')

def error(v1, v2):
    '''percentage relative error'''
    return abs(v1 - v2) / max(abs(v1), abs(v2)) * 100

def calcGmxPotEnerDiff(res):
    def getEnergies(template):
        cmd = template % cmdDict
        pdebug(cmd)
        out = _getoutput(cmd)
        out = out.split('\n')[-nEner:]
        dictEner = {}
        for item in out:
            k, v = item.replace(' Dih.', '_Dih.').replace(' (SR)', '_(SR)').replace('c En.', 'c_En.').replace(' (bar)', '_(bar)').replace('l E', 'l_E').split()[:2]
            v = eval(v)
            dictEner[k] = v
            if 'Dih.' in k or 'Bell.' in k:
                if 'Dihedral' in list(dictEner.keys()):
                    dictEner['Dihedral'] = dictEner['Dihedral'] + v
                else:
                    dictEner['Dihedral'] = v
            if 'Coulomb' in k or 'LJ' in k:
                if 'Non_Bonded' in list(dictEner.keys()):
                    dictEner['Non_Bonded'] = dictEner['Non_Bonded'] + v
                else:
                    dictEner['Non_Bonded'] = v
        dictEner['Bonded'] = dictEner.get('Potential') - dictEner.get('Non_Bonded')
        return dictEner

    os.chdir(tmpDir)
    nEner = 9 # number of energy entries from g_energy
    tEner = ' '.join([str(x) for x in range(1, nEner + 1)])
    open('SPE.mdp', 'w').writelines(spe_mdp)
    prefix = 'a'
    if doOpls:
        prefix = 'aog'
    cmdDict = {'pdb2gmx':pdb2gmx, 'grompp':grompp, 'mdrun':mdrun, 'res':res,
               'g_energy':g_energy, 'tEner':tEner, 'water':water, 'prefix':prefix}

    # calc Pot Ener for agaXXX.pdb or aogXXX.pdb
#    if not useGmx45:
#        shutil.copy('ag%s.pdb' % res, 'ag_%s.pdb' % res)
#        createAmberPdb3('ag_%s.pdb' % res)
#    else:
#        shutil.copy('ag%s.pdb' % res, 'aga%s.pdb' % res)
    template = '''%(pdb2gmx)s -ff amber99sb -f %(prefix)s%(res)s.pdb -o %(prefix)s%(res)s_.pdb -p %(prefix)s%(res)s.top %(water)s
    %(grompp)s -c %(prefix)s%(res)s_.pdb -p %(prefix)s%(res)s.top -f SPE.mdp -o %(prefix)s%(res)s.tpr
    %(mdrun)s -v -deffnm %(prefix)s%(res)s
    echo %(tEner)s | %(g_energy)s -f %(prefix)s%(res)s.edr
    '''
    dictEner1 = getEnergies(template)
    #print(dictEner1)

    #calc Pot Ener for agXXX.acpype/agXXX.pdb
    template = '''%(grompp)s -c ag%(res)s.acpype/ag%(res)s_NEW.pdb -p ag%(res)s.acpype/ag%(res)s_GMX.top -f SPE.mdp -o ag%(res)s.tpr
    %(mdrun)s -v -deffnm ag%(res)s
    echo %(tEner)s | %(g_energy)s -f ag%(res)s.edr
    '''
    dictEner2 = getEnergies(template)
    #print(dictEner2)

    for k, v in dictEner1.items():
        v2 = dictEner2.get(k)
        rerror = error(v2, v)
        if rerror > 0.1:
            print("%15s   %.3f   %5.6f x %5.6f" % (k, rerror, v2, v))
        else:
            print("%15s   %.3f" % (k, rerror))

if __name__ == '__main__':

    '''order: (tleap/EM or pymol) AAA.pdb -f-> oAAA.pdb --> (pdb2gmx) ogAAA.pdb -f->
         -f-> aogAAA.pdb --> (pdb2gmx) agAAA.pdb -f-> acpype
    '''
    os.system('rm -fr \#*')
    global amber2oplsDict, a2oD, ac2opls
    amber2oplsDict = {}
    a2oD = {}
    ac2opls = {}
    if useGmx45:
        aNb = nbDict(file(gmxTopDir + '/gromacs/top/amber99sb.ff/ffnonbonded.itp').readlines())
        oNb = nbDict(file(gmxTopDir + '/gromacs/top/oplsaa.ff/ffnonbonded.itp').readlines())
        aBon = parseTopFile(file(gmxTopDir + '/gromacs/top/amber99sb.ff/ffbonded.itp').readlines())
        oBon = parseTopFile(file(gmxTopDir + '/gromacs/top/oplsaa.ff/ffbonded.itp').readlines())
    else:
        aNb = nbDict(file(gmxTopDir + '/gromacs/top/ffamber99sbnb.itp').readlines())
        oNb = nbDict(file(gmxTopDir + '/gromacs/top/ffoplsaanb.itp').readlines())
        aBon = parseTopFile(file(gmxTopDir + '/gromacs/top/ffamber99sbbon.itp').readlines())
        oBon = parseTopFile(file(gmxTopDir + '/gromacs/top/ffoplsaabon.itp').readlines())

    tmpFile = 'tempScript.py'
    if not os.path.exists(tmpDir):
        os.mkdir(tmpDir)
    os.chdir(tmpDir)
    #create res.pdb
    if usePymol:
        ff = open(tmpFile, 'w')
        ff.writelines(pymolScript)
        ff.close()
        cmd = "%s -qc %s" % (exePymol, tmpFile)
        os.system(cmd)
    else:
        build_residues_tleap()
    listRes = os.listdir('.') # list all res.pdb
    listRes = glob('???.pdb')
    listRes.sort()
    for resFile in listRes:
        res, ext = os.path.splitext(resFile) # eg. res = 'AAA'
        if len(resFile) == 7 and ext == ".pdb" and resFile[:3].isupper():
            #if not resFile == 'HHH.pdb': continue
            print("\nFile %s : residue %s" % (resFile, aa_dict[res[0]].upper()))

            if doOpls:
                if usePymol:
                    # from pdb pymol to opls pdb and top
                    opdb = createOplsPdb(resFile)
                else:
                    # from pdb org to opls pdb and top
                    opdb = createOplsPdb2(resFile)
                ogpdb = 'og%s' % resFile
                ogtop = 'og%s.top' % res

                if res in hisDictConv.keys():
                    hisConv(opdb, res, 0) # fix HIE etc. from pymol/tleap pdb to pdb opls
                #cmd = ' % s - f % s - o % s - p % s - ff oplsaa - water none' % (pdb2gmx, opdb, ogpdb, ogtop) # coords for O term changed# -water none for gmx 4.5
                cmd = ' %s -f %s -o %s -p %s -ff oplsaa %s' % (pdb2gmx, opdb, ogpdb, ogtop, water)
                #print(1111, cmd)
                #out = commands.getstatusoutput(cmd)
                out = _getoutput(cmd)
                #parseOut(out)
                apdb = createAmberPdb(ogpdb, 1) # create file aogAAA.pdb
            else:
                _pdb = createOldPdb2(resFile) # using my own dict
                #_pdb = createOldPdb(resFile) # using pdb2gmx 4.5
                #if res in hisDictConv.keys():
                #    hisConv(resFile, res, 2)
                apdb = createAmberPdb3(_pdb) # create file aAAA.pdb with NXXX, CXXX

            # from ogpdb to amber pdb and top
            agpdb = 'ag%s.pdb' % res # output name
            agtop = 'ag%s.top' % res
            #if res in hisDictConv.keys():
                    #hisConv(apdb, res, 2)
            #        hisFlag = ''
            cmd = ' %s -f %s -o %s -p %s -ff amber99sb %s' % (pdb2gmx, apdb, agpdb, agtop, water)
            pdebug(cmd)
            out = _getoutput(cmd)
            #print(out)
            #parseOut(out)

            # acpype on agpdb file
            fixRes4Acpype(agpdb)
            #if doOpls:
            #    fixRes4Acpype(agpdb)
            #else:
            #    shutil.copy(resFile, agpdb)
            cv = None
            #cmd = "%s -dfi %s -c %s -a %s" % (acpypeExe, agpdb, cType, ffType)
            if res in ['JJJ', 'RRR'] and cType == 'bcc': cv = 3 #cmd += ' -n 3' # acpype failed to get correct charge
            mol = ACTopol(agpdb, chargeType = cType, atomType = ffType, chargeVal = cv,
                          debug = False, verbose = debug, gmx45 = useGmx45)
            mol.createACTopol()
            mol.createMolTopol()

            #out = commands.getstatusoutput(cmd)
            #parseOut(out[1])
            print("Compare ACPYPE x GMX AMBER99SB topol & param")
            checkTopAcpype(res, 'a')
            # calc Pot Energy
            print("Compare ACPYPE x GMX AMBER99SB Pot. Energy (|ERROR|%)")
            # calc gmx amber energies
            calcGmxPotEnerDiff(res)
            # calc gmx acpype energies

            # acpype on ogpdb file
            if doOpls:
                fixRes4Acpype(ogpdb)
                #cmd = "%s -dfi %s -c %s -a %s" % (acpypeExe, ogpdb, cType, ffType)
                if res == 'JJJ' and cType == 'bcc': cmd += ' -n 3' # acpype failed to get correct charge
                mol = ACTopol(ogpdb, chargeType = cType, atomType = ffType, verbose = False)
                mol.createACTopol()
                mol.createMolTopol()

                #out = commands.getstatusoutput(cmd)
                #parseOut(out[1])
                print("Compare ACPYPE x GMX OPLS/AA only topol")
                checkTopAcpype(res, 'o')

#            if out[0] :
#                print("pdb2gmx for %s FAILED... aborting." % res)
#                print(cmd)
#            #break
    os.system('rm -f %s/\#* posre.itp tempScript.py' % tmpDir)
    os.system("find . -name 'ag*GMX*.itp' | xargs grep -v 'created by acpype on' > standard_ag_itp.txt")
    os.system("find . -name 'og*GMX*.itp' | xargs grep -v 'created by acpype on' > standard_og_itp.txt")

    #print ac2opls
    if doOpls:
        tmp = {}
        for k, v in ac2opls.items():
            vOpls, vMass = v
            vOpls = [x.replace('opls_', '') for x in vOpls]
            vMass = map(str, vMass)
            vOpls.sort()
            vMass.sort()
            if len(vMass) > 1:
                print(k, v)
            tmp[k] = vOpls + vMass
        #print tmp

        #otopFiles = commands.getoutput('ls o*.top').splitlines()
    #    print amber2oplsDict
    #    print a2oD
    #    ll = []
    #    for k in amber2oplsDict:
    #        v = amber2oplsDict[k]
    #        if len(v) == 1:
    #            ll.append([1,k,list(v)])
    #        elif len(v) == 2:
    #            ll.append([2,k,list(v)])
    #        elif len(v) == 3:
    #            ll.append([3,k,list(v)])
    #        elif len(v) == 4:
    #            ll.append([4,k,list(v)])
    #
    #    aNbRev = {}
    #    for k,v in aNb.items():
    #        if aNbRev.has_key(v):
    #            aNbRev[v].append(k)
    #        else:
    #            aNbRev[v] = [k]
    #    for k,v in aNbRev.items():
    #        v.sort()
    #        v0 = v[0]
    #        aNbRev[k] = v0
    #    from operator import itemgetter
    #    print sorted(aNbRev.items(), key=itemgetter(1))
    #
    #    o1 = file('opls_sorted.txt').readlines()
    #    dictOplsAtomType2OplsGmxCode = {}
    #    dictOplsMass = {}
    #    for line in o1:
    #        line = line[:-1]
    #        atOpls, gmxCodeOpls, mass = line.split()
    #        if dictOplsAtomType2OplsGmxCode.has_key(atOpls):
    #            dictOplsAtomType2OplsGmxCode[atOpls].append(gmxCodeOpls)
    #            dictOplsMass[atOpls].add(mass)
    #        else:
    #            dictOplsAtomType2OplsGmxCode[atOpls] = [gmxCodeOpls]
    #            dictOplsMass[atOpls] = set([mass])
    #
    #
    #    notOpls = []
    #    for ambKey, ambVal in dictAmbAtomType2AmbGmxCode.items():
    #        if dictOplsAtomType2OplsGmxCode.has_key(ambKey):
    #            print ambKey,ambVal,dictOplsAtomType2OplsGmxCode[ambKey]
    #        else:
    #            notOpls.append(ambKey)
    #    print "No opls AT for %s" % str(notOpls)
