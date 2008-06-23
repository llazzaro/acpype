#!/usr/bin/env python

from commands import getoutput
from shutil import copy2
import os
import cPickle as pickle
import sys

# input parameters: file.pdb [file.mol2] [-c gas|bcc] [-nc 0]

#call antechamber
USAGE = \
"""
    acpypi -i _file_ [-c _string_] [-n _int_] [-m _int_] [-a _string_] [-f]
    -i    input file name with either extension '.pdb' or '.mol2' (mandatory)
    -c    charge method: gas, bcc (default)
    -n    net molecular charge (int), for gas default is 0
    -m    multiplicity (2S+1), default is 1
    -a    atom type, can be gaff, amber, bcc and sybyl, default is gaff
    -f    force recalculate topologies
"""

SLEAP_TEMPLATE = \
"""
source leaprc.ff03
source leaprc.gaff
set default fastbld on
%(base)s = loadpdb %(base)s.pdb
saveamberparm %(base)s %(base)s.top %(base)s.xyz
quit
"""

TLEAP_TEMPLATE = \
"""
verbosity 1
source leaprc.gaff
mods = loadamberparams frcmod
%(base)s = loadmol2 %(acMol2File)s
saveamberparm %(base)s %(base)s.top %(base)s.xyz
quit
"""

def invalidArgs():

    print USAGE
    sys.exit(1)

def parseArgs(args):

    import getopt

    options = 'hi:c:n:m:a:f'

    ctList = [ 'gas', 'bcc']
    atList = [ 'gaff', 'amber', 'bcc', 'sybyl']

    try:
        opt_list, args = getopt.getopt(args, options) #, long_options)
    except:
        invalidArgs()

    d = {}

    for key, value in opt_list:
        if key in d:
            invalidArgs()

        if value == '':
            value = None

        d[key] = value

    if not d and not args:
        invalidArgs()

    not_none = ('-i', '-c', '-n', '-m','-a')

    for option in not_none:
        if option in d:
            if d[option] is None:
                invalidArgs()

    if '-c' in d.keys():
        if d['-c'] not in ctList:
            invalidArgs()

    if '-a' in d.keys():
        if d['-a'] not in atList:
            invalidArgs()

    if '-h' in d.keys():
        invalidArgs()

    if not '-i' in d.keys():
        invalidArgs()

    if not os.path.exists(d['-i']):
        print "ERROR: input file '%s' doesn't exist" % d['-i']
        invalidArgs()

    if args:
        invalidArgs()

    return d

class ACTopol:
    """
        Class to build the AC topologies (Antechamber AmberTools 1.0)
    """

    def __init__(self, inputFile, chargeType = 'bcc', chargeVal = None,
                 multiplicity = '1', atomType = 'gaff', force = False):

        self.inputFile = inputFile
        self.rootDir = os.path.abspath('.')
        self.absInputFile = os.path.abspath(inputFile)
        if not os.path.exists(self.absInputFile):
            print "WARNING: input file doesn't exist"
        base, ext = os.path.splitext(inputFile)
        self.baseName = base
        self.ext = ext
        self.homeDir = self.baseName + '.acpypi'
        self.chargeType = chargeType
        self.chargeVal = chargeVal
        self.multiplicity = multiplicity
        self.atomType = atomType
        self.force = force
        self.acExe = getoutput('which antechamber') or None
        self.tleapExe = getoutput('which tleap') or None
        self.sleapExe = getoutput('which sleap') or None
        self.parmchkExe = getoutput('which parmchk') or None
        self.babelExe = getoutput('which babel') or None
        self.acXyz = None
        self.acTop = None
        if self.chargeVal == None:
            print "WARNING: no charge value given, trying to guess one..."
            if self.guessCharge():
                print "ERROR: failed to guess charge"
        acMol2File = '%s_%s_%s.mol2' % (base, chargeType, atomType)
        self.acMol2File = acMol2File
        self.acParDict = {'base':base,'ext':ext[1:], 'acMol2File':acMol2File}

    def guessCharge(self):
        """
            Guess the charge of a system based on antechamber
            Returns None i case of error
        """
        charge = 0
        localDir = os.path.abspath('.')
        tmpDir = '/tmp/acpypi'
        if not os.path.exists(tmpDir):
            os.mkdir(tmpDir)
        if not os.path.exists(os.path.join(tmpDir, self.inputFile)):
            copy2(self.absInputFile, tmpDir)
        os.chdir(tmpDir)
        if self.ext == ".pdb":
            cmd = '%s -ipdb %s -omol2 %s.mol2' % (self.babelExe, self.inputFile,
                                                  self.baseName)
            _out = getoutput(cmd)
            # print _out
        cmd = '%s -i %s.mol2 -fi mol2 -o tmp.mol2 -fo mol2 -c gas' % \
                                                    (self.acExe, self.baseName)
        log = getoutput(cmd)
        #print log
        m = log.split()
        if len(m) >= 21:
            charge = int(float(m[14].replace('(','').replace(')','')))
        else:
            print "ERROR: guessCharge failed"
            os.chdir(localDir)
            return True
        self.chargeVal = str(charge)
        print "... charge set to", charge
        os.chdir(localDir)

    def execAntechamber(self, chargeType = None, atomType = None):
        """
            To call Antechamber and execute it

Usage: antechamber -i   input file name
                   -fi  input file format
                   -o   output file name
                   -fo  output file format
                   -c   charge method
                   -cf  charge file name
                   -nc  net molecular charge (int)
                   -a   additional file name
                   -fa  additional file format
                   -ao  additional file operation
                        crd : only read in coordinate
                        crg: only read in charge
                        name  : only read in atom name
                        type  : only read in atom type
                        bond  : only read in bond type
                   -m   multiplicity (2S+1), default is 1
                   -rn  residue name, if not available in the input file, default is MOL
                   -rf  residue toplogy file name in prep input file, default is molecule.res
                   -ch  check file name in gaussian input file, default is molecule
                   -ek  empirical calculation (mopac or divcon) keyword in a pair of quotation marks
                   -gk  gaussian keyword in a pair of quotation marks
                   -df  use divcon flag, 1 - use divcon; 0 - use mopac (default)
                   -at  atom type, can be gaff, amber, bcc and sybyl, default is gaff
                   -du  check atom name duplications, can be yes(y) or no(n), default is yes
                   -j   atom type and bond type prediction index, default is 4
                        0    : no assignment
                        1    : atom type
                        2    : full  bond types
                        3    : part  bond types
                        4    : atom and full bond type
                        5    : atom and part bond type
                   -s   status information, can be 0 (brief), 1 (the default) and 2 (verbose)
                   -pf  remove the intermediate files: can be yes (y) and no (n), default is no
                   -i -o -fi and -fo must appear in command lines and the others are optional

                             List of the File Formats

                file format type  abbre. index | file format type abbre. index
                ---------------------------------------------------------------
                Antechamber        ac       1  | Sybyl Mol2         mol2    2
                PDB                pdb      3  | Modified PDB       mpdb    4
                AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6
                Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8
                Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10
                Gaussian Output    gout    11  | Mopac Output       mopout 12
                Alchemy            alc     13  | CSD                csd    14
                MDL                mdl     15  | Hyper              hin    16
                AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18
                Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20
                Divcon Input       divcrt  21  | Divcon Output      divout 22
                Charmm             charmm  23
                --------------------------------------------------------------

                AMBER restart file can only be read in as additional file.

                             List of the Charge Methods

                charge method     abbre.  index | charge method      abbre. index
                ----------------------------------------------------------------
                RESP               resp     1  |  AM1-BCC            bcc     2
                CM1                cm1      3  |  CM2                cm2     4
                ESP (Kollman)      esp      5  |  Mulliken           mul     6
                Gasteiger          gas      7  |  Read in charge     rc      8
                Write out charge   wc       9  |  Delete Charge      dc     10
                ----------------------------------------------------------------
        """

        print "Executing Antechamber..."

        self.makeDir()

        ct = chargeType or self.chargeType
        at = atomType or self.atomType

        if self.chargeVal:
            #print "1111111111"
            cmd = '%s -i %s -fi %s -o %s -fo mol2 -c %s -nc %s -m %s -s 2 -df 0 -at %s' % \
            (self.acExe, self.inputFile, self.ext[1:], self.acMol2File, ct,
             self.chargeVal, self.multiplicity, at)
        else:
            #print "2222222222"
            cmd = '%s -i %s -fi %s -o %s -fo mol2 -c %s -s 2 -df 0 -at %s' % \
            (self.acExe, self.inputFile, self.ext[1:], self.acMol2File, ct, at)


        if os.path.exists(self.acMol2File) and not self.force:
            print "AC output file present... doing nothing"
        else:
            try: os.remove(self.acMol2File)
            except: pass
            self.acLog = getoutput(cmd)

        if os.path.exists(self.acMol2File):
            print "==> Antechamber OK"
        else:
            print self.acLog
            return True

    def checkXyzAndTopFiles(self):
        fileXyz = self.baseName + '.xyz'
        fileTop = self.baseName + '.top'
        if os.path.exists(fileXyz) and os.path.exists(fileTop):
            self.acXyz = fileXyz
            self.acTop = fileTop
            return True
        return False

    def execSleap(self):

        self.makeDir()

        if self.ext == '.mol2':
            print "WARNING: Sleap doesn't work with mol2 files yet..."
            return True

        if self.chargeType != 'bcc':
            print "WARNING: Sleap works only with bcc charge method"
            return True

        sleapScpt = SLEAP_TEMPLATE % self.acParDict

        fp = open('sleap.in','w')
        fp.write(sleapScpt)
        fp.close()

        cmd = '%s -f sleap.in' % self.sleapExe

        if self.checkXyzAndTopFiles() and not self.force:
            print "Topologies files already present... doing nothing"
        else:
            try: os.remove(self.acTop) ; os.remove(self.acXyz)
            except: pass
            self.sleapLog = getoutput(cmd)

            if self.checkXyzAndTopFiles():
                print "==> Sleap OK"
            else:
                print self.sleapLog
                return True

    def execTleap(self):

        self.makeDir()

        if self.ext == ".pdb":
            print '... converting pdb input file to mol2 input file'
            if self.convertPdbToMol2():
                print "ERROR: convertPdbToMol2 failed"

        #print self.chargeVal

        if self.execAntechamber():
            print "ERROR: Antechamber failed"

        if self.execParmchk():
            print "ERROR: Parmchk failed"

        tleapScpt = TLEAP_TEMPLATE % self.acParDict

        fp = open('tleap.in','w')
        fp.write(tleapScpt)
        fp.close()

        cmd = '%s -f tleap.in' % self.tleapExe

        if self.checkXyzAndTopFiles() and not self.force:
            print "Topologies files already present... doing nothing"
        else:
            try: os.remove(self.acTop) ; os.remove(self.acXyz)
            except: pass
            self.tleapLog = getoutput(cmd)

        if self.checkXyzAndTopFiles():
            print "==> Tleap OK"
        else:
            print self.tleapLog
            return True

    def execParmchk(self):

        self.makeDir()
        cmd = '%s -i %s -f mol2 -o frcmod' % (self.parmchkExe, self.acMol2File)
        self.parmchkLog = getoutput(cmd)

        if os.path.exists('frcmod'):
            print "==> Parmchk OK"
        else:
            print self.parmchkLog
            return True

    def convertPdbToMol2(self):
        if self.ext == '.pdb':
            if self.execBabel():
                print "ERROR: convert pdb to mol2 via babel failed"
                return True

    def execBabel(self):

        self.makeDir()

        cmd = '%s -ipdb %s.pdb -omol2 %s.mol2' % (self.babelExe, self.baseName,
                                                  self.baseName)
        self.babelLog = getoutput(cmd)
        self.ext = '.mol2'
        self.inputFile = self.baseName+self.ext
        self.acParDict['ext'] = 'mol2'
        if os.path.exists(self.inputFile):
            print "==> Babel OK"
        else:
            print self.babelLog
            return True

    def makeDir(self):

        os.chdir(self.rootDir)
        self.absHomeDir = os.path.abspath(self.homeDir)
        if not os.path.exists(self.homeDir):
            os.mkdir(self.homeDir)
        os.chdir(self.homeDir)
        if not os.path.exists(self.inputFile):
            copy2(self.absInputFile, '.')

        return True

    def createACTopol(self):
        """
            If successful, a Amber Top and Xyz file will be generated
        """
        if self.execSleap():
            print "ERROR: Sleap failed"
            print "... trying Tleap"
            if self.execTleap():
                print "ERROR: Tleap failed"
        self.pickleSave()

    def pickleSave(self):
        """
            To restore:
                from acpypi import *
                import cPickle as pickle
                o = pickle.load(open('DDD.pkl'))
            It's fails with ipython
        """
        pklFile = self.baseName+".pkl"
        if not os.path.exists(pklFile):
            print "Writing pickle file %s" % pklFile
            pickle.dump(self, open(pklFile,"w"))
        elif self.force:
            print "Overwriting pickle file %s" % pklFile
            pickle.dump(self, open(pklFile,"w"))
        else:
            print "Pickle file %s already present... doing nothing" % pklFile

    def createMolTopol(self):
        """
            Create molTop obj
        """
        self.molTopol = MolTopol(self)

class MolTopol:
    """"
        http://amber.scripps.edu/formats.html (not updated to amber 10 yet)
        Parser, take information in AC xyz and top files and convert to objects
        INPUTS: acFileXyz and acFileTop
        RETURN: molTopol obj or None
    """
    def __init__(self, acTopolObj = None, acFileXyz = None, acFileTop = None):

        if acTopolObj:
            if not acFileXyz: acFileXyz = acTopolObj.acXyz
            if not acFileTop: acFileTop = acTopolObj.acTop
        if not os.path.exists(acFileXyz) and not os.path.exists(acFileTop):
            print "ERROR: Files '%s' and '%s' don't exist"
            print "       molTopol object won't be created"
            return None

        self.xyzFile = open(acFileXyz, 'r').readlines()
        self.topFile = open(acFileTop, 'r').readlines()

        self.pointers = self.getFlagData('POINTERS')

        self.getResidueLabel()

        self.getAtoms()

        self.getBonds()

        self.getAngles()

        self.getDihedrals()

        # a list of FLAGS from acTopFile that matter
        self.flags = ( 'POINTERS', 'ATOM_NAME', 'CHARGE', 'MASS', 'ATOM_TYPE_INDEX',
                  'NUMBER_EXCLUDED_ATOMS', 'NONBONDED_PARM_INDEX',
                  'RESIDUE_LABEL', 'BOND_FORCE_CONSTANT', 'BOND_EQUIL_VALUE',
                  'ANGLE_FORCE_CONSTANT', 'ANGLE_EQUIL_VALUE',
                  'DIHEDRAL_FORCE_CONSTANT', 'DIHEDRAL_PERIODICITY',
                  'DIHEDRAL_PHASE', 'AMBER_ATOM_TYPE' )

    def getFlagData(self, flag):
        """
            For a given acFileTop flag, return a list of the data related
        """
        block = False
        tFlag = '%FLAG ' + flag
        data = ''

        for rawLine in self.topFile:
            line = rawLine[:-1]
            if '%FORMAT' in line: continue
            if tFlag in line:
                block = True
                continue
            if block and '%FLAG ' in line: break
            if block: data += line
        sdata = data.split()
        if '+' and '.' in data: # it's a float
            ndata = map(float, sdata)
        else:
            try: # try if it's integer
                ndata = map(int, sdata)
            except: # it's string
                ndata = sdata
        return ndata # a list

    def getCoords(self):
        """
            For a given acFileXyz file, return a list of coords as:
            [[x1,y1,z1],[x2,y2,z2], etc.]
        """
        data = ''
        for rawLine in self.xyzFile[2:]:
            line = rawLine[:-1]
            data += line
        sdata = data.split()
        ndata = map(float, sdata)

        gdata = []
        for i in range(0, len(ndata), 3):
            gdata.append([ndata[i], ndata[i+1], ndata[i+2]])

        return gdata

    def getResidueLabel(self):
        """
            Get a 3 capital letters code from acFileTop
        """
        residueLabel = self.getFlagData('RESIDUE_LABEL')
        self.residueLabel = residueLabel[0]

    def getAtoms(self):
        """
            Set a list with all atoms objects build from dat in acFileTop
            Set also if molTopol atom type system is gaff or amber
            Set also list atomTypes
            Set also molTopol total charge
        """
        atomNameList = self.getFlagData('ATOM_NAME')
        atomTypeNameList = self.getFlagData('AMBER_ATOM_TYPE')
        self._atomTypeNameList = atomTypeNameList
        massList = self.getFlagData('MASS')
        chargeList = self.getFlagData('CHARGE')
        #uniqAtomTypeId = self.getFlagData('ATOM_TYPE_INDEX') # for LJ
        balanceChargeList = self.balanceCharges(chargeList)
        coords = self.getCoords()
        ACOEFs, BCOEFs = self.getABCOEFs()

        atoms = []
        atomTypes = []
        tmpList = [] # a list with unique atom types
        totalCharge = 0.0
        for atomName in atomNameList:
            id = atomNameList.index(atomName)
            atomTypeName = atomTypeNameList[id]
            mass = massList[id]
            charge = balanceChargeList[id]
            totalCharge += charge
            coord = coords[id]
            ACOEF = ACOEFs[id]
            BCOEF = BCOEFs[id]
            if atomTypeName not in tmpList:
                tmpList.append(atomTypeName)
                atomType = AtomType(atomTypeName, mass, ACOEF, BCOEF)
                atomTypes.append(atomType)
            atom = Atom(atomName, atomType, mass, charge, coord, ACOEF, BCOEF)
            atoms.append(atom)

        if atomTypeName[0].islower:
            self.atomTypeSystem = 'gaff'
        else:
            self.atomTypeSystem = 'amber'

        self.totalCharge = int(totalCharge)

        self.atoms = atoms
        self.atomTypes = atomTypes

    def getBonds(self):
        uniqKbList = self.getFlagData('BOND_FORCE_CONSTANT')
        uniqReqList = self.getFlagData('BOND_EQUIL_VALUE')
        #numbnd = len(uniqKbList)
        # for list below, true atom number = index/3 + 1
        bondCodeHList = self.getFlagData('BONDS_INC_HYDROGEN')
        bondCodeNonHList = self.getFlagData('BONDS_WITHOUT_HYDROGEN')
        bondCodeList = bondCodeHList + bondCodeNonHList
        bonds = []
        for i in range(0, len(bondCodeList), 3):
            idAtom1 = bondCodeList[i] / 3 # remember python starts with id 0
            idAtom2 = bondCodeList[i+1] / 3
            bondTypeId = bondCodeList[i+2] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            kb = uniqKbList[bondTypeId]
            req = uniqReqList[bondTypeId]
            bond = Bond(atom1, atom2, kb, req)
            bonds.append(bond)
        self.bonds = bonds

    def getAngles(self):
        uniqKtList = self.getFlagData('ANGLE_FORCE_CONSTANT')
        uniqTeqList = self.getFlagData('ANGLE_EQUIL_VALUE')
        # for list below, true atom number = index/3 + 1
        angleCodeHList = self.getFlagData('ANGLES_INC_HYDROGEN')
        angleCodeNonHList = self.getFlagData('ANGLES_WITHOUT_HYDROGEN')
        angleCodeList = angleCodeHList + angleCodeNonHList
        angles = []
        for i in range(0, len(angleCodeList), 4):
            idAtom1 = angleCodeList[i] / 3 # remember python starts with id 0
            idAtom2 = angleCodeList[i+1] / 3
            idAtom3 = angleCodeList[i+2] / 3
            angleTypeId = angleCodeList[i+3] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            atom3 = self.atoms[idAtom3]
            kt = uniqKtList[angleTypeId]
            teq = uniqTeqList[angleTypeId] # angle given in rad in prmtop
            angle = Angle(atom1, atom2, atom3, kt, teq)
            angles.append(angle)
        self.angles = angles

    def getDihedrals(self):
        uniqKpList = self.getFlagData('DIHEDRAL_FORCE_CONSTANT')
        uniqPeriodList = self.getFlagData('DIHEDRAL_PERIODICITY')
        uniqPhaseList = self.getFlagData('DIHEDRAL_PHASE')
        # for list below, true atom number = abs(index)/3 + 1
        dihCodeHList = self.getFlagData('DIHEDRALS_INC_HYDROGEN')
        dihCodeNonHList = self.getFlagData('DIHEDRALS_WITHOUT_HYDROGEN')
        dihCodeList = dihCodeHList + dihCodeNonHList
        properDih = []
        improperDih = []
        for i in range(0, len(dihCodeList), 5):
            idAtom1 = dihCodeList[i] / 3 # remember python starts with id 0
            idAtom2 = dihCodeList[i+1] / 3
            # 3 and 4 indexes can be negative: if id3 < 0, end group interations
            # in amber are to be ignored; if id4 < 0, dihedral is improper
            # TODO: investigate if id3 signal is important for CNS and GMX.
            # using only id4 signal for the moment.
            idAtom3 = dihCodeList[i+2] / 3 # can be negative
            idAtom4raw = dihCodeList[i+3] / 3 # can be negative
            idAtom4 = idAtom4raw
            dihTypeId = dihCodeList[i+4] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            atom3 = self.atoms[idAtom3]
            atom4 = self.atoms[idAtom4]
            kPhi = uniqKpList[dihTypeId] # already divided by IDIVF
            period = int(uniqPeriodList[dihTypeId]) # integer
            phase = uniqPhaseList[dihTypeId]# angle given in rad in prmtop
            dihedral = Dihedral(atom1, atom2, atom3, atom4, kPhi, period, phase)
            if idAtom4raw > 0:
                properDih.append(dihedral)
            else:
                improperDih.append(dihedral)

        self.properDihedrals = properDih
        self.improperDihedrals = improperDih

    def balanceCharges(self, chargeList):
        """
            Note that python is very annoying about floating points.
            Even after balance, there will always be some residue of
            order e-12 to e-16
        """
        total = sum(chargeList)
        maxVal = max(chargeList)
        minVal = min(chargeList)
        if abs(maxVal) >= abs(minVal): lim = maxVal
        else: lim = minVal
        limId = chargeList.index(lim)
        diff = total - int(total)
        fix = lim - diff
        chargeList[limId] = fix
        return chargeList

    def getABCOEFs(self):
        uniqAtomTypeIdList = self.getFlagData('ATOM_TYPE_INDEX')
        nonBonIdList = self.getFlagData('NONBONDED_PARM_INDEX')
        rawACOEFs = self.getFlagData('LENNARD_JONES_ACOEF')
        rawBCOEFs = self.getFlagData('LENNARD_JONES_BCOEF')
        #print nonBonIdList, len(nonBonIdList), rawACOEFs, len(rawACOEFs)
        ACOEFs = []
        BCOEFs = []
        ntypes = max(uniqAtomTypeIdList)
        for atName in self._atomTypeNameList:
            id = self._atomTypeNameList.index(atName)
            atomTypeId = uniqAtomTypeIdList[id]
            index = ntypes * (atomTypeId - 1) + atomTypeId
            nonBondId = nonBonIdList[index - 1]
            #print "*****", index, ntypes, atName, id, atomTypeId, nonBondId
            ACOEFs.append(rawACOEFs[nonBondId - 1])
            BCOEFs.append(rawBCOEFs[nonBondId - 1])
        #print ACOEFs
        return ACOEFs, BCOEFs

    def writeGromacsTopolFiles(self):
        """
            # from ~/Programmes/amber10/dat/leap/parm/gaff.dat
            #atom type        atomic mass        atomic polarizability        comments
            ca                12.01                 0.360                    Sp2 C in pure aromatic systems
            ha                1.008                 0.135                    H bonded to aromatic carbon

            #bonded atoms        harmonic force kcal/mol/A^2       eq. dist. Ang.  comments
            ca-ha                  344.3*                           1.087**         SOURCE3  1496    0.0024    0.0045
            * for gmx: 344.3 * 4.184 * 100 * 2 = 288110 kJ/mol/nm^2 (why factor 2?)
            ** convert Ang to nm ( div by 10) for gmx: 1.087 A = 0.1087 nm
            # CA HA         1    0.10800   307105.6 ; ged from 340. bsd on C6H6 nmodes; PHE,TRP,TYR (from ffamber99bon.itp)
            # CA-HA  367.0    1.080       changed from 340. bsd on C6H6 nmodes; PHE,TRP,TYR (from parm99.dat)

            # angle        HF kcal/mol/rad^2    eq angle degrees     comments
            ca-ca-ha        48.5*             120.01                SOURCE3 2980   0.1509   0.2511
            * to convert to gmx: 48.5 * 4.184 * 2 = 405.848 kJ/mol/rad^2 (why factor 2?)
            # CA  CA  HA           1   120.000    418.400 ; new99 (from ffamber99bon.itp)
            # CA-CA-HA    50.0      120.00 (from parm99.dat)

            # dihedral    idivf        barrier hight/2 kcal/mol  phase degrees       periodicity     comments
            X -ca-ca-X    4           14.500*                     180.000                2.000             intrpol.bsd.on C6H6
            * to convert to gmx: 14.5/4 * 4.184 * 2 (?) (yes in amb2gmx, no in topolbuild, why?) = 30.334 or 15.167 kJ/mol
            # X -CA-CA-X    4   14.50        180.0             2.         intrpol.bsd.on C6H6 (from parm99.dat)
            # X   CA  CA  X     3    30.33400     0.00000   -30.33400     0.00000     0.00000     0.00000   ; intrpol.bsd.on C6H6
            ;propers treated as RBs in GROMACS to use combine multiple AMBER torsions per quartet (from ffamber99bon.itp)

            # impr. dihedral        barrier hight/2      phase degrees       periodicity     comments
            X -X -ca-ha             1.1*                  180.                      2.                   bsd.on C6H6 nmodes
            * to convert to gmx: 1.1 * 4.184 = 4.6024 kJ/mol/rad^2
            # X -X -CA-HA         1.1          180.          2.           bsd.on C6H6 nmodes (from parm99.dat)
            # X   X   CA  HA       1      180.00     4.60240     2      ; bsd.on C6H6 nmodes
            ;impropers treated as propers in GROMACS to use correct AMBER analytical function (from ffamber99bon.itp)

            # 6-12 parms     sigma = 2 * r * 2^(-1/6)    epsilon
            # atomtype        radius Ang.                    pot. well depth kcal/mol      comments
              ha                  1.4590*                      0.0150**                         Spellmeyer
              ca                  1.9080                    0.0860                            OPLS
            * to convert to gmx:
                sigma = 1.4590 * 2^(-1/6) * 2 = 2 * 1.29982 Ang. = 2 * 0.129982 nm  = 1.4590 * 2^(5/6)/10 =  0.259964 nm
            ** to convert to gmx: 0.0150 * 4.184 = 0.06276 kJ/mol
            # amber99_3    CA     0.0000  0.0000  A   3.39967e-01  3.59824e-01 (from ffamber99nb.itp)
            # amber99_22   HA     0.0000  0.0000  A   2.59964e-01  6.27600e-02 (from ffamber99nb.itp)
            # C*          1.9080  0.0860             Spellmeyer
            # HA          1.4590  0.0150             Spellmeyer (from parm99.dat)
            # to convert r and epsilon to ACOEF and BCOEF
            # ACOEF = sqrt(e1*e2) * (r1 + r2)^12 ; BCOEF = 2 * sqrt(e1*e2) * (r1 + r2)^6 = 2 * ACOEF/(r1+r2)^6
            #   ca   ca       819971.66        531.10
            #   ca   ha        76245.15        104.66
            #   ha   ha         5716.30         18.52
        """

class Atom:
    """
        Charges in prmtop file has to be divide by 18.2223 to convert to charge
        in units of the electron charge.
        To convert ACOEF and BCOEF to r0 (Ang.) and epsilon (kcal/mol), as seen
        in gaff.dat for example; same atom type (i = j):
            r0 = 1/2 * (2 * ACOEF/BCOEF)^(1/6)
            epsilon = 1/(4 * A) * BCOEF^2
        To convert r0 and epsilon to ACOEF and BCOEF
            ACOEF = sqrt(ep_i * ep_j) * (r0_i + r0_j)^12
            BCOEF = 2 * sqrt(ep_i * ep_j) * (r0_i + r0_j)^6
                  = 2 * ACOEF/(r0_i + r0_j)^6
        where index i and j for atom types.
        Coord is given in Ang. and mass in Atomic Mass Unit.
    """
    def __init__(self, atomName, atomType, mass, charge, coord, ACOEF, BCOEF):
        self.atomName = atomName
        self.atomType = atomType
        self.mass = mass
        self.charge = charge / 18.2223
        self.coords = coord

class AtomType:
    """
        AtomType per atom in gaff or amber.
    """
    def __init__(self, atomTypeName, mass, ACOEF, BCOEF):
        self.atomTypeName = atomTypeName
        self.mass = mass
        self.ACOEF = ACOEF
        self.BCOEF = BCOEF

class Bond:
    """
        attributes: pair of Atoms, spring constant (kcal/mol), dist. eq. (Ang)
    """
    def __init__(self, atom1, atom2, kBond, rEq):
        self.atom1 = atom1
        self.atom2 = atom2
        self.kBond = kBond
        self.rEq = rEq

class Angle:
    """
        attributes: 3 Atoms, spring constant (kcal/mol/rad^2), angle eq. (rad)
    """
    def __init__(self, atom1, atom2, atom3, kTheta, thetaEq):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.kTheta = kTheta
        self.thetaEq = thetaEq # rad, to convert to degree: thetaEq * 180/Pi

class Dihedral:
    """
        attributes: 4 Atoms, spring constant (kcal/mol), periodicity,
        phase (rad)
    """
    def __init__(self, atom1, atom2, atom3, atom4, kPhi, period, phase):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.kPhi = kPhi
        self.period = period
        self.phase = phase # rad, to convert to degree: thetaEq * 180/Pi

class CnsWriter:
    pass

if __name__ == '__main__':
    argsDict = parseArgs(sys.argv[1:])
    iF = argsDict.get('-i')
    cT = argsDict.get('-c', 'bcc')
    cV = argsDict.get('-n', None)
    mt = argsDict.get('-m', '1')
    at = argsDict.get('-F', 'gaff')
    fs = False
    if '-f' in argsDict.keys(): fs = True

    molecule = ACTopol(iF, chargeType = cT, chargeVal = cV,
                               multiplicity = mt, atomType = at, force = fs)

    #molecule.convertPdbToMol2()
    #print molecule.babelLog

    #molecule.execAntechamber()
    #print molecule.acLog

    #molecule.execSleap()
    #print molecule.sleapLog

    #molecule.execTleap()
    #print molecule.tleapLog, molecule.parmchkLog

    molecule.createACTopol()
    molecule.createMolTopol()
