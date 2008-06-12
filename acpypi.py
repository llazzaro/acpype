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
        cmd = '%s -i %s.mol2 -fi mol2 -o tmp.mol2 -fo mol2 -c gas' % \
                                                    (self.acExe, self.baseName)
        log = getoutput(cmd)
        m = log.split()
        if len(m) >= 21:
            charge = int(float(m[14].replace('(','').replace(')','')))
        else:
            print "ERROR: guessCharge failed"
            os.chdir(localDir)
            return True
        self.chargeVal = charge
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

        self.makeDir()

        ct = chargeType or self.chargeType
        at = atomType or self.atomType

        if self.chargeVal:
            cmd = '%s -i %s -fi %s -o %s -fo mol2 -c %s -nc %s -m %s -s 2 -df 0 -at %s' % \
            (self.acExe, self.inputFile, self.ext[1:], self.acMol2File, ct,
             self.chargeVal, self.multiplicity, at)
        else:
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

        self.getResidueLabel()

        self.getAtoms()

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

    def getResidueLabel(self):
        """
            Get a 3 capital letters code from acFileTop
        """
        residueLabel = self.getFlagData('RESIDUE_LABEL')
        self.residueLabel = residueLabel[0]

    def getAtoms(self):
        """
            Returns a tuple with all atoms objects build from dat in acFileTop
            Set also if molTopol atom type system is gaff or amber
            Set also molTopol total charge
        """
        atomNameList = self.getFlagData('ATOM_NAME')
        atomTypeList = self.getFlagData('AMBER_ATOM_TYPE')
        massList = self.getFlagData('MASS')
        chargeList = self.getFlagData('CHARGE')

        atoms = []
        totalCharge = 0.0
        for atomName in atomNameList:
            id = atomNameList.index(atomName)
            atomType = atomTypeList[id]
            mass = massList[id]
            charge = chargeList[id]
            totalCharge += charge
            atom = Atom(atomName, atomType, mass, charge)
            atoms.append(atom)
        if atomType[0].islower:
            self.atomType = 'gaff'
        else:
            self.atomType = 'amber'

        self.totalCharge = int(totalCharge)

        self.atoms = atoms

class Atom:
    """
        charges in prmtop file has to be divide by 18.2223 to convert to charge
        in units of the electron charge
    """
    def __init__(self, atomName, atomType, mass, charge):
        self.atomName = atomName
        self.atomType = atomType
        self.mass = mass
        self.charge = charge / 18.2223

class Bond:
    """
        attributes: id, pair of Atoms
    """
    pass

class Angle:
    pass

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
