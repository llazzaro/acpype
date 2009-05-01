#!/usr/bin/env python
import commands
import sys #@UnusedImport
import os

aa_dict = {'A':'ala', 'C':'cys', 'D':'asp', 'E':'glu', 'F':'phe', 'G':'gly',
           'H':'hie', 'J':'hip', 'O':'hid', 'I':'ile', 'K':'lys', 'L':'leu',
           'M':'met', 'N':'asn', 'P':'pro', 'Q':'gln', 'R':'arg', 'S':'ser',
           'T':'thr', 'V':'val', 'W':'trp', 'Y':'tyr'}

hisDictConv = {'HHH':['HIE', 'HISB'], 'JJJ':['HIP', 'HISH'], 'OOO':['HID', 'HISA']}

pymolScript = 'aa_dict = %s' % str(aa_dict) + \
'''
def build3(object_name, sequence, first_residue = "1"):
    if len(sequence):
        codeNT, code, codeCT = sequence
        if code in ['J', 'R']:
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

def hisConv(fName, his, opt):
    ''' opt = 0 : pymol pdb to opls
        opt = 1 : apdb (from opdb) to amber pdb
    '''
    hisA, hisO = hisDictConv[his]
    if opt == 0:
        cmd = "sed s/%s\ /%s/g %s > tMpFile; mv tMpFile %s" % (hisA,hisO,fName,fName)
    elif opt == 1:
        cmd = "sed s/%s/%s/g %s > tMpFile; mv tMpFile %s" % ('HIS',hisA,fName,fName)
    os.system(cmd)

def createOplsPdb(fpdb):
    '''pdb -> opdb (opls)'''
    fLines = file(fpdb).readlines()
    fName = 'o'+fpdb
    fopls = open(fName,'w')
    for line in fLines:
        if 'ATOM  ' in line:
            if 'AAA' not in fpdb:
                line = line.replace(' 3HB ',' 1HB ')
            line = line.replace(' 3HG ',' 1HG ')
            line = line.replace('3HA  GLY','1HA  GLY')
            line = line.replace(' HA  GLY','2HA  GLY')
            line = line.replace(' 3HD ',' 1HD ')
            line = line.replace('3HE  LYS','1HE  LYS')
            line = line.replace('3HG1 ILE','1HG1 ILE')
            line = line.replace('3H   PRO','1H   PRO')
        fopls.write(line)
    fopls.close()
    return fName

def createAmberPdb(fpdb, opt):
    '''pdb -> apdb (amber)
       opt = 0 : pymol pdb to apdb
       opt = 1 : ogpdb to apdb
    '''
    fLines = file(fpdb).readlines()
    fName = 'a'+fpdb
    famb = open(fName,'w')
    for line in fLines:
        if 'ATOM  ' in line:
            line = line.replace('CYS','CYN')
            line = line.replace('LYS','LYP')
            if opt == 0:
                line = line.replace('CYS','CYN')
                if 'AAA' not in fpdb:
                    line = line.replace(' 3HB ',' 1HB ')
                line = line.replace(' 3HG ',' 1HG ')
                line = line.replace(' 3HA ',' 1HA ')
                line = line.replace(' HA  GLY','2HA  GLY')
                line = line.replace('CD1 ILE','CD  ILE')
                line = line.replace('3HG1 ILE','1HG1 ILE')
                line = line.replace('1HD1 ILE',' HD1 ILE')
                line = line.replace('2HD1 ILE',' HD2 ILE')
                line = line.replace('3HD1 ILE',' HD3 ILE')
                line = line.replace('LYS','LYP')
                line = line.replace(' 3HD ',' 1HD ')
                line = line.replace('3HE  LYP','1HE  LYP')
                line = line.replace('3H   PRO','1H   PRO')
                if line[22:26] == '   3':
                    line = line.replace('  O  ', '  OC1').replace('  OXT','  OC2')
            elif opt == 1:
                if line[22:26] == '   3':
                    line = line.replace(' O1 ', ' OC1').replace(' O2 ',' OC2')
            if line[22:26] == '   1':
                res = line[17:20]
                line = line.replace('%s ' % res, 'N%s' % res)
            if line[22:26] == '   3':
                res = line[17:20]
                line = line.replace('%s ' % res, 'C%s' % res)
        famb.write(line)
    famb.close()
    return fName

def parseTopFile(lista):
    ''' lista = top file in list format
        return a dict with fields'''
    defDict = {}
    parDict = {}
    for line in lista:
        line = line.split(';')[0]
        line = line.strip()
        if line.startswith('#def'):
            vals = line.split()
            defDict[vals[1]] = map(eval,vals[2:])
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
#    mRes = parseTopFile(file('g%s.acpypi/g%s_GMX.itp' % (res,res)).readlines())

#    flags = ['pairs', 'bonds', 'angles', 'dihedrals'], ['dihedraltypes', 'angletypes', 'bondtypes']

    agPairs = agRes[0]['pairs']
    ogPairs = ogRes[0]['pairs']
#    mPairs = mRes[0]['pairs']

    if agPairs != ogPairs:
        print "!!!!!!", res

def roundAllFloats(lista,l):
    nlista = []
    for ii in lista:
        tt = ii[:l+1]
        for jj in ii[l+1:]:
            if jj > 100.0:
                jj = round(jj,-1)
            nn = round(jj,3)
            tt.append(nn)
        nlista.append(tt)
    return nlista

def checkTopAcpypi(res, ff):
    '''Compare acpypi gmx itp against amber or opls pdb2gmx results'''
    global amber2oplsDict, ac2opls
    def addParam(l,item):
        dict = {2:'bondtypes', 3:'angletypes', 4:'dihedraltypes'}
        dType = {} # dict
        for type in ffBon[0][dict[l]]:
            i = type[:l+1]
            j = type[l+1:]
            dType[str(i)] = j
        entries = []
        lNum = item[l]
        ent = [ffDictAtom[x] for x in item[:l]]
        rent = ent[:]
        rent.reverse()
        entries.append(ent+[lNum])
        entries.append(rent+[lNum])
        if l == 4:
            if len(item) == 6:
                par = ffBon[1][item[5]]
                return par
            tent = ent[:]
            rtent = rent[:]
            if lNum == 3: # dih proper
                tent[0] = 'X'
                entries.append(tent+[lNum])
                tent.reverse()
                entries.append(tent+[lNum])
                tent[0] = 'X'
                entries.append(tent+[lNum])
                tent.reverse()
                entries.append(tent+[lNum])
                rtent[0] = 'X'
                entries.append(rtent+[lNum])
                rtent.reverse()
                entries.append(rtent+[lNum])
            if lNum == 1: # dih improp
                tent[0] = 'X'
                entries.append(tent+[lNum])
                tent[1] = 'X'
                entries.append(tent+[lNum])
                rtent[0] = 'X'
                entries.append(rtent+[lNum])
                rtent[1] = 'X'
                entries.append(rtent+[lNum])
        found = False
        for e in entries:
            try:
                par = dType[str(e)]
                found = True
                break
            except: pass
        if not found:
            print "%s %s %s not found in %s Bon" % (dict[l], ent, item, ffType)
        return par

    compareParameters = False
    if ff == 'a' : compareParameters = True

    agRes = parseTopFile(file('ag%s.top' % (res)).readlines())
    ogRes = parseTopFile(file('og%s.top' % (res)).readlines())
    #ffgRes = parseTopFile(file('%sg%s.top' % (ff,res)).readlines())
    acRes = parseTopFile(file('%sg%s.acpypi/%sg%s_GMX.itp' % (ff,res,ff,res)).readlines())
    if ff == 'a':
        ffNb = aNb
        ffBon = aBon
        ffgRes = agRes
    if ff == 'o':
        ffNb = oNb
        ffBon = oBon
        ffgRes = ogRes

    ffDictAtom = {} # dict link res atom numbers to amber or opls atom types
    for item in ffgRes[0]['atoms']:
        i,j = item[:2]
        ffDictAtom[i] = ffNb[j]

    acDictAtom = {}
    for item in acRes[0]['atoms']:
        i,j = item[:2]
        acDictAtom[i] = j

    for k in acDictAtom.keys():
        if amber2oplsDict.has_key(acDictAtom[k]):
            amber2oplsDict[acDictAtom[k]].add(ffDictAtom[k])
        else:
            amber2oplsDict[acDictAtom[k]] = set([ffDictAtom[k]])

    # to build a dict acpypi atom types (amber or gaff) to [[opls_code],[opls_mass]]
    #ac2opls
    if ff == 'o':
        for item in acRes[0]['atoms']:
            id, at = item[:2]
            it = ogRes[0]['atoms'][id-1]
            ocode, omass = it[1],it[-1]
            #print item, it, ocode, omass
            if ac2opls.has_key(at):
                ac2opls[at][0].add(ocode)
                ac2opls[at][1].add(omass)
            else:
                ac2opls[at] = [set([ocode]),set([omass])]

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

    flags = [('pairs',2), ('bonds',2), ('angles',3), ('dihedrals',4)] #, ['dihedraltypes', 'angletypes', 'bondtypes']

    for flag,l in flags:
        print "    ==> Comparing %s" % flag
        if flag != flags[0][0] and compareParameters:
            agres = []
            tAgRes = ffgRes[0][flag]
            for item in tAgRes:
                if flag == flags[1][0]:
                    par = addParam(l,item)
                elif flag == flags[2][0]:
                    par = addParam(l,item)
                elif flag == flags[3][0]:
                    par = addParam(l,item)
                    if len(item) == 6:
                        item.pop()
                agres.append(item + par)
            if compareParameters: acres = acRes[0][flag]
        else:
            if flag == flags[3][0]:
                agres = [x[:l+1] for x in ffgRes[0][flag]]
            else: agres = ffgRes[0][flag]
            if compareParameters: acres = acRes[0][flag]
            else: acres = [x[:l+1] for x in acRes[0][flag]]

        act = acres[:]
        agt = agres[:]

        if compareParameters:
            if flag != flags[0][0]:
                # round parameters
                act = roundAllFloats(act,l)
                agt = roundAllFloats(agt,l)
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
            print "    ac: ", act3
        if agt3:
            agt3.sort()
            print "    %sg: " % ff, agt3

        if not act3 and not agt3:
            print "        %s OK" % flag

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
                print line
        elif 'Total charge' in line:
            print "    ",line
        elif 'ERROR' in line.upper():
            if 'Fatal error' in line:
                print line, lines[count+1]
            else:
                print line
        count += 1

def fixRes4Acpypi(fpdb):
    code = fpdb[2]
    fLines = file(fpdb).readlines()
    famb = open(fpdb,'w')
    for line in fLines:
        if 'ATOM  ' in line:
            line = line[:17]+aa_dict[code].upper()+line[20:]
        famb.write(line)
    famb.close()

if __name__ == '__main__':

    '''order: (pymol) AAA.pdb -f-> oAAA.pdb --> (pdb2gmx) ogAAA.pdb -f->
         -f-> aogAAA.pdb --> (pdb2gmx) agAAA.pdb -f-> acpypi
    '''
    global amber2oplsDict, a2oD, ac2opls
    amber2oplsDict = {}
    a2oD = {}
    ac2opls = {}
    aNb = nbDict(file('/sw/share/gromacs/top/ffamber99sbnb.itp').readlines())
    oNb = nbDict(file('/sw/share/gromacs/top/ffoplsaanb.itp').readlines())

    aBon = parseTopFile(file('/sw/share/gromacs/top/ffamber99sbbon.itp').readlines())
    oBon = parseTopFile(file('/sw/share/gromacs/top/ffoplsaabon.itp').readlines())

    tmpDir = '/tmp/testAcpypi'
    tmpFile = 'tempScript.py'
    os.putenv('PATH', '/usr/bin:/bin:/usr/sbin:/sbin:/sw/bin:/usr/local/bin:/Users/alan/Programmes/amber10/exe')
    exePymol = 'pymol' # '/sw/bin/pymol'
#    print os.environ['PATH']
    if not os.path.exists(tmpDir):
        os.mkdir(tmpDir)
    os.chdir(tmpDir)
    ff = open(tmpFile,'w')
    ff.writelines(pymolScript)
    ff.close()
    cmd = "%s -qc %s" % (exePymol,tmpFile)
    os.system(cmd)
    listRes = os.listdir('.')
    listRes.sort()
    for resFile in listRes:
        res, ext = os.path.splitext(resFile) # eg. res = 'AAA'
        if len(resFile) == 7 and ext == ".pdb":
            print "File %s : residue %s" % (resFile, aa_dict[res[0]])

            # from pdb pymol to opls pdb and top
            opdb = createOplsPdb(resFile)
            ogpdb = 'og%s' % resFile
            ogtop = 'og%s.top' % res
            if res in hisDictConv.keys():
                hisConv(opdb, res, 0) # fix HIE etc. from pymol pdb to pdb opls
            cmd = 'pdb2gmx -f %s -o %s -p %s -ff oplsaa' % \
                (opdb, ogpdb, ogtop) # coords for O term changed
            out = commands.getstatusoutput(cmd)
            #parseOut(out[1])

            # from ogpdb to amber pdb and top
            apdb = createAmberPdb(ogpdb, 1)
            agpdb = 'ag%s.pdb' % res
            agtop = 'ag%s.top' % res
            if res in hisDictConv.keys():
                hisConv(apdb, res, 1)
            cmd = 'pdb2gmx -f %s -o %s -p %s -ff amber99sb' % \
                (apdb, agpdb, agtop)
            out = commands.getstatusoutput(cmd)
            #parseOut(out[1])

            # acpypi on agpdb file
#            fixRes4Acpypi(agpdb)
#            ffType = 'amber' # gaff
#            cType = 'gas' # bcc
#            cmd = "acpypi -dfi %s -c %s -a %s" % (agpdb, cType, ffType)
#            if res == 'JJJ' and cType == 'bcc': cmd += ' -n 1' # acpypi failed to get correct charge
#            out = commands.getstatusoutput(cmd)
#            #parseOut(out[1])
#            checkTopAcpypi(res)

            # acpypi on ogpdb file
            fixRes4Acpypi(ogpdb)
            ffType = 'gaff' #'amber' # gaff
            cType = 'gas' # bcc
            cmd = "acpypi -dfi %s -c %s -a %s" % (ogpdb, cType, ffType)
            if res == 'JJJ' and cType == 'bcc': cmd += ' -n 1' # acpypi failed to get correct charge
            out = commands.getstatusoutput(cmd)
            parseOut(out[1])
            checkTopAcpypi(res, 'o')

            if out[0] :
                print "pdb2gmx for %s FAILED... aborting." % res
                print cmd
                #sys.exit(1)
            #break
    os.system('rm -f %s/\#* posre.itp tempScript.py' % tmpDir)
    os.system("find . -name '*.itp' | xargs grep -v 'created by acpypi on' > standard_itp.txt")
    #print ac2opls
    tmp = {}
    for k,v in ac2opls.items():
        vOpls, vMass = v
        vOpls = [x.replace('opls_','') for x in vOpls]
        vMass = map(str,vMass)
        vOpls.sort()
        vMass.sort()
        if len(vMass) > 1:
            print k,v
        tmp[k] = vOpls + vMass
    print tmp

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