#!/usr/bin/env python

import sys
import os

hisDictConv = {'HHH':['HIE', 'HISB'], 'JJJ':['HIP', 'HISH'], 'OOO':['HID', 'HISA']}

pymolScript = \
'''
aa_dict = {'A' : 'ala','C' : 'cys','D' : 'asp','E' : 'glu','F' : 'phe',
           'G' : 'gly','H' : 'hie','J' : 'hip','O' : 'hid','I' : 'ile',
           'K' : 'lys','L' : 'leu','M' : 'met','N' : 'asn','P' : 'pro',
           'Q' : 'gln','R' : 'arg','S' : 'ser','T' : 'thr','V' : 'val',
           'W' : 'trp','Y' : 'tyr'}

def build(object_name, sequence, first_residue = "1"):
    if len(sequence):
        code = sequence[0]
        cmd.fragment(aa_dict[code],object_name)
        cmd.alter(object_name,'resi="%s"'%first_residue)
        cmd.edit(object_name+" and name C")
        for code in sequence[1:]:
            editor.attach_amino_acid("pk1",aa_dict[code])
        cmd.edit()

def build3(object_name, sequence, first_residue = "1"):
    if len(sequence):
        code = sequence[0]
        cmd.fragment('nt_'+aa_dict[code],object_name)
        cmd.alter(object_name,'resi="%s"'%first_residue)
        cmd.edit(object_name+" and name C")
        for code in sequence[1:-1]:
            editor.attach_amino_acid("pk1",aa_dict[code])
        editor.attach_amino_acid("pk1",'ct_'+aa_dict[code])
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

def hisConv(his,opt):
    ''' opt = 0 : pdb to opls
        opt = 1 : apdb to amber
    '''
    his1, his2 = hisDictConv[his]
    if opt:
        cmd = "sed s/%s/%s/g a%s.pdb > tMpFile; mv tMpFile a%s.pdb" % ('HIS',his1,his,his)
    else:
        cmd = "sed s/%s\ /%s/g %s.pdb > tMpFile; mv tMpFile %s.pdb" % (his1,his2,his,his)
    os.system(cmd)

def createAmberPdb(fpdb):
#    pdb = os.path.splitext(fpdb)[0]
    fLines = file(fpdb).readlines()
    famb = open(fpdb.replace('g','a'),'w')
    for line in fLines:
        if 'ATOM  ' in line:
            line = line.replace('CYS','CYN')
            line = line.replace('LYS','LYP')
#            line = line.replace('3HB','1HB')
#            line = line.replace('3HG','1HG')
#            line = line.replace('3HA','1HA')
#            line = line.replace('3HD','1HD')
#            line = line.replace('3HE','1HE')
#            #line = line.replace('3HZ','1HZ')
#            line = line.replace(' HA  GLY','2HA  GLY')
#            line = line.replace('CD1 ILE','CD  ILE')
#            line = line.replace('1HD1 ILE',' HD1 ILE')
#            line = line.replace('2HD1 ILE',' HD2 ILE')
#            line = line.replace('3HD1 ILE',' HD3 ILE')
            if line[22:26] == '   1':
                res = line[17:20]
                line = line.replace('%s ' % res, 'N%s' % res)
            if line[22:26] == '   3':
                res = line[17:20]
                line = line.replace('%s ' % res, 'C%s' % res).replace(' O1 ', ' OC1').replace(' O2 ',' OC2')
        famb.write(line)
    famb.close()

if __name__ == '__main__':
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
    listPdb = os.listdir('.')
    listPdb.sort()
    for pdbFile in listPdb[:]:
        pdb, ext = os.path.splitext(pdbFile)
        if len(pdbFile) == 7 and ext == ".pdb":
            if pdb in hisDictConv.keys():
                hisConv(pdb, 0)
            cmd = 'pdb2gmx -f %s.pdb -o g%s.pdb -p o%s.top -ff oplsaa -ignh' % \
                (pdb, pdb, pdb)
            out = os.system(cmd)
            createAmberPdb("g%s.pdb" % pdb)
            if pdb in hisDictConv.keys():
                hisConv(pdb, 1)
            cmd = 'pdb2gmx -f a%s.pdb -o a%s.pdb -p a%s.top -ff amber99sb' % \
                (pdb, pdb, pdb)
            out = os.system(cmd)
            cmd = "acpypi -di g%s.pdb -c gas" % pdb
            out = os.system(cmd)
            if out :
                print "pdb2gmx for %s FAILED... aborting." % pdb
                print cmd
                sys.exit(1)
    os.system('rm -f %s/\#*' % tmpDir)
