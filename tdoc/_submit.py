# -*- coding: utf-8 -*-
from ase import Atom
from math import ceil
from glob import glob
from rdkit import Chem
from shutil import copy
from random import sample
from re import search, findall
from collections import Counter
from os import chdir, listdir, makedirs

from tdoc._check import Check
from tdoc._conformer import Conformer



class Submit(object):
    
    """ To initiate parameters for Submit. """
    def __init__(self,para):
        self.para = para



    """ To get submitted file of M062X. """
    def get_M062X_out(self):
        
        # Go to working directory and determine absent molecules.
        chdir('{}/gaussian_out/M062X'.format(self.para['work_path']))
        smiles = []
        for smi in self.para['species']:
            if smi + '.dat' not in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], 'M062X')):
                if smi + '.xyz' not in listdir('../GFNFF'):
                    smiles.append(smi)
        
        # Check whether calculations are completed.
        Check(self.para, smiles, 'M062X').check_method_out()
        
        # Construct input files of M062X for Gaussian.
        for smi in smiles:
            if smi + '.gjf' not in listdir('.'):
                self.build_structure_from_smi(smi, self.para, self.para['conformer_method'])
                print('This molecule has been constructed: {}'.format(smi))
        
        # Generate submitted scripts.
        submitted_smiles = [smi for smi in smiles if smi + '.gjf' in listdir('.') and smi + '.out' not in listdir('.')]
        if submitted_smiles:
            self.write_submitted_script(submitted_smiles, 'M062X')


    
    """ To get submitted file of GFNFF. """
    def get_GFNFF_out(self):
        
        # Go to working directory and determine absent molecules.
        chdir('{}/gaussian_out/GFNFF'.format(self.para['work_path']))
        smiles = []
        for smi in self.para['species']:
            if smi + '.dat' not in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], 'M062X')):
                smiles.append(smi)
        
        # Check whether calculations are completed.
        Check(self.para, smiles, 'GFNFF').check_method_out()

        # Construct input files of GFNFF for CREST.
        for smi in smiles:
            if smi + '.xyz' not in listdir('.') and smi + '.out' in listdir('../M062X'):
                with open('../M062X/{}.out'.format(smi)) as f, open(smi + '.xyz', 'w') as p:
                    res = search('temp=200 [\s\S]+?Charge =  (\d) Multiplicity = (\d)[\s\S]+?Standard[\s\S]+?Z\n -+\n((.{70}\n)+?) -+', f.read())
                    coord = []
                    for x in res.group(3).strip().split('\n'):
                        for y in [list(map(float, x.split()))]:
                            coord.append('%s%16.6f%16.6f%16.6f\n'%(Atom(int(y[1])).symbol, y[3], y[4], y[5]))
                    p.write('{}\n\n{}'.format(len(coord), ''.join(coord)))
                    print('This molecule has been converted: {}'.format(smi))       
        
        # Generate submitted scripts.
        submitted_smiles = [smi for smi in smiles if smi + '.xyz' in listdir('.') and smi + '.out' not in listdir('.')]
        if submitted_smiles:
            self.write_submitted_script(submitted_smiles, 'GFNFF')


    
    """ To get submitted file of CCSD(T). """
    def get_CCSDT_out(self):

        # Go to working directory and determine absent molecules.
        chdir(self.para['work_path'] + '/gaussian_out/CCSDT')
        smiles = []
        for smi in self.para['micro_mol']:
            if smi + '.dat' not in listdir('{}/chemkin_dat/CCSDT'.format(self.para['work_path'])):
                smiles.append(smi)
 
        # Check whether calculations are completed.
        Check(self.para, smiles, 'CCSDT').check_method_out()

        # Construct input files of CCSDT for Molpro.
        for smi in smiles:
            if smi + '.out' in listdir('../M062X') and smi + '.gjf' not in listdir('.'):
                with open('../M062X/{}.out'.format(smi)) as f, open('../M062X/{}.gjf'.format(smi)) as p, open(smi + '.gjf', 'w') as q:
                    content = f.read()
                    newsmi = smi[:-2] if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
                    res = search('temp=200 [\s\S]+?Charge =  (\d) Multiplicity = (\d)[\s\S]+?Standard[\s\S]+?Z\n -+\n((.{70}\n)+?) -+', content)
                    spin = int(res.group(2))
                    coord = []
                    for x in res.group(3).strip().split('\n'):
                        for y in [list(map(float, x.split()))]:
                            coord.append('%s%16.6f%16.6f%16.6f\n'%(Atom(int(y[1])).symbol, y[3], y[4], y[5]))
                    atoms = Counter(x.GetSymbol() for x in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAtoms())
                    electrons = sum([self.para['electrons_atoms'][k] * atoms[k] for k in atoms])
                    if spin > 1:
                        q.write('***, {}\n\ngeometry={{\n{}}}\n\nbasis=vdz\n{{rhf;wf,{},spin,{},state,1}}\nuccsd(t)\n\nbasis=vtz\nrhf\nuccsd(t)\n\nbasis=vqz\nrhf\n'.format(smi, ''.join(coord), electrons, spin - 1))
                    else:
                        q.write('***, {}\n\ngeometry={{\n{}}}\n\nbasis=vdz\n{{hf;wf,{},spin,{},state,1}}\nccsd(t)\n\nbasis=vtz\nhf\nccsd(t)\n\nbasis=vqz\nhf\n'.format(smi, ''.join(coord), electrons, spin - 1))
                    print('This molecule has been converted: {}'.format(smi))

        # Generate submitted scripts.
        submitted_smiles = [smi for smi in smiles if smi + '.gjf' in glob('*.gjf') and smi + '.out' not in glob('*.out')]
        if submitted_smiles:
            self.write_submitted_script(submitted_smiles, 'CCSDT')



    """ To get submitted file of G4. """
    def get_G4_out(self, micro_mol):


        # Go to working directory and determine absent molecules.
        chdir(self.para['work_path'] + '/gaussian_out/G4')
        smiles = []
        for smi in self.para['micro_mol']:
            if smi + '.dat' not in listdir('{}/chemkin_dat/G4'.format(self.para['work_path'])):
                smiles.append(smi)
 
        # Check whether calculations are completed.
        Check(self.para, smiles, 'G4').check_method_out()

  
        # Construct input files of G4 for Gaussian.
        for smi in smiles:
            if smi + '.out' in listdir('../M062X') and smi + '.gjf' not in listdir('.'):
                with open('../M062X/{}.out'.format(smi)) as f, open('../M062X/{}.gjf'.format(smi)) as p, open(smi + '.gjf', 'w') as q:
                    res = search('temp=200 [\s\S]+?Charge =  (\d) Multiplicity = (\d)[\s\S]+?Standard[\s\S]+?Z\n -+\n((.{70}\n)+?) -+', f.read())
                    symm = ''.join(findall('\s*(\w*symm=*\w*)', p.read()))
                    coord = []
                    for x in res.group(3).strip().split('\n'):
                        for y in [list(map(float, x.split()))]:
                            coord.append('%s%16.6f%16.6f%16.6f\n'%(Atom(int(y[1])).symbol, y[3], y[4], y[5]))
                    q.write('%nprocshared={}\n# G4=sp scf=fermi {}\n\n{}\n\n0 {}\n{}\n'.format(self.para['number_process'], symm, smi, res.group(2)), ''.join(coord))
                    print('This molecule has been converted: {}'.format(smi))
        
        # Generate submitted scripts.
        submitted_smiles = [smi for smi in smiles if smi + '.gjf' in glob('*.gjf') and smi + '.out' not in glob('*.out')]
        if submitted_smiles:
            self.write_submitted_script(submitted_smiles, 'G4')


    
    """ To get submitted file of CBSQB3. """
    def get_CBSQB3_out(self, micro_mol):
        
         # Go to working directory and determine absent molecules.
        chdir(self.para['work_path'] + '/gaussian_out/CBSQB3')
        smiles = []
        for smi in self.para['micro_mol']:
            if smi + '.dat' not in listdir('{}/chemkin_dat/CBSQB3'.format(self.para['work_path'])):
                smiles.append(smi)
 
        # Check whether calculations are completed.
        Check(self.para, smiles, 'CBSQB3').check_method_out()

  
        # Construct input files of CBSQB3 for Gaussian.
        for smi in smiles:
            if smi + '.out' in listdir('../M062X') and smi + '.gjf' not in listdir('.'):
                with open('../M062X/{}.out'.format(smi)) as f, open('../M062X/{}.gjf'.format(smi)) as p, open(smi + '.gjf', 'w') as q:
                    res = search('temp=200 [\s\S]+?Charge =  (\d) Multiplicity = (\d)[\s\S]+?Standard[\s\S]+?Z\n -+\n((.{70}\n)+?) -+', f.read())
                    symm = ''.join(findall('\s*(\w*symm=*\w*)', p.read()))
                    coord = []
                    for x in res.group(3).strip().split('\n'):
                        for y in [list(map(float, x.split()))]:
                            coord.append('%s%16.6f%16.6f%16.6f\n'%(Atom(int(y[1])).symbol, y[3], y[4], y[5]))
                    q.write('%nprocshared={}\n# G4=sp scf=fermi {}\n\n{}\n\n0 {}\n{}\n'.format(self.para['number_process'], symm, smi, res.group(2)), ''.join(coord))
                    print('This molecule has been converted: {}'.format(smi))

        # Generate submitted scripts.
        submitted_smiles = [smi for smi in species if smi + '.gjf' in glob('*.gjf') and smi + '.out' not in glob('*.out')]
        if submitted_smiles:
            self.write_submitted_script(submitted_smiles, 'CBSQB3')


 
    """ To construct 3D coordinates from SMILES. """
    def build_structure_from_smi(self, smi, para, structure_method = 1):

        # Build conformational structure.
        print('\nThis molecule of structure is being constructed: {}'.format(smi))
        Conformer(smi, self.para, structure_method).output_standard_gjf()

        # Check whether bug files exist.
        with open(smi + '.gjf') as f: content = f.read()
        if '*' in content or len(findall('0.00000         0.00000         0.00000', content)) > 1:
            Conformer(smi, self.para, 2).output_standard_gjf()
    
    
    """ To write submitted scripts for linux systems. """
    def write_submitted_script(self, smiles, method):

        # Generate the submitted directory.
        if self.para['submitted_type'] not in [1, 2]:
            print('!!! This type has not been considered. Default type of nohup has been applied. !!!')
        gjf=["'{}.gjf'".format(smi) if set(smi) & set(['(', '[', '=']) else smi + '.gjf' for smi in smiles]
        makedirs('{}/submitted_inp/{}'.format(self.para['work_path'], method), True)
        
        # Generate sumbmitted scripts for CREST.
        if method == 'GFNFF':
            chunks = lambda a, b: [a[i:i+b] for i in range(0, len(a), b)]
            xyz=["'{}.xyz'".format(smi) if set(smi) & set(['(', '[', '=']) else smi + '.xyz' for smi in smiles]
            with open('{}/submitted_inp/GFNFF/{}'.format(self.para['work_path'], self.para['submitted_shell']), 'w', newline = '\n') as f:
                f.write('#!/bin/bash\n\n')
                for i, v in enumerate(chunks(sample(xyz, len(xyz)), ceil(len(xyz) / int(self.para['number_task'])))):
                    with open('{}/submitted_inp/GFNFF/job{}.sh'.format(self.para['work_path'], i + 1), 'w', newline = '\n') as p:
                        p.write('#!/bin/bash\n\n')
                        if self.para['submitted_shell'] == 2:
                            f.write('mkdir {0}\ncp job{0}.sh {0}\ncp {1} {0}\ncd {0}\nqsub *.sh\ncd ..\n\n'.format(i + 1, ' '.join(v)))
                            p.write('#PBS -N {}-job{}\n#PBS -l nodes=1:ppn={}\n#PBS -q medium\n#PBS -j eo\ncd $PBS_O_WORKDIR\n\n'.format(method, i + 1, self.para['number_process']))
                            for x in v:
                                smi = x[:-4] if x[0] != '\'' else x[1:-5]
                                with open('../M062X/' + smi + '.gjf') as q: chrg, spin = map(int, findall('\n\n(\d)\s+(\d)\s*\n', q.read())[0])
                                p.write('crest {0} >{1} -gfnff -v4 -quick -ewin 50 -T {2} -chrg {3} -uhf {4}\nwait\ncp -f crest_best.xyz ../{0}\ncp {1} ..\n'.format(x, x.replace('xyz', 'out'), self.para['number_process'], chrg, spin - 1))
                        else:
                            f.write('mkdir {0}\ncp job{0}.sh {0}\ncp {1} {0}\ncd {0}\nchmod +x *.sh\nnohup ./*.sh >nohup.out &\ncd ..\n\n'.format(i + 1, ' '.join(v)))
                            for x in v:
                                smi = x[:-4] if x[0] != '\'' else x[1:-5]
                                with open('../M062X/' + smi + '.gjf') as q: chrg, spin = map(int, findall('\n\n(\d)\s+(\d)\s*\n', q.read())[0])
                                p.write('nohup crest {0} >{1} -gfnff -v4 -quick -ewin 50 -T {2} -chrg {3} -uhf {4} &\nwait\ncp -f crest_best.xyz ../{0}\ncp {1} ..\n'.format(x, x.replace('xyz', 'out'), self.para['number_process'], chrg, spin - 1))
                if self.para['submitted_type'] != 2: f.write('wait\nrm -rf {}\n\n'.format(' '.join(map(str, range(1, i + 2)))))
        
        # Generate sumbmitted scripts for Molpro.
        elif method == 'CCSDT':
            chunks = lambda a, b: [a[i:i+b] for i in range(0, len(a), b)]
            for i, v in enumerate(chunks(sample(gjf, len(gjf)), ceil(len(gjf) / int(self.para['number_task'])))):
                with open('{}/submitted_inp/{}/job{}.sh'.format(self.para['work_path'], method, i + 1), 'w', newline = '\n') as f:
                    f.write('#!/bin/bash\n\n')
                    if self.para['submitted_type'] == 2:
                        f.write('#PBS -N {}-job{}\n#PBS -l nodes=1:ppn={}\n#PBS -q medium\n#PBS -j eo\ncd $PBS_O_WORKDIR\n\n'.format(method, i + 1, self.para['number_process']))
                        for x in v: f.write('molpro -n {} -m 300m <{} >{}\nwait\nsleep 10'.format(self.para['number_process'], x, x.replace('gjf', 'out')))
                    else:
                        for x in v: f.write('nohup molpro -n {} -m 300m <{} >{} &\nwait\nsleep 10\n'.format(self.para['number_process'], x, x.replace('gjf', 'out')))
        
        # Generate sumbmitted scripts for Gaussian.
        else:
            chunks = lambda a, b: [a[i:i+b] for i in range(0, len(a), b)]
            for i, v in enumerate(chunks(sample(gjf, len(gjf)), ceil(len(gjf) / int(self.para['number_process'])))):
                with open('{}/submitted_inp/{}/job{}.sh'.format(self.para['work_path'], method, i + 1), 'w', newline = '\n') as f:
                    f.write('#!/bin/bash\n\n')
                    if self.para['submitted_type'] == '2':
                        f.write('#PBS -N {}-job{}\n#PBS -l nodes=1:ppn={}\n#PBS -q medium\n#PBS -j eo\ncd $PBS_O_WORKDIR\n\n'.format(method, i + 1, self.para['number_process']))
                        for x in v: f.write('g09 <{} >{}\nwait\n'.format(x, x.replace('gjf', 'out')))
                    else:
                        for x in v: f.write('nohup g09 <{} >{} &\nwait\n'.format(x, x.replace('gjf', 'out')))
        
        # Prepare related files to the submitted directory.
        if method == 'GFNFF':
            for smi in smiles:
                copy(smi + '.xyz', '{}/submitted_inp/{}'.format(self.para['work_path'], method))
        else:
            with open('{}/submitted_inp/{}/{}'.format(self.para['work_path'], method, self.para['submitted_shell']), 'w', newline = '\n') as f:
                if self.para['submitted_type'] == 2:
                    write_content = '  qsub $s'
                else:
                    write_content = '  chmod +x $s\n  nohup ./$s &'
                f.write('#!/bin/bash\n\nfor s in *.sh;\ndo\n  dos2unix $s\n{}\ndone\n'.format(write_content))
            for smi in smiles:
                copy(smi + '.gjf', '{}/submitted_inp/{}'.format(self.para['work_path'], method))



    """ To construct submitted files for Linux system. """
    def get_submitted_out(self):
        print('\nProcess submited files ...\n')
        self.get_M062X_out()
        self.get_GFNFF_out()
        exec('self.get_{}_out()'.format(self.para['high_level']))
