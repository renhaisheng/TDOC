# -*- coding: utf-8 -*-
from os import remove, chdir, getcwd, listdir, makedirs
from subprocess import Popen, PIPE
from collections import Counter
from rdkit.Chem import AllChem
from shutil import copy, rmtree 
from os.path import exists
from rdkit import Chem
from re import findall
from time import sleep
from glob import glob
from ._CBH import CBH


class PROJECT(object):

    def __init__(self):
        self.submit_scripts = {'1': 'job-nohup', '2': 'job-qsub', '3': 'job-bsub'}    # The submited scriptes of linux. 
        self.multiplicities = {'S': '1', 'D': '2', 'T': '3', 'Q': '4', 'P': '5'}    # The Symbol of spin multiplicities.
        self.high_precision = {'1': 'CCSDT', '2': 'G4', '3': 'CBSQB3'}    # The calculated method of high precision.
        self.path = getcwd()
        self.number_task = '4'
        self.number_process = '8'
        self.submit_type = '1'
        self.calculated_method = '1'
        self.conformer_method = '2'
        self.spin_method = '1'
        self.key_words = '# opt freq=hinderedrotor M062X/6-311++g(d,p)'
        self.Parameters()
            
    # To get SMILES.
    def Get_smiles(self):
        chdir(self.path)
        smiles, smilesdict = [], {}
        if 'smiles.txt' in listdir(self.path):
            if 'program_out' not in listdir(self.path): makedirs('program_out')
            with open('smiles.txt') as f, open('program_out/canonical_smiles.txt', 'w') as p:
                p.write('{:>30}{:>30}{:>30}\n'.format('species', 'smiles', 'canonical smiles'))
                for x in f:
                    if x.strip():
                        species, smi = x.strip().split()
                        newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.multiplicities else smi
                        mol = Chem.MolFromSmiles(newsmi)
                        Chem.Kekulize(mol, clearAromaticFlags = True)
                        cal_smi = Chem.MolToSmiles(mol) + smi[-2:] if '-' in smi and smi[-1] in self.multiplicities else Chem.MolToSmiles(mol)
                        p.write('{:>30}{:>30}{:>30}\n'.format(species, smi, cal_smi))
                        smilesdict.update({species: cal_smi})
        else: print('!!! Input files are missing. !!!'), sleep(5), exit()
        if 'chem.inp' in listdir(self.path):
            with open('chem.inp') as f:
                species = ''.join(findall('\n\s*\w.+', findall('SPECIES([\s\S]+?)REACTIONS', f.read())[0])).split()
                for mol in species:
                    if mol in smilesdict: smiles.append(smilesdict[mol])
                    else: print('!!! This molecule is not in smiles.dict: {}. !!!'.format(mol)), sleep(5), exit() 
        else: smiles = set(smilesdict.values())
        return smiles, smilesdict
    
    # To identify SMILES and classify them into macro and micro molecules.  
    def Classify_smiles(self):
        print('\nClassify smiles ...\n')
        micro_mol, macro_mol, CBH_equation = set(), set(), {}
        for smi in self.species:
            newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.multiplicities else smi
            reac, prod = CBH().Get_CBH3(newsmi)
            if Chem.MolFromSmiles(newsmi).GetAromaticAtoms(): micro_mol.add(smi), CBH_equation.update({smi:[Counter([smi]), Counter([smi])]})
            else:
                reac, prod = CBH().Get_CBH3(newsmi)
                if reac != prod:
                    if findall('=|#', smi):
                        cmd = 'obabel -:"{}" -O obabel.sdf --confab --gen3D'.format(newsmi)
                        p = Popen(cmd, shell = True, stdout = PIPE, stderr = PIPE)
                        res = findall('tot conformations = (\d+)', p.communicate()[0].decode())
                        if not res and len(findall('=|#',  findall('1(.+)1', smi if '1' in smi else '1=1')[0])) > 1:
                            micro_mol.add(smi), CBH_equation.update({smi:[Counter([smi]), Counter([smi])]})
                        elif len(findall('=|#', smi)) > 1:
                            atoms = Counter(x.GetSymbol() for x in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAtoms())
                            if 0.5 * (atoms.get('C', 0) * 2 + 2 - atoms.get('H', 0)) > atoms.get('C', 0) + atoms.get('O', 0) - 3:
                                micro_mol.add(smi), CBH_equation.update({smi:[Counter([smi]), Counter([smi])]})
                            else:
                                macro_mol.add(smi), micro_mol.update(set(reac + prod - Counter([newsmi]))), CBH_equation.update({smi: [reac, prod]})
                        else:
                            macro_mol.add(smi), micro_mol.update(set(reac + prod - Counter([newsmi]))), CBH_equation.update({smi: [reac, prod]})
                    else:
                        macro_mol.add(smi), micro_mol.update(set(reac + prod - Counter([newsmi]))), CBH_equation.update({smi: [reac, prod]})
                else: micro_mol.add(smi), CBH_equation.update({smi: [Counter([smi]), Counter([smi])]})
        if 'obabel.sdf' in listdir('.'): remove('obabel.sdf')
        return micro_mol, macro_mol, CBH_equation

    def Parameters(self):
        if 'input.ini' in listdir(self.path):
            with open(self.path + '/input.ini') as f: self.number_task, self.number_process, self.submit_type, self.calculated_method, self.conformer_method, spin_method, key_words = findall('=\s*(.+?)\s*@', f.read())
        else: print('!!! The input.ini is missing and the defualt values will be used. !!!'), sleep(5), exit()
        self.high_level = self.high_precision[self.calculated_method]
        for x in ['GFNFF', 'M062X', self.high_level]:
            for y in ['gaussian_out', 'chemkin_dat']:
                if not exists('{}/{}/{}'.format(self.path, y, x)): makedirs('{}/{}/{}'.format(self.path, y, x))
                rmtree('{}/chemkin_dat/GFNFF'.format(self.path), True)
            if exists('{}/preexisted_dat/{}'.format(self.path, x)):
                for z in glob('{}/preexisted_dat/{}/*.dat'.format(self.path, x)): copy(z, '{}/chemkin_dat/{}'.format(self.path, x))
        if exists('{}/submitted_inp'.format(self.path)): rmtree('{}/submitted_inp'.format(self.path),  True)
        self.species = set([smi for smi in self.Get_smiles()[0] if smi + '.dat' not in listdir('{}/chemkin_dat/{}'.format(self.path, self.high_level))])
        self.micro_mol, self.macro_mol, _ = self.Classify_smiles()
        self.species = self.micro_mol | self.macro_mol
