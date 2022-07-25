# -*- coding: utf-8 -*-
from glob import glob
from rdkit import Chem
from re import findall
from os.path import exists
from shutil import copy, rmtree 
from collections import Counter
from os import chdir, listdir, makedirs

from tdoc._CBH import CBH



class Parameters():
    
    """ To initiate parameters for Parameters. """
    def __init__(self,input_file,para_file, work_path):
        self.para={'input_file':input_file,'para_file':para_file, 'work_path':work_path}
        self.get_input_parameters()
        
    """ To get SMILES. """
    def get_smiles(self):
        smiles, smilesdict = [], {}
        if 'program_out' not in listdir(self.para['work_path']):
            makedirs('{}/program_out'.format(self.para['work_path']))
        with open('{}/{}'.format(self.para['work_path'], self.para['input_file'])) as f, open('{}/program_out/canonical_smiles.txt'.format(self.para['work_path']), 'w') as p:
            p.write('{:>30}{:>30}{:>30}\n'.format('species', 'smiles', 'canonical smiles'))
            for x in f:
                if x.strip():
                    species, smi = x.strip().split()
                    newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
                    mol = Chem.MolFromSmiles(newsmi)
                    Chem.Kekulize(mol, clearAromaticFlags = True)
                    cal_smi = Chem.MolToSmiles(mol) + smi[-2:] if '-' in smi and smi[-1] in self.para['multiplicities'] else Chem.MolToSmiles(mol)
                    p.write('{:>30}{:>30}{:>30}\n'.format(species, smi, cal_smi))
                    smilesdict.update({species: cal_smi})
        smiles = set(smilesdict.values())
        return smiles, smilesdict
    
    """ To identify SMILES and classify them into macro and micro molecules. """
    def classify_smiles(self,species):
        print('\nClassify smiles ...\n')
        micro_mol, macro_mol, CBH_equation = set(), set(), {}
        for smi in species:
            newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
            reac, prod = CBH(self.para['aro_rings']).get_CBH3(newsmi)
            if reac != prod:
                if findall('=|#', smi):
                    mol = Chem.MolFromSmiles(newsmi)
                    if not mol.GetAromaticAtoms():
                        atoms = Counter(x.GetSymbol() for x in Chem.AddHs(mol).GetAtoms())
                        if 0.5 * (atoms.get('C', 0) * 2 + 2 - atoms.get('H', 0)) > atoms.get('C', 0) + atoms.get('O', 0) - 2:
                            micro_mol.add(smi), CBH_equation.update({smi:[Counter([smi]), Counter([smi])]})
                        else:
                            macro_mol.add(smi), micro_mol.update(set(reac + prod - Counter([newsmi]))), CBH_equation.update({smi: [reac, prod]})
                    else:
                        macro_mol.add(smi), micro_mol.update(set(reac + prod - Counter([newsmi]))), CBH_equation.update({smi: [reac, prod]})
                else:
                    macro_mol.add(smi), micro_mol.update(set(reac + prod - Counter([newsmi]))), CBH_equation.update({smi: [reac, prod]})
            else: micro_mol.add(smi), CBH_equation.update({smi: [Counter([smi]), Counter([smi])]})
        return micro_mol, macro_mol, CBH_equation
  
    """ To get input parameters. """
    def get_input_parameters(self):
        
        def input_parameters(**args):
            self.para.update(args)

        def default_parameters(**args):
            self.para.update(args)

        with open('{}/{}'.format(self.para['work_path'], self.para['para_file'])) as f:
            exec(f.read())
    
    """ To build working directories. """
    def build_work_dir(self):
        for x in ['GFNFF', 'M062X', self.para['high_level']]:
            for y in ['gaussian_out', 'chemkin_dat']:
                if not exists('{}/{}/{}'.format(self.para['work_path'], y, x)):
                    makedirs('{}/{}/{}'.format(self.para['work_path'], y, x))
                rmtree('{}/chemkin_dat/GFNFF'.format(self.para['work_path']), True)
            if exists('{}/preexisted_dat/{}'.format(self.para['work_path'], x)):
                for z in glob('{}/preexisted_dat/{}/*.dat'.format(self.para['work_path'], x)):
                    copy(z, '{}/chemkin_dat/{}'.format(self.para['work_path'], x))
        if exists('{}/submitted_inp'.format(self.para['work_path'])):
            rmtree('{}/submitted_inp'.format(self.para['work_path']),  True)
    
    """ To get all species. """
    def get_all_species(self):
        species = set()
        for smi in self.get_smiles()[0]:
            if smi + '.dat' not in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'],self.para['high_level'])):
                species.add(smi)
        micro_mol, macro_mol, _ = self.classify_smiles(species)
        species = micro_mol | macro_mol
        return species, micro_mol, macro_mol
    
    """ To get all parameters. """     
    def get_all_parameters(self):
        self.para['high_level']=self.para['calcualted_levels'][self.para['calculated_method']]
        self.para['submitted_shell']=self.para['submitted_scripts'][self.para['submitted_type']]
        self.build_work_dir()
        self.para.update(dict(zip(['species','micro_mol','macro_mol'],self.get_all_species())))
        return self.para
