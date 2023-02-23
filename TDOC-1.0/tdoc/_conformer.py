# -*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel



class Conformer(object):

    """ To initiate parameters for Conformer. """
    def __init__(self, smi, para, structure_method):
        self.smi = smi
        self.para = para
        self.structure_method = structure_method



    """ To generate the gjf file from rdkit. """
    def rdkit_to_gjf(self):
        newsmi = self.smi[:-2] if '-' in self.smi and self.smi[-1] in self.para['multiplicities'] else self.smi
        mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        molblock = Chem.MolToMolBlock(mol)
        mol = pybel.readstring('mol',molblock)
        gjf = mol.write('gjf').split('\n')
        gjf[3] = self.smi + '_' + str(self.structure_method)
        if '-' in self.smi and self.smi[-1] in self.para['multiplicities']:
            gjf[5] = gjf[5][:-1] + str(self.structure_method)
        gjf = '\n'.join(gjf)
        return gjf



    """ To generate the gjf file from openbabel. """
    def pybel_to_gjf(self):
        spin = self.rdkit_to_gjf().split('\n')[5][-1]
        newsmi = self.smi[:-2] if '-' in self.smi and self.smi[-1] in self.para['multiplicities'] else self.smi
        mol = pybel.readstring('smi', newsmi)
        mol.make3D()
        gjf = mol.write('gjf').split('\n')
        gjf[3] = self.smi + '_' + str(self.structure_method)
        gjf[5] = gjf[5][:-1] + spin
        gjf = '\n'.join(gjf)
        return gjf

    

    """ To generate the standard gjf file for Gaussian. """
    def output_standard_gjf(self):
        
        # Get 3D structure.
        if self.structure_method == 2:
            gjf = self.rdkit_to_gjf()
        else:
            gjf = self.pybel_to_gjf()
        
        # Produce standard input gjf file.
        if 'press' not in self.para['key_words']:
            key_words = self.para['key_words'] + ' press=0.987'
        link_words = '# freq=(readfc,hinderedrotor)' if 'hinderedrotor' in key_words else '# freq=readfc'
        with open(self.smi + '.gjf', 'w') as f:
            f.write('%nprocshared={}\n%chk={}.chk\n{}\n'.format(self.para['number_process'], self.smi, key_words))
            f.write('\n'.join(gjf.split('\n')[2:]))
            for T in list(range(200, 1000, 50)) + list(range(1000, 5050, 100)):
                f.write('\n--Link1--\n%chk={}.chk\n{} temp={} press=0.987 geom=allcheck guess=read\n'.format(self.smi, link_words, T))
