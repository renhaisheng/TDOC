# -*- coding: utf-8 -*-
from rdkit import Chem
from openbabel import pybel
from rdkit.Chem import AllChem, rdMolTransforms


class Conformer(object):

    """ To initiate parameters for Conformer. """
    def __init__(self, smi, para, structure_method):
        self.smi = smi
        self.para = para
        self.structure_method = structure_method



    def get_atom_info(self,mol,atom):
        p_hybrid_degree=str(atom.GetHybridization())[-1] if str(atom.GetHybridization())[-1] in ["2","3"] else "1"
        adjacent_heavy_atoms=[x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol()!='H']
        adjacent_H_atoms=[x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol()=='H']
        atom_info='{}{}{}{}'.format(atom.GetSymbol(),len(adjacent_heavy_atoms),len(adjacent_H_atoms),p_hybrid_degree) if atom.GetSymbol()!='H' else 'H'
        return atom_info



    def get_cis_mol(self):
        bondmark = {1.0:'-', 2.0:'=', 3.0:'#', 1.5:'$'}
        newsmi = self.smi[4:]
        mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))
        AllChem.Compute2DCoords(mol)
        AllChem.EmbedMolecule(mol)
        for bond in mol.GetBonds():
            leftatom = bond.GetBeginAtom()
            rightatom = bond.GetEndAtom()
            newmol = Chem.RWMol(Chem.MolFromSmiles(''))
            if leftatom.GetDegree() == 1 or rightatom.GetDegree() == 1: continue
            leftheavyatoms = [x.GetIdx() for x in leftatom.GetNeighbors() if x.GetSymbol() != 'H' and x.GetIdx() != rightatom.GetIdx()]
            rightheavyatoms = [x.GetIdx() for x in rightatom.GetNeighbors() if x.GetSymbol() != 'H' and x.GetIdx() != leftatom.GetIdx()]
            leftbonds = [mol.GetBondBetweenAtoms(leftatom.GetIdx(),atomidx) for atomidx in leftheavyatoms]
            rightbonds = [mol.GetBondBetweenAtoms(rightatom.GetIdx(),atomidx) for atomidx in rightheavyatoms]
            leftoutside = sorted([(self.get_atom_info(mol,x.GetOtherAtom(leftatom)), bondmark[x.GetBondTypeAsDouble()], x.GetOtherAtom(leftatom).GetIdx()) for x in leftbonds])
            rightoutside = sorted([(self.get_atom_info(mol,x.GetOtherAtom(rightatom)), bondmark[x.GetBondTypeAsDouble()], x.GetOtherAtom(rightatom).GetIdx()) for x in rightbonds])
            if leftheavyatoms and rightheavyatoms:
                if [x.GetIsConjugated() for x in leftbonds+rightbonds].count(True) >= 1 and bond.GetIsConjugated() or bond.GetBondTypeAsDouble() > 1.0:
                    rdMolTransforms.SetDihedralDeg(mol.GetConformer(), leftoutside[0][2], leftatom.GetIdx(),rightatom.GetIdx(),rightoutside[0][2], 0)
                    molblock = Chem.MolToMolBlock(mol)
                    return molblock



    """ To generate the gjf file from rdkit. """
    def rdkit_to_gjf(self):
        if self.smi.startswith('cis-'):
            molblock = self.get_cis_mol()
        else:
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
        if self.smi.startswith('cis-'):
            gjf = self.rdkit_to_gjf()
        else:
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
            for T in range(200, 3050, 50):
                f.write('\n--Link1--\n%chk={}.chk\n{} temp={} press=0.987 geom=allcheck guess=read\n'.format(self.smi, link_words, T))
