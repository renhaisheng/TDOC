# -*- coding: utf-8 -*-
from rdkit import Chem
from collections import Counter
from rdkit.Chem import AllChem, rdMolTransforms



class CBH(object):



    """ To initiate parameters for CBH. """
    def __init__(self, para):
        self.para = para



    """ To get atomic information. """
    def get_atom_info(self,mol,atom):
        p_hybrid_degree = str(atom.GetHybridization())[-1] if str(atom.GetHybridization())[-1] in ["2", "3"] else "1"
        adjacent_heavy_atoms = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() != 'H']
        adjacent_H_atoms = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() == 'H']
        atom_info = '{}{}{}{}'.format(atom.GetSymbol(), len(adjacent_heavy_atoms), len(adjacent_H_atoms), p_hybrid_degree) if atom.GetSymbol() != 'H' else 'H'
        return atom_info



    """ To get right term of CBH0. """
    def CBH0_Counter(self, smi):
        newsmi = smi[4:] if smi.startswith('cis-') else smi
        smiles, mol = [], Chem.MolFromSmiles(newsmi)
        if mol.GetNumHeavyAtoms() <= 1: return Counter([smi])
        for x in mol.GetAtoms():
            smiles.append(x.GetSymbol())
        return Counter(smiles)
    


    """ To get right term of CBH1. """
    def CBH1_Counter(self, smi):
        newsmi = smi[4:] if smi.startswith('cis-') else smi
        smiles, mol = [], Chem.MolFromSmiles(newsmi)
        if mol.GetNumHeavyAtoms() <= 2: return Counter([smi])
        if mol.GetAromaticAtoms(): Chem.Kekulize(mol, clearAromaticFlags = True)
        for bond in mol.GetBonds():
            leftatom = bond.GetBeginAtom()
            rightatom = bond.GetEndAtom()
            newmol = Chem.RWMol(Chem.MolFromSmiles(''))
            newmol.AddAtom(leftatom)
            newmol.AddAtom(rightatom)
            newmol.AddBond(0, 1, bond.GetBondType())
            smiles.append(Chem.MolToSmiles(newmol))
        if not smiles: smiles = [smi]
        return Counter(smiles)
   


    """ To get right term of CBH2. """
    def CBH2_Counter(self, smi):
        newsmi = smi[4:] if smi.startswith('cis-') else smi
        smiles, mol = [], Chem.MolFromSmiles(newsmi)
        if mol.GetNumHeavyAtoms() <= 3: return Counter([smi])
        if mol.GetAromaticAtoms(): Chem.Kekulize(mol, clearAromaticFlags = True)
        for x in mol.GetAtoms():
            newmol = Chem.RWMol(Chem.MolFromSmiles(''))
            if x.GetDegree() == 1: continue
            newmol.AddAtom(x)
            for i, y in enumerate(x.GetBonds()):
                newmol.AddAtom(y.GetOtherAtom(x))
                newmol.AddBond(0, i + 1, y.GetBondType())
            smiles.append(Chem.MolToSmiles(newmol))
        if not smiles: smiles = [smi]
        return Counter(smiles)



    """ To get right term of CBH3. """
    def CBH3_Counter(self, smi):
        newsmi=smi[4:] if smi.startswith('cis-') else smi
        smiles, mol = [], Chem.MolFromSmiles(newsmi)
        AllChem.Compute2DCoords(mol)
        if mol.GetNumHeavyAtoms() <= 4: return Counter([smi])
        if mol.GetAromaticAtoms(): Chem.Kekulize(mol, clearAromaticFlags = True)
        for bond in mol.GetBonds():
            cis_tag=''
            leftatom = bond.GetBeginAtom()
            rightatom = bond.GetEndAtom()
            newmol = Chem.RWMol(Chem.MolFromSmiles(''))
            if leftatom.GetDegree() == 1 or rightatom.GetDegree() == 1: continue
            leftheavyatoms=[x.GetIdx() for x in leftatom.GetNeighbors() if x.GetIdx()!=rightatom.GetIdx()]
            rightheavyatoms=[x.GetIdx() for x in rightatom.GetNeighbors() if x.GetIdx()!=leftatom.GetIdx()]
            leftbonds=[mol.GetBondBetweenAtoms(leftatom.GetIdx(), atomidx) for atomidx in leftheavyatoms]
            rightbonds=[mol.GetBondBetweenAtoms(rightatom.GetIdx(), atomidx) for atomidx in rightheavyatoms]
            leftoutside = sorted([(self.get_atom_info(mol, x.GetOtherAtom(leftatom)), self.para['bond_marks'][x.GetBondTypeAsDouble()], x.GetOtherAtom(leftatom).GetIdx()) for x in leftbonds])
            rightoutside=sorted([(self.get_atom_info(mol, x.GetOtherAtom(rightatom)), self.para['bond_marks'][x.GetBondTypeAsDouble()], x.GetOtherAtom(rightatom).GetIdx()) for x in rightbonds])
            if leftheavyatoms and rightheavyatoms:
                if [x.GetIsConjugated() for x in leftbonds+rightbonds].count(True)>=1 and bond.GetIsConjugated() or bond.GetBondTypeAsDouble()==2.0:
                    dihedral_angle=abs(rdMolTransforms.GetDihedralDeg(mol.GetConformer(), leftoutside[0][2], leftatom.GetIdx(),rightatom.GetIdx(),rightoutside[0][2]))
                    cis_tag='cis-' if dihedral_angle<90 else ''
            newmol.AddAtom(leftatom)
            newmol.AddAtom(rightatom)
            newmol.AddBond(0, 1, bond.GetBondType())
            for i, y in enumerate(m for m in leftatom.GetBonds() if m.GetIdx() != bond.GetIdx()):
                newmol.AddAtom(y.GetOtherAtom(leftatom))
                newmol.AddBond(0, i + 2, y.GetBondType())
            for j, z in enumerate(m for m in rightatom.GetBonds() if m.GetIdx() != bond.GetIdx()):
                newmol.AddAtom(z.GetOtherAtom(rightatom))
                newmol.AddBond(1, i + j + 3, z.GetBondType())
            smiles.append(cis_tag+Chem.MolToSmiles(newmol))
        if not smiles: smiles = [smi]
        return Counter(smiles)



    """ To get all species of CBH1. """   
    def get_CBH1(self, smi):
        left, right, CBH1 = self.CBH0_Counter(smi), Counter({}), self.CBH1_Counter(smi)
        for k, v in CBH1.items(): right += Counter({x: y * v for x, y in self.CBH0_Counter(k).items()})
        reac = Counter([smi]) + right - left
        return reac, CBH1
    


    """ To get all species of CBH2. """ 
    def get_CBH2(self, smi):
        left, right, CBH2 = self.CBH1_Counter(smi), Counter({}), self.CBH2_Counter(smi)
        for k, v in CBH2.items(): right += Counter({x: y * v for x, y in self.CBH1_Counter(k).items()})
        reac = Counter([smi]) + right - left
        return reac, CBH2
  


    """ To get all species of CBH3. """ 
    def get_CBH3(self, smi):
        left, right, CBH3 = self.CBH2_Counter(smi), Counter({}), self.CBH3_Counter(smi)
        for k, v in CBH3.items(): right += Counter({x: y * v for x, y in self.CBH2_Counter(k).items()})
        reac = Counter([smi]) + right - left
        return reac, CBH3
