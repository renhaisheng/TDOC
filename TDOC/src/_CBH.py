# -*- coding: utf-8 -*-
from collections import Counter
from rdkit import Chem
from itertools import chain
from rdkit.Chem import AllChem

class CBH(object):

    def __init__(self,aro_rings=['C12=CC=CC=C1C=CC=C2','C1=CC=CC=C1','C1=COC=C1']):
        self.replaced_unit=['S','P','As']
        self.aro_rings=aro_rings
    
    def replace_mol(self, mol):
        Chem.Kekulize(mol, clearAromaticFlags = True)
        for i, v in enumerate(self.aro_rings):
            smarts = Chem.MolFromSmarts(v)
            units = Chem.MolFromSmarts(self.replaced_unit[i])
            mol = AllChem.ReplaceSubstructs(mol, smarts , units, replaceAll = True)[0]
        return mol

    def recover_mol(self, mol):
        Chem.Kekulize(mol, clearAromaticFlags = True)
        for i, v in enumerate(self.aro_rings):
            smarts = Chem.MolFromSmarts(v)
            units = Chem.MolFromSmarts(self.replaced_unit[i])
            mol = AllChem.ReplaceSubstructs(mol, units , smarts, replaceAll = True)[0]
        return mol

    def CBH0_Counter(self, smi):
        smiles, mol = [], Chem.MolFromSmiles(smi)
        if mol.GetNumHeavyAtoms() <= 1:
            return Counter([smi])
        mol = self.replace_mol(mol)
        for x in mol.GetAtoms():
            smiles.append(x.GetSymbol())
        for x in set(smiles)&set(self.replaced_unit):
            mol = Chem.MolFromSmiles(x)
            mol = self.recover_mol(mol)
            newsmi = Chem.MolToSmiles(mol)
            smiles = '\n'.join(smiles).replace(x, newsmi).split('\n')
        return Counter(smiles)

    def CBH1_Counter(self, smi, ):
        smiles, mol = [], Chem.RWMol(Chem.MolFromSmiles(smi))
        if mol.GetNumHeavyAtoms() <= 2:
            return Counter([smi])
        mol = self.replace_mol(mol)
        for bond in mol.GetBonds():
            leftatom = bond.GetBeginAtom()
            rightatom = bond.GetEndAtom()
            newmol = Chem.RWMol(Chem.MolFromSmiles(''))
            newmol.AddAtom(leftatom)
            newmol.AddAtom(rightatom)
            newmol.AddBond(0, 1, bond.GetBondType())
            newmol = self.recover_mol(newmol)
            smiles.append(Chem.MolToSmiles(newmol))
        if not smiles: smiles = [smi]
        return Counter(smiles)

    def CBH2_Counter(self, smi):
        smiles, mol = [], Chem.RWMol(Chem.MolFromSmiles(smi))
        if mol.GetNumHeavyAtoms() <= 3: return Counter([smi])
        mol = self.replace_mol(mol)
        for x in mol.GetAtoms():
            newmol = Chem.RWMol(Chem.MolFromSmiles(''))
            if x.GetDegree() == 1: continue
            newmol.AddAtom(x)
            for i, y in enumerate(x.GetBonds()):
                newmol.AddAtom(y.GetOtherAtom(x))
                newmol.AddBond(0, i + 1, y.GetBondType())
            newmol = self.recover_mol(newmol)
            smiles.append(Chem.MolToSmiles(newmol))
        return Counter(smiles)

    def CBH3_Counter(self, smi):
        smiles, mol = [], Chem.RWMol(Chem.MolFromSmiles(smi))
        if mol.GetNumHeavyAtoms() <= 4:
            return Counter([smi])
        mol = self.replace_mol(mol)
        for bond in mol.GetBonds():
            leftatom = bond.GetBeginAtom()
            rightatom = bond.GetEndAtom()
            newmol = Chem.RWMol(Chem.MolFromSmiles(''))
            if leftatom.GetDegree() == 1 or rightatom.GetDegree() == 1: continue
            newmol.AddAtom(leftatom)
            newmol.AddAtom(rightatom)
            newmol.AddBond(0, 1, bond.GetBondType())
            for i, y in enumerate(m for m in leftatom.GetBonds() if m.GetIdx() != bond.GetIdx()):
                newmol.AddAtom(y.GetOtherAtom(leftatom))
                newmol.AddBond(0, i + 2, y.GetBondType())
            for j, z in enumerate(m for m in rightatom.GetBonds() if m.GetIdx() != bond.GetIdx()):
                newmol.AddAtom(z.GetOtherAtom(rightatom))
                newmol.AddBond(1, i + j + 3, z.GetBondType())
            newmol = self.recover_mol(newmol)
            smiles.append(Chem.MolToSmiles(newmol))
        if not smiles: smiles = [smi]
        return Counter(smiles)

        
    def Get_CBH1(self, smi):
        left, right, CBH1 = self.CBH0_Counter(smi), Counter({}), self.CBH1_Counter(smi)
        for k, v in CBH1.items(): right += Counter({x: y * v for x, y in self.CBH0_Counter(k).items()})
        reac = Counter([smi]) + right - left
        return reac, CBH1

    def Get_CBH2(self, smi):
        left, right, CBH2 = self.CBH1_Counter(smi), Counter({}), self.CBH2_Counter(smi)
        for k, v in CBH2.items(): right += Counter({x: y * v for x, y in self.CBH1_Counter(k).items()})
        reac = Counter([smi]) + right - left
        return reac, CBH2

    def Get_CBH3(self, smi):
        left, right, CBH3 = self.CBH2_Counter(smi), Counter({}), self.CBH3_Counter(smi)
        for k, v in CBH3.items(): right += Counter({x: y * v for x, y in self.CBH2_Counter(k).items()})
        reac = Counter([smi]) + right - left
        return reac, CBH3

