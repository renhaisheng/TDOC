# -*- coding: utf-8 -*-
from collections import Counter
from rdkit import Chem


class CBH(object):
    def CBH0_Counter(self, smi):
        smiles, mol = [], Chem.MolFromSmiles(smi)
        if mol.GetNumHeavyAtoms() <= 1: return Counter([smi])
        for x in mol.GetAtoms(): smiles.append(x.GetSymbol())
        return Counter(smiles)

    def CBH1_Counter(self, smi):
        smiles, mol = [], Chem.RWMol(Chem.MolFromSmiles(smi))
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

    def CBH2_Counter(self, smi):
        smiles, mol = [], Chem.RWMol(Chem.MolFromSmiles(smi))
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
        return Counter(smiles)

    def CBH3_Counter(self, smi):
        smiles, mol = [], Chem.RWMol(Chem.MolFromSmiles(smi))
        if mol.GetNumHeavyAtoms() <= 4: return Counter([smi])
        if mol.GetAromaticAtoms(): Chem.Kekulize(mol, clearAromaticFlags = True)
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
            smiles.append(Chem.MolToSmiles(newmol))
        if not smiles: smiles = [smi]
        return Counter(smiles)

    def CBH4_Counter(self, smi):
        smiles, mol = [], Chem.RWMol(Chem.MolFromSmiles(smi))
        if mol.GetAromaticAtoms(): Chem.Kekulize(mol, clearAromaticFlags = True)
        for x in mol.GetAtoms():
            newmol = Chem.RWMol(Chem.MolFromSmiles(''))
            if x.GetDegree() == 1 or all(m.GetDegree() == 1 for m in x.GetNeighbors()): continue
            newmol.AddAtom(x)
            for y in x.GetBonds():
                Natom = newmol.GetNumAtoms()
                newmol.AddAtom(y.GetOtherAtom(x))
                newmol.AddBond(0, Natom, y.GetBondType())
                for i, z in enumerate(m for m in y.GetOtherAtom(x).GetBonds() if m.GetIdx() != y.GetIdx()):
                    newmol.AddAtom(z.GetOtherAtom(y.GetOtherAtom(x)))
                    newmol.AddBond(Natom, Natom + i + 1, z.GetBondType())
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

    def Get_CBH4(self, smi):
        left, right, CBH4 = self.CBH3_Counter(smi), Counter({}), self.CBH4_Counter(smi)
        for k, v in CBH4.items(): right += Counter({x: y * v for x, y in self.CBH3_Counter(k).items()})
        reac = Counter([smi]) + right - left
        return reac, CBH4
