# -*- coding: utf-8 -*-
import numpy as np
from math import ceil
from rdkit import Chem
from re import findall
from itertools import chain
from os import chdir, listdir
from collections import Counter
from rdkit.Chem import AllChem, rdMolTransforms


from tdoc._CBH import CBH
from tdoc._thermofit import Thermofit
from tdoc._parameters import Parameters



class Thermo(object):

    """ To initiate parameters for Thermo. """
    def __init__(self, para):
        self.para = para



    """ To calculate thermodymamic data for fitting thermodynamic paramters. """
    def get_data_from_out(self, smi, method):

        # Bypass the molecules that do not need thermal fitting. 
        if smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], method)):
            return None, None, None
        elif smi + '.out' not in listdir('{}/gaussian_out/GFNFF'.format(self.para['work_path'])):
            return None, None, None
        elif smi + '.out' not in listdir('{}/gaussian_out/{}'.format(self.para['work_path'], method)):
            if smi in self.para['macro_mol'] and method == 'CCSDT':
                if smi + '.out' not in listdir('{}/gaussian_out/MP2'.format(self.para['work_path'])):
                    return None, None, None
            else:
                return None, None, None
        
        # Get basic ZPE, H_corr, H_mol from low-level calculations.
        newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
        newsmi = newsmi[4:] if smi.startswith('cis-') else newsmi
        atoms = Counter(x.GetSymbol() for x in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAtoms())
        with open('{}/gaussian_out/B3LYP/{}.out'.format(self.para['work_path'], smi)) as f: content = f.read()
        ZPE, H_corr, H_mol = map(float, findall('- Thermochemistry -\n.+\n.+?298.150 Kelvin[\s\S]+?correction=\s+(\S+?)\s[\s\S]+?Enthalpy=\s+(\S+?)\n[\s\S]+?Enthalpies=\s+(\S+?)\n', content)[0])
        
        # Get E, S, Cv from Rigid-Rotor-Harmonic-Oscillator (RRHO) approximation.
        if ' - Thermochemistry For Hindered Internal Rotation -' not in content:
            E, Cv, S = zip(*findall('\n Total\s+(-*\d+\.\d+)\s+(\S+)\s+(\S+)\n', content))
            deltaH, Cp, S = 0, np.round(np.array(Cv, float) + 1.987, 3), float(S[0])
        
        # Get E, S, Cv from Hindered-Rotor (HR) approximation.
        else:
            res = findall('\n Total\s+(-*\d+\.\d+)\s+(\S+)\s+(\S+)\n', content)
            E, Cv, S = zip(*res[4::5])
            deltaH, S, (_, corr_Cv, _) = float(res[1][1]) / 627.5095, float(res[2][1]) + float(S[0]), zip(*res[3::5])
            Cp = np.round(np.array(Cv, float) + np.array(corr_Cv, float) + 1.987, 3)
        
        temp=np.array(findall('- Thermochemistry -\n.+\n Temperature\s+(\S+) Kelvin', content), float)

       
        # Check BAC and BDC values. 
        BAC = self.get_BAC_from_smi(smi)
        BDC = self.get_BDC_from_smi(smi)
        if type(BAC) == dict:
            return BAC
        if type(BDC) == dict:
            return BDC

        # Get Hf_mol from CCSDT calculations.
        if method == 'CCSDT' and smi in self.para['micro_mol']:
            with open('{}/gaussian_out/{}/{}.out'.format(self.para['work_path'], method, smi)) as f: content = f.read()
            
            # Get Hf_mol from CCSD(T)/CBS with Bond Additivity Corrections (BACs).
            H_mol = self.get_CBS_energy_from_out(smi) + H_corr + deltaH - (1 - self.para['zpe_scale']) * ZPE
            Hf_mol = 627.5095 * (H_mol - sum(v * self.para['H_atoms']['CCSDT'][k] for k, v in atoms.items()) + sum(v * self.para['Hf_atoms'][k] for k, v in atoms.items()))
            Hf_mol = Hf_mol + self.get_soc_from_smi(smi) + BAC

        # Get Hf_mol from MP2 calculations.
        elif method == 'MP2':
            H_mol = self.get_MP2_energy_from_out(smi) + H_corr + deltaH - (1 - self.para['zpe_scale']) * ZPE
            Hf_mol = 627.5095 * (H_mol - sum(v * self.para['H_atoms'][method][k] for k, v in atoms.items()) + sum(v * self.para['Hf_atoms'][k] for k, v in atoms.items()))
            Hf_mol = Hf_mol + self.get_soc_from_smi(smi) + BAC
        
        # Get Hf from CBH-3 extrapolation.
        if method == 'CCSDT' and smi not in self.para['micro_mol']:
            reac, prod = CBH(self.para).get_CBH3(newsmi)
            thermodata = Thermofit(self.para)
            try:
                H_mol = self.get_MP2_energy_from_out(smi) + H_corr + deltaH - (1 - self.para['zpe_scale']) * ZPE
                Hf_mol_MP2 = 627.5095 * (H_mol - sum(v * self.para['H_atoms']['MP2'][k] for k, v in atoms.items()) + sum(v * self.para['Hf_atoms'][k] for k, v in atoms.items()))
                Hf_reac = sum(v * (thermodata.get_thermo_from_data(k, method)[0] - thermodata.get_thermo_from_data(k, 'MP2')[0]) for k, v in (reac - Counter([newsmi])).items())
                Hf_prod = sum(v * (thermodata.get_thermo_from_data(k, method)[0] - thermodata.get_thermo_from_data(k, 'MP2')[0]) for k, v in prod.items())
                Hf_mol = Hf_mol_MP2 + Hf_prod - Hf_reac + self.get_soc_from_smi(smi) + BAC + BDC
            except:
                Hf_mol, S, Cp = None, None, None
        
        # Get corrections of Conformational Sampling (CS).
        Hf_conf, S_conf, Cp_conf = self.get_conf_correction_from_GFNFF(smi, temp)

        if Hf_mol:
            if Hf_conf != None:
                Hf_mol, S, Cp = round(Hf_mol, 3) + Hf_conf, S + S_conf, dict(tuple(zip(temp, np.round(Cp + Cp_conf, 3))))
            else: 
                Hf_mol, S, Cp = round(Hf_mol, 3), S, dict(tuple(zip(temp, np.round(Cp, 3))))

        return Hf_mol, S, Cp



    """ To MP2 energy. """
    def get_MP2_energy_from_out(self, smi):
        with open('{}/gaussian_out/MP2/{}.out'.format(self.para['work_path'], smi)) as f: content = f.read()
        MP2 = float(content.strip().split('\n')[-3].split()[0])
        MP2 = round(MP2, 7)
        return MP2



    """ To derived CCSD(T)/CBS energy. """
    def get_CBS_energy_from_out(self, smi):
        with open('{}/gaussian_out/CCSDT/{}.out'.format(self.para['work_path'], smi)) as f: content = f.read()
        HF = list(map(float, content.strip().split('\n')[-3].split()[::2][::-1]))
        CC = list(map(float, content.strip().split('\n')[-3].split()[1::2][::-1]))
        HF_CBS = (HF[1] ** 2 - HF[0] * HF[2]) / (2 * HF[1] - HF[0] - HF[2])
        corr_CBS = (343 * (CC[1] - HF[1]) - 125 * (CC[0] - HF[0])) / 218
        CC_CBS = round(HF_CBS + corr_CBS, 7)
        return CC_CBS



    """ To obtain conformational corrections. """
    def get_conf_correction_from_GFNFF(self, smi, temp):
        with open('{}/gaussian_out/GFNFF/{}.out'.format(self.para['work_path'], smi), encoding='utf-8') as f:
            content, E, g = f.read(), [], [] 
        try:
            for x in [x.split() for x in findall('origin\n([\s\S]+?)\nT', content)[-1].split('\n')]:
                if len(x) > 5:
                    E.append(x[2]), g.append(x[-2])
        except:
            return None, None, None
        E, g, Cp = np.array(E, float) * 627.5095 * 1000, np.array(g, int), {}
        E, T = E[0] * np.ones(len(E)) - E, 1 / (1.987 * temp)
        Hf = 1.987 * 298.15 / 1000 * np.dot(g,  - E * T[2]*np.exp(E * T[2])) / np.dot(g, np.exp(E * T[2]))
        S = - 1.987 * np.dot(g * np.exp(E * T[2]) / (np.dot(g, np.exp(E * T[2]))), np.log(g * np.exp(E * T[2]) / (np.dot(g, np.exp(E * T[2])))))
        E, T = E.reshape(1, len(E)), T.reshape(len(T), 1)
        Cp = 1.987 * np.dot((E * T) ** 2 * np.exp(E * T), g) / np.dot(np.exp(E * T), g) - 1.987 * (np.dot(E * T * np.exp(E * T), g) / np.dot(np.exp(E * T), g)) ** 2
        return round(Hf, 3), round(S, 3), np.round(Cp, 3)


    def get_soc_from_smi(self, smi):
        soc = 0
        newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
        newsmi = newsmi[4:] if smi.startswith('cis-') else newsmi
        atoms = Counter(x.GetSymbol() for x in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAtoms())
        if smi not in ['[C]-T', '[O]-T']:
            soc = atoms.get('C', 0) * 627.5095 / 1000 * (-0.14) + atoms.get('O', 0) * 627.5095 / 1000 * (-0.36)
        return soc


    """ To get atom information. """
    def get_BAC_atom_info(self, mol, atom):
        p_hybrid_degree = str(atom.GetHybridization())[-1] if str(atom.GetHybridization())[-1] in ["2", "3"] else "1"
        electrons = atom.GetNumRadicalElectrons()
        
        # To get atom info by EP method.
        atom_info = '{}{}{}'.format(atom.GetSymbol(), electrons, p_hybrid_degree) if atom.GetSymbol() != 'H' else 'H'
        return atom_info



    """ To get BAC bonds. """
    def get_BAC_bonds_from_smi(self, smi):
        all_bonds = []
        newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
        newsmi = newsmi[4:] if smi.startswith('cis-') else newsmi
        mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))
        
        # To get all bonds.
        for bond in mol.GetBonds():
            left_atom = bond.GetBeginAtom()
            right_atom = bond.GetEndAtom()
            bond_info = self.para['bond_marks'][bond.GetBondTypeAsDouble()].join(sorted([self.get_BAC_atom_info(mol,left_atom), self.get_BAC_atom_info(mol,right_atom)]))
            if bond_info[0] == 'H':
                bond_info = bond_info[2:] + bond_info[1] + bond_info[0]
            all_bonds.append(bond_info)
        all_bonds = Counter(all_bonds)
        return all_bonds



    """ To get BACs. """
    def get_BAC_from_smi(self, smi):
        BAC_bonds = self.get_BAC_bonds_from_smi(smi)
        try:
            BAC = sum([self.para['BAC_parameters'][k] * v for k, v in BAC_bonds.items()])
        except:
            return {'BAC': {smi: ', '.join(set(BAC_bonds) - set(self.para['BAC_parameters']))}}
        return BAC



    """ To BDC atom information. """
    def get_BDC_atom_info(self, mol, atom):
        p_hybrid_degree = str(atom.GetHybridization())[-1] if str(atom.GetHybridization())[-1] in ["2", "3"] else "1"
        adjacent_heavy_atoms = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() != 'H']
        adjacent_H_atoms = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() == 'H']
        
        # To get atom info by HHP method.
        atom_info = '{}{}{}{}'.format(atom.GetSymbol(), len(adjacent_heavy_atoms), len(adjacent_H_atoms), p_hybrid_degree) if atom.GetSymbol() != 'H' else 'H'
        return atom_info


    
    """ To get singlet bonds. """
    def get_singlet_bonds(self, smi):
        all_bonds = []
        newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
        newsmi = newsmi[4:] if smi.startswith('cis-') else newsmi
        mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))

        # To get all bonds.
        for bond in mol.GetBonds():
            left_atom = bond.GetBeginAtom()
            right_atom = bond.GetEndAtom()
            bond_info = self.para['bond_marks'][bond.GetBondTypeAsDouble()].join(sorted([self.get_BDC_atom_info(mol, left_atom), self.get_BDC_atom_info(mol, right_atom)]))
            if bond_info[0] == 'H':
                bond_info = bond_info[2:] + bond_info[1] + bond_info[0]
            all_bonds.append(bond_info)
        all_bonds = Counter(all_bonds)
        return all_bonds



    """ To get doublet bonds. """
    def get_doublet_bonds(self, smi):
        all_bonds = []
        newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
        newsmi = newsmi[4:] if smi.startswith('cis-') else newsmi
        mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))
        
        # To get all bonds.
        for atom in mol.GetAtoms():
            if len([True for x in atom.GetNeighbors() if x.GetSymbol() != 'H']) <2 : continue
            info = sorted([(self.get_BDC_atom_info(mol, x), self.para['bond_marks'][mol.GetBondBetweenAtoms(x.GetIdx(), atom.GetIdx()).GetBondTypeAsDouble()]) for x in atom.GetNeighbors() if x.GetSymbol() != 'H'])
            all_atoms = [x for x in atom.GetNeighbors() if x.GetSymbol() != 'H'] + [atom]
            ring_tag = 'ring-' if all([x.IsInRing() for x in all_atoms]) else ''
            if len(info[:ceil(len(info) / 2)]) > 1:
                left_outside_info = '[{}]'.format('|'.join([x[0] + x[1] for x in info[:ceil(len(info) / 2)]]))
            else:
                left_outside_info = info[0][0] + info[0][1]
            if len(info[ceil(len(info) / 2):]) > 1:
                right_outside_info = '[{}]'.format('|'.join([x[1] + x[0] for x in info[ceil(len(info) / 2):]]))
            else:
                right_outside_info = info[-1][-1] + info[-1][0]
            bond_info = '{}{}{}{}'.format(ring_tag, left_outside_info, self.get_BDC_atom_info(mol, atom), right_outside_info)
            all_bonds.append(bond_info)
        all_bonds = Counter(all_bonds)
        return all_bonds



    """ To get cis-trans bonds. """
    def get_cis_trans_bonds(self, smi):
        all_bonds = []
        newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
        newsmi = newsmi[4:] if smi.startswith('cis-') else newsmi
        mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))
        AllChem.Compute2DCoords(mol)

        # To get all bonds.
        for bond in mol.GetBonds():
            left_atom = bond.GetBeginAtom()
            right_atom = bond.GetEndAtom()
            left_heavy_atoms = [x.GetIdx() for x in left_atom.GetNeighbors() if x.GetSymbol() != 'H' and x.GetIdx() != right_atom.GetIdx()]
            right_heavy_atoms = [x.GetIdx() for x in right_atom.GetNeighbors() if x.GetSymbol() != 'H' and x.GetIdx() != left_atom.GetIdx()]
            left_bonds = [mol.GetBondBetweenAtoms(left_atom.GetIdx(), atom_idx) for atom_idx in left_heavy_atoms]
            right_bonds = [mol.GetBondBetweenAtoms(right_atom.GetIdx(), atom_idx) for atom_idx in right_heavy_atoms]
            left_central_atom_info = self.get_BDC_atom_info(mol, left_atom)
            right_central_atom_info = self.get_BDC_atom_info(mol, right_atom)
            left_outside = sorted([(self.get_BDC_atom_info(mol, x.GetOtherAtom(left_atom)), self.para['bond_marks'][x.GetBondTypeAsDouble()], x.GetOtherAtom(left_atom).GetIdx()) for x in left_bonds])
            right_outside = sorted([(self.get_BDC_atom_info(mol, x.GetOtherAtom(right_atom)), self.para['bond_marks'][x.GetBondTypeAsDouble()], x.GetOtherAtom(right_atom).GetIdx()) for x in right_bonds])
            if left_heavy_atoms and right_heavy_atoms:
                info = sorted([[left_outside, left_central_atom_info], [right_outside, right_central_atom_info]])
                left_outside_info = '|'.join([x[0] + x[1] for x in info[0][0]])
                right_outside_info = '|'.join([x[1] + x[0] for x in info[1][0]])
                left_outside_info = '[{}]'.format(left_outside_info) if '|' in left_outside_info else left_outside_info
                right_outside_info = '[{}]'.format(right_outside_info) if '|' in right_outside_info else right_outside_info
                if ([x.GetIsConjugated() for x in left_bonds + right_bonds].count(True) >= 2 and bond.GetIsConjugated()) or bond.GetBondTypeAsDouble() > 1.0:
                    dihedral_angle = abs(rdMolTransforms.GetDihedralDeg(mol.GetConformer(), left_outside[0][2], left_atom.GetIdx(), right_atom.GetIdx(), right_outside[0][2]))
                    cis_tag = 'cis-' if dihedral_angle < 90 else ''
                    cis_tag = 'cis-' if smi.startswith('cis-') else cis_tag
                    bond_info = '{}{}{}{}{}{}'.format(cis_tag, left_outside_info, info[0][1], self.para['bond_marks'][bond.GetBondTypeAsDouble()], info[1][1], right_outside_info)
                    all_bonds.append(bond_info) 
        all_bonds = Counter(all_bonds)
        return all_bonds



    """ To get mol bonds from smi. """
    def get_mol_bonds_from_smi(self, smi):
        all_bonds = self.get_singlet_bonds(smi) + self.get_doublet_bonds(smi) + self.get_cis_trans_bonds(smi)
        return all_bonds


    """ To get delta bonds from smi. """
    def get_delta_bonds_from_smi(self, smi):
        reac, prod = CBH(self.para).get_CBH3(smi)
        if reac == prod:
            return Counter({})
        bond_types = []
        for k, v in (reac + prod).items():
            bonds = self.get_mol_bonds_from_smi(k)
            bond_types.extend(bonds.keys())
        bond_types = sorted(Counter(bond_types))
        L, R= [0] * len(bond_types), [0] * len(bond_types)
        for x, y in reac.items():
            bonds = self.get_mol_bonds_from_smi(x)
            for i, v in enumerate(bond_types):
                L[i] = L[i] + bonds.get(v, 0) * y
        for x, y in prod.items():
            bonds = self.get_mol_bonds_from_smi(x)
            for i, v in enumerate(bond_types):
                R[i] = R[i] + bonds.get(v, 0) * y
        delta_coefficients = np.array(R) - np.array(L)
        all_bonds = Counter(dict(zip(bond_types, delta_coefficients)))

        # To remove vacant bonds.
        for x in list(all_bonds):
            if all_bonds[x] == 0:
                del all_bonds[x]
        return all_bonds



    """ To get BDCs. """
    def get_BDC_from_smi(self, smi):
        BDC_bonds = self.get_delta_bonds_from_smi(smi)
        try:
            BDC = sum([self.para['BDC_parameters'][k] * v for k, v in BDC_bonds.items()])
        except:
            return {'BDC': {smi: ', '.join(set(BDC_bonds)-set(self.para['BDC_parameters']))}}
        return BDC



    """ To write thermodynamic files. """
    def write_thermo_dat(self):
        
        # Go to working directory.
        chdir('{}/program_out'.format(self.para['work_path']))
        smiles, smilesdict = Parameters(self.para['input_file'], self.para['para_file'], self.para['work_path']).get_smiles()
        
        # Write thermodynamic parameters for all species to a dat file for Chemkin.
        with open('thermo.dat', 'w') as f:
            f.write('THERM ALL\n   300.000  1000.000  5000.000\n')
            for species, smi in smilesdict.items():
                if smi + '.dat' in listdir('{}/chemkin_dat/CCSDT'.format(self.para['work_path'])):
                    with open('{}/chemkin_dat/CCSDT/{}.dat'.format(self.para['work_path'], smi)) as p:
                        dat = findall('.{24}((.+\n){4})', p.read())[0][0]
                        f.write('{:<24}{}'.format(species, dat))
            f.write('END\n')
        
        # Write thermodyamic parameters for macro molecules to a dat file for Chemkin.
        with open('macro_mol.dat', 'w') as f:
            f.write('THERM ALL\n   300.000  1000.000  5000.000\n')
            for species, smi in smilesdict.items():
                if smi + '.dat' in listdir('{}/chemkin_dat/CCSDT'.format(self.para['work_path'])):
                    if smi in self.para['macro_mol']:
                        with open('{}/chemkin_dat/CCSDT/{}.dat'.format(self.para['work_path'], smi)) as p:
                            dat = findall('.{24}((.+\n){4})', p.read())[0][0]
                            f.write('{:<24}{}'.format(species, dat))
            f.write('END\n')
        
        # Write all thermodyamic data at 298.15 K.
        with open('thermo_data.txt', 'w') as f:
            f.write('T: 298.15K,  Hf: kcal/mol,  S: cal/mol/K,  Cp: cal/mol/K\n')
            f.write('%30s\t%60s\t%20s%20s%20s%20s\n'%('species', 'smiles', 'Hf_MP2', 'Hf_CCSDT', 'S', 'Cp'))
            for species, smi in smilesdict.items():
                if all(smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], method)) for method in ['MP2', 'CCSDT']):
                    Hf_MP2, S, Cp = Thermofit(self.para).get_thermo_from_data(smi, 'MP2')
                    Hf_CCSDT, S, Cp = Thermofit(self.para).get_thermo_from_data(smi, 'CCSDT')
                    f.write('%30s\t%60s\t%20.2f%20.2f%20.2f%20.2f\n'%(species, smi, Hf_MP2, Hf_CCSDT, S, Cp))
        
        # Write all thermodyamic data for all smiles at 298.15 K.
        with open('smiles_data.txt', 'w') as f:
            f.write('T: 298.15K,  Hf: kcal/mol,  S: cal/mol/K,  Cp: cal/mol/K\n')
            f.write('%60s\t%20s%20s%20s%20s\n'%('smiles', 'Hf_MP2', 'Hf_CCSDT', 'S', 'Cp'))
            for smi in sorted(findall('(\S+?)\.dat', ' '.join(set(chain(*[listdir('{}/chemkin_dat/MP2'.format(self.para['work_path']))])))), key = len):
                Hf = locals()
                for method in ['CCSDT', 'MP2']:
                    if smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], method)):
                        Hf[method], S, Cp = Thermofit(self.para).get_thermo_from_data(smi, method)
                    else:
                        break
                else:
                    f.write('%60s\t%20.2f%20.2f%20.2f%20.2f\n'%(smi, Hf['MP2'], Hf['CCSDT'], S, Cp))
                
        # Check whether all thermodynamic parameters are completed.
        if all(x + '.dat' in listdir('{}/chemkin_dat/CCSDT'.format(self.para['work_path'])) for x in smiles):
            print('\n{:*^30}\n'.format(' Thermo completed! '))
            print('\n{:*^30}\n'.format(''))
            print('Exit this run with normal termination. \n')
        else:
            print('\n{:*^30}\n'.format(' Thermo incompleted! '))
            print('\n{:*^30}\n'.format(''))
            print('Exit this run and wait for the next run before Linux files are finished. \n')



    """ To get all thermodynamic parameters. """
    def get_all_thermodat(self):
        print('\nProcess thermodynamic data ...\n')
        missing_parameters = {'BAC': {}, 'BDC': {}}
        for mols, method in zip([self.para['species'], self.para['micro_mol'], self.para['macro_mol']], ['MP2', 'CCSDT', 'CCSDT']):
            for smi in mols:
                res = self.get_data_from_out(smi, method)
                if len(res) == 1:
                    for k, v in res.items():
                        missing_parameters[k].update(v)
                elif res[0] != None:
                    Hf, S, Cp = res
                    Thermofit(self.para).output_dat_from_data(smi, Hf, S, Cp, method)
                    print('Completed {}: {}'.format(method, smi))
        
        # Output all missing BAC and BDC parameters.
        for smi in missing_parameters['BAC']:
            print('{:30} missing BAC parameters: {}'.format(smi, missing_parameters['BAC'][smi]))
        if missing_parameters['BAC'] and missing_parameters['BDC']:
            print('\n{:*^30}\n'.format(''))
        for smi in missing_parameters['BDC']:
            print('{:30} missing BDC parameters: {}'.format(smi, missing_parameters['BDC'][smi]))
        
        # Write all thermodynamic data.
        self.write_thermo_dat()
         