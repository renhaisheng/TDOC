# -*- coding: utf-8 -*-
import numpy as np
from rdkit import Chem
from re import findall
from shutil import copy
from itertools import chain
from os import chdir, listdir
from collections import Counter

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
            if smi in self.para['macro_mol'] and method == self.para['high_level']:
                if smi + '.out' not in listdir('{}/gaussian_out/M062X'.format(self.para['work_path'])):
                    return None, None, None
            else:
                return None, None, None      
        
        # Get basic ZPE, H_corr, H_mol from low-level calculations.
        newsmi = smi[:-2] if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
        atoms = Counter(x.GetSymbol() for x in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAtoms())
        with open('{}/gaussian_out/M062X/{}.out'.format(self.para['work_path'], smi)) as f: content = f.read()
        ZPE, H_corr, H_mol = map(float, findall('- Thermochemistry -\n.+\n.+?298.150 Kelvin[\s\S]+?correction=\s+(\S+?)\s[\s\S]+?Enthalpy=\s+(\S+?)\n[\s\S]+?Enthalpies=\s+(\S+?)\n', content)[0])
        
        # Get E, S, Cv from Rigid-Rotor-Harmonic-Oscillator (RRHO) approximation.
        if ' - Thermochemistry For Hindered Internal Rotation -' not in content:
            E, Cv, S = zip(*findall('\n Total\s+(-*\d+\.\d+)\s+(\S+)\s+(\S+)\n', content))
            deltaH, Cp, S = 0, np.round(np.array(Cv[1:3] + Cv[0:1] + Cv[3:], float) + 1.987, 3), float(S[0])
        
        # Get E, S, Cv from Hindered-Rotor (HR) approximation.
        else:
            res = findall('\n Total\s+(-*\d+\.\d+)\s+(\S+)\s+(\S+)\n', content)
            E, Cv, S = zip(*res[4::5])
            deltaH, S, (_, corr_Cv, _) = float(res[1][1]) / 627.5095, float(res[2][1]) + float(S[0]), zip(*res[3::5])
            Cp = np.round(np.array(Cv[1: 3] + Cv[0:1] + Cv[3:], float) + np.array(corr_Cv[1:3] + corr_Cv[0:1] + corr_Cv[3:], float) + 1.987, 3)
        
        # Get Hf_mol from high-level calculations.
        if method == self.para['high_level'] and smi in self.para['micro_mol']:
            with open('{}/gaussian_out/{}/{}.out'.format(self.para['work_path'], method, smi)) as f: content = f.read()
            
            # Get Hf_mol from CCSD(T)/CBS with Bond Additivity Corrections (BACs).
            if method == 'CCSDT':
                H_mol = self.get_CBS_energy_from_out(smi) + H_corr + deltaH - (1 - 0.97) * ZPE + self.get_EBAC_from_smi(smi) / 627.5095
                Hf_mol = 627.5095 * (H_mol - sum(v * self.para['H_atoms']['CCSDT'][k] for k, v in atoms.items()) + sum(v * self.para['Hf_atoms'][k] for k, v in atoms.items()))
            else:
                H_mol = float(findall(self.para['high_level'] + '=(-\d*\.\d*)', content.replace('\n ',  ''))[-1]) + H_corr + deltaH - (1 - 0.97) * ZPE
                Hf_mol = 627.5095 * (H_mol - sum(v * self.para['H_atoms'][method][k] for k, v in atoms.items()) + sum(v * self.para['Hf_atoms'][k] for k, v in atoms.items()))
        else:
            H_mol = H_mol + deltaH - (1 - 0.97) * ZPE
            Hf_mol = 627.5095 * (H_mol - sum(v * self.para['H_atoms'][method][k] for k, v in atoms.items()) + sum(v * self.para['Hf_atoms'][k] for k, v in atoms.items()))
        
        # Get Hf from CBH-3 extrapolation.
        if method == self.para['high_level'] and smi not in self.para['micro_mol']:
            reac, prod = CBH(self.para['aro_rings']).get_CBH3(newsmi)
            thermodata = Thermofit(self.para)
            try:
                Hf_mol_M062X = 627.5095 * (H_mol - sum(v * self.para['H_atoms']['M062X'][k] for k, v in atoms.items()) + sum(v * self.para['Hf_atoms'][k] for k, v in atoms.items()))
                Hf_reac = sum(v * (thermodata.get_thermo_from_data(k, method)[0] - thermodata.get_thermo_from_data(k, 'M062X')[0]) for k, v in (reac - Counter([newsmi])).items())
                Hf_prod = sum(v * (thermodata.get_thermo_from_data(k, method)[0] - thermodata.get_thermo_from_data(k, 'M062X')[0]) for k, v in prod.items())
                Hf_mol = Hf_mol_M062X + Hf_prod - Hf_reac
            except: Hf_mol, S, Cp = None, None, None
        
        # Get corrections of Conformational Sampling (CS).
        Hf_conf, S_conf, Cp_conf = self.get_conf_correction_from_GFNFF(smi)
        if Hf_mol:
            if Hf_conf != None:
                Hf_mol, S, Cp = round(Hf_mol, 3) + Hf_conf, S + S_conf, dict(tuple(zip([200, 250, 298.15] + list(range(300, 5050, 50)), np.round(Cp + Cp_conf, 3))))
            else: 
                Hf_mol, S, Cp = round(Hf_mol, 3), S, dict(tuple(zip([200, 250, 298.15] + list(range(300, 5050, 50)), np.round(Cp, 3))))
        return Hf_mol, S, Cp



    """ To derive CCSD(T)/CBS energy. """
    def get_CBS_energy_from_out(self, smi):
        with open('{}/gaussian_out/CCSDT/{}.out'.format(self.para['work_path'], smi)) as f: content = f.read()
        HF = list(map(float, content.strip().split('\n')[-3].split()[::2][::-1]))
        CC = list(map(float, content.strip().split('\n')[-3].split()[1::2][::-1]))
        HF_CBS = (HF[1] ** 2 - HF[0] * HF[2]) / (2 * HF[1] - HF[0] - HF[2])
        corr_CBS = (343 * (CC[1] - HF[1]) - 125 * (CC[0] - HF[0])) / 218
        CC_CBS = round(HF_CBS + corr_CBS, 7)
        return CC_CBS



    """ To obtain conformational corrections. """
    def get_conf_correction_from_GFNFF(self, smi):
        with open('{}/gaussian_out/GFNFF/{}.out'.format(self.para['work_path'], smi), encoding='utf-8') as f:
            content, E, g = f.read(), [], [] 
        try:
            for x in [x.split() for x in findall('origin\n([\s\S]+?)\nT', content)[-1].split('\n')]:
                if len(x) > 5:
                    E.append(x[2]), g.append(x[-2])
        except:
            return None, None, None
        E, g, Cp = np.array(E, dtype = float) * 627.5095 * 1000, np.array(g, dtype = int), {}
        E, T = E[0] * np.ones(len(E)) - E, 1 / (1.987 * np.array([200, 250, 298.15] + list(range(300, 5050, 50))))
        Hf = 1.987 * 298.15 / 1000 * np.dot(g,  - E * T[2]*np.exp(E * T[2])) / np.dot(g, np.exp(E * T[2]))
        S = - 1.987 * np.dot(g * np.exp(E * T[2]) / (np.dot(g, np.exp(E * T[2]))), np.log(g * np.exp(E * T[2]) / (np.dot(g, np.exp(E * T[2])))))
        E, T = E.reshape(1, len(E)), T.reshape(len(T), 1)
        Cp = 1.987 * np.dot((E * T) ** 2 * np.exp(E * T), g) / np.dot(np.exp(E * T), g) - 1.987 * (np.dot(E * T * np.exp(E * T), g) / np.dot(np.exp(E * T), g)) ** 2
        return round(Hf, 3), round(S, 3), np.round(Cp, 3)


    
    """ To get BACs. """
    def get_EBAC_from_smi(self, smi):
        newsmi = smi[:-2] if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
        if Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAromaticAtoms():
            BAC = {'CH1.0': 0.57,  'CC1.0': -1.6,  'CC1.5': -1.86,  'CC2.0': -2.47,  'CC3.0': -3.09,  'CO1.0': -2.51,  'CO2.0': -2.49,  'CO1.5': -2.09,  'HO1.0': 1.34}
        else:
            BAC = {'CH1.0': 0.14,  'CC1.0': -0.82,  'CC2.0': -2.21,  'CC3.0': -3.64,  'CO1.0': -1.46,  'CO2.0': -3.25,  'CO3.0': -3.9,  'OO1.0': -2.87,  'OO2.0': -4.15,  'HO1.0': -0.14,  'HH1.0': 0.67}
        bonds = []
        for bond in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetBonds():
            bonds.append(''.join(sorted([bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()], key = ['C', 'H', 'O'].index)) + str(bond.GetBondTypeAsDouble()))
        all_bonds = Counter(bonds)
        if newsmi == '[O][O]':
            all_bonds = Counter({'OO2.0': 1})
        EBAC = sum(BAC[k] * v for k, v in all_bonds.items())
        return EBAC



    """ To write thermodynamic files. """
    def write_thermo_dat(self):
        
        # Go to working directory.
        chdir('{}/program_out'.format(self.para['work_path']))
        smiles, smilesdict = Parameters(self.para['input_file'], self.para['para_file'], self.para['work_path']).get_smiles()
        
        # Write thermodynamic parameters for all species to a dat file for Chemkin.
        with open('thermo.dat', 'w') as f:
            f.write('THERM ALL\n   300.000  1000.000  5000.000\n')
            for species, smi in smilesdict.items():
                if smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], self.para['high_level'])):
                    with open('{}/chemkin_dat/{}/{}.dat'.format(self.para['work_path'], self.para['high_level'], smi)) as p:
                        dat = findall('.{24}((.+\n){4})', p.read())[0][0]
                        f.write('{:<24}{}'.format(species, dat))
            f.write('END\n')
        
        # Write thermodyamic parameters for macro molecules to a dat file for Chemkin.
        with open('macro_mol.dat', 'w') as f:
            f.write('THERM ALL\n   300.000  1000.000  5000.000\n')
            for species, smi in smilesdict.items():
                if smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], self.para['high_level'])):
                    if smi in self.para['macro_mol']:
                        with open('{}/chemkin_dat/{}/{}.dat'.format(self.para['work_path'], self.para['high_level'], smi)) as p:
                            dat = findall('.{24}((.+\n){4})', p.read())[0][0]
                            f.write('{:<24}{}'.format(species, dat))
            f.write('END\n')
        
        # Write all thermodyamic data at 298.15 K.
        with open('thermo_data.txt', 'w') as f:
            f.write('T: 298.15K,  Hf: kcal/mol,  S: cal/mol/K,  Cp: cal/mol/K\n')
            f.write('%16s\t%24s\t%20s%20s%20s%20s\n'%('species', 'smiles', 'Hf_M062X', 'Hf_' + self.para['high_level'], 'S', 'Cp'))
            for species, smi in smilesdict.items():
                if all(smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], method)) for method in ['M062X', self.para['high_level']]):
                    Hf_M062X, S, Cp = Thermofit(self.para).get_thermo_from_data(smi, 'M062X')
                    Hf_high_level, S, Cp = Thermofit(self.para).get_thermo_from_data(smi, self.para['high_level'])
                    f.write('%16s\t%24s\t%20.2f%20.2f%20.2f%20.2f\n'%(species, smi, Hf_M062X, Hf_high_level, S, Cp))
        
        # Write all thermodyamic data for all smiles at 298.15 K.
        with open('smiles_data.txt', 'w') as f:
            f.write('T: 298.15K,  Hf: kcal/mol,  S: cal/mol/K,  Cp: cal/mol/K\n')
            f.write('%24s\t%20s%20s%20s%20s\n'%('smiles', 'Hf_M062X', 'Hf_' + self.para['high_level'], 'S', 'Cp'))
            for smi in sorted(findall('(\S+?)\.dat', ' '.join(set(chain(*[listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], 'M062X'))])))), key = len):
                Hf = locals()
                for method in [self.para['high_level'], 'M062X']:
                    if smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], method)): Hf[method], S, Cp = Thermofit(self.para).get_thermo_from_data(smi, method)
                    else: Hf[method], S, Cp = 0, 0, 0
                f.write('%24s\t%20.2f%20.2f%20.2f%20.2f\n'%(smi, Hf['M062X'], Hf[self.para['high_level']], S, Cp))
        
        # Check whether all thermodynamic parameters are completed.
        if all(x + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.para['work_path'], self.para['high_level'])) for x in smiles):
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
        for mols, method in zip([self.para['species'], self.para['micro_mol'], self.para['macro_mol']], ['M062X', self.para['high_level'], self.para['high_level']]):
            for smi in mols:
                Hf, S, Cp = self.get_data_from_out(smi, method)
                if Hf:
                    Thermofit(self.para).output_dat_from_data(smi, Hf, S, Cp, method)
                    print('Completed {}: {}'.format(method, smi))
        self.write_thermo_dat()
         