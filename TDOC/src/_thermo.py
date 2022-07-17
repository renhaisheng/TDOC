# -*- coding: utf-8 -*-
from collections import Counter, defaultdict
from os import chdir, listdir, system, remove
from itertools import chain
from shutil import copy 
from rdkit import Chem
from re import findall 
from ._CBH import CBH
from ._submit import SUBMIT
import numpy as np


class THERMO(SUBMIT):
	
    # To calculate thermodymamic data for fitting thermodynamic paramters.
    def Get_data_from_out(self, smi, micro_mol, method):
        HfX, HX = {'C': 0.2729652, 'H': 0.0830318, 'O': 0.0949038}, defaultdict(dict)
        HX['CCSDT']['M062X'] = {'C': -37.8380477, 'H': -0.4958344, 'O': -75.0586831}
        HX['CCSDT']['CCSDT'] = {'C': -37.7900384, 'H': -0.4976322, 'O': -75.0050381}
        HX['G4']['M062X'] = {'C': -37.8380477, 'H': -0.4958344, 'O': -75.0586831}
        HX['G4']['G4'] = {'C': -37.8318079, 'H': -0.4990600, 'O': -75.0431406}
        HX['CBSQB3']['M062X'] = {'C': -37.8380477, 'H': -0.4958344, 'O': -75.0586831}
        HX['CBSQB3']['CBSQB3'] = {'C': -37.7830150, 'H': -0.4974574, 'O': -74.9852586}
        if smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.path, method)): return None, None, None, None
        elif smi + '.out' not in listdir('{}/gaussian_out/GFNFF'.format(self.path)): return None, None, None, None
        elif smi + '.out' not in listdir('{}/gaussian_out/{}'.format(self.path, method)):
            if smi not in self.micro_mol and method == self.high_level and smi + '.out' in listdir('{}/gaussian_out/M062X'.format(self.path)): pass
            else:
                return None, None, None, None
        Hf_conf, S_conf, Cp_conf = self.Get_conf_correction_from_GFNFF(smi)
        newsmi = smi[:-2] if '-' in smi and smi[-1] in self.multiplicities else smi
        atoms = Counter(x.GetSymbol() for x in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAtoms())
        with open('{}/gaussian_out/M062X/{}.out'.format(self.path, smi)) as f: content = f.read()
        ZPE, Hcorr, Hmol = map(float, findall('- Thermochemistry -\n.+\n.+?298.150 Kelvin[\s\S]+?correction=\s+(\S+?)\s[\s\S]+?Enthalpy=\s+(\S+?)\n[\s\S]+?Enthalpies=\s+(\S+?)\n', content)[0])
        if ' - Thermochemistry For Hindered Internal Rotation -' not in content:
            E, Cv, S = zip(*findall('\n Total\s+(-*\d+\.\d+)\s+(\S+)\s+(\S+)\n', content))
            deltaH, Cp, S = 0, np.round(np.array(Cv[1:3] + Cv[0:1] + Cv[3:], float) + 1.987, 3), float(S[0])
        else:
            res = findall('\n Total\s+(-*\d+\.\d+)\s+(\S+)\s+(\S+)\n', content)
            E, Cv, S = zip(*res[4::5])
            deltaH, S, (_, corr_Cv, _) = float(res[1][1]) / 627.5095, float(res[2][1]) + float(S[0]), zip(*res[3::5])
            Cp = np.round(np.array(Cv[1: 3] + Cv[0:1] + Cv[3:], float) + np.array(corr_Cv[1:3] + corr_Cv[0:1] + corr_Cv[3:], float) + 1.987, 3)
        if method == self.high_level and smi in self.micro_mol:
            with open('{}/gaussian_out/{}/{}.out'.format(self.path, method, smi)) as f: content = f.read()
            if method == 'CCSDT':
                Hmol = self.Get_CBS_energy_from_out(smi) + Hcorr + deltaH - (1 - 0.97) * ZPE + self.Get_EBAC_from_smi(smi) / 627.5095
                Hfmol = 627.5095 * (Hmol - sum(v * HX[self.high_level][method][k] for k, v in atoms.items()) + sum(v * HfX[k] for k, v in atoms.items()))
            else:
                Hmol = float(findall(self.high_level + '=(-\d*\.\d*)', content.replace('\n ',  ''))[-1]) + Hcorr + deltaH - (1 - 0.97) * ZPE
                Hfmol = 627.5095 * (Hmol - sum(v * HX[self.high_level][method][k] for k, v in atoms.items()) + sum(v * HfX[k] for k, v in atoms.items()))
        else:
            Hmol = Hmol + deltaH - (1 - 0.97) * ZPE
            Hfmol = 627.5095 * (Hmol - sum(v * HX[self.high_level][method][k] for k, v in atoms.items()) + sum(v * HfX[k] for k, v in atoms.items()))
        if method == self.high_level and smi not in self.micro_mol:
            reac, prod = CBH().Get_CBH3(newsmi)
            try:
                Hfmol_M062X = 627.5095 * (Hmol - sum(v * HX[self.high_level]['M062X'][k] for k, v in atoms.items()) + sum(v * HfX[k] for k, v in atoms.items()))
                Hf_reac = sum(v * (self.Get_thermo_from_data(k, 298.15, method)[0] - self.Get_thermo_from_data(k, 298.15, 'M062X')[0] - self.Get_EBAC_from_smi(k)) for k, v in (reac - Counter([newsmi])).items())
                Hf_prod = sum(v * (self.Get_thermo_from_data(k, 298.15, method)[0] - self.Get_thermo_from_data(k, 298.15, 'M062X')[0] - self.Get_EBAC_from_smi(k)) for k, v in prod.items())
                Hfmol = Hfmol_M062X + Hf_prod - Hf_reac + self.Get_EBAC_from_smi(smi)
            except: smi, Hfmol, S, Cp = None, None, None, None
        if Hfmol:
            if Hf_conf != None: Hfmol, S, Cp = round(Hfmol, 3) + Hf_conf, S + S_conf, dict(tuple(zip([200, 250, 298.15] + list(range(300, 5050, 50)), np.round(Cp + Cp_conf, 3))))
            else: Hfmol, S, Cp = round(Hfmol, 3), S, dict(tuple(zip([200, 250, 298.15] + list(range(300, 5050, 50)), np.round(Cp, 3))))
        return smi, Hfmol, S, Cp

    # To derive CCSD(T)/CBS energy.
    def Get_CBS_energy_from_out(self, smi):
        with open('{}/gaussian_out/CCSDT/{}.out'.format(self.path, smi)) as f: content = f.read()
        HF = list(map(float, content.strip().split('\n')[-3].split()[::2][::-1]))
        CC = list(map(float, content.strip().split('\n')[-3].split()[1::2][::-1]))
        HF_CBS = (HF[1] ** 2 - HF[0] * HF[2]) / (2 * HF[1] - HF[0] - HF[2])
        corr_CBS = (343 * (CC[1] - HF[1]) - 125 * (CC[0] - HF[0])) / 218
        CC_CBS = round(HF_CBS + corr_CBS, 7)
        return CC_CBS

    # To obtain conformational corrections.
    def Get_conf_correction_from_GFNFF(self, smi):
        with open('{}/gaussian_out/GFNFF/{}.out'.format(self.path, smi), encoding='utf-8') as f: content, E, g = f.read(), [], [] 
        try:
            for x in [x.split() for x in findall('origin\n([\s\S]+?)\nT', content)[-1].split('\n')]:
                if len(x) > 5: E.append(x[2]), g.append(x[-2])
        except:
            return None, None, None
        E, g, Cp = np.array(E, dtype = float) * 627.5095 * 1000, np.array(g, dtype = int), {}
        E, T = E[0] * np.ones(len(E)) - E, 1 / (1.987 * np.array([200, 250, 298.15] + list(range(300, 5050, 50))))
        Hf = 1.987 * 298.15 / 1000 * np.dot(g,  - E * T[2]*np.exp(E * T[2])) / np.dot(g, np.exp(E * T[2]))
        S= - 1.987 * np.dot(g * np.exp(E * T[2]) / (np.dot(g, np.exp(E * T[2]))), np.log(g * np.exp(E * T[2]) / (np.dot(g, np.exp(E * T[2])))))
        E, T = E.reshape(1, len(E)), T.reshape(len(T), 1)
        Cp = 1.987 * np.dot((E * T) ** 2 * np.exp(E * T), g) / np.dot(np.exp(E * T), g) - 1.987 * (np.dot(E * T * np.exp(E * T), g) / np.dot(np.exp(E * T), g)) ** 2
        return round(Hf, 3), round(S, 3), np.round(Cp, 3)
    
    # To get BAC corrections.
    def Get_EBAC_from_smi(self, smi):
        newsmi = smi[:-2] if '-' in smi and smi[-1] in self.multiplicities else smi
        if Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAromaticAtoms(): BAC = {'CH1.0': 0.57,  'CC1.0': -1.6,  'CC1.5': -1.86,  'CC2.0': -2.47,  'CC3.0': -3.09,  'CO1.0': -2.51,  'CO2.0': -2.49,  'CO1.5': -2.09,  'HO1.0': 1.34}
        else: BAC = {'CH1.0': 0.14,  'CC1.0': -0.82,  'CC2.0': -2.21,  'CC3.0': -3.64,  'CO1.0': -1.46,  'CO2.0': -3.25,  'CO3.0': -3.9,  'OO1.0': -2.87,  'OO2.0': -4.15,  'HO1.0': -0.14,  'HH1.0': 0.67}
        bonds = Counter(''.join(sorted([x.GetBeginAtom().GetSymbol(), x.GetEndAtom().GetSymbol()], key = ['C', 'H', 'O'].index)) + str(x.GetBondTypeAsDouble()) for x in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetBonds())
        if newsmi == '[O][O]': bonds = Counter({'OO2.0': 1})
        EBAC = sum(BAC[k] * v for k, v in bonds.items())
        return EBAC

    # To read thermodynamic data at target temperature. 
    def Get_thermo_from_data(self, smi, T, method):
        func_Cp = lambda T, a1, a2, a3, a4, a5, :(a1 + a2 * T + a3 * T ** 2 + a4 * T ** 3 + a5 * T ** 4) * 1.987
        func_Hf = lambda T, a1, a2, a3, a4, a5, a6:(a1 * T + a2/2 * T ** 2 + a3 / 3 * T ** 3 + a4 / 4 * T ** 4 + a5/5 * T ** 5 + a6) * 1.987 / 1000
        func_S = lambda T, a1, a2, a3, a4, a5, a7:(a1 * np.log(T) + a2 * T + a3 / 2 * T ** 2 + a4 / 3 * T ** 3 + a5/4 * T ** 4 + a7) * 1.987
        with open('{}/chemkin_dat/{}/{}.dat'.format(self.path, method, smi)) as f: content = f.read()
        res = list(map(float, content[45:75].split()))
        Tmin, Tmax, Tbreak = res if len(res) == 3 else res + [1000]
        para = list(map(float, findall('.{15}', findall('\n([\s\S]+)\s{10}', content)[0])))
        if Tmax >= T >= Tbreak: Hf, S, Cp = func_Hf(T, *para[0:6]), func_S(T, *(para[0:5] + para[6:7])), func_Cp(T, *para[0:5])
        else: Hf, S, Cp = func_Hf(T, *para[7: 13]), func_S(T, *(para[7: 12] + para[13:])), func_Cp(T, *para[7: 12])
        return Hf, S, Cp
    
    # To output fitted themodyamic file.
    def Output_chemkin_format_data(self, smi, Hfmol, S, Cp, method):
        chdir('{}/src'.format(self.path))
        with open('input.dat', 'w') as f:
            f.write('H: kcal/mol,  S: cal/(mol k),  Cp: cal/(mol k)\n[name] {}\n[Hf]   {:.6f}\n[S]   {}\n[Cp]  {}\n'.format(smi, Hfmol, S, len(Cp)))
            for x in sorted(Cp): f.write('%12s\t\t%12s\n'%(x, Cp[x]))
        system('python thermofit.py')
        copy('out.dat', '{}/chemkin_dat/{}/{}.dat'.format(self.path, method, smi))
        print('Completed {}: {}'.format(method, smi)), remove('input.dat'), remove('out.dat')

    def Get_thermo_out(self):
        print('\nProcess thermodynamic data ...\n')
        for mols, method in zip([self.species, self.micro_mol, self.macro_mol], ['M062X', self.high_level, self.high_level]):
            for smi in mols:
                smi, Hfmol, S, Cp = self.Get_data_from_out(smi, self.micro_mol, method)
                if Hfmol: self.Output_chemkin_format_data(smi, Hfmol, S, Cp, method)

    # To write substituted thermodynamic file for mechanism.
    def Write_thermo_dat(self):
        smiles, smilesdict = self.Get_smiles()
        chdir('{}/program_out'.format(self.path))
        with open('thermo.dat', 'w') as f:
            f.write('THERM ALL\n   300.000  1000.000  5000.000\n')
            for species, smi in smilesdict.items():
                if smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.path, self.high_level)):
                    with open('{}/chemkin_dat/{}/{}.dat'.format(self.path, self.high_level, smi)) as p:
                        dat = findall('.{24}((.+\n){4})', p.read())[0][0]
                        f.write('{:<24}{}'.format(species, dat))
            f.write('END\n')
        with open('macro_mol.dat', 'w') as f:
            f.write('THERM ALL\n   300.000  1000.000  5000.000\n')
            for species, smi in smilesdict.items():
                if smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.path, self.high_level)):
                    if smi in self.macro_mol:
                        with open('{}/chemkin_dat/{}/{}.dat'.format(self.path, self.high_level, smi)) as p:
                            dat = findall('.{24}((.+\n){4})', p.read())[0][0]
                            f.write('{:<24}{}'.format(species, dat))
            f.write('END\n')
        with open('thermo_data.txt', 'w') as f:
            f.write('T: 298.15K,  Hf: kcal/mol,  S: cal/mol/K,  Cp: cal/mol/K\n')
            f.write('%16s\t%24s\t%20s%20s%20s%20s\n'%('species', 'smiles', 'Hf_M062X', 'Hf_' + self.high_level, 'S', 'Cp'))
            for species, smi in smilesdict.items():
                if all(smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.path, method)) for method in ['M062X', self.high_level]):
                    Hf_M062X, S, Cp = self.Get_thermo_from_data(smi, 298.15, 'M062X')
                    Hf_high_level, S, Cp = self.Get_thermo_from_data(smi, 298.15, self.high_level)
                    f.write('%16s\t%24s\t%20.2f%20.2f%20.2f%20.2f\n'%(species, smi, Hf_M062X, Hf_high_level, S, Cp))
        with open('smiles_data.txt', 'w') as f:
            f.write('T: 298.15K,  Hf: kcal/mol,  S: cal/mol/K,  Cp: cal/mol/K\n')
            f.write('%24s\t%20s%20s%20s%20s\n'%('smiles', 'Hf_M062X', 'Hf_' + self.high_level, 'S', 'Cp'))
            for smi in sorted(findall('(\S+?)\.dat', ' '.join(set(chain(*[listdir('{}/chemkin_dat/{}'.format(self.path, 'M062X'))])))), key = len):
                Hf = locals()
                for method in [self.high_level, 'M062X']:
                    if smi + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.path, method)): Hf[method], S, Cp = self.Get_thermo_from_data(smi, 298.15, method)
                    else: Hf[method], S, Cp = 0, 0, 0
                f.write('%24s\t%20.2f%20.2f%20.2f%20.2f\n'%(smi, Hf['M062X'], Hf[self.high_level], S, Cp))
        if all(x + '.dat' in listdir('{}/chemkin_dat/{}'.format(self.path, self.high_level)) for x in smiles):
            print('\n{:*^30}\n'.format(' Thermo completed! '))
            print('\n{:*^30}\n'.format(''))
            print('Exit this run with normal termination. \n')
        else:
            print('\n{:*^30}\n'.format(' Thermo incompleted! '))
            print('\n{:*^30}\n'.format(''))
            print('Exit this run and wait for the next run before Linux files are finished. \n')
            