# -*- coding: utf-8 -*-
from ase.formula import Formula
from collections import Counter
from scipy import optimize
from rdkit import Chem
import numpy as np
import matplotlib.pyplot as plt


#Basic parameters.
multiplicities = {'S':'1', 'D':'2', 'T':'3', 'Q':'4', 'P':'5'}
with open('input.dat') as f: content = f.readlines()
species = content[1].split()[1]
Hf298 = float(content[2].split()[1]) * 1000
S298 = float(content[3].split()[1])
Tlist, Cp = list(zip(*[map(float, x.split()) for x in content[5:]]))

#NASA polynomials.
Cp_fit = lambda T, p0, p1, p2, p3, p4, p7, p8, p9, p10, p11: np.piecewise(T, [T >= 1000, T < 1000], [lambda T:(p0 + p1 * T + p2 * T ** 2 + p3 * T ** 3 + p4 * T ** 4) * 1.987, lambda T:(p7 + p8 * T + p9 * T ** 2 + p10 * T ** 3 + p11 * T ** 4) *1.987])
funcCp = lambda T, a0, a1, a2, a3, a4, :(a0 + a1 * T + a2 * T ** 2 + a3 * T ** 3 + a4 * T ** 4) * 1.987
func_Hf = lambda T, a0, a1, a2, a3, a4, a5:(a0 * T + a1 / 2 * T ** 2 + a2 / 3 * T ** 3 + a3 / 4 * T ** 4 + a4 / 5 * T ** 5 + a5) * 1.987
func_S = lambda T, a0, a1, a2, a3, a4, a6:(a0 * np.log(T) + a1 * T + a2 / 2 * T ** 2 + a3 / 3 * T ** 3 + a4 / 4 * T ** 4 + a6) * 1.987
p0, p1, p2, p3, p4, p7, p8, p9, p10, p11 = optimize.curve_fit(Cp_fit, Tlist, Cp)[0]
para = [p0, p1, p2, p3, p4, 0, 0, p7, p8, p9, p10, p11, 0, 0]
para[12] = Hf298 / 1.987 - func_Hf(298.15, *para[7:13]) / 1.987
para[5] = func_Hf(1000, *para[7:13]) / 1.987 - func_Hf(1000, *para[:6]) / 1.987
para[13] = S298 / 1.987 - func_S(298.15, *para[7:12], para[13]) / 1.987
para[6] = func_S(1000, *para[7:12], para[13]) / 1.987 - func_S(1000, *para[:5], para[6]) / 1.987

#Thermodynamic parameters.
try:
    smi = species[:-2] if '-' in species and species[-1] in multiplicities else species
    atoms = Counter(x.GetSymbol() for x in Chem.AddHs(Chem.MolFromSmiles(smi)).GetAtoms())
except:
    atoms = Formula(species).count()
atomlist = ''.join('{:<2}{:>3}'.format(x, atoms[x]) for x in atoms)
species = ''.join(x + str(atoms[x]) if atoms[x] != 1 else x for x in atoms)
with open('out.dat', 'w') as f:
    f.write('{:<24}{:<20}G{:>10.3f}{:>10.3f}{:>8.2f}{:>7}\n'.format(species, atomlist, Tlist[0], Tlist[-1], 1000.00, 1))
    f.write('{:>15.8E}{:>15.8E}{:>15.8E}{:>15.8E}{:>15.8E}{:>5}\n'.format(*para[:5], 2))
    f.write('{:>15.8E}{:>15.8E}{:>15.8E}{:>15.8E}{:>15.8E}{:>5}\n'.format(*para[5:10], 3))
    f.write('{:>15.8E}{:>15.8E}{:>15.8E}{:>15.8E}{:>20}\n'.format(*para[10:], 4))
