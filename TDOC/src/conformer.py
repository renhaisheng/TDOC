# -*- coding: utf-8 -*-
try:
    from ChemScript14 import *
    NoChemScript = False
except:
    NoChemScript = True
from os import system, remove
from re import findall
from glob import glob


multiplicities = {'S':'1', 'D':'2', 'T':'3', 'Q':'4', 'P':'5'}
with open('input.txt') as f: newsmi, smi, number_process, conformer_method, spin_method, key_words = [x.strip() for x in f.readlines()]
if 'press' not in key_words: key_words = key_words + ' press=0.987'
if conformer_method == '1' and NoChemScript: conformer_method = '2'
link_words = '# freq=(readfc,hinderedrotor)' if 'hinderedrotor' in key_words else '# freq=readfc'
with open('rdkit.gjf') as f: benchmark = f.readlines()
if conformer_method == '1':
    mol = StructureData().LoadData(newsmi, 'smiles')
    mol.ConvertTo3DStructure()
    mol.WriteFile('chem3d.mol')
    system('obabel chem3d.mol -O chem3d.gjf')
    with open('chem3d.gjf') as f: content = f.readlines()
    remove('chem3d.mol'), remove('chem3d.gjf')
    if len(benchmark) != len(content): exit()
if conformer_method == '2':
    system('obabel -:"{}" -O obabel.mol --gen3D'.format(newsmi))
    system('obabel obabel.mol -O obabel.gjf')
    with open('obabel.gjf') as f: content = f.readlines()
    if len(benchmark) != len(content): conformer_method = '3'
    remove('obabel.mol'), remove('obabel.gjf')
if conformer_method == '3':
    with open('rdkit.gjf') as f: content = f.readlines()
if newsmi != smi:
    spin = multiplicities[smi[-1]]
elif spin_method == '1':
    with open('rdkit.gjf') as f: spin = f.readlines()[5].strip()[-1]
else:
    system('obabel -:"' + newsmi + '" -O obabel.mol --gen3D')
    system('obabel obabel.mol -O obabel.gjf')
    with open('obabel.gjf') as f: spin = f.readlines()[5].strip()[-1]
    remove('obabel.mol'), remove('obabel.gjf')
if spin != '1': key_words = key_words.replace('M062X', 'UM062X')
with open(smi + '.gjf', 'w') as f:
    f.write('%nprocshared={}\n%chk={}.chk\n{}\n\n{}\n\n0 {}\n{}'.format(number_process, smi, key_words, smi + '_' + conformer_method, spin, ''.join(content[6:])))
    for T in range(200, 5050, 50):
        f.write('\n--Link1--\n%chk={}.chk\n{} temp={} press=0.987 geom=allcheck guess=read\n'.format(smi, link_words, T))
