# -*- coding: utf-8 -*-
from os import system, chdir, listdir, remove, makedirs
from utilspie.iterutils import get_chunks
from re import search, findall, sub
from collections import Counter
from rdkit.Chem import AllChem
from shutil import copy, which
from itertools import chain
from ase import Atom, Atoms
from random import sample
from rdkit import Chem
from math import ceil
from glob import glob
from openbabel import *
from ._project import PROJECT


class SUBMIT(PROJECT):

    # To get submitted file of M062X.
    def Get_M062X_out(self):
        chdir('{}/gaussian_out/M062X'.format(self.path))
        smiles = [smi for smi in self.species if smi + '.dat' not in  listdir('{}/chemkin_dat/{}'.format(self.path, 'M062X')) and smi + '.xyz' not in listdir('../GFNFF')]
        self.Check_method_out(smiles, 'M062X')
        if not smiles:
            print('\n{:*^30}\n'.format(' M062X completed! '))
            return None
        elif all(smi + '.out' in listdir('.') for smi in smiles):
            print('\n{:*^30}\n'.format(' M062X completed! '))
            return None
        else: print('\n{:*^30}\n'.format(' M062X incompleted! '))
        for smi in smiles:
            if smi + '.gjf' not in listdir('.'):
                self.Build_structure_from_smi(smi, self.conformer_method), print('This molecule has been constructed: {}'.format(smi))
                if 'conformer.py' in listdir('.'): remove('conformer.py')
        submit_smiles = [smi for smi in smiles if smi + '.gjf' in listdir('.') and smi + '.out' not in listdir('.')]
        if submit_smiles: self.Write_submit_script(submit_smiles, 'M062X')
    
    # To get submitted file of GFNFF.
    def Get_GFNFF_out(self):
        chdir('{}/gaussian_out/GFNFF'.format(self.path))
        smiles = [smi for smi in self.species if smi + '.dat' not in listdir('{}/chemkin_dat/{}'.format(self.path, 'M062X'))]
        self.Check_method_out(smiles, 'GFNFF')
        if not smiles:
            print('\n{:*^30}\n'.format(' GFNFF completed! '))
            return None
        elif all(smi + '.out' in listdir('.') for smi in smiles):
            print('\n{:*^30}\n'.format(' GFNFF completed! '))
            return None
        else: print('\n{:*^30}\n'.format(' GFNFF incompleted! '))
        for smi in smiles:
            if smi + '.xyz' not in listdir('.') and smi + '.out' in listdir('../M062X'):
                with open('../M062X/{}.out'.format(smi)) as f, open(smi + '.xyz', 'w') as p:
                    res = search('temp=200 [\s\S]+?Charge =  (\d) Multiplicity = (\d)[\s\S]+?Standard[\s\S]+?Z\n -+\n((.{70}\n)+?) -+', f.read())
                    coord = ['%s%16.6f%16.6f%16.6f\n'%(Atom(int(x[1])).symbol, x[3], x[4], x[5]) for x in [list(map(float, x.split())) for x in res.group(3).strip().split('\n')]]
                    p.write('{}\n\n{}'.format(len(coord), ''.join(coord))), print('This molecule has been converted: {}'.format(smi))       
        submit_smiles = [smi for smi in smiles if smi + '.xyz' in listdir('.') and smi + '.out' not in listdir('.')]
        if submit_smiles: self.Write_submit_script(submit_smiles, 'GFNFF')
    
    # To get submitted file of CCSD(T).
    def Get_CCSDT_out(self):
        chdir(self.path + '/gaussian_out/CCSDT')
        ele = {'C': 6, 'H': 1, 'O': 8}
        smiles = [smi for smi in self.micro_mol if smi + '.dat' not in listdir('{}/chemkin_dat/CCSDT'.format(self.path))]
        self.Check_method_out(smiles, 'CCSDT')
        if not smiles:
            print('\n{:*^30}\n'.format(' CCSDT completed! '))
            return None
        elif all(smi + '.out' in listdir('.') for smi in smiles):
            print('\n{:*^30}\n'.format(' CCSDT completed! '))
            return None
        else:
            print('\n{:*^30}\n'.format(' CCSDT incompleted! '))
        fore_smiles = [smi for smi in smiles if smi + '.gjf' not in listdir('.')]
        for smi in fore_smiles:
            if smi + '.out' in listdir('../M062X'):
                with open('../M062X/{}.out'.format(smi)) as f, open('../M062X/{}.gjf'.format(smi)) as p, open(smi + '.gjf', 'w') as q:
                    content = f.read()
                    newsmi = smi[:-2] if '-' in smi and smi[-1] in self.multiplicities else smi
                    res = search('temp=200 [\s\S]+?Charge =  (\d) Multiplicity = (\d)[\s\S]+?Standard[\s\S]+?Z\n -+\n((.{70}\n)+?) -+', content)
                    spin = int(res.group(2))
                    coord = '\n'.join('%s%16.6f%16.6f%16.6f'%(Atom(int(x[1])).symbol, x[3], x[4], x[5]) for x in [list(map(float, x.split())) for x in res.group(3).strip().split('\n')])
                    atoms = Counter(x.GetSymbol() for x in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAtoms())
                    electronic = sum([ele[k] * atoms[k] for k in atoms])
                    if spin > 1: q.write('***,{}\n\ngeometry={{\n{}\n}}\n\nbasis=vdz\n{{rhf;wf,{},spin,{},state,1}}\nuccsd(t)\n\nbasis=vtz\nrhf\nuccsd(t)\n\nbasis=vqz\nrhf\n'.format(smi, coord, electronic, spin - 1))
                    else: q.write('***,{}\n\ngeometry={{\n{}\n}}\n\nbasis=vdz\n{{hf;wf,{},spin,{},state,1}}\nccsd(t)\n\nbasis=vtz\nhf\nccsd(t)\n\nbasis=vqz\nhf\n'.format(smi, coord, electronic, spin - 1))
                    print('This molecule has been converted: {}'.format(smi))
        submit_smiles = [smi for smi in smiles if smi + '.gjf' in glob('*.gjf') and smi + '.out' not in glob('*.out')]
        if submit_smiles: self.Write_submit_script(submit_smiles, 'CCSDT')

    # To get submitted file of G4.
    def Get_G4_out(self, micro_mol):
        chdir(path + '/gaussian_out/G4')
        smiles = [smi for smi in micro_mol if smi + '.dat' not in listdir('{}/chemkin_dat/{}'.format(path, 'G4'))]
        self.Check_method_out(smiles, 'G4')
        if not smiles:
            print('\n{:*^30}\n'.format(' G4 completed! '))
            return None
        elif all(smi + '.out' in listdir('.') for smi in smiles):
            print('\n{:*^30}\n'.format(' G4 completed! '))
            return None
        else: print('\n{:*^30}\n'.format(' G4 incompleted! '))
        fore_smiles = [smi for smi in smiles if smi + '.gjf' not in listdir('.')]
        for smi in fore_smiles:
            if smi + '.out' in listdir('../M062X'):
                with open('../M062X/{}.out'.format(smi)) as f, open('../M062X/{}.gjf'.format(smi)) as p, open(smi + '.gjf', 'w') as q:
                    res = search('temp=200 [\s\S]+?Charge =  (\d) Multiplicity = (\d)[\s\S]+?Standard[\s\S]+?Z\n -+\n((.{70}\n)+?) -+', f.read())
                    symm = ''.join(findall('\s*(\w*symm=*\w*)', p.read()))
                    q.write('%nprocshared={}\n# G4=sp scf=fermi {}\n\n{}\n\n{} {}\n'.format(number_process, symm, smi, res.group(1), res.group(2)))
                    for x in [list(map(float, x.split())) for x in res.group(3).strip().split('\n')]: q.write('%s%16.6f%16.6f%16.6f\n'%(Atom(int(x[1])).symbol, x[3], x[4], x[5]))
                    q.write('\n'), print('This molecule has been converted: {}'.format(smi))
        submit_smiles = [smi for smi in smiles if smi + '.gjf' in glob('*.gjf') and smi + '.out' not in glob('*.out')]
        if submit_smiles: self.Write_submit_script(submit_smiles, 'G4')
    
    # To get submitted file of CBSQB3.
    def Get_CBSQB3_out(self, micro_mol):
        chdir(path + '/gaussian_out/CBSQB3')
        smiles = [smi for smi in micro_mol if smi + '.dat' not in listdir('{}/chemkin_dat/{}'.format(path, 'G4'))]
        self.Check_method_out(smiles, 'CBSQB3')
        if not smiles:
            print('\n{:*^30}\n'.format(' CBSQB3 completed! '))
            return None
        elif all(smi + '.out' in glob('*.out') for smi in smiles):
            print('\n{:*^30}\n'.format(' CBSQB3 completed! '))
            return None
        else: print('\n{:*^30}\n'.format(' CBSQB3 incompleted! '))
        fore_smiles = [smi for smi in smiles if smi + '.gjf' not in glob('*.gjf')]
        for smi in fore_smiles:
            if smi + '.out' in listdir('{}/gaussian_out/M062X/'.format(path)):
                with open('{}/gaussian_out/M062X/{}.out'.format(path, smi)) as f, open('{}/gaussian_out/M062X/{}.gjf'.format(path, smi)) as p, open(smi + '.gjf', 'w') as q:
                    res = search('temp=200 [\s\S]+?Charge =  (\d) Multiplicity = (\d)[\s\S]+?Standard[\s\S]+?Z\n -+\n((.{70}\n)+?) -+', f.read())
                    symm = ''.join(findall('\s*(\w*symm=*\w*)', p.read()))
                    q.write('%nprocshared={}\n# CBSQB3=sp scf=fermi {}\n\n{}\n\n{} {}\n'.format(number_process, symm, smi, res.group(1), res.group(2)))
                    for x in [list(map(float, x.split())) for x in res.group(3).strip().split('\n')]: q.write('%s%16.6f%16.6f%16.6f\n'%(Atom(int(x[1])).symbol, x[3], x[4], x[5]))
                    q.write('\n'), print('This molecule has been converted: {}'.format(smi))
        submit_smiles = [smi for smi in species if smi + '.gjf' in glob('*.gjf') and smi + '.out' not in glob('*.out')]
        if submit_smiles: self.Write_submit_script(submit_smiles, 'CBSQB3')
    
    # To build 3D coordinates from SMILES.
    def Build_structure_from_smi(self, smi, structure_method):
        print('\nThis molecule of structure is being constructed: {}'.format(smi))
        if 'conformer.py' not in listdir('.'): copy('{}/src/conformer.py'.format(self.path), '.')
        newsmi = smi[:-2] if '-' in smi and smi[-1] in self.multiplicities else smi
        mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))
        AllChem.EmbedMolecule(mol), AllChem.MMFFOptimizeMolecule(mol), Chem.MolToMolFile(mol, 'rdkit.mol'), system('obabel rdkit.mol -O rdkit.gjf')
        exit()
        with open('input.txt', 'w') as f: f.write('{}\n{}\n{}\n{}\n{}\n{}'.format(newsmi, smi, self.number_process, structure_method, self.spin_method, self.key_words))
        system('python conformer.py >conformer.txt')
        remove('rdkit.mol'), remove('rdkit.gjf'), remove('input.txt'), remove('conformer.txt')
        with open(smi + '.gjf') as f: content = f.readlines()
        if '*' in ''.join(content): self.Build_structure_from_smi(smi, structure_method)
        else:
            try:
                lines = findall('\n\d\s+\d\n([\s\S]+?)\n\n', ''.join(content))[0].split('\n')
                distances = Atoms([Atom(e, (float(x), float(y), float(z))) for e, (x, y, z) in [(x.split()[0], x.split()[-3:]) for x in lines]]).get_all_distances()
                if sum(1 if 0.1 < x <0.7 else 0 for x in chain(*distances)) >= 2:
                    for i, v in enumerate(lines): content[7 + i] = v.split()[0] + ' '.join('%16.5f'%(float(x) * 1.8) for x in v.split()[-3:]) + '\n'
                    with open(smi + '.gjf', 'w') as f: f.write(''.join(content))
            except:
                self.Build_structure_from_smi(smi, structure_method)
    
    # To check whether the output file is complete and give corresponding solutions.
    def Check_method_out(self, smiles, method):
        print('\nCheck calculations of {} ...\n'.format(method))
        for x in chain(*[glob('*.chk'), glob('job*'), glob('nohup.out'), glob('fort.7')], glob('molpro*')): remove(x)
        post_smiles = [smi for smi in smiles if smi + '.out' in listdir('.')]
        manually_process = []
        if method == 'GFNFF':
            for smi in post_smiles:
                with open(smi + '.out', 'rb') as f:
                    f.seek(-500, 2)
                    endinfo = f.read().decode()
                if 'normally' not in endinfo: print('\nNot over file: {}.out\nProcessed: Perform calculation again.\n'.format(smi)), remove(smi + '.out')
        elif method == 'M062X':
            copy('{}/src/Shermo.exe'.format(self.path), '.'), copy('{}/src/settings.ini'.format(self.path), '.')
            scale = findall('scale=(\s*\d+?\.\d+)', self.key_words)[0] if findall('scale=(\s*\d+?\.\d*?)', self.key_words) else '1'
            with open('settings.ini') as f: lines = f.read()
            with open('settings.ini', 'w') as f: f.write(sub('(sclS=)(\s*\S+)', '\\1 {}'.format(scale), lines))
            for smi in post_smiles:
                with open(smi + '.out', 'rb') as f:
                    content = f.read().decode()
                    f.seek(-500, 2)
                    endinfo = f.read().decode()
                with open(smi + '.gjf') as f: lines = f.readlines()
                spin = int(lines[6][-2])
                if 'Error' in endinfo:
                    print('\nError end file: {}.out\nError messages:\n{}'.format(smi, '\n'.join(endinfo.split('\n')[-6:-4]))), remove(smi + '.out')
                    if 'Problem with the number of degrees of freedom' in endinfo:
                        if any(x in smi for x in ['=C=', '#', ']C=', ']=', '=[']):
                            print('Processed: Remove hinderedrotor.\n')
                            coord = findall('(\S+)\s+(\S+)\s+(\S+)\n', search('[\s\S]+Input orientation[\s\S]+?Z\n.+\n([\s\S]+?\n) --+[\s\S]+?$', content).group(1))
                            for i, v in enumerate(coord): lines[7 + i]=sub('(\s*[a-zA-Z]+)(\s+\S.+)\n', '\\1%16s%16s%16s\n'%v, lines[7 + i])
                            lines[2] = lines[2].replace('=hinderedrotor', '').replace(', hinderedrotor', '').replace('hinderedrotor, ', '')
                            with open(smi + '.gjf', 'w') as f: f.write(''.join(lines).replace('(readfc, hinderedrotor)', 'readfc'))
                        else:
                            structure_method = 1 if '_' not in lines[4] else int(lines[4][-2])
                            if structure_method < 3: print('Processed: Change conformer construction.\n'), self.Build_structure_from_smi(smi, str(structure_method + 1))
                            else: manually_process.append(smi), print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
                    elif 'Linear angle in Tors' in endinfo or 'FormBX had a problem.' in endinfo:
                        res = findall('opt\S*\s', lines[2])[0]
                        if 'cartesian' not in res:
                            print('Processed: Change cartesian optimazition.\n')
                            if '=' in res and '(' not in res: lines[2] = lines[2].replace(res, res.replace('=', '=(cartesian, ').replace(' ', ') '))
                            else: lines[2] = lines[2].replace(res, res.replace('opt ', 'opt=cartesian ').replace('(', '(cartesian, '))
                            with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                        else: manually_process.append(smi), print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
                    elif 'l9999' in endinfo:
                        max_force = float(search('[\s\S]+\n Maximum Force\s+(\S+)\s+[\s\S]+?$', content).group(1))
                        if max_force < 0.000200:
                            coord = findall('(\S+)\s+(\S+)\s+(\S+)\n', search('[\s\S]+Input orientation[\s\S]+?Z\n.+\n([\s\S]+?\n) --+[\s\S]+?$', content).group(1))
                            for i, v in enumerate(coord): lines[7 + i] = sub('(\s*[a-zA-Z]+)(\s+\S.+)\n', '\\1%16s%16s%16s\n'%v, lines[7 + i])
                            print('Processed: Increase grid accuracy.\n')
                            if 'int=ultrafine' not in lines[2] and 'int=superfine' not in lines[2]:
                                version = search('Normal termination of Gaussian (\d\d)', content).group(1)
                                if version == '09': lines[2] = lines[2].replace('\n', ' int=ultrafine\n')
                                else: lines[2] = lines[2].replace('\n', ' int=superfine\n')
                            with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                        else:
                            if 'gdiis' not in lines[2]:
                                print('Processed: Change optimazition algorithm.\n')
                                lines[2] = lines[2].replace('opt ', 'opt=gdiis ').replace('opt=(', 'opt=(gdiis, ').replace('scf=fermi', '')
                                if 'gdiis' not in lines[2]: lines[2] = sub('opt=(\S+)', 'opt=(gdiis, \\1)', lines[2])
                                with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                            else:
                                structure_method = 1 if '_' not in lines else int(lines[4][-2])
                                if structure_method < 3: print('Processed: Change conformer construction.\n'), self.Build_structure_from_smi(smi, str(structure_method + 1))
                                else: manually_process.append(smi), print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
                    else: manually_process.append(smi), print('Please manually choose applicable key words to elimate error for {}'.format(smi))
                elif 'Normal' not in endinfo: print('\nNot over file: {}.out\nProcessed: Perform calculation again.\n'.format(smi)), remove(smi + '.out')
                elif 'imaginary frequencies' in content:
                    coord = findall('(\S+)\s+(\S+)\s+(\S+)\n', search('temp=200[\s\S]+? Standard orientation[\s\S]+?Z\n.+\n([\s\S]+?\n) --+', content).group(1))
                    res = search(' Frequencies --(\s.+\n)[\s\S]+?Z(\n[\s\S]+?)\n\s{15}', content)
                    imag, offset = [x for x in res.group(1).split() if float(x) < 0], findall('\n\s+\d+\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)', res.group(2))
                    print('\nImaginary frequencies file: {}.out\nImaginary frequencies: {}\nProcessed: Apply vibration displacement.\n'.format(smi, ', '.join(imag))), remove(smi + '.out')
                    for i, (x, y) in enumerate(zip(coord, offset)): lines[7 + i] = sub('(\s*[a-zA-Z]+)(\s.+)\n', '\\1 %16.6f%16.6f%16.6f\n'%tuple(float(m) + 0.15 * float(n) for m, n in zip(x, y)), lines[7 + i])
                    if 'int=ultrafine' not in lines[2] and 'int=superfine' not in lines[2]:
                        version = search('Normal termination of Gaussian (\d\d)', content).group(1)
                        if version == '09': lines[2] = lines[2].replace('\n', ' int=ultrafine\n')
                        else: lines[2] = lines[2].replace('\n', ' int=superfine\n')
                    with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                elif 'NO' in [x.split()[-1] for x in findall('Converged\?\n([\s\S]+)\n Predicted', search(' Thermochemistry ([\s\S]+)Normal termination', content).group(1))[0].split('\n')][:2]:
                    print('\nFreq no convergence: {}.out\nProcessed: Continue to optimize.\n'.format(smi)), remove(smi + '.out')
                    coord = findall('(\S+)\s+(\S+)\s+(\S+)\n', search('temp=200[\s\S]+? Standard orientation[\s\S]+?Z\n.+\n([\s\S]+?\n) --+', content).group(1))
                    if 'int=superfine' not in lines[2]: lines[2] = lines[2].replace('\n', ' int=superfine\n')
                    for i, v in enumerate(coord): lines[7 + i] = sub('(\s*[a-zA-Z]+)(\s+\S.+)\n', '\\1%16s%16s%16s\n'%v, lines[7 + i])
                    with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                elif spin>1:
                    res = list(map(float, findall('S\*\*2 before annihilation\s+(\S+?),.+?\s+(\S+)\n', content)[-1]))
                    if abs(res[0] - res[1]) > res[1]* 0.1:
                        print('\nSpin contamination exceeds 10% of S**2: {}.out\nProcessed: Change UM062X to ROM062X.\n'.format(smi)), remove(smi + '.out')
                        lines=''.join(lines).replace('UM062X','ROM062X')
                        with open(smi + '.gjf') as f: f.write(lines)
                else:
                    S_gaussian = float(search('Mol-Kelvin\n Total\s+(\S+)\s+(\S+)\s+(\S+)\n', content).group(3))
                    system('Shermo.exe {}.out >temp.txt'.format(smi))
                    with open('scan_SCq.txt') as f: S_shermo = float(f.readlines()[3].split()[2])
                    if abs(S_gaussian - S_shermo) > 0.5:
                        print('\nNot corrected point group: {}.out\nProcessed: Decrease structure symmetry.\n'.format(smi)), remove(smi + '.out')
                        coord = findall('(\S+)\s+(\S+)\s+(\S+)\n', search('temp=200[\s\S]+? Standard orientation[\s\S]+?Z\n.+\n([\s\S]+?\n) --+', content).group(1))
                        if 'loose' not in lines[2]: lines[2] = lines[2].replace('\n', ' symm=loose\n')
                        else: lines[2] = lines[2].replace('=loose', '=veryloose').replace('nosymm', 'symm=veryloose').replace('(loose)', '(veryloose)')
                        for i, v in enumerate(coord): lines[7 + i] = sub('(\s*[a-zA-Z]+)(\s+\S.+)\n', '\\1%16s%16s%16s\n'%v, lines[7 + i])
                        with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                    remove('scan_SCq.txt'), remove('scan_UHG.txt'), remove('temp.txt')
            remove('Shermo.exe'), remove('settings.ini')
            if manually_process: print('\n\n!!! Manually Processed smiles. !!!\n{}\n'.format('\t'.join(manually_process)))
        elif method == 'CCSDT':
            for smi in post_smiles:
                with open(smi + '.out') as f: content = f.read()
                spin = int(search(',(\d),state', content)[1])
                res = findall('Spin contamination.+?\s+?(\S+?)\n',content)
                spin_contamination = float(res[-1]) if res else 0
                if 'Molpro calculation terminated' not in content[-100:] and 'fail' not in content[-200:]:
                    print('\nNot over file: {}.out\nProcessed: Perform calculation again.\n'.format(smi)), remove(smi + '.out')
                elif 'fail' in content[-200:]:
                    manually_process.append(smi), print('Please choose applicable key words to elimate error for {}\n'.format(smi)), remove(smi + '.out')
                elif spin > 0 and spin_contamination > 0.1 * 0.25 * spin * (spin + 2):
                    print('\nSpin contamination exceeds 10% of S**2: {}.out\nProcessed: Change UCCSD(T) to RCCSD(T).\n'.format(smi)), remove(smi + '.out')
                    with open(smi + '.gjf') as f: content = f.read()
                    content = content.replace('uccsd(t)', 'rccsd(t)')
                    with open(smi + '.gjf', 'w') as f: f.write(content)
                else:
                    T1 = float(findall('T1 diagnostic:\s*(\S+)\n', content)[-1])
                    if T1 > 0.04 and 'rccsd(t)' not in content:
                        print('\nT1 Diagnostic exceeds 0.04: {}.out\nProcessed: Change UCCSD(T) to RCCSD(T).\n'.format(smi)), remove(smi + '.out')
                        with open(smi + '.gjf') as f: content = f.read()
                        content = content.replace('uccsd(t)', 'rccsd(t)')
                        with open(smi + '.gjf', 'w') as f: f.write(content)
                    elif T1 > 0.04 and 'rccsd(t)' in content:
                        print('\nT1 Diagnostic exceeds 0.04: {}.out\n\n'.format(smi)), remove(smi + '.out')
                        manually_process.append(smi), print('Please choose applicable solution to elimate this error for {}\n'.format(smi))
        else:
            for smi in post_smiles:
                with open(smi + '.out', 'rb') as f: endinfo = f.read().decode()[-500:]
                if 'Error' in endinfo:
                    print('\nError end file: {}.out\nError messages:\n{}'.format(smi, '\n'.join(endinfo.split('\n')[-6:-4]))), remove(smi + '.out')
                    with open(smi + '.gjf') as f: lines = f.readlines()
                    if 'l502' in endinfo and 'maxcyc=300' not in lines[1]:
                        print('Processed: Increase iteration steps.\n')
                        lines[1] = lines[1].replace('scf=fermi', 'scf=(fermi,maxcyc=300)')
                        with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                    elif 'xqc' not in lines[1]:
                        print('Processed: Change to quadratic convergence.\n')
                        lines[1] = lines[1].replace('scf=(fermi,maxcyc=300)', 'scf=(fermi,xqc,maxcyc=300)')
                        with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                    else: manually_process.append(smi), print('Please choose applicable key words to elimate error for {}\n'.format(smi))
                elif 'Normal' not in endinfo: print('\nNot over file: {}.out\nProcessed: Perform calculation again.\n'.format(smi)), remove(smi + '.out')
    
    # To write submitted scripts for linux systems.
    def Write_submit_script(self, smiles, method):
        if self.submit_type not in ['1', '2']: print('!!! This type has not been considered. Default type of nohup has been applied. !!!')
        gjf=["'{}.gjf'".format(smi) if set(smi) & set(['(', '[', '=']) else smi + '.gjf' for smi in smiles]
        makedirs('{}/submitted_inp/{}'.format(self.path, method), True)
        if method == 'GFNFF':
            xyz=["'{}.xyz'".format(smi) if set(smi) & set(['(', '[', '=']) else smi + '.xyz' for smi in smiles]
            with open('{}/submitted_inp/GFNFF/{}'.format(self.path, self.submit_scripts[self.submit_type]), 'w', newline = '\n') as f:
                f.write('#!/bin/bash\n\n')
                for i, v in enumerate(list(get_chunks(sample(xyz, len(xyz)), ceil(len(xyz) / int(self.number_task))))):
                    with open('{}/submitted_inp/GFNFF/job{}.sh'.format(self.path, i + 1), 'w', newline = '\n') as p:
                        p.write('#!/bin/bash\n\n')
                        if self.submit_type == '2':
                            f.write('mkdir {0}\ncp job{0}.sh {0}\ncp {1} {0}\ncd {0}\nqsub *.sh\ncd ..\n\n'.format(i + 1, ' '.join(v)))
                            p.write('#PBS -N {}-job{}\n#PBS -l nodes=1:ppn={}\n#PBS -q medium\n#PBS -j eo\ncd $PBS_O_WORKDIR\n\n'.format(method, i + 1, self.number_process))
                            for x in v:
                                smi = x[:-4] if x[0] != '\'' else x[1:-5]
                                with open('../M062X/' + smi + '.gjf') as q: chrg, spin = map(int, findall('\n\n(\d)\s+(\d)\s*\n', q.read())[0])
                                p.write('crest {0} >{1} -gfnff -v4 -quick -ewin {2} -T {3} -chrg {4} -uhf {5}\nwait\ncp -f crest_best.xyz ../{0}\ncp {1} ..\n'.format(x, x.replace('xyz', 'out'), 50, self.number_process, chrg, spin - 1))
                        else:
                            f.write('mkdir {0}\ncp job{0}.sh {0}\ncp {1} {0}\ncd {0}\nchmod +x *.sh\nnohup ./*.sh >nohup.out &\ncd ..\n\n'.format(i + 1, ' '.join(v)))
                            for x in v:
                                smi = x[:-4] if x[0] != '\'' else x[1:-5]
                                with open('../M062X/' + smi + '.gjf') as q: chrg, spin = map(int, findall('\n\n(\d)\s+(\d)\s*\n', q.read())[0])
                                p.write('nohup crest {0} >{1} -gfnff -v4 -quick -ewin {2} -T {3} -chrg {4} -uhf {5} &\nwait\ncp -f crest_best.xyz ../{0}\ncp {1} ..\n'.format(x, x.replace('xyz', 'out'), 50, self.number_process, chrg, spin - 1))
                if self.submit_type != '2': f.write('wait\nrm -rf {}\n\n'.format(' '.join(map(str, range(1, i + 2)))))
            for smi in smiles: copy(smi + '.xyz', '{}/submitted_inp/{}'.format(self.path, method))
        elif method == 'CCSDT':
            for i, v in enumerate(list(get_chunks(sample(gjf, len(gjf)), ceil(len(gjf) / int(self.number_task))))):
                with open('{}/submitted_inp/{}/job{}.sh'.format(self.path, method, i + 1), 'w', newline = '\n') as f:
                    f.write('#!/bin/bash\n\n')
                    if self.submit_type == '2':
                        f.write('#PBS -N {}-job{}\n#PBS -l nodes=1:ppn={}\n#PBS -q medium\n#PBS -j eo\ncd $PBS_O_WORKDIR\n\n'.format(method, i + 1, self.number_process))
                        for x in v: f.write('molpro -n {} -m 300m <{} >{}\nwait\nsleep 10'.format(self.number_process, x, x.replace('gjf', 'out')))
                    else:
                        for x in v: f.write('nohup molpro -n {} -m 300m <{} >{} &\nwait\nsleep 10\n'.format(self.number_process, x, x.replace('gjf', 'out')))
        else:
            for i, v in enumerate(list(get_chunks(sample(gjf, len(gjf)), ceil(len(gjf) / int(self.number_task))))):
                with open('{}/submitted_inp/{}/job{}.sh'.format(self.path, method, i + 1), 'w', newline = '\n') as f:
                    f.write('#!/bin/bash\n\n')
                    if self.submit_type == '2':
                        f.write('#PBS -N {}-job{}\n#PBS -l nodes=1:ppn={}\n#PBS -q medium\n#PBS -j eo\ncd $PBS_O_WORKDIR\n\n'.format(method, i + 1, self.number_process))
                        for x in v: f.write('g09 <{} >{}\nwait\n'.format(x, x.replace('gjf', 'out')))
                    else:
                        for x in v: f.write('nohup g09 <{} >{} &\nwait\n'.format(x, x.replace('gjf', 'out')))
        if method != 'GFNFF':
            copy('{}/src/{}'.format(self.path, self.submit_scripts[self.submit_type]), '{}/submitted_inp/{}'.format(self.path, method))
            for x in glob('job*'): move(x, '{}/submitted_inp/{}'.format(self.path, method))
            for smi in smiles: copy(smi + '.gjf', '{}/submitted_inp/{}'.format(self.path, method))

    def Get_submit_out(self):
        print('\nProcess submited files ...\n')
        self.Get_M062X_out()
        self.Get_GFNFF_out()
        eval('self.Get_{}_out()'.format(self.high_level))
