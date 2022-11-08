# -*- coding: utf-8 -*-
from glob import glob
from numpy import array
from os import listdir, remove
from re import findall, search, sub

import tdoc._submit


class Check(object):

    """ To initiate parameters for Check. """
    def __init__(self, para, smiles, method):
        self.para = para
        self.smiles = smiles
        self.method = method


    
    """ To check M062X calculations. """
    def check_M062X_out(self, finished_smiles):
        manually_process = []
        for smi in finished_smiles:
            with open(smi + '.out', 'rb') as f:
                content = f.read().decode()
                f.seek(-500, 2)
                endinfo = f.read().decode()
            with open(smi + '.gjf') as f: lines = f.readlines()
            spin = int(lines[6][-2])

            # Deal with unfinished calculations.
            if 'Normal' not in endinfo and 'Error' not in endinfo:
                remove(smi + '.out')
                print('\nNot over file: {}\nProcessed: Perform calculation again.\n'.format(smi))

            # Deal with runtime errors.
            elif 'Error' in endinfo:
                remove(smi + '.out')
                print('\nError end file: {}\nError messages:\n{}'.format(smi, '\n'.join(endinfo.split('\n')[-6:-4])))

                # Process the freedom of hindered rotors.
                if 'Problem with the number of degrees of freedom' in endinfo:
                    print('Processed: Remove hinderedrotor.\n')
                    aim_content = search('[\s\S]+Input orientation[\s\S]+?Z\n.+\n([\s\S]+?\n) --+[\s\S]+?$', content).group(1)
                    coord = findall('(\S+)\s+(\S+)\s+(\S+)\n', aim_content)
                    for i, v in enumerate(coord):
                        lines[7 + i] = sub('(\s*[a-zA-Z]+)(\s+\S.+)\n', '\\1%16s%16s%16s\n'%v, lines[7 + i])
                    lines[2] = lines[2].replace('=hinderedrotor', '')
                    with open(smi + '.gjf', 'w') as f:
                        f.write(''.join(lines).replace('(readfc,hinderedrotor)', 'readfc'))

                # Process the problem of internal coordinates.
                elif 'Linear angle in Tors' in endinfo or 'FormBX had a problem' in endinfo:
                    if 'cartesian' not in lines[2]:
                        print('Processed: Change to cartesian optimazition.\n')
                        lines[2] = lines[2].replace('opt ', 'opt=cartesian ').replace('opt=tight', 'opt=(tight,cartesian)')
                        with open(smi + '.gjf', 'w') as f:
                            f.write(''.join(lines))
                    else:
                        manually_process.append(smi)
                        print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))

                # Process the problem of optimization convergence.
                elif 'l9999' in endinfo:
                    aim_content = search('[\s\S]+Input orientation[\s\S]+?Z\n.+\n([\s\S]+?\n) --+[\s\S]+?$', content).group(1)
                    coord = findall('(\S+)\s+(\S+)\s+(\S+)\n', aim_content)
                    for i, v in enumerate(coord):
                        lines[7 + i] = sub('(\s*[a-zA-Z]+)(\s+\S.+)\n', '\\1%16s%16s%16s\n'%v, lines[7 + i])
                    max_force = float(search('[\s\S]+\n Maximum Force\s+(\S+)\s+[\s\S]+?$', content).group(1))
                    if self.para['gaussian_version'] == 'g09' and max_force < 0.000500:
                        if 'int=ultrafine' not in lines[2]:
                            print('Processed: Increase grid accuracy.\n')
                            lines[2] = lines[2].replace('\n', ' int=ultrafine\n')
                        elif 'opt=tight' in lines[2]:
                            print('Processed: Decrease convergence accuracy.\n')
                            lines[2] = lines[2].replace('opt=tight', 'opt').replace('opt=(tight,cartesian)', 'opt=cartesian')
                        else:
                            print('Processed: Continue to optimize.\n')
                        with open(smi + '.gjf', 'w') as f:
                            f.write(''.join(lines))
                    elif self.para['gaussian_version'] == 'g16' and max_force < 0.000500:
                        if 'int=superfine' not in lines[2]:
                            print('Processed: Increase grid accuracy.\n')
                            lines[2] = lines[2].replace('\n', ' int=superfine\n')
                        elif 'opt=tight' in lines[2]:
                            print('Processed: Decrease convergence accuracy.\n')
                            lines[2] = lines[2].replace('opt=tight', 'opt').replace('opt=(tight,cartesian)', 'opt=cartesian')
                        else:
                            print('Processed: Continue to optimize.\n')
                        with open(smi + '.gjf', 'w') as f:
                            f.write(''.join(lines))
                    else:
                        structure_method = 1 if '_' not in lines[4] else int(lines[4][-2])
                        if structure_method == 1: 
                            print('Processed: Change to conformer construction.\n')
                            tdoc._submit.Submit(self.para).build_structure_from_smi(smi, 2)
                            print('\nThis molecule has been constructed: {}'.format(smi))
                        else: 
                            manually_process.append(smi)
                            print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
                    
                # Process the problem of bad 3D structure.
                elif 'l202' in endinfo:
                    structure_method = 1 if '_' not in lines[4] else int(lines[4][-2])
                    if structure_method == 1:
                        print('Processed: Change to conformer construction.\n')
                        tdoc._submit.Submit(self.para).build_structure_from_smi(smi, 2)
                    else:
                        manually_process.append(smi)
                        print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
                    
                # Process the problem of SCF convergence.
                elif 'l502' in endinfo:
                    if 'maxcyc=300' not in lines[2]:
                        print('Processed: Increase iteration steps.\n')
                        lines[2] = lines[2].replace('scf=fermi', 'scf=(fermi,maxcyc=300)')
                        with open(smi + '.gjf', 'w') as f:
                            f.write(''.join(lines))
                    elif 'vshift=300' not in lines[2]:
                        print('Processed: Apply shift level.\n')
                        lines[2] = lines[2].replace('scf=(fermi,maxcyc=300)', 'scf=(fermi,maxcyc=300,vshift=300)')
                        with open(smi + '.gjf', 'w') as f:
                            f.write(''.join(lines))
                    elif 'xqc' not in lines[2]:
                        print('Processed: Change to quadratic convergence.\n')
                        lines[2] = lines[2].replace('scf=(fermi,maxcyc=300,vshift=300)', 'scf=(fermi,xqc,maxcyc=300,vshift=300)')
                        with open(smi + '.gjf', 'w') as f:
                            f.write(''.join(lines))
                    else:
                        manually_process.append(smi)
                        print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
                    
                # Process the unknown problem by users manually.
                else:
                    manually_process.append(smi)
                    print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
                
            # Deal with imaginary frequencies.
            elif 'imaginary frequencies' in content:
                aim_content = search('temp=200[\s\S]+? Standard orientation[\s\S]+?Z\n.+\n([\s\S]+?\n) --+', content).group(1)
                coord = findall('(\S+)\s+(\S+)\s+(\S+)\n', aim_content)
                res = search(' Frequencies --(\s.+\n)[\s\S]+?Z(\n[\s\S]+?)\n\s{15}', content)
                imag = [x for x in res.group(1).split() if float(x) < 0]
                offset = findall('\n\s+\d+\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)', res.group(2))
                remove(smi + '.out')
                print('\nImaginary frequencies file: {}\nImaginary frequencies: {}\nProcessed: Apply vibration displacement.\n'.format(smi, ', '.join(imag)))
                for i, (x, y) in enumerate(zip(coord, offset)):
                    lines[7 + i] = sub('(\s*[a-zA-Z]+)(\s.+)\n', '\\1 %16.6f%16.6f%16.6f\n'%tuple(float(m) + 0.3 * float(n) for m, n in zip(x, y)), lines[7 + i])
                if self.para['gaussian_version'] == 'g16':
                    if 'int=superfine' not in lines[2]:
                        print('Processed: Increase grid accuracy.\n')
                        lines[2] = lines[2].replace('\n', ' int=superfine\n')
                    else:
                        print('Processed: Increase convergence accuracy.\n')
                        lines[2] = lines[2].replace(' opt ', ' opt=tight ').replace('opt=cartesian', 'opt=(tight,cartesian)')

                else:
                    if 'int=ultrafine' not in lines[2]:
                        print('Processed: Increase grid accuracy.\n')
                        lines[2] = lines[2].replace('\n', ' int=ultrafine\n')
                    else:
                        print('Processed: Increase convergence accuracy.\n')
                        lines[2] = lines[2].replace(' opt ', ' opt=tight ').replace('opt=cartesian', 'opt=(tight,cartesian)')
                with open(smi + '.gjf', 'w') as f:
                    f.write(''.join(lines))
                
            # Deal with the spin contamination.
            elif spin > 1:
                res = list(map(float, findall('before annihilation\s+(\S+?),\s+after\s+(\S+?)\n', content)[-1]))
                if abs(res[0] - res[1]) > res[1] * 0.1:
                    remove(smi + '.out')
                    opt_level = findall(' (\w+?)/',lines[2])[0]
                    print('\nSpin contamination exceeds 10% of S**2: {0}\nProcessed: Change {1} to RO{1}.\n'.format(smi, opt_level))
                    lines=''.join(lines).replace('{}'.format(opt_level), 'RO{}'.format(opt_level))
                    with open(smi + '.gjf') as f:
                        f.write(lines)
            
        # Show manually processed molecules.
        if manually_process:
            print('\n\n!!! Manually Processed smiles. !!!\n{}\n'.format('\t'.join(manually_process)))



    """ To check GFNFF calculations. """
    def check_GFNFF_out(self, finished_smiles):
        for smi in finished_smiles:
            with open(smi + '.out', 'rb') as f:
                content = f.read().decode()

            # Deal with unfinished calculations.
            if 'normally' not in content:
                remove(smi + '.out')
                print('\nNot over file: {}\nProcessed: Perform calculation again.\n'.format(smi))


    
    """ To check whether the calculations of related methods are completed. """
    def check_method_completed(self):
        unfinish_smiles = [smi for smi in self.smiles if smi + '.out' not in listdir('.')]
        if unfinish_smiles:
            print('\n{:*^30}\n'.format(' {} incompleted! '.format(self.method)))
        else:
            print('\n{:*^30}\n'.format(' {} completed! '.format(self.method)))
            return None


    """ To check CCSDT calculations. """
    def check_CCSDT_out(self, finished_smiles):
        manually_process = []
        for smi in finished_smiles:
            with open(smi + '.out') as f:
                content = f.read()

            # Deal with unfinished calculations.
            if 'Molpro calculation terminated' not in content and 'fail' not in content:
                remove(smi + '.out')
                print('\nNot over file: {}\nProcessed: Perform calculation again.\n'.format(smi))

            # Deal with insufficient memory.
            elif 'fail' in content and 'Insufficient memory to allocate' in content:
                remove(smi + '.out')
                print('\nError end file: {}\nError messages:\nInsufficient memory.'.format(smi))
                with open(smi +'.gjf') as f:
                    input_content = f.read()
                newsmi = smi[:-2] if '-' in smi and smi[-1] in self.para['multiplicities'] else smi
                atoms = Counter(x.GetSymbol() for x in Chem.AddHs(Chem.MolFromSmiles(newsmi)).GetAtoms())
                mem = int(findall('memory,(\d+),m', input_content)[0])
                added_mem = (atoms.get('C', 0) + atoms.get('O', 0) + atoms.get('N', 0)) * 100
                mem = mem + added_mem
                input_content = sub('(?<=memory,)(\d+)', str(mem), input_content)
                with open(smi + '.gjf', 'w') as f:
                    f.write(input_content)
                print('Processed: Increase input memory.\n')

            # Deal with the problem of SCF convergence.
            elif 'fail' in content and 'No convergence' in content:
                remove(smi + '.out')
                print('\nError end file: {}\nError messages:\nNo convergence.'.format(smi))
                with open(smi +'.gjf') as f:
                    input_content = f.read()
                if 'shift' not in input_content:
                    input_content = input_content.replace('maxit,80','maxit,80;shift,-0.3,0.3')
                    with open(smi + '.gjf', 'w') as f:
                        f.write(input_content)
                    print('Processed: Apply shift level.\n')
                elif 'shift,-0.3,0.3' in input_content:
                    input_content = input_content.replace('shift,-0.3,0.3','shift,0.3,-0.3')
                    with open(smi + '.gjf', 'w') as f:
                        f.write(input_content)
                    print('Processed: Apply shift level.\n')
                elif 'shift,0.3,-0.3' in input_content:
                    input_content = input_content.replace('shift,0.3,-0.3','shift,0.3,-0.5')
                    with open(smi + '.gjf', 'w') as f:
                        f.write(input_content)
                    print('Processed: Apply shift level.\n') 
                elif 'shift,0.3,-0.5' in input_content:
                    input_content = input_content.replace('shift,0.3,-0.5','shift,-0.7,-0.4')
                    with open(smi + '.gjf', 'w') as f:
                        f.write(input_content)
                    print('Processed: Apply shift level.\n')
                else:
                    manually_process.append(smi)
                    print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
                
            # Deal with the unknown problem by users manually.
            elif 'fail' in content:
                remove(smi + '.out')
                manually_process.append(smi)
                print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
                
            # Deal with the problem of multireference effects.
            else:
                T1 = round(float(findall('T1 diagnostic:\s*(\S+)\n', content)[-1]),3)
                if T1 > 0.05:
                    remove(smi + '.out')
                    print('\nT1 diagnostic exceeds 0.05: {}\n\n'.format(smi))
                    manually_process.append(smi)
                    print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))

        # Show manually processed molecules.
        if manually_process:
            print('\n\n!!! Manually Processed smiles. !!!\n{}\n'.format('\t'.join(manually_process)))


    
    """ To check G4 calculations. """
    def check_G4_out(self, finished_smiles):
        manually_process = []
        for smi in finished_smiles:
            with open(smi + '.out', 'rb') as f:
                endinfo = f.read().decode()[-500:]

            # Deal with unfinished calculations.
            if 'Normal' not in endinfo:
                remove(smi + '.out')
                print('\nNot over file: {}\nProcessed: Perform calculation again.\n'.format(smi))

            # Deal with runtime errors.
            elif 'Error' in endinfo:
                remove(smi + '.out')
                print('\nError end file: {}\nError messages:\n{2}'.format(smi, '\n'.join(endinfo.split('\n')[-6:-4])))
                with open(smi + '.gjf') as f: lines = f.readlines()

                # Process the problem of SCF convergence.
                if 'l502' in endinfo and 'maxcyc=300' not in lines[1]:
                    print('Processed: Increase iteration steps.\n')
                    lines[1] = lines[1].replace('scf=fermi', 'scf=(fermi,maxcyc=300)')
                    with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                elif 'l502' in endinfo and 'vshift=300' not in lines[1]:
                    print('Processed: Apply shift level.\n')
                    lines[1] = lines[1].replace('scf=(fermi,maxcyc=300)', 'scf=(fermi,maxcyc=300,vshift=300)')
                    with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                elif 'l502' in endinfo and 'xqc' not in lines[1]:
                    print('Processed: Change to quadratic convergence.\n')
                    lines[1] = lines[1].replace('scf=(fermi,maxcyc=300,vshift=300)', 'scf=(fermi,xqc,maxcyc=300,vshift=300)')
                    with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                    
                # Process the unknown problem by users manually.
                else:
                    manually_process.append(smi)
                    print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
            
        # Show manually processed molecules.   
        if manually_process: print('\n\n!!! Manually Processed smiles. !!!\n{}\n'.format('\t'.join(manually_process)))



    """ To check CBSQB3 calculations. """
    def check_CBSQB3_out(self, finished_smiles):
        manually_process = []
        for smi in finished_smiles:
            with open(smi + '.out', 'rb') as f:
                endinfo = f.read().decode()[-500:]

            # Deal with unfinished calculations.
            if 'Normal' not in endinfo:
                remove(smi + '.out')
                print('\nNot over file: {}\nProcessed: Perform calculation again.\n'.format(smi))

            # Deal with runtime errors.
            elif 'Error' in endinfo:
                remove(smi + '.out')
                print('\nError end file: {}\nError messages:\n{2}'.format(smi, '\n'.join(endinfo.split('\n')[-6:-4])))
                with open(smi + '.gjf') as f: lines = f.readlines()

                # Process the problem of SCF convergence.
                if 'l502' in endinfo and 'maxcyc=300' not in lines[1]:
                    print('Processed: Increase iteration steps.\n')
                    lines[1] = lines[1].replace('scf=fermi', 'scf=(fermi,maxcyc=300)')
                    with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                elif 'l502' in endinfo and 'vshift=300' not in lines[1]:
                    print('Processed: Apply shift level.\n')
                    lines[1] = lines[1].replace('scf=(fermi,maxcyc=300)', 'scf=(fermi,maxcyc=300,vshift=300)')
                    with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
                elif 'l502' in endinfo and 'xqc' not in lines[1]:
                    print('Processed: Change to quadratic convergence.\n')
                    lines[1] = lines[1].replace('scf=(fermi,maxcyc=300,vshift=300)', 'scf=(fermi,xqc,maxcyc=300,vshift=300)')
                    with open(smi + '.gjf', 'w') as f: f.write(''.join(lines))
   
                # Process the unknown problem by users manually.
                else:
                    manually_process.append(smi)
                    print('Please manually choose applicable key words to elimate error for {}\n'.format(smi))
            
        # Show manually processed molecules.   
        if manually_process: print('\n\n!!! Manually Processed smiles. !!!\n{}\n'.format('\t'.join(manually_process)))
    
    

    """ To check the output results and provide solutions to solve incorrect output. """
    def check_method_out(self):
        
        # Remove redundant files and determine finished molecules.
        print('\nCheck calculations of {} ...\n'.format(self.method))
        redundant_files = array([glob('*.chk'), glob('job*'), glob('nohup.out'), glob('fort.7'), glob('molpro*')]).flatten()
        for file in redundant_files: remove(file)
        finished_smiles = [smi for smi in self.smiles if smi + '.out' in listdir('.')]
        
        # Check calculations of each method whether completed.
        exec('self.check_{}_out({})'.format(self.method, finished_smiles))
        self.check_method_completed()
        
