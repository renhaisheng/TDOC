import numpy as np
from  rdkit import Chem
from numpy.linalg import lstsq
from collections import Counter


bond_mark = {1.0: '-', 2.0:'=', 3.0:'#', 1.5:'$'}

multiplicities = {'S':'1', 'D':'2', 'T':'3', 'Q':'4', 'P':'5'}


def get_training_data(file):
    smiles = []
    all_ben = []
    all_raw = []

    with open(file) as f:
        for line in f.readlines()[1:]:
            if not line.strip(): break
            smi, ben, raw = line.split()
            smiles.append(smi)
            all_ben.append(ben)
            all_raw.append(raw)

    return smiles, all_ben, all_raw



def get_BAC_atom_info(mol, atom):
    p_hybrid_degree = str(atom.GetHybridization())[-1] if str(atom.GetHybridization())[-1] in ["2", "3"] else "1"
    electrons=atom.GetNumRadicalElectrons()

    adjacent_heavy_atoms = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() != 'H']
    adjacent_H_atoms = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() == 'H']

    # BAC method.
    #atom_info = '{}'.format(atom.GetSymbol()) if atom.GetSymbol() != 'H' else 'H'
    
    # BAC-P method. 
    #atom_info = '{}{}'.format(atom.GetSymbol(), p_hybrid_degree) if atom.GetSymbol() != 'H' else 'H'
    
    # BAC-EP method. 
    atom_info = '{}{}{}'.format(atom.GetSymbol(), electrons, p_hybrid_degree) if atom.GetSymbol() != 'H' else 'H'
    
    # BAC-HHP method.
    #atom_info = '{}{}{}{}'.format(atom.GetSymbol(), len(adjacent_heavy_atoms), len(adjacent_H_atoms), p_hybrid_degree) if atom.GetSymbol() != 'H' else 'H'  #
    
    return atom_info



def get_BAC_bonds(smi):
    all_bonds = []
    newsmi = ''.join(smi[:-2]) if '-' in smi and smi[-1] in multiplicities else smi
    mol=Chem.AddHs(Chem.MolFromSmiles(newsmi))

    for bond in mol.GetBonds():
        left_atom = bond.GetBeginAtom()
        right_atom = bond.GetEndAtom()

        bond_info = bond_mark[bond.GetBondTypeAsDouble()].join(sorted([get_BAC_atom_info(mol,left_atom), get_BAC_atom_info(mol, right_atom)]))
        
        if bond_info[0] == 'H':
            bond_info = bond_info[2:] + bond_info[1] + bond_info[0]
        all_bonds.append(bond_info)

    all_bonds = Counter(all_bonds)

    return all_bonds



def get_all_mol_bonds(smiles):
    all_bond_types=[]
    all_bonds=[]

    for smi in smiles:
        bonds = get_BAC_bonds(smi)
        for x in bonds.keys():
            if x not in all_bond_types:
                all_bond_types.append(x)

    for smi in smiles:
        bonds = get_BAC_bonds(smi)
        all_bonds.append([bonds.get(x, 0) for x in all_bond_types])
    
    return all_bond_types,all_bonds



def get_fitted_parameters(smiles, ben, raw):
    ben, raw = np.array(ben, dtype = float), np.array(raw, dtype = float)
    bond_parameters = {}
    all_bonds = []
    bond_types, all_bonds = get_all_mol_bonds(smiles)

    xdata = np.array(all_bonds)
    ydata = ben - raw
    
    parameters = np.linalg.pinv(xdata).dot(ydata)
    parameters = np.round(np.array(parameters), 3)
    bond_parameters = dict(zip(bond_types, parameters))

    fitted_data = raw + np.dot(xdata, parameters)
    residual_error = fitted_data - ben

    with open('fitted_data.txt','w') as f:
        print('%25s%6.2f kcal/mol'%('Largest deviation: ', np.absolute(residual_error).max()), file = f)
        print('%25s%6.2f kcal/mol'%('Mean absolute deviation: ', np.absolute(residual_error).mean()), file = f)
        print('%60s%24s%24s%24s%24s'%('SMILES', 'benchmark (kcal/mol)', 'raw data (kcal/mol)', 'fitted data (kcal/mol)', 'residual error'), file = f)
        
        file_content = ''

        for i, smi in enumerate(smiles):
            file_content += '%60s%24.2f%24.2f%24.2f%24.2f\n'%(smi, ben[i], raw[i], fitted_data[i], residual_error[i])
        print(file_content, file = f)

    return bond_parameters






if __name__ == '__main__':

    #################  Fitting  #################

    print('Fitting bond parameters...')


    smiles, all_ben, all_raw = get_training_data('training_data.txt')


    bond_parameters = get_fitted_parameters(smiles, all_ben, all_raw)


    with open('BAC_parameters.txt','w') as f:
        f.write(str(bond_parameters))
