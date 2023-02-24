import numpy as np
from math import ceil
from rdkit import Chem
from collections import Counter
from numpy.linalg import lstsq
from rdkit.Chem import AllChem, rdMolTransforms
from _CBH import CBH


bond_mark = {1.0: '-', 2.0:'=', 3.0:'#', 1.5:'$'}


def get_training_data(file):
    smiles = []
    ori_RED = []

    with open(file) as f:
        for line in f.readlines()[1:]:
            if not line.strip(): break
            smi, RED = line.split()
            smiles.append(smi)
            ori_RED.append(RED)

    return smiles, ori_RED


def get_atom_info(mol, atom):
    p_hybrid_degree = str(atom.GetHybridization())[-1] if str(atom.GetHybridization())[-1] in ["2", "3"] else "1"
    adjacent_heavy_atoms = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() != 'H']
    adjacent_H_atoms = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() == 'H']

    atom_info = '{}{}{}{}'.format(atom.GetSymbol(), len(adjacent_heavy_atoms), len(adjacent_H_atoms), p_hybrid_degree) if atom.GetSymbol() != 'H' else 'H'

    return atom_info



def get_single_bonds(smi):
    all_bonds = []
    newsmi = smi[4:] if smi.startswith('cis-') else smi
    mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))

    for bond in mol.GetBonds():
        left_atom = bond.GetBeginAtom()
        right_atom = bond.GetEndAtom()
        bond_info = bond_mark[bond.GetBondTypeAsDouble()].join(sorted([get_atom_info(mol, left_atom), get_atom_info(mol, right_atom)]))
        
        if bond_info[0] == 'H':
            bond_info = bond_info[2:] + bond_info[1] + bond_info[0]
        all_bonds.append(bond_info)

    all_bonds = Counter(all_bonds)

    return all_bonds


def get_binary_bonds(smi):
    all_bonds = []
    newsmi = smi[4:] if smi.startswith('cis-') else smi
    mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))

    for atom in mol.GetAtoms():
        if len([True for x in atom.GetNeighbors() if x.GetSymbol() != 'H']) < 2: continue
        info = sorted([(get_atom_info(mol,x),bond_mark[mol.GetBondBetweenAtoms(x.GetIdx(),atom.GetIdx()).GetBondTypeAsDouble()]) for x in atom.GetNeighbors() if x.GetSymbol() != 'H'])
        all_atoms = [x for x in atom.GetNeighbors() if x.GetSymbol() != 'H'] + [atom]
        ring_tag = 'ring-' if all([x.IsInRing() for x in all_atoms]) else ''
        
        if len(info[:ceil(len(info) / 2)]) > 1:
            left_outside_info = '[{}]'.format('|'.join([x[0]+x[1] for x in info[:ceil(len(info) / 2)]]))
        else:
            left_outside_info = info[0][0] + info[0][1]
        if len(info[ceil(len(info) / 2):]) > 1:
            right_outside_info='[{}]'.format('|'.join([x[1] + x[0] for x in info[ceil(len(info) / 2):]]))
        else:
            right_outside_info = info[-1][-1] + info[-1][0]
        bond_info = '{}{}{}{}'.format(ring_tag, left_outside_info, get_atom_info(mol,atom), right_outside_info)
        all_bonds.append(bond_info)

    all_bonds = Counter(all_bonds)

    return all_bonds



def get_conjugated_bonds(smi):
    all_bonds = []
    newsmi = smi[4:] if smi.startswith('cis-') else smi
    mol = Chem.AddHs(Chem.MolFromSmiles(newsmi))
    AllChem.Compute2DCoords(mol)

    for bond in mol.GetBonds():
        left_atom = bond.GetBeginAtom()
        right_atom = bond.GetEndAtom()
        left_heavy_atoms = [x.GetIdx() for x in left_atom.GetNeighbors() if x.GetSymbol() != 'H' and x.GetIdx() != right_atom.GetIdx()]
        right_heavy_atoms = [x.GetIdx() for x in right_atom.GetNeighbors() if x.GetSymbol() != 'H' and x.GetIdx() != left_atom.GetIdx()]
        left_bonds = [mol.GetBondBetweenAtoms(left_atom.GetIdx(), atom_idx) for atom_idx in left_heavy_atoms]
        right_bonds = [mol.GetBondBetweenAtoms(right_atom.GetIdx(), atom_idx) for atom_idx in right_heavy_atoms]
        left_central_atom_info = get_atom_info(mol, left_atom)
        right_central_atom_info = get_atom_info(mol, right_atom)
        left_outside = sorted([(get_atom_info(mol, x.GetOtherAtom(left_atom)), bond_mark[x.GetBondTypeAsDouble()], x.GetOtherAtom(left_atom).GetIdx()) for x in left_bonds])
        right_outside = sorted([(get_atom_info(mol, x.GetOtherAtom(right_atom)), bond_mark[x.GetBondTypeAsDouble()], x.GetOtherAtom(right_atom).GetIdx()) for x in right_bonds])

        if left_heavy_atoms and right_heavy_atoms:

            info = sorted([[left_outside, left_central_atom_info], [right_outside, right_central_atom_info]])

            left_outside_info = '|'.join([x[0] + x[1] for x in info[0][0]])
            right_outside_info = '|'.join([x[1] + x[0] for x in info[1][0]])

            left_outside_info = '[{}]'.format(left_outside_info) if '|' in left_outside_info else left_outside_info
            right_outside_info = '[{}]'.format(right_outside_info) if '|' in right_outside_info else right_outside_info

            if ([x.GetIsConjugated() for x in left_bonds + right_bonds].count(True) >= 2 and bond.GetIsConjugated()) or bond.GetBondTypeAsDouble() > 1.0:
                dihedral_angle = abs(rdMolTransforms.GetDihedralDeg(mol.GetConformer(), left_outside[0][2], left_atom.GetIdx(), right_atom.GetIdx(), right_outside[0][2]))
                cis_tag = 'cis-' if dihedral_angle < 90 else ''
                if smi.startswith('cis-'):
                    cis_tag = 'cis-'
                
                bond_info = '{}{}{}{}{}{}'.format(cis_tag, left_outside_info, info[0][1], bond_mark[bond.GetBondTypeAsDouble()], info[1][1], right_outside_info)
                all_bonds.append(bond_info)

    all_bonds = Counter(all_bonds)

    return all_bonds



def get_single_mol_bonds(smi):
    all_bonds = get_single_bonds(smi) + get_binary_bonds(smi) + get_conjugated_bonds(smi)

    return all_bonds



def get_reac_delta_bonds(smi):
    reac, prod=CBH().get_CBH3(smi)
    if reac == prod:
        return Counter({})
    bond_types = []

    for k, v in (reac + prod).items():
        bonds = get_single_mol_bonds(k)
        bond_types.extend(bonds.keys())
    bond_types = sorted(Counter(bond_types))

    L = [0] * len(bond_types)
    R = [0] * len(bond_types)
    for x, y in reac.items():
        bonds = get_single_mol_bonds(x)
        
        for i, v in enumerate(bond_types):
            L[i] = L[i] + bonds.get(v,0) * y
    
    for x, y in prod.items():
        bonds = get_single_mol_bonds(x)
        
        for i, v in enumerate(bond_types):
            R[i] = R[i] + bonds.get(v, 0) * y

    delta_coefficients = np.array(R) - np.array(L)
    all_bonds = Counter(dict(zip(bond_types, delta_coefficients)))

    for x in list(all_bonds):
        if all_bonds[x] == 0:
            del all_bonds[x]

    return all_bonds



def get_all_delta_bonds(smiles):
    all_bonds = []
    all_bond_types = []

    for smi in smiles:
        bonds = get_reac_delta_bonds(smi)
        for x in bonds.keys():
            if x not in all_bond_types:
                all_bond_types.append(x)

    for smi in smiles:
        bonds = get_reac_delta_bonds(smi)
        all_bonds.append([bonds.get(x, 0) for x in all_bond_types])

    return all_bond_types, all_bonds



def get_fitted_parameters(smiles, ori_RED):
    ori_RED = np.array(ori_RED, dtype = float)
    bond_parameters = {}
    all_bonds = []

    types, bonds = get_all_delta_bonds(smiles)
    bond_types = sorted(set(types))


    for i, v in enumerate(smiles):
        bond = dict(zip(types, bonds[i]))
        all_bonds.append([bond.get(x, 0) for x in bond_types])

    xdata = np.array(all_bonds)
    ydata = ori_RED

    parameters = np.linalg.pinv(xdata).dot(ydata)
    parameters = np.round(np.array(parameters), 3)
    bond_parameters = dict(zip(bond_types, - parameters)) # The BDC parameter is the negative value of bond energy difference.

    fitted_RED = np.dot(xdata, parameters)
    residual_RED = np.array(fitted_RED) - np.array(ydata) 
    
    with open('fitted_data.txt', 'w') as f:
        print('%25s%6.2f kcal/mol'%('Largest deviation: ', np.absolute(residual_RED).max()), file=f)
        print('%25s%6.2f kcal/mol'%('Mean absolute deviation: ',np.absolute(residual_RED).mean()), file=f)
        print('%60s%24s%24s%24s\n'%('SMILES', 'original RED', 'fitted RED', 'residue RED'), file = f)
        
        file_content = ''

        for i, smi in enumerate(smiles):
            file_content += '%60s%24.2f%24.2f%24.2f\n'%(smi, ori_RED[i], fitted_RED[i], residual_RED[i])
        print(file_content, file = f)

    return bond_parameters






if __name__ == '__main__':

    #################  Fitting  #################

    print('Fitting bond parameters...')


    smiles, ori_RED= get_training_data('training_data.txt')


    bond_parameters = get_fitted_parameters(smiles, ori_RED)


    with open('BDC_parameters.txt','w') as f:
        f.write(str(bond_parameters))
