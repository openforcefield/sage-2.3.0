import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from openff.toolkit import Molecule, ForceField, Topology
import json
import tqdm
import sys
from yammbs import MoleculeStore
from openff.units import unit
import click
import logging
logging.getLogger("openff").setLevel(logging.ERROR)

# def add_mol_to_dict(qm_file,mm_file,ff,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter=False):
def add_mol_to_dict(qm_mol_rd,mm_mol_rd,smiles,recordid,molecule_force_list,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter=None):
    '''
    Function to enumerate FF parameters assigned to a molecule (e.g. bond length, angle, dihedral angle),
    calculate the actual value from both QM and MM, and compile a set of parameter dictionaries for plotting.

    Inputs:
        qm_file: SDF file with the optimized QM geometry
        mm_file: SDF file with the optimized MM geometry
        ff: ForceField object to assign parameters with
        bond_data_dict: dictionary where bond lengths are enumerated. Can be empty
        angle_data_dict: dictionary where angles are enumerated. Can be empty
        proper_data_dict: dictionary where dihedral angles corresponding to
                          proper torsion parameters are enumerated. Can be empty
        improper_data_dict: dictionary where dihedral angles corresponding to
                          improper torsion parameters are enumerated. Can be empty

    Returns:
        Nothing is returned, but data dictionaries are modified in place.
    '''
    # print('Adding mol to dict')
    qm_conf = qm_mol_rd.GetConformer()
    mm_conf = mm_mol_rd.GetConformer()
    # print(mm_conf)

    # these have the form {(atom1 idx, atom2 idx,...): FF parameter} for all bonds/angles/torsions in the molecule
    bond_dict = dict(molecule_force_list[0]['Bonds'])
    angle_dict = dict(molecule_force_list[0]['Angles'])
    proper_dict = dict(molecule_force_list[0]['ProperTorsions'])
    improper_dict = dict(molecule_force_list[0]['ImproperTorsions'])
    # print(angle_dict)
    # print('Loaded dicts')

    # qm_data_dict[recordid] = {'Bonds':{},'Angles':{},'ProperTorsions':{},'ImproperTorsions':{}} # record ID: {'Bonds': (idx):qm_bl} --> can be read in to avoid recalc
    # mm_data_dict[recordid] = {'Bonds':{},'Angles':{},'ProperTorsions':{},'ImproperTorsions':{}} # same but for MM

    if filter:
        qm_mol = Molecule.from_rdkit(qm_mol_rd,allow_undefined_stereo=True)
        matches = qm_mol.chemical_environment_matches(filter)

    #Bonds
    for idx in bond_dict:
        b = bond_dict[idx]

        if filter:
            filter_condition = (idx[0],) in matches or (idx[1],) in matches

        if (filter and filter_condition) or not filter:
            # This part could potentially be eliminated by reading in QM geom data
            qm_bl = Chem.rdMolTransforms.GetBondLength(qm_conf,idx[0],idx[1])
            # qm_data_dict[recordid]['Bonds'][idx] = qm_bl

            mm_bl = Chem.rdMolTransforms.GetBondLength(mm_conf,idx[0],idx[1]) # A
            # mm_data_dict[recordid]['Bonds'][idx] = mm_bl

            if b.smirks in bond_data_dict:
                bond_data_dict[b.smirks]['molecules'].append(smiles)
                bond_data_dict[b.smirks]['envs'].append(idx)
                bond_data_dict[b.smirks]['qm_values'].append(qm_bl)
                bond_data_dict[b.smirks]['mm_values'].append(mm_bl)
            else:
                bond_data_dict[b.smirks] = {'ident':b.id,
                                            'sage_value': b.length.magnitude,
                                            'qm_values': [qm_bl],
                                            'mm_values': [mm_bl],
                                            'molecules': [smiles],
                                            'envs': [idx]}

    # print(len(angle_dict.keys()))
    # Angles
    for idx in angle_dict.keys():
        b = angle_dict[idx]

        if filter:
            filter_condition = ((idx[0],) in matches or (idx[1],) in matches) or (idx[2],) in matches

        if (filter and filter_condition) or not filter:

            # This part could potentially be eliminated by reading in QM geom data
            qm_ang = Chem.rdMolTransforms.GetAngleDeg(qm_conf,idx[0],idx[1],idx[2])
            # qm_data_dict[recordid]['Angles'][idx] = qm_ang
            mm_ang = Chem.rdMolTransforms.GetAngleDeg(mm_conf,idx[0],idx[1],idx[2])
            # mm_data_dict[recordid]['Angles'][idx] = mm_ang

            if b.smirks in angle_data_dict:
                angle_data_dict[b.smirks]['molecules'].append(smiles)
                angle_data_dict[b.smirks]['envs'].append(idx)
                angle_data_dict[b.smirks]['qm_values'].append(qm_ang)
                angle_data_dict[b.smirks]['mm_values'].append(mm_ang)
            else:
                angle_data_dict[b.smirks] = {'ident':b.id,
                                            'sage_value': b.angle.magnitude,
                                            'qm_values': [qm_ang],
                                            'mm_values': [mm_ang],
                                            'molecules': [smiles],
                                            'envs': [idx]}

    # Proper Torsions
    for idx in proper_dict:
        b = proper_dict[idx]

        if filter:
            filter_condition = ((idx[0],) in matches or (idx[1],) in matches) or ((idx[2],) in matches or (idx[3],) in matches)

        if (filter and filter_condition) or not filter:

            # This part could potentially be eliminated by reading in QM geom data
            qm_tor = Chem.rdMolTransforms.GetDihedralDeg(qm_conf,idx[0],idx[1],idx[2],idx[3])
            # qm_data_dict[recordid]['ProperTorsions'][idx] = qm_tor
            mm_tor = Chem.rdMolTransforms.GetDihedralDeg(mm_conf,idx[0],idx[1],idx[2],idx[3])
            # mm_data_dict[recordid]['ProperTorsions'][idx] = mm_tor

            if b.smirks in proper_data_dict:
                proper_data_dict[b.smirks]['molecules'].append(smiles)
                proper_data_dict[b.smirks]['envs'].append(idx)
                proper_data_dict[b.smirks]['qm_values'].append(qm_tor)
                proper_data_dict[b.smirks]['mm_values'].append(mm_tor)
            else:
                proper_data_dict[b.smirks] = {'ident':b.id,
                                            # 'sage_value': b.angle.magnitude,
                                            'qm_values': [qm_tor],
                                            'mm_values': [mm_tor],
                                            'molecules': [smiles],
                                            'envs': [idx]}

    # Improper Torsions
    for idx in improper_dict:
        b = improper_dict[idx]

        if filter:
            filter_condition = ((idx[0],) in matches or (idx[1],) in matches) or ((idx[2],) in matches or (idx[3],) in matches)

        if (filter and filter_condition) or not filter:

            # This part could potentially be eliminated by reading in QM geom data
            qm_tor = Chem.rdMolTransforms.GetDihedralDeg(qm_conf,idx[0],idx[1],idx[2],idx[3])
            # qm_data_dict[recordid]['ProperTorsions'][idx] = qm_tor
            mm_tor = Chem.rdMolTransforms.GetDihedralDeg(mm_conf,idx[0],idx[1],idx[2],idx[3])
            # mm_data_dict[recordid]['ProperTorsions'][idx] = mm_tor

            if b.smirks in improper_data_dict:
                improper_data_dict[b.smirks]['molecules'].append(smiles)
                improper_data_dict[b.smirks]['envs'].append(idx)
                improper_data_dict[b.smirks]['qm_values'].append(qm_tor)
                improper_data_dict[b.smirks]['mm_values'].append(mm_tor)
            else:
                improper_data_dict[b.smirks] = {'ident':b.id,
                                            # 'sage_value': b.angle.magnitude,
                                            'qm_values': [qm_tor],
                                            'mm_values': [mm_tor],
                                            'molecules': [smiles],
                                            'envs': [idx]}

    return

def get_mols_from_sqlite(sqlite_store,ff_file,all_data_dicts,filter_pattern=None,conformers=False,outliers=None):
    print("Loading molecules from database")
    # print(ff)
    store = MoleculeStore(sqlite_store)
    mol_ids = store.get_molecule_ids() # each molecule has many conformers
    # print(mol_ids)
    all_smiles = store.get_smiles() # these are the unique mapped smiles
    qca_ids = [store.get_qcarchive_ids_by_molecule_id(i) for i in mol_ids] # grouped by conformer


    print(store.get_force_fields())
    # mm_confs = []
    print('Loading FF')
    ff = ForceField(ff_file,allow_cosmetic_attributes=True)
    ff_file_name = store.get_force_fields()[0].split('/')[-1]
    smiles_no_errs = []
    qca_ids_no_errs = []
    qm_confs = []
    mm_confs = []
    errors = []
    bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,qm_data_dict,mm_data_dict = all_data_dicts

    mol_ids = mol_ids#[:10]
    total = 0
    if not conformers:
        # print(mol_ids)
        for idx in tqdm.tqdm(range(0,len(mol_ids)),desc='Calculating geometric parameters'):
            i = mol_ids[idx]
            # print(i)
            total += 1
            try:
                qca_id = qca_ids[idx][0]
                if qca_id not in outliers:
                    # print(qca_ids[idx])
                    # print(qca_id)
                    # mapped_smiles = all_smiles[idx]
                    mapped_smiles = store.get_smiles_by_molecule_id(i)
                    # print(mapped_smiles)

                    mm_conf = store.get_mm_conformers_by_molecule_id(i,force_field = ff_file_name)[0]
                    # mm_conf = store.get_mm_conformer_by_qcarchive_id(qca_id,force_field = ff_file_name)#[0]
                    # print(mm_conf)
                    qm_conf = store.get_qm_conformers_by_molecule_id(i)[0]
                    #print(qm_conf)
                    # qm_conf = store.get_qm_conformer_by_qcarchive_id(qca_id)
                    # print(qm_conf)

                    mm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
                    mm_mol.add_conformer(unit.Quantity(mm_conf,'angstrom'))#*unit.angstrom)
                    mm_mol_rdkit = mm_mol.to_rdkit()
                    # print(mm_mol_rdkit)

                    qm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
                    qm_mol.add_conformer(unit.Quantity(qm_conf,'angstrom'))#*unit.angstrom)
                    qm_mol_rdkit = qm_mol.to_rdkit()
                    # print(qm_mol_rdkit)

                    topol = Topology.from_molecules([qm_mol])
                    mol_force_list = ff.label_molecules(topol)
                    # print(mol_force_list)

                    add_mol_to_dict(qm_mol_rdkit,mm_mol_rdkit, mapped_smiles,qca_id,mol_force_list,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter=filter_pattern)

            except IndexError:
                # print('Error with molecule ID {}'.format(i))
                errors.append(i)
                # break

    else:
        for idx in tqdm.tqdm(range(0,len(mol_ids)),desc='Calculating geometric parameters'):
        # for i,mol in enuemrate(qca_ids): # this is grouped by _molecule_
            # mm_confs_i = store.get_mm_conformers_by_molecule_id(mol_ids[idx],force_field = ff_file_name)
            # qm_confs_i = store.get_qm_conformers_by_molecule_id(mol_ids[idx])
            for j,qca_id in enumerate(qca_ids[idx]):
                total += 1
                try:
                    if qca_id not in outliers:
                        # mm_conf = mm_confs_i[j] # flatten these
                        # qm_conf = qm_confs_i[j]
                        mm_conf = store.get_mm_conformer_by_qcarchive_id(qca_id,force_field = ff_file_name)
                        qm_conf = store.get_qm_conformer_by_qcarchive_id(qca_id)
                        # mapped_smiles = all_smiles[idx] # pull by the index of the _molecule_ and repeat nconf number of times
                        mapped_smiles = store.get_smiles_by_molecule_id(mol_ids[idx])
                        # qca_id = end(qca_ids[i][j]) # flatten qca_ids

                        mm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
                        mm_mol.add_conformer(unit.Quantity(mm_conf,unit.angstrom))
                        mm_mol_rdkit = mm_mol.to_rdkit()

                        qm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
                        qm_mol.add_conformer(unit.Quantity(qm_conf,unit.angstrom))
                        qm_mol_rdkit = qm_mol.to_rdkit()

                        topol = Topology.from_molecules([qm_mol])
                        mol_force_list = ff.label_molecules(topol)

                        add_mol_to_dict(qm_mol_rdkit,mm_mol_rdkit, mapped_smiles,qca_id,mol_force_list,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter=filter_pattern)

                except IndexError:
                    # print('Error with molecule ID {} conformer {}'.format(i,j))
                    errors.append(mol_ids[idx])

    print('Number of errors: {} out of {}'.format(len(errors),total))


def get_mols_from_files(mm_dir0,ff_file,qm_dir0,all_data_dicts,filter_pattern=None,conformers=False,dir = '/Users/lexiemcisaac/Documents/OpenFF/conformer_energy_ordering/swope_scripts/benchmarking/'):

    if conformers:
        compound_list = dir + 'all_molecules.txt'
    else:
        compound_list = dir + 'compound.list'
    qm_dir = dir + qm_dir0
    mm_dir = dir + mm_dir0
    ff = ForceField(ff_file,allow_cosmetic_attributes=True)

    all_mols = np.loadtxt(compound_list,dtype='str')#[:10]

    bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,qm_data_dict,mm_data_dict = all_data_dicts
    for mol in tqdm.tqdm(all_mols,desc='Calculating geometric parameters'):
        if not conformers:
            mol += '-00.sdf'
        qm_file = qm_dir +'/'+ mol
        mm_file = mm_dir +'/'+ mol

        qm_mol = Molecule(qm_file,allow_undefined_stereo=True)
        qm_mol_rdkit = Chem.SDMolSupplier(qm_file,removeHs=False)[0]


        mm_mol = Molecule(mm_file,allow_undefined_stereo=True)
        mm_mol_rdkit = Chem.SDMolSupplier(mm_file,removeHs=False)[0]


        qca_id = qm_mol_rdkit.GetPropsAsDict()['Record QCArchive']
        mapped_smiles = qm_mol.to_smiles(mapped=True)

        topol = Topology.from_molecules([qm_mol])
        mol_force_list = ff.label_molecules(topol) # dictionary of forces for the molecule, keys are type of force ('Bonds', 'Angles', etc)
        try:
            add_mol_to_dict(qm_mol_rdkit, mm_mol_rdkit, mapped_smiles,qca_id,mol_force_list,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter_pattern=filter_pattern)
            # add_mol_to_dict(qm_file,mm_file,sage,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict)#,filter='[r4:1]')
        except OSError:
            print('Error with molecule ',mol)
            pass


@click.command()
@click.option('--db',default=None,help='SQlite database to read from')
@click.option('--mm_dir',default=None,help='directory with MM optimization SDF files')
@click.option('--qm_dir',default='b3lyp-d3bj_dzvp',help='directory with QM optimization SDF files')
@click.option('--ff_file',default='openff-2.1.0.offxml',help='Force field to group parameters by')
@click.option('--conformers',default=False,help='Whether to use all conformers. If False, just use the first conformer for each molecule')
@click.option('--dir',default = '/Users/lexiemcisaac/Documents/OpenFF/conformer_energy_ordering/swope_scripts/benchmarking/',help='Directory where QM and MM directories are located')
@click.option('--label',default=None,help='Label to save data with')
@click.option('--filter',default=None,help='SMARTS pattern to filter for. Must have at least one tagged atom')
@click.option('--problem_file',default=[],multiple=True)
def main(db,mm_dir,qm_dir,ff_file,conformers,dir,label,filter,problem_file):
    # Collecting data for all molecules
    bond_data_dict = {}
    angle_data_dict = {}
    proper_data_dict = {}
    improper_data_dict = {}
    qm_data_dict = {}
    mm_data_dict = {}

    all_data_dicts = [bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,qm_data_dict,mm_data_dict]

    problem_ids = [] # will be empty if there are no outliers specified
    if len(problem_file) > 0:

        for filter_id_file in problem_file:
            filter_id = np.loadtxt(filter_id_file)
            problem_ids.extend(filter_id)

    if db == None:
        get_mols_from_files(mm_dir,ff_file,qm_dir,all_data_dicts,filter_pattern=filter,conformers=conformers,dir=dir)
    else:
        get_mols_from_sqlite(db,ff_file,all_data_dicts,filter_pattern=filter,conformers=conformers,outliers=problem_ids)

    if label == None:
        label = '.'.join(ff_file.split('/')[-1].split('.')[:-1])
    # print(ff_name)

    with open('bonds_qmv{}.json'.format(label),'w') as jsonfile:
        json.dump(bond_data_dict,jsonfile,indent=4)

    with open('angles_qmv{}.json'.format(label),'w') as jsonfile:
        json.dump(angle_data_dict,jsonfile,indent=4)

    with open('propers_qmv{}.json'.format(label),'w') as jsonfile:
        json.dump(proper_data_dict,jsonfile,indent=4)

    with open('impropers_qmv{}.json'.format(label),'w') as jsonfile:
        json.dump(improper_data_dict,jsonfile,indent=4)


if __name__ == '__main__':
    main()
