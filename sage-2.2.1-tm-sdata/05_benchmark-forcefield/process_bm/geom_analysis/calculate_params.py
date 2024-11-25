import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from openff.toolkit import Molecule, ForceField, Topology
import json
import tqdm
from yammbs import MoleculeStore
from openff.units import unit
import click
import os
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
    qm_conf = qm_mol_rd.GetConformer()
    mm_conf = mm_mol_rd.GetConformer()

    # these have the form {(atom1 idx, atom2 idx,...): FF parameter} for all bonds/angles/torsions in the molecule
    bond_dict = dict(molecule_force_list[0]['Bonds'])
    angle_dict = dict(molecule_force_list[0]['Angles'])
    proper_dict = dict(molecule_force_list[0]['ProperTorsions'])
    improper_dict = dict(molecule_force_list[0]['ImproperTorsions'])

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


    # Angles
    for idx in angle_dict.keys():
        b = angle_dict[idx]

        if filter: # To do: change this so that the match is exact, e.g. all 3 indices are present. Complicated if only one or two atoms in pattern.
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
    # To do: standardize boundary condition
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

def get_mols_from_sqlite(sqlite_store,ff_file,ff_yammbs,all_data_dicts,filter_pattern=None,conformers=False,outliers=None):
    print("Loading molecules from database ",sqlite_store)
    store = MoleculeStore(sqlite_store)
    mol_ids = store.get_molecule_ids() # each molecule has many conformers

    all_smiles = store.get_smiles() # these are the unique mapped smiles
    qca_ids = [store.get_qcarchive_ids_by_molecule_id(i) for i in mol_ids] # grouped by conformer

    # Need to insert some checks here. E.g. what if no force fields available, what if no database available, etc
    if ff_yammbs == None:
        all_ffs = store.get_force_fields()
        print('No YAMMBS force field provided. Available force fields in database: ',all_ffs)
        print('Using geometries from', all_ffs[0])
        ff_yammbs = all_ffs[0]

    print('Grouping parameters based on ',ff_file)
    ff = ForceField(ff_file,allow_cosmetic_attributes=True)

    smiles_no_errs = []
    qca_ids_no_errs = []
    qm_confs = []
    mm_confs = []
    errors = []
    bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,qm_data_dict,mm_data_dict = all_data_dicts

    mol_ids = mol_ids
    total = 0
    if not conformers: # If the user doesn't want to use all conformers, we just pick the first conformer.
        for idx in tqdm.tqdm(range(0,len(mol_ids)),desc='Calculating geometric parameters'):
            i = mol_ids[idx]
            total += 1
            try:
                qca_id = qca_ids[idx][0]
                if qca_id not in outliers:
                    mapped_smiles = store.get_smiles_by_molecule_id(i)

                    mm_conf = store.get_mm_conformers_by_molecule_id(i,force_field = ff_yammbs)[0]
                    qm_conf = store.get_qm_conformers_by_molecule_id(i)[0]

                    mm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
                    mm_mol.add_conformer(unit.Quantity(mm_conf,'angstrom'))
                    mm_mol_rdkit = mm_mol.to_rdkit()

                    qm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
                    qm_mol.add_conformer(unit.Quantity(qm_conf,'angstrom'))
                    qm_mol_rdkit = qm_mol.to_rdkit()

                    topol = Topology.from_molecules([qm_mol])
                    mol_force_list = ff.label_molecules(topol)

                    add_mol_to_dict(qm_mol_rdkit,mm_mol_rdkit, mapped_smiles,qca_id,mol_force_list,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter=filter_pattern)

            except IndexError: # If no conformers stored for this molecule
                errors.append(i)

    else: # If the user requested all conformers be analyzed. To do: parallelize
        for idx in tqdm.tqdm(range(0,len(mol_ids)),desc='Calculating geometric parameters'):
            for j,qca_id in enumerate(qca_ids[idx]):
                total += 1
                try:
                    if qca_id not in outliers:
                        mm_conf = store.get_mm_conformer_by_qcarchive_id(qca_id,force_field = ff_yammbs)
                        qm_conf = store.get_qm_conformer_by_qcarchive_id(qca_id)

                        mapped_smiles = store.get_smiles_by_molecule_id(mol_ids[idx])

                        mm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
                        mm_mol.add_conformer(unit.Quantity(mm_conf,'angstrom'))
                        mm_mol_rdkit = mm_mol.to_rdkit()

                        qm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
                        qm_mol.add_conformer(unit.Quantity(qm_conf,'angstrom'))
                        qm_mol_rdkit = qm_mol.to_rdkit()

                        topol = Topology.from_molecules([qm_mol])
                        mol_force_list = ff.label_molecules(topol)

                        add_mol_to_dict(qm_mol_rdkit,mm_mol_rdkit, mapped_smiles,qca_id,mol_force_list,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter=filter_pattern)

                except IndexError:
                    errors.append(mol_ids[idx])

    print('Number of errors: {} out of {}'.format(len(errors),total))


def get_mols_from_files(mm_dir0,ff_file,qm_dir0,all_data_dicts,filter_pattern=None,conformers=False,dir = '.',outliers=None):
    # This function is written using my conformer/molecule naming convention.
    # Make sure to modify it if you use a different convention.
    # Structures should be named mol-i-conf-j.sdf for this script to work

    qm_dir = dir + qm_dir0
    mm_dir = dir + mm_dir0
    print('Loading QM structures from ',qm_dir)
    print('Loading MM structures from ',mm_dir)

    # Obtain list of structures from directory
    if conformers:
        all_mols = sorted(os.listdir(qm_dir))
    else:
        all_mols = sorted(np.unique(['-'.join(filename.split('-')[:-2]) for filename in os.listdir(qm_dir)]))

    print('Grouping parameters based on provided force field ',ff_file)
    ff = ForceField(ff_file,allow_cosmetic_attributes=True)

    bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,qm_data_dict,mm_data_dict = all_data_dicts
    for mol in tqdm.tqdm(all_mols,desc='Calculating geometric parameters'): # Note this progress bar will index number of conformers, while SQLite indexes number of molecules

        if not conformers:
            mol += '-conf-00.sdf'

        qm_file = qm_dir +'/'+ mol
        mm_file = mm_dir +'/'+ mol

        qm_mol = Molecule(qm_file,allow_undefined_stereo=True)
        qm_mol_rdkit = Chem.SDMolSupplier(qm_file,removeHs=False)[0]

        try:
            qca_id = qm_mol_rdkit.GetPropsAsDict()['Record QCArchive'] # This would have to be manually modified if you use a different naming convention
        except KeyError:
            qca_id = 0
            if len(outliers)>0:
                print('WARNING: QCArchive ID not present in SDF file. All files will be included.')
        if qca_id not in outliers:
            try:
                mapped_smiles = qm_mol_rdkit.GetPropsAsDict()['Mapped SMILES']
            except KeyError:
                mapped_smiles = qm_mol.to_smiles(mapped=True)

            topol = Topology.from_molecules([qm_mol])
            mol_force_list = ff.label_molecules(topol) # dictionary of forces for the molecule, keys are type of force ('Bonds', 'Angles', etc)

            mm_mol = Molecule(mm_file,allow_undefined_stereo=True)
            mm_mol_rdkit = Chem.SDMolSupplier(mm_file,removeHs=False)[0]


            try:
                add_mol_to_dict(qm_mol_rdkit, mm_mol_rdkit, mapped_smiles,qca_id,mol_force_list,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter=filter_pattern)

            except OSError: # Usually file not found
                print('Error with molecule ',mol)
                pass


@click.command()
@click.option('--db',default=None,help='SQlite database to read from (if using)')
@click.option('--mm_dir',default=None,help='Directory with MM optimization SDF files (if not reading from database)')
@click.option('--qm_dir',default=None,help='Directory with QM optimization SDF files (if not reading from database)')
@click.option('--ff_yammbs',default=None,help='Force field used to calculate geometries with yammbs (if reading from database). If not provided, it will use the first stored force field in the database.')
@click.option('--ff_file',default='openff-2.1.0.offxml',help='Force field to group parameters by.')
@click.option('--conformers',default=False,help='Whether to use all conformers. If False, just use the first conformer for each molecule. Note that setting this to True can lead to very large files and long run times.')
@click.option('--dir',default = './',help='Directory where QM and MM directories are located')
@click.option('--label',default=None,help='Label to save data files with. Data will be stored in a directory with this name.')
@click.option('--filter',default=None,help='SMARTS pattern to filter for. Must have at least one tagged atom.')
@click.option('--problem_file',default=[],multiple=True,help='File(s) listing QCAIDs of conformers to exclude from analysis.')
def main(db,mm_dir,qm_dir,ff_yammbs,ff_file,conformers,dir,label,filter,problem_file):
    # Collecting data for all molecules
    bond_data_dict = {}
    angle_data_dict = {}
    proper_data_dict = {}
    improper_data_dict = {}
    qm_data_dict = {}
    mm_data_dict = {}

    all_data_dicts = [bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,qm_data_dict,mm_data_dict]

    # Collect QCAIDs of any problematic entries
    problem_ids = [] # will be empty if there are no outliers specified
    if len(problem_file) > 0:

        for filter_id_file in problem_file:
            filter_id = np.loadtxt(filter_id_file)
            problem_ids.extend(filter_id)

    # If no yammbs database is provided, attempt to read from qm and mm directories
    if db == None:
        get_mols_from_files(mm_dir,ff_file,qm_dir,all_data_dicts,filter_pattern=filter,conformers=conformers,dir=dir,outliers=problem_ids)

    # If yammbs database is provided, calculate parameters based on stored geometries
    else:
        get_mols_from_sqlite(db,ff_file,ff_yammbs,all_data_dicts,filter_pattern=filter,conformers=conformers,outliers=problem_ids)

    if label == None:
        label = '.'.join(ff_file.split('/')[-1].split('.')[:-1])
        print('WARNING: No label provided; saving data using force field name ',label)

    # Make data directory
    try:
        os.mkdir(label)
    except FileExistsError:
        pass

    with open('{}/bonds_qmv{}.json'.format(label,label),'w') as jsonfile:
        json.dump(bond_data_dict,jsonfile,indent=4)

    with open('{}/angles_qmv{}.json'.format(label,label),'w') as jsonfile:
        json.dump(angle_data_dict,jsonfile,indent=4)

    with open('{}/propers_qmv{}.json'.format(label,label),'w') as jsonfile:
        json.dump(proper_data_dict,jsonfile,indent=4)

    with open('{}/impropers_qmv{}.json'.format(label,label),'w') as jsonfile:
        json.dump(improper_data_dict,jsonfile,indent=4)


if __name__ == '__main__':
    main()
