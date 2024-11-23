from collections import defaultdict
import logging
import os
import json
import typing

import click
import numpy as np
from openff.units import unit
import tqdm
# suppress stereochemistry warnings
logging.getLogger("openff").setLevel(logging.ERROR)

if typing.TYPE_CHECKING:
    from openff.toolkit import Molecule, ForceField
    from qcportal.models import OptimizationRecord, ResultRecord


def calculate_parameters(
    qc_record: "ResultRecord",
    molecule: "Molecule",
    forcefield: "ForceField",
) -> typing.Dict[str, typing.Dict[str, typing.List[unit.Quantity]]]:
    """
    Calculate the modified seminario parameters for the given input molecule
    and store them by OFF SMIRKS.
    """
    from qubekit.molecules import Ligand
    from qubekit.bonded.mod_seminario import ModSeminario

    mod_sem = ModSeminario()

    # create the qube molecule, this should be in the same order as the off_mol
    qube_mol = Ligand.from_rdkit(molecule.to_rdkit(), name="offmol")
    qube_mol.hessian = qc_record.return_result
    # calculate the modified seminario parameters and store in the molecule
    qube_mol = mod_sem.run(qube_mol)
    # label the openff molecule
    labels = forcefield.label_molecules(molecule.to_topology())[0]
    # loop over all bonds and angles and collect the results in nm/ kj/mol / radians(openMM units)
    all_parameters = {
        "bond_eq": defaultdict(list),
        "bond_k": defaultdict(list),
        "angle_eq": defaultdict(list),
        "angle_k": defaultdict(list),
    }


    for bond, parameter in labels["Bonds"].items():
        # bond is a tuple of the atom index the parameter is applied to
        qube_param = qube_mol.BondForce[bond]
        all_parameters["bond_eq"][parameter.smirks].append(qube_param.length)
        all_parameters["bond_k"][parameter.smirks].append(qube_param.k)
        
    for angle, parameter in labels["Angles"].items():
        qube_param = qube_mol.AngleForce[angle]
        all_parameters["angle_eq"][parameter.smirks].append(qube_param.angle)
        all_parameters["angle_k"][parameter.smirks].append(qube_param.k)
        
    return all_parameters


@click.command()
@click.option(
    "--initial-force-field",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the initial force field file (OFFXML).",
)
@click.option(
    "--output",
    "output_force_field",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the output force field file (OFFXML).",
)
@click.option(
    "--optimization-dataset",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the optimization dataset.",
)
@click.option(
    "--working-directory",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    required=False,
    help=(
        "The path to the working directory. "
        "Intermediate files are saved here if provided"
    ),
)
@click.option(
    "--hessian-datasets",
    type = str,required=False,multiple=True,default = [], help = ('Name of Hessian dataset(s) corresponding to single points for the molecules in the Optimization dataset')
)
@click.option(
    "--frozen-angle-file",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=False,
    help="The path to a JSON file containing frozen angle smirks (e.g. for linear angles)",
)
@click.option(
    "--verbose/--no-verbose",
    default=False,
    help="Enable verbose logging.",
)
def main(
    initial_force_field: str,
    output_force_field: str,
    optimization_dataset: str,
    working_directory: typing.Optional[str] = None,
    hessian_datasets: typing.Optional[list] = [],
    frozen_angle_file: typing.Optional[str]=None,
    verbose: bool = False
):
    from openff.toolkit import ForceField
    from openff.qcsubmit.results import (
        BasicResultCollection,
        OptimizationResultCollection
    )
    from openff.qcsubmit.results.filters import LowestEnergyFilter

    dataset = OptimizationResultCollection.parse_file(optimization_dataset)

    # filter for lowest energy results
    filtered = dataset.filter(LowestEnergyFilter())

    # filter to only keep entries with hessians calculated
    hessian_set = filtered.to_basic_result_collection(driver="hessian")


    if len(hessian_datasets) > 0:
        if verbose:
            print(filtered.n_results,' optimization records after LowestEnergyFilter')
            print(hessian_set.n_results,' hessians after converting to basic dataset')

        # Download Hessian datasets
        from qcportal import PortalClient
        client = PortalClient(address="https://api.qcarchive.molssi.org:443/")
        hess_dataset = BasicResultCollection.from_server(
            client=client,
            datasets=hessian_datasets,
            spec_name="default",
        )

        # Match ids between (filtered) opt dataset and (unfiltered) hess dataset, to avoid costly filtering step
        print('Converting Hessian-only dataset to records')
        opt_records = filtered.to_records()
        hess_records = hess_dataset.to_records()

        print('Matching extra Hessian records to Optimization records')
        record_ids_to_add = []
        n_no_match = 0
        for i,rec in enumerate(hess_records):
            any_rmsd_match = False
            for j,rec2 in enumerate(opt_records): # only loop over the filtered opt results
                if rec[1].to_smiles() == rec2[1].to_smiles(): # make sure smiles match
                    if not np.any(np.abs(rec[1].conformers[0] - rec2[1].conformers[0]).magnitude > 0.00000001): # QCA uses 1e-9, need to use 1e-8
                        # Only keep entries where the MolID isn't the same, since if it is the same it'll already be in the DS
                        if rec2[0].final_molecule_id != rec[0].molecule_id:
                            record_ids_to_add.append(rec[0].id)
                            any_rmsd_match = True
                            break
            if not any_rmsd_match:
                n_no_match += 1
        print('Number of Hessian records found on QCA:',hess_dataset.n_results) # Should be 899
        print('Number of records matched:   ',len(record_ids_to_add)) 
        print('Number of records filtered out: ',n_no_match) # Anything here would mean that that entry was filtered out of Opt, or removed for not being the lowest energy conf

        # Add new hessian entries to our dataset, if the corresponding opt molecule id is in the original opt set, e.g. if it wasn't filtered out
        hessian_set.entries[client.address] += [entry for entry in hess_dataset.entries[client.address] if entry.record_id in record_ids_to_add]

    if working_directory is not None:
        hessian_file = os.path.join(working_directory, "hessian_set.json")
        with open(hessian_file, "w") as f:
            f.write(hessian_set.json(indent=2))
        if verbose:
            print(f"Hessian set written to: {hessian_file}")
    
    if verbose:
        print(f"Found {hessian_set.n_results} hessian calculations")
        print(f"Found {hessian_set.n_molecules} hessian molecules")

    ff = ForceField(initial_force_field, allow_cosmetic_attributes=True)

    # load potential frozen angles
    frozen_angle_smirks = []
    if frozen_angle_file is not None:
        with open(frozen_angle_file, "r") as f:
            frozen_angle_smirks = json.load(f)["smirks"]

    # calculate MSM parameters for the dataset
    records_and_molecules = list(hessian_set.to_records())
    if verbose:
        records_and_molecules = tqdm.tqdm(
            records_and_molecules,
            desc="Calculating parameters",
        )

    all_parameters = {
        "bond_eq": defaultdict(list),
        "bond_k": defaultdict(list),
        "angle_eq": defaultdict(list),
        "angle_k": defaultdict(list),
    }
    errored_records_and_molecules = []
    for record, molecule in records_and_molecules:
        try:
            parameters = calculate_parameters(record, molecule, ff)
        except BaseException:
            errored_records_and_molecules.append((record, molecule))
            continue
        else:
            for key, values in parameters.items():
                for smirks, value in values.items():
                    if key == "angle_eq" and smirks in frozen_angle_smirks: # don't set frozen angles
                        continue
                    all_parameters[key][smirks].extend(value)
    
    if working_directory is not None:
        seminario_file = os.path.join(working_directory, "seminario_parameters.json")
        with open(seminario_file, "w") as file:
            json.dump(all_parameters, file, indent=2)

    if verbose:
        print(f"Found {len(errored_records_and_molecules)} errored calculations")
    if working_directory is not None:
        if len(errored_records_and_molecules):
            key = list(dataset.entries.keys())[0]
            opt_records_by_id = {
                record.record_id: record
                for record in hessian_set.entries[key]
            }
            records, _ = zip(*errored_records_and_molecules)
            errored_records = [
                opt_records_by_id[record.id]
                for record in records
            ]
            errored_dataset = BasicResultCollection(
                entries={
                    key: errored_records
                }
            )
            error_file = os.path.join(working_directory, "errored_dataset.json")
            with open(error_file, "w") as f:
                f.write(errored_dataset.json(indent=2))
            if verbose:
                print(f"Errored dataset written to: {error_file}")
    

    # now we need to update the initial FF parameters to be the
    # mean of the MSM parameters for each molecule covered by that parameter
    kj_per_mol_per_nm2 = unit.kilojoule_per_mole / unit.nanometer ** 2
    bond_handler = ff.get_parameter_handler("Bonds")
    for smirks in all_parameters["bond_eq"]:
        bond = bond_handler.parameters[smirks]

        bond_length = np.mean(all_parameters["bond_eq"][smirks]) * unit.nanometer
        bond.length = bond_length.to(unit.angstrom)

        bond_k = np.mean(all_parameters["bond_k"][smirks]) * kj_per_mol_per_nm2
        bond.k = bond_k.to(unit.kilocalorie_per_mole / (unit.angstrom ** 2))

    kj_per_mol_per_rad2 = unit.kilojoule_per_mole / (unit.radian ** 2)
    angle_handler = ff.get_parameter_handler("Angles")
    for smirks in all_parameters["angle_k"]:
        angle = angle_handler.parameters[smirks]

        # allow for empty smirks                                                            
        if all_parameters["angle_eq"][smirks]:
            angle_eq = np.mean(all_parameters["angle_eq"][smirks]) * unit.radian
            angle.angle = angle_eq.to(unit.degree)
                                                                                             
        if all_parameters["angle_k"][smirks]:
            angle_k = np.mean(all_parameters["angle_k"][smirks]) * kj_per_mol_per_rad2       
            angle.k = angle_k.to(unit.kilocalorie_per_mole / unit.radian ** 2)               
    
    ff.to_file(output_force_field)



if __name__ == "__main__":
    main()
