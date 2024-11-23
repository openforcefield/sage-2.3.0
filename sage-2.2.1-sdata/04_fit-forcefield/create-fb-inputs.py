import os
import json
import pathlib
import typing

import click


def filter_for_smarts_or_smiles(
    entry,
    smarts_to_exclude: typing.Optional[typing.List[str]] = None,
    inchi_keys_to_exclude: typing.Optional[typing.List[str]] = None,
):
    """
    Normal filtering with the to_records() call is incredibly slow;
    copy out the actual filtering function to speed it up
    """
    from openff.toolkit import Molecule

    mol = Molecule.from_mapped_smiles(entry.cmiles, allow_undefined_stereo=True)

    if smarts_to_exclude:
        for smarts in smarts_to_exclude:
            if mol.chemical_environment_matches(smarts):
                return False
    
    if inchi_keys_to_exclude:
        if mol.to_inchikey(fixed_hydrogens=True) in inchi_keys_to_exclude:
            return False
    return True
        

def filter_dataset(
    dataset,
    smarts_to_exclude: typing.Optional[typing.List[str]] = None,
    smiles_to_exclude: typing.Optional[typing.List[str]] = None,
):
    from openff.toolkit import Molecule

    inchi_keys_to_exclude = []
    if smiles_to_exclude:
        inchi_keys_to_exclude = [
            Molecule.from_smiles(smiles).to_inchikey(fixed_hydrogens=True)
            for smiles in smiles_to_exclude
        ]
    
    key = list(dataset.entries)[0]
    original_dataset = dataset.entries[key]
    filtered = [
        entry
        for entry in original_dataset
        if filter_for_smarts_or_smiles(
            entry,
            smarts_to_exclude=smarts_to_exclude,
            inchi_keys_to_exclude=inchi_keys_to_exclude,
        )
    ]
    dataset.entries[key] = filtered



def load_training_data(
    optimization_dataset: str,
    torsion_dataset: str,
    smarts_to_exclude: typing.Optional[str] = None,
    smiles_to_exclude: typing.Optional[str] = None,
    verbose: bool = False
):
    from openff.qcsubmit.results import (
        OptimizationResultCollection,
        TorsionDriveResultCollection,
    )
    from openff.qcsubmit.results.filters import (
        SMARTSFilter,
        SMILESFilter,
    )

    if smarts_to_exclude is not None:
        exclude_smarts = pathlib.Path(smarts_to_exclude).read_text().splitlines()
    else:
        exclude_smarts = []

    if smiles_to_exclude is not None:
        exclude_smiles = pathlib.Path(smiles_to_exclude).read_text().splitlines()
    else:
        exclude_smiles = []

    torsion_training_set = TorsionDriveResultCollection.parse_file(torsion_dataset)
    if verbose:
        print(f"Loaded torsion training set with {torsion_training_set.n_results} entries.")

    # torsion_training_set = torsion_training_set.filter(
    #     SMARTSFilter(smarts_to_exclude=exclude_smarts),
    #     SMILESFilter(smiles_to_exclude=exclude_smiles),
    # )
    filter_dataset(
        torsion_training_set,
        smarts_to_exclude=exclude_smarts,
        smiles_to_exclude=exclude_smiles,
    )

    if verbose:
        print(f"Filtered torsion training set to {torsion_training_set.n_results} entries.")

    optimization_training_set = OptimizationResultCollection.parse_file(optimization_dataset)
    if verbose:
        print(f"Loaded optimization training set with {optimization_training_set.n_results} entries.")
    # optimization_training_set = optimization_training_set.filter(
    #     SMARTSFilter(smarts_to_exclude=exclude_smarts),
    #     SMILESFilter(smiles_to_exclude=exclude_smiles),
    # )
    filter_dataset(
        optimization_training_set,
        smarts_to_exclude=exclude_smarts,
        smiles_to_exclude=exclude_smiles,
    )
    if verbose:
        print(f"Filtered optimization training set to {optimization_training_set.n_results} entries.")


    return torsion_training_set, optimization_training_set


@click.command()
@click.option(
    "--tag",
    type=str,
    default="fb-fit",
    help="The tag to use for the fitting run.",
)
@click.option(
    "--optimization-dataset",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the optimization dataset to use. (JSON)",
)
@click.option(
    "--torsion-dataset",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the torsion dataset to use. (JSON)",
)
@click.option(
    "--forcefield",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the force field to use. (offxml)",
)
@click.option(
    "--valence-to-optimize",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the valence parameters to optimize (JSON).",
)
@click.option(
    "--torsions-to-optimize",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the torsions to optimize (JSON).",
)
@click.option(
    "--output-directory",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    required=True,
    help="The directory to write the ForceBalance inputs to.",
)
@click.option(
    "--frozen-angle-file",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=False,
    help="The path to a JSON file containing frozen angle smirks (e.g. for linear angles)",
)
@click.option(
    "--smarts-to-exclude",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default=None,
    help=(
        "The path to a file containing a list of SMARTS patterns "
        "to exclude from the training set. "
        "The patterns should be separated by new lines."
    ),
)
@click.option(
    "--smiles-to-exclude",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default=None,
    help=(
        "The path to a file containing a list of SMILES patterns "
        "to exclude from the training set. "
        "The patterns should be separated by new lines."
    ),
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Whether to print verbose logging messages.",
)
@click.option(
    "--max-iterations",
    type=int,
    default=50,
    show_default=True,
    help="The maximum number of iterations to run the fitting for.",
)
@click.option(
    "--port",
    type=int,
    default=55387,
    show_default=True,
    help="The port to run the server on.",
)
def generate(
    tag: str,
    optimization_dataset: str,
    torsion_dataset: str,
    forcefield: str,
    valence_to_optimize: str,
    torsions_to_optimize: str,
    output_directory: str,
    frozen_angle_file: typing.Optional[str] = None,
    smarts_to_exclude: typing.Optional[str] = None,
    smiles_to_exclude: typing.Optional[str] = None,
    verbose: bool = False,
    max_iterations: int = 50,
    port: int = 55387,

):
    from openff.toolkit import ForceField
    from openff.bespokefit.optimizers.forcebalance import ForceBalanceInputFactory
    from openff.bespokefit.schema.fitting import OptimizationSchema, OptimizationStageSchema
    from openff.bespokefit.schema.optimizers import ForceBalanceSchema
    from openff.bespokefit.schema.targets import (
        OptGeoTargetSchema,
        TorsionProfileTargetSchema,
    )
    from openff.bespokefit.schema.smirnoff import (
        AngleHyperparameters,
        AngleSMIRKS,
        BondHyperparameters,
        ImproperTorsionHyperparameters,
        ProperTorsionHyperparameters,
        BondSMIRKS,
        ProperTorsionSMIRKS,
        ImproperTorsionSMIRKS,
    )

    torsion_training_set, optimization_training_set = load_training_data(
        optimization_dataset=optimization_dataset,
        torsion_dataset=torsion_dataset,
        smarts_to_exclude=smarts_to_exclude,
        smiles_to_exclude=smiles_to_exclude,
        verbose=verbose
    )

    optimizer = ForceBalanceSchema(
        max_iterations=max_iterations,
        step_convergence_threshold=0.01,
        objective_convergence_threshold=0.1,
        gradient_convergence_threshold=0.1,
        n_criteria=2,
        initial_trust_radius=-1.0,
        finite_difference_h=0.01,
        extras={
            "wq_port": str(port),
            "asynchronous": "True",
            "search_tolerance": "0.1",
            "backup": "0",
            "retain_micro_outputs": "0",
        },
    )

    targets = [
        TorsionProfileTargetSchema(
            reference_data=torsion_training_set,
            energy_denominator=1.0,
            energy_cutoff=8.0,
            extras={"remote": "1"},
        ),
        OptGeoTargetSchema(
            reference_data=optimization_training_set,
            weight=0.01,
            extras={"batch_size": 30, "remote": "1"},
            bond_denominator=0.05,
            angle_denominator=5.0,
            dihedral_denominator=10.0,
            improper_denominator=10.0,
        ),
    ]

    linear_angle_smirks = []
    if frozen_angle_file:
        with open(frozen_angle_file, "r") as f:
            linear_angle_smirks = json.load(f)["smirks"]
    print(f"Frozen angle SMIRKS: {linear_angle_smirks}")
        

    with open(valence_to_optimize, "r") as f:
        valence_smirks = json.load(f)
    with open(torsions_to_optimize, "r") as f:
        torsion_smirks = json.load(f)

    target_parameters = []
    for smirks in valence_smirks["Angles"]:
        if smirks in linear_angle_smirks:
            parameter = AngleSMIRKS(smirks=smirks, attributes={"k"})
        else:
            parameter = AngleSMIRKS(smirks=smirks, attributes={"k", "angle"})
        target_parameters.append(parameter)

    for smirks in valence_smirks["Bonds"]:
        target_parameters.append(BondSMIRKS(smirks=smirks, attributes={"k", "length"}))

    ff = ForceField(forcefield)

    torsion_handler = ff.get_parameter_handler("ProperTorsions")
    for smirks in torsion_smirks["ProperTorsions"]:
        original_k = torsion_handler.parameters[smirks].k
        attributes = {f"k{i + 1}" for i in range(len(original_k))}
        target_parameters.append(
            ProperTorsionSMIRKS(
            smirks=smirks,
            attributes=attributes
            )
        )
    # re-fit impropers directly from force field
    improper_torsion_handler = ff.get_parameter_handler("ImproperTorsions")
    for improper in improper_torsion_handler.parameters:
        original_k = improper.k
        attributes = {f"k{i + 1}" for i in range(len(original_k))}
        target_parameters.append(
            ImproperTorsionSMIRKS(
                smirks=improper.smirks,
                attributes=attributes
            )
        )


    optimization_schema = OptimizationSchema(
        id=tag,
        initial_force_field=os.path.abspath(forcefield),
        stages=[
            OptimizationStageSchema(
                optimizer=optimizer,
                targets=targets,
                parameters=target_parameters,
                parameter_hyperparameters=[
                    AngleHyperparameters(priors={"k": 100.0, "angle": 5.0}),
                    BondHyperparameters(priors={"k": 100.0, "length": 0.1}),
                    ProperTorsionHyperparameters(priors={"k": 5}),
                    ImproperTorsionHyperparameters(priors={"k": 5}),
                ],
            )
        ]
    )

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    optdir = output_directory / "schemas" / "optimizations"
    optdir.mkdir(parents=True, exist_ok=True)

    optfile = optdir / f"{optimization_schema.id}.json"
    with optfile.open("w") as f:
        f.write(optimization_schema.json(indent=2))

    # Generate the ForceBalance inputs
    ForceBalanceInputFactory.generate(
        os.path.join(optimization_schema.id),
        optimization_schema.stages[0],
        ForceField(optimization_schema.initial_force_field),
    )



if __name__ == "__main__":
    generate()
