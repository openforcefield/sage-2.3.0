from qcportal import PortalClient
from openff.qcsubmit.results import OptimizationResultCollection, TorsionDriveResultCollection
from openff.qcsubmit.results.filters import (
    ConnectivityFilter,
    RecordStatusFilter,)

from openff.toolkit import ForceField

from qcportal.record_models import RecordStatusEnum
import contextlib
import json
import itertools
import math
import pathlib
import typing

import click
from click_option_group import optgroup
import tqdm

import dask
from dask import distributed
dask.config.set({
"distributed.comm.timeouts.tcp": "50s",
"distributed.scheduler.allowed-failures": 999
})




# The majority of this code is from Lily

# ==== helper functions to submit batch workers, copied from openff-nagl ==== #
# +   if NAGL is already installed, you can bypass this by just importing   + #
# +   batch_distributed from openff.nagl.utils._parallelization             + #
# +   although installing NAGL comes with heavy dependencies, e.g. pytorch  + #
# =========================================================================== #

def batch_entries(
    entries: typing.Iterable,
    batch_size: int = 1,
):
    """Convert a flat list of entries into batches to submit"""
    size = batch_size - 1
    entries = iter(entries)
    for x in entries:
        yield list(itertools.chain([x], itertools.islice(entries, size)))


def reconcile_batch_workers(
    entries: typing.Iterable,
    n_entries: typing.Optional[int] = None,
    batch_size: int = -1,
    n_workers: int = -1,
):
    """Calculate how many workers to submit with constraints"""
    if n_entries is None:
        n_entries = len(entries)
    if batch_size < 0:
        n_batches = 1
        batch_size = n_entries
    else:
        n_batches = int(math.ceil(n_entries / batch_size))

    if n_workers < 0:
        n_workers = n_batches
    if n_workers > n_batches:
        n_workers = n_batches

    return n_workers, n_batches, batch_size


@contextlib.contextmanager
def batch_distributed(
    entries: typing.Iterable[typing.Any],
    n_entries: typing.Optional[int] = None,
    batch_size: int = -1,
    n_workers: int = -1,
    worker_type: typing.Literal["lsf", "slurm", "local"] = "slurm",
    queue: str = "free",
    account: typing.Optional[str] = None,
    conda_environment: str = "openff-nagl",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    package_manager: typing.Literal["conda", "micromamba"] = "conda",
    **kwargs
):
    from distributed import LocalCluster
    from dask_jobqueue import LSFCluster, SLURMCluster

    n_workers, n_batches, batch_size = reconcile_batch_workers(
        entries, n_entries, batch_size, n_workers
    )

    print(
        f"Setting n_workers={n_workers} for {n_batches} batches"
    )

    env_extra = dask.config.get(f"jobqueue.{worker_type}.job-script-prologue", default=[])
    if not env_extra:
        env_extra = []
    env_extra.append("source ~/.bashrc")
    env_extra.append(f"{package_manager} activate {conda_environment}")

    if worker_type == "local":
        cluster = LocalCluster(n_workers=n_workers)
    else:
        CLUSTER_CLASSES = {
            "lsf": LSFCluster,
            "slurm": SLURMCluster,
        }
        cluster = CLUSTER_CLASSES[worker_type](
            queue=queue,
            project=account,
            cores=1,
            memory=f"{memory * 1e9}B",
            walltime=f"{walltime}:00:00",
            job_script_prologue=env_extra,
        )
        cluster.scale(n=n_workers)
    
    client = distributed.Client(cluster)

    def wrapper(func, **kwargs):
        for batch in batch_entries(entries, batch_size):
            future = client.submit(func, batch, **kwargs)
            yield future

    try:
        yield wrapper
    finally:
        if worker_type != "local":
            cluster.scale(n=0)


# ==== actual script part ==== #

def check_molecule_can_assign_charges(
    smiles: str,
    charge_backend: typing.Literal["openeye", "ambertools"] = "openeye",
    ff: str = 'openff_unconstrained-2.0.0.offxml'
):
    """
    Quick check to assign partial charges
    """
    from openff.toolkit import Molecule,ForceField
    from openff.toolkit.utils.toolkit_registry import ToolkitRegistry
    from openff.toolkit.utils import AmberToolsToolkitWrapper, RDKitToolkitWrapper, OpenEyeToolkitWrapper

    REGISTRIES = {
        "openeye": OpenEyeToolkitWrapper(),
        "ambertools": ToolkitRegistry([RDKitToolkitWrapper(), AmberToolsToolkitWrapper()]),
    }
    registry = REGISTRIES[charge_backend.lower()]

    # the actual charge check filter converts to SDF first --
    # the code here can be replaced with that if necessary
    # but requires having qcsubmit available, I think!

    # First, assign charges with AM1BCC. Anything that can be charged with AM1BCC should be good for ELF10, too
    offmol = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
    offmol.assign_partial_charges("am1bcc", toolkit_registry=registry)

    # Check for FF coverage by creating an interchange, feed these charges in so that it doesn't do expensive ELF10
    force_field = ForceField(ff)
    system = force_field.create_interchange(offmol.to_topology(),charge_from_molecules=[offmol]).to_openmm(
                     combine_nonbonded_forces=False
                     )



def batch_label(
    entries: list[dict[str, str]],
    charge_backend: typing.Literal["openeye", "ambertools"] = "openeye",
    ff:str = 'openff_unconstrained-2.0.0.offxml'
):
    new_entries = []
    errors = []
    for entry in tqdm.tqdm(entries):
        try:
            check_molecule_can_assign_charges(
                entry["cmiles"],
                charge_backend=charge_backend,
                ff = ff
            )
        except BaseException as e:
            errors.append((entry["cmiles"], e))
        else:
            new_entries.append(entry)


    return new_entries, errors


@click.command()
@click.option(
    "--input",
    "input_dataset",
    multiple=True,
    type=str,
    required=True,
    help="The names of the input dataset(s) to filter.",
)
@click.option(
    "--charge-backend",
    type=click.Choice(["openeye", "ambertools"]),
    required=True,
    help="The charge backend to use",
)
@click.option(
    "--output",
    "output_dataset",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to save the filtered dataset to (*.json).",
)
@click.option(
    "--forcefield",
    "forcefield",
    type=str,
    help='force field to check coverage for'
)   
@optgroup.group("Parallelization configuration")
@optgroup.option(
    "--n-workers",
    help="The number of workers to distribute the labelling across. Use -1 to request "
    "one worker per batch.",
    type=int,
    default=1,
    show_default=True,
)
@optgroup.option(
    "--worker-type",
    help="The type of worker to distribute the labelling across.",
    type=click.Choice(["lsf", "local", "slurm"]),
    default="local",
    show_default=True,
)
@optgroup.option(
    "--batch-size",
    help="The number of molecules to processes at once on a particular worker.",
    type=int,
    default=500,
    show_default=True,
)
@optgroup.group("LSF configuration", help="Options to configure LSF workers.")
@optgroup.option(
    "--memory",
    help="The amount of memory (GB) to request per LSF queue worker.",
    type=int,
    default=3,
    show_default=True,
)
@optgroup.option(
    "--walltime",
    help="The maximum wall-clock hours to request per LSF queue worker.",
    type=int,
    default=2,
    show_default=True,
)
@optgroup.option(
    "--queue",
    help="The LSF queue to submit workers to.",
    type=str,
    default="cpuqueue",
    show_default=True,
)
@optgroup.option(
    "--conda-environment",
    help="The conda environment that LSF workers should run using.",
    type=str,
)
def main(
    input_dataset: str,
    output_dataset: str,
    charge_backend: typing.Literal["openeye", "ambertools"] = "openeye",
    forcefield: str='openff_unconstrained-2.0.0.offxml',
    worker_type: typing.Literal["lsf", "local"] = "local",
    queue: str = "cpuqueue",
    conda_environment: str = "openff-nagl",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):


    client = PortalClient("https://api.qcarchive.molssi.org:443") 
    contents_nofilter = OptimizationResultCollection.from_server( client = client, datasets = input_dataset, spec_name = 'default')

    print(contents_nofilter.n_molecules , ' molecules and ', contents_nofilter.n_results, ' conformers before filtering for record status and connectivity')

    contents = contents_nofilter.filter( RecordStatusFilter(status=RecordStatusEnum.complete),
                                           ConnectivityFilter(tolerance=1.2)
                                         )

    print(contents.n_molecules, ' molecules and ', contents.n_results, ' conformers after filtering for record status and connectivity')

    with open('datasets/filter_record_connec.json','w') as file:
        file.write(contents.json(indent=2))

    with open('datasets/filter_record_connec.json','r') as file:
        contents = json.load(file) 

    all_entries = contents['entries']["https://api.qcarchive.molssi.org:443/"] 


    filtered_entries = []
    all_errors = []
    with batch_distributed(
        all_entries,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(batch_label, charge_backend=charge_backend,ff=forcefield))
        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures),
            desc="Checking charge assignment and coverage",
        ):
            batch, errors = future.result()
            if errors:
                all_errors.extend(errors)
            filtered_entries.extend(batch)

            temp_ds = {'entries': {'https://api.qcarchive.molssi.org:443/': filtered_entries }}
            with open('datasets/temp.json','w') as f:
                json.dump(temp_ds,f,indent=2)

            with open('datasets/temp.err','w') as f:
                for smiles, error in all_errors:
                    f.write(f"SMILES: {smiles}\n{error}\n")
                
    
    print(f"Filtered {len(filtered_entries)} from {len(all_entries)}")
    new_dataset = {
        "entries": {
            "https://api.qcarchive.molssi.org:443/": filtered_entries
        }
    }
    with open(output_dataset, "w") as f:
        json.dump(new_dataset, f, indent=2)
    print(f"Wrote filtered dataset to {output_dataset}")

    output_dataset = pathlib.Path(output_dataset)
    error_file = output_dataset.parent / (output_dataset.stem + "_errors.dat")
    with open(error_file, "w") as file:
        for smiles, error in all_errors:
            file.write(f"SMILES: {smiles}\n{error}\n")


if __name__ == "__main__":
    main()
