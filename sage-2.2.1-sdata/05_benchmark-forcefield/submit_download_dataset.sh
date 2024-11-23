#!/bin/bash
#SBATCH -J download_filter_dataset 
#SBATCH -p standard
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6GB
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o download_filter_dataset.out 
#SBATCH -e download_filter_dataset.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env-dask

python filter_dataset_parallel.py \
    --input                         "OpenFF Industry Benchmark Season 1 v1.1"        \
    --input                         "OpenFF Sulfur Optimization Benchmarking Coverage Supplement v1.0" \
    --output                        "datasets/OpenFF-Industry-v1.1-Sulfur-v1.0-chargecoverage.json"         \
    --charge-backend                "openeye"            \
    --forcefield                    "openff_unconstrained-2.0.0.offxml" \
    --n-workers                     300                     \
    --worker-type                   "slurm"                 \
    --batch-size                    10                      \
    --memory                        6                       \
    --walltime                      48                      \
    --queue                         "free"                  \
    --conda-environment             "sep-2024-env-dask"
