#!/bin/bash
#SBATCH -J openff_unconstrained-2.2.1-tm-sdata
#SBATCH -p standard
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o openff_unconstrained-2.2.1-tm-sdata.out
#SBATCH -e openff_unconstrained-2.2.1-tm-sdata.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env

savedir="openff_unconstrained-2.2.1-tm-sdata"
mkdir $savedir

python -c "from openff.toolkit.utils import *; assert OpenEyeToolkitWrapper().is_available"

python -u  benchmark.py -f "openff_unconstrained-2.2.1-tm-sdata.offxml" -d "datasets/OpenFF-Industry-v1.1-Sulfur-v1.0-2TM-v1.0-chargecoverage-cache.json" -s "openff_unconstrained-2.2.1-tm-sdata.sqlite" -o $savedir --procs 16

date
