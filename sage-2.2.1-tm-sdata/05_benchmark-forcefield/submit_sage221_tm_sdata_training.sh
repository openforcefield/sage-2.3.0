#!/bin/bash
#SBATCH -J openff_unconstrained-2.2.1-sdata_training
#SBATCH -p standard
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o openff_unconstrained-2.2.1-sdata_training.out
#SBATCH -e openff_unconstrained-2.2.1-sdata_training.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env

savedir="openff_unconstrained-2.2.1-sdata_training"
mkdir $savedir

python -c "from openff.toolkit.utils import *; assert OpenEyeToolkitWrapper().is_available"

python -u  benchmark.py -f "openff_unconstrained-2.2.1-sdata.offxml" -d "datasets/optimization-training-set-cache.json" -s "openff_unconstrained-2.2.1-sdata_training.sqlite" -o $savedir --procs 16

date
