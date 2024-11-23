#!/bin/bash
#SBATCH -J filter_bm_smarts
#SBATCH -p standard
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o filter_bm_smarts.out
#SBATCH -e filter_bm_smarts.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env

python filter_bm_smarts.py

date
