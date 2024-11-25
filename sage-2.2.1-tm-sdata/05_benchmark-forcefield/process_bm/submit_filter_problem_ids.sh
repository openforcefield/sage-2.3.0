#!/bin/bash
#SBATCH -J filter_sx4
#SBATCH -p standard
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o filter_sx4.out
#SBATCH -e filter_sx4.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env 

python filter_sx4.py

date
