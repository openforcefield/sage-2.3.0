#!/bin/bash
#SBATCH -J plot_hist
#SBATCH -p standard
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem 24G
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o plot_hist.out
#SBATCH -e plot_hist.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env-sns

python plot_param_hist.py

date
