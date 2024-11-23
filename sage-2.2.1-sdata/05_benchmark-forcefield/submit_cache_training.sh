#!/bin/bash
#SBATCH -J cache_training 
#SBATCH -p standard
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6GB
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o cache_training.out 
#SBATCH -e cache_training.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env 

python cache_training.py 8
