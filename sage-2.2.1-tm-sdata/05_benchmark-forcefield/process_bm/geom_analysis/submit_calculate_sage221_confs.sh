#!/bin/bash
#SBATCH -J calc_params_221_confs
#SBATCH -p standard
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 6G
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o calc_params_221_confs.out
#SBATCH -e calc_params_221_confs.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env-dask

python calculate_params.py --db ../../openff_unconstrained-2.2.1.sqlite --ff_file openff_unconstrained-2.2.1.offxml --label openff_unconstrained-2.2.1_confs --problem_file ../problem_ids/all_r7_outliers.txt --problem_file ../problem_ids/all_r8_outliers.txt --problem_file ../problem_ids/sx4_outliers.txt --problem_file ../problem_ids/proton_rearr.txt --conformers True

date
