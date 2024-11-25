#!/bin/bash
#SBATCH -J msm
#SBATCH -p free
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o msm.out
#SBATCH -e msm.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env 

python create-msm-ff.py                                                                                     \
    --initial-force-field       "../01_generate-forcefield/output/initial-force-field.offxml"  \
    --optimization-dataset      "../02_curate-data/output/optimization-training-set.json"                   \
    --hessian-datasets          "OpenFF Gen 2 Opt Set 1 Roche"                      \
    --hessian-datasets          "OpenFF Gen 2 Opt Set 2 Coverage"                   \
    --hessian-datasets          "OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy"         \
    --hessian-datasets          "OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy"     \
    --hessian-datasets          "OpenFF Gen 2 Opt Set 5 Bayer"                      \
    --hessian-datasets          "OpenFF Optimization Set 1"                         \
    --hessian-datasets          "SMIRNOFF Coverage Set 1"                           \
    --hessian-datasets-2          "OpenFF Sulfur Hessian Training Coverage Supplement v1.0" \
    --working-directory         "working-directory"                                                         \
    --output                    "output/initial-force-field-msm.offxml" \
    --verbose True

date
