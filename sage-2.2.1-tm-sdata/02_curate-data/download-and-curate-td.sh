#!/bin/bash
#SBATCH -J download_and_filter_td
#SBATCH -p standard
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o download_and_filter_td.out
#SBATCH -e download_and_filter_td.err

date
hostname

source ~/.bashrc
conda activate sep-2024-env

python curate-dataset.py download-td                                                \
    --core-td-dataset       "OpenFF Gen 2 Torsion Set 1 Roche 2"                    \
    --core-td-dataset       "OpenFF Gen 2 Torsion Set 2 Coverage 2"                 \
    --core-td-dataset       "OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2"       \
    --core-td-dataset       "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2"   \
    --core-td-dataset       "OpenFF Gen 2 Torsion Set 5 Bayer 2"                    \
    --core-td-dataset       "OpenFF Gen 2 Torsion Set 6 supplemental 2"             \
    --aux-td-dataset        "SMIRNOFF Coverage Torsion Set 1"                       \
    --aux-td-dataset        "OpenFF Group1 Torsions"                                \
    --aux-td-dataset        "OpenFF Group1 Torsions 2"                              \
    --aux-td-dataset        "OpenFF Group1 Torsions 3"                              \
    --aux-td-dataset        "Pfizer discrepancy torsion dataset 1"                  \
    --aux-td-dataset        "OpenFF Gen3 Torsion Set v1.0"                          \
    --aux-td-dataset        "OpenFF Amide Torsion Set v1.0"                         \
    --aux-td-dataset        "OpenFF WBO Conjugated Series v1.0"                     \
    --aux-td-dataset        "OpenFF DANCE 1 eMolecules t142 v1.0"                   \
    --initial-forcefield    "../01_generate-forcefield/output/initial-force-field.offxml" \
    --td-records-to-remove  "td_records_to_remove.dat"                              \
    --additional-td-records "additional_td_records.json"                            \
    --cap-size              5                                                       \
    --cap-method            "pick_random"                                           \
    --n-processes           8                                                       \
    --output                "output/torsion-training-set.json"                      \
    --output-parameter-smirks "output/training-torsion-smirks.json"                 \
    --verbose

date
