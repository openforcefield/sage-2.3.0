#!/bin/bash

# Only first conformer used, nothing filtered out
python ../calculate_params.py --db OpenFF_Gen2_Coverage_sage220.sqlite --conformers False --ff_yammbs openff_unconstrained-2.2.0.offxml --ff_file openff_unconstrained-2.2.0.offxml --label sage220-fromsqlite-noconfs-noprobs

## Only first conformer used, random QCAIDs filtered out
#python ../calculate_params.py --db OpenFF_Gen2_Coverage_sage220.sqlite --conformers False --ff_yammbs openff_unconstrained-2.2.0.offxml --ff_file openff_unconstrained-2.2.0.offxml --problem_file problem_ids/random_qcaids.txt --label sage220-fromsqlite-noconfs-probs
#
## All conformers used, nothing filtered out                                
#python ../calculate_params.py --db OpenFF_Gen2_Coverage_sage220.sqlite --conformers True --ff_yammbs openff_unconstrained-2.2.0.offxml --ff_file openff_unconstrained-2.2.0.offxml --problem_file problem_ids/random_qcaids.txt --label sage220-fromsqlite-confs-noprobs 
# 
## All conformers used, random QCAIDs filtered out
#python ../calculate_params.py --db OpenFF_Gen2_Coverage_sage220.sqlite --conformers True  --ff_yammbs openff_unconstrained-2.2.0.offxml --ff_file openff_unconstrained-2.2.0.offxml --problem_file problem_ids/random_qcaids.txt --label sage220-fromsqlite-confs-probs  
