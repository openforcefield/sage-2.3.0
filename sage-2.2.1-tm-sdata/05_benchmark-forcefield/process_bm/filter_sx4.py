import numpy as np
from matplotlib import pyplot as plt
import copy
# suppress stereochemistry warnings
import logging
logging.getLogger("openff").setLevel(logging.ERROR)
import sys
from openff.qcsubmit.results.filters import SMARTSFilter
from openff.qcsubmit.results import OptimizationResultCollection
from openff.units import unit
import pandas as pd
from openff.toolkit import Molecule, ForceField




smarts_to_filter = [['[*:1]-[#16X4:2](-[*])(-[*])-[*:3]']]

industry_dataset = OptimizationResultCollection.parse_file('../datasets/OpenFF-Industry-v1.1-Sulfur-v1.0-2TM-v1.0-chargecoverage.json')
industry_entries = industry_dataset.entries['https://api.qcarchive.molssi.org:443/']
industry_ids =[rec.record_id for rec in industry_entries]
# industry_ids = np.loadtxt(sys.argv[1])

for smarts_filter_pattern in smarts_to_filter:
    try:
        print('Loaded industry dataset, number of entries: ',industry_dataset.n_results)
        industry_filtered = industry_dataset.filter(SMARTSFilter(smarts_to_include=smarts_filter_pattern))
        print('Filtered dataset. Number of entries: ',industry_filtered.n_results)

        # First, find molecules that were filtered out by looking for record ID's that
        # appear only in the full dataset
        #print(industry_filtered.entries)
        filtered_entries = industry_filtered.entries['https://api.qcarchive.molssi.org:443/']


        print('Number of entries in unfiltered dataset: ',len(industry_ids))
        filtered_ids = [rec.record_id for rec in filtered_entries]
        print('Number of entries in filtered dataset:   ',len(filtered_ids))
        filtered_out_ids = [rec for rec in industry_ids if rec not in filtered_ids]
        print('Number of entries filtered out:          ',len(filtered_out_ids))

        np.savetxt('problem_ids/sx4_outliers.txt'.format(''.join(smarts_filter_pattern)),filtered_ids,fmt='%s')
    except KeyError:
        print('No matches for pattern ',smarts_filter_pattern)
