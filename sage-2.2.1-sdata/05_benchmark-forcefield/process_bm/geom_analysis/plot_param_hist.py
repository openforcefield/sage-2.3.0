import pandas as pd
import json
import numpy as np
import tqdm
# For drawing:
from openff.toolkit.topology import FrozenMolecule
from openff.toolkit import Molecule,ForceField

import sys
import matplotlib.pyplot as plt

def mae(array):
    return np.abs(array).mean()


def plot_hist(datasets,labels,filename,title='Industry dataset',xlim=None,ylim=False,legend=True,lw=1,xlab = ''):
    plt.figure()
    plt.axvline(x=0,linestyle='--',color='k',linewidth=0.5)
    for i,data in enumerate(datasets):
        plt.hist(data,histtype='step',bins=50,label=labels[i],linewidth=lw)

    if xlim != None: plt.xlim(xlim[0],xlim[1])
    plt.title(title)
    plt.ylabel('Count')
    plt.xlabel(xlab)
    if legend: plt.legend()
    if np.any(ylim): plt.ylim(ylim[0],ylim[1])
    plt.savefig(filename)
    plt.close()

def boxplot(datasets,labels,filename,title='',ylab = ''):
    # data_dict = {label:datasets[i] for i,label in enumerate(labels)}

    fig,ax = plt.subplots(nrows=1,ncols=1)
    ax.boxplot(datasets,labels=labels)
    left,right=ax.get_xlim()
    ax.hlines(y=0,linestyle='--',color='k',linewidth=0.5,xmin=left,xmax=right)
    # sns.boxplot(data=pd.DataFrame(data_dict))
    ax.set_title(title)
    ax.set_ylabel(ylab)
    plt.savefig(filename)
    plt.close()

suffixs = ['openff_unconstrained-2.1.0_grp_confs','openff_unconstrained-2.2.1_confs' , 'openff_unconstrained-2.2.1-sdata_confs']
angles = [param.id for param in ForceField('openff_unconstrained-2.2.1.offxml').get_parameter_handler("Angles").parameters]
bonds = [param.id for param in ForceField('openff_unconstrained-2.2.1.offxml').get_parameter_handler("Bonds").parameters]

# angles = ['a31','a40']
# bonds = []

# Explicitly changed--might be hard to compare since things leave/arrive
# bonds = ['b13','b13a','b43','b53','b57','b57a']
# angles = ['a13', 'a13a','a13b','a14', 'a14a','a14b','a14c','a15','a15a','a28','a28a','a29','a29a']

bond_errs = {}
angle_errs = {}
proper_errs = {}
improper_errs = {}
for l,suffix in enumerate(suffixs):
    print(suffix)
    bond_errs[suffix]={param : {'errors': []} for param in bonds}
    angle_errs[suffix]={param : {'errors': []} for param in angles}
    # proper_errs[suffix]={param : {'errors': []} for param in propers}
    # improper_errs[suffix]={param : {'errors': []} for param in impropers}

    # First load in all data so it's available easily for toggling
    BOND_DATA_FILE='bonds_qmv{}.json'.format(suffix)
    ANGLE_DATA_FILE='angles_qmv{}.json'.format(suffix)
    # PROPER_DATA_FILE='propers_qmv{}.json'.format(suffix)
    # IMPROPER_DATA_FILE='impropers_qmv{}.json'.format(suffix)

    with open(BOND_DATA_FILE,'r') as jsonfile:
        BOND_JSON = dict(json.load(jsonfile))

    with open(ANGLE_DATA_FILE,'r') as jsonfile:
        ANGLE_JSON = dict(json.load(jsonfile))

    # with open(PROPER_DATA_FILE,'r') as jsonfile:
    #     PROPER_JSON = dict(json.load(jsonfile))
    #
    # with open(IMPROPER_DATA_FILE,'r') as jsonfile:
    #     IMPROPER_JSON = dict(json.load(jsonfile))

    # Would be better to just loop these once, but the smirks changes for each FF so idt it will work
    for smirks ,value in BOND_JSON.items():
        # print(value)
        param_id = value['ident']
        if param_id in bonds:
            bond_errs[suffix][param_id]['errors']=np.array(value['mm_values'] )-np.array(value['qm_values'])
            bond_errs[suffix][param_id]['smirks'] = smirks

    for smirks ,value in ANGLE_JSON.items():
        param_id = value['ident']
        print(param_id)
        if param_id in angles:
            # these are lists/arrays

            # make sure to pick out only the good molecules
            angle_errs[suffix][param_id]['errors']=np.array(value['mm_values'] )-np.array(value['qm_values'])
            angle_errs[suffix][param_id]['smirks'] = smirks

    # for smirks ,value in PROPER_JSON.items():
    #     param_id = value['ident']
    #     if param_id in propers:
    #         # print(param_id)
    #         mm_prop = np.array(value['mm_values'] )
    #         qm_prop = np.array(value['qm_values'] )
    #         mm_min_qm = mm_prop - qm_prop
    #
    #         # try shifting with minimum image convention
    #         # angles are already between -180 to 180
    #         x_size = 360
    #
    #         mm_min_qm[mm_min_qm < -x_size*0.5] += x_size
    #         mm_min_qm[mm_min_qm >= x_size*0.5] -= x_size
    #
    #         proper_errs[suffix][param_id]['errors']=mm_min_qm
    #         proper_errs[suffix][param_id]['smirks'] = smirks
    #         # print(mm_prop[np.abs(mm_min_qm) > 100])
    #         # print(qm_prop[np.abs(mm_min_qm) > 100])
    #         # print uncorrected
    #         # print(np.array(value['mm_values'] )[np.abs(mm_min_qm) > 100])
    #         # print(np.array(value['qm_values'] )[np.abs(mm_min_qm) > 100])
    #
    # for smirks ,value in IMPROPER_JSON.items():
    #     param_id = value['ident']
    #     if param_id in impropers:
    #         # print(param_id)
    #         mm_prop = np.array(value['mm_values'] )
    #         qm_prop = np.array(value['qm_values'] )
    #         mm_min_qm = mm_prop - qm_prop
    #
    #         # try shifting with minimum image convention
    #         # angles are already between -180 to 180
    #         x_size = 360
    #
    #         mm_min_qm[mm_min_qm < -x_size*0.5] += x_size
    #         mm_min_qm[mm_min_qm >= x_size*0.5] -= x_size
    #
    #         improper_errs[suffix][param_id]['errors']=mm_min_qm
    #         improper_errs[suffix][param_id]['smirks'] = smirks

labels = ['Sage 2.1.0 (regroup)','Sage 2.2.1','Sage 2.2.1-sdata' ]
# labels = suffixs
for i,bond in enumerate(bonds):
    try:
        bond_errs_x =[bond_errs[x][bond]['errors'] for x in suffixs]
        try:
            smirks = bond_errs[suffixs[1]][bond]['smirks']
        except KeyError:
            smirks = bond_errs[suffixs[0]][bond]['smirks']
        plot_hist(bond_errs_x,labels,filename='params_incl210_grp/'+bond+ '_hist.pdf',title=bond + ' ' + smirks,xlab='Bond length error (A)')
        boxplot(bond_errs_x,labels,filename='params_incl210_grp/'+bond+ '_boxplot.pdf',title=bond + ' ' + smirks ,ylab='Bond length error (A)' )
    except KeyError:
        print(bond)
        pass

for i,bond in enumerate(angles):
    bond_errs_x =[angle_errs[x][bond]['errors'] for x in suffixs]
    try:
        smirks = angle_errs[suffixs[1]][bond]['smirks']

    except KeyError:
        try:
            smirks = angle_errs[suffixs[0]][bond]['smirks']
        except KeyError:
            print(bond)
            pass
            # try:
            #     smirks = angle_errs[suffixs[3]][bond]['smirks']
            # except KeyError:
            #     smirks = angle_errs[suffixs[4]][bond]['smirks']

    plot_hist(bond_errs_x,labels,filename='params_incl210_grp/'+bond+ '_hist.pdf',title=bond + ' ' + smirks,xlab='Angle error (deg)',xlim=[-25,25])
    boxplot(bond_errs_x,labels,filename='params_incl210_grp/'+bond+ '_boxplot.pdf',title=bond + ' ' + smirks ,ylab = 'Angle error (deg)')

# for i,bond in enumerate(propers):
#     bond_errs_x =[proper_errs[x][bond]['errors'] for x in suffixs]
#     try:
#         smirks = proper_errs[suffixs[1]][bond]['smirks']
#     except KeyError:
#         try:
#             smirks = proper_errs[suffixs[2]][bond]['smirks']
#         except KeyError:
#             smirks = ''
#
#     plot_hist(bond_errs_x,labels,filename='params/'+bond+ '_hist_pbc.pdf',title=bond + ' ' + smirks,xlab='Dihedral angle error (deg)')
#     boxplot(bond_errs_x,labels,filename='params/'+bond+ '_boxplot_pbc.pdf',title=bond + ' ' + smirks  ,ylab= 'Dihedral error (deg)')
#
# for i,bond in enumerate(impropers):
#     bond_errs_x =[improper_errs[x][bond]['errors'] for x in suffixs]
#     try:
#         smirks = improper_errs[suffixs[1]][bond]['smirks']
#     except KeyError:
#         smirks = improper_errs[suffixs[2]][bond]['smirks']
#
#     plot_hist(bond_errs_x,labels,filename='params/'+bond+ '_hist_pbc.pdf',title=bond + ' ' + smirks,xlab='Improper error (deg)')
#     boxplot(bond_errs_x,labels,filename='params/'+bond+ '_boxplot_pbc.pdf',title=bond + ' ' + smirks ,ylab = 'Improper error (deg)' )
