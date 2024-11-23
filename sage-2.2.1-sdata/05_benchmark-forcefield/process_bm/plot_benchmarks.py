import matplotlib.pyplot as plt
import numpy as np
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import pandas as pd

import seaborn as sns
import sys
import click
import os
# import warnings
warnings.filterwarnings(
    "ignore", message="divide by zero", category=RuntimeWarning
)

def plot_hist(datasets,labels,filename,title='Industry dataset',xlim=[-2,2],ylim=False,legend=True,lw=1):
    plt.figure()
    plt.axvline(x=0,linestyle='--',color='k',linewidth=0.5)
    for i,data in enumerate(datasets):
        plt.hist(data,histtype='step',bins=400,range=(-100.25,99.75),label=labels[i],linewidth=lw)

    plt.xlim(xlim[0],xlim[1])
    plt.title(title)
    plt.ylabel('Count')
    plt.xlabel('DDE (kcal/mol)')
    if legend: plt.legend()
    if np.any(ylim): plt.ylim(ylim[0],ylim[1])
    plt.savefig(filename)

def plot_kde(datasets,labels,filename,title='Industry dataset',bw=1,xlim=[-2,2],cumulative=False):
    if bw > .5:
        gs = 1
    elif bw > .2:
        gs = 2
    else:
        gs = 3

    # Industry
    plt.figure()
    plt.axvline(x=0,linestyle='--',color='k',linewidth=0.5)
    for i,data in enumerate(datasets):
        sns.kdeplot(data,bw_adjust=bw,gridsize=10000*gs,label=labels[i],cumulative=cumulative)

    plt.xlim(xlim[0],xlim[1])
    plt.title(title)
    plt.ylabel('Density')
    plt.xlabel('DDE (kcal/mol)')
    plt.legend()
    plt.savefig(filename)


def plot_cdf(datasets,labels,filename,title='Industry dataset',bw=1,xlim=[-2,2],xaxis='RMSD (A)'):

    # Industry
    fig,ax = plt.subplots(1)
    ax.axvline(x=0,linestyle='--',color='k',linewidth=0.5)
    ax.axhline(y=1,linestyle = '--',color='k',linewidth=0.5)
    ax.axhline(y=0.75,linestyle = '--',color='k',linewidth=0.5)
    for i,data in enumerate(datasets):

        ax.ecdf(np.abs(data[np.logical_not(np.isnan(data))]),label=labels[i])
    ax.set_xlim(xlim[0],xlim[1])
    ax.set_title(title)
    ax.set_ylabel('Density')
    ax.set_xlabel(xaxis)
    ax.legend()
    plt.savefig(filename)


def mae(array):
    return np.abs(array).mean()

def get_outliers(all_data,all_data_names,artificial_low=False,artificial_high=False):
    # Get actual outliers--molecules that are outside the IQR
    outlier_ns = []
    for i,data in enumerate(all_data):
        q1, q3 = np.percentile(data, [25, 75])
        whisker_low = q1 - (q3 - q1) * 1.5
        whisker_high = q3 + (q3 - q1) * 1.5
        if artificial_low: whisker_low = artificial_low
        if artificial_high: whisker_high = artificial_high
        outliers = data[(data > whisker_high) | (data < whisker_low)]
        # outlier_ids = all_ids[i][(data > whisker_high) | (data < whisker_low)]
        # if save: np.savetxt(all_data_names[i]+'_outlierids.txt',outlier_ids,fmt='%i')
        outlier_ns.append(outliers.shape[0])
    return outlier_ns

def print_stats(datasets,labels,lines,outlier_low=False,outlier_high=False):
    # print(datasets)
    # print('Stats:')
    lines.append('{:<25s}   {:>8s} {:>10s} {:>10s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}'.format('Dataset','Size','Min','Max','MSE','MAE','Median','STD','# outliers'))
    longest_label = max([len(label) for label in labels])
    # print(longest_label)
    outliers = get_outliers(datasets,labels,artificial_low=outlier_low,artificial_high=outlier_high)
    for i,ds in enumerate(datasets):
        # print(data)
        data = ds[~np.isnan(ds)]

        lines.append('{:<25s}   {:>-8d} {:>-10.4f} {:>-10.4f} {:>-8.5f} {:>-8.5f} {:>-8.5f} {:>-8.5f} {:>-8d}'.format(labels[i],data.shape[0],data.min(),data.max(),data.mean(),mae(data),np.median(data),data.std(),outliers[i]))


def boxplot(datasets,labels,filename,title=''):
    data_dict = {label:datasets[i] for i,label in enumerate(labels)}
    plt.figure()
    sns.boxplot(data=pd.DataFrame(data_dict))
    plt.ylim(0,6)
    plt.title(title)
    plt.savefig(filename)

def plot_log_kde(datasets,labels,filename,title='Industry dataset',bw=1,xlim=[-2,1],xlabel='RMSD (Log A)'):
    if bw > .5:
        gs = 1
    elif bw > .2:
        gs = 2
    else:
        gs = 3

    # Industry
    plt.figure()
    plt.axvline(x=0,linestyle='--',color='k',linewidth=0.5)
    for i,data in enumerate(datasets):
        sns.kdeplot(np.log10(data),bw_adjust=bw,gridsize=10000*gs,label=labels[i])

    plt.xlim(xlim[0],xlim[1])
    plt.title(title)
    plt.ylabel('Density')
    plt.xlabel('Log10 '+xlabel)
    plt.legend()
    plt.savefig(filename)
    plt.close()

def plot_conf(all_data,all_data_names,dir,type='rmsd'):
    if type == 'tfd':
        xlab = 'TFD'
        xlim1 = [0,0.5]
        xlim2 = [-5,0.2]
    else:
        xlab = "RMSD (Log A)"
        xlim1 = [0,1]
        xlim2 = [-2,1]
    plot_log_kde(all_data,all_data_names,dir+'{}_kde.pdf'.format(type),xlabel=xlab,xlim=xlim2)
    plt.close()
    plot_cdf(all_data,all_data_names,dir+'{}_cdf.pdf'.format(type),xlim=xlim1,xaxis=xlab.replace('Log ',''))
    plt.close()
    try:
        boxplot(all_data,all_data_names,dir+'{}_boxplot.pdf'.format(type))
        plt.close()
    except ValueError:
        pass

def plot_dde(all_data,all_data_names,dir,type,ic_type,xlim=[-2,2]):
    if type == 'dde':
        plttype = 'dde'
    elif type == 'icrmsd':
        plttype = ic_type
        xlim = [0.9*min([dat.min() for dat in all_data]),max([dat.max() for dat in all_data])*1.1]
    try:
        plot_hist(all_data,all_data_names,dir+'{}_hist.pdf'.format(plttype),xlim=xlim)
        plt.close()
        plot_hist(all_data,all_data_names,dir+'{}_hist_zoom.pdf'.format(plttype),xlim=[-0.5,0.5],lw=2)
        plt.close()
        plot_kde(all_data,all_data_names,dir+'{}_kde.pdf'.format(plttype),xlim=xlim)
        plt.close()
        if plttype == 'dde':
            plot_cdf(all_data,all_data_names,dir+'{}_cdf.pdf'.format(type),xlim=[-0.5,5],xaxis='|ddE| (kcal/mol)')
            plt.close()
    except ValueError:
        pass
    try:
        boxplot(all_data,all_data_names,dir+'{}_boxplot.pdf'.format(plttype))
        plt.close()
    except ValueError as e:
        print(e)
        print([d.shape for d in all_data])
        pass

def plot_conf_filter(kept_by_filter,filtered_out,all_data_names,dir,filter_pattern,type='rmsd'):
    if type == 'tfd':
        xlab = 'TFD'
        xlim1 = [0,0.5]
        xlim2 = [-5,0.2]
    else:
        xlab = "RMSD (Log A)"
        xlim1 = [0,1]
        xlim2 = [-2,1]



    title = 'Molecules with {}'.format(filter_pattern)
    plot_log_kde(kept_by_filter,all_data_names,dir+'{}_{}_kde.pdf'.format(type,filter_pattern),title=title,xlabel=xlab,xlim=xlim2)
    plt.close()
    boxplot(kept_by_filter,all_data_names,dir+'{}_{}_boxplot.pdf'.format(type,filter_pattern),title=title)
    plt.close()
    plot_cdf(kept_by_filter,all_data_names,dir+'{}_{}_cdf.pdf'.format(type,filter_pattern),xlim=xlim1,xaxis=xlab.replace('Log ',''))
    plt.close()

    title = 'Molecules without {}'.format(filter_pattern)
    plot_log_kde(filtered_out,all_data_names,dir+'{}_{}fo_kde.pdf'.format(type,filter_pattern),title=title,xlabel=xlab,xlim=xlim2)
    plt.close()
    boxplot(filtered_out,all_data_names,dir+'{}_{}fo_boxplot.pdf'.format(type,filter_pattern),title=title)
    plt.close()
    plot_cdf(filtered_out,all_data_names,dir+'{}_{}fo_cdf.pdf'.format(type,filter_pattern),xlim=xlim1,xaxis=xlab.replace('Log ',''))
    plt.close()

def plot_dde_filter(kept_by_filter,filtered_out,all_data_names,dir,filter_pattern,type,ic_type):
    if type == 'dde':
        plttype = 'dde'
    elif type == 'icrmsd':
        plttype = ic_type

    # kept by filter
    title = 'Molecules with {}'.format(filter_pattern)
    plot_hist(kept_by_filter,all_data_names,dir+'{}_{}_hist.pdf'.format(plttype,filter_pattern),title=title)
    plt.close()
    plot_hist(kept_by_filter,all_data_names,dir+'{}_{}_hist_zoom.pdf'.format(plttype,filter_pattern),xlim=[-0.5,0.5],lw=2,title=title)
    plt.close()
    plot_kde(kept_by_filter,all_data_names,dir+'{}_{}_kde.pdf'.format(plttype,filter_pattern),title=title)
    plt.close()
    if plttype == 'dde':
        plot_cdf(kept_by_filter,all_data_names,dir+'{}_{}_cdf.pdf'.format(plttype,filter_pattern),xlim=[-0.5,5],xaxis='|ddE| (kcal/mol)')
        plt.close()
    try:
        boxplot(kept_by_filter,all_data_names,dir+'{}_{}_boxplot.pdf'.format(plttype,filter_pattern),title=title)
        plt.close()
    except ValueError:
        pass

    title = 'Molecules without {}'.format(filter_pattern)
    plot_hist(filtered_out,all_data_names,dir+'{}_{}fo_hist.pdf'.format(plttype,filter_pattern),title=title)
    plt.close()
    plot_hist(filtered_out,all_data_names,dir+'{}_{}fo_hist_zoom.pdf'.format(plttype,filter_pattern),xlim=[-0.5,0.5],lw=2,title=title)
    plt.close()
    plot_kde(filtered_out,all_data_names,dir+'{}_{}fo_kde.pdf'.format(plttype,filter_pattern),title=title)
    plt.close()
    if plttype == 'dde':
        plot_cdf(filtered_out,all_data_names,dir+'{}_{}fo_cdf.pdf'.format(plttype,filter_pattern),xlim=[-0.5,5],xaxis='|ddE| (kcal/mol)')
        plt.close()
    try:
        boxplot(filtered_out,all_data_names,dir+'{}_{}fo_boxplot.pdf'.format(plttype,filter_pattern),title=title)
        plt.close()
    except ValueError:
        pass

# def plot_dde()

def filter_data_difsize(datasets,dataset_ids,filter_ids,lines):#,rmsd_filter_idx):
    lines.append("Data has different sizes, filtering datasets individually")
    # loop through datasets
    all_data_covered = []
    all_data_notcovered = []
    all_ids_covered = []
    all_ids_notcovered = []
    filter_id_files = [np.loadtxt(file) for file in filter_ids]
    for i,data in enumerate(datasets):
        did = np.array(dataset_ids[i])

        # print(data.shape)
        # print(did.shape)
        # print(rmsd_idx.shape)
        # print(np.count_nonzero(rmsd_idx))

        # loop through filter files--maybe not the most efficient
        n_data_flag = False
        for filter_id in filter_id_files:
            # filter_id = np.loadtxt(filter_id_file)

            # print(filter_id)
            if n_data_flag: # if we've already filtered once, need to add these in
                covered_idx = np.logical_or(covered_idx,np.array([i in filter_id for i in did]))
            else:
                covered_idx = np.array([i in filter_id for i in did])
                n_data_flag = True

        # print(covered_idx.shape)
        # these will just be all "True" if no rmsd filter, so should be OK
        not_covered_idx = np.logical_not(covered_idx)#[rmsd_idx]
        covered_idx = covered_idx#[rmsd_idx]
        # print(covered_idx.shape)

        # don't need to do any asserting since we are doing each one in its own order
        all_data_covered.append(data[covered_idx]) # data that includes the filtered IDS
        all_data_notcovered.append(data[not_covered_idx]) # data that does not include the filtere IDs
        all_ids_covered.append(did[covered_idx])#[rmsd_idx][covered_idx])
        all_ids_notcovered.append(did[not_covered_idx])#[rmsd_idx][not_covered_idx])

    return all_data_covered,all_data_notcovered,all_ids_covered,all_ids_notcovered

def filter_data(datasets,dataset_ids,filter_ids,lines): #,rmsd_filter_idx):
    try:
        sage_ids = dataset_ids[0] # Assume Sage is first

        n_data_flag = False
        # loop through all filter ID files
        for filter_id_file in filter_ids:
            filter_id = np.loadtxt(filter_id_file)

            # print(filter_id)
            if n_data_flag: # if we've already filtered once, need to add these in
                covered_idx = np.logical_or(covered_idx,np.array([i in filter_id for i in sage_ids]))
            else:
                covered_idx = np.array([i in filter_id for i in sage_ids])
                n_data_flag = True

        not_covered_idx = np.logical_not(covered_idx)

        # all_covered_idx = [covered_idx[rmsd_idx] for rmsd_idx in rmsd_filter_idx]
        # all_not_covered_idx = [not_covered_idx[rmsd_idx] for rmsd_idx in rmsd_filter_idx]
        # for x in all_covered_idx:
        #     print(x.shape)

        all_data_covered = []
        all_data_notcovered = []
        all_ids_covered = []
        all_ids_notcovered = []
        for i,data in enumerate(datasets):
            # covered_idx = all_covered_idx[i]
            # not_covered_idx = all_not_covered_idx[i]
            assert np.all(dataset_ids[i]==sage_ids) # Make sure they're all in the same order

            all_data_covered.append(data[covered_idx]) # data that includes the filtered IDS
            all_data_notcovered.append(data[not_covered_idx]) # data that does not include the filtere IDs
            all_ids_covered.append(dataset_ids[i][covered_idx])
            all_ids_notcovered.append(dataset_ids[i][not_covered_idx])

        return all_data_covered,all_data_notcovered,all_ids_covered,all_ids_notcovered
    except (ValueError,IndexError):
        return filter_data_difsize(datasets,dataset_ids,filter_ids,lines)#,rmsd_filter_idx)



def main_filter(dir,all_data,all_data_names,all_ids,filter_file,filter_pattern,type,ic_type,lines):#,rmsd_idx): # outliers have already been excluded
    dir = dir+filter_pattern+'/'
    if not os.path.isdir(dir):
    # except FileNotFoundError:
        os.mkdir(dir)
    dir += type+'/'
    if not os.path.isdir(dir): # make parent dir if necessary
    # except FileNotFoundError:
        os.mkdir(dir)


    lines.append('Filtering dataset')
    kept_by_filter,filtered_out,kept_ids,filtered_out_ids = filter_data(all_data,all_ids,filter_file,lines)#,rmsd_idx)
    lines.append('Number of molecules kept by filter: ' + str(len(kept_ids[0])))
    lines.append('Number of molecules filtered out: ' + str(len(filtered_out_ids[0])))
    lines.append('')

    lines.append('')
    lines.append('Stats for molecules that include {}'.format(filter_pattern))
    print_stats(kept_by_filter,all_data_names,lines)
    lines.append('')
    lines.append("Stats for molecules that don't include {}".format(filter_pattern))
    print_stats(filtered_out,all_data_names,lines)

    lines.append('Saving plots to '+ dir)

    if type == 'icrmsd':
        with open(dir+'/'+type+'_'+ic_type+'.log','w') as writefile:
            for line in lines:
                writefile.write(line+'\n')
    else:
        with open(dir+'/'+type+'.log','w') as writefile:
            for line in lines:
                writefile.write(line+'\n')

    if type == 'dde' or type == 'icrmsd':
        plot_dde_filter(kept_by_filter,filtered_out,all_data_names,dir,filter_pattern,type,ic_type)
    else:
        plot_conf_filter(kept_by_filter,filtered_out,all_data_names,dir,filter_pattern,type=type)

def main_nofilter(dir,all_data,all_data_names,type,ic_type,lines):
    dir += type+'/'
    if not os.path.isdir(dir): # make parent dir if necessary
    # except FileNotFoundError:
        os.mkdir(dir)

    lines.append('Saving plots to '+ dir)

    if type == 'icrmsd':
        with open(dir+'/'+type+'_'+ic_type+'.log','w') as writefile:
            for line in lines:
                writefile.write(line+'\n')
    else:
        with open(dir+'/'+type+'.log','w') as writefile:
            for line in lines:
                writefile.write(line+'\n')

    if type == 'dde' or type == 'icrmsd':
        plot_dde(all_data,all_data_names,dir,type,ic_type)
    else:
        plot_conf(all_data,all_data_names,dir,type=type)

@click.command()
@click.option('--data_dirs',multiple=True)
@click.option('--save_dir',default='./',help='Directory to save plots to')
@click.option('--problem_files',default=[],multiple=True,help='File with problem IDs to filter out')
@click.option('--filter_file',default=[],multiple=True,help='File with IDs of records to filter. Should contain molecules kept by filter. ')
@click.option('--filter_pattern',default=None,help='Description of filter to be used for titles etc')
@click.option('--type',type=click.Choice(['dde','rmsd','tfd','icrmsd']),default='dde',help='Property to plot')
@click.option('--ic_type',type=click.Choice(['bond','angle','dihedral','improper']),default='bond',help='which ICRMSD to plot')
def main(data_dirs,save_dir,problem_files,filter_file,filter_pattern,type,ic_type):
    dir = save_dir
    # start keeping track of lines for a log file
    lines = []
    # make sure we can use dir by adding a slash and making the directory if necessary
    if dir[-1] != '/':
            dir += '/'

    if not os.path.isdir(dir): # make parent dir if necessary
    # except FileNotFoundError:
        os.mkdir(dir)

    if type == 'icrmsd':
        icrmsd_idx = {'bond':0,'angle':1,'dihedral':2,'improper':3}[ic_type] + 1
    else:
        icrmsd_idx = 1


    all_data = []
    all_ids = []
    all_data_names = []
    for ff in data_dirs:
        ids,data = np.genfromtxt('{}/{}.csv'.format(ff,type),delimiter = ',',skip_header=1,missing_values='',unpack=True,usecols=(0,icrmsd_idx))
        all_data.append(data)
        all_ids.append(ids)
        all_data_names.append(ff.replace('../',''))

    if len(problem_files)>0:
        dir += 'problems_removed/'
        if not os.path.isdir(dir): # make parent dir if necessary
        # except FileNotFoundError:
            os.mkdir(dir)

        lines.append('Excluding problem molecules: '+str(problem_files))
        all_outliers,all_data,all_outlier_ids,all_ids = filter_data(all_data,all_ids,problem_files,lines)#,rmsd_idx)
        lines.append('Stats for all data, with problem molecules excluded:')
        # dir += 'outliers'
    else:
        lines.append('Stats for all data:')
    print_stats(all_data,all_data_names,lines)

    # print(filter_file,len(filter_file))
    if len(filter_file)>0:
        main_filter(dir,all_data,all_data_names,all_ids,filter_file,filter_pattern,type,ic_type,lines)#,rmsd_idx)
    else:
        main_nofilter(dir,all_data,all_data_names,type,ic_type,lines)



#dde
if __name__ == '__main__':

    main()

