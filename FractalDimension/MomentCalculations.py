import pathlib
cwd = pathlib.Path.cwd()
import os
import re
import numpy as np
import glob
from plotly import graph_objects as go
import timeit
import DistanceClass as distance
# from scipy.spatial import distance
# from statsmodels.graphics.tsaplots import plot_acf
import pandas
import pickle
import itertools
import matplotlib.pyplot as plt
import GeneClass as Gene
import random
import Heatmaps as hm
plt.rcParams['text.usetex'] = True


def main():
    '''
    Look i knew I keep doing this but never import the things... the project got rather intense and I never was able to keep up with all
    the demands while keeping everything clean. There should be an import script so I don't recrete the thing over and over again, and probably
    a few others, but oh well.

    Compare exon k-mer to k-1 -mer.

    do log transform on y axis for moment plots (makes it easier to view).
    '''
    # recover_data(cwd / "Dicts")
    # exit()
    k = 4

    exon_file = cwd / "Dicts" / f"Exon_{k}mer.pkl"
    intron_file = cwd / "Dicts" / f"Intron_{k}mer.pkl"

    # data = hm._import_data(exon_file, just_import=True)
    # data = pandas.DataFrame(data)
    # data = data.pow(2.1)
    # print(data)
    # value = data.to_numpy().sum()
    # print(value**(1/2.1))
    # # print(data)
    # print(type(data))
    step = 100
    max_n = 2
    min_n = 0.5
    ms = [n / step for n in range(int(step*(min_n)), int(step*max_n) + 1)]
    title = f"E v I log w/N divided"
    moments_v2(cwd / "Dicts", ms, min_k = 1, max_k = 6, convergence = False, logy = False, N_value = True, title = title, x_ticks={0.5: 0.5, 1.0: 1.0, 1.5: 1.5, 2.0: 2.0})
    # evalue = moment(exon_file, unlog = True, k = 2*k, n = 1)
    # ivalue = moment(intron_file, unlog = True, k = 2*k, n = 1)
    # print(f"k = {k}\te value = {evalue}\ti = value = {ivalue}")


def moments_v2(file_dir: pathlib.Path, 
               ms: list, 
               min_k: int = 1, max_k: int = 6, 
               unlog: bool = False, convergence: bool = False, logy: bool = False, N_value: bool = False, 
               title: str = None, 
               x_ticks: dict = None, y_ticks: dict = None,
               *args, **kwargs):
    '''
    Moment calculations do NOT need the 4**s part multiplied, everything just needs to be scaled from 0 to 1 in a traditional normalization.*

    You have to feed in the list of moments, because I don't want to deal with a for float loop... I have that code somewhere I just can't find it. It's a for loop with floats, and it's pretty easy but there's a trick and I just don't want to do it... should be called the run on sentence bandit...


    Does a for loop of this. Each file needs to be a pickle file and have the same basic file name. Just give it the path and the filename and it will iterate over it.

    Exon files should be f'Exon_{k}mer.pkl' and Intron files should be f'Intron_{k}mer.pkl'

    The for loop for the k iterations is range(min_k, max_k + 1) so you do not need to artifically increase the max_k. That's done for you.

    *Do this as well with a N**(1/m - 1) term where N = 4**k (or in this case 4**2k) which should set all moment calculations for an evenly distributed heatmap to 1 (no matter what the moment is) and will
    adjust the other moment calculations to be either above or below 1. Sort of like subtracting the mean.
    '''

    if convergence and (max_k - min_k) < 1:
        print("Cannot continue, convergence can only be checked with more then 1 k value")
        return None 

    me = {}
    mi = {}
    pd = {}

    uni = {}

    for k in range(min_k, max_k + 1):
        print(f"K = {k}")

        me[k]: list = []
        mi[k]: list = []
        pd[k]: list = []
        uni[k]: list = []

        exon_file = file_dir / f"Exon_{k}mer.pkl"
        intron_file = file_dir / f"Intron_{k}mer.pkl"

        exon_data = hm._import_data(exon_file, just_import = True)
        exon_data = pandas.DataFrame(exon_data)
        if unlog:
            exon_data = _unrenormalize(exon_data, 2*k)
        print(f"\tImported Exon File")

        intron_data = hm._import_data(intron_file, just_import = True)
        intron_data = pandas.DataFrame(intron_data)
        if unlog:
            intron_data = _unrenormalize(intron_data, 2*k)
        print(f"\tImported Intron File")

        for m in ms:

            if N_value:
                N = 4**(2*k)
                N = N**((1/m) - 1)
                print(f"\t\tm = {m}\tN = {N}")
            else:
                N = 1
                print(f"\t\tm = {m}")

            exon_v = moment(exon_data, m = m, unlog = False, N = N, *args, **kwargs)
            # exon_u = uniform_density(exon_data, m = m, N = N)
            intron_v = moment(intron_data, m = m, unlog = False, N = N, *args, **kwargs)
            # intron_u = uniform_density(intron_data, m = m, N = N, *args, **kwargs)

            # if N_value:
            #     N = 4**(2*k)
            #     N = N**((1/m) - 1)
            #     exon_v = N * exon_v
            #     intron_v = N * intron_v

            if logy:
                exon_v = np.log2(exon_v)
                intron_v = np.log2(intron_v)
            
            try:
                per_diff = _percent_difference(exon_v, intron_v)
            except ZeroDivisionError:
                per_diff = 0
            
            me[k].append(exon_v)
            mi[k].append(intron_v)
            pd[k].append(per_diff)

            uni[k].append(1)

    print(f"Finished data, printing to plots")
    plt.rc("font", size = 20)
    plt.rc('axes', linewidth = 2)
    for key, value in me.items():
        print(f"plotting for {2*key}-mer")
        fig, axs = plt.subplots()
        fig.set_size_inches(8, 8)
        fig_pd, axs_pd = plt.subplots()

        axs.plot(ms, value, label = f"Exon")
        axs.plot(ms, mi[key], label = f"Intron")
        axs.plot(ms, uni[key], linestyle = "dotted")
        axs.set_title(f"{2*key}-mer")
        axs.set_ylabel(r"$(\sum\rho^{m})^{1/m}$") #, rotation='horizontal')
        axs.set_xlabel(r"$m$")
        axs.legend()

        if isinstance(x_ticks, dict):
            axs.set_xticks(list(x_ticks.keys()))
            axs.set_xticklabels(list(x_ticks.values()))
        else:
            axs.set_xticks([])
        if isinstance(y_ticks, dict):
            axs.set_yticks(list(y_ticks.keys()))
            axs.set_yticklabels(list(y_ticks.values()))
        else:
            axs.set_yticks([])

        # axs.yaxis.set_label_coords(-0.0, 0.5)

        axs_pd.plot(ms, pd[key])
        axs_pd.set_title(f"Per Diff for {2*key}-mer")
        axs_pd.set_xlabel("m")

        fig.savefig(cwd / "TE_Images_ForPaper" / "Moments" / f"{2*key}-mer")
        fig_pd.savefig(cwd / "TE_Images_ForPaper" / "Moments" / f"{2*key}-mer_PD")


def moments(file_dir: pathlib.Path, ms: list, min_k: int = 1, max_k: int = 6, unlog: bool = False, convergence: bool = False, logy: bool = False, N_value: bool = False, title: str = None, *args, **kwargs):
    '''
    Moment calculations do NOT need the 4**s part multiplied, everything just needs to be scaled from 0 to 1 in a traditional normalization.*

    You have to feed in the list of moments, because I don't want to deal with a for float loop... I have that code somewhere I just can't find it. It's a for loop with floats, and it's pretty easy but there's a trick and I just don't want to do it... should be called the run on sentence bandit...


    Does a for loop of this. Each file needs to be a pickle file and have the same basic file name. Just give it the path and the filename and it will iterate over it.

    Exon files should be f'Exon_{k}mer.pkl' and Intron files should be f'Intron_{k}mer.pkl'

    The for loop for the k iterations is range(min_k, max_k + 1) so you do not need to artifically increase the max_k. That's done for you.

    *Do this as well with a N**(1/m - 1) term where N = 4**k (or in this case 4**2k) which should set all moment calculations for an evenly distributed heatmap to 1 (no matter what the moment is) and will
    adjust the other moment calculations to be either above or below 1. Sort of like subtracting the mean.
    '''

    if convergence and (max_k - min_k) < 1:
        print("Cannot continue, convergence can only be checked with more then 1 k value")
        return None 

    me = {}
    mi = {}
    pd = {}

    uni = {}

    for k in range(min_k, max_k + 1):
        print(f"K = {k}")

        me[k]: list = []
        mi[k]: list = []
        pd[k]: list = []
        uni[k]: list = []

        exon_file = file_dir / f"Exon_{k}mer.pkl"
        intron_file = file_dir / f"Intron_{k}mer.pkl"

        exon_data = hm._import_data(exon_file, just_import = True)
        exon_data = pandas.DataFrame(exon_data)
        if unlog:
            exon_data = _unrenormalize(exon_data, 2*k)
        print(f"\tImported Exon File")

        intron_data = hm._import_data(intron_file, just_import = True)
        intron_data = pandas.DataFrame(intron_data)
        if unlog:
            intron_data = _unrenormalize(intron_data, 2*k)
        print(f"\tImported Intron File")

        for m in ms:

            if N_value:
                N = 4**(2*k)
                N = N**((1/m) - 1)
                print(f"\t\tm = {m}\tN = {N}")
            else:
                N = 1
                print(f"\t\tm = {m}")

            exon_v = moment(exon_data, m = m, unlog = False, N = N, *args, **kwargs)
            # exon_u = uniform_density(exon_data, m = m, N = N)
            intron_v = moment(intron_data, m = m, unlog = False, N = N, *args, **kwargs)
            # intron_u = uniform_density(intron_data, m = m, N = N, *args, **kwargs)

            # if N_value:
            #     N = 4**(2*k)
            #     N = N**((1/m) - 1)
            #     exon_v = N * exon_v
            #     intron_v = N * intron_v

            if logy:
                exon_v = np.log2(exon_v)
                intron_v = np.log2(intron_v)
            
            try:
                per_diff = _percent_difference(exon_v, intron_v)
            except ZeroDivisionError:
                per_diff = 0
            
            me[k].append(exon_v)
            mi[k].append(intron_v)
            pd[k].append(per_diff)

            uni[k].append(1)

    print(f"Finished data, printing to plots")
    if convergence:
        fig, axs = plt.subplots(nrows = max_k - min_k, ncols = 4)

        keys = list(me.keys())
        # keys = keys[0: len(keys) - 2] # removing last entry, so that way I don't have to 

        for i in range(len(keys) - 1):

            axs[i, 0].plot(ms, mi[keys[i]], label = f"{2*keys[i]}-mer")
            axs[i, 0].plot(ms, mi[keys[i + 1]], label = f"{2*keys[i + 1]}-mer")
            axs[i, 0].legend()

            ipd = perdiff_list(mi[keys[i]], mi[keys[i + 1]])

            axs[i, 1].plot(ms, ipd)
            axs[i, 1].set_ylabel("PD")

            axs[i, 2].plot(ms, me[keys[i]], label = f"{2*keys[i]}-mer")
            axs[i, 2].plot(ms, me[keys[i + 1]], label = f"{2*keys[i + 1]}-mer")
            axs[i, 2].legend()

            epd = perdiff_list(me[keys[i]], me[keys[i + 1]])

            axs[i, 3].plot(ms, epd)
            axs[i, 3].set_ylabel("PD")

    else:
        fig, axs = plt.subplots(nrows = max_k - min_k + 1, ncols = 2)

        if len(me) > 1:
            for i, (key, value) in enumerate(me.items()):
                axs[i, 0].plot(ms, value, label = f"Exon")
                axs[i, 0].plot(ms, mi[key], label = f"Intron")
                axs[i, 0].plot(ms, uni[key], linestyle = "dotted")
                # axs[i, 0].plot(ms, ued[key], label = f"Exon U Density", linestyle = "dotted")
                axs[i, 0].set_ylabel(f"{2*key}-mer")
                axs[i, 0].legend()

                axs[i, 1].plot(ms, pd[key])
                axs[i, 1].set_ylabel("PD (0 to 1)")

        else:
            for key, value in me.items():
                axs[0].plot(ms, value, label = f"Exon")
                axs[0].plot(ms, mi[key], label = f"Intron")
                axs[0].plot(ms, uni[key], linestyle = "dotted")
                # axs[0].plot(ms, ued[key], label = f"Exon U Density", linestyle = "dotted")
                axs[0].set_ylabel(f"{2*key}-mer")
                axs[0].legend()

                axs[1].plot(ms, pd[key])
                axs[1].set_ylabel("Percent Difference")

    if isinstance(title, str):
        plt.suptitle(title)

    plt.show()


def perdiff_list(list1: list, list2: list):
    '''
    Because I don't want the moment script to be to long, I'm putting this part here.
    '''

    pd = []
    for i, ei in enumerate(list1):
        try:
            diff = _percent_difference(ei, list2[i])
        except ZeroDivisionError as e:
            diff = 0
        pd.append(diff)

    return pd



def _percent_difference(a, b):
    '''
    Because there's a possibility of a divide by zero warning, I have to make this its own function.
    '''
    demoninator = (a + b) / 2
    numerator = a - b

    if demoninator == 0 and numerator == 0:
        diff = 0
    elif demoninator == 0 and abs(numerator) > 0:
        diff = np.inf
    else:
        diff = numerator / demoninator
    
    return diff



def moment(data: pandas.DataFrame or dict or pathlib, m: float = 1, N: float = 1, unlog: bool = False, *args, **kwargs):
    '''
    Takes a region and finds the moment of that region.

    m = [sum(region)**n / sum(area)]**(1/n)

    Acutally I never did this next part so basically you can ignore it. Don't know if it would have been useful:
    It will have to import the whole data, then isolate a single area and count that. Maybe I shoudl count by region, add them, then calculate each moment.
    This would divide the space into R regions
    '''
    w: float = 1 / m  # think of this as an upside down m

    if isinstance(data, pathlib.Path) or isinstance(data, str):
        data = hm._import_data(data, just_import = True)

    if isinstance(data, dict):
        data: pandas.DataFrame = pandas.DataFrame(data)

    if unlog:
        data = _unrenormalize(data, *args, **kwargs)

    data = data.pow(m)
    values = data.to_numpy().sum()
    values = (1 / N) * (values**w)

    return values


def _unrenormalize(data: pandas.DataFrame, k: int, *args, **kwargs):
    '''
    Yes, that says un-re-normalize. This is because I saved only the log2(x + 1) data where x = 4**2 * (y/sum(y)). This undoes that
    because I'm lazy.

    log2( (4**s * x) + 1) = y
    4**s * x + 1 = 2**y
    4**s * x = 2**y - 1
    x = ((2**y) - 1) / 4**s

    s = 2*k -> I've always done this as k forward = k backward even though the scripts can handle the two being different.
    '''

    rows, cols = data.shape
    k = 4**k

    for row in range(rows):
        for col in range(cols):
            y = data.iloc[row, col]

            x = ((2**y) - 1) / k
            data.iloc[row, col] = x

    return data


def recover_data(dir: pathlib.Path, kmer: int = 6):
    '''
    Because I saved the data with a log2 transform with an extra term that totally screws with the moment calculations.

    This takes a file, saves the original file with a logtransformed extra in its filename, and unrenormalizes the data, and then saves that as a pkl file.

    I'm writing this as a general function in case I need to do it again.
    '''

    for k in range(1, kmer + 1):
        exon_file = dir / f"Exon_{k}mer.pkl"
        log_exon_file = dir / f"Exon_{k}mer_LT.pkl"
        intron_file = dir / f"Intron_{k}mer.pkl"
        log_intron_file = dir / f"Intron_{k}mer_LT.pkl"

        print(f"Importing from {exon_file}")
        exon_data = hm._import_data(exon_file, just_import = True)
        exon_data = pandas.DataFrame(exon_data)
        write_pickel(exon_data, log_exon_file)

        print(f"Importing from {exon_file}")
        intron_data = hm._import_data(intron_file, just_import = True)
        intron_data = pandas.DataFrame(intron_data)
        write_pickel(intron_data, log_intron_file)

        print(f"\tUnrenormalizing Exon Data")
        exon_data = _unrenormalize(exon_data, 2*k)
        write_pickel(exon_data, exon_file)

        print(f"\tUnrenormalizing Intron Data")
        intron_data = _unrenormalize(intron_data, 2*k)
        write_pickel(intron_data, intron_file)


def write_pickel(file_data: pandas.DataFrame, file_path: pathlib.Path):
    '''
    '''
    print(f"Exporting to {file_path}")
    with open(file_path, "wb") as file:
        file_data.to_pickle(file)


def uniform_density(data: pandas.DataFrame, m: int = 1, *args, **kwargs):
    '''
    Finds the moment for a uniform density for comparison to the rest of the moment plots.
    '''

    rows, cols = data.shape

    ave = data.to_numpy().mean()

    data = pandas.DataFrame(0, index = np.arange(rows), columns = np.arange(cols))
    data = data + ave

    uoment = moment(data, m = m, *args, **kwargs)

    return uoment



if __name__ in "__main__":
    main()