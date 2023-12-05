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
    ns = [n / step for n in range(int(step*(min_n)), int(step*max_n) + 1)]
    title = f"E v I log"
    moments(cwd / "Dicts", ns, min_k = 1, max_k = 6, convergence = False, logy = True, N_value = False, title = title)
    # evalue = moment(exon_file, unlog = True, k = 2*k, n = 1)
    # ivalue = moment(intron_file, unlog = True, k = 2*k, n = 1)
    # print(f"k = {k}\te value = {evalue}\ti = value = {ivalue}")



def moments(file_dir: pathlib.Path, ns: list, min_k: int = 1, max_k: int = 6, unlog: bool = False, convergence: bool = False, logy: bool = False, N_value: bool = False, title: str = None, *args, **kwargs):
    '''
    Moment calculations do NOT need the 4**s part multiplied, everything just needs to be scaled from 0 to 1 in a traditional normalization.

    You have to feed in the list of n, because I don't want to deal with a for float loop... I have that code somewhere I just can't find it.


    Does a for loop of this. Each file needs to be a pickle file and have the same basic file name. Just give it the path and the filename and it will iterate over it.

    Exon files should be f'Exon_{k}mer.pkl' and Intron files should be f'Intron_{k}mer.pkl'

    The for loop for the k iterations is range(min_k, max_k + 1) so you do not need to artifically increase the max_k. That's done for you.

    Do this as well with a N**(1/m - 1) term where N = 4**k (or in this case 4**2k) which should set all moment calculations for an evenly distributed heatmap to 1 (no matter what the moment is) and will
    adjust the other moment calculations to be either above or below 1. Sort of like subtracting the mean.
    '''

    if convergence and (max_k - min_k) < 1:
        print("Cannot continue, convergence can only be checked with more then 1 k value")
        return None 

    ne = {}
    ni = {}
    pd = {}

    for k in range(min_k, max_k + 1):
        print(f"K = {k}")

        ne[k]: list = []
        ni[k]: list = []
        pd[k]: list = []

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

        for n in ns:

            print(f"\t\tn = {n}")
            
            exon_v = moment(exon_data, n = n, k = 2*k, unlog = False, *args, **kwargs)
            intron_v = moment(intron_data, n = n, k = 2*k, unlog = False, *args, **kwargs)

            if N_value:
                N = 4**(2*k)
                N = N**((1/n) - 1)
                exon_v = N * exon_v
                intron_v = N * intron_v

            if logy:
                exon_v = np.log2(exon_v)
                intron_v = np.log2(intron_v)
            
            try:
                per_diff = _percent_difference(exon_v, intron_v)
            except ZeroDivisionError:
                per_diff = 0
            
            ne[k].append(exon_v)
            ni[k].append(intron_v)
            pd[k].append(per_diff)

    print(f"Finished data, printing to plots")
    if convergence:
        fig, axs = plt.subplots(nrows = max_k - min_k, ncols = 4)

        keys = list(ne.keys())
        # keys = keys[0: len(keys) - 2] # removing last entry, so that way I don't have to 

        for i in range(len(keys) - 1):

            axs[i, 0].plot(ns, ni[keys[i]], label = f"{2*keys[i]}-mer")
            axs[i, 0].plot(ns, ni[keys[i + 1]], label = f"{2*keys[i + 1]}-mer")
            axs[i, 0].legend()

            ipd = perdiff_list(ni[keys[i]], ni[keys[i + 1]])

            axs[i, 1].plot(ns, ipd)
            axs[i, 1].set_ylabel("PD")

            axs[i, 2].plot(ns, ne[keys[i]], label = f"{2*keys[i]}-mer")
            axs[i, 2].plot(ns, ne[keys[i + 1]], label = f"{2*keys[i + 1]}-mer")
            axs[i, 2].legend()

            epd = perdiff_list(ne[keys[i]], ne[keys[i + 1]])

            axs[i, 3].plot(ns, epd)
            axs[i, 3].set_ylabel("PD")

    else:
        fig, axs = plt.subplots(nrows = max_k - min_k + 1, ncols = 2)

        if len(ne) > 1:
            for i, (key, value) in enumerate(ne.items()):
                axs[i, 0].plot(ns, value, label = f"Exon")
                axs[i, 0].plot(ns, ni[key], label = f"Intron")
                axs[i, 0].set_ylabel(f"{2*key}-mer")
                axs[i, 0].legend()

                axs[i, 1].plot(ns, pd[key])
                axs[i, 1].set_ylabel("PD (0 to 1)")

        else:
            for key, value in ne.items():
                axs[0].plot(ns, value, label = f"Exon")
                axs[0].plot(ns, ni[key], label = f"Intron")
                axs[0].set_ylabel(f"{2*key}-mer")
                axs[0].legend()

                axs[1].plot(ns, pd[key])
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



def moment(data: pandas.DataFrame or dict or pathlib, n: float = 1, regions = None, unlog: bool = False, *args, **kwargs):
    '''
    Takes a region and finds the moment of that region.

    m = [sum(region)**n / sum(area)]**(1/n)

    It will have to import the whole data, then isolate a single area and count that. Maybe I shoudl count by region, add them, then calculate each moment.
    This would divide the space into R regions
    '''
    u: float = 1 / n  # think of this as an upside down n

    if isinstance(data, pathlib.Path) or isinstance(data, str):
        data = hm._import_data(data, just_import = True)

    if isinstance(data, dict):
        data: pandas.DataFrame = pandas.DataFrame(data)

    if unlog:
        data = _unrenormalize(data, *args, **kwargs)

    data = data.pow(n)
    values = data.to_numpy().sum()
    values = values**u

    return values


def _unrenormalize(data: pandas.DataFrame, k: int, *args, **kwargs):
    '''
    Yes, that says un-re-normalize. This is because I saved only the log2(x + 1) data where x = 4**2 * (y/sum(y)). This undoes that
    because I'm lazy.

    log2( (4**s * x) + 1) = y
    4**s * x + 1 = 2**y
    4**s * x = 2**y - 1
    x = ((2**y) - 1) / 4**s
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




if __name__ in "__main__":
    main()