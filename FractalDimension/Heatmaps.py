import pathlib
cwd = pathlib.Path.cwd()
import os
import re
import numpy as np
import glob
from plotly import graph_objects as go
import timeit
import DistanceClass as distance
import pandas
import pickle
import itertools
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import seaborn as sns
import GeneClass as Gene
import random
from icecream import ic
from typing import Callable


def main():
    '''
    Create a heat map: X = past, Y = future, Z = count. So count the occurences of a k-mer - to k-mer +.
    '''

    linux_path = f"/media/ethanspeakman/Elements/"
    windows_path = f"F:/"
    data_path = windows_path
    data_set = 2
    n = 100_000
    transform = np.log2
    colors_map = "hot"
    nucsequence: str = "AGTC"
    cbins = 10
    # cticks = [0, 0.6, 0.65, 0.75, 0.8, 0.85, 0.9, 0.91, 0.91, 0.91, 0.91, 1.2, 1.3, 1.4]
    # title_details = f"Independent Scaling Log2 Transform {nucsequence} Order Cbins = {cbins}"
    # file_details = title_details.replace(" ", "")

    # k = 1
    # colors: list = ['darkred', 'red', 'orange', 'gold', 'yellow'] 
    # bounds: list = [0, 0.8, 1.0, 1.1, 2.0]
    # # bounds: list = [0, 0.8, 1.0, 1.1, 2.0]
    # # bounds: list = [0, 0.8, 1.0, 1.1, 1.2]

    # k = 2
    # colors: list = ['darkred', 'firebrick', 'red', 'orange', 'gold', 'yellow']
    # # bounds: list = [0, 0.65, 0.87, 0.97, 1.1, 1.2] # 0.92
    # # bounds: list = [0, 0.65, 0.87, 0.97, 1.1, 1.38]  # 0.95
    # bounds: list = [0, 0.65, 0.87, 0.97, 1.1, 1.49]  # 0.98

    # k = 3
    # colors: list = ['darkred', 'firebrick', 'red', 'orange', 'gold', 'yellow']
    # # bounds: list = [0, 0.45, 0.85, 0.95, 1.2, 1.75]  # 0.98
    # # bounds: list = [0, 0.45, 0.85, 0.95, 1.2, 1.5]  # 0.95
    # bounds: list = [0, 0.45, 0.85, 0.95, 1.2, 1.4]  # 0.92

    k = 4
    colors: list = ['darkred', 'firebrick', 'red', 'orange', 'gold', 'yellow']
    # bounds: list = [0, 0.35, 0.65, 0.95, 1.2, 2.0]  # 0.98
    # bounds: list = [0, 0.35, 0.65, 0.95, 1.2, 1.76]  # 0.95
    bounds: list = [0, 0.35, 0.65, 0.95, 1.2, 1.56]  # 0.92

    k = 5
    colors: list = ['darkred', 'firebrick', 'red', 'orange', 'gold', 'yellow']
    bounds: list = [0, 0.15, 0.45, 0.85, 1.26, 3.5]
    bounds: list = [0, 0.15, 0.45, 0.85, 1.26, 3.5]
    bounds: list = [0, 0.15, 0.45, 0.85, 1.26, 3.5]

    # k = 6
    # colors: list = ['black', 'darkred', 'firebrick', 'red', 'orange', 'gold', 'yellow']
    # bounds: list = [0, 0.754, 1.056, 1.457, 1.859, 2.462, 3.165]
    # bounds: list = [0, 0.754, 1.056, 1.457, 1.859, 2.462, 3.165]
    # bounds: list = [0, 0.754, 1.056, 1.457, 1.859, 2.462, 3.165]

    if 1 == k:
        label_points = True
        label_axis = True
    elif ((3 > k) and (k > 1)):
        label_points = False
        label_axis = True
    else:
        label_points = False
        label_axis = False

    exon_file = cwd / "Dicts" / f"Exon_{k}mer.pkl"
    intron_file = cwd / "Dicts" / f"Intron_{k}mer.pkl"

    output_files = [cwd / "TE_Images_ForPaper" / "Dict" / f"Heat_DS{data_set}_Master.pkl", cwd / "TE_Images_ForPaper" / "Dict" / f"Heat_DS{data_set}_Exon.pkl", cwd / "TE_Images_ForPaper" / "Dict" / f"Heat_DS{data_set}_Intron.pkl"]

    # _, exon, intron, _, _, exon_max, exon_min, intron_max, intron_min = heat_embedding(cwd / "TE_Images_ForPaper" / "Dict" / "Seq_For_Images_n100000_minLength12.pkl",
    #                                                                                    n = n,
    #                                                                                    k_m = k, k_p = k,
    #                                                                                    just_import = True,
    #                                                                                    log_transform = transform, 
    #                                                                                    nucsequence = nucsequence)  # Just import is turned on because I have one dictionary with all the data. I'll use that one set for all the images from now on. Just makes things a little faster and easier.

    exon = _import_data(exon_file, just_import=True)
    intron = _import_data(intron_file, just_import=True)
    intron_max, intron_min = min_max(intron)

    # cticks = [((intron_max - 0) / cbins)*i for i in range(cbins)]
    # title_details = f"Independent Scaling Log2 Transform {nucsequence} Order Cbins = {len(bounds)}"
    title_details = f"Same Area Heatmap bounds v3"
    file_details = title_details.replace(" ", "")


    # pdf_cdf(exon_file, intron_file,
    #         # print_targets = [0.7, 0.8, 0.85, 0.9, 0.95, 0.98], 
    #         print_targets = [0.2, 0.4, 0.6, 0.8, 0.85, 0.9, 0.92, 0.95, 0.98], 
    #         pdf_title = f"PDF for {2*k}mer", pdf_outpath = cwd / "TE_Images_ForPaper" / "Histograms" / f"PDF_{2*k}mer.png",
    #         cdf_title = f"CDF for {2*k}mer", cdf_outpath = cwd / "TE_Images_ForPaper" / "Histograms" / f"CDF_{2*k}mer.png")
    # exit()

    # # master_max, master_min = min_max(master)
    # print(f"Master Max = {minmax[0]}\tMaster Min = {minmax[1]}")
    # # exon_max, exon_min = min_max(exon)
    # print(f"Exon Max = {minmax[2]}\Exon Min = {minmax[3]}")
    # # intron_max, intron_min = min_max(intron)
    # print(f"Intron Max = {minmax[4]}\tIntron Min = {minmax[5]}")

    # These values came out of the heat_embedding script above for k = 6 (12mer)
    # master_rescale = [7758, 0.25]
    # exon_rescale = [3102, 0.84]
    # intron_rescale = [9792, 0.37]
    # minmax = [7758, 0, 3102, 0, 9792, 0]

    # k_max = intron_max if intron_max > exon_max else exon_max
    # k_min = intron_min if intron_min < exon_min else exon_min

    # # master = _rescale_data(master, minmax[0], minmax[1])
    # exon = _rescale_data(exon, k_max, k_min)
    # intron = _rescale_data(intron, k_max, k_min)

    # heatmap(master, 
    #         label_points = label_points, label_axis = label_axis,
    #         title = f"Master Heatmap {2*k}mer") #vmin = master_rescale[0], vmax = master_rescale[1]) #, file_output = cwd / "TE_Images_ForPaper" / "Heatmaps" / f"Master_DS{data_set}_n{n}_{k*2}mer.png")
    
    heatmapv2(exon, colors = colors, bounds = bounds, title = f"Exon Equal Area Heatmap\n{2*k}mer", fileoutput = cwd / "TE_Images_ForPaper" / "Heatmaps" / "EA" / f"Exon_HEA_DS{data_set}_n{n}_{k*2}mer_{file_details}.png")
    heatmapv2(intron, colors = colors, bounds = bounds, title = f"Intron Equal Area Heatmap\n{2*k}mer", fileoutput = cwd / "TE_Images_ForPaper" / "Heatmaps" / "EA" / f"Intron_HEA_DS{data_set}_n{n}_{k*2}mer_{file_details}.png")
    exit()


    # heatmap(exon,
    #         label_points = label_points, label_axis = label_axis,
    #         title = f"Exon Heatmap {2*k}mer\n{title_details}",
    #         nucsequence = nucsequence,
    #         vmin = intron_min, vmax = intron_max, cticks = cticks,
    #         colors_map = colors_map,
    #         file_output = cwd / "TE_Images_ForPaper" / "Heatmaps" / f"Exon_DS{data_set}_n{n}_{k*2}mer_{file_details}.png") #vmin = exon_rescale[0], vmax = exon_rescale[1]) #, file_output = cwd / "TE_Images_ForPaper" / "Heatmaps" / f"Exon_DS{data_set}_n{n}_{k*2}mer.png")
    # heatmap(intron,
    #         label_points = label_points, label_axis = label_axis,
    #         title = f"Intron Heatmap {2*k}mer\n{title_details}",
    #         nucsequence = nucsequence,
    #         vmin = intron_min, vmax = intron_max, cticks = cticks,
    #         colors_map = colors_map,
    #         file_output = cwd / "TE_Images_ForPaper" / "Heatmaps" / f"Intron_DS{data_set}_n{n}_{k*2}mer_{file_details}.png") #vmin = intron_rescale[0], vmax = intron_rescale[1]) #, file_output = cwd / "TE_Images_ForPaper" / "Heatmaps" / f"Intron_DS{data_set}_n{n}_{k*2}mer.png")


def pdf_cdf(data1: dict or pathlib.Path, data2: dict or pandas.DataFrame, 
            pdf_title: str = None, cdf_title: str = None, 
            pdf_outpath: pathlib.Path = None, cdf_outpath: pathlib.Path = None,
            print_targets: list = None,
            *args, **kwargs):
    '''
    Does the PDF and CDF of the data all at once.

    It's just going to be easier to do it like this.
    '''
    if isinstance(data1, pathlib.Path):
        data1 = _import_data(data1, just_import = True)
    if isinstance(data2, pathlib.Path):
        data2 = _import_data(data2, just_import = True)

    data1 = _count_data(data1)
    data2 = _count_data(data2)

    data12 = {}
    data12 = _dict_deep_copy(data12, data1, data2)

    max_key = max(list(data12.keys()))
    bins = int((max_key - 0)*10)
    delta = (max_key - 0) / bins

    pdf, cdf, xticks = _init_pdf_cdf(bins, delta)

    for x, c in data12.items():
        for i in range(bins):
            if (pdf[i, 1] <= x) and (x < pdf[i, 2]):
                pdf[i, 3] += c
                break
        else:
            pdf[i, 3] += c

    cumm = 0
    for i in range(bins):
        cumm += pdf[i, 3]
        cdf[i, 3] = cumm

    if isinstance(print_targets, list):
        for target in print_targets:
            idx = np.asarray(cdf[:, 3])
            idx = (np.abs(idx - target)).argmin()

            print(f"Closest to {target}: {xticks[idx]}")

    bar_chart(pdf, xticks = xticks, title = pdf_title, outpath = pdf_outpath)
    bar_chart(cdf, xticks = xticks, title = cdf_title, outpath = cdf_outpath)


def bar_chart(data: np.array, xticks: list = None, title: str = None, outpath: pathlib.Path = None, *args, **kwargs):
    '''
    Does the PDF/CDF chart for input data.
    '''
    fig, ax = plt.subplots()
    ax.bar(data[:, 0], data[:, 3])

    if isinstance(xticks, list):
        ax.set_xticks(data[:, 0], xticks, rotation = 90)
    
    if isinstance(title, str):
        plt.title(title)

    if isinstance(outpath, pathlib.Path):
        plt.savefig(outpath)
        plt.close()
        print(f"Output image to\n\t{outpath}")
    else:
        plt.show()
    # exit()



def _init_pdf_cdf(bins: int, delta: float):
    '''
    Numpy arrays just feel a little better for this then a dictionary.

    column 0 = an index
    column 1 = x min for the bin
    column 2 = x max for the bin
    column 3 = midpoint between x min and x max - for placement of the actual bin, mostly the label. This also is returned from the function as a list of strs.
    column 4 = bin count (normalized preferably, right now this will be zero)

    Actually returns two: one for the pdf and one for the cdf.

    OK this looks a lot like a pandas Dataframe... sue me.
    '''
    pdf = np.zeros(shape = (bins, 4))
    cdf = np.zeros(shape = (bins, 4))
    xticks = []
    for i in range(bins):
        pdf[i, 0], cdf[i, 0] = i, i

        pdf[i, 1], pdf[i, 2] = (i)*delta, (i + 1)*delta
        cdf[i, 1], cdf[i, 2] = (i)*delta, (i + 1)*delta

        xticks.append(str(round((pdf[i, 1] + pdf[i, 2])/2, 3)))

    return pdf, cdf, xticks



def _count_data(data: pandas.DataFrame or dict, *args, **kwargs):
    '''
    '''
    if isinstance(data, dict):
        data: pandas.DataFrame = pandas.DataFrame(data)
    data = _reorder_frame(data, *args, **kwargs)
    rows, cols = data.shape

    hist_data = {}

    for row in range(rows):
        for col in range(cols):
            key = data.iloc[row, col]
            if key not in hist_data.keys():
                hist_data[key] = 1
            else:
                hist_data[key] += 1

    return hist_data


def heatmapv2(data: dict or pandas.DataFrame or pathlib.Path,
              colors: list = ['darkred', 'red', 'orange', 'gold', 'yellow'], 
              bounds: list = [0, 0.8, 1.0, 1.1, 2.0],
              title: str = None,
              fileoutput: pathlib.Path = None,
              *args, **kwargs):
    '''
    Uses seaborn and does an 'equal area' heatmap.
    '''
    if isinstance(data, pathlib.Path) or isinstance(data, str):
        data = _import_data(data, just_import = True)

    if isinstance(data, dict):
        data: pandas.DataFrame = pandas.DataFrame(data)

    data = _reorder_frame(data, *args, **kwargs)


    cus_cmap = ListedColormap(colors)
    cus_norm = BoundaryNorm(bounds, ncolors = len(bounds))

    fig, ax = plt.subplots()
    sns.heatmap(data,
                ax = ax,
                cmap = cus_cmap,
                norm = cus_norm)
    
    plt.xlabel("Back")
    plt.ylabel("Forward")
    if isinstance(title, str):
        plt.title(title)
    else:
        plt.title("Heatmap of Equal area")

    if isinstance(fileoutput, pathlib.Path):
        plt.savefig(fileoutput)
        plt.close
        print(f"Output file to\n\t{fileoutput}")
    else:
        plt.show()



def heatmap(data: dict or pandas.DataFrame or pathlib.Path,
            label_points: bool = False,
            label_axis: bool = True,
            file_output: str or pathlib.Path = None,
            title: str = None,
            vmin: float = None, vmax: float = None, cbins: int = None,
            colors_map: str = "gray",
            cticks: list = None,
            *args, **kwargs):
    '''
    Creates a heatmap. The input can either be a pathlib (in which case it will import the data), a dictionary (in which case it will convert it to a dataframe), or a dataframe (in which case it will create a heatmap).
    '''
    if isinstance(data, pathlib.Path) or isinstance(data, str):
        data = _import_data(data, just_import = True)

    if isinstance(data, dict):
        data: pandas.DataFrame = pandas.DataFrame(data)

    # data = data.fillna(0)
    # print("###Entering Reorder frame###")
    # print(data)
    data = _reorder_frame(data, *args, **kwargs)
    # print("###Returning from Reorder Frame###")
    # print(data)
    # exit()
    # # data = data.sort_index()
    # row_names, col_names = list(data.index), list(data.columns)
    # print(row_names)
    # print(type(row_names))
    # exit()
    rows, cols = data.shape
    
    fig, ax = plt.subplots()
    if (isinstance(vmin, float) and isinstance(vmax, float) and (isinstance(cticks, list))):
        cmap = plt.cm.get_cmap(colors_map, len(cticks))
        print(type(cmap))
        im = ax.imshow(data, cmap = cmap, vmin = vmin, vmax = vmax)
    else:
        im = ax.imshow(data)
    # print(list(data.index))
    
    if label_axis:
        ax.set_xticks(np.arange(rows), labels = list(data.columns))
        ax.set_yticks(np.arange(cols), labels = list(data.index))
        plt.setp(ax.get_xticklabels(), ha = "right", rotation_mode = "anchor")

    if label_points:
        for i in range(rows):
            for j in range(cols):
                text = ax.text(j, i, round(data.iloc[i, j], 2), ha = "center", va = "center", color = "w")

    if (isinstance(vmin, float) and isinstance(vmax, float) and (isinstance(cticks, list))):
        plt.colorbar(im, ticks = cticks)
    else:
        plt.colorbar(im)
    plt.xlabel("Back")
    plt.ylabel("Forward")
    if isinstance(title, str):
        plt.title(title)
    else:
        plt.title("Heatmap of Back vs Forward")

    if isinstance(file_output, str) or isinstance(file_output, pathlib.Path):
        print(f"Outputing file to\n\t{file_output}")
        plt.savefig(file_output)
    else:   
        plt.show()


def heat_embedding(data: pandas.DataFrame,
                   n: int = 50_000,
                   k_p: int = 6, k_m: int = 6, 
                   backwards: bool = True, # I'm leaving this in here in case I want to use it later. Probably not though
                   nucsequence: str = "AGTC", sequence_name: str = "Seq", classification_name: str = "Classificaion",
                   dict_output_files: list = None,
                   *args, **kwargs):
    '''
    Somewhat similar to the time embedding, this looks at the past and the future and creates a heatmap.

    This script will look at the previous nucleotides and counts the occurrences of the next nucleotides, returning a dictionary in the form:
    return_dict[past_nucleotide] = sub_dict[future_nucleotide]: occurrences

    If the dict_output_files is a list then the paths/strs inside will be used to save the dictionaries to a pickel file. The order of them is [master, exon, intron]. Sorry, I'm not going to think of a way to predict what order you want the files to be saved as, so that's the order.
    If it's left as None, the dictionaries are returned as a tuple.
    '''
    if isinstance(data, str) or isinstance(data, pathlib.Path):
        data = _import_data(data, n, k_p, k_m, *args, **kwargs)

    rows, _ = data.shape

    master_dict, exon_dict, intron_dict = _init_dict(k_p, nucsequence = nucsequence)
    # ic(exon_dict)

    for row in range(rows):
        try:
            sequence = data.loc[row, sequence_name].upper()
        except Exception as e:
            print(type(e))
            print(e)
            print(row)
            continue

        region = data.loc[row, classification_name]

        local_dict = _heat_data(sequence, k_p = k_p, k_m = k_m, nucsequence = nucsequence)
        master_dict = _dict_deep_merge(master_dict, local_dict)

        if (region in "Exon") or (region in "exon") or (region in "e"):  # because I realized I was doing this like 3 different ways... I probably should have been more precise.
            exon_dict = _dict_deep_merge(exon_dict, local_dict)

        elif (region in "Intron") or (region in "intron") or (region in "i"):
            intron_dict = _dict_deep_merge(intron_dict, local_dict)

    master_dict, master_max, master_min = _normalize(master_dict, int(4**(k_p + k_m)), *args, **kwargs)
    exon_dict, exon_max, exon_min = _normalize(exon_dict, int(4**(k_p + k_m)), *args, **kwargs)
    intron_dict, intron_max, intron_min = _normalize(intron_dict, int(4**(k_p + k_m)), *args, **kwargs)

    # if k_m + k_p < 5:
    #     ic(master_dict)

    if isinstance(dict_output_files, list):
        try:
            master_filepath: str or pathlib.Path = dict_output_files[0]
        except IndexError as i:
            print("Cannot output master dict")
        else:
            _write_pickle(master_dict, master_filepath)

        try:
            exon_filepath: str or pathlib.Path = dict_output_files[1]
        except IndexError as i:
            print("Cannot output exon dict")
        else:
            _write_pickle(exon_dict, exon_filepath)

        try:
            intron_filepath: str or pathlib.Path = dict_output_files[2]
        except IndexError as i:
            print("Cannot output intron dict")
        else:
            _write_pickle(intron_dict, intron_filepath)

    else:
        return master_dict, exon_dict, intron_dict, master_max, master_min, exon_max, exon_min, intron_max, intron_min
    

def _reorder_frame(dataframe: pandas.DataFrame, *args, **kwargs):
    '''
    finds the current index and columns, then reorders them in an alternative alphabet i.e. ACGT or AGCT or AGTC.

    This is probably redundent since I'm using Perms earlier to set the order, but I'm going to keep this here in case I want to change the order during the script.
    '''
    rows = list(dataframe.index)
    new_rows = {row: _digitize_seq(row, *args, **kwargs) for row in rows}

    # print(dataframe)

    dataframe = dataframe.rename(index = new_rows)
    dataframe = dataframe.rename(columns = new_rows)
    # print("###Renamed###")
    # print(dataframe)

    # print(dataframe)
    # print(list(dataframe.columns)[0], type(list(dataframe.columns)[0]))

    # sorted_rows = sorted(list(new_rows.values()))
    # columns = list(dataframe.columns)
    # print(columns[0], type(columns[0]))
    # print(sorted_rows, type(sorted_rows))

    # dataframe = dataframe[[sorted_rows]]
    dataframe = dataframe.reindex(sorted(dataframe.columns), axis = 1)
    # dataframe = dataframe.reindex(sorted(dataframe.columns).reverse(), axis = 0)
    dataframe = dataframe.sort_index(ascending = False)
    # print("###Reindexed###")
    # print(dataframe)
    # print(dataframe)

    return dataframe




def _digitize_seq(seq: str, nucseq = "AGTC", *args, **kwargs):
    '''
    Converts a sequence into a digial version - OK not digital but a numerical version, but I like this name better.
    '''
    digital_version = ""
    for n in seq:
        dn = "0" if n in nucseq[0] else "1" if n in nucseq[1] else "2" if n in nucseq[2] else "3" if n in nucseq[3] else "5"
        digital_version = f"{digital_version}{dn}"

    return digital_version


def _dict_deep_merge(dict_to_merge: dict, dict_to_add: dict):
    '''
    dict_to_add copies the data into dict_to_merge

    This assumes that there are sub dictionaries in the merging dictionaries.
    '''
    # ic(dict_to_add)
    forward_dict: dict

    for back_key, forward_dict in dict_to_add.items():
        for forward_key, forward_value in forward_dict.items():
            dict_to_merge[back_key][forward_key] += forward_value

    return dict_to_merge


def _dict_deep_copy(dictionary: dict, *args, **kwargs):
    '''
    Merges keys and values of multiple dictionaries into the main one.

    Also normalizes them.
    '''
    d: dict
    norm_value = 0
    # print(args)

    for d in args:
        for key, value in d.items():
            if key not in dictionary.keys():
                dictionary[key] = value
            else:
                dictionary[key] += value

            norm_value += value

    for key, value in dictionary.items():
        dictionary[key] = value / norm_value

    return dictionary


def _dict_normalize(dictionary: dict):
    '''
    Normalizes the values in the dictionary. No sub dictionaries.
    '''
    norm_value = 0
    for value in dictionary.values():
        norm_value += value

    for key, value in dictionary.items():
        dictionary[key] = value / norm_value

    return dictionary



def _heat_data(sequence: str, 
              k_p: int = 6, k_m: int = 6, 
              nucsequence: str = "AGTC") -> dict:
    '''
    This part does the actual counting of the occurrences, and will return a local dictionary from teh sequence.
    '''
    local_dict = dict()

    sequence = sequence.upper()
    nucsequence = nucsequence.upper()
    seq_length = len(sequence)

    if seq_length < (k_m + k_p):
        print("Cannont find Trajectory for this gene: to small")
        return None

    k_minus = [sequence[k_prime:k_prime + k_m] for k_prime in range(0, seq_length - (k_p + k_m))]
    k_plus = [sequence[k_prime:k_prime + k_p] for k_prime in range(k_m, seq_length - k_p)]

    for i, k_prime in enumerate(k_minus):
        k_prime = k_prime[::-1]
        
        if k_prime not in local_dict.keys():
            local_dict[k_prime] = {}
            local_dict[k_prime][k_plus[i]] = 1
        else:
            if k_plus[i] not in local_dict[k_prime].keys():
                local_dict[k_prime][k_plus[i]] = 1
            else:
                local_dict[k_prime][k_plus[i]] += 1

    return local_dict


def _import_data(data: pathlib.Path,
                 n: int = 50_000,
                 k_p: int = 6, k_m: int = 6,
                 just_import: bool = False,
                 *args, **kwargs):
    '''
    Generates the heat map data. Feeds in a pandas Dataframe (or pathlib to import the data), then outputs the overall data.
    '''
    try:
        if isinstance(data, pathlib.Path):
            with open(data, "rb") as p:
                data = pickle.load(p)
    except FileNotFoundError:
        print(f"File was not found, did you mean to try to import this file?\n{data}")
        exit()

    if just_import:
        pass
    
    else:
        try:
            data = data.sample(n = n).reset_index()
        except ValueError as e:
            data = data.reset_index()
        except Exception as e:
            print(type(e))
            print(e)
            print(data.shape)
            print(n)
            exit()

        if "Length" in data.columns:
            data = data[data["Length"] > (k_m + k_p)]
        else:
            data["Length"] = data.Seq.str.len()
            data = data[data["Length"] > (k_m + k_p)]

    return data


def _init_dict(k: int, nucsequence = "AGTC"):
    '''
    '''
    master, exon, intron = dict(), dict(), dict()
    perms = tuple(itertools.product([nucsequence[0], nucsequence[1], nucsequence[2], nucsequence[3]], repeat = k))
    for p in perms:
        p = "".join(p)
        master[p], exon[p], intron[p] = dict(), dict(), dict()
        for pp in perms:
            pp = "".join(pp)
            master[p][pp], exon[p][pp], intron[p][pp] = 0, 0, 0

    return master, exon, intron



def min_max(data: dict):
    '''
    Finds the minimum and maximum values in the dictionaries. Because they're embedded dictionaries this will have to use some embedded for loops.
    Even though Dictionaries are fast I'm doing both at once since I don't want to iterate across the same datasets twice.
    '''
    max_value, min_value = None, None
    for _, sub_dict in data.items():
        for _, _ in sub_dict.items():

            if max_value is None:
                max_key: str = max(sub_dict, key = sub_dict.get)
                max_value: float = sub_dict[max_key]

                min_key: str = min(sub_dict, key = sub_dict.get)
                min_value: float = sub_dict[min_key]

            else:
                max_key: str = max(sub_dict, key = sub_dict.get)
                max_local: float = sub_dict[max_key]
                if max_local > max_value:
                    max_value = max_local

                min_key: str = min(sub_dict, key = sub_dict.get)
                min_local: float = sub_dict[min_key]
                if min_value > min_local:
                    min_value = min_local

    return max_value, min_value



def _normalize(data: dict, k: int, norm_constant: int = None, log_transform: Callable = None, *args, **kwargs):
    '''
    Normalizes the data by: counting and summing all values for a normalization constant, then doing: value = 4**(kp + km) * value / constant

    4**(kp + km) = k. I don't want to figure it out here so... you do it. It's a necessary input.

    Add in a part to get the min/max out of here. Find it as you move through the dict.
    '''
    max_value, min_value = 0, None
    if isinstance(norm_constant, int) or isinstance(norm_constant, float):
        pass
    else:
        norm_constant = 0
        for sub_dict in data.values():
            for value in sub_dict.values():
                norm_constant += value

    for master_key, sub_dict in data.items():
        for key, value in sub_dict.items():
            new_value = (k * value) / norm_constant
            if isinstance(log_transform, np.ufunc):
                new_value = log_transform(new_value + 1)
            data[master_key][key] = new_value

            if new_value > max_value:
                max_value = new_value

            if min_value is None:
                min_value = new_value
            elif min_value > new_value:
                min_value = new_value

    return data, max_value, min_value


def _rescale_data(data: dict, rescale_max: float, rescale_min: float):
    '''
    '''
    old_max, _ = min_max(data)
    M = ((rescale_max - rescale_min) / old_max)

    for master_key, sub_dict in data.items():
        for key, value in sub_dict.items():

            rescale_value = M * value + rescale_min
            data[master_key][key] = rescale_value

    return data



def _write_pickle(data, filepath: str):
    '''
    Writes a pickle file
    '''
    try:
        with open(filepath, "wb") as file:
            print(f"Writing to {filepath}")
            pickle.dump(data, file)
    except Exception as e:
        print("Unable to save to {filepath}")
        print("Reason:")
        print(f"\t{type(e)}")
        print(e)


if __name__ in "__main__":
    main()