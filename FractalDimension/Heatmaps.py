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
import GeneClass as Gene
import random
from icecream import ic


def main():
    '''
    Create a heat map: X = past, Y = future, Z = count. So count the occurences of a k-mer - to k-mer +.
    '''
    linux_path = f"/media/ethanspeakman/Elements/"
    windows_path = f"G:/"

    data_path = windows_path
    data_set = 2
    output_files = [cwd / "TE_Images_ForPaper" / "Dict" / f"Heat_DS{data_set}_Master.pkl", cwd / "TE_Images_ForPaper" / "Dict" / f"Heat_DS{data_set}_Exon.pkl", cwd / "TE_Images_ForPaper" / "Dict" / f"Heat_DS{data_set}_Intron.pkl"]
    
    master, _, _ = heat_embedding(pathlib.Path(f"{data_path}/Gene_Data_Sets/Data_Set_{data_set}_histogram.pkl"),
                   n = 1_000,
                   k_m = 1, k_p = 1)
                #    dict_output_files = output_files)

    heatmap(master, label_points=True)


def heatmap(data: dict,
            label_points: bool = False):
    '''
    This feels like it's doing this backwards and it's driving me nuts.
    '''
    if isinstance(data, pathlib.Path) or isinstance(data, str):
        data = _import_data(data, just_import = True)

    data["A_p"] = data.pop("A")

    data: pandas.DataFrame = pandas.DataFrame(data).T
    data = data.reindex(sorted(data.columns), axis = 1)
    data = data.sort_index()
    print(data)
    rows, cols = data.shape
    
    fig, ax = plt.subplots()
    im = ax.imshow(data)
    # print(list(data.index))

    ax.set_xticks(np.arange(rows), labels = list(data.index))
    ax.set_yticks(np.arange(cols), labels = list(data.columns))
    plt.setp(ax.get_xticklabels(), ha = "right", rotation_mode = "anchor")

    if label_points:
        for i in range(rows):
            for j in range(cols):
                text = ax.text(j, i, round(data.iloc[i, j], 2), ha = "center", va = "center", color = "w")

    plt.colorbar(im)
    plt.xlabel("Back")
    plt.ylabel("Forward")
    plt.title("Heatmap of Back vs Forward")
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

    master_dict = dict()
    exon_dict = dict()
    intron_dict = dict()

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

    master_dict = _normalize(master_dict, int(4**(k_p + k_m)))
    exon_dict = _normalize(exon_dict, int(4**(k_p + k_m)))
    intron_dict = _normalize(intron_dict, int(4**(k_p + k_m)))

    ic(master_dict)

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
        return master_dict, exon_dict, intron_dict


def _dict_deep_merge(dict_to_merge: dict, dict_to_add: dict):
    '''
    dict_to_add copies the data into dict_to_merge
    '''
    # ic(dict_to_add)
    forward_dict: dict

    for back_key, forward_dict in dict_to_add.items():
        # key does not exist in master_key, therefore dict does not exist
        # key does exist in master_key
            # sub_key not in master_sub_dict, therefore no value
            # sub_key does exist in master_sub_dict, update value
        
        if back_key not in dict_to_merge.keys():
            dict_to_merge[back_key] = dict()
            for forward_key, forward_value in forward_dict.items():
                dict_to_merge[back_key][forward_key] = forward_value

        else:
            for forward_key, forward_value in forward_dict.items():
                if forward_key not in dict_to_merge[back_key].keys():
                    dict_to_merge[back_key][forward_key] = forward_value
                else:
                    dict_to_merge[back_key][forward_key] += forward_value

    return dict_to_merge


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
        
        if k_prime not in local_dict.keys():
            local_dict[k_prime] = {}
            local_dict[k_prime][k_plus[i]] = 1
        else:
            if k_plus[i] not in local_dict[k_prime].keys():
                local_dict[k_prime][k_plus[i]] = 1
            else:
                local_dict[k_prime][k_plus[i]] += 1

    return local_dict


def _import_data(data: pandas.DataFrame,
                 n: int = 50_000,
                 k_p: int = 6, k_m: int = 6,
                 just_import: bool = False,
                 *args, **kwargs):
    '''
    Generates the heat map data. Feeds in a pandas Dataframe (or pathlib to import the data), then outputs the overall data.
    '''
    if isinstance(data, pathlib.Path):
        with open(data, "rb") as p:
            data = pickle.load(p)

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


def _normalize(data: dict, k: int):
    '''
    Normalizes the data by: counting and summing all values for a normalization constant, then doing: value = 4**(kp + km) * value / constant

    4**(kp + km) = k. I don't want to figure it out here so... you do it. It's a necessary input.
    '''
    norm_constant = 0
    for sub_dict in data.values():
        for value in sub_dict.values():
            norm_constant += value

    for master_key, sub_dict in data.items():
        for key, value in sub_dict.items():
            data[master_key][key] = (k * value) / norm_constant

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