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
import pandas
import pickle
import GeneClass as Gene
import TimeEmbedding as TE


def main():
    '''
    To whoever reads this next: I'm sorry. Things got out of hand. I've got folders and folders of scripts and things are diverging and it's hard to keep track of what is supposed to be what.
    It made sense at the start, and now I'm loosing track of what is going on.

    This script is meant to scan through all the different regions and take note of what they are composed of.
    '''
    kmer = 2

    data = pathlib.Path("G:\Data_Set_1.pkl")  # I know the method says it was a dataframe, but I also coded it where it can just take a pathlib and load the data. I got lazy
    output_file = cwd / "GeneSurvey.pkl"
    survey(data, output_file, kmer)


def count_occurences(sequence: str, permutations: list):
    '''
    Takes a string that represents the sequence and outputs how many times a given occurance 
    '''
    occurences = dict()

    for perm in permutations:
        counts = len(re.findall(perm, sequence))

        occurences[perm] = counts

    return occurences


def nucleotide_counter(sequence: str, window_size: int) -> dict:
    '''
    '''
    keys: set = set()
    counter = dict()
    master_count = 0


    for i in range(len(sequence) - window_size):
        seq = sequence[i: i + window_size].upper()

        if seq not in keys:
            keys.add(seq)
            counter[seq] = 1
            master_count += 1

        else:
            counter[seq] += 1
            master_count += 1

    # for key, value in counter.items():
    #     counter[key] = value / master_count

    return counter


def regional_survey(data: pandas.DataFrame, kmer: int, sequence_name: str = "Seq") -> dict:
    '''
    '''
    rows, cols = data.shape
    regional_mer = dict()

    for row in range(rows):
        seq_of_interest = data.loc[row, sequence_name]

        nuc_count = nucleotide_counter(seq_of_interest, kmer)

        for key, value in nuc_count.items():
            if key in regional_mer.keys():
                regional_mer[key] += value
            else:
                regional_mer[key] = value

    return regional_mer



def statistics(global_mer, exon_mer, intron_mer, kmer: int) -> pandas.DataFrame:
    '''
    '''
    average = 0.25**kmer

    counts_frame = pandas.DataFrame([global_mer, exon_mer, intron_mer]).T
    counts_frame.columns = ["Global", "Exon", "Intron"]
    counts_frame = counts_frame.fillna(0)

    motiffs = list(counts_frame.index)
    for motiff in motiffs:
        if re.search("N", motiff):
            counts_frame = counts_frame.drop(index = (motiff))

    count_global_sum = counts_frame["Global"].sum()
    count_exon_sum = counts_frame["Exon"].sum()
    count_intron_sum = counts_frame["Intron"].sum()

    counts_frame["G%"] = (counts_frame["Global"] / count_global_sum) - average
    counts_frame["E%"] = (counts_frame["Exon"] / count_exon_sum) - average
    counts_frame["I%"] = (counts_frame["Intron"] / count_intron_sum) - average

    counts_frame = counts_frame.sort_index()

    return counts_frame


def survey(data: pandas.DataFrame or pathlib.Path, output_file: pathlib.Path, kmer: int, min_length: int = 10, 
           sequence_column: str = "Seq", classification_name: str = "Classificaion", exon_name: str = "exon", intron_name: str = "intron"):
    '''
    '''
    if isinstance(data, pathlib.Path):
        with open(data, "rb") as p:
            data = pickle.load(p)

    keep = np.where(data[sequence_column].str.len() >= min_length)[0]
    data = data.iloc[keep, :]
    data = data.reset_index()

    exon_data = data[data[classification_name] == exon_name]
    exon_data = exon_data.reset_index()
    intron_data = data[data[classification_name] == intron_name]
    intron_data = intron_data.reset_index()

    del data  # deletes the data frame to free up memory. It's going to be big, and while my home computer has 48 Gigs, the school computer only has 16

    exon_mer: dict = regional_survey(exon_data, kmer)
    intron_mer: dict = regional_survey(intron_data, kmer)

    global_mer: dict = exon_mer

    for key, value in intron_mer.items():
        if key in global_mer.keys():
            global_mer[key] += value
        else:
            global_mer[key] = value

    counts_frame = statistics(global_mer, exon_mer, intron_mer, kmer)

    counts_frame.to_pickle(output_file)


if __name__ in "__main__":
    main()