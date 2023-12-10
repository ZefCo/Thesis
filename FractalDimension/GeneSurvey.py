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
    kmer = 3
    linux_path = f"/media/ethanspeakman/Elements/"
    windows_path = f"G:/"

    data_path = windows_path

    # This needs to use the same data for the histogram
    data_1 = pathlib.Path(f"{data_path}/Gene_Data_Sets/Data_Set_1_histogram.pkl")  # I know the method says it was a dataframe, but I also coded it where it can just take a pathlib and load the data. I got lazy
    output_file_1 = cwd / f"GeneSurvey_1_{kmer}mer.pkl"
    data_2 = pathlib.Path(f"{data_path}/Gene_Data_Sets/Data_Set_2_histogram.pkl")  # I know the method says it was a dataframe, but I also coded it where it can just take a pathlib and load the data. I got lazy
    output_file_2 = cwd / f"GeneSurvey_2_{kmer}mer.pkl"

    # print(data_1)
    # survey(data_1, output_file_1, kmer, reverse = True)
    # survey(data_2, output_file_2, kmer, reverse = True)
    recreate(output_file_1, kmer, title = "Dataset 1 Histogram Method Seqeunces", labels = True, output_file = pathlib.Path(cwd / f"GeneSurvey_Dataset1_{kmer}mer.html"))
    recreate(output_file_2, kmer, title = "Dataset 2 Histogram Method Seqeunces", labels = True, output_file = pathlib.Path(cwd / f"GeneSurvey_Dataset2_{kmer}mer.html"))


def count_occurences(sequence: str, permutations: list):
    '''
    Takes a string that represents the sequence and outputs how many times a given occurance 
    '''
    occurences = dict()

    for perm in permutations:
        counts = len(re.findall(perm, sequence))

        occurences[perm] = counts

    return occurences


def nucleotide_counter(sequence: str, window_size: int, reverse: bool = False) -> dict:
    '''
    '''
    keys: set = set()
    counter = dict()
    master_count = 0


    for i in range(len(sequence) - window_size):
        seq = sequence[i: i + window_size].upper()
        if reverse:
            seq = seq[::-1]

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


def regional_survey(data: pandas.DataFrame, kmer: int, sequence_name: str = "Seq", *args, **kwargs) -> dict:
    '''
    '''
    rows, cols = data.shape
    regional_mer = dict()

    for row in range(rows):
        seq_of_interest = data.loc[row, sequence_name]

        nuc_count = nucleotide_counter(seq_of_interest, kmer, **kwargs)

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
           sequence_column: str = "Seq", classification_name: str = "Classificaion", exon_name: str = "exon", intron_name: str = "intron",
           *args, **kwargs):
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
    print(counts_frame)

    counts_frame.to_pickle(output_file, )


def recreate(filepath: pathlib.Path, kmer, output_file: pathlib.Path = None, 
             shaded = False, labels = False, title = None, 
             read_type = "pkl"):
    '''
    creates the plots. It's called recreate as a holdover from when I first did this, which is the WindowCounts.py script.

    Note the read_type: that can be pkl or csv right now. By default it's pkl
    '''
    if read_type in "pkl":
        data: pandas.DataFrame = pandas.read_pickle(filepath)
        x_column = data.index
    elif read_type in "csv":
        data: pandas.DataFrame = pandas.read_csv(filepath, header = 0)
        x_column = data["Unnamed: 0"]
    # print(data)

    average = 0.25**kmer

    fig = go.Figure()
    fig.add_trace(go.Bar(x = x_column, y = data["G%"], name = "Full Seq"))
    fig.add_trace(go.Bar(x = x_column, y = data["E%"], name = "Exon"))
    fig.add_trace(go.Bar(x = x_column, y = data["I%"], name = "Intron"))

    # fig = make_subplots(rows = 3, cols = 1)
    # fig.add_trace(go.Bar(x = list(data.index), y = data["G%"], name = "Global"), row = 1, col = 1)
    # fig.add_trace(go.Bar(x = list(data.index), y = data["E%"], name = "Exon"), row = 2, col = 1)
    # fig.add_trace(go.Bar(x = list(data.index), y = data["I%"], name = "Intron"), row = 3, col = 1)

    if shaded:
        if kmer > 1:
            sections = (4**kmer)/4
            fig.add_vrect(x0 =            - 0.5, x1 =   sections - 0.5, col = "all", fillcolor = "red",    opacity = 0.25, annotation_text = "A")
            fig.add_vrect(x0 =   sections - 0.5, x1 = 2*sections - 0.5, col = "all", fillcolor = "green",  opacity = 0.25, annotation_text = "C")
            fig.add_vrect(x0 = 2*sections - 0.5, x1 = 3*sections - 0.5, col = "all", fillcolor = "blue",   opacity = 0.25, annotation_text = "G")
            fig.add_vrect(x0 = 3*sections - 0.5, x1 = 4*sections - 0.5, col = "all", fillcolor = "yellow", opacity = 0.25, annotation_text = "T")

    if title is None:
        title = f"Occurences of {kmer}-Mer window Sequences <br> Subtracted by a mean of {average}"
    else:
        title = f"{title}<br>Occurences of {kmer}-Mer window Sequences <br> Subtracted by a mean of {average}"
    fig.update_layout(barmode = "overlay", title = title)
    fig.update_traces(opacity = 0.9)
    fig.update_xaxes(showticklabels = labels)
    fig.show()
    if output_file is not None:
        print(f"Output to \t{output_file}")
        fig.write_html(str(output_file))



if __name__ in "__main__":
    main()