import pandas
import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import re

kmer = 6

def nucleotide_counter(sequence: str, window_size: int):
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

global_data: pandas.DataFrame = pandas.read_pickle(cwd / "TrainingGeneData_v5.pkl")
# print(global_data.shape)

keep = np.where(global_data["Seq"].str.len() >= 100)[0]
# print(keep)
global_data = global_data.iloc[keep, :]

# rows, cols = global_data.shape

global_data = global_data.reset_index()

# print(global_data["Type"].unique())

# print(global_data.shape)
# print(global_data)

# print(global_data)

intron, intname = 0, "Intron"
exon, exname = 0, "Exon"

global_mer = dict()
exon_mer = dict()
intron_mer = dict()

exon_data = global_data[global_data["Type"] == exname]
exon_data = exon_data.reset_index()
eows, cols = exon_data.shape

print(f"Starting Exon run\nRows = {eows}")
for row in range(eows):
    seq_of_interest = exon_data.loc[row, "Seq"]
    nuc_count = nucleotide_counter(seq_of_interest, kmer)

    for key, value in nuc_count.items():
        if key in exon_mer.keys():
            exon_mer[key] += value
            global_mer[key] += value
        else:
            exon_mer[key] = value
            global_mer[key] = value

    # if row < 10:
        # break


intron_data = global_data[global_data["Type"] == intname]
intron_data = intron_data.reset_index()
iows, cols = intron_data.shape

print(f"Starting Intron run\nRows = {iows}")
for row in range(iows):
    seq_of_interest = intron_data.loc[row, "Seq"]
    nuc_count = nucleotide_counter(seq_of_interest, kmer)

    for key, value in nuc_count.items():
        if key in intron_mer.keys():
            intron_mer[key] += value
        else:
            intron_mer[key] = value

        if key in global_mer.keys():
            global_mer[key] += value
        else:
            global_mer[key] = value

    # if row < 10:
        # break

counts_frame = pandas.DataFrame([global_mer, exon_mer, intron_mer]).T
counts_frame.columns = ["Global", "Exon", "Intron"]
counts_frame = counts_frame.fillna(0)
print(counts_frame)
counts_frame.to_pickle(cwd / "KMerCounts.pkl")

# print(exon_mer)
# print(intron_mer)
# print(global_mer)
