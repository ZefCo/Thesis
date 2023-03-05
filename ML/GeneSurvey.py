import pandas
import pathlib
import re
# from tensorflow.keras.models import Sequential, load_module
# from tensorflow.keras.layers import Conv1D
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import accuracy_score
# import tensorflow as tf
import random
import numpy as np
import pandas as pd
# import keras
# import Permutations
import itertools
import re
import plotly.graph_objects as go
import copy

cwd = pathlib.Path.cwd()
data = np.load("TrainingGeneData_v3.pkl", allow_pickle=True)

print(data)

data["Type"] = pandas.Categorical(data["Type"])
seq_types = pandas.unique(data["Type"])


type_counts = dict()
for t in seq_types:
    type_rows = data[data["Type"] == t].shape[0]
    # print(f"Type = {t}\t{type_rows}")
    type_counts[t] = type_rows

# Find how many gene parts there are
# Find which chromosomes they're taken from
    # Check how the set was built, it's probably there
# Find breakdown of all possible nucleotide permutations (of 2)
    # Exon, Intron, UTR5, UTR3

def nucleotide_permutations() -> list:

    nucs = ["A", "C", "G", "T"]
    nucz = nucs.copy()

    nuc_perm = list()

    for s in nucs:
        for z in nucz:
            nuc_perm.append(f"{s}{z}")

    return nuc_perm


def count_occurences(n_data, permutations: list):
    '''
    '''
    occurences = dict()

    for perm in permutations:
        counts = len(re.findall(perm, n_data))

        occurences[perm] = counts

    return occurences


nuc_perm = nucleotide_permutations()
master_occurences = {perm:0 for perm in nuc_perm}
exon_occurences = copy.deepcopy(master_occurences)
intron_occurences = copy.deepcopy(master_occurences)
utr5_occurences = copy.deepcopy(master_occurences)
utr3_occurences = copy.deepcopy(master_occurences)
missing_no = copy.deepcopy(master_occurences)

master_norm = 0
exon_norm, intron_norm, utr5_norm, utr3_norm = master_norm, master_norm, master_norm, master_norm

# for perm in nuc_perm:
#     counts = len(re.findall(perm, "ACACACACACACACACACAC"))

#     occurences[perm] = counts

# print(occurences)

rows, _ = data.shape

for row in range(rows):
    row_of_interest = data.iloc[row, :].copy()
    sequences = row_of_interest["Seq"].upper()

    local_occurances = count_occurences(sequences, nuc_perm)

    for p, c in local_occurances.items():
        master_occurences[p] += c
        master_norm += c

    if row_of_interest["Type"] in "Intron":
        for p, c in local_occurances.items():
            intron_occurences[p] += c
            intron_norm += c
    if row_of_interest["Type"] in "CDS":
        for p, c in local_occurances.items():
            exon_occurences[p] += c
            exon_norm += c
    if row_of_interest["Type"] in "UTR5":
        for p, c in local_occurances.items():
            utr5_occurences[p] += c
            utr5_norm += c
    if row_of_interest["Type"] in "UTR3":
        for p, c in local_occurances.items():
            utr3_occurences[p] += c
            utr3_norm += c


for p in master_occurences.keys():
    master_occurences[p] = master_occurences[p] / master_norm
    exon_occurences[p] = exon_occurences[p] / exon_norm
    intron_occurences[p] = intron_occurences[p] / intron_norm
    utr5_occurences[p] = utr5_occurences[p] / utr5_norm
    utr3_occurences[p] = utr3_occurences[p] / utr3_norm


# print(master_occurences)
# master_occurences = pandas.DataFrame(master_occurences)
# print(master_occurences)
# print(exon_occurences)
# print(intron_occurences)
# print(utr5_occurences)
# print(utr3_occurences)
# print(type_counts)

hist = go.Figure()

hist.add_trace(go.Bar(x = list(master_occurences.keys()), y = list(master_occurences.values()), name = "All"))
hist.add_trace(go.Bar(x = list(exon_occurences.keys()), y = list(exon_occurences.values()), name = "CDS"))
hist.add_trace(go.Bar(x = list(intron_occurences.keys()), y = list(intron_occurences.values()), name = "Intron"))
hist.add_trace(go.Bar(x = list(utr5_occurences.keys()), y = list(utr5_occurences.values()), name = "UTR5"))
hist.add_trace(go.Bar(x = list(utr3_occurences.keys()), y = list(utr3_occurences.values()), name = "UTR3"))

# hist.update_traces(opacity = 0.5)
hist.update_layout(title = "All occurances")

hist.show()