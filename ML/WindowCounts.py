import pandas
import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import re

mer = 6

global_data: pandas.DataFrame = pandas.read_pickle(cwd / "TrainingGeneData_v5.pkl")
# print(global_data.shape)

keep = np.where(global_data["Seq"].str.len() >= 100)[0]
# print(keep)
global_data = global_data.iloc[keep, :]

rows, cols = global_data.shape

global_data = global_data.reset_index()

# print(global_data["Type"].unique())

# print(global_data.shape)
# print(global_data)

print(global_data)

# intron, intname = 0, "Intron"
# exon, exname = 0, "Exon"

# global_mer = dict()
# exon_mer = dict()
# intron_mer = dict()

# exon_data = global_data[global_data[""]]

# for row in range(rows):
#     row_of_interest = global_data.iloc[row, :]

