import pandas
import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import re
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import itertools


def main():
    '''
    '''
    histogram()

def sameLength(threshold: float = 0.15):
    '''
    '''
    percent_diff = lambda x1, x2: (abs(x1 - x2)) / (0.5*(x1 + x2))

    data = cwd / "TrainingGeneData_v5.pkl"
    new_data = pandas.DataFrame()

    data: pandas.DataFrame = pandas.read_pickle(data)


    data["Length"] = data["Seq"].str.len()
    data = data[data["Length"] > 100]

    exon_subset = data[data["Type"] == "Exon"]
    exon_subset = exon_subset.reset_index(drop = True)
    intron_subset = data[data["Type"] == "Intron"]
    intron_subset = intron_subset.reset_index(drop = True)
    print(intron_subset)

    eows, eols = exon_subset.shape

    for eow in range(eows):
        target_length = exon_subset.loc[eow, "Length"]

        intron_keep = np.where(percent_diff(target_length, intron_subset["Length"]) <= threshold)[0]
        intron_subsub = intron_subset.iloc[intron_keep, :]
        # print(intron_subsub)

        if (intron_subsub.shape[0] > 0):

            rand_samp = intron_subsub.sample(n = 1)

            new_data = pandas.concat([new_data, exon_subset.iloc[eow, :].to_frame().T, rand_samp], ignore_index = True).reset_index(drop = True)

            intron_subset = intron_subset.drop([rand_samp.index[0]])
            intron_subset = intron_subset.reset_index(drop = True)

            # intron_subset = intron_subset.drop([])
            # exit()


        # for iow in range(intron_subset.shape[0]):
        #     intron_length = intron_subset.loc[iow, "Length"]
        #     if percent_diff(target_length, intron_length) <= threshold:
        #         new_data = pandas.concat([new_data, exon_subset.iloc[eow, :].to_frame().T, intron_subset.iloc[iow, :].to_frame().T], ignore_index=True)
        #         new_data = new_data.reset_index(drop = True)

        #         # print(new_data)

        #         intron_subset = intron_subset.drop([iow])
        #         intron_subset = intron_subset.reset_index(drop = True)
        #         break

        if (eow % 1000) == 0:
            print(f"Completed row {eow}")

    # # new_data = new_data.drop(["index"])
    # new_data["Type"] = pandas.Categorical(new_data["Type"])
    # print(new_data["Type"].unique())
    # # print(data)
    # print(new_data[new_data["Type"] == "Exon"])
    # print(new_data[new_data["Type"] == "Intron"])
    new_data.to_pickle(cwd / "TrainingData_SameSize_11.pkl")


def histogram(data_file: pathlib.Path, min_length: int = 100, max_length: int = 10_000, classification_col: str = "Type"):
    '''
    '''
    # data = cwd / "TrainingGeneData_v5.pkl"
    new_data = pandas.DataFrame()

    data: pandas.DataFrame = pandas.read_pickle(data_file)


    data["Length"] = data["Seq"].str.len()
    # print(data)
    # data = data[data["Length"] > 100]

    data = data[(data["Length"] >= min_length) & (data["Length"] <= max_length)]


    # data.to_pickle(cwd / "TrainingData_SameSize_His.pkl")

    exon_subset = data[(data[classification_col] == "Exon") | (data[classification_col] == "exon") | (data[classification_col] == "e")]
    exon_subset = exon_subset.reset_index(drop = True)
    intron_subset = data[(data[classification_col] == "Intron") | (data[classification_col] == "intron") | (data[classification_col] == "i")]
    intron_subset = intron_subset.reset_index(drop = True)
    # print(intron_subset)

    print(exon_subset["Length"].mean())
    print(intron_subset["Length"].mean())

    exon_subset["LogLength"] = np.log10(exon_subset["Length"])
    intron_subset["LogLength"] = np.log10(intron_subset["Length"])

    # exon_len_values = exon_subset["Length"].unique()
    # print(min(exon_len_values), max(exon_len_values))
    # intron_len_values = intron_subset["Length"].unique()
    # print(min(intron_len_values), max(intron_len_values))

    figRaw = go.Figure()
    figRaw.add_trace(go.Histogram(x = intron_subset["Length"], name = "Intron", nbinsx = len(intron_subset["Length"].unique())))
    figRaw.add_trace(go.Histogram(x = exon_subset["Length"], name = "Exon", nbinsx = len(intron_subset["Length"].unique())))
    figRaw.update_layout(barmode = "overlay")
    figRaw.show()

    figLog = go.Figure()
    figLog.add_trace(go.Histogram(x = intron_subset["LogLength"], name = "Intron", nbinsx = len(intron_subset["LogLength"].unique())))
    figLog.add_trace(go.Histogram(x = exon_subset["LogLength"], name = "Exon", nbinsx = len(intron_subset["LogLength"].unique())))
    figLog.update_layout(barmode = "overlay")
    figLog.show()






if __name__ in "__main__":
    main()
