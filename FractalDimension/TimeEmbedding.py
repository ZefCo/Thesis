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
from statsmodels.graphics.tsaplots import plot_acf
import pandas
import pickle
import itertools

def main():
    '''
    '''
    time_embedding(max_rows = 5000)
    # score_keys()



def score_keys(k = 9):
    '''
    '''
    perms = tuple(itertools.product(["A", "C", "G", "T"], repeat = k))

    w_p = [(0.25)**n for n in range(1, k + 1)]
    w_m = [(0.25)**n for n in range(1, k + 1)]
    w_m.reverse()

    scores = {}
    for p in perms:
        key = ""
        for n in p:
            key = f"{key}{n}"

        scores[key] = {"history": 0, "future": 0}

        digital = [0 if n in "A" else 1 if n in "C" else 2 if n in "G" else 3 if n in "T" else 100 for n in key]

        score_p = np.dot(w_p, digital)
        score_m = np.dot(w_m, digital)

        scores[key]["future"] = score_p
        scores[key]["history"] = score_m

    scores = pandas.DataFrame(data = scores).T
    print(scores)

    scores.to_csv(cwd / "ScoreKey.csv")



def time_embedding(k_p = 9, k_m = 9, max_rows = 200):
    '''
    Will use the G number idea to find the forward and backward number.

    A = 0
    C = 1
    G = 2
    T = 3

    N = ?? <- I hope there are no Ns

    load the sequence data, take a window of 6, find the numbers, plot the numbers

    Switch from using the Training Data to a single gene from the Fractal Di mension stuff.
    '''

    w_p = [(0.25)**n for n in range(1, k_p + 1)]
    w_m = [(0.25)**n for n in range(1, k_m + 1)]
    w_m.reverse()

    ML_Folder = cwd.parent / "ML"
    
    with open(str(ML_Folder / "TrainingData_SameSize.pkl"), "rb") as p:
        data: pandas.DataFrame = pickle.load(p)

    fig = go.Figure()
    eig = go.Figure()
    iig = go.Figure()
    e_count, i_count = 0, 0

    # print(data.loc[4, "Seq"])
    # exit()

    rows, cols = data.shape

    data = data[data["Length"] > (k_m + k_p)]

    for row in range(rows):
        sequence = data.loc[row, "Seq"].upper()
        length = data.loc[row, "Length"]
        region = data.loc[row, "Type"]
        x, y = [], []
        # ex, ey = [], []
        # ix, iy = [], []

        k_minus = [sequence[k_prime:k_prime + k_m] for k_prime in range(0, length - (k_p + k_m))]
        k_plus = [sequence[k_prime:k_prime + k_p] for k_prime in range(k_m, length - k_p)]

        for i, k_prime in enumerate(k_minus):
            n = [0 if n in "A" else 1 if n in "C" else 2 if n in "G" else 3 if n in "T" else 100 for n in k_prime]
            k_x = np.dot(w_m, n)

            n = [0 if n in "A" else 1 if n in "C" else 2 if n in "G" else 3 if n in "T" else 100 for n in k_plus[i]]
            k_y = np.dot(w_p, n)

            x.append(k_x), y.append(k_y)

        if region in "Exon":
            name = f"Exon_{e_count}"
            legend_group = f"exon{e_count}"

            eig.add_trace(go.Scatter(x = x, y = y, name = name, legendgroup = legend_group, mode = "markers", marker = dict(size = 1)))

            e_count += 1
        else:
            name = f"Intron_{i_count}"
            legend_group = f"intron{i_count}"

            iig.add_trace(go.Scatter(x = x, y = y, name = name, legendgroup = legend_group, mode = "markers", marker = dict(size = 1)))

            i_count += 1


        fig.add_trace(go.Scatter(x = x, y = y, name = name, legendgroup = legend_group, mode = "markers", marker = dict(size = 1)))
        # fig.add_trace(go.Scatter(x = x, y = y, name = name, legendgroup = legend_group, mode = "lines", line = dict(width = 0.25)))

        if row > max_rows:
            break

    fig.update_layout(title = f"Exons and Introns, Time Embedding<br>{e_count + i_count} Total Regions - Unknown how many genes are represented", xaxis_title = "History", yaxis_title = "Future")
    # fig.show()
    fig.update_layout(showlegend = False)
    fig.write_image(str(cwd / "ExonIntron_HistoryvsFuture.png"))

    eig.update_layout(title = f"Exons Only<br>{e_count} total Exons", xaxis_title = "History", yaxis_title = "Future")
    # eig.show()
    eig.update_layout(showlegend = False)
    eig.write_image(str(cwd / "Exon_HistoryvsFuture.png"))

    iig.update_layout(title = f"Introns Only<br>{i_count} total Introns", xaxis_title = "History", yaxis_title = "Future")
    # iig.show()
    iig.update_layout(showlegend = False)
    iig.write_image(str(cwd / "Intron_HistoryvsFuture.png"))

        # exit()



def scoring(N: list):
    '''
    '''
    n = [0 if n in "A" else 1 if n in "C" else 2 if n in "G" else 3 if n in "T" else 100 for n in N]
    w = [(0.25)**n for n in range(len(N))]

    print(n)
    print(w)

    return np.dot(np.array(w), np.array(n))



if __name__ in "__main__":
    main()