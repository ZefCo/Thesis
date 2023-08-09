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

def main():
    '''
    I really need a new name for this, Time Embedding doens't work.

    k-Mer Time Series. That's a better name.

    Time embedding v3: make one that goes through a gene (all its different forms) and plots the trajectory of the k windows. Does something happen for the introns and the exons?
    You'll have to also do the exons and introns seperatly, but we want to see how the trajectory can "jump" from exon to intron: maybe there is something of interest there?
    Time embedding v3: make one that goes through a gene (all its different forms) and plots the trajectory of the k windows. Does something happen for the introns and the exons?
    You'll have to also do the exons and introns seperatly, but we want to see how the trajectory can "jump" from exon to intron: maybe there is something of interest there?
    '''
    time_embedding_v2(str(cwd.parent / "ML" / "TrainingData_SameSize.pkl"), 
                      max_rows = 5000, 
                      gap = 0, 
                      k_p = 6, 
                      k_m = 6, 
                      backwards = True,
                      compliment = False)
    # time_embedding_v2(str(cwd.parent / "ML" / "TrainingData_SameSize.pkl"), 
    #                   max_rows = 5000, 
    #                   gap = 0, 
    #                   k_p = 6, 
    #                   k_m = 6, 
    #                   backwards = False,
    #                   compliment = False)
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



def time_embedding(k_p = 9, k_m = 9, gap = 0, max_rows = 200):
    '''
    Will use the G number idea to find the forward and backward number.

    A = 0
    C = 1
    G = 2
    T = 3

    N = ?? <- I hope there are no Ns

    load the sequence data, take a window of 6, find the numbers, plot the numbers

    Switch from using the Training Data to a single gene from the Fractal Di mension stuff.

    This uses plotly, which is great, but has issues with large number of points.
    '''
    image_dir = cwd / "TE_Images"
    image_dir.mkdir(parents = True, exist_ok = True)
    exon_dir = image_dir / "Exon"
    exon_dir.mkdir(parents = True, exist_ok = True)
    intron_dir = image_dir / "Intron"
    intron_dir.mkdir(parents = True, exist_ok = True)
    both_dir = image_dir / "Both"
    both_dir.mkdir(parents = True, exist_ok = True)

    w_p = [(0.25)**n for n in range(1, k_p + 1)]
    w_m = [(0.25)**n for n in range(1, k_m + 1)]
    w_m.reverse()

    ML_Folder = cwd.parent / "ML"
    
    with open(str(ML_Folder / "TrainingGeneData_SLSGHis.pkl"), "rb") as p:
        data: pandas.DataFrame = pickle.load(p)

    fig = go.Figure()
    eig = go.Figure()
    iig = go.Figure()
    e_count, i_count = 0, 0

    # print(data.loc[4, "Seq"])
    # exit()


    if "Length" in data.columns:
        data = data[data["Length"] > (k_m + k_p + gap)]
    else:
        data["Length"] = data.Seq.str.len()
        data = data[data["Length"] > (k_m + k_p + gap)]

    data = data.sample(n = max_rows).reset_index()
    rows, cols = data.shape

    # print(data)


    for row in range(rows):
        try:
            sequence = data.loc[row, "Seq"].upper()
        except Exception as e:
            print(type(e))
            print(e)
            print(row)
            exit()

        length = data.loc[row, "Length"]

        region = data.loc[row, "Type"]
        x, y = [], []
        # ex, ey = [], []
        # ix, iy = [], []

        k_minus = [sequence[k_prime:k_prime + k_m] for k_prime in range(0, length - (k_p + k_m + gap))]
        k_plus = [sequence[k_prime:k_prime + k_p] for k_prime in range(gap + k_m, length - k_p)]

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

        if ((row % 1000) == 0):
            print(f"Finished row {row}")

        # if row > max_rows:
        #     break

    fig.update_layout(title = f"Exons and Introns, Time Embedding w/ {gap}-mer Gap between + and -<br>{e_count + i_count} Total Regions - Unknown how many genes are represented", xaxis_title = "History, -", yaxis_title = "Future, +")
    # fig.show()
    fig.update_layout(showlegend = False, autosize = False, width = 1500, height = 1500)
    fig.write_image(str(both_dir / f"both_gap_{gap}.png"))

    eig.update_layout(title = f"Exons Only w/ {gap}-mer Gap<br>{e_count} total Exons", xaxis_title = "History, -", yaxis_title = "Future, +")
    # eig.show()
    eig.update_layout(showlegend = False, autosize = False, width = 1500, height = 1500)
    eig.write_image(str(exon_dir / f"exon_gap_{gap}.png"))

    iig.update_layout(title = f"Introns Only w/ {gap}-mer Gap<br>{i_count} total Introns", xaxis_title = "History, -", yaxis_title = "Future, +")
    # iig.show()
    iig.update_layout(showlegend = False, autosize = False, width = 1500, height = 1500)
    iig.write_image(str(intron_dir / f"intron_gap_{gap}.png"))

        # exit()




def time_embedding_v2(pickle_file, k_p = 9, k_m = 9, gap = 0, max_rows = 200, backwards = True, compliment = False):
    '''
    Almost identical to the above, but uses matplotlib to output the images. Doesn't look as slick but does output the images much faster.

    References Time Embedding v3 so I can make this more consistent with various changes I make.
    '''
    image_dir = cwd / "TE_Images"
    image_dir.mkdir(parents = True, exist_ok = True)
    exon_dir = image_dir / "Exon"
    exon_dir.mkdir(parents = True, exist_ok = True)
    intron_dir = image_dir / "Intron"
    intron_dir.mkdir(parents = True, exist_ok = True)
    both_dir = image_dir / "Both"
    both_dir.mkdir(parents = True, exist_ok = True)

    with open(str(pickle_file), "rb") as p:
        data: pandas.DataFrame = pickle.load(p)

    if "Length" in data.columns:
        data = data[data["Length"] > (k_m + k_p + gap)]
    else:
        data["Length"] = data.Seq.str.len()
        data = data[data["Length"] > (k_m + k_p + gap)]

    try:
        data = data.sample(n = max_rows).reset_index()
    except ValueError as e:
        data = data.reset_index()
    except Exception as e:
        print(type(e))
        print(e)
        print(data.shape)
        print(max_rows)
        exit()
    rows, cols = data.shape

    # print(data)

    # b_frame, e_frame, i_frame = np.zeros(shape = (rows, 1, 2)), np.zeros(shape = (data[data["Type"] == "Exon"].shape[0], 1, 2)), np.zeros(shape = (data[data["Type"] == "Intron"].shape[0], 1, 2))
    b_frame, e_frame, i_frame = [], [], []
    e_count, i_count = 0, 0

    for row in range(rows):
        try:
            sequence = data.loc[row, "Seq"].upper()
        except Exception as e:
            print(type(e))
            print(e)
            print(row)
            exit()

        region = data.loc[row, "Type"]

        xy = time_embedding_v3(sequence, k_p = k_p, k_m = k_m, gap = gap, m_backwards = backwards, compliment = compliment)
        
        if region in "Exon":
            e_frame.append(xy)
            e_count += 1
        else:
            i_frame.append(xy)
            i_count += 1

        b_frame.append(xy)

        if ((row % 1000) == 0):
            print(f"Finished row {row}")

    # Both Plot
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 20)
    for points in b_frame:
        ax.scatter(points[:, 0], points[:, 1], s = 0.1)

    if backwards:
        b_title = f"Exons and Introns, Time Embedding w/ {gap}-mer Gap between + and -\n{e_count + i_count} Total Regions - Unknown how many genes are represented: weights are forwards and backwards\nCompliment = {compliment}"
        e_title = f"Exons, Time Embedding w/ {gap}-mer Gap\n{e_count} Total Regions: weights are forwards and backwards\nCompliment = {compliment}"
        i_title = f"Introns, Time Embedding w/ {gap}-mer Gap\n{i_count} Total Regions: weights are forwards and backwards\nCompliment = {compliment}"
    else:
        b_title = f"Exons and Introns, Time Embedding w/ {gap}-mer Gap between + and -\n{e_count + i_count} Total Regions - Unknown how many genes are represented: weights are forwards\nCompliment = {compliment}"
        e_title = f"Exons, Time Embedding w/ {gap}-mer Gap\n{e_count} Total Regions: weights are forwards\nCompliment = {compliment}"
        i_title = f"Introns, Time Embedding w/ {gap}-mer Gap\n{i_count} Total Regions: weights are forwards\nCompliment = {compliment}"

    x_title = f"History: {k_m}-Mer"
    y_title = f"Future: {k_p}-Mer"

    both_file = str(both_dir / f"both_gap_{gap}_{k_m}v{k_p}_Back_{backwards}_Comp_{compliment}.png")
    plt.title(b_title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.savefig(both_file)
    print(f"Output image to {both_file}")
    plt.close()

    # Exon Plot
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 20)
    for points in e_frame:
        ax.scatter(points[:, 0], points[:, 1], s = 0.1)
    exon_file = str(exon_dir / f"exon_gap_{gap}_{k_m}v{k_p}_Back_{backwards}_Comp_{compliment}.png")
    plt.title(e_title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.savefig(exon_file)
    print(f"Output image to {exon_file}")
    plt.close()

    # Intron Plot
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 20)
    for points in i_frame:
        ax.scatter(points[:, 0], points[:, 1], s = 0.1)
    intron_file = str(intron_dir / f"intron_gap_{gap}_{k_m}v{k_p}_Back_{backwards}_Comp_{compliment}.png")
    plt.title(i_title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.savefig(intron_file)
    print(f"Output image to {intron_file}")
    plt.close()


def time_embedding_v3(sequence: str, k_p = 9, k_m = 9, gap = 0, m_backwards = True, p_backwards = False, compliment = False):
    '''
    Feeds in a sequence, and it finds the xy coordinates for that sequence.

    There is an option for the compliment strand: probably should never be used.
    '''
    seq_length = len(sequence)

    if seq_length < (k_m + k_p + abs(gap)):  # I'm making this an |gap| becuase I don't want to think about how it should be done if g < 0. It has to be a certain length, and that length needs to be long.
        print("Cannont find Trajectory for this gene: to small")
        return None

    w_p = [(0.25)**n for n in range(1, k_p + 1)]
    w_m = [(0.25)**n for n in range(1, k_m + 1)]

    if m_backwards:
        w_m.reverse()
    if p_backwards:
        w_p.reverse()

    sequence = sequence.upper()

    xy = np.zeros(shape=(seq_length - (k_p + k_m + gap), 2))

    k_minus = [sequence[k_prime:k_prime + k_m] for k_prime in range(0, seq_length - (k_p + k_m + gap))]
    k_plus = [sequence[k_prime:k_prime + k_p] for k_prime in range(gap + k_m, seq_length - k_p)]

    if compliment:  # probably should never be used.
        for i, k_prime in enumerate(k_minus):
            n = [0 if n in "A" else 1 if n in "C" else 2 if n in "G" else 3 if n in "T" else 100 for n in k_prime]
            k_x = np.dot(w_m, n)

            n = [3 if n in "A" else 2 if n in "C" else 1 if n in "G" else 0 if n in "T" else 100 for n in k_plus[i]]
            k_y = np.dot(w_p, n)

            xy[i][0], xy[i][1] = k_x, k_y

    else:
        for i, k_prime in enumerate(k_minus):
            n = [0 if n in "A" else 1 if n in "C" else 2 if n in "G" else 3 if n in "T" else 100 for n in k_prime]
            k_x = np.dot(w_m, n)

            n = [0 if n in "A" else 1 if n in "C" else 2 if n in "G" else 3 if n in "T" else 100 for n in k_plus[i]]
            k_y = np.dot(w_p, n)

            xy[i][0], xy[i][1] = k_x, k_y

    return xy




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