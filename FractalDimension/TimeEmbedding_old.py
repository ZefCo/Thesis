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
import GeneClass as Gene
import random

def main():
    '''
    I really need a new name for this, Time Embedding doens't work.

    k-Mer Time Series. That's a better name.

    Time embedding v3: make one that goes through a gene (all its different forms) and plots the trajectory of the k windows. Does something happen for the introns and the exons?
    You'll have to also do the exons and introns seperatly, but we want to see how the trajectory can "jump" from exon to intron: maybe there is something of interest there?
    Time embedding v3: make one that goes through a gene (all its different forms) and plots the trajectory of the k windows. Does something happen for the introns and the exons?
    You'll have to also do the exons and introns seperatly, but we want to see how the trajectory can "jump" from exon to intron: maybe there is something of interest there?
    '''
    # score_keys(k = 6)
    # exit()
    # test_seq = "GGCGGACCGGGCGTCCCTACCAAT"  # this sequence was created in the iterative Map becuase I'm having trouble understanding what the hell Dr G is purposing with his
    #                                        # multiple by 4 and add the y and take the fractional etc. I'm just going to study it using this.
    
    # xy = time_embedding(test_seq)

    # print(xy)

    linux_path = f"/media/ethanspeakman/Elements/"
    windows_path = f"F:/"

    data_path = windows_path
    
    x_lim = [0.75, 1.0]
    y_lim = [0.25, 0.5]
    x_ticks = {x/100:x/100 for x in range(75, 100 + 1, 1)}
    y_ticks = {y/100:y/100 for y in range(25, 50 + 1, 1)}
    # x_lim = None
    # y_lim = None
    # x_ticks = None
    # y_ticks = None

    # s = 1
    # s = 0.7
    s = 0.01
    # data_set = f"1&2"
    data_set = 1
    # box1 = [[0.75, 0.25], 0.25, 0.25, "red"]  # major box
    # box2 = [[0.4375, 0.0], 0.5 - 0.4375, 1.0, "red"]  # forward box
    # box3 = [[0, 0.8125], 1.0, 0.875 - 0.8125, "red"]  # backward box
    # box4 = [[0, 0.25], 0.25, 0.25, "green"] # comparison box
    # boxes = [box1, box2, box3, box4]
    boxes = None
    # title = "Still working on placement of boxes"
    title = None

    n = 100_000

    reload_matplotlib(cwd / "TE_Images_ForPaper" / "Intron" / "BF_DS1_boxed_n100000_intron_gap_0_6v6_Back_True_NucSeq_AGTC.pkl",
                      cwd / "TE_Images_ForPaper" / "Intron" / "Test_zoom.png",
                      x_lim = x_lim, y_lim = y_lim,
                      x_tick_marks = x_ticks, y_tick_marks = y_ticks)
    exit()

    time_embedding_v2(pathlib.Path(f"{data_path}/Gene_Data_Sets/Data_Set_{data_set}_histogram.pkl"), 
                      output_file = f"BF_DS{data_set}_boxed_n{n}", 
                      n = n, 
                      backwards = True, 
                      x_lim = x_lim, y_lim = y_lim, 
                      dot_size = s,
                      ioxes = boxes,
                      eoxes = boxes,
                      title = title,
                      x_tick_marks = x_ticks, y_tick_marks = y_ticks,
                      file_extension = ".pkl")

    # # back_forward_trajectories(cwd / "2mer_occ.csv", cwd / "2mer_occ_wScores.csv")
    
    # time_embedding_v2(pathlib.Path(f"G:\Gene_Data_Sets\Data_Set_{data_set}_histogram.pkl"), 
    #                   output_file = f"BF_DS{data_set}_zoom", 
    #                   n = 10_000, 
    #                   backwards = True, 
    #                   x_lim = x_lim, y_lim = y_lim, 
    #                   dot_size = s)

    # back_forward_trajectories(cwd / "2mer_occ.csv", cwd / "2mer_occ_wScores.csv")
    # time_embedding_v3(pathlib.Path(f"F:\Gene_Data_Sets\Data_Set_{data_set}_histogram.pkl"), 
    #                  output_file = f"BF_DS{data_set}_zoom", 
    #                  n = 1000,
    #                  N = 10)

    # thing = np.array([[1, 1], [1, 2], [3, 4], [5, 1], [1, 1]])
    # print(thing[1, 1])
    # something = np.where(thing[:, 0] == 1)
    # print(something)


def generate_sequence(k = 12, nucsequence: str = "AGTC"):
    '''
    Generates a random sequence.
    '''

    sequence = ""
    for _ in range(k):
        rand_n = random.randint(0, 3)
        sequence = f"{sequence}{nucsequence[rand_n]}"

    return sequence


def back_forward_trajectories(occ_data: pandas.DataFrame, output_file: pathlib.Path, kmer: int = 2):
    '''
    '''
    w_p = [(0.25)**n for n in range(1, kmer + 1)]
    w_m = [(0.25)**n for n in range(1, kmer + 1)]
    w_m.reverse()

    forwards, backwards = [], []

    x_ave = (0.25**kmer)*100

    if isinstance(occ_data, pathlib.Path) or isinstance(occ_data, str):
        occ_data = pandas.read_csv(occ_data, header = 0, index_col = "Seq")
    
    occ_data["Occ"] = (occ_data["Occ"]*100) + x_ave

    seq_mer = tuple(occ_data.index)

    for seq in seq_mer:
        digital_seq = [0 if n in "A" else 1 if n in "G" else 2 if n in "T" else 3 for n in seq]
        forwards.append(np.dot(digital_seq, w_p))
        backwards.append(np.dot(digital_seq, w_m))

    occ_data["Foward"] = forwards
    occ_data["Backward"] = backwards

    print(occ_data["Occ"].sum())

    occ_data.to_csv(output_file)


def score_keys(k = 9, nucsequence: str = "AGTC"):
    '''
    Generates a score key. The sequence used by default is AGTC, which can be changed by inputing a new string as CGTA (for example).
    '''
    perms = tuple(itertools.product([nucsequence[0], nucsequence[1], nucsequence[2], nucsequence[3]], repeat = k))

    w_p = [(0.25)**n for n in range(1, k + 1)]
    w_m = [(0.25)**n for n in range(1, k + 1)]
    w_m.reverse()

    scores = {}
    for p in perms:
        key = ""
        for n in p:
            key = f"{key}{n}"

        scores[key] = {"history": 0, "future": 0}

        digital = [0 if n in nucsequence[0] else 1 if n in nucsequence[1] else 2 if n in nucsequence[2] else 3 if n in nucsequence[3] else 100 for n in key]

        score_p = np.dot(w_p, digital)
        score_m = np.dot(w_m, digital)

        scores[key]["future"] = score_p
        scores[key]["history"] = score_m

    scores = pandas.DataFrame(data = scores).T
    print(scores)

    scores.to_csv(cwd / f"ScoreKey_kmer_{k}_{nucsequence}.csv")



def matplotfigure(frame: dict,
                  dir: str = None,
                  count: int = 0,
                  k_p: int = 6, k_m: int = 6, gap: int = 0, 
                  backwards: bool = True, 
                  nucsequence: str = "AGTC",
                  title: str = None,
                  file_name: str = None,
                  x_lim: list = None,
                  y_lim: list = None,
                  x_tick_marks: dict = None,
                  y_tick_marks: dict = None,
                  dot_size: float = 0.1,
                  boxes: list = None,
                  file_extension: str = ".png",
                  *args, **kwargs):
    '''
    The figure thing is redundent, so I'm going to try a single function to do everything.
    '''

    if isinstance(title, str):
        title = f"{title}\nTime Embedding w/ {gap}-mer Gap\n{count} Total Regions"
    else:
        title = f"Time Embedding w/ {gap}-mer Gap\n{count} Total Regions"

    if backwards:
        weights = ": weights are forwards and backwards"
    else:
        weights = ": weights are forwards"


    # b_title = f"{b_title}{weights}{NS}"
    title = f"{title}{weights}{nucsequence}"

    x_title = f"History: {k_m}-Mer"
    y_title = f"Future: {k_p}-Mer"

    # Plot
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 20)
    for i, points in enumerate(frame.values()):
        # ax.scatter(points[:, 0], points[:, 1], s = dot_size, marker = "s", color = "k", linewidths = 0)
        plt.plot(points[:, 0], points[:, 1], "s", markersize = dot_size, marker = "s", color = "k")
    plt.title(title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    if isinstance(x_lim, list):
        plt.xlim(x_lim[0], x_lim[1])
        file_name = f"{file_name}_x_{x_lim[0]}v{x_lim[1]}"
    if isinstance(y_lim, list):
        plt.ylim(y_lim[0], y_lim[1])
        file_name = f"{file_name}_y_{y_lim[0]}v{y_lim[1]}"
    file = str(dir / f"{file_name}")

    if isinstance(boxes, list):
        for box in boxes:
            if isinstance(box[3], str):
                color = box[3]
            else:
                color = None
            p = plt.Rectangle(box[0], box[1], box[2], edgecolor = color, linewidth = 2, fill = False)  #set_fill = False, 

            ax.add_patch(p)
    
    if isinstance(x_tick_marks, dict):
        ax.set_xticks(list(x_tick_marks.keys()))
        ax.set_xticklabels(list(x_tick_marks.values()))
    if isinstance(y_tick_marks, dict):
        ax.set_yticks(list(y_tick_marks.keys()))
        ax.set_yticklabels(list(y_tick_marks.values()))

    if (file_extension in ".png") or (file_extension in ".pdf"):
        plt.savefig(f"{file}{file_extension}")
    elif file_extension in ".pkl":
        pickle.dump(fig, open(f"{file}{file_extension}", "wb"))

    print(f"Output image to {file}{file_extension}")
    plt.close()


def reload_matplotlib(file_path: str or pathlib.Path,
                      output_file: str or pathlib.Path,
                      x_lim: list = None, y_lim: list = None,
                      x_tick_marks: dict = None, y_tick_marks: dict = None,
                      dot_size: float = 0.1,
                      boxes: list = None):
    '''
    Reloads a pickeled image. If the image was pickled this will reload it and allow the image to be zoomed in or out, or to put boxes on it. It's beacuse this beats
    having to do the same thing over and over again with N = very large number. Now I can make one plot with a very large N and then reload that same image to put boxes
    on them as I find them.
    '''
    with open(file_path, "rb") as ax:
        re_plt: plt  = pickle.load(ax)

    if isinstance(x_lim, list):
        plt.xlim(x_lim[0], x_lim[1])
    if isinstance(y_lim, list):
        plt.ylim(y_lim[0], y_lim[1])

    # # print(type(fig))
    # # fig.show()
    # if isinstance(boxes, list):
    #     for box in boxes:
    #         if isinstance(box[3], str):
    #             color = box[3]
    #         else:
    #             color = None
    #         p = plt.Rectangle(box[0], box[1], box[2], edgecolor = color, linewidth = 2, fill = False)  #set_fill = False, 

    #         ax.add_patch(p)
    
    # if isinstance(x_tick_marks, dict):
    #     ax.set_xticks(list(x_tick_marks.keys()))
    #     ax.set_xticklabels(list(x_tick_marks.values()))
    # if isinstance(y_tick_marks, dict):
    #     ax.set_yticks(list(y_tick_marks.keys()))
    #     ax.set_yticklabels(list(y_tick_marks.values()))

    re_plt.savefig(output_file)
    print(f"Output image to {output_file}")



def old_time_embedding_plotly(k_p = 9, k_m = 9, gap = 0, max_rows = 200):
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




def old_time_embedding_matplot(pickle_file, k_p = 9, k_m = 9, gap = 0, max_rows = 200, backwards = True, PyPu = False, nucsequence: str = "AGTC"):
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

        if PyPu:
            xy = time_embedding_PyPu(sequence, k_p = k_p, k_m = k_m, gap = gap, m_backwards = backwards)
        else:
            xy = time_embedding(sequence, k_p = k_p, k_m = k_m, gap = gap, m_backwards = backwards, nucsequence = nucsequence)
        
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

    b_title = f"Exons and Introns, Time Embedding w/ {gap}-mer Gap between + and -\n{e_count + i_count} Total Regions - Unknown number of genes are represented"
    e_title = f"Exons, Time Embedding w/ {gap}-mer Gap\n{e_count} Total Regions: weights are forwards and backwards"
    i_title = f"Introns, Time Embedding w/ {gap}-mer Gap\n{i_count} Total Regions: weights are forwards and backwards"

    if backwards:
        weights = ": weights are forwards and backwards"
    else:
        weights = ": weights are forwards"

    if PyPu:
        NS = f"\nPy = 0, Pu = 1"
        file_NS = f"_PyPu_{PyPu}"
    else:
        NS = f"\n{nucsequence}"
        file_NS = f"_NucSeq_{nucsequence}"

    b_title = f"{b_title}{weights}{NS}"
    e_title = f"{e_title}{weights}{NS}"
    i_title = f"{i_title}{weights}{NS}"


    x_title = f"History: {k_m}-Mer"
    y_title = f"Future: {k_p}-Mer"

    both_file = str(both_dir / f"both_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}.png")
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
    exon_file = str(exon_dir / f"exon_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}.png")
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
    intron_file = str(intron_dir / f"intron_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}.png")
    plt.title(i_title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.savefig(intron_file)
    print(f"Output image to {intron_file}")
    plt.close()


def time_embedding(sequence: str, 
                      k_p: int = 6, k_m: int = 6, gap: int = 0, 
                      m_backwards: bool = True, p_backwards: bool = False, 
                    #   compliment: bool = False, 
                      nucsequence: str = "AGTC"):
    '''
    Feeds in a sequence, and it finds the xy coordinates for that sequence.

    The nucleotide to number order can be altered. By default it is A = 0, G = 1, T = 2, C = 3. To alter it just feed in a new str with your prefered order. The first index is 0, the next
    index is 1, and so on.

    There is an option for the compliment strand: probably should never be used.
    '''
    sequence = sequence.upper()
    nucsequence = nucsequence.upper() # Just in case someone puts in a different order and forgets to capitalize everything
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


    xy = np.zeros(shape=(seq_length - (k_p + k_m + gap), 2))

    k_minus = [sequence[k_prime:k_prime + k_m] for k_prime in range(0, seq_length - (k_p + k_m + gap))]
    k_plus = [sequence[k_prime:k_prime + k_p] for k_prime in range(gap + k_m, seq_length - k_p)]

    # if compliment:  # probably should never be used.
    #     for i, k_prime in enumerate(k_minus):
    #         n = [0 if n in nucsequence[0] else 1 if n in nucsequence[1] else 2 if n in nucsequence[2] else 3 if n in nucsequence[3] else 100 for n in k_prime]
    #         k_x = np.dot(w_m, n)

    #         n = [3 if n in nucsequence[0] else 2 if n in nucsequence[1] else 1 if n in nucsequence[2] else 0 if n in nucsequence[3] else 100 for n in k_plus[i]]
    #         k_y = np.dot(w_p, n)

    #         xy[i][0], xy[i][1] = k_x, k_y

    # else:
    for i, k_prime in enumerate(k_minus):
        n = [0 if n in nucsequence[0] else 1 if n in nucsequence[1] else 2 if n in nucsequence[2] else 3 if n in nucsequence[3] else 100 for n in k_prime]
        k_x = np.dot(w_m, n)

        n = [0 if n in nucsequence[0] else 1 if n in nucsequence[1] else 2 if n in nucsequence[2] else 3 if n in nucsequence[3] else 100 for n in k_plus[i]]
        k_y = np.dot(w_p, n)

        xy[i][0], xy[i][1] = k_x, k_y

    return xy


def time_embedding_PyPu(sequence: str, 
                        k_p: int = 6, k_m: int = 6, gap: int = 0, 
                        m_backwards: bool = True, p_backwards: bool = False):
    '''
    Does the time embedding, but this one is based off of Pyrimdines and Purines.
    '''
    sequence = sequence.upper()
    seq_length = len(sequence)

    if seq_length < (k_m + k_p + abs(gap)):  # I'm making this an |gap| becuase I don't want to think about how it should be done if g < 0. It has to be a certain length, and that length needs to be long.
        print("Cannont find Trajectory for this gene: to small")
        return None

    w_p = [(0.5)**n for n in range(1, k_p + 1)]
    w_m = [(0.5)**n for n in range(1, k_m + 1)]

    if m_backwards:
        w_m.reverse()
    if p_backwards:
        w_p.reverse()


    xy = np.zeros(shape=(seq_length - (k_p + k_m + gap), 2))

    k_minus = [sequence[k_prime:k_prime + k_m] for k_prime in range(0, seq_length - (k_p + k_m + gap))]
    k_plus = [sequence[k_prime:k_prime + k_p] for k_prime in range(gap + k_m, seq_length - k_p)]

    for i, k_prime in enumerate(k_minus):
        n = [0 if ((n in "A") or (n in "G")) else 1 if ((n in "C") or (n in "T")) else 100 for n in k_prime]
        k_x = np.dot(w_m, n)

        n = [0 if ((n in "A") or (n in "G")) else 1 if ((n in "C") or (n in "T")) else 100 for n in k_plus[i]]
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


def time_embedding_v2a(data: pandas.DataFrame, 
                      n: int,
                      k_p: int = 6, k_m: int = 6, gap: int = 0, 
                      backwards: bool = True, 
                      nucsequence: str = "AGTC", PyPu: bool = False,
                      sequence_name: str = "Seq",
                      classification_name: str = "Classificaion",
                      title: str = None,
                      output_file: str = None,
                      x_lim: list = None,
                      y_lim: list = None,
                      x_tick_marks: dict = None,
                      y_tick_marks: dict = None,
                      dot_size: float = 0.1,
                      eoxes: list = None, ioxes: list = None,
                      file_extension: str = ".png",
                      *args, **kwargs):
    '''
    The new way of doing things with the updated data. You can put a pathlib in place of a Dataframe which then opens that file, but that data better be a Dataframe. I'm not going to code
    other ways of handeling that data. It's a Dataframe. Deal with it. WE'RE DEALING WITH THINGS TED!

    n represents the number of samples to be taken. I can't use all of them because some datasets have as much as 150,000, which basically freezes matplotlib.

    The classification & sequence name is because I've done this a few different times with different data standards... and I probably screwed myself for that. But now I get to try to salvage that...
    And yes, that says Classificaion, because I misspelled something.

    Boxes allows you to draw on the plot. Boxes is a list of lists, and every sublist has an anchor point, a height, and a width. The anchor point itself should be a tuple. So it's boxes = [[[x0, y0], h0, w0, color0], [[x1, y1], h1, w1, color1], ..., [[xn, yn], hn, wn, colorn]]
    '''

    if isinstance(data, pathlib.Path):
        with open(data, "rb") as p:
            data = pickle.load(p)

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


    print(data.columns)

    image_dir = cwd / "TE_Images_ForPaper"
    image_dir.mkdir(parents = True, exist_ok = True)
    exon_dir = image_dir / "Exon"
    exon_dir.mkdir(parents = True, exist_ok = True)
    intron_dir = image_dir / "Intron"
    intron_dir.mkdir(parents = True, exist_ok = True)
    both_dir = image_dir / "Both"
    both_dir.mkdir(parents = True, exist_ok = True)

    # gene: Gene.Gene
    # ncib: str

    rows, cols = data.shape

    if "Length" in data.columns:
        data = data[data["Length"] > (k_m + k_p + gap)]
    else:
        data["Length"] = data.Seq.str.len()
        data = data[data["Length"] > (k_m + k_p + gap)]

    b_frame, e_frame, i_frame = {}, {}, {}  # OK this may seem weird, but hear me out: lists are slow, especially when you have a very large N in a list. But a dictionary is hashable, so I can store the lists in the dictionary, who cares about the key, and it should be faster for large N
    b_count, e_count, i_count = 0, 0, 0

    for row in range(rows):
        try:
            sequence = data.loc[row, sequence_name].upper()
        except Exception as e:
            print(type(e))
            print(e)
            print(row)
            continue
            exit()

        region = data.loc[row, classification_name]

        if PyPu:
            xy = time_embedding_PyPu(sequence, k_p = k_p, k_m = k_m, gap = gap, m_backwards = backwards)
        else:
            xy = time_embedding(sequence, k_p = k_p, k_m = k_m, gap = gap, m_backwards = backwards, nucsequence = nucsequence)

        if isinstance(x_lim, list):
            xy = xy[(x_lim[1] >= xy[:, 0]) & (xy[:, 0] >= x_lim[0])]
            # xy = np.where((x_lim[1] >= xy[:, 0]) & (xy[:, 0] >= x_lim[0]))  # Filter the number of points down to a zoomed in area. Should help speed up the process of making the actual graph
        if isinstance(y_lim, list):
            xy = xy[(y_lim[1] >= xy[:, 1]) & (xy[:, 1] >= y_lim[0])]
            # xy = np.where((y_lim[1] >= xy[:, 1]) & (xy[:, 1] >= y_lim[0]))
        
        if (region in "Exon") or (region in "exon") or (region in "e"):  # because I realized I was doing this like 3 different ways... I probably should have been more precise.
            e_frame[e_count] = xy
            e_count += 1
        elif (region in "Intron") or (region in "intron") or (region in "i"):
            i_frame[i_count] = xy
            i_count += 1

        # b_frame[b_count] = xy
        # b_count += 1

        if ((row % 1000) == 0):
            print(f"Finished row {row} of {rows}")


    # Both Plot
    # fig, ax = plt.subplots()
    # fig.set_size_inches(20, 20)
    # for points in b_frame:
    #     ax.scatter(points[:, 0], points[:, 1], s = 0.1)

    # b_title = f"Exons and Introns, Time Embedding w/ {gap}-mer Gap between + and -\n{e_count + i_count} Total Regions - Unknown number of genes are represented"
    if isinstance(title, str):
        e_title = f"{title}\nExons, Time Embedding w/ {gap}-mer Gap\n{e_count} Total Regions"
        i_title = f"{title}\nIntrons, Time Embedding w/ {gap}-mer Gap\n{i_count} Total Regions"
    else:
        e_title = f"Exons, Time Embedding w/ {gap}-mer Gap\n{e_count} Total Regions"
        i_title = f"Introns, Time Embedding w/ {gap}-mer Gap\n{i_count} Total Regions"

    if backwards:
        weights = ": weights are forwards and backwards"
    else:
        weights = ": weights are forwards"

    if PyPu:
        NS = f"\nPy = 0, Pu = 1"
        file_NS = f"_PyPu_{PyPu}"
    else:
        NS = f"\n{nucsequence}"
        file_NS = f"_NucSeq_{nucsequence}"

    # b_title = f"{b_title}{weights}{NS}"
    e_title = f"{e_title}{weights}{NS}"
    i_title = f"{i_title}{weights}{NS}"


    x_title = f"History: {k_m}-Mer"
    y_title = f"Future: {k_p}-Mer"

    # both_file = str(both_dir / f"both_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}.png")
    # plt.title(b_title)
    # plt.xlabel(x_title)
    # plt.ylabel(y_title)
    # plt.savefig(both_file)
    # print(f"Output image to {both_file}")
    # plt.close()

    # Exon Plot
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 20)
    for points in e_frame.values():
        ax.scatter(points[:, 0], points[:, 1], s = dot_size, marker = "s", color = "k")
    plt.title(e_title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    e_file_name: str = f"{output_file}_exon_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}"
    if isinstance(x_lim, list):
        plt.xlim(x_lim[0], x_lim[1])
        e_file_name = f"{e_file_name}_x_{x_lim[0]}v{x_lim[1]}"
    if isinstance(y_lim, list):
        plt.ylim(y_lim[0], y_lim[1])
        e_file_name = f"{e_file_name}_y_{y_lim[0]}v{y_lim[1]}"
    exon_file = str(exon_dir / f"{e_file_name}")

    if isinstance(eoxes, list):
        for box in eoxes:
            if isinstance(box[3], str):
                color = box[3]
            else:
                color = None
            p = plt.Rectangle(box[0], box[1], box[2], edgecolor = color, linewidth = 2, fill = False)  #set_fill = False, 

            ax.add_patch(p)
    
    if isinstance(x_tick_marks, dict):
        ax.set_xticks(list(x_tick_marks.keys()))
        ax.set_xticklabels(list(x_tick_marks.values()))
    if isinstance(y_tick_marks, dict):
        ax.set_yticks(list(y_tick_marks.keys()))
        ax.set_yticklabels(list(y_tick_marks.values()))

    plt.savefig(exon_file)
    print(f"Output image to {exon_file}")
    plt.close()

    # Intron Plot
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 20)
    for points in i_frame.values():
        ax.scatter(points[:, 0], points[:, 1], s = dot_size, marker = "s", color = "k")
    plt.title(i_title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    i_file_name: str = f"{output_file}_intron_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}"
    if isinstance(x_lim, list):
        plt.xlim(x_lim[0], x_lim[1])
        i_file_name = f"{i_file_name}_x_{x_lim[0]}v{x_lim[1]}"
    if isinstance(y_lim, list):
        plt.ylim(y_lim[0], y_lim[1])
        i_file_name = f"{i_file_name}_y_{y_lim[0]}v{y_lim[1]}"
    intron_file = str(intron_dir / f"{i_file_name}")

    if isinstance(ioxes, list):
        for box in ioxes:
            if isinstance(box[3], str):
                color = box[3]
            else:
                color = None
            p = plt.Rectangle(box[0], box[1], box[2], edgecolor = color, linewidth = 2, fill = False)  #set_fill = False, 

            ax.add_patch(p)
    
    if isinstance(x_tick_marks, dict):
        ax.set_xticks(list(x_tick_marks.keys()))
        ax.set_xticklabels(list(x_tick_marks.values()))
    if isinstance(y_tick_marks, dict):
        ax.set_yticks(list(y_tick_marks.keys()))
        ax.set_yticklabels(list(y_tick_marks.values()))

    plt.savefig(intron_file)
    print(f"Output image to {intron_file}")
    plt.close()



def time_embedding_v2(data: pandas.DataFrame, 
                      n: int,
                      k_p: int = 6, k_m: int = 6, gap: int = 0, 
                      backwards: bool = True, 
                      nucsequence: str = "AGTC", PyPu: bool = False,
                      sequence_name: str = "Seq",
                      classification_name: str = "Classificaion",
                      title: str = None,
                      output_file: str = None,
                      x_lim: list = None,
                      y_lim: list = None,
                      x_tick_marks: dict = None,
                      y_tick_marks: dict = None,
                      dot_size: float = 0.1,
                      eoxes: list = None, ioxes: list = None,
                      file_extension: str = ".png",
                      *args, **kwargs):
    '''
    The new way of doing things with the updated data. You can put a pathlib in place of a Dataframe which then opens that file, but that data better be a Dataframe. I'm not going to code
    other ways of handeling that data. It's a Dataframe. Deal with it. WE'RE DEALING WITH THINGS TED!

    n represents the number of samples to be taken. I can't use all of them because some datasets have as much as 150,000, which basically freezes matplotlib.

    The classification & sequence name is because I've done this a few different times with different data standards... and I probably screwed myself for that. But now I get to try to salvage that...
    And yes, that says Classificaion, because I misspelled something.

    Boxes allows you to draw on the plot. Boxes is a list of lists, and every sublist has an anchor point, a height, and a width. The anchor point itself should be a tuple. So it's boxes = [[[x0, y0], h0, w0, color0], [[x1, y1], h1, w1, color1], ..., [[xn, yn], hn, wn, colorn]]
    '''
    if PyPu:
        NS = f"\nPy = 0, Pu = 1"
        file_NS = f"_PyPu_{PyPu}"
    else:
        NS = f"\n{nucsequence}"
        file_NS = f"_NucSeq_{nucsequence}"

    if isinstance(data, pathlib.Path):
        with open(data, "rb") as p:
            data = pickle.load(p)

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


    print(data.columns)

    image_dir = cwd / "TE_Images_ForPaper"
    image_dir.mkdir(parents = True, exist_ok = True)
    exon_dir = image_dir / "Exon"
    exon_dir.mkdir(parents = True, exist_ok = True)
    intron_dir = image_dir / "Intron"
    intron_dir.mkdir(parents = True, exist_ok = True)
    both_dir = image_dir / "Both"
    both_dir.mkdir(parents = True, exist_ok = True)

    # gene: Gene.Gene
    # ncib: str

    rows, cols = data.shape

    if "Length" in data.columns:
        data = data[data["Length"] > (k_m + k_p + gap)]
    else:
        data["Length"] = data.Seq.str.len()
        data = data[data["Length"] > (k_m + k_p + gap)]

    b_frame, e_frame, i_frame = {}, {}, {}  # OK this may seem weird, but hear me out: lists are slow, especially when you have a very large N in a list. But a dictionary is hashable, so I can store the lists in the dictionary, who cares about the key, and it should be faster for large N
    b_count, e_count, i_count = 0, 0, 0

    for row in range(rows):
        try:
            sequence = data.loc[row, sequence_name].upper()
        except Exception as e:
            print(type(e))
            print(e)
            print(row)
            continue
            exit()

        region = data.loc[row, classification_name]

        if PyPu:
            xy = time_embedding_PyPu(sequence, k_p = k_p, k_m = k_m, gap = gap, m_backwards = backwards)
        else:
            xy = time_embedding(sequence, k_p = k_p, k_m = k_m, gap = gap, m_backwards = backwards, nucsequence = nucsequence)

        if isinstance(x_lim, list):
            xy = xy[(x_lim[1] >= xy[:, 0]) & (xy[:, 0] >= x_lim[0])]
            # xy = np.where((x_lim[1] >= xy[:, 0]) & (xy[:, 0] >= x_lim[0]))  # Filter the number of points down to a zoomed in area. Should help speed up the process of making the actual graph
        if isinstance(y_lim, list):
            xy = xy[(y_lim[1] >= xy[:, 1]) & (xy[:, 1] >= y_lim[0])]
            # xy = np.where((y_lim[1] >= xy[:, 1]) & (xy[:, 1] >= y_lim[0]))
        
        if (region in "Exon") or (region in "exon") or (region in "e"):  # because I realized I was doing this like 3 different ways... I probably should have been more precise.
            e_frame[e_count] = xy
            e_count += 1
        elif (region in "Intron") or (region in "intron") or (region in "i"):
            i_frame[i_count] = xy
            i_count += 1

        # b_frame[b_count] = xy
        # b_count += 1

        if ((row % 1000) == 0):
            print(f"Finished row {row} of {rows}")

    e_file_name: str = f"{output_file}_exon_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}"
    matplotfigure(e_frame, exon_dir, e_count, 
                  k_p = k_p, k_m = k_m, gap = gap, backwards = backwards, nucsequence = nucsequence, PyPu = PyPu, 
                  title = title, file_name = e_file_name,
                  x_lim = x_lim, y_lim = y_lim,
                  x_tick_marks = x_tick_marks, y_tick_marks = y_tick_marks,
                  dot_size = dot_size, boxes = eoxes, file_extension = file_extension)

    i_file_name: str = f"{output_file}_intron_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}"
    matplotfigure(i_frame, intron_dir, i_count, 
                  k_p = k_p, k_m = k_m, gap = gap, backwards = backwards, nucsequence = nucsequence, PyPu = PyPu, 
                  title = title, file_name = i_file_name,
                  x_lim = x_lim, y_lim = y_lim,
                  x_tick_marks = x_tick_marks, y_tick_marks = y_tick_marks,
                  dot_size = dot_size, boxes = eoxes, file_extension = file_extension)





def time_embedding_v3(data: pandas.DataFrame,
                      n: int = 10_000,
                      k_m: int = 6, k_p: int = 6, gap: int = 0,
                      nucsequence: str = "AGTC", sequence_name: str = "Seq", classification_name: str = "Classificaion",
                      N: int = 10,
                      *args, **kwargs):
    '''
    Maybe I should rename this to "Danish Pastry" method.

    
    I have this idea that I want to write down before I forget it: create a numpy vector that is n**2 x 3 for the points in the time embedding map (danish pastrey?). it would be [[x y z]] where x and y are the points and z is the count.
    That way if I want a specific area I can just filter by either x vs y or by z.
    '''
    # master_xy = np.zeros(shape = (1, 3))
    backwards: bool = True
    if isinstance(data, pathlib.Path):
        with open(data, "rb") as p:
            data = pickle.load(p)

    # print(data.iloc[0, :])
    # exit()

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

    # data = data.iloc[0:2, :]
    # data.loc[0, "Seq"] = "AGTCGGAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTC"

    # print(data.columns)

    image_dir = cwd / "TE_Images_ForPaper"
    image_dir.mkdir(parents = True, exist_ok = True)
    exon_dir = image_dir / "Exon"
    exon_dir.mkdir(parents = True, exist_ok = True)
    intron_dir = image_dir / "Intron"
    intron_dir.mkdir(parents = True, exist_ok = True)
    both_dir = image_dir / "Both"
    both_dir.mkdir(parents = True, exist_ok = True)

    # gene: Gene.Gene
    # ncib: str

    rows, cols = data.shape
    # rows = 1

    if "Length" in data.columns:
        data = data[data["Length"] > (k_m + k_p + gap)]
    else:
        data["Length"] = data.Seq.str.len()
        data = data[data["Length"] > (k_m + k_p + gap)]

    b_frame, e_frame, i_frame = {}, {}, {}  # OK this may seem weird, but hear me out: lists are slow, especially when you have a very large N in a list. But a dictionary is hashable, so I can store the lists in the dictionary, who cares about the key, and it should be faster for large N
    b_count, e_count, i_count = 0, 0, 0

    for row in range(rows):
        try:
            sequence = data.loc[row, sequence_name].upper()
        except Exception as e:
            print(type(e))
            print(e)
            print(row)
            continue
            exit()

        region = data.loc[row, classification_name]

        xy = time_embedding(sequence, k_p = k_p, k_m = k_m, gap = gap, m_backwards = backwards, nucsequence = nucsequence)
        # print(xy)
        # exit()

        # Need to figure out a faster way instead of iterating and adding it one at a time.

        # This is where things are going to get crazy. I have to find the xy points, check if they're in the master thing, add them if they are not, and add to the z position if they are.
        for j, r in enumerate(xy):
            index = f"{r[0]}_{r[1]}"
            if index not in b_frame.keys():
                b_frame[index] = {"x": r[0], "y": r[1], "z": 1}
            else:
                b_frame[index]["z"] += 1

            if (region in "Exon") or (region in "exon") or (region in "e"):  # because I realized I was doing this like 3 different ways... I probably should have been more precise.
                if index not in e_frame.keys():
                    e_frame[index] = {"x": r[0], "y": r[1], "z": 1}
                else:
                    e_frame[index]["z"] += 1
            elif (region in "Intron") or (region in "intron") or (region in "i"):
                if index not in i_frame.keys():
                    i_frame[index] = {"x": r[0], "y": r[1], "z": 1}
                else:
                    i_frame[index]["z"] += 1

        if ((row % 1000) == 0):
            print(f"Finished row {row} of {rows}")
    # print(b_frame)
    # print(e_frame)
    # print(i_frame)
        
    # print(master_xy)

    # sorted_e = dict(sorted(e_frame.items(), key = lambda x:x[1], reverse = True)[:N])
    # print(sorted_e)
    master_xy = np.zeros(shape = (len(e_frame), 3))
    for i, (_, value) in enumerate(e_frame.items()):
        master_xy[i, 0] = value["x"]
        master_xy[i, 1] = value["y"]
        master_xy[i, 2] = value["z"]

    print(np.unique(master_xy[:, 2]))
        

if __name__ in "__main__":
    main()