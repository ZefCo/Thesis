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

    linux_path = f"/media/ethanspeakman/Elements/"
    windows_path = f"F:/"

    data_path = windows_path

    # # Zoomed in Picture    
    # x_lim = [0.75, 1.0]
    # y_lim = [0.25, 0.5]
    # x_ticks = {x/100:x/100 for x in range(75, 100 + 1, 1)}
    # y_ticks = {y/100:y/100 for y in range(25, 50 + 1, 1)}
    # s = 0.5

    # Non Zoomed Picture
    x_lim = None
    y_lim = None
    x_ticks = None
    y_ticks = None
    s = 0.01

    # s = 1
    # data_set = f"1&2"
    data_set = 1

    box1 = [[0, 0.8125], 1.0, 0.875 - 0.8125, ["blue", 4]]  # backward box
    
    box2 = [[0.75, 0.25], 0.25, 0.25, ["red", 4]]  # major box
    arrow = arrows([[0, 0.8125], [1.0, 0.875 - 0.8125]], [[0.75, 0.25], [0.25, 0.25]])
    # print(arrow)
    box2.append([arrow[0][0], arrow[0][1], arrow[1][0], arrow[1][1], 0.008, "blue"])
    
    box3 = [[0.4375, 0.0], 0.5 - 0.4375, 1.0, ["green", 4]]  # forward box
    arrow = arrows([[0.75, 0.25], [0.25, 0.25]], [[0.4375, 0.0], [0.5 - 0.4375, 1.0]])
    box3.append([arrow[0][0], arrow[0][1], arrow[1][0], arrow[1][1], 0.008, "red"])

    # print(box1)
    # print(box2)
    # print(box3)
    # exit()
    
    # box4 = [[0, 0.25], 0.25, 0.25, "green"] # comparison box
    
    boxes = [box1, box2, box3]
    # boxes = None
    # title = "Still working on placement of boxes"
    title = None

    # arrow_box_expl = [0 + 1.0/2, 0.8125 + (0.875 - 0.8125)/2, (0.75 + (0.25/2)) - (0 + 1.0/2), (0.25 + (0.25/2)) - (0.8125 + ((0.875 - 0.8125) / 2)) , 0.008, "blue"]
    # arrow_box_func = arrows([[0, 0.8125], [1.0, 0.875 - 0.8125], "blue"], [[0.75, 0.25], [0.25, 0.25], "red"])

    # print(arrow_box_expl)
    # print(arrow_box_func)
    # exit()

    n = 100_000
    if isinstance(x_lim, list):
        exon_dict_file = f"ExonData_n{n}_DS{data_set}_kp6_km6_zoomed_x{x_lim[0]}by{x_lim[1]}_y{y_lim[0]}by{y_lim[1]}"
    else:
        exon_dict_file = f"ExonData_n{n}_DS{data_set}_kp6_km6"
    if isinstance(y_lim, list):
        intron_dict_file = f"IntronData_n{n}_DS{data_set}_kp6_km6_zoomed_x{x_lim[0]}by{x_lim[1]}_y{y_lim[0]}by{y_lim[1]}"
    else:
        intron_dict_file = f"IntronData_n{n}_DS{data_set}_kp6_km6"

    # time_embedding_v2(pathlib.Path(f"{data_path}/Gene_Data_Sets/Data_Set_{data_set}_histogram.pkl"), 
    #                   n = n, 
    #                   backwards = True, 
    #                   x_lim = x_lim, y_lim = y_lim, 
    #                   ioxes = boxes,
    #                   eoxes = boxes,
    #                   exon_outfile = cwd / "TE_Images_ForPaper" / "Dict" / f"{exon_dict_file}.pkl",
    #                   intron_outfile = cwd / "TE_Images_ForPaper" / "Dict" / f"{intron_dict_file}.pkl")
    
    matplotfigure(cwd / "TE_Images_ForPaper" / "Dict" / f"{exon_dict_file}.pkl",
                  cwd / "TE_Images_ForPaper" / "Exon",
                  f"{exon_dict_file}.png",
                  x_lim=x_lim, y_lim=y_lim,
                  x_tick_marks=x_ticks, y_tick_marks=y_ticks,
                  boxes = boxes,
                  dot_size=s)
    
    matplotfigure(cwd / "TE_Images_ForPaper" / "Dict" / f"{intron_dict_file}.pkl",
                  cwd / "TE_Images_ForPaper" / "Intron",
                  f"{intron_dict_file}.png",
                  x_lim=x_lim, y_lim=y_lim,
                  x_tick_marks=x_ticks, y_tick_marks=y_ticks,
                  boxes = boxes,
                  dot_size=s)

    # reload_matplotlib(cwd / "TE_Images_ForPaper" / "Exon" / f"{exon_dict_file}_pkltest.pkl", cwd / "TE_Images_ForPaper" / "Exon" / f"{exon_dict_file}_pkltest.png", boxes = boxes)
    # reload_matplotlib(cwd / "TE_Images_ForPaper" / "Intron" / f"{intron_dict_file}_pkltest.pkl", cwd / "TE_Images_ForPaper" / "Intron" / f"{intron_dict_file}_pkltest.png", boxes = boxes)




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



def _clean_filepath(filepath: pathlib.Path) -> tuple:
    '''
    Just cleans the filepath and removes unnecessary .
    '''
    file_extension = filepath.suffix
    filepath = filepath.with_suffix("")
    filepath = pathlib.Path(str(filepath).replace(".", ""))
    filepath = filepath.with_suffix(file_extension)

    return filepath, file_extension



def arrows(xy0: np.ndarray or list, xy1: np.ndarray or list) -> list:
    '''
    Getting the arrows is a bit of a challenge, so I'm going to make a single function that will do the work for me.

    The arrows need x, y, dx, dy. I want to get all of those form x0, y0, dx0, dy0, x1, y1, dx1, dy1 which generate the boxes.

    The inputs should be a 2x2 numpy array, the first row being x, y, and the second row being dx, dy, or something that can be indexed twice like embedded lists.
    '''
    R = X = 0
    D = Y = 1
    r_half = lambda r, dr: r + (dr / 2)

    start = [r_half(xy0[R][X], xy0[D][X]), r_half(xy0[R][Y], xy0[D][Y])]
    delta = [r_half(xy1[R][X], xy1[D][X]) - r_half(xy0[R][X], xy0[D][X]), r_half(xy1[R][Y], xy1[D][Y]) - r_half(xy0[R][Y], xy0[D][Y])]

    return [start, delta]



def matplotfigure(frame: dict or pathlib.Path or str,
                  dir: str or pathlib.Path,
                  file_name: str or pathlib.Path,
                  k_p: int = 6, k_m: int = 6, gap: int = 0, 
                  backwards: bool = True, 
                  nucsequence: str = "AGTC",
                  title: str = None,
                  x_tick_marks: dict = None,
                  y_tick_marks: dict = None,
                  dot_size: float = 0.1,
                  boxes: list = None,
                  *args, **kwargs):
    '''
    The figure thing is redundent, so I'm going to try a single function to do everything.

    If frame == pathlib.Path or str then it loads that pickle file into the frame. dir is the output dir.

    Boxes allows you to draw on the plot. Boxes is a list of lists, and every sublist has an anchor point, a height, and a width. 
    The anchor point itself should be a tuple. So it's boxes = [[[x0, y0], h0, w0, [color0, linewidth0]], [[x1, y1], h1, w1, [color1, linewidth1]], ..., [[xn, yn], hn, wn, [colorn, linewidthn]]]
    A fifth element can be added which includes arrow information: [x, y, dx, dy, width, color]
    '''
    if isinstance(file_name, str):
        file_name: pathlib.Path = pathlib.Path(file_name)
    file_name, file_extension = _clean_filepath(file_name)
    file = dir / f"{file_name}"

    # print(file_name)
    # # print(file_extension)
    # exit()

    if isinstance(frame, pathlib.Path) or isinstance(frame, str):
        with open(frame, "rb") as f:
            frame = pickle.load(f)

    count = list(frame.keys())[len(list(frame.keys())) - 1]

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
        plt.plot(points[:, 0], points[:, 1], markersize = dot_size, marker = "s", linestyle = "", color = "k")
    plt.title(title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)

    if isinstance(boxes, list):
        for box in boxes:
            if isinstance(box[3], list):
                color = box[3][0]
                linewidth = box[3][1]
            else:
                color = None
                linewidth = 2
            p = plt.Rectangle(box[0], box[1], box[2], edgecolor = color, linewidth = linewidth, fill = False)  #set_fill = False, 
            ax.add_patch(p)

            if len(box) == 5:
                arrow = box[4]
                plt.arrow(arrow[0], arrow[1], arrow[2], arrow[3], width = arrow[4], facecolor = arrow[5], edgecolor = "none", length_includes_head = True)
    
    if isinstance(x_tick_marks, dict):
        ax.set_xticks(list(x_tick_marks.keys()))
        ax.set_xticklabels(list(x_tick_marks.values()))
    if isinstance(y_tick_marks, dict):
        ax.set_yticks(list(y_tick_marks.keys()))
        ax.set_yticklabels(list(y_tick_marks.values()))

    if (file_extension in ".png") or (file_extension in ".pdf"):
        plt.savefig(f"{file}")
    elif file_extension in ".pkl":
        pickle.dump(fig, open(f"{file}", "wb"))

    print(f"Output image to {file}")
    plt.close()


def reload_matplotlib(file_path: str or pathlib.Path,
                      output_file: str or pathlib.Path,
                      boxes: list = None):
    '''
    This is overall not working. There is some difficulty in trying to get the plots to actually redo things, and I dont' want to waste time on it. I'm going to make dictionaries and pickle
    those to save the data, then put it into a matplotlib and modify that. It will be a little slower but at least I can do things that way.


    Reloads a pickeled image. If the image was pickled this will reload it and allow the image to be zoomed in or out, or to put boxes on it. It's beacuse this beats
    having to do the same thing over and over again with N = very large number. Now I can make one plot with a very large N and then reload that same image to put boxes
    on them as I find them.
    '''
    file_path, _ = _clean_filepath(file_path)
    re_plt, ax = _steal_manager()

    with open(file_path, "rb") as file:
        try:
            re_plt: plt  = pickle.load(file)
        except Exception as e:
            print(type(e))
            print(e)
            exit()

    # print(type(fig))
    # fig.show()
    if isinstance(boxes, list):
        for box in boxes:
            if isinstance(box[3], str):
                color = box[3]
            else:
                color = None
            p = plt.Rectangle(box[0], box[1], box[2], edgecolor = color, linewidth = 2, fill = False)  #set_fill = False, 

            ax.add_patch(p)
    
    re_plt.savefig(output_file)
    print(f"Output image to {output_file}")



def _steal_manager():
    '''
    You can't just load the pickle file from matplot lib, you have to steal the canvas manager. This creates a dummy figure to hijack and use when loading a pickle file.
    '''
    fig, ax = plt.subplots()

    return fig, ax


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



def time_embedding_v2(data: pandas.DataFrame, 
                      n: int,
                      k_p: int = 6, k_m: int = 6, gap: int = 0, 
                      backwards: bool = True, 
                      nucsequence: str = "AGTC", PyPu: bool = False,
                      sequence_name: str = "Seq",
                      classification_name: str = "Classificaion",
                      x_lim: list = None,
                      y_lim: list = None,
                      exon_outfile: pathlib.Path or str = None,
                      intron_outfile: pathlib.Path or str = None,
                      *args, **kwargs):
    '''
    Rewrite this so that it pickles the data. Then we'll save that to a figure seperately.


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

    e_frame, i_frame = {}, {}  # OK this may seem weird, but hear me out: lists are slow, especially when you have a very large N in a list. But a dictionary is hashable, so I can store the lists in the dictionary, who cares about the key, and it should be faster for large N
    e_count, i_count = 0, 0

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

    with open(exon_outfile, "wb") as e_file:
        print(f"Writing to {exon_outfile}.pkl")
        pickle.dump(e_frame, e_file)
    with open(intron_outfile, "wb") as i_file:
        print(f"Writing to {intron_outfile}.pkl")
        pickle.dump(i_frame, i_file)

    # save the dictionary here

    # e_file_name: str = f"{output_file}_exon_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}"
    # matplotfigure(e_frame, exon_dir, e_count, 
    #               k_p = k_p, k_m = k_m, gap = gap, backwards = backwards, nucsequence = nucsequence, PyPu = PyPu, 
    #               title = title, file_name = e_file_name,
    #               x_lim = x_lim, y_lim = y_lim,
    #               x_tick_marks = x_tick_marks, y_tick_marks = y_tick_marks,
    #               dot_size = dot_size, boxes = eoxes, file_extension = file_extension)

    # i_file_name: str = f"{output_file}_intron_gap_{gap}_{k_m}v{k_p}_Back_{backwards}{file_NS}"
    # matplotfigure(i_frame, intron_dir, i_count, 
    #               k_p = k_p, k_m = k_m, gap = gap, backwards = backwards, nucsequence = nucsequence, PyPu = PyPu, 
    #               title = title, file_name = i_file_name,
    #               x_lim = x_lim, y_lim = y_lim,
    #               x_tick_marks = x_tick_marks, y_tick_marks = y_tick_marks,
    #               dot_size = dot_size, boxes = eoxes, file_extension = file_extension)


if __name__ in "__main__":
    main()