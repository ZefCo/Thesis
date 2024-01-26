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
    Zoom in on CGA|CTA

    x = CTA
    y = AGC
    box1 = A|G
    box2 = GA|CT
    box3 = CGA|CTA

    0.25 -> 1/4
    0.0625 -> 1/8
    0.015625 -> 1/16
    0.00390625 -> 1/32
    0.0009765625 -> 1/64

    So traditionally I did the backwards on the x axis and forwards on the y axis, and that's getting me into trouble. The Kneeding transform is best done with it flipped. So I introduced a simple
    transform bool, but... that's getting confusing. The point is: be careful how you use this.


    Do this with different permutations of nucelotides.

    Break this out and create a new script for playing with the boxes and arrows. Use something like: Green = 0th, Red = -N iterations, Blue = +N iterations. Have it fade or try lighter shades every iteration.
    Also need to figure out a way to get the density of the points in these regions. Want to compare the Exon to Intron, so target low density regions. Then we want to count the points. Maybe consider a correlation dimension (my idea, not Dr. Gs)

    I really need a new name for this, Time Embedding doens't work.

    k-Mer Time Series. That's a better name.

    Time embedding v3: make one that goes through a gene (all its different forms) and plots the trajectory of the k windows. Does something happen for the introns and the exons?
    You'll have to also do the exons and introns seperatly, but we want to see how the trajectory can "jump" from exon to intron: maybe there is something of interest there?
    Time embedding v3: make one that goes through a gene (all its different forms) and plots the trajectory of the k windows. Does something happen for the introns and the exons?
    You'll have to also do the exons and introns seperatly, but we want to see how the trajectory can "jump" from exon to intron: maybe there is something of interest there?
    '''
    # score_keys(k = 6)
    # exit()
    species = "Mouse"
    # species = "Fly"

    linux_path = f"/media/ethanspeakman/Elements/"
    windows_path = f"F:/"

    data_path = windows_path

    # Non Zoomed Picture
    exon_dict_file = f"{species}_ExonData_n100000_DS1_kp6_km6"
    intron_dict_file = f"{species}_IntronData_n100000_DS1_kp6_km6"
    x_lim = None
    y_lim = None
    x_ticks = {0.0: 0.0, 1.0: 1.0}
    y_ticks = {0.0: 0.0, 1.0: 1.0}
    s = 0.01
    bfa = None
    inches = 20
    exon_png_file = f"{species}_Exon_KT"
    intron_png_file = f"{species}_Intron_KT"
    box_1 = [0.75, 0.00, 0.25, 0.25, "red", 4, True, 0.2]
    bfa = [box_1]
    bfa = None
    border = 8
    # bfa = None
    textsize = 40
    n = 100_000
    data_set = 1
    x_title, y_title = "X", "Y"
    label_coords = [0.5, -0.02, -0.02, 0.5]  # [None, None, -0.1, 0.5] for "Data_zoomed_x0875y00625T"

    # # Zoomed in Picture
    # general_dict_file = "Data_zoomed_x025y075T"
    # exon_dict_file = f"Exon{general_dict_file}"
    # intron_dict_file = f"Intron{general_dict_file}"
    # x_lim = [0.109375, 0.125]
    # y_lim = [0.875, 0.890625]
    # y_ticks = {0.75: 0.75, 1.0: 1.0}
    # x_ticks = {0.25: 0.25, 0.5: 0.5}
    # border = 8
    # s = 0.5
    # exon_png_file = "ExonZoomed_x025_050_y075_100"
    # intron_png_file = "IntronZoomed_x025_050_y075_100"
    # inches = 20
    # textsize = 40
    # n = 100_000
    # data_set = 1
    # x_title, y_title = "X", "Y"
    # label_coords = [0.5, -0.02, -0.02, 0.5]  # [None, None, -0.1, 0.5] for "Data_zoomed_x0875y00625T"
    # box_1 = [0.875, 0.109375, 0.015625, 0.015625, "red", 4, True, 0.2]

    # box_2 = [0.3125, 0.9375, 0.375 - 0.3125, 1.0 - 0.9375, "blue", 4, True, 0.4]

    # box_3 = [0.25, 0.75, 0.0625, 0.0625, "green", 4, True, 0.4]

    # bfa = [box_1]
    # bfa = None

    # time_embedding_v2(pathlib.Path(f"{data_path}/Gene_Data_Sets/Data_Set_{data_set}_histogram.pkl"), 
    time_embedding_v2(pathlib.Path(f"D:\Downloads\GeneData\{species}\Selected{species}Dict/Data_Set_{data_set}_histogram.pkl"), 
                      n = n, 
                      x_lim = x_lim, y_lim = y_lim, 
                      exon_outfile = cwd / "TE_Images_ForPaper" / "Dict" / f"{exon_dict_file}.pkl",
                      intron_outfile = cwd / "TE_Images_ForPaper" / "Dict" / f"{intron_dict_file}.pkl")
    
    matplotfigure(cwd / "TE_Images_ForPaper" / "Dict" / f"{exon_dict_file}.pkl",
                  cwd / "TE_Images_ForPaper" / "Exon",
                  f"{exon_png_file}",
                  default_title = False,
                #   transpose=False,
                #   x_lim=x_lim, y_lim=y_lim,
                  x_title = x_title, y_title = y_title,
                  border= border,
                  x_tick_marks=x_ticks, y_tick_marks=y_ticks,
                  inches = inches,
                  title = f"{species} Exon",
                  box_fill_arrow = bfa,
                  textsize = textsize,
                  dot_size = s, label_coords = label_coords)
    
    matplotfigure(cwd / "TE_Images_ForPaper" / "Dict" / f"{intron_dict_file}.pkl",
                  cwd / "TE_Images_ForPaper" / "Intron",
                  f"{intron_png_file}",
                  default_title = False,
                #   transpose=False,
                #   x_lim=x_lim, y_lim=y_lim,
                  x_title = x_title, y_title = y_title,
                  border= border,
                  x_tick_marks=x_ticks, y_tick_marks=y_ticks,
                  inches = inches,
                  title = f"{species} Intron",
                  box_fill_arrow = bfa,
                  textsize = textsize,
                  dot_size = s, label_coords = label_coords)



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
    # print(scores)

    cheat_sheet = {nucsequence[0]: np.dot(0, w_p), nucsequence[1]: np.dot(1, w_p), nucsequence[2]: np.dot(2, w_p), nucsequence[3]: np.dot(3, w_p)}
    cheat_sheet = pandas.DataFrame(data = cheat_sheet).T
    cheat_sheet.columns = [f"n/{x + 1}" for x in range(k)]

    scores.to_csv(cwd / f"ScoreKey_kmer_{k}_{nucsequence}.csv")
    cheat_sheet.to_csv(cwd / f"CheatSheet_kmer_{k}_{nucsequence}.csv")
    print("Finished Score and Cheat Sheet")



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
                  default_title: bool = True,
                  x_tick_marks: dict = None,
                  y_tick_marks: dict = None,
                  dot_size: float = 0.1,
                  box_fill_arrow: list = None,
                  transpose: bool = True,
                  inches: float = 20,
                  textsize: int = None,
                  border: int = 1,
                  x_title: str = "Forward", y_title: str = "Backward",
                  label_coords: list = None,
                  *args, **kwargs):
    '''
    The figure thing is redundent, so I'm going to try a single function to do everything.

    If frame == pathlib.Path or str then it loads that pickle file into the frame. dir is the output dir.

    box_fill_arrow allows you to draw on the plot. 
    Boxes is a list of the anchor point (xy), the displacement (dxdy), and an optional color (str), linewidth (float), fill (bool), and transparency (float), so box = [x, y, dx, dy, color, width, fill, alpha]. The first 4 are required, the last four are optional and default to None, 2, False, and 0.5.
    
    This section is currenlty not working: you cannot put arrows on the plots.
    Arrows is a list for placing an arrow, with an anchor point (xy), displacement (dxdy), color and linewidth (optional), so arrow = [x, y, dx, dy, color, width].
    Putting them all together, box_fill_arrow should look like box_fill_arrow = [box, fill, arrow].
    End of section that is not working. At some point it might work, but not today.

    Transpose is because I originally put the backwards on the x axis and the forwards on the y axis. Instead of going back and swapping those in other places, it's up to the user to know which ones are which. This script is designed such that the backwards is on the x axis and the forwards is on the y: 
    if you want to do it the other way use Transpose = True (which is true by default), else do false. This flips the coordinates for you.

    ax.patch.set_edgecolor('pink') to adjust the actual edge of the box: make this thicker and darker.

    Remove outer frame of the image.
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

    if transpose:
        for i, (key, points) in enumerate(frame.items()):
            old_x, old_y = np.copy(points[:, 0]), np.copy(points[:, 1])
            frame[key][:, 0] = old_y
            frame[key][:, 1] = old_x

    count = list(frame.keys())[len(list(frame.keys())) - 1]

    if default_title:
        base_title = f"Time Embedding w/ {gap}-mer Gap\n{count} Total Regions"

        if isinstance(title, str):
            title = f"{title}\n{base_title}"
        
        if backwards:
            weights = ": weights are forwards and backwards"
        else:
            weights = ": weights are forwards"

        title = f"{title}{weights}{nucsequence}"

    # Plot
    if isinstance(textsize, int):
        plt.rc("font", size = textsize)

    fig, ax = plt.subplots()
    fig.set_size_inches(inches, inches)
    for i, points in enumerate(frame.values()):
        # ax.scatter(points[:, 0], points[:, 1], s = dot_size, marker = "s", color = "k", linewidths = 0)
        plt.plot(points[:, 0], points[:, 1], markersize = dot_size, marker = "s", linestyle = "", color = "k")
    plt.title(title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    
    if isinstance(box_fill_arrow, list):
        boxes = box_fill_arrow
        arrows = None

        # boxes = box_fill_arrow[0]
        # try:
        #     arrows = box_fill_arrow[2]
        # except IndexError:
        #     arrows = None

        for box in boxes:
            box_color, box_width, box_fill, box_alpha = _colorwidth(box)
            p = plt.Rectangle([box[0], box[1]], box[2], box[3], color = box_color, linewidth = box_width, fill = box_fill, alpha = box_alpha)  #set_fill = False, 
            ax.add_patch(p)

        if isinstance(arrows, list):
            for arrow in arrows:
                arrow_color, arrow_width, _ = _colorwidth(arrow, width = 0.08)
                plt.arrow(arrow[0], arrow[1], arrow[2], arrow[3], width = arrow_color, facecolor = arrow_width, edgecolor = "none", length_includes_head = True)
    
    if isinstance(x_tick_marks, dict):
        ax.set_xticks(list(x_tick_marks.keys()))
        ax.set_xticklabels(list(x_tick_marks.values()))
    if isinstance(y_tick_marks, dict):
        ax.set_yticks(list(y_tick_marks.keys()))
        ax.set_yticklabels(list(y_tick_marks.values()))

    if isinstance(label_coords, list):
        if isinstance(label_coords[0], float) and isinstance(label_coords[1], float):
            ax.xaxis.set_label_coords(label_coords[0], label_coords[1])
        if isinstance(label_coords[2], float) and isinstance(label_coords[3], float):
            ax.yaxis.set_label_coords(label_coords[2], label_coords[3])


    if isinstance(x_tick_marks, dict) and isinstance(y_tick_marks, dict):
        p = plt.Rectangle([list(x_tick_marks.values())[0], list(y_tick_marks.values())[0]], list(x_tick_marks.values())[1] - list(x_tick_marks.values())[0], list(y_tick_marks.values())[1] - list(y_tick_marks.values())[0], color = "black", linewidth = border, fill = False)  #set_fill = False, 
        ax.add_patch(p)

    # axes = ax.axes.flat

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0)

    if (file_extension in ".png") or (file_extension in ".pdf"):
        plt.savefig(f"{file}")
    elif file_extension in ".pkl":
        pickle.dump(fig, open(f"{file}", "wb"))

    print(f"Output image to {file}")
    plt.close()



def _colorwidth(thing: list, width: float = 2, fill: bool = False, alpha: float = 0.5):
    '''
    4 == color, 5 == width, 6 == fill, 7 == alpha/transparency
    '''
    try:
        color = thing[4]
    except IndexError:
        color = None
    try:
        width = thing[5]
    except IndexError:
        width = width
    try:
        fill = thing[6]
    except IndexError:
        fill = fill
    try:
        alpha = thing[7]
    except IndexError:
        alpha = alpha
    
    return color, width, fill, alpha


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


    xy = np.zeros(shape=(seq_length - (k_p + k_m + gap) + 1, 2))

    k_minus = [sequence[k_prime:k_prime + k_m] for k_prime in range(0, seq_length - (k_p + k_m + gap) + 1)]
    k_plus = [sequence[k_prime:k_prime + k_p] for k_prime in range(gap + k_m, seq_length - k_p + 1)]

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

    _check_file(exon_outfile)
    _check_file(intron_outfile)

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
        print(f"Writing to {exon_outfile}")
        pickle.dump(e_frame, e_file)
    with open(intron_outfile, "wb") as i_file:
        print(f"Writing to {intron_outfile}")
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


def _check_file(file_path: pathlib.Path):
    '''
    Checks if the file path is free: if it is not it gives the option to abandon the script. If the it is free then a temp empty dictionary is put there
    in case multiple scripts are running at once.
    '''
    if file_path.is_file():
        print(f"{file_path} already exists!")
        overwrite = input(f"Do you wish to continue [y/n]? ")
        if (overwrite in "n") or (overwrite in "No") or (overwrite in "no") or (overwrite in "NO"):
            print("Ending script")
            exit()
        
        else:
            temp = {}
            with open(file_path, "wb") as file:
                pickle.dump(temp, file)



def _junk():
    '''
    So this stuff at one point was useful, and I don't think it is anymore, but I don't have the guts to delete it. So it goes here in case I do need it.
    '''
    # box_1 = [0.7634, 0.26523, 0.87456 - 0.7634, 0.375 - 0.26523, "red", 4, True, 0.4]

    # box_2 = [0.82389, 0.32356, 0.843 - 0.82389, 0.34467 - 0.32356, "blue", 4, True, 0.4]

    # box_3 = [0.87478, 0.375, 0.985278 - 0.87478, 0.4377 - 0.375, "green", 4, True, 0.4]

    # box_4 = [0.87478, 0.25, 1.0 - 0.87478, 0.4377 - 0.25, "yellow", 4, True, 0.4]

    # box_5 = [0.75, 0.4377, 0.25, 0.5 - 0.4377, "teal", 4, True, 0.4]

    # box_6 = [0.75, 0.375, 0.25, 0.3912 - 0.375, "brown", 4, True, 0.4]

    # box_7 = [0.75, 0.25, 0.25, 0.26523 - 0.25, "violet", 4, True, 0.4]

    # boxes = [box_1, box_2, box_3, box_4, box_5, box_6, box_7]
    # exon_png_file = "ExonZoomed_x075_10_y025_05_7Regions"
    # intron_png_file = "IntronZoomed_x075_10_y025_05_7Regions"



    # box_1 = [0.8125, 0.0, 0.875 - 0.8125, 1.0, "red", 4, True, 0.4]  #[x, y, dx, dy, color, width, fill, alpha]

    # box_2 = [0.25, 0.75, 0.25, 0.25, "green", 4, True, 0.4]

    # box_3 = [0.0, 0.4375, 1.0, 0.5 - 0.4375, "blue", 4, True, 0.4]

    # # box_4 = [0.87478, 0.25, 1.0 - 0.87478, 0.4377 - 0.25, "yellow", 4, True, 0.4]

    # # box_5 = [0.75, 0.4377, 0.25, 0.5 - 0.4377, "teal", 4, True, 0.4]

    # # box_6 = [0.75, 0.375, 0.25, 0.3912 - 0.375, "brown", 4, True, 0.4]

    # # box_7 = [0.75, 0.25, 0.25, 0.26523 - 0.25, "violet", 4, True, 0.4]

    # boxes = [box_1, box_2, box_3]
    # # boxes = None
    # bfa = boxes
    # exon_png_file = "Exon_IterativeMapRGB_x025y075"
    # intron_png_file = "Intron_IterativeMapRGB_x025y075"

    # n = 100_000
    # data_set = 1
    # s = 1
    # data_set = f"1&2"
    # data_set = 1

    # box = [xB, yB, dxB, dyB, "color of box", box thickness, fill (true/false), alpha]
    # fill = [x1, x2, y1, y2, color, alpha]
    # arrow[xA, yA, dxA, dxY, arrow thickness, "color of arrow"]

    
    # box_m = [0.75, 0.25, 0.25, 0.25, "green", 4, True, 0.5]  # major box
    
    # boxp1 = [0.0, 0.8125, 1.0, 0.875 - 0.8125, "blue", 4, True, 0.5]  # forward box
    # boxp2 = [0.0, 0.203125, 1.0, 0.21875 - 0.203125, "blue", 4, True, 0.4]  # forward box
    # boxp3 = [0.0, 0.05078125, 1.0, 0.0546875 - 0.05078125, "blue", 4, True, 0.3]  # forward box
    # boxp4 = [0.0, 0.0126953125, 1.0, 0.013671875 - 0.0126953125, "blue", 4, True, 0.2]  # forward box

    # boxes = [box_1]
    # exon_png_file = "ExonZoomed_x075_10_y025_05_R1"
    # intron_png_file = "IntronZoomed_x075_10_y025_05_R1"

    # boxes = [box_2]
    # exon_png_file = "ExonZoomed_x075_10_y025_05_R2"
    # intron_png_file = "IntronZoomed_x075_10_y025_05_R2"

    # boxes = [box_3]
    # exon_png_file = "ExonZoomed_x075_10_y025_05_R3"
    # intron_png_file = "IntronZoomed_x075_10_y025_05_R3"

    # boxes = [box_4]
    # exon_png_file = "ExonZoomed_x075_10_y025_05_R4"
    # intron_png_file = "IntronZoomed_x075_10_y025_05_R4"

    # boxes = [box_5]
    # exon_png_file = "ExonZoomed_x075_10_y025_05_R5"
    # intron_png_file = "IntronZoomed_x075_10_y025_05_R5"

    # boxes = [box_6]
    # exon_png_file = "ExonZoomed_x075_10_y025_05_R6"
    # intron_png_file = "IntronZoomed_x075_10_y025_05_R6"

    # boxes = [box_7]
    # exon_png_file = "ExonZoomed_x075_10_y025_05_R7"
    # intron_png_file = "IntronZoomed_x075_10_y025_05_R7"

    # bfa = [boxes]
    # bfa = None

    # if isinstance(x_lim, list):
    #     exon_dict_file = f"ExonData_n{n}_DS{data_set}_kp6_km6_zoomed_x{x_lim[0]}by{x_lim[1]}_y{y_lim[0]}by{y_lim[1]}"
    #     intron_dict_file = f"IntronData_n{n}_DS{data_set}_kp6_km6_zoomed_x{x_lim[0]}by{x_lim[1]}_y{y_lim[0]}by{y_lim[1]}"
    #     exon_png_file = f"{exon_dict_file}_IterativeMap.png"
    #     intron_png_file = f"{intron_dict_file}_IterativeMap.png"
    # else:
    #     exon_dict_file = f"ExonData_n{n}_DS{data_set}_kp6_km6"
    #     intron_dict_file = f"IntronData_n{n}_DS{data_set}_kp6_km6"
    #     exon_png_file = f"{exon_dict_file}_IterativeMap_rgb_filled.png"
    #     intron_png_file = f"{intron_dict_file}_IterativeMap_rgb_filled.png"



if __name__ in "__main__":
    main()