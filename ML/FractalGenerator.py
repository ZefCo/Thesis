import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import itertools
import os
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image
import pickle
import GeneClass as Gene
import random
import shutil


def main():
    '''
    This is going to be the main script for generating all fractals from now on.
    '''

    kmer = 6
    min_length = 40
    method = "KTA"
    data_set = 1
    ex_col: str = "exon"
    int_col: str = "intron"
    classification_col: str = "Classificaion"
    gene_name_col: str = "NCIBName"

    generator(f"G:/Gene_Data_Sets/Data_Set_{data_set}_histogram.pkl", kmer, target_dir = f"G:/Gene_Data_Sets/Data_Set_{data_set}_histogram_{kmer}mer_{method}", ex_col = ex_col, int_col = int_col, classification_col = classification_col, method = method, min_length = min_length, gene_name_col = gene_name_col)


def generator(data_file: pathlib.Path, kmer: int,
              target_dir: pathlib.Path,
              method: str = "CGR",
              min_length: int = 100, 
              *args, **kwargs):
    '''
    This function calls the other functions to generate the Fractals.

    data_file should be a pickled dataframe. I could make a csv option, but I wont.
    '''
    if isinstance(data_file, str):
        data_file = pathlib.Path(data_file)
    if isinstance(target_dir, str):
        target_dir = pathlib.Path(target_dir)

    try:
        train_data: pandas.DataFrame = pandas.read_pickle(data_file)
    except Exception as e:
        print(f"Unable to open file {data_file}\nError Type: {type(e)}\nError {e}\nAre you sure it's a file?")
        exit()

    keep = np.where(train_data["Seq"].str.len() >= min_length)[0]
    train_data = train_data.iloc[keep, :]
    train_data = train_data.reset_index()

    print(f"Writing training images to {target_dir}")
    try:
        target_dir.mkdir(parents = True, exist_ok = False)
    except FileExistsError as e:
        print("Data Already exists: continuing will delete that data and start over.")
        response: str = input("Do you wish to continue? [y/n] ")

        if response in "y":
            print("Deleting files")
            shutil.rmtree(target_dir)
            # exit()
            target_dir.mkdir(parents = True, exist_ok = False)
            print("Files Deleted, resuming script")
        else:
            print("Exiting script")
            exit()

    except Exception as e:
        print(f"New Exception")
        print(f"Type: {type(e)}")
        print(f"Error: {e}")
        exit()

    exon_images: pathlib.Path = target_dir / "EXON"
    intron_images: pathlib.Path = target_dir / "INTRON"
    exon_images.mkdir(parents = True, exist_ok = True)
    intron_images.mkdir(parents = True, exist_ok = True)
    meta_data = target_dir / f"Meta_Data.txt"
    kwargs["meta_data"] = meta_data

    if method in "CGR":
        kwargs["k"] = kmer
        kwargs["target_dir"] = target_dir
        # print(kwargs)
        cgr_generator(train_data, *args, **kwargs)
    elif method in "KTA":
        # print("Unable to handle time embedding yet, exiting")
        # exit()
        # It's going to take a bit more work. Time Embedding returns a xy vector that needs to be turned into a picture.
        matplotlib.rc('figure', figsize=(1, 1), dpi = np.sqrt(4**kmer))
        kwargs["k_p"], kwargs["k_m"] = kmer, kmer
        kwargs["target_dir"] = target_dir
        # print(kwargs)
        kta_generator(train_data, *args, **kwargs)




def manhattan_position(nuc: int, x0: int):
    '''
    '''
    x1: np.array = (nuc - x0) / 2

    return x1.astype(int)


def nucleotide_permutations(sequence: str = "ACGT", length: int = 3) -> dict:
    nuc_perm = dict()

    if len(sequence) < length:
        return None

    perms = itertools.permutations(sequence, length)
    for p in perms:

        key = ""
        for n in p:
            key = f"{key}{n}"

        nuc_perm[key] = 0

    return nuc_perm


def nucleotide_counter(sequence: str, window_size: int):
    '''
    '''
    keys: set = set()
    counter = dict()
    master_count = 0


    for i in range(len(sequence) - window_size):
        seq = sequence[i: i + window_size]

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

def kta_generator(train_data: pandas.DataFrame,
                  target_dir: pathlib.Path = None,
                  int_col: str = None, ex_col: str = None,
                  classification_col: str = None, gene_name_col: str = None,
                  *args, **kwargs):
    '''
    for when the method == KTA
    '''

    exon, intron = 0, 0
    exon_images = target_dir / "EXON"
    intron_images = target_dir / "INTRON"

    meta_data = target_dir / f"Meta_Data.txt"

    rows, cols = train_data.shape
    
    for row in range(rows):
        seq = train_data.loc[row, "Seq"]
        typ = train_data.loc[row, classification_col]
        gene_name = train_data.loc[row, gene_name_col]
        length = len(seq)

        xy = time_embedding(seq, gap = 0, *args, **kwargs)
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_axis_off()
        plt.scatter(xy[:, 0], xy[:, 1], color = "black", s = 0.01)
        
        if typ in int_col:
            filepath = intron_images / f"Intron_{intron}.png"
            intron += 1
            try:
                plt.savefig(filepath)
                plt.close()
                with open(meta_data, "at") as mf:
                    mf.write(f"{filepath.name}\t{gene_name}\Intron {intron}\tLength = {length}\n")

            except Exception as e:
                print(type(e))
                print(e)

        if typ in ex_col:
            filepath = exon_images / f"Exon_{exon}.png"
            exon += 1
            try:
                plt.savefig(filepath)
                plt.close()
                with open(meta_data, "at") as mf:
                    mf.write(f"{filepath.name}\t{gene_name}\tExon {intron}\tLength = {length}\n")

            except Exception as e:
                print(type(e))
                print(e)

    
def cgr_generator(train_data: pandas.DataFrame, 
                  classification_col: str = None, gene_name_col: str = None, 
                  int_col: str = None, ex_col: str = None, 
                  target_dir: pathlib.Path = None,
                  meta_data: pathlib.Path = None, 
                  *args, **kwargs):
    '''
    For when the method == CGR
    '''
    rows, cols = train_data.shape
    exon, intron = 0, 0
    intron_images = target_dir / "INTRON"
    exon_images = target_dir / "EXON"

    for row in range(rows):
        seq = train_data.loc[row, "Seq"]
        typ = train_data.loc[row, classification_col]
        gene_name = train_data.loc[row, gene_name_col]
        length = len(seq)

        cgr = chaos_game_representation(seq, **kwargs)

        if typ in int_col:
            filepath = intron_images / f"Intron_{intron}.png"
            try:
                plt.imsave(filepath, cgr, cmap = "gray")
                intron += 1
                with open(meta_data, "at") as mf:
                    mf.write(f"{filepath.name}\t{gene_name}\Intron {intron}\tLength = {length}\n")

            except Exception as e:
                print(type(e))
                print(e)

        elif typ in ex_col:
            filepath = exon_images / f"Exon_{exon}.png"
            try:
                plt.imsave(filepath, cgr, cmap = "gray")
                exon += 1
                with open(meta_data, "at") as mf:
                    mf.write(f"{filepath.name}\t{gene_name}\tExon {intron}\tLength = {length}\n")

            except Exception as e:
                print(type(e))
                print(e)

    # print(filepath)



def chaos_game_representation(seq: str, k: int, *args, **kwargs):
    array_size = int(np.sqrt(4**k))

    cgr = np.zeros(shape = (array_size, array_size))

    probabilities: dict = nucleotide_counter(seq, k)

    for key, value in probabilities.items():
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
 
        for char in key:
            if char == "T":
                posx += maxx / 2
            elif char == "C":
                posy += maxy / 2
            elif char == "G":
                posx += maxx / 2
                posy += maxy / 2
            maxx /=  2
            maxy /= 2

        # print(int(posy), int(posx))
        cgr[int(posy - 1)][int(posx - 1)] = value


    return cgr


def generate_trajectories(data_file: pathlib.Path, choices: int = 5000, kp: int = 6, km: int = 6):
    '''
    Want to grab genes based on the histogram model (find range of exons, get introns with same range, do this by gene) and generate the images from that.
    '''


    gene: Gene.Gene
    exon_count, intron_count = 0, 0

    fractal_dir = cwd / "FractalTrajectory"
    exon_dir = fractal_dir / "EXON"
    intron_dir = fractal_dir / "INTRON"

    try:
        fractal_dir.mkdir(parents = True, exist_ok = False)
    except FileExistsError as e:
        print("Data Already exists: continuing will delete that data and start over.")
        response: str = input("Do you wish to continue? [y/n] ")

        if response in "y":
            print("Deleting files")
            shutil.rmtree(fractal_dir)
            # exit()
            fractal_dir.mkdir(parents = True, exist_ok = False)

        else:
            print("Exiting script")
    except Exception as e:
        print(type(e))
        print(e)
        exit()

    exon_dir.mkdir(parents = True, exist_ok = True)
    intron_dir.mkdir(parents = True, exist_ok = True)

    meta_data = fractal_dir / f"Meta_Data.txt"

    with open(data_file, "rb") as p:
        data: dict = pickle.load(p)

    rand_selection = set(data.keys())
    rand_selection = random.sample(list(rand_selection), k = choices)
    # print(len(rand_selection))
    # print(len(set(rand_selection)))
    # # print(rand_selection)
    # exit()

    for choice in rand_selection:
        gene = data[choice]
        # print(ncibname, gene.exon_seq[0])
        min_len, max_len = 0, 0

        # print(type(gene.exon_seq), type(gene.intron_seq))
        if (gene.exon_seq is not None) and (gene.intron_seq is not None):

            if len(gene.exon_seq) > 0:
                min_len, max_len, local_exon = 20, 100, 0
                for exon in gene.exon_seq:
                    local_exon += 1

                    length = len(exon)
                    if length > min_len:
                        if length > max_len:
                            max_len = length

                        xy = time_embedding(exon, kp, km, 0)
                        # print(xy)
                        # exit()
                        filepath = exon_dir / f"Exon_{exon_count}.png"
                        # plt.imshow(xy, cmap = "gray")
                        plt.scatter(xy[:, 0], xy[:, 1], color = "black")
                        # ax = plt.axes()
                        # ax.set_facecolor("black")
                        # plt.imsave(filepath, xy, cmap = "gray")
                        plt.savefig(filepath)
                        exon_count += 1
                        plt.close()

                        with open(meta_data, "at") as mf:
                            mf.write(f"{filepath.name}\t{choice}\tExon {local_exon}\tLength = {length}\n")
                            
            if ((max_len > 0) and (min_len > 0) and (not None)):
                local_intron = 0
                # print(gene.intron_seq)
                for intron in gene.intron_seq:
                    local_intron += 1
                    length = len(intron)

                    if (max_len >= length) and (length >= min_len):

                        xy = time_embedding(intron, kp, km, 0)
                        filepath = intron_dir / f"Intron_{intron_count}.png"
                        # plt.imshow(xy, cmap = "gray")
                        plt.scatter(xy[:, 0], xy[:, 1], color = "black")
                        # ax = plt.axes()
                        # ax.set_facecolor("black")
                        # plt.imsave(filepath, xy, cmap = "gray")
                        plt.savefig(filepath)
                        intron_count += 1
                        plt.close()

                        with open(meta_data, "at") as mf:
                            mf.write(f"{filepath.name}\t{choice}\tIntron {local_intron}\tLength = {length}\n")




def time_embedding(sequence: str, 
                      k_p: int = 6, k_m: int = 6, gap: int = 0, 
                      m_backwards: bool = True, p_backwards: bool = False, 
                      nucsequence: str = "AGTC",
                      *args, **kwargs):
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

    for i, k_prime in enumerate(k_minus):
        n = [0 if n in nucsequence[0] else 1 if n in nucsequence[1] else 2 if n in nucsequence[2] else 3 if n in nucsequence[3] else 100 for n in k_prime]
        k_x = np.dot(w_m, n)

        n = [0 if n in nucsequence[0] else 1 if n in nucsequence[1] else 2 if n in nucsequence[2] else 3 if n in nucsequence[3] else 100 for n in k_plus[i]]
        k_y = np.dot(w_p, n)

        xy[i][0], xy[i][1] = k_x, k_y

    return xy


if __name__ in "__main__":
    main()