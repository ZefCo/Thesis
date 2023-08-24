import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import itertools
import os
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
    generator(f"/media/ethanspeakman/Elements/Gene_Data_Sets/Data_Set_2_histogram.pkl", kmer, ex_col = "exon", int_col="intron", classification_col = "Classificaion")


def generator(data_file: pathlib.Path, kmer: int, 
              method: str = "CGR",
              min_length: int = 100, 
              ex_col: str = "Exon", int_col: str = "Intron", classification_col: str = "Type", gene_name_col: str = "NCIBName",
              *args, **kwargs):
    '''
    This function calls the other functions to generate the Fractals.

    data_file should be a pickled dataframe. I could make a csv option, but I wont.
    '''
    try:
        train_data: pandas.DataFrame = pandas.read_pickle(cwd / data_file)
    except Exception as e:
        print(f"Unable to open file {data_file}\nError Type: {type(e)}\nError {e}\nAre you sure it's a file?")
        exit()

    keep = np.where(train_data["Seq"].str.len() >= min_length)[0]
    train_data = train_data.iloc[keep, :]
    train_data = train_data.reset_index()

    rows, cols = train_data.shape

    exon = 0
    intron = 0

    new_folder: pathlib.Path = pathlib.Path( cwd / str(os.path.splitext(data_file)[0] + f"_{kmer}mer"))
    print(f"Writing training images to {new_folder}")
    try:
        new_folder.mkdir(parents = True, exist_ok = False)
    except FileExistsError as e:
        print("Data Already exists: continuing will delete that data and start over.")
        response: str = input("Do you wish to continue? [y/n] ")

        if response in "y":
            print("Deleting files")
            shutil.rmtree(new_folder)
            # exit()
            new_folder.mkdir(parents = True, exist_ok = False)
        else:
            print("Exiting script")

    except Exception as e:
        print(f"New Exception")
        print(f"Type: {type(e)}")
        print(f"Error: {e}")
        exit()

    exon_images: pathlib.Path = new_folder / "EXON"
    intron_images: pathlib.Path = new_folder / "INTRON"
    exon_images.mkdir(parents = True, exist_ok = True)
    intron_images.mkdir(parents = True, exist_ok = True)
    meta_data = new_folder / f"Meta_Data.txt"


    if method in "CGR":
        fractal_method = chaos_game_representation
        kwargs["k"] = kmer
    elif method in "KTA":
        print("Unable to handle time embedding yet, exiting")
        exit()
        # It's going to take a bit more work. Time Embedding returns a xy vector that needs to be turned into a picture.
        fractal_method = time_embedding_v3
        kwargs["k_p"], kwargs["k_m"] = kmer, kmer

    for row in range(rows):
        seq = train_data.loc[row, "Seq"]
        typ = train_data.loc[row, classification_col]
        gene_name = train_data.loc[row, gene_name_col]
        length = len(seq)

        cgr = fractal_method(seq, **kwargs)

        if typ in int_col:
            filepath = intron_images / f"Intron_{intron}.png"
            try:
                plt.imsave(filepath, cgr, cmap = "gray")
                intron += 1
                with open(meta_data, "at") as mf:
                    mf.write(f"{filepath.name}\t{gene_name}\tExon {intron}\tLength = {length}\n")

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

                        xy = time_embedding_v3(exon, kp, km, 0)
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

                        xy = time_embedding_v3(intron, kp, km, 0)
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


def gene_trajectory(sequence: str, k_p = 9, k_m = 9, *args, **kwargs):
    '''
    Does the time embedding image one at a time.
    '''



def time_embedding_v3(sequence: str, k_p = 9, k_m = 9, gap = 0, backwards = True, compliment = False, *args, **kwargs):
    '''
    Feeds in a sequence, and it finds the xy coordinates for that sequence.

    This is out of date. Do not use
    '''
    seq_length = len(sequence)

    if seq_length < (k_m + k_p + gap):
        print("Cannont find Trajectory for this gene: to small")
        return None

    w_p = [(0.25)**n for n in range(1, k_p + 1)]
    w_m = [(0.25)**n for n in range(1, k_m + 1)]

    if backwards:
        w_m.reverse()

    # b_frame, e_frame, i_frame = [], [], []
    # e_count, i_count = 0, 0

    sequence = sequence.upper()


    xy = np.zeros(shape=(seq_length - (k_p + k_m + gap), 2))

    k_minus = [sequence[k_prime:k_prime + k_m] for k_prime in range(0, seq_length - (k_p + k_m + gap))]
    k_plus = [sequence[k_prime:k_prime + k_p] for k_prime in range(gap + k_m, seq_length - k_p)]

    if compliment:
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


if __name__ in "__main__":
    main()