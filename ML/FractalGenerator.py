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
    '''
    kmer = 6
    generate_trajectories(pathlib.Path(cwd / "Known_Genes_hg19_NCBIGene_DICT.pkl"), choices = 12)


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



def time_embedding_v3(sequence: str, k_p = 9, k_m = 9, gap = 0, backwards = True, compliment = False):
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