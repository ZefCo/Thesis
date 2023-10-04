import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import itertools
import os
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from PIL import Image
import NC
import pickle
from pptx import Presentation
from pptx.util import Inches, Pt
import GeneClass as Gene
import TimeEmbedding as KTA
# import pywt
import random
import TimeEmbedding as TE


def main():
    '''
    Start with x = C and y = G. Then we want to see where it came from and where it's going to. Find this for all possible combinations of n. The formula for the x and y coordinates is:

    n_6 n_5 n_4 n_3 n_2 n_1 | n^1 n^2 n^3 n^4 n^5 n^6

    x = sum 1/4^i n_i
    y = sum 1/4^i n^i

    so we would put C and G in the n_6 and n_5 position then do 8 iterations for all possible combinations of nucleotides

    Niether of these are working.

    Trying to find all possible sequences is filling up my memory, and trying to use the equation is yeilding static numbers.
    '''
    # xy = np.array([0.9, 0.3])
    # something = scores(xy, k = 100)
    # print(something)
    # seq = sequences("")
    # print(seq)

    test_sequence = TE.generate_sequence(k = 100)
    print(test_sequence)
    xy = TE.time_embedding(test_sequence)
    # print(xy)
    # print(len(xy))
    # print(type(xy))
    xy_0 = xy[0]
    xy_p = scores(xy_0, k = len(xy))
    e = error(xy, xy_p)

    data = pandas.DataFrame(data = {"X": xy[:, 0], "Y": xy[:, 1], "X'": xy_p[:, 0], "Y'": xy_p[:, 1], "Ex": e[:, 0], "Ey": e[:, 1]})
    print(data)
    data.to_csv(f"Iterative_Process_{test_sequence}.csv")

    # print(e)




def scores(score: np.ndarray, k: int = 4):
    '''
    Does k iterations backwards and 4 iterations forwards

    {} -> decimal/fractional portion
    [] -> integer portion

    y_k+1 = {4y_k}
    x_k+1 = 0.25 * ([4y_k] + x_k)
    '''
    # x = score
    xy = np.zeros(shape = (k, 2))

    xy[0][0] = score[0]
    xy[0][1] = score[1]

    for i in range(1, k):
        x, y = xy[i - 1][0], xy[i - 1][1]

        y_p = (4*y) % 1
        x_p = 0.25 * (int(4*y) + x)
        
        xy[i][0], xy[i][1] = x_p, y_p

    return xy


def sequences(target: str, nucsequence: str = "AGTC", k: int = 11):
    '''
    This is going to find all possible combinations of n_6 n_5 n_4 n_3 n_2 C | G n^2 n^3 n^4 n^5 n^6. Break this up into two seperate chunks of 6 each, and them combine them.

    Ultimatly I want to iterate from n_11 n_10 n_9 n_8 n_7 n_6 | n_5 n_4 n_3 n_2 C G to C G n^2 n^3 n^4 n^5 | n^6 n^7 n^8 n^9 n^10 n^11

    This creates a bunch of sequences with the target thing at the middle. That way I can just iterate through them.

    This creates a memory error because it's something like 10^15 possible sequences... yeah we're just going to create one random one instead.
    '''
    perms = tuple(itertools.product([nucsequence[0], nucsequence[1], nucsequence[2], nucsequence[3]], repeat = k))
    # seq = set()

    berms = _seqs(perms)
    ferms = _seqs(perms)
    totals: int = len(berms) * len(ferms)

    random_b = random.randint(0, len(berms))
    random_f = random.randint(0, len(ferms))

    seq = f"{berms[random_b]}{target}{ferms[random_f]}"

    # for berm in berms:
    #     for ferm in ferms:
    #         sequence = f"{berm}{target}{ferm}"
    #         seq.add(sequence)
    #         # print(len(seq))

    #         if (len(seq) % 1000000) == 0:
    #             print(f"Finished {len(seq)} of {totals}")

    # print(seq[0])
    # print(seq[len(seq) - 1])

    return seq

def _seqs(permutations: list):
    '''
    '''
    perms = []
    for perm in permutations:
        p = ""
        for n in perm:
            p = f"{p}{n}"
        
        perms.append(p)

    return perms



def error(xy_t, xy_m):
    '''
    finds the error between two arrays:

    e = (m - t)/t
    '''

    e = (abs(xy_m - xy_t) / xy_t) * 100

    return e





if __name__ in "__main__":
    main()