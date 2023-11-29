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
    I'll keep the below notes, but the idea behind it is that you are tracking the lines, not the points. Really your x,y points are defined from a sum function going from i = 1 to infinity, so we really
    only care about the lines. 


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

    # test_sequence = TE.generate_sequence(k = 100)
    # print(test_sequence)
    # xy = TE.time_embedding(test_sequence)
    # # print(xy)
    # # print(len(xy))
    # # print(type(xy))
    # xy_0 = xy[0]
    # xy_p = scores(xy_0, k = len(xy))
    # e = error(xy, xy_p)

    # data = pandas.DataFrame(data = {"X": xy[:, 0], "Y": xy[:, 1], "X'": xy_p[:, 0], "Y'": xy_p[:, 1], "Ex": e[:, 0], "Ey": e[:, 1]})
    # print(data)
    # data.to_csv(f"Iterative_Process_{test_sequence}.csv")

    # print(e)

    x1 = float(input("Input a value for x1: "))
    x2 = float(input("Input a value for x2: "))
    y1 = float(input("Input a value for y1: "))
    y2 = float(input("Input a value for y2: "))
    k = int(input("Input number of iterations: "))
    X, Y = 0, 1

    # print("Testing with [0.75, 1.0], [0.25. 0.5]")
    # x1 = 0.75
    # x2 = 1.0
    # y1 = 0.25
    # y2 = 0.5

    xy = np.array([[x1, x2], [y1, y2]])
    score_file = f"Scores_x_{x1}_{x2}_y_{y1}_{y2}"
    score_file = score_file.replace(".", "")
    score_file = cwd / "TE_Images_ForPaper" / "Scores" / f"{score_file}.csv"

    # xy = scores(xy)
    # print(f"Backwards propogation:\n\tx1 = {xy[0][X][0]}\tx2 = {xy[0][X][1]}\n\ty1 = {xy[0][Y][0]}\ty2 = {xy[0][Y][1]}")
    # print(f"Original position:\n\tx1 = {xy[1][X][0]}\tx2 = {xy[1][X][1]}\n\ty1 = {xy[1][Y][0]}\ty2 = {xy[1][Y][1]}")
    # print(f"Forwards propogation:\n\tx1 = {xy[2][X][0]}\tx2 = {xy[2][X][1]}\n\ty1 = {xy[2][Y][0]}\ty2 = {xy[2][Y][1]}")

    xy = kscores(xy, k)
    xy = convert2frame(xy)
    xy.to_csv(score_file)
    print(f"Printed to {score_file}")


def score_propogation(x1: float, x2: float):
    '''
    '''
    x1p = (4*x1) % 1

    x2p = 0.25 * (int(4*x1) + x2)

    return x1p, x2p


def scores(score: np.ndarray) -> np.ndarray:
    '''
    Finds the forward and backward propogation of a given input. Expects the input to be a 2x2 array where the first row is the x and the second row is the y. Will output a 3x2x2 matrix representing the -1, 0, +1 iterations. 0 is the input.

    {} -> decimal/fractional portion
    [] -> integer portion

    x_k+1 = {4x_k}
    y_k+1 = 0.25 * ([4x_k] + y_k)

    y_k-1 = {4y_k}
    x_k-1 = 0.25 * ([4y_k] + x_k)
    '''
    X, Y  = 0, 1
    xy = np.zeros(shape = (3, 2, 2))

    xy[1][X] = score[X]
    xy[1][Y] = score[Y]


    # # backwards propogation
    y, x = _propogation(xy[1][Y], xy[1][X])
    xy[0][X][0], xy[0][X][1] = x[0], x[1]
    xy[0][Y][0], xy[0][Y][1] = y[0], y[1]

    # forwards propogation
    x, y = _propogation(xy[1][X], xy[1][Y])
    xy[2][X][0], xy[2][X][1] = x[0], x[1]
    xy[2][Y][0], xy[2][Y][1] = y[0], y[1]

    return xy

def kscores(score: np.ndarray, k: int) -> np.ndarray:
    '''
    Finds 2*k scores from an inital score, returning 2*k+1 scores.
    '''
    # print(f"Starting at {k}")

    X, Y  = 0, 1
    xy = np.zeros(shape = (2*k + 1, 2, 2))

    xy[k][X] = score[X]
    xy[k][Y] = score[Y]

    # print(f"Back propogation")
    for i in range(k - 1, -1, -1):
        y, x = _propogation(xy[i + 1][Y], xy[i + 1][X])
        xy[i][X][0], xy[i][X][1] = x[0], x[1]
        xy[i][Y][0], xy[i][Y][1] = y[0], y[1]


    # print(f"Froward propogation")
    for j in range(k + 1, 2*k + 1, 1):
        x, y = _propogation(xy[j - 1][X], xy[j - 1][Y])
        xy[j][X][0], xy[j][X][1] = x[0], x[1]
        xy[j][Y][0], xy[j][Y][1] = y[0], y[1]

    return xy


def convert2frame(array: np.ndarray) -> pandas.DataFrame:
    '''
    Input an array that is kx2x2 and output a kx4 dataframe
    '''
    X, Y = 0, 1
    i, _, _ = array.shape
    
    zeros: np.ndarray = np.zeros(shape = (i, 4))
    index = [j - int(i/2) for j in range(i)]
    frame: pandas.DataFrame = pandas.DataFrame(zeros, columns=["X1", "X2", "Y1", "Y2"], index=index)

    for j in range(i):
        xy = array[j]
        frame.iloc[j, 0], frame.iloc[j, 1] = xy[X][0], xy[X][1]
        frame.iloc[j, 2], frame.iloc[j, 3] = xy[Y][0], xy[Y][1]

    return frame


def _propogation(x: np.ndarray, y: np.ndarray) -> tuple:
    '''
    The actual calculation. I felt it was better as a function and a bit more complicated then a lambda function. To do the backwards propogation, simply switch the x and y inputs.
    Each x and y are a 1x2 array.
    '''
    x_k = np.zeros(shape = (2, 1))
    y_k = np.zeros(shape = (2, 1))

    x_k[0] = 4*x[0] % 1
    x_k[1] = 4*x[1] % 1

    y_k[0] = 0.25 * (int(4*x[0]) + y[0])
    y_k[1] = 0.25 * (int(4*x[0]) + y[1])
    # y_k[1] = 0.25 * (int(4*x[1]) + y[1]) if (x_k[1] > 0) else 0.25 * (int(4*x[0]) + y[1])

    if (x_k[0] == 0) and (x_k[1] == 0):
        x_k[1] = 1.0

    return x_k, y_k


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