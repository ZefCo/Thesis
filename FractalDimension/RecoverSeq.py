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
    '''
    score: float = float(input("Enter a score: "))
    k: int = int(input("Enter the length of the sequence: "))
    recovered_seq, remainder, score = recover_seq(score, k)

    print(f"Recovered Sequenc: {recovered_seq}\tRemainder = {remainder}\tOriginal Score: {score}")

    # scores = [0.7634, 
    #           0.87456, 
    #           0.26523,  
    #           0.82389, 
    #           0.843, 
    #           0.32356, 
    #           0.34467, 
    #           0.87478, 
    #           0.985278, 
    #           0.375, 
    #           0.4377,
    #           0.26523,
    #           0.3912]
    
    # for score in scores:
    #     seq, remain, score = recover_seq(score, 6)

    #     print(f"Score = {score}\tSeq = {seq}\t\tRemainder = {remain}")





def recover_seq(score: float, k: int, nucsequence: str = "AGTC") -> tuple:
    '''
    takes an input score and "recovers" the original sequence. Iterativly subtracts the possible scores from it and gives back the probably sequence.
    Returns the sequence, the remainder, and the original score.
    '''
    recovred_seq = ""

    w_p = [(0.25)**n for n in range(1, k + 1)]
    cheat_sheet = {nucsequence[0]: np.dot(0, w_p), nucsequence[1]: np.dot(1, w_p), nucsequence[2]: np.dot(2, w_p), nucsequence[3]: np.dot(3, w_p)}
    cheat_sheet = pandas.DataFrame(data = cheat_sheet).T
    cheat_sheet.columns = [x for x in range(k)]
    ecneuqesucn = nucsequence[::-1] # varialbe name is nucsequence backwards

    remainder = score

    for ki in range(k):
        for n in ecneuqesucn:
            nscore = cheat_sheet.loc[n, ki]

            if remainder > nscore:
                remainder -= nscore
                recovred_seq = f"{recovred_seq}{n}"
                break
    
    return recovred_seq, remainder, score



if __name__ in "__main__":
    main()