import pandas
import numpy as np
import os
import pathlib
cwd = pathlib.Path.cwd()
import MomentCalculations as MC
import Heatmaps as heat
import TimeEmbedding as TE
import pickle
import GeneClass as Gene
from typing import Tuple
from typing import Generator


def main():
    '''
    '''

    # text = read_file(cwd.parent / "Data_Files/Intergenomic/Seq2Check/Good/NM_000141.5_NM_001001976.3")
    # # print(type(text))
    # # print(len(text))
    # # for line in text:
    # #     print(line)
    # print(type(text[6]))
    # for file in get_files(cwd.parent / "Data_Files/Intergenomic/Seq2Check/Good"):
    #     print(file)

    step = 100
    max_n = 2
    min_n = 0.5
    ms = [n / step for n in range(int(step*(min_n)), int(step*max_n) + 1)]

    data = collect_sequences(directory = cwd.parent / "Data_Files/Intergenomic/Seq2Check/Good")
    n, _ = data.shape

    TE.time_embedding_v2(data = data, n = n, k_p = 6, k_m = 6, inter_outfile = cwd.parent / "Data_Files/Intergenomic/Dicts/InterTE.pkl")
    TE.matplotfigure(cwd.parent / "Data_Files/Intergenomic/Dicts/InterTE.pkl", cwd / "TE_Images_ForPaper" / "InterGenome", "Inter_v1.png", title = "InterGenomic")

    master, *_ = heat.heat_embedding(data, n = n, skip_ExonIntron = True)
    with open(cwd / "HS_Dicts_EF" / "Inter_6mer.pkl", "wb") as file:
        pickle.dump(master, file)

    me, mi, mn, uni, _ = MC.moments_v3(cwd / "HS_Dicts_EF", ms, min_k=6, N_value = True)
    
    me = {"Human": {"data": me, "marker": None}}
    mi = {"Human": {"data": mi, "marker": None}}
    mn = {"Human": {"data": mn, "marker": None}}

    MC.multiple_species_plots(ms, me, mi, mn, uni, cwd / "TE_Images_ForPaper" / "InterGenome", x_ticks={0.5: 0.5, 1.0: 1.0, 1.5: 1.5, 2.0: 2.0},  y_ticks = {0: 0, 7: 7})




def collect_sequences(*args, **kwargs) -> pandas.DataFrame:
    '''
    Returns a Dataframe of: NCIBName (which is none), Classificaion (which is "Inter"), and the sequences.
    '''
    sequences = pandas.DataFrame(columns = ["NCIBName", "Classificaion", "Seq"])
    for f, file in enumerate(get_files(*args, **kwargs)):
        text = read_file(file, *args, **kwargs)
        sequence = text[6]
        sequences.loc[f, "NCIBName"] = "None"
        sequences.loc[f, "Classificaion"] = "Inter"
        sequences.loc[f, "Seq"] = sequence

    return sequences



def get_files(directory: pathlib.Path, *args, **kwargs):
    '''
    Iterates through a directory and gets the names and paths of all the files in that directory.

    It actually yeilds them, this is a generator. Make sure all the files in the path are intended to be used.
    '''
    for file in directory.rglob("*"):
        yield file



def read_file(file: pathlib.Path, *args, **kwargs) -> str:
    '''
    Reads the file. Neede the filepath, after that it reads and finds the 7th line: that is the start of the sequence. Returns a string, the sequence.
    '''
    with open(file, "r") as text:
        text = text.readlines()

    return text


if __name__ in "__main__":
    main()