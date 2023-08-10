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
import CorrelationDimension as CD



def main():
    '''
    '''
    # plotlog = True
    # plotboth = False
    # autocorplot = False
    # dot_size = 2
    # line_width = 0.25

    dict_list = walk_directory(str(cwd / "Distance_Dict"), extension=".pkl")
    load_plot(dict_list)
    # TESTFILE = cwd / "Distance_Dict" / "EXON" / "TUBGCP5_E_23_587_4096x4096.pkl"

    # with open(TESTFILE, "rb") as p:
    #     test_file = pickle.load(p)

    # print(test_file)


def load_plot(list_of_files: list):
    '''
    Loads the pickle file and puts it into a plot.
    '''

    for file in list_of_files:
        with open(file, "rb") as p:
            data = pickle.load(p)

        CD.correlation_dimension_figure_v2(data, title = str(file))

        # break




def walk_directory(folder: str or pathlib.Path, extension = ".npy") -> list:
    '''
    Simply walks thorugh the directory.

    I should really build a common script file for things like this. I use this type of thing over and over again and am tired of rewriting it.
    '''

    file: str
    files = list()
    # print(folder)
    # exit()

    for r, _, f in os.walk(str(folder)):
        for file in f:
            if file.endswith(extension):
                files.append(pathlib.Path(os.path.join(r, file)))

    return files


def least_squares(x_points: np.array, y_points: np.array):
    '''
    '''
    if isinstance(x_points, list):
        x_points = np.array(x_points)
    if isinstance(y_points, list):
        y_points = np.array(y_points)

    print(y_points)

    ave_x = np.mean(x_points)
    ave_y = np.mean(y_points)

    # print(x_points)
    # print(x_points.max())
    # print(y_points)

    m = np.dot(x_points - ave_x, y_points - ave_y) / np.dot(x_points - ave_x, x_points - ave_x)
    b = ave_y - (m * ave_x)
    print(ave_x, ave_y, m, b)

    ls = {e: m*e + b for e in x_points}

    return ls, m, b


if __name__ in "__main__":
    main()