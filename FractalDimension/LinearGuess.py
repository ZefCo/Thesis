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
from statsmodels.graphics.tsaplots import plot_acf
import pandas



def main():
    '''
    '''
    plotlog = True
    plotboth = False
    autocorplot = False
    dot_size = 2
    line_width = 0.25

    random_list = walk_directory(str(cwd / "GenePerCrhome"))

    



def walk_directory(folder: str or pathlib.Path) -> list:
    '''
    '''

    file: str
    files = list()
    # print(folder)
    # exit()

    for r, _, f in os.walk(str(folder)):
        for file in f:
            if file.endswith(".npy"):
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