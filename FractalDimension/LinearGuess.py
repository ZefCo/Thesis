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
import DistanceClass as Dist



def main():
    '''
    '''
    # plotlog = True
    # plotboth = False
    # autocorplot = False
    # dot_size = 2
    # line_width = 0.25

    dict_list = walk_directory(str(cwd / "Distance_Dict"), extension=".pkl")
    # load_plot(dict_list)
    # TESTFILE = cwd / "Distance_Dict" / "EXON" / "TUBGCP5_E_23_587_4096x4096.pkl"

    # with open(TESTFILE, "rb") as p:
    #     test_file = pickle.load(p)

    # print(test_file)
    linear_approximation(dict_list)


def load_plot(list_of_files: list):
    '''
    Loads the pickle file and puts it into a plot.
    '''

    for file in list_of_files:
        with open(file, "rb") as p:
            data = pickle.load(p)

        CD.correlation_dimension_figure_v2(data, title = str(file))

        # break



def linear_approximation(list_of_files: list, std_max: float = 0.4, std_min: float = 0.15, slope_thresh = 0.002):
    '''
    Generates a linear approximation to the log plot. First it compuetes the distance between all points on the plot. Then it discards all points above and below a
    threshold. From those points a linear approximation is made based on the least squares method of the remaining points.

    The threshold for triming is based on the std of the distance. By default it is set to +0.4*std & -0.15*std. These can be adjusted manually. This does a fairly good job of triming
    the wider spots in the data, but not so much the later points at saturation.
    '''
    per_thresh = lambda x1, x2: abs(x1 - x2) # / (0.5*(x1 + x2))

    for file in list_of_files:
        with open(file, "rb") as p:
            data: dict = pickle.load(p)

        x_values, y_values = tuple(data.keys()), tuple(data.values())

        data = [[x, y_values[i]] for i, x in enumerate(x_values)]

        data: np.ndarray
        data_og: np.ndarray
        data = data_og = np.array(data)

        DData = Dist.DistanceMatrix(data, nearest_neighbor = True, method = "manhattan")

        # min_d = DData.D1.min()
        # max_d = DData.D1.max()
        mean_d = DData.D1.mean()
        std_d = DData.D1.std()

        min_thresh = mean_d - (std_min * std_d)
        max_thresh = mean_d + (std_max * std_d)
        something = DData.D1  # I need to call this something other then someting.
        points = (max_thresh <= something) | (something <= min_thresh)  # note the direction of the shapes: we want the points that are outside of our desired threshold
        something[points] = 0  # setting those points to be zero. Also I should call this variable something other then something.

        something = np.transpose(np.transpose(np.where(something > 0)))[0] # grabing the points that are greater then 0: that's where the important distances are
        data = data[something]  # actually pulling those indicies out

        something = np.flip(data, axis = 0)

        for i, xy in enumerate(something):
            if ((per_thresh(xy[1], something[i + 1][1])) < slope_thresh):
                data = np.delete(data, data.shape[0] - 1, axis = 0)
            else:
                break

        ls, m, b = least_squares(data[:, 0], data[:, 1])

        fig = go.Figure()
        fig.add_trace(go.Scatter(x = data_og[:,0], y = data_og[:,1], name = "Original Data"))
        fig.add_trace(go.Scatter(x = data[:,0], y = data[:,1], name = "Selected Data"))
        fig.add_trace(go.Scatter(x = tuple(ls.keys()), y = tuple(ls.values()), name = f"Interpolated Lines<br>y = {round(m, 3)}x + {round(b, 3)}"))
        fig.update_layout(title = f"File:<br>{file}")
        fig.show()

        # exit()



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
    Returns a y = mx + b plot as a tuple. (interpolated x and y as a dict, slope, y intercept)
    '''
    if isinstance(x_points, list):
        x_points = np.array(x_points)
    if isinstance(y_points, list):
        y_points = np.array(y_points)

    # print(y_points)

    ave_x = np.mean(x_points)
    ave_y = np.mean(y_points)

    # print(x_points)
    # print(x_points.max())
    # print(y_points)

    m = np.dot(x_points - ave_x, y_points - ave_y) / np.dot(x_points - ave_x, x_points - ave_x)
    b = ave_y - (m * ave_x)
    # print(ave_x, ave_y, m, b)

    ls = {e: m*e + b for e in x_points}

    return ls, m, b


if __name__ in "__main__":
    main()