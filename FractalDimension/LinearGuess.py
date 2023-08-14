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
    linear_approximation(dict_list, show_figs = False, method = "manhattan")
    # F_numbers(dict_list)


def load_plot(list_of_files: list):
    '''
    Loads the pickle file and puts it into a plot.
    '''

    for file in list_of_files:
        with open(file, "rb") as p:
            data = pickle.load(p)

        CD.correlation_dimension_figure_v2(data, title = str(file))

        # break



def linear_approximation(list_of_files: list, std_max: float = 0.4, std_min: float = 0.15, slope_thresh = 0.002, show_figs: bool = False, *args, **kwargs):
    '''
    Generates a linear approximation to the log plot. First it compuetes the distance between all points on the plot. Then it discards all points above and below a
    threshold. From those points a linear approximation is made based on the least squares method of the remaining points.

    The threshold for triming is based on the std of the distance. By default it is set to +0.4*std & -0.15*std. These can be adjusted manually. This does a fairly good job of triming
    the wider spots in the data, but not so much the later points at saturation.
    '''
    file: pathlib.Path
    per_thresh = lambda x1, x2: abs(x1 - x2) # / (0.5*(x1 + x2))

    data_keys = ["Gene", "Region Type", "Region Index", "Length", "Slope", "Y Intercept", "R**2", "Image Size", "F Number"]
    meta_data = {data_keys[0]: [], 
                 data_keys[1]: [], 
                 data_keys[2]: [], 
                 data_keys[3]: [], 
                 data_keys[4]: [], 
                 data_keys[5]: [], 
                 data_keys[5]: [], 
                 data_keys[6]: [], 
                 data_keys[7]: [], 
                 data_keys[8]: []}

    for file in list_of_files:
        print(file.name)
        with open(file, "rb") as p:
            data: dict = pickle.load(p)

        x_values, y_values = tuple(data.keys()), tuple(data.values())

        data = [[x, y_values[i]] for i, x in enumerate(x_values)]

        data: np.array
        data_og: np.array
        data = data_og = np.array(data)

        DData = Dist.DistanceMatrix(data, nearest_neighbor = True, *args, **kwargs)
        DData.F_number()

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
        r = r_squared(data, ls)
        file_info = re.split("_", file.name)  # 0 = gene, 1 = region type, 2 = region index, 3 = length, 4 = image size + file extension (.pkl)

        meta_data[data_keys[0]].append(file_info[0])
        meta_data[data_keys[1]].append(file_info[1]) 
        meta_data[data_keys[2]].append(file_info[2]) 
        meta_data[data_keys[3]].append(file_info[3]) 
        meta_data[data_keys[4]].append(m)
        meta_data[data_keys[5]].append(b)
        meta_data[data_keys[6]].append(r)
        meta_data[data_keys[7]].append(re.split(r"\.", file_info[4])[0])
        meta_data[data_keys[8]].append(DData.F[DData.F.shape[0] - 1])

        if show_figs:
            fig = go.Figure()
            fig.add_trace(go.Scatter(x = data_og[:,0], y = data_og[:,1], name = "Original Data"))
            fig.add_trace(go.Scatter(x = data[:,0], y = data[:,1], name = "Selected Data"))
            fig.add_trace(go.Scatter(x = tuple(ls.keys()), y = tuple(ls.values()), name = f"Interpolated Lines<br>y = {round(m, 3)}x + {round(b, 3)}<br>R**2 = {r}"))
            fig.update_layout(title = f"File:<br>{file}")
            fig.show()

        # exit()

    meta_data = pandas.DataFrame(meta_data)

    meta_data.to_csv(cwd / "GenePerChrom_META.csv")
    meta_data.to_pickle(cwd / "GenePerChrom_META.pkl")



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


def dict_to_numpy(input_dict: dict) -> np.array:
    '''
    Converst a dictionary to a 2D numpy array. Assumes the dict is structured as {X: Y} and will return the array as [[X Y]] 
    '''
    x_values, y_values = tuple(input_dict.keys()), tuple(input_dict.values())
    input_dict = [[x, y_values[i]] for i, x in enumerate(x_values)]
    
    return_numpy = np.array(input_dict)

    return return_numpy


def total_variation(y_values: np.array):
    '''
    Probably could be a lambda function. Finds the mean of the values, subtracts that from the actual values, squrars it, and returns that value.
    '''
    var = y_values - y_values.mean()
    var = np.dot(var, var)
    
    return var


def error(y_values: np.array, y_predict: np.array):
    '''
    This could also probably be a lambda function. Finds the difference between the actual values and the predicted values and squares them.
    '''
    e = (y_values - y_predict)
    e = np.dot(e, e)

    return e




def r_squared(actual_data: np.array or dict, predicted_data: np.array or dict) -> float:
    '''
    finds the R Square of the data.
    '''
    if isinstance(actual_data, dict):        
        actual_data: np.array = dict_to_numpy(actual_data)

    if isinstance(predicted_data, dict):
        predicted_data: np.array = dict_to_numpy(predicted_data)

    if actual_data.shape == predicted_data.shape:
        # print("\tYep")
        e = error(actual_data[:,1], predicted_data[:,1])
        v = total_variation(actual_data[:,1])

        r = 1 - (e / v)
        
        return r

    else:
        # print("\tNope")

        return None
    

def F_numbers(list_of_files: list):
    '''
    '''
    file: pathlib.Path
    per_thresh = lambda x1, x2: abs(x1 - x2) # / (0.5*(x1 + x2))

    data_keys = ["Gene", "Region Type", "Region Index", "Length", "Slope", "Y Intercept", "R**2", "Image Size"]
    meta_data = {data_keys[0]: [], data_keys[1]: [], data_keys[2]: [], data_keys[3]: [], data_keys[4]: [], data_keys[5]: [], data_keys[5]: [], data_keys[6]: [], data_keys[7]: []}

    for file in list_of_files:
        print(file.name)
        with open(file, "rb") as p:
            data: dict = pickle.load(p)

        x_values, y_values = tuple(data.keys()), tuple(data.values())

        data = [[x, y_values[i]] for i, x in enumerate(x_values)]
        data = np.array(data)

        DData = Dist.DistanceMatrix(data, nearest_neighbor = True, method = "euclidean")
        DData.F_number()

        print(f"\t{DData.F[DData.F.shape[0] - 1]}")

        # exit()





if __name__ in "__main__":
    main()