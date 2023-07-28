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


def main():
    '''
    '''
    # doeverything(plotlog = True, plotboth = False, autocorplot = False, dot_size = 2, line_width = 0.25)
    doExonIntron()



def doExonIntron():
    '''
    '''
    gene_folder = cwd / "GenePerChrome"
    cgr: np.ndarray

    exon_files, intron_files = walk_exon_intron(gene_folder)

    D_folder = cwd / "Distance_Dict"
    E_folder = D_folder / "EXON"
    I_folder = D_folder / "INTRON"

    E_folder.mkdir(parents = True, exist_ok = True)
    I_folder.mkdir(parents = True, exist_ok = True)

    for file in exon_files:
        cgr = load_cgr(file)

        if isinstance(cgr, np.ndarray):
            L, W = cgr.shape
            new_file = f"{file.stem}_{L}x{W}"
            new_file = pathlib.Path(new_file).with_suffix(".pkl")

            start = timeit.default_timer()
            _, local_c = correlation_dimension(cgr)
            log_c = log_plot(local_c, normalize = False)

            with open(str(E_folder / new_file), "wb") as f:
                pickle.dump(log_c, f)

            stop = timeit.default_timer()

            print(f"Time = {stop - start}")
            # exit()

    for file in intron_files:
        cgr = load_cgr(file)

        if isinstance(cgr, np.ndarray):
            L, W = cgr.shape
            new_file = f"{file.stem}_{L}x{W}"
            new_file = pathlib.Path(new_file).with_suffix(".pkl")

            start = timeit.default_timer()
            _, local_c = correlation_dimension(cgr)
            log_c = log_plot(local_c, normalize = False)

            with open(str(I_folder / new_file), "wb") as f:
                pickle.dump(log_c, f)

            stop = timeit.default_timer()

            print(f"Time = {stop - start}")
            # exit()




def doeverything(plotlog = True, plotboth = False, autocorplot = False, dot_size = 2, line_width = 0.25):
    '''
    '''
    gene_folder = cwd / "GenePerChrome"
    cgr: np.ndarray

    cgr_files = walk_directory(gene_folder)
    # cgr_files = [pathlib.Path("D:\Coding\Thesis\FractalDimension\GenePerChrome\FAM177B\EXON\FAM177B_E_3.npy"), 
    #              pathlib.Path("D:\Coding\Thesis\FractalDimension\GenePerChrome\FAM177B\INTRON\FAM177B_I_1.npy")]


    for file in cgr_files:
        # print(file)
        # print(type(file))
        try:
            cgr = np.load(str(file), allow_pickle=True)
        except Exception as e:
            print(type(e))
            print(e)
            print(file)
            cgr = None

        try:
            L, W, = cgr.shape
        except Exception as e:
            # Looks like this is caused by a .png.npy file being saved. Not sure how that happened, but I could probably
            # do some sort of filtering earlier to get around that. Until then this should suffice. Also if something went
            # weird when the file was saved this could cover me.
            print(f"Error after loading\n\t{file}")
            print(type(e))
            print(e)

            cgr = None


        if isinstance(cgr, np.ndarray):
            title = f"{file.name} - {L}x{W} - {int(np.emath.logn(4, L*W))} mer window"

            start = timeit.default_timer()
            correlation_dimension_figure(cgr, title = title, plotlog = plotlog, plotboth = plotboth, autocorplot = autocorplot, dot_size = dot_size, line_width = line_width)
            stop = timeit.default_timer()

            print(f"Time = {stop - start}")
            # exit()



def load_cgr(file):
    '''
    '''
    try:
        cgr = np.load(str(file), allow_pickle = True)
    except Exception as e:
        print(type(e))
        print(e)
        print(file)
        cgr = None

    try:
        L, W, = cgr.shape
    except Exception as e:
        # Looks like this is caused by a .png.npy file being saved. Not sure how that happened, but I could probably
        # do some sort of filtering earlier to get around that. Until then this should suffice. Also if something went
        # weird when the file was saved this could cover me.
        print(f"Error after loading\n\t{file}")
        print(type(e))
        print(e)

        cgr = None

    return cgr



def walk_exon_intron(folder: str or pathlib.Path) -> list:
    '''
    Will only go over the introns and the exons
    '''

    file: str
    iiles = list()
    eiles = list()
    # print(folder)
    # exit()

    for r, _, f in os.walk(str(folder)):
        if pathlib.Path(r).name in "EXON":
            for file in f:
                if file.endswith(".npy"):
                    eiles.append(pathlib.Path(os.path.join(r, file)))
        elif pathlib.Path(r).name in "INTRON":
            for file in f:
                if file.endswith(".npy"):
                    iiles.append(pathlib.Path(os.path.join(r, file)))

    return eiles, iiles


def walk_directory(folder: str or pathlib.Path) -> list:
    '''
    Will walk the direcotry and return everything that ends in the .npy
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
    


    
def correlation_dimension_figure(cgr, title = None, plotlog = True, plotboth = False, autocorplot = False, dot_size = 2, line_width = 0.25):
    '''
    Plot both will plot both the log and the non log plots. This will superseed the logplot.
    '''
    local_c: dict
    local_n: dict

    local_n, local_c = correlation_dimension(cgr)
    log_c = log_plot(local_c, normalize=False)

    if autocorplot:
        auto_x = autocorrelation(log_c)
    # print(type(auto_x))
    # print(auto_x)

    if isinstance(local_n, dict):

        x_values = tuple(local_n.keys())
        y_values = tuple(local_n.values())

        x_logs = tuple(log_c.keys())
        y_logs = tuple(log_c.values())

        if autocorplot:
            x_auto = tuple(auto_x.keys())
            y_auto = tuple(auto_x.values())

        if plotboth:
            # log_ls, m, _ = least_squares(x_logs, y_logs)
            # title = f"{title}<br>slope = {m}"

            fig = go.Figure()
            fig.add_trace(go.Scatter(x = x_values, y = y_values, name = "C(e)", mode = "lines+markers", marker = dict(size = dot_size), line = dict(width = line_width)))
            fig.add_trace(go.Scatter(x = x_logs, y = y_logs, name = "Log", mode = "lines+markers", marker = dict(size = dot_size), line = dict(width = line_width)))
            if autocorplot:
                fig.add_trace(go.Scatter(x = x_auto, y = y_auto, name = "Auto Correlation", mode = "lines+markers", marker = dict(size = dot_size), line = dict(width = line_width)))
            # fig.add_trace(go.Scatter(x = tuple(log_ls.keys()), y = tuple(log_ls.values()), name = "Log 2 LoBF"))
            fig.update_layout(title = title)
            fig.show()

        elif plotlog:
            fig = go.Figure()
            fig.add_trace(go.Scatter(x = x_logs, y = y_logs, name = "Log Transform", mode = "lines+markers", marker = dict(size = dot_size), line = dict(width = line_width)))
            if autocorplot:
                fig.add_trace(go.Scatter(x = x_auto, y = y_auto, name = "Auto Correlation", mode = "lines+markers", marker = dict(size = dot_size), line = dict(width = line_width)))
            # fig.add_trace(go.Scatter(x = tuple(log_ls.keys()), y = tuple(log_ls.values()), name = "Log 2 LoBF"))
            fig.update_layout(title = title)
            fig.show()

        else:
            fig = go.Figure()
            fig.add_trace(go.Scatter(x = x_values, y = y_values, mode = "lines+markers", marker = dict(size = dot_size), line = dict(width = line_width)))
            fig.update_layout(title = title)
            fig.show()
    
    # exit()



def log_plot(C: dict, normalize = True) -> dict:
    '''
    Creates a log-log plot
    '''

    x_values = tuple(C.keys())
    y_values = tuple(C.values())

    x_values = [np.log(x) if x !=0 else 0 for x in x_values]
    y_values = [np.log(y) if y !=0 else 0 for y in y_values]

    if normalize:
        x_max = max(x_values)
        y_max = max(y_values)

        x_values = [x / x_max for x in x_values]
        y_values = [y / y_max for y in y_values]

    return dict(zip(x_values, y_values))

    

    # x_values = tuple(x.keys())
    # y_values = tuple(x.values())

    # x_min, y_min = min(x_values), min(y_values)
    # print(x_min, y_min)

    # x_values = [np.log2(v) if (v != 0) else np.log2(x_min) for v in x_values]
    # y_values = [np.log2(v) if (v != 0) else np.log2(y_min) for v in y_values]
    
    # # x_values = [0 if (v == 0) else np.log2(v) for v in x_values]
    # # y_values = [0 if (v == 0) else np.log2(v) for v in y_values]
    
    # x_values = normalize(x_values)
    # y_values = normalize(y_values)

    # return dict(zip(x_values, y_values))


def normalize(x_list: list):
    '''
    A lazy way of doing this: find the maximum value and divides that value into all other values. The reason for the True/False of max is that
    if the values are negative then it actually needs the minimum value as that has the largest magnitude.
    '''
    index_max = [0 if x == float("-inf") else abs(x) for x in x_list]
    index_max = index_max.index(max(index_max))

    norm = x_list[index_max]
    # print(norm, type(norm))

    return [x / norm for x in x_list]


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



# def correlation_dimension_old(cgr: np.ndarray) -> dict:
#     '''
#     Take the euclidean distance between every point in the matrix. Start with a small radius then increase this to the whole thing.

#     Will need to find the size of the input matrix. Assume it's always a square (in this case it is).
#     Also don't need to calculate every single point over and over again. Can find a way to calculate only half the values. Also once the values
#     are found, save them so I don't need to recalculate them.

#     Also using the Heaviside function, so if (xi - xj) is a point: 1. Else 0.
#     '''
#     non_zero_indices = np.transpose(np.where(cgr > 0))
#     hash_index = lambda x: f"{x[0]}, {x[1]}"

#     C = dict()

#     distances = dict()
#     for i in non_zero_indices:
#         hash_i = hash_index(i)
#         distances[hash_i] = dict()
#         for j in non_zero_indices:
#             hash_j = hash_index(j)
#             if hash_i not in hash_j:
#                 distances[hash_i][hash_j] = 0

#     l, w = cgr.shape
#     # N = (l*w)**(-2)
#     N = non_zero_indices.shape[0]**(-2)
#     size = int(np.sqrt(l**2 + w**2)) + 1

#     for e in range(size):
#         C[e] = 0
        
#         for i in non_zero_indices:
#             hash_i = hash_index(i)

#             for j in non_zero_indices:
#                 hash_j = hash_index(j)

#                 # This part is working
#                 if hash_i not in hash_j:
                
#                     if distances[hash_i][hash_j] == 0:
#                         d = heaviside_function(i, j, e)
#                         distances[hash_i][hash_j] = d
#                         # distances[hash_j][hash_i] = d

#                     C[e] += distances[hash_i][hash_j]
                

#         C[e] = C[e] * N  # note the multiplication: all the work for making N be 1/N**2 was taken care of earlier

#     # exit()

#     return C

#     # for r in range(size):
#     #     for index in non_zero_indices:
#     #         c = heaviside_function()



def correlation_dimension(cgr: np.ndarray) -> tuple:
    '''
    Take the euclidean distance between every point in the matrix. Start with a small radius then increase this to the whole thing.

    Will need to find the size of the input matrix. Assume it's always a square (in this case it is).
    Also don't need to calculate every single point over and over again. Can find a way to calculate only half the values. Also once the values
    are found, save them so I don't need to recalculate them.

    Also using the Heaviside function, so if (xi - xj) is a point: 1. Else 0.
    '''
    C, Cn = list(), list()
    radius, nadius = list(), list()

    L, W = cgr.shape
    size = int(np.sqrt(L**2 + W**2))

    x = np.transpose(np.where(cgr > 0))

    D = distance.DistanceMatrix(x)

    if isinstance(D.D1, np.ndarray):
        N = (D.l)**(-2)

        for e in range(size):
            nadius.append((e + 1) / size), radius.append((e + 1))
            # C[epsilon] = np.count_nonzero((epsilon >= D) & (D > 0))
            Cn.append((2 * np.count_nonzero((e + 1 >= D.D1) & (D.D1 > 0), axis = 0)) * N), C.append((2 * np.count_nonzero((e + 1 >= D.D1) & (D.D1 > 0), axis = 0)))
            # print(C[epsilon])

        C = dict(zip(radius, C))
        Cn = dict(zip(nadius, Cn))

        return Cn, C
    
    else:
        return None, None
    


def derivative(f: dict) -> dict:
    '''
    '''
    forward_d = lambda f, x, h: (f[x + h] - f[x]) / h
    backward_d = lambda f, x, h: (f[x] - (f[x - h])) / h
    central_d = lambda f, x, h: (f[x + h] - (f[x - h])) / (2*h)
    slope: lambda x1, y1, x2, y2: (y1 - y2) / (x1 - x2)

    df = dict()
    norm = 0

    for i, (x, y) in enumerate(f.items()):


        print(i, x, y)
        if i == 0:
            dx = forward_d(f, x, 1)

        elif (i + 1) == len(f):
            dx = backward_d(f, x, 1)

        else:
            dx = central_d(f, x, 1)

        df[x] = dx
        norm += dx

    for x, y in df.items():
        df[x] = y / norm

    return df




def heaviside_function(xyi: float, xyj: float, r: float) -> int:
    '''
    Actually computes the manhattan distance but that's because I'm lazy and don't want to deal with sqrt right now.
    '''

    # d = abs(xyi[0] - xyj[0]) + abs(xyi[1] - xyj[1])
    d = np.sqrt(np.power(xyi[0] - xyj[0], 2) + np.power(xyi[1] - xyj[1], 2))

    if ((r - d) >= 0):
        return 1
    else:
        return 0
    


def autocorrelation(data: dict):
    '''
    '''
    # data = np.array(data)
    # print(data)

    # try:
    #     something = plot_acf(data)
    # except Exception as e:
    #     print(type(e))
    #     print(e)
    X = dict()
    M = tuple(data.values())
    # print(M)

    T = tuple(data.keys())
    # print(T)

    tmax = T[-1]

    for i, t in enumerate(T):
        if i == len(T) - 1:
            break
        alpha = (tmax - t)**(-1)
        # try:
        #     alpha = (tmax - t)**(-1)
        # except ZeroDivisionError as e:
        #     pass
        # except Exception as e:
        #     print(tmax, t)
        #     print(type(e))
        #     print(e)
        #     exit()

        first = 0
        second_a = 0
        second_b = 0

        for tp in range(0, len(T) - i):
            first += M[tp]*M[tp + i]

            second_a += M[tp]
            second_b += M[tp + i]

        X[t] = (alpha*first) - (alpha*second_a*alpha*second_b)

    x_min = abs(min(tuple(X.values())))
    x_max = abs(max(tuple(X.values())))
    norm = x_max if x_min > x_min else x_min

    for t, x in X.items():
        X[t] = x / norm


    return X


if __name__ in "__main__":
    main()
