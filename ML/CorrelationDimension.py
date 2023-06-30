import pathlib
cwd = pathlib.Path.cwd()
import os
import re
import numpy as np
import glob
from plotly import graph_objects as go
import timeit
from scipy.spatial import distance

logplot = False

def main():
    '''
    '''
    gene_folder = cwd / "GenePerChrome"

    cgr_files = []

    for root, _, files in os.walk(gene_folder):
        for file in files:
            if file.endswith(".npy"):
                cgr_files.append(os.path.join(root, file))

    for file in cgr_files:
        file = pathlib.Path(file)
        # print(file)
        # print(type(file))
        try:
            cgr = np.load(str(file), allow_pickle=True)
        except Exception as e:
            print(type(e))
            print(e)
            print(file)
            cgr = None

        if isinstance(cgr, np.ndarray):
            title = file.name

            start = timeit.default_timer()
            correlation_dimension_figure(cgr, title = title, logplot = logplot)
            stop = timeit.default_timer()

            print(f"Time = {stop - start}")
            # exit()


    
def correlation_dimension_figure(cgr, title = None, logplot = True):
    '''
    '''

    local_c: dict = correlation_dimension(cgr)

    if isinstance(local_c, dict):

        if logplot:
            x_values = np.log(list(local_c.keys()), where = list(local_c.keys()) !=0)
            y_values = np.log(list(local_c.values()), where = list(local_c.values()) !=0)

        else:
            x_values = list(local_c.keys())
            y_values = list(local_c.values())

        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x = x_values, y = y_values))
        fig.update_layout(title = title)
        fig.show()




def correlation_dimension_old(cgr: np.ndarray) -> dict:
    '''
    Take the euclidean distance between every point in the matrix. Start with a small radius then increase this to the whole thing.

    Will need to find the size of the input matrix. Assume it's always a square (in this case it is).
    Also don't need to calculate every single point over and over again. Can find a way to calculate only half the values. Also once the values
    are found, save them so I don't need to recalculate them.

    Also using the Heaviside function, so if (xi - xj) is a point: 1. Else 0.
    '''
    non_zero_indices = np.transpose(np.where(cgr > 0))
    hash_index = lambda x: f"{x[0]}, {x[1]}"

    C = dict()

    distances = dict()
    for i in non_zero_indices:
        hash_i = hash_index(i)
        distances[hash_i] = dict()
        for j in non_zero_indices:
            hash_j = hash_index(j)
            if hash_i not in hash_j:
                distances[hash_i][hash_j] = 0

    l, w = cgr.shape
    # N = (l*w)**(-2)
    N = non_zero_indices.shape[0]**(-2)
    size = int(np.sqrt(l**2 + w**2)) + 1

    for e in range(size):
        C[e] = 0
        
        for i in non_zero_indices:
            hash_i = hash_index(i)

            for j in non_zero_indices:
                hash_j = hash_index(j)

                # This part is working
                if hash_i not in hash_j:
                
                    if distances[hash_i][hash_j] == 0:
                        d = heaviside_function(i, j, e)
                        distances[hash_i][hash_j] = d
                        # distances[hash_j][hash_i] = d

                    C[e] += distances[hash_i][hash_j]
                

        C[e] = C[e] * N  # note the multiplication: all the work for making N be 1/N**2 was taken care of earlier

    # exit()

    return C

    # for r in range(size):
    #     for index in non_zero_indices:
    #         c = heaviside_function()



def correlation_dimension(cgr: np.ndarray) -> dict:
    '''
    Take the euclidean distance between every point in the matrix. Start with a small radius then increase this to the whole thing.

    Will need to find the size of the input matrix. Assume it's always a square (in this case it is).
    Also don't need to calculate every single point over and over again. Can find a way to calculate only half the values. Also once the values
    are found, save them so I don't need to recalculate them.

    Also using the Heaviside function, so if (xi - xj) is a point: 1. Else 0.
    '''
    C = dict()

    L, W = cgr.shape
    size = int(np.sqrt(L**2 + W**2)) + 1

    D = euc_dis_matrix(cgr)
    # print(type(D))
    # exit()
    # D = distance.c

    if isinstance(D, np.ndarray):
        N = (D.size)**(-2)

        for e in range(size):
            epsilon = e + 1
            # C[epsilon] = np.count_nonzero((epsilon >= D) & (D > 0))
            C[epsilon] = np.count_nonzero(epsilon >= D)
            
        #     for i in non_zero_indices:
        #         hash_i = hash_index(i)

        #         for j in non_zero_indices:
        #             hash_j = hash_index(j)

        #             # This part is working
        #             if hash_i not in hash_j:
                    
        #                 if distances[hash_i][hash_j] == 0:
        #                     d = heaviside_function(i, j, e)
        #                     distances[hash_i][hash_j] = d
        #                     # distances[hash_j][hash_i] = d

        #                 C[e] += distances[hash_i][hash_j]
                    

            C[epsilon] = C[epsilon] * N  # note the multiplication: all the work for making N be 1/N**2 was taken care of earlier

        # # exit()

        return C
    
    else:
        return None




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
    


def euc_dis_matrix(matrix) -> np.array:
    '''
    '''
    euc = lambda x1, x2: np.sqrt(((x1[0] - x2[0])**2) + ((x1[1] - x2[1])**2))
    points = lambda x: (0.5*x**2) - (0.5*x)

    x = np.transpose(np.where(matrix > 0))

    L, W = matrix.shape
    # N = (l*w)**(-2)
    l, _ = x.shape

    # slower and has same problem
    # D = distance.cdist(x, x, "euclidean")
    # return D

    # fast, but has an issue with the size of the matrix. Takes up to much memory
    # try:
    #     D = np.zeros(shape=(l, l))
    # except np.core._exceptions._ArrayMemoryError as e:
    #     print("memory error... stupid non root privledge...")
    #     print(f"Size = {l*l}\tItem Size = 8?\tGB = {l*l*8*1E-9}")

    #     return None

    # except Exception as e:
    #     print(type(e))
    #     print(e)

    #     print(f"Size = {l*l}\tItem Size = 8?\tGB = {l*l*8*1E-9}")
    #     return None
    # else:

    #     print(f"Size = {D.size}\tItem Size = {D.itemsize}\tGB = {D.size*D.itemsize*1E-9}")

    #     w = len(x)
    #     for k in range(w):
    #         step = 1 + k
    #         for i in range(0, l - k - 1):
    #             d = euc(x[i], x[i + step])

    #             D[i, i + step] = d
    #             D[i + step, i] = d

    #     return D

    p = points(l)
    try:
        D = np.zeros(int(p))

    except np.core._exceptions._ArrayMemoryError as e:
        print("Still have a memory error")
        print(f"Size = {p}\tItem Size = 8?\tGB = {p*8*1E-9}")
        return None

    except Exception as e:
        print(type(e))
        print(e)

        print(f"Size = {p}\tItem Size = 8?\tGB = {p*8*1E-9}")
        return None

    else:
        print(f"Size = {D.size}\tItem Size = {D.itemsize}\tGB = {D.size*D.itemsize*1E-9}")

        w = len(x)
        for k in range(w):
            step = 1 + k
            for i in range(0, l - k - 1):
                d = euc(x[i], x[i + step])

                D[i] = d

        return D






if __name__ in "__main__":
    main()
