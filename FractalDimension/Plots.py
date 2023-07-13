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
import random
import pickle
import shutil



def main():
    '''
    '''
    # r_files = random_selection()
    # plotfolder(str(cwd / "Single_Gene"))
    eyeball_method()




def random_selection(k = 10) -> list:
    '''
    '''
    file: pathlib.Path

    D_Folder = cwd / "Distance_Dict"
    E_Folder = D_Folder / "EXON"
    I_Folder = D_Folder / "INTRON"
    R_Folder = cwd / "Random_Selection"
    R_Folder.mkdir(parents=True, exist_ok=True)

    e_files = walk_directory(E_Folder, extension = ".pkl")
    i_files = walk_directory(I_Folder, extension = ".pkl")

    e_files = random.choices(e_files, k = k)
    i_files = random.choices(i_files, k = k)

    r_files = e_files + i_files

    for file in r_files:
        new_path = R_Folder / file.name
        shutil.copy(str(file), str(new_path))



def eyeball_method():
    '''
    '''

    m_g = []
    m_e = []
    m_i = []

    fig = go.Figure()

    # slices = {"D:\Coding\Thesis\FractalDimension\Random_Selection\GATAD2A_E_9_280_512x512.pkl": [2, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\EPHX2_I_17_356_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\BPIFB6_I_12_799_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\AMPD3_I_4_2195_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\ABAT_I_12_3806_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\GATAD2A_I_1_49992_512x512.pkl": [1, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\GPR107_I_10_1096_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\GPR107_I_1_21562_512x512.pkl": [1, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\KLC1_I_11_1658_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\KLC1_I_6_583_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\LCLAT1_I_4_5767_512x512.pkl": [1, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\EPHX2_E_11_86_512x512.pkl": [2.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\GPR107_E_15_50_512x512.pkl": [3, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\LRCH1_E_4_106_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\MGAM_E_52_97_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\MGAM_E_63_117_512x512.pkl": [1.5, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\NCAPH2_E_2_102_512x512.pkl": [2.4, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\NCAPH2_E_4_85_512x512.pkl": [1, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\SLC37A1_E_15_49_512x512.pkl": [3, 4],
    #           "D:\Coding\Thesis\FractalDimension\Random_Selection\TMEM236_E_2_73_512x512.pkl": [1, 4]}

    slices = {"D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_10_63_512x512.pkl": [2, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_11_64_512x512.pkl": [1.8, 4.2],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_12_149_512x512.pkl": [1.8, 4.2],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_13_138_512x512.pkl": [1.8, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_14_168_512x512.pkl": [1, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_15_147_512x512.pkl": [2, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_16_112_512x512.pkl": [2.4, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_17_3262_512x512.pkl": [1, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_1_95_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_2_111_512x512.pkl": [1, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_3_98_512x512.pkl": [1.5, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_4_30_512x512.pkl": [3, 4.9],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_5_118_512x512.pkl": [2, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_6_96_512x512.pkl": [1.5, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_7_50_512x512.pkl": [2, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_8_81_512x512.pkl": [0.5, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\EXON\ABAT_E_9_93_512x512.pkl": [2, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_10_1922_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_11_568_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_12_3806_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_13_1972_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_14_1286_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_15_2988_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_16_1718_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_1_60987_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_2_10191_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_3_2009_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_4_2284_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_5_730_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_6_6391_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_7_6262_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_8_588_512x512.pkl": [0, 4],
              "D:\Coding\Thesis\FractalDimension\Single_Gene\INTRON\ABAT_I_9_1377_512x512.pkl": [0, 4]}

    for i, (file, s) in enumerate(slices.items()):

        legend_title = file.split("\\")
        legend_title = legend_title[len(legend_title) - 1]

        x_min, x_max = s[0], s[1]
        data: dict = load_pickle(file)

        x_points, y_points = [], []

        for x, y in data.items():
            if (x >= x_min) and (x <= x_max):
                x_points.append(x)
                y_points.append(y)


        x_points = np.array(x_points)
        y_points = np.array(y_points)

        f, m, b = least_squares(x_points, y_points)
        mx_b = f"{round(m, 6)}*x + {round(b, 6)}"
        m_g.append(m)

        region = re.split("_", legend_title)[1]
        if region in "E":
            m_e.append(m)
        elif region in "I":
            m_i.append(m)

        fig.add_trace(go.Scatter(x = x_points, y = y_points, mode = "markers", name = "Selected Data", legendgroup = i, legendgrouptitle_text = legend_title))
        fig.add_trace(go.Scatter(x = tuple(f.keys()), y = tuple(f.values()), name = mx_b, mode = "lines", legendgroup = i))

        print(f"{legend_title}:\n\tslope = {m}\t\tb = {b}\n")


    m_g, m_e, m_i = np.array(m_g), np.array(m_e), np.array(m_i)
    print(f"GLOBAL: slope mean: {m_g.mean()}\tslope max: {m_g.max()}\tslope min: {m_g.min()}\tslope std: {m_g.std()}")
    print(f"EXON:   slope mean: {m_e.mean()}\tslope max: {m_e.max()}\tslope min: {m_e.min()}\tslope std: {m_e.std()}")
    print(f"INTRON: slope mean: {m_i.mean()}\tslope max: {m_i.max()}\tslope min: {m_i.min()}\tslope std: {m_i.std()}")

    fig.show()

def plotfolder(path):
    '''
    '''

    files = walk_directory(path, extension = ".pkl")

    fig = go.Figure()

    for i, file in enumerate(files):
        print(file)
        d: dict = load_pickle(file)

        x_values = tuple(d.keys())
        y_values = tuple(d.values())

        fig.add_trace(go.Scatter(x = x_values, y = y_values, name = file.name, mode = "lines+markers"))

    fig.show()


def everything():
    '''
    '''
    file: pathlib.Path

    D_Folder = cwd / "Distance_Dict"
    E_Folder = D_Folder / "EXON"
    I_Folder = D_Folder / "INTRON"

    e_files = walk_directory(E_Folder, extension = ".pkl")
    i_files = walk_directory(I_Folder, extension = ".pkl")

    fig = go.Figure()

    for i, file in enumerate(e_files):
        d: dict = load_pickle(file)

        x_values = tuple(d.keys())
        y_values = tuple(d.values())

        if i == 0:
            fig.add_trace(go.Scatter(x = x_values, y = y_values, name = file.name, mode = "lines+markers", legendgroup = "exon", legendgrouptitle_text = "Exons"))
        else:
            fig.add_trace(go.Scatter(x = x_values, y = y_values, name = file.name, mode = "lines+markers", legendgroup = "exon"))

    for i, file in enumerate(i_files):
        d: dict = load_pickle(file)

        x_values = tuple(d.keys())
        y_values = tuple(d.values())

        if i == 0:
            fig.add_trace(go.Scatter(x = x_values, y = y_values, name = file.name, mode = "lines+markers", legendgroup = "intron", legendgrouptitle_text = "Introns"))
        else:
            fig.add_trace(go.Scatter(x = x_values, y = y_values, name = file.name, mode = "lines+markers", legendgroup = "intron"))
    
    fig.show()



def load_pickle(filepath):
    '''
    '''

    with open(str(filepath), "rb") as f:
        pickle_dict = pickle.load(f)

    return pickle_dict


def least_squares(x_points: np.array, y_points: np.array):
    '''
    '''
    if isinstance(x_points, list):
        x_points = np.array(x_points)
    if isinstance(y_points, list):
        y_points = np.array(y_points)

    ave_x = np.mean(x_points)
    ave_y = np.mean(y_points)

    m = np.dot(x_points - ave_x, y_points - ave_y) / np.dot(x_points - ave_x, x_points - ave_x)
    b = ave_y - (m * ave_x)

    ls = {e: m*e + b for e in x_points}

    return ls, m, b




def walk_directory(folder: str or pathlib.Path, extension: str = None) -> list:
    '''
    Will walk the direcotry and return everything that ends in the .npy
    '''

    file: str
    files = list()
    # print(folder)
    # exit()

    if isinstance(extension, str):
        for r, _, f in os.walk(str(folder)):
            for file in f:
                if file.endswith(extension):
                    files.append(pathlib.Path(os.path.join(r, file)))

    else:
        for r, _, f in os.walk(str(folder)):
            for file in f:
                files.append(pathlib.Path(os.path.join(r, file)))


    return files




if __name__ in "__main__":
    main()