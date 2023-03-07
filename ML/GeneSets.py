import pandas
import pathlib
import re
# from tensorflow.keras.models import Sequential, load_module
# from tensorflow.keras.layers import Conv1D
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import tensorflow as tf
import random
import numpy as np
import pandas as pd
# import plotly.graph_objects as go
from matplotlib import pyplot as plt



def fileter_frame(indata: pandas.DataFrame, condition):
    '''
    '''

    data: pandas.DataFrame = indata[condition]
    data = data.reset_index()
    data.pop("index")

    return data



def sets(indata: pandas.Series, data_len: int, window_size: int):
    '''
    '''
    if not isinstance(data_len, int):
        data_len = int(data_len)
    
    data_sets = []
    
    for window in range(data_len - window_size + 1):
        data_sets.append(indata[window: window + window_size])
    
    data_sets = set(data_sets)
    
    return data_sets



def main(seq_path: str, window_size: int = 10):
    '''
    '''
    cds_data: pandas.DataFrame = pandas.read_pickle(seq_path)

    cds_data = fileter_frame(cds_data, cds_data["Type"] == "CDS")
    # # print(cds_data.shape)
    rows, _ = cds_data.shape

    cds_data["Len"] = 0
    cds_data["Sets"] = None
    cds_data["Intersections"] = None

    for row in range(rows):
        row_of_interest = cds_data.iloc[row, :].copy()
        length_of_interest = len(row_of_interest["Seq"])
        # cds_data.loc[row, "Len"] = len(row_of_interest["Seq"])

        if length_of_interest >= window_size:
            cds_data.loc[row, "Len"] = length_of_interest
            
            try:
                cds_sets = sets(row_of_interest["Seq"], length_of_interest, window_size)
                cds_data.at[row, "Sets"] = cds_sets
            
            except ValueError as v:
                print(row)
                print(type(cds_sets))
                print(cds_sets)
                print(row_of_interest["Seq"])
                print(length_of_interest)
                print(window_size)

                print(f"\n\n")
                print(type(v))
                print(v)
                exit()

            except Exception as e:
                print(type(e))
                print(e)

                exit()

    cds_data = fileter_frame(cds_data, cds_data["Len"] >= window_size)

    # print(cds_data)

    rows, _ = cds_data.shape

    for i in range(rows):
        i_set: set = cds_data.loc[i, "Sets"]

        i_intersection = set()

        for j in range(rows):
            if i != j:
                j_set: set = cds_data.loc[j, "Sets"]

                intersection = i_set.intersection(j_set)

                i_intersection = i_intersection | intersection

        cds_data.at[i, "Intersections"] = i_intersection

    print(cds_data)



if __name__ in "__main__":
    cwd = pathlib.Path.cwd()
    main(str(cwd / "TrainingGeneData_v3.pkl"))