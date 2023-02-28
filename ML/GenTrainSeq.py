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
# import keras


def import_file(file_path: str or pathlib.Path, seq_name: str, header = None) -> pandas.DataFrame:
    '''
    '''
    if header is not None:
        return_data = pandas.read_csv(str(file_path), header = header)
        print(return_data)
    
    else:
        return_data = pandas.read_csv(str(file_path), header = None)

        print(return_data)
        quit()

    column_names = tuple(return_data.columns)
    new_col_names = dict()

    for name in column_names:
        if re.match("Cord", name):
            return_data.pop(name)
        elif re.match("Seq", name):
            new_col_names[name] = f"{seq_name}_{name}"
        # else:
        #     new_col_names[name] = name

    return_data.pop("GName")
    return_data.pop("EName")

    return_data = return_data.rename(columns=new_col_names)

    return return_data




cwd = pathlib.Path.cwd()
TraningData = cwd.parent / "Data_Files" / "TrainingData"
# data_set_1: pandas.DataFrame = pandas.read_pickle(str(TraningData / "TrainingData_SeqComponent.pkl"))
# print(data_set_1)
# print(data_set_1.shape)
# print(f"\n\n")
# data_set_2: pandas.DataFrame = pandas.read_pickle(str(TraningData / "TrainingGeneData_v2_Sample.pkl"))
# print(data_set_2)
# print(data_set_2.shape)
# print(f"\n\n")
# data_set_3: pandas.DataFrame = pandas.read_pickle(str(TraningData / "TrainingGeneData.pkl"))
# print(data_set_3)
# print(data_set_3.shape)
# print(f"\n\n")
# data_set_4: pandas.DataFrame = pandas.read_pickle(str(TraningData / "TrainingGeneDataSample.pkl"))
# print(data_set_4)
# print(data_set_4.shape)
# print(f"\n\n")
data_set_5: pandas.DataFrame = pandas.read_pickle(str(TraningData / "TrainingGeneData_v2.pkl"))
print(data_set_5)
print(data_set_5.shape)
print(f"\n\n")


# utr5_data = import_file(str(TraningData / "UTR5_seq.csv"), "UTR5", header = 0)
# utr3_data = import_file(str(TraningData / "UTR3_seq.csv"), "UTR3", header = 0)
# cds_data = import_file(str(TraningData / "CDS_seq.csv"), "CDS", header = 0)
# intron_data = import_file(str(TraningData / "Intron_seq.csv"), "INT", header = None)

# sequence_types = {"utr5": utr5_data, "utr3": utr3_data, "cds": cds_data, "intron": intron_data}

# print(data_set_3)
# data_set_3["Chr"] = pandas.Categorical(data_set_3["Chr"])

# chromes = data_set_3["Chr"].unique()

# for chr in chromes:
#     print(chr)

# chrs = pandas.Categorical(data_set_3["Chr"])

# rows, _ = data_set_3.shape[0]

# length = 100_000

# for row in range(rows):
#     row_of_interest = data_set_3.iloc[row, :].copy()

# for stype, data in sequence_types.items():
#     print(data.shape)

# for title in data_set_2:
#     print(title)


