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
from PIL import Image
import matplotlib.image
# import keras

cwd = pathlib.Path.cwd()

TraningData = cwd.parent / "Data_Files" / "TrainingData"

data_set_1: pandas.DataFrame = pandas.read_pickle(str(TraningData / "TrainingGeneData_v3.pkl"))
print(data_set_1)

data_set_1["Type"] = pandas.Categorical(data_set_1["Type"])
types = data_set_1["Type"].unique()

for t in types:
    dump = data_set_1[data_set_1["Type"] == t]
    print(t, dump.shape)

# data_set_1 = data_set_1.rename(columns = {"Strand": "Chrm", "Chrm": "Strand"})
# data_set_1 = data_set_1.reset_index()

rows, cols = data_set_1.shape

length = 100_000


# # Padding Value
# for row in range(rows):
#     row_of_interest = data_set_1.iloc[row, :].copy()
#     roi_len = len(row_of_interest["Seq"])
#     if length < roi_len:
#         length = roi_len

# print(data_set_1)

data = np.ndarray(shape = (rows, length, 4))
labels = np.ndarray(shape = (rows, 1))


for row in range(rows):
    row_of_interest = data_set_1.iloc[row, :].copy()

    seq = row_of_interest["Seq"]

    for x, n in enumerate(range(20)):
        seq = f"X{seq}"

    if len(seq) > length:
        pass

    else:

        if len(seq) < length:
            for _ in range(length - len(seq)):
                seq = f"{seq}X"
        seq_type = row_of_interest["Type"]
        # if re.search(r"\d", seq_type):
        #     seq_type = re.sub(r"\d", "", seq_type)
        # chr = row_of_interest["Chrm"]
        # strand = "sense" if row_of_interest["Strand"] in "+" else "antisense"

        # print(f"Len = {len(seq)}\tType = {seq_type}\tChrm = {chr}\tStrand = {strand}")
        image_seq = np.ndarray(shape = (len(seq), 4), dtype = np.uint8)

        # image_seq = list()

        for x, n in enumerate(seq):
            if n in "A":
                image_seq[x, 0] = 1
                image_seq[x, 1] = 0
                image_seq[x, 2] = 0
                image_seq[x, 3] = 0
            elif n in "C":
                image_seq[x, 0] = 0
                image_seq[x, 1] = 1
                image_seq[x, 2] = 0
                image_seq[x, 3] = 0
            elif n in "G":
                image_seq[x, 0] = 0
                image_seq[x, 1] = 0
                image_seq[x, 2] = 1
                image_seq[x, 3] = 0
            elif n in "T":
                image_seq[x, 0] = 0
                image_seq[x, 1] = 0
                image_seq[x, 2] = 0
                image_seq[x, 3] = 1
            else:
                for i in range(4):
                    image_seq[x, i] = 0
        
        data[row] = image_seq
        labels[row] = 1 if seq_type in "CDS" else 2 if seq_type in "Intron" else 3 if seq_type in "UTR5" else 4
        # print(image_seq.shape)
        # break

        # data = np.concatenate([data, image_seq])

# print(labels[0])

np.save(cwd / "Data", data)
np.save(cwd / "Labels", labels)


# print(data.shape)
# print(data[0])

    # # print(image_seq.shape)
    # # print(image_seq.dtype)

    # # file_path = TraningData / strand / chr / seq_type
    # file_path = TraningData / "SeqImages" / seq_type
    # file_path.mkdir(parents = True, exist_ok = True)
    # file_name = f"{seq_type}_image_{row}.png"
    # # np.save(str(file_path / file_name), image_seq)

    # try:
    #     matplotlib.image.imsave(str(file_path / file_name), image_seq)
    # except Exception as e:
    #     print(type(e))
    #     print(e)
    #     print(f"Len = {len(seq)}\tChrm = {chr}\t\tStrand = {strand}\tType = {seq_type}")

    # # image = Image.fromarray(image_seq)
    # # image.save(str(file_path / file_name))

    # # if row == 0:
    # #     break

    

