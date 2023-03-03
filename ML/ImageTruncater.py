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


def set_num_seq(in_seq: str, length = 100_000) -> np.array:
    '''
    '''
    sength = len(in_seq)
    image_seq = np.ndarray(shape = (length, 4), dtype = "uint8")

    if sength > length:
        in_seq = in_seq[:length]

    if sength <= length:


        for _ in range(length - sength):
            in_seq = f"{in_seq}X"

        # print(in_seq)

        for x, n in enumerate(in_seq):
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

            # for i in range(4):
            #     print(f"x = {x} i = {i} value = {image_seq[x, i]}")
            #     # if (image_seq[x, i] != 0) or (image_seq[x, i] != 1):
            #     #     print(x, i, image_seq[x, i])

    uvalues = set(np.unique(image_seq))

    if uvalues != {0, 1}:
        print(image_seq)
        print(uvalues)
        raise ValueError("Values should be 0 or 1, got values other then that")
    
    else:
        return image_seq



def set_num_lab(in_tepe: str) -> np.uint8:
    '''
    '''

    label = 1 if in_tepe in "CDS" else 2 if in_tepe in "Intron" else 3 if in_tepe in "UTR5" else 4

    if label < 1:
        raise ValueError("Value should be between 1 and 4 inclusive and value is below 1")
    
    else:
        return label




def digitize_sequences(input_path, out_labels, out_data, length = 100_000):
    # cwd = pathlib.Path.cwd()

    # # length = 100_000

    # TraningData = cwd.parent / "Data_Files" / "TrainingData"
    input_data: pandas.DataFrame = pandas.read_pickle(str(input_path))

    print(input_data)

    rows = input_data.shape[0]

    data, labels = np.ndarray(shape = (rows, length, 4), dtype="uint8"), np.ndarray(shape = (rows, 1), dtype="uint8")

    for row in range(rows):
        row_of_interest = input_data.loc[row, :].copy()

        seq = row_of_interest["Seq"]
        stype = row_of_interest["Type"]

        try:
            label = set_num_lab(stype)
        except ValueError as v:
            print("Inappropriate value for Label")
            print(f"Type = {stype} at row = {row}")
            label = None
        except Exception as e:
            print("Something else went wrong when doing labels")
            print(f"Error type = {type(e)}")
            print(f"Row = {row}, type = {stype}")


        try:
            digital_seq = set_num_seq(seq)
        except ValueError as v:
            print("Inapprorpiate value(s) for digital seqeunce")
            print(f"Type = {stype} at row = {row}")
            digital_seq = None
        except Exception as e:
            print("Something else went wrong when doing Seqs")
            print(f"Error type = {type(e)}")
            print(f"Row = {row}, type = {stype}")


        # print(digital_seq)
        # print(seq)

        if (digital_seq is not None) and (label is not None):
            labels[row] = label
            data[row] = digital_seq


    np.save(str(out_labels), labels), np.save(str(out_data), data)



if __name__ in "__main__":
    digitize_sequences()