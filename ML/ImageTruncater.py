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


def set_num_seq(in_seqa: str, length = 100_000) -> np.array:
    '''
    '''
    in_seq = in_seqa
    in_seq = in_seq.upper()
    sength = len(in_seq)
    image_seq = np.zeros(shape = (length, 4), dtype = "int32")
    # print(image_seq)
    # print(image_seq.shape)
    # exit()

    if sength > length:
        in_seq = in_seq[:length]
        sength = len(in_seq)

    # print(sength, length)

    if sength <= length:

        for _ in range(length - sength):
            in_seq = f"{in_seq}X"

        # print(in_seq)

        for x, n in enumerate(in_seq):
            if n == "A":
                image_seq[x, 0] = 1
            elif n == "C":
                image_seq[x, 1] = 1
            elif n == "G":
                image_seq[x, 2] = 1
            elif n == "T":
                image_seq[x, 3] = 1

    uvalues = set(np.unique(image_seq))

    if uvalues == {0, 1}:    
        return image_seq
    
    elif uvalues == {0}:
        # print(image_seq)
        # print(uvalues)
        # print(in_seqa)
        # print(f"\n\n")
        # print(in_seq)
        # print("\n\n")
        # print(image_seq)
        # print(f"\n\n")
        # # exit()
        raise ValueError("Values should be 0 or 1, this sequences is Null")
    
    else:
        raise ValueError("Values should be 0 or 1, this is a new value error")




def set_num_lab(in_type: str) -> np.int32:
    '''
    '''

    label = 0 if in_type in "CDS" else 1 if in_type in "Intron" else 2 if in_type in "UTR5" else 3

    if (label < 0) or (label > 3):
        raise ValueError("Value should be between 1 and 4 inclusive and value is below 1")
    
    else:
        return label
    

def set_one_hot(in_type) -> np.int32:
    '''
    '''

    label = np.array([1, 0, 0, 0], dtype='int32') if in_type in "CDS" else np.array([0, 1, 0, 0], dtype='int32') if in_type in "Intron" else np.array([0, 0, 1, 0], dtype='int32') if in_type in "UTR5" else np.array([0, 0, 0, 1], dtype='int32')

    if (set(np.unique(label))) == {0, 1}:
        return label
    
    else:
        raise ValueError("Values should be 0 or 1, got values other then that")




def digitize_sequences(input_path, out_labels, out_data, length = 100_000):
    total_seq_error, total_lab_error = 0, 0
    input_data: pandas.DataFrame = pandas.read_pickle(str(input_path))

    out_data = f"{out_data}_L{length}"
    out_labels = f"{out_labels}_L{length}"
    out_hots = f"{out_labels}_OneHots_L{length}"

    print(input_data)

    rows = input_data.shape[0]

    data, labels, one_hots = np.zeros(shape = (rows, length, 4), dtype = "int32"), np.zeros(shape = (rows, 1), dtype = "int32"), np.zeros(shape = (rows, 4), dtype = "int32")
    print(data.shape)

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
            total_lab_error += 1
        except Exception as e:
            print("Something else went wrong when doing labels")
            print(f"Error type = {type(e)}")
            print(f"Row = {row}, type = {stype}")
            label = None
            total_lab_error += 1

        try:
            one_hot = set_one_hot(stype)
        except ValueError as v:
            print("Inappropriate values for One Hot Label")
            print(f"Type = {stype} at row = {row}")
            one_hot = None
            total_lab_error += 1
        except Exception as e:
            print("something else went wrong doing one hot labels")
            print(f"Error type = {type(e)}")
            print(f"Row = {row}, type = {stype}")
            one_hot = None
            total_lab_error += 1

        try:
            digital_seq = set_num_seq(seq, length=length)
        except ValueError as v:
            print("Inapprorpiate value(s) for digital seqeunce")
            print(f"Type = {stype} at row = {row}")
            digital_seq = None
            total_seq_error += 1
        except Exception as e:
            print("Something else went wrong when doing Seqs")
            print(f"Error type = {type(e)}")
            print(f"Row = {row}, type = {stype}")
            digital_seq = None
            total_seq_error += 1
        
        # exit()


        # print(digital_seq)
        # print(seq)

        if (digital_seq is not None) and (label is not None) and (one_hot is not None):
            # if label == 0:
            #     print(f"### {row} ###")
            labels[row] = label
            data[row] = digital_seq
            one_hots[row] = one_hot
        
        else:
            labels[row] = -1
            data[row] = -1
            one_hots[row] = -1
            

    # print(f"\n\n")
    # for row in range(rows):
    #     if labels[row] == 0:
    #         print(f"Row = {row}")
    
    print(np.unique(labels))
    print(np.unique(data))
    print(np.unique(one_hots))
    print(total_seq_error, total_lab_error)

    # print(labels)
    # print(one_hots)


    np.save(str(out_labels), labels), np.save(str(out_data), data), np.save(str(out_hots), one_hots)



if __name__ in "__main__":
    cwd = pathlib.Path.cwd()
    # length = 100_000
    TraningData = cwd.parent / "Data_Files" / "TrainingData"

    digitize_sequences(str(TraningData / "TrainingGeneData_v3.pkl"), str("Labels"), str("Data"), length=1000)