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

# Split the sequence
# Classify the splits


# https://stackoverflow.com/questions/41711190/keras-how-to-get-the-output-of-each-layer

def sample_data():
    '''
    '''

    cwd = pathlib.Path.cwd()

    train_data: pandas.DataFrame = pandas.read_pickle(cwd.parent / "Data_Files" / "TrainingData" / "TrainingGeneDataSample.pkl")
    train_data.pop("ename"), train_data.pop("gname")
    train_data = train_data.reset_index()

    train_data_cols = tuple(train_data.columns)

    # print(train_data.shape)

    # utr3_list, cds_list, intron_list, utr5_list = list(), list(), list(), list()

    for col in train_data_cols:
        if re.search("Cord", col):
        # if re.search("Seq", col):
            train_data.pop(col)

    train_data_2 = pandas.DataFrame(columns=["Strand", "Chrm", "Seq", "Type"])

    rows, cols = train_data.shape

    for row in range(rows):
        row_of_interest = train_data.iloc[row, :]

        row_of_interest = row_of_interest.dropna()
        strand, chrm = row_of_interest["Chr"], row_of_interest["Strand"]

        indicies = list(row_of_interest.index)
        indicies.remove("index"), indicies.remove("Name"), indicies.remove("ncibname"), indicies.remove("Chr"), indicies.remove("Strand")
        ran_index = random.randint(0, len(indicies) - 1)
        # print(len(indicies), ran_index)
        # print(indicies)
        try:
            index = indicies[ran_index]
        except Exception as e:
            print(e)
            print(type(e))
            print(indicies)
            print(len(indicies), ran_index)

            break

        seq = row_of_interest[index]

        seq_type = re.sub(r"Seq\d", "", index)

        # print(seq_type, type(seq_type))

        new_row = pandas.Series(data = [strand, chrm, seq, seq_type], index = ["Strand", "Chrm", "Seq", "Type"])
        train_data_2 = pandas.concat([train_data_2, new_row.to_frame().T])

    # print(train_data_2)
    train_data_2.to_pickle(cwd.parent / "Data_Files" / "TrainingData" / "TrainingData_SeqComponent.pkl")


def convert_to_numpy(df: pandas.DataFrame) -> np.ndarray:
    '''
    Does more then just convert to a numpy ndarray: also adds in a dummy dimension for use with a 1D Convolution.
    '''

    if isinstance(df, pandas.DataFrame):
        # print("Initial DF")
        # print(df)
        df = df.to_numpy()
        # print(f"\n\nConverted to Numpy")
        # print(df)
        df = df[..., None]
        # print(f"\n\nPadded with Dummy Dimension")
        # print(df)

        return df

    else:
        return None


def old_main():
    '''
    Keeping this here because I might need some of it later.
    '''

    cwd = pathlib.Path.cwd()
    train_data: pandas.DataFrame = pandas.read_pickle(cwd.parent / "Data_Files" / "TrainingData" / "TrainingData_SeqComponent.pkl")
    train_data = train_data.rename(columns = {"Strand": "Chrm", "Chrm": "Strand"})
    train_data = train_data.reset_index()
    # print(train_data)

    # train_data = train_data.to_numpy()
    # print(train_data)

    # https://stackoverflow.com/questions/64631701/valueerror-expected-min-ndim-3-found-ndim-2

    # train_data = train_data[..., "None"]

    # X = pandas.get_dummies(train_data.drop(["Strand", "Chrm"], axis = 1))
    X = train_data[["Seq", "Chrm", "Strand"]]
    # print(X)
    # print(f"\n\n")
    X = convert_to_numpy(X)

    y = train_data["Type"].apply(lambda x: 1 if x == "CDS" else 2 if x == "Intron" else 3 if x == "UTR5" else 4)
    # print(y)
    # print(f"\n\n")
    y = convert_to_numpy(y)

    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = .2)
    # print(type(X_train))

    # # print(X_train.head())

    model = tf.keras.models.Sequential()  # init model

    model.add(tf.keras.layers.Conv1D(64, 20))  # add layers

    model.add(tf.keras.layers.Dense(units = 1))  # final output layer
    model.compile(loss = "binaryCrossentropy", optimizer = "sgd", metrics = "accuracy")

    try:
        model.fit(X, y, epochs = 200, batch_size = 32)
    except Exception as e:
        print("Error encountered")
        print(type(e))
        print(e)


def main():
    '''
    '''
    cwd = pathlib.Path.cwd()
    # train_dir = cwd.parent / "Data_Files" / "TrainingData" / "SeqImages"

    batch_size = 10
    seed = 123

    train_ds = np.load(cwd / "Data.npy")
    train_lb = np.load(cwd / "Labels.npy")
    # print(train_ds.shape)

    valid_ds = train_ds[300:]
    valid_lb = train_lb[300:]

    index2delete = [i for i in range(300, 349)]

    train_ds = np.delete(train_ds, index2delete, 0)
    train_lb = np.delete(train_lb, index2delete, 0)

    input_shape = (None, train_ds.shape[1], train_ds.shape[2], train_ds.shape[3])
    # print(input_shape)

    print(train_ds.shape, train_lb.shape)
    print(valid_ds.shape, valid_lb.shape)

    # exit()

    # valid_ds = train_ds[]

    # try:
    #     train_ds = tf.keras.utils.image_dataset_from_directory(train_dir, validation_split=.2, batch_size=batch_size, subset="training", seed = seed)
    # except Exception as e:
    #     print("Error in training dataset")
    #     print(type(e))
    #     print(e)

    # try:
    #     valid_ds = tf.keras.utils.image_dataset_from_directory(train_dir, validation_split=.2, batch_size=batch_size, subset="validation", seed = seed)
    # except Exception as e:
    #     print("Error in valid dataset")
    #     print(type(e))
    #     print(e)

    # sanity = train_ds.unbatch()
    # print(sanity)

    # print(type(train_ds))

    model = tf.keras.models.Sequential()  # init model

    model.add(tf.keras.layers.Conv1D(1, 1, activation = "relu"))  # add layers
    model.build(input_shape)

    model.summary()

    # # model.add(tf.keras.layers.Conv1D(4, 4))  # add layers
    model.compile(loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=False), optimizer = "adam", metrics = ["accuracy"])


    try:
        model.fit(train_ds, train_lb, epochs = 10, validation_data = (valid_ds, valid_lb))
    except Exception as e:
        print("Error in model fit")
        print(f"\n### {type(e)} ###\n\n")
        print(e)
        print(f"\n\n")

    print("### Finished ###")




if __name__ in "__main__":
    main()