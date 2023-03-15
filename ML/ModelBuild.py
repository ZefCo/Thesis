import pandas
import pathlib
import re
# from tensorflow.keras.models import Sequential, load_module
# from tensorflow.keras.layers import Conv1D
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import accuracy_score
import tensorflow as tf
# print(tf.__version__)
# exit()
import random
import numpy as np
import pandas as pd
# import plotly.graph_objects as go
import datetime
from matplotlib import pyplot as plt
from contextlib import redirect_stdout
# import keras
# https://stackoverflow.com/questions/53066762/understanding-1d-convolution-of-dna-sequences-encoded-as-a-one-hot-vector
# https://stackoverflow.com/questions/53514495/what-does-batch-repeat-and-shuffle-do-with-tensorflow-dataset
# https://stackoverflow.com/questions/53066762/understanding-1d-convolution-of-dna-sequences-encoded-as-a-one-hot-vector


def make_dir(dir_path: pathlib.Path):
    '''
    '''
    if not dir_path.is_dir():
        dir_path.mkdir()

    return dir_path





def model_summery(summery, now):
    '''
    '''
    with open(f"Model_Loss_{now.year}-{now.month}-{now.day}_{now.hour}-{now.minute}-{now.second}.txt", "a") as f:
        print(summery, file = f)



def import_data(data_filepath, label_filepath, onehot_filepath):
    '''
    '''
    train_lb = np.load(label_filepath)
    keep = np.where(train_lb >= 0)[0]
    train_lb = train_lb[keep]

    train_on = np.load(onehot_filepath)
    train_on = train_on[keep]

    train_ds = np.load(data_filepath)
    train_ds = train_ds[keep]

    print(f"Imported DS Shape = {train_ds.shape}")
    print(f"Imported Labels Shape = {train_lb.shape}")
    print(f"Imported One Hot Shape = {train_on.shape}")

    print(f"Unique Train Labels = {np.unique(train_lb)}")
    print(f"That should be {0, 1, 2, 3}")

    # all_samples = train_lb.shape[0]

    print(f"Keep dimension = {keep.shape}")
    print(f"Hopefully that is the same size as the data")

    return train_ds, train_lb, train_on


def data_groups(data, labels, onehot):
    '''
    '''
    all_samples = labels.shape[0]

    random_sample = int(0.2 * all_samples)
    tv_samples = np.random.choice(all_samples, size = random_sample, replace=False)
    tv_sample = tv_samples.shape[0]

    test_sample = int(tv_sample / 2)
    valid_sample = tv_sample - test_sample

    test_samples = tv_samples[0:test_sample]
    valid_samples = tv_samples[valid_sample - 1:tv_samples.shape[0]]

    test_ds = data[test_samples]
    print(f"Testing Data set = {test_ds.shape}")
    test_lb = labels[test_samples]
    print(f"Testing Labels = {test_lb.shape}")
    test_on = onehot[test_samples]
    print(f"Testing One Hots = {test_on.shape}")

    valid_ds = data[valid_samples]
    print(f"Valid Data set = {valid_ds.shape}")
    valid_lb = labels[valid_samples]
    print(f"Valid Labels = {valid_lb.shape}")
    valid_on = onehot[valid_samples]
    print(f"Valid One Hots = {valid_on.shape}")

    train_ds = np.delete(data, tv_samples, axis = 0)
    print(f"Training Data set = {train_ds.shape}")
    train_lb = np.delete(labels, tv_samples, axis = 0)
    print(f"Training Labels = {train_lb.shape}")
    train_on = np.delete(onehot, tv_samples, axis = 0)
    print(f"Training One Hots = {train_on.shape}")

    print("Datasets, labels, and one hots should all match in size")

    return train_ds, train_lb, train_on, test_ds, test_lb, test_on, valid_ds, valid_lb, valid_on



def main(data_filepath, label_filepath, onehot_filepath,
         epochs = 100, batch_size = 250):
    '''
    '''

    model_summery_dir: pathlib.Path = make_dir(cwd / "Model_Summeries")
    model_history_dir: pathlib.Path = make_dir(cwd / "Model_History")
    model_lp_dir: pathlib.Path = make_dir(cwd / "Model_Loss_Plots")
    model_ap_dir: pathlib.Path = make_dir(cwd / "Model_Accuracy_Plots")

    train_ds, train_lb, train_on = import_data(data_filepath, label_filepath, onehot_filepath)
    train_ds, train_lb, train_on, test_ds, test_lb, test_on, valid_ds, vlaid_lb, valid_on = data_groups(train_ds, train_lb, train_on)

    nuc_length = train_ds.shape[1]


    data_size = train_ds.shape[0]

    steps_per_epoch = int(data_size / batch_size) + 1
    print(f"Model will use {batch_size} batches per epoch")
    print(f"Steps per epoch = {steps_per_epoch}\nMeaning {batch_size * steps_per_epoch} Seqs will be processed per Epoch\nFor {data_size} Seqs")

    # https://stackoverflow.com/questions/64787511/tensorflow-dataset-batch-size-and-steps-per-epoch

    dataset = tf.data.Dataset.from_tensor_slices((train_ds, train_on))
    dataset = dataset.shuffle(buffer_size = data_size).batch(batch_size)
    dataset = dataset.repeat()

    # dataset = dataset.batch(batch_size).repeat()

    testset = tf.data.Dataset.from_tensor_slices((test_ds, test_on))
    validset = tf.data.Dataset.from_tensor_slices((valid_ds, valid_on))
    validset = validset.batch(batch_size)
    # validset = validset.repeat()

    input_layer = tf.keras.Input(shape = (nuc_length, 4))
    x = tf.keras.layers.Conv1D(filters = 3, kernel_size = 3)(input_layer)
    # x = tf.keras.layers.Dropout(.25)(x)
    # x = tf.keras.layers.Conv1D(filters = 2, kernel_size = 1, activation="tanh")(x)
    # x = tf.keras.layers.Dropout(.25)(x)
    # x = tf.keras.layers.Conv1D(filters = 3, kernel_size = 1, activation="tanh")(x)
    # x = tf.keras.layers.Dropout(.25)(x)
    x = tf.keras.layers.Conv1D(filters = 20, kernel_size = 3, activation="relu")(x)
    # x = tf.keras.layers.Dropout(.4)(x)
    x = tf.keras.layers.Conv1D(filters = 50, kernel_size = 3, activation="relu")(x)
    # x = tf.keras.layers.Dropout(.4)(x)
    x = tf.keras.layers.Conv1D(filters = 200, kernel_size = 3, activation="relu")(x)
    x = tf.keras.layers.Conv1D(filters = 500, kernel_size = 3, activation="relu")(x)
    x = tf.keras.layers.Flatten()(x)
    output_layer = tf.keras.layers.Dense(4, activation = "softmax")(x)
    # output_layer = tf.keras.layers.Dense(3, activation = "softmax")(x)

    # # 2/28/2023 Attempt

    # input_layer = tf.keras.Input(shape = (100000, 4))
    # x = tf.keras.layers.Conv1D(4, 1)(input_layer)
    # x = tf.keras.layers.Conv1D(4, 1, strides = 20)(x)
    # # x = tf.keras.layers.Conv1D(4, 1)(x)
    # # x = tf.keras.layers.Conv1D(4, 1)(x)
    # x = tf.keras.layers.Flatten()(x)
    # output_layer = tf.keras.layers.Dense(4, activation = "softmax")(x)
    # # output_layer = tf.keras.layers.Dense(3, activation = "softmax")(x)

    now = datetime.datetime.now()

    model = tf.keras.Model(inputs = input_layer, outputs = output_layer)

    model.compile(optimizer = 'adam',
                loss = 'categorical_crossentropy',
                #   metrics = ['accuracy'])
                metrics = [tf.keras.metrics.CategoricalAccuracy()])

    
    with open(str(model_summery_dir / f"Model_Summary_{now.year}-{now.month}-{now.day}_{now.hour}-{now.minute}-{now.second}.txt"), "w") as f:
        with redirect_stdout(f):
            model.summary()

    # model.fit(train_ds, train_lb, batch_size=50, epochs=3, steps_per_epoch=10, validation_split=.2)
    history = model.fit(dataset, 
                        batch_size = batch_size, 
                        epochs = epochs, 
                        steps_per_epoch = steps_per_epoch,
                        validation_data = validset)

    # def eval_input_fn():
    #     batched_dataset = dataset.test(flags_obj.data_dir).batch(flags_obj.batch_size)
    #     return batched_dataset.__iter__()

    # File "c:\Users\tokyo\anaconda3\envs\thesis\lib\site-packages\keras\backend.py", line 5532, in categorical_crossentropy
            # target.shape.assert_is_compatible_with(output.shape)
        # ValueError: Shapes (None, 1) and (None, 2) are incompatible

    # Is the shape of the cross entropy wrong?

    results = model.evaluate(test_ds, test_on, batch_size=batch_size)
    print(f"Results are: loss = {results[0]}, accuracy = {results[1]}")
    # print(history)
    history_data = history.history

    pandas.DataFrame(history_data).to_csv(str(model_history_dir / f"Model_History_{now.year}-{now.month}-{now.day}_{now.hour}-{now.minute}-{now.second}.csv"))

    # print(history_data.keys())

    plt.plot(history_data["loss"])
    plt.plot(history_data["val_loss"])
    plt.xlabel("Epochs")
    plt.legend(["Loss", "Valid Loss"])
    plt.savefig(str(model_lp_dir / f"Model_Loss_{now.year}-{now.month}-{now.day}_{now.hour}-{now.minute}-{now.second}.png"))
    plt.close()

    plt.plot(history_data["categorical_accuracy"])
    plt.plot(history_data["val_categorical_accuracy"])
    plt.xlabel("Epochs")
    plt.legend(["Accuracy", "Valid Accuracy"])
    plt.savefig(str(model_ap_dir / f"Model_Accuracy_{now.year}-{now.month}-{now.day}_{now.hour}-{now.minute}-{now.second}.png"))
    plt.close()


if __name__ in "__main__":
    cwd = pathlib.Path.cwd()
    # print(datetime.datetime.now().year)
    # now = datetime.datetime.now()
    # print(f"Model_Preformance_{now.year}-{now.month}-{now.day}_{now.hour}-{now.minute}-{now.second}.html")
    # print(now.hour)
    # print(now.minute)
    # print(now.day)

    main(str(cwd / "Data_L1000.npy"), str(cwd / "Labels_L1000.npy"), str(cwd / "Labels_OneHots_L1000.npy"),
        epochs = 100)