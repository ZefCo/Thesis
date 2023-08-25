# https://towardsdatascience.com/convolutional-neural-networks-understanding-and-organizing-your-data-set-ba3e8b4086cb
# https://datagy.io/python-decorators/
# https://gitpython.readthedocs.io/en/stable/intro.html

# https://stackoverflow.com/questions/26134026/how-to-get-the-current-checked-out-git-branch-name-through-pygit2

# https://gitpython.readthedocs.io/en/stable/intro.html

import os
# os.environ["GET_PYTHON_REFRESH"] = "quiet"
# import git
import pandas
import pathlib
cwd = pathlib.Path.cwd()
import shutil
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
from PIL import Image
# import keras
# https://stackoverflow.com/questions/53066762/understanding-1d-convolution-of-dna-sequences-encoded-as-a-one-hot-vector
# https://stackoverflow.com/questions/53514495/what-does-batch-repeat-and-shuffle-do-with-tensorflow-dataset
# https://stackoverflow.com/questions/53066762/understanding-1d-convolution-of-dna-sequences-encoded-as-a-one-hot-vector


def main():
    '''
    This isn't going to be a one size fits all kind of script: the model has to be rebuilt. However, having saved the scripts used and the weights from the scripts,
    this shouldn't be to hard, just a bunch of copy and pasting.

    I need to think about this. I can't get the history without redoing the model (that data is not saved) but could get it via a custom callback. How much is it
    worth it to actually do all that work?
    '''
    version_dir: pathlib.Path = cwd / "FractalModels" / "IvE" / "version_slurm_22"  
    image_dir: pathlib.Path = pathlib.Path("/media/ethanspeakman/Elements/Gene_Data_Sets/FractalImageEvI_DS1_6mer_hist")
    model_file: pathlib.Path = version_dir / "Model"
    # reload_model(version_dir, image_dir)
    reload_saved_model(model_file)


def reload_saved_model(model_dir: pathlib.Path):
    '''
    '''
    imported_model = tf.saved_model.load(str(model_dir))

    print(imported_model.history)

    # print(type(imported_model))


def reload_model(version_dir: pathlib.Path, image_dir: pathlib.Path, output_classes: int = 2, *args, **kwargs):
    '''
    What goes in here has to be copy and pasted from the apropriate model that you are trying to resume.
    '''
    def get_num_pixels(filepath):
        width, height = Image.open(filepath).size
        return width, height

    def latest_checkpoint(model_dir: pathlib.Path):
        latest = tf.train.latest_checkpoint(model_dir)
        return latest

    model_dir = version_dir / "Model"

    w, h = get_num_pixels(image_dir / "EXON" / "Exon_0.png")  # probably shouldn't hard code it to the first image, but I'm feeling lazy this morning
    
    model = create_model(model_dir, version_dir, output_classes, w, h)
    latest = latest_checkpoint(model_dir)

    model.load_weights(latest)


def create_model(model_dir: pathlib.Path, version_dir: pathlib.Path, output_classes: int, w: int, h: int):
    '''
    This should be copy and pasted over from the original model script.
    '''
    input_layer = tf.keras.Input(shape = (w, h, 1))
    a = tf.keras.layers.Conv2D(w, (1, 1), activation = "gelu", kernel_regularizer = tf.keras.regularizers.l2(l = 0.001))(input_layer)
    a = tf.keras.layers.Dropout(.5)(a)
    a = tf.keras.layers.BatchNormalization()(a)
    b = tf.keras.layers.Conv2D(h, (6, 6), activation = "gelu", kernel_regularizer = tf.keras.regularizers.l2(l = 0.001))(a)
    b = tf.keras.layers.Dropout(.5)(b)
    b = tf.keras.layers.BatchNormalization()(b)
    # c = tf.keras.layers.Conv2D(64, (3, 3), activation = "gelu", kernel_regularizer = tf.keras.regularizers.l2(l = 0.001))(b)
    # c = tf.keras.layers.Dropout(.5)(c)
    # c = tf.keras.layers.BatchNormalization()(c)
    # d = tf.keras.layers.Conv2D(64, (1, 1), activation = "gelu", kernel_regularizer = tf.keras.regularizers.l2(l = 0.001))(c)
    # d = tf.keras.layers.Dropout(.5)(d)
    # d = tf.keras.layers.BatchNormalization()(d)
    flatten = tf.keras.layers.Flatten()(b)
    final = tf.keras.layers.Dropout(.5)(flatten)
    output_layer = tf.keras.layers.Dense(output_classes, activation = "softmax")(final)
    model = tf.keras.Model(inputs = input_layer, outputs = output_layer)

    model.compile(optimizer = 'adam',
                loss = 'categorical_crossentropy',
                #   metrics = ['accuracy'])
                metrics = [tf.keras.metrics.CategoricalAccuracy()])
    
    return model

if __name__ in "__main__":
    main()