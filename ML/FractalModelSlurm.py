print("Starting script")
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
import shutil
# import re
# from tensorflow.keras.models import Sequential, load_module
# from tensorflow.keras.layers import Conv1D
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import accuracy_score
import tensorflow as tf
# print(tf.__version__)
# exit()
import random
# import numpy as np
# import pandas as pd
# import plotly.graph_objects as go
# import datetime
# from matplotlib import pyplot as plt
from contextlib import redirect_stdout
# from PIL import Image
# import sys
# import keras
# https://stackoverflow.com/questions/53066762/understanding-1d-convolution-of-dna-sequences-encoded-as-a-one-hot-vector
# https://stackoverflow.com/questions/53514495/what-does-batch-repeat-and-shuffle-do-with-tensorflow-dataset
# https://stackoverflow.com/questions/53066762/understanding-1d-convolution-of-dna-sequences-encoded-as-a-one-hot-vector
cwd = pathlib.Path.cwd()
import datetime
# exit()
width, height, depth = 64, 64, 1
# def get_num_pixels(filepath):
#     '''
#     So this used to be the way that the height and width were determined: it just needed a filepath and would open the file, determine the pixel size
#     and use that for the width and height input of tensorflow. It was a great system that allowed me to ignore the image size and just run the code. But
#     as of 3/21/24, noon, when I load the Anaconda and TensorFlow module they somehow conflict and I loose access to matplotlib and PIL and tensorflow.
#     The work-around: this is now hard coded. You HAVE to manually tell it what the image size it. Until UH fixes itself (HA!) you're fucked.
#     '''
#     width, height = Image.open(filepath).size
#     return width, height


# report = git.Repo(search_parent_directories=True)
# branch = report.active_branch
branch = "slurm"

image_dir = cwd / "GMSubsample"
ive_dir = cwd / "FractalModels" / "IvE"
date = datetime.datetime.now()
version_num = f"{date.year}-{date.month}-{date.day}-{date.hour}-{date.minute}-{date.second}-{date.microsecond}"# len(next(os.walk(ive_dir))[1]) + 60
version_dir = ive_dir / f"version_{branch}_{version_num}"
version_dir.mkdir(parents = True, exist_ok = True)
model_dir = version_dir / "Model"
model_dir.mkdir(parents = True, exist_ok = True)
output_classes = 2

# image_dir = cwd / "FractalImageCvIvU"
# model_dir = cwd / "FractalModels" / "CvIvU"
# version_num = len(next(os.walk(model_dir))[1])
# version_dir = model_dir / f"version_{version_num}"
# output_classes = 4


# w, h = get_num_pixels(image_dir / "EXON" / "Exon_0.png")  # probably shouldn't hard code it to the first image, but I'm feeling lazy this morning
# The reason you don't hard code this is because there might not be an exon_0. Since there might not be you need to grab a random image and check that size.
# w, h = get_num_pixels(image_dir / "EXON" / random.choice(os.listdir(image_dir / "EXON")))

seed = random.randint(1000000, 9000000)


batch_size = 521
epochs = 100

pngs = 0
for folder in os.listdir(image_dir):
    for file in os.listdir(image_dir / folder):
        if file.endswith("png"):
            pngs += 1


def normalization(image, label):
    image = tf.cast(image / 255, tf.float32)
    return image, label

train_set = tf.keras.preprocessing.image_dataset_from_directory(str(image_dir), image_size = (w, h), shuffle = True, seed = seed, subset = "training", validation_split = 0.2, batch_size = batch_size, label_mode = "categorical", color_mode = "grayscale")
valid_set = tf.keras.preprocessing.image_dataset_from_directory(str(image_dir), image_size = (w, h), shuffle = True, seed = seed, subset = "validation", validation_split = 0.2, batch_size = batch_size, label_mode = "categorical", color_mode = "grayscale")

train_set = train_set.map(normalization)
valid_set = valid_set.map(normalization)

steps_per_epoch = int(pngs / batch_size) + 1
# print(steps_per_epoch)

cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath = str(model_dir))
python_script = pathlib.Path(__file__)
shutil.copy(str(python_script), str(version_dir / python_script.name))

input_layer = tf.keras.Input(shape = (width, height, depth))
a = tf.keras.layers.Conv2D(3, (3, 3), padding = "same", activation="gelu")(input_layer) #), activation = "gelu", kernel_regularizer = tf.keras.regularizers.l2(l = 0.001))(input_layer)
# a = tf.keras.layers.ReLU()(a)
a = tf.keras.layers.BatchNormalization()(a)
a = tf.keras.layers.MaxPooling2D(pool_size = (3, 3))(a)
b = tf.keras.layers.Conv2D(5, (3, 3), padding = "same", activation="gelu")(a)
# b = tf.keras.layers.ReLU()(b)
b = tf.keras.layers.BatchNormalization()(b)
b = tf.keras.layers.MaxPooling2D(pool_size = (2, 2))(b)
c = tf.keras.layers.Conv2D(9, (3, 3), padding = "same", activation="gelu")(b)
# c = tf.keras.layers.ReLU()(c)
c = tf.keras.layers.BatchNormalization()(c)
flatten = tf.keras.layers.Flatten()(c)
dense = tf.keras.layers.Dense(50, activation = "gelu", kernel_regularizer = tf.keras.regularizers.l2(l = 0.01))(flatten)
final = tf.keras.layers.Dropout(.25)(dense)
output_layer = tf.keras.layers.Dense(output_classes, activation = "softmax")(final)
model = tf.keras.Model(inputs = input_layer, outputs = output_layer)

model.compile(optimizer = 'adam',
              loss = 'categorical_crossentropy',
              #   metrics = ['accuracy'])
              metrics = [tf.keras.metrics.CategoricalAccuracy()])



with open(str(version_dir / "Summery.txt"), "w") as f:
    f.write(f"Seed: {seed}\n\n")
    with redirect_stdout(f):
        model.summary()

# steps_per_epoch = int(data_size / batch_size) + 1
try:
    history = model.fit(train_set, 
                        batch_size = batch_size, 
                        epochs = epochs,
                        # steps_per_epoch = steps_per_epoch,
                        validation_data = valid_set,
                        callbacks = [cp_callback])
except Exception as E:
    print(E)
    print(f"Error type {type(E)}")
    exit()


history_data = history.history
pandas.DataFrame(history_data).to_csv(str(version_dir / "ModelSteps.csv"))


# plt.plot(history_data["loss"])
# plt.plot(history_data["val_loss"])
# plt.xlabel("Epochs")
# plt.legend(["Loss", "Valid Loss"])
# plt.ylim(top = 10, bottom = 0)
# plt.savefig(str(version_dir / "Loss_Graph.png"))

# plt.close()


# plt.plot(history_data["categorical_accuracy"])
# plt.plot(history_data["val_categorical_accuracy"])
# plt.xlabel("Epochs")
# plt.legend(["Accuracy", "Valid Accuracy"])
# # plt.show()
# plt.savefig(str(version_dir / "Accuracy_Graph.png"))

print(f"\nOutput files to\n{version_dir}")
