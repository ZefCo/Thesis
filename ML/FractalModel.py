# https://towardsdatascience.com/convolutional-neural-networks-understanding-and-organizing-your-data-set-ba3e8b4086cb
# https://datagy.io/python-decorators/
# https://gitpython.readthedocs.io/en/stable/intro.html

# https://stackoverflow.com/questions/26134026/how-to-get-the-current-checked-out-git-branch-name-through-pygit2

# https://gitpython.readthedocs.io/en/stable/intro.html

import os
# os.environ["GET_PYTHON_REFRESH"] = "quiet"
import git
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

cwd = pathlib.Path.cwd()

report = git.Repo(search_parent_directories=True)
branch = report.active_branch

image_dir = cwd / "FractalImageEvI"
model_dir = cwd / "FractalModels" / "IvE"
version_num = len(next(os.walk(model_dir))[1]) + 1
version_dir = model_dir / f"version_{branch}_{version_num}"
version_dir.mkdir(parents = True, exist_ok = True)
output_classes = 2

# image_dir = cwd / "FractalImageCvIvU"
# model_dir = cwd / "FractalModels" / "CvIvU"
# version_num = len(next(os.walk(model_dir))[1])
# version_dir = model_dir / f"version_{version_num}"
# output_classes = 4

seed = random.randint(1000000, 9000000)


batch_size = 200
epochs = 100

pngs = 0
for folder in os.listdir(image_dir):
    for file in os.listdir(image_dir / folder):
        if file.endswith("png"):
            pngs += 1


train_set = tf.keras.preprocessing.image_dataset_from_directory(str(image_dir), image_size = (64, 64), seed = seed, subset = "training", validation_split = 0.2, batch_size = batch_size, label_mode = "categorical", color_mode = "rgba")
valid_set = tf.keras.preprocessing.image_dataset_from_directory(str(image_dir), image_size = (64, 64), seed = seed, subset = "validation", validation_split = 0.2, batch_size = batch_size, label_mode = "categorical", color_mode = "rgba")

steps_per_epoch = int(pngs / batch_size) + 1
# print(steps_per_epoch)


input_layer = tf.keras.Input(shape = (64, 64, 4))
x = tf.keras.layers.Conv2D(64, (1, 1), activation = "gelu", kernel_regularizer = tf.keras.regularizers.l2(l = 0.001))(input_layer)
x = tf.keras.layers.Dropout(.25)(x)
x = tf.keras.layers.BatchNormalization()(x)
x = tf.keras.layers.Conv2D(64, (1, 1), activation = "gelu", kernel_regularizer = tf.keras.regularizers.l2(l = 0.001))(x)
x = tf.keras.layers.Dropout(.25)(x)
x = tf.keras.layers.BatchNormalization()(x)
x = tf.keras.layers.Conv2D(64, (1, 1), activation = "gelu", kernel_regularizer = tf.keras.regularizers.l2(l = 0.001))(x)
x = tf.keras.layers.Dropout(.25)(x)
x = tf.keras.layers.BatchNormalization()(x)
x = tf.keras.layers.Flatten()(x)
output_layer = tf.keras.layers.Dense(output_classes, activation = "softmax")(x)

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
                        validation_data = valid_set)
except Exception as E:
    print(type(E))
    print(E)


history_data = history.history
pandas.DataFrame(history_data).to_csv(str(version_dir / "ModelSteps.csv"))


plt.plot(history_data["loss"])
plt.plot(history_data["val_loss"])
plt.xlabel("Epochs")
plt.legend(["Loss", "Valid Loss"])
plt.ylim(top = 10, bottom = 0)
plt.savefig(str(version_dir / "Loss_Graph.png"))

plt.close()


plt.plot(history_data["categorical_accuracy"])
plt.plot(history_data["val_categorical_accuracy"])
plt.xlabel("Epochs")
plt.legend(["Accuracy", "Valid Accuracy"])
# plt.show()
plt.savefig(str(version_dir / "Accuracy_Graph.png"))

print(f"\nOutput files to {version_dir}")
