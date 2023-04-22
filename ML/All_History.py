import pathlib
import os
import re
# from tensorflow.keras.models import Sequential, load_module
# from tensorflow.keras.layers import Conv1D
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import accuracy_score
# import tensorflow as tf
# print(tf.__version__)
# exit()
import random
import numpy as np
import pandas
import plotly.graph_objects as go
import datetime
# from matplotlib import pyplot as plt
# from contextlib import redirect_stdout
# import keras


cwd = pathlib.Path.cwd()

model_histories = cwd / "Model_History"


models = pandas.DataFrame()
model_files = list(model_histories.rglob("*.csv"))

for model in model_files:
    # print(model.name)
    # print(type(model.name))
    local_model = pandas.read_csv(model, header = 0)
    local_model["Date"] = str(model.name)
    local_model["Date"] = pandas.Categorical(local_model["Date"])
    local_model["Epoch"] = local_model.index

    models = pandas.concat([models, local_model])

models.pop("Unnamed: 0")
models = models.reset_index()
models.pop("index")
# models["Epoch"]
print(models)
# print(models.shape)


cat_acc = go.Figure()
cat_los = go.Figure()
val_acc = go.Figure()
val_los = go.Figure()
for model in pandas.unique(models["Date"]):
    local_model = models[models["Date"] == model]
    cat_acc.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["categorical_accuracy"], name = model))
    val_acc.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_categorical_accuracy"], name = model))
    cat_los.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["loss"], name = model))
    val_los.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_loss"], name = model))
    # hist.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_loss"], name = model))
    # hist.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_categorical_accuracy"], name = model))
    # hist.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["loss"], name = model))

cat_acc.update_xaxes(title = "Epoch")
val_acc.update_xaxes(title = "Epoch")
cat_los.update_xaxes(title = "Epoch")
val_los.update_xaxes(title = "Epoch")
cat_acc.update_layout(title = "Categorical Accuracy")
val_acc.update_layout(title = "Validation Accuracy")
cat_los.update_layout(title = "Categorical Loss")
val_los.update_layout(title = "Validation Loss")
cat_acc.show()
cat_los.show()
val_acc.show()
val_los.show()


# loss = go.Figure()
# for model in pan