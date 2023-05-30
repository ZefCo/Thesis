import pathlib
cwd = pathlib.Path.cwd()
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


def History_1D():
    '''
    '''
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

    all_acc = go.Figure()
    all_los = go.Figure()

    for m, model in enumerate(pandas.unique(models["Date"])):
        local_model = models[models["Date"] == model]
        # cat_acc.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["categorical_accuracy"])) #, name = model))
        # val_acc.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_categorical_accuracy"])) #, name = model))
        # cat_los.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["loss"])) #, name = model))
        # val_los.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_loss"])) #, name = model))

        all_acc.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["categorical_accuracy"], name = f"Training Acc {m}"))
        all_acc.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_categorical_accuracy"], name = f"Validation Acc {m}"))

        all_los.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["loss"], name = f"Training Loss {m}"))
        all_los.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_loss"], name = f"Validation Loss {m}"))

        # hist.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_loss"], name = model))
        # hist.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_categorical_accuracy"], name = model))
        # hist.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["loss"], name = model))

    # cat_acc.update_xaxes(title = "Epoch")
    # val_acc.update_xaxes(title = "Epoch")
    # cat_los.update_xaxes(title = "Epoch")
    # val_los.update_xaxes(title = "Epoch")
    all_acc.update_xaxes(title = "Epoch")
    all_los.update_xaxes(title = "Epoch")

    # cat_acc.update_layout(title = "Categorical Accuracy")
    # val_acc.update_layout(title = "Validation Accuracy")
    # cat_los.update_layout(title = "Categorical Loss")
    # val_los.update_layout(title = "Validation Loss")

    all_acc.update_layout(title = "Accuracy")
    all_los.update_layout(title = "Loss")

    # cat_acc.show()
    # cat_los.show()
    # val_acc.show()
    # val_los.show()

    all_acc.show()
    all_los.show()
    # loss = go.Figure()
    # for model in pan



def History_F():
    '''
    '''
    fractal_models = cwd / "FractalModels" / "IvE"
    models = pandas.DataFrame()

    model_count = 0

    for parent, folder, files in os.walk(fractal_models):
        for file in files:
            if file in "ModelSteps.csv":
                model_count += 1
                local_model = pandas.read_csv(pathlib.Path(parent, file), header = 0)
                local_model["ModelIndex"] = f"Model_{model_count}"
                local_model["ModelIndex"] = pandas.Categorical(local_model["ModelIndex"])
                local_model["Epoch"] = local_model.index
        
                models = pandas.concat([models, local_model])
                # print(pathlib.Path(parent, file))
    
    models.pop("Unnamed: 0")
    models = models.reset_index()
    models.pop("index")
    # # models["Epoch"]
    print(models)
    # # print(models.shape)

    all_acc = go.Figure()
    all_los = go.Figure()

    for m, model in enumerate(pandas.unique(models["ModelIndex"])):
        local_model = models[models["ModelIndex"] == model]

        all_acc.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["categorical_accuracy"], name = f"Training Acc {m}"))
        all_acc.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_categorical_accuracy"], name = f"Validation Acc {m}"))

        all_los.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["loss"], name = f"Training Loss {m}"))
        all_los.add_trace(go.Scatter(x = local_model["Epoch"], y = local_model["val_loss"], name = f"Validation Loss {m}"))

    all_acc.update_xaxes(title = "Epoch")
    all_los.update_xaxes(title = "Epoch")

    all_acc.update_layout(title = "Accuracy")
    all_los.update_layout(title = "Loss")

    all_acc.show()
    all_los.show()



if __name__ in "__main__":
    History_1D()
    # History_F()