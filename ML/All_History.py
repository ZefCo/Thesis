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

epochs = None



def main():
    '''
    '''
    # History_1D()
    # History_F()
    SelectedHistory(GMD_320_183921 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-20-18-39-21",
                    GMD_320_18437 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-20-18-43-7",
                    GMD_320_184338 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-20-18-43-38",
                    GMD_321_8557 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-21-8-5-57-994627",
                    GMD_321_21467 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-21-21-46-7-886939",
                    GMD_322_74153 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-7-41-53-606828",
                    GMD_322_74917 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-7-49-17-995433",
                    GMD_322_13120 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-13-12-0-637765",
                    GMD_322_131356 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-13-13-56-583206",
                    GMD_322_133756 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-13-37-56-620653",
                    GMD_322_14720 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-14-7-20-431707",
                    GMD_322_17955 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-17-9-55-837319",
                    GMD_322_172844 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-17-28-44-920291",
                    GMD_322_181335 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-18-13-35-948525",
                    GMD_322_181532 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-18-15-32-995204",
                    GMD_322_194111 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-22-19-41-11-248177",
                    GMD_325_123957 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-25-12-39-57-34508",
                    GMD_325_124528 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-25-12-45-28-551325",
                    GMD_326_11223 = cwd / "FractalModels" / "IvE" / "version_slurm_2024-3-26-11-22-3-83685",
        
                    # GMD_94 = cwd / "FractalModels" / "IvE" / "version_slurm_94",
                    # GMD_95 = cwd / "FractalModels" / "IvE" / "version_slurm_95",
                    # GMD_96 = cwd / "FractalModels" / "IvE" / "version_slurm_96",
                    # GMD_97 = cwd / "FractalModels" / "IvE" / "version_slurm_97",
                    # GMD_98 = cwd / "FractalModels" / "IvE" / "version_slurm_98",
                    # GMD_99 = cwd / "FractalModels" / "IvE" / "version_slurm_99",
                    # SLSGHis = cwd / "FractalModels" / "IvE" / "version_slurm_47",
                    # DS1_hist = cwd / "FractalModels" / "IvE" / "version_slurm_23",
                    # DS1 = cwd / "FractalModels" / "IvE" / "version_slurm_22",
                    # DS2 = cwd / "FractalModels" / "IvE" / "version_slurm_24",
                    # mer9_v1 = cwd / "FractalModels" / "IvE" / "version_slurm_50",
                    # mer8_v1 = cwd / "FractalModels" / "IvE" / "version_slurm_51",
                    # mer8_v2 = cwd / "FractalModels" / "IvE" / "version_slurm_52",
                    # mer9_v2 = cwd / "FractalModels" / "IvE" / "version_slurm_55",
                    # new_v70 = cwd / "FractalModels" / "IvE" / "version_slurm_70",
                    # new_71 = cwd / "FractalModels" / "IvE" / "version_slurm_71",
                    # new_73 = cwd / "FractalModels" / "IvE" / "version_slurm_73",
                    # new_74 = cwd / "FractalModels" / "IvE" / "version_slurm_74",
                    # new_75 = cwd / "FractalModels" / "IvE" / "version_slurm_75",
                    # new_76 = cwd / "FractalModels" / "IvE" / "version_slurm_76",
                    # new_77 = cwd / "FractalModels" / "IvE" / "version_slurm_77",
                    # new_78 = cwd / "FractalModels" / "IvE" / "version_slurm_78",
                    # new_79 = cwd / "FractalModels" / "IvE" / "version_slurm_79",
                    # new_80 = cwd / "FractalModels" / "IvE" / "version_slurm_80",
                    # new_81 = cwd / "FractalModels" / "IvE" / "version_slurm_81",
                    # new_82 = cwd / "FractalModels" / "IvE" / "version_slurm_82",
                    # new_83 = cwd / "FractalModels" / "IvE" / "version_slurm_83",
                    # new_85 = cwd / "FractalModels" / "IvE" / "version_slurm_85",
                    # new_86 = cwd / "FractalModels" / "IvE" / "version_slurm_86",
                    # new_87 = cwd / "FractalModels" / "IvE" / "version_slurm_87",
                    # new_88 = cwd / "FractalModels" / "IvE" / "version_slurm_88",
                    # new_89 = cwd / "FractalModels" / "IvE" / "version_slurm_89",

                    line_of_best_epochs = epochs)


def least_squares(x_points: np.array, y_points: np.array):
    '''
    '''
    ave_x = np.mean(x_points)
    ave_y = np.mean(y_points)

    # print(x_points)
    # print(y_points)

    m = np.dot(x_points - ave_x, y_points - ave_y) / np.dot(x_points - ave_x, x_points - ave_x)
    b = ave_y - (m * ave_x)

    ls = {e: m*e + b for e in range(0, x_points.max() + 1)}

    return ls, m, b


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


def SelectedHistory(line_of_best_epochs = None, *args, **kwargs):
    '''
    Right now this cannot combine args and kwargs, so use only one at a time.
    '''
    models = pandas.DataFrame()

    if (len(args) > 0) and (len(kwargs) < 1):

        model_count = 0

        for path in args:
            for _, _, files in os.walk(path):
                # print(file)
                for file in files:
                    if file in "ModelSteps.csv":
                        model_count += 1
                        local_model = pandas.read_csv(pathlib.Path(path, file), header = 0)
                        local_model["ModelIndex"] = f"Model_{model_count}"
                        local_model["ModelIndex"] = pandas.Categorical(local_model["ModelIndex"])
                        local_model["Epoch"] = local_model.index + 1
                
                        models = pandas.concat([models, local_model])
                        # print(pathlib.Path(parent, file))

    elif (len(args) < 1) and (len(kwargs) > 0):
        for model, path in kwargs.items():
            for _, _, files in os.walk(path):
                # print(file)
                for file in files:
                    if file in "ModelSteps.csv":
                        local_model = pandas.read_csv(pathlib.Path(path, file), header = 0)
                        local_model["ModelIndex"] = model
                        local_model["ModelIndex"] = pandas.Categorical(local_model["ModelIndex"])
                        local_model["Epoch"] = local_model.index + 1
                
                        models = pandas.concat([models, local_model])
                        # print(pathlib.Path(parent, file))

    
    models.pop("Unnamed: 0")
    models = models.reset_index()
    models.pop("index")
    # # models["Epoch"]
    # print(models)
    # # print(models.shape)
    # epochs = models["Epoch"].unique()
    # print(epochs)

    all_acc = go.Figure()
    all_los = go.Figure()

    for m, model in enumerate(pandas.unique(models["ModelIndex"])):
        # print(model)
        local_model: pandas.DataFrame = models[models["ModelIndex"] == model]
        local_model = local_model.set_index("Epoch")
        # print(local_model)

        all_acc.add_trace(go.Scatter(x = local_model.index, y = local_model["categorical_accuracy"], name = f"{model}: TA", legendgroup = model))
        all_acc.add_trace(go.Scatter(x = local_model.index, y = local_model["val_categorical_accuracy"], name = f"{model}: VA", legendgroup = model))

        all_los.add_trace(go.Scatter(x = local_model.index, y = local_model["loss"], name = f"{model}: TL", legendgroup = model))
        all_los.add_trace(go.Scatter(x = local_model.index, y = local_model["val_loss"], name = f"{model}: VL", legendgroup = model))

        if line_of_best_epochs is not None:
            epochs = local_model.index

            acc_range = local_model.loc[line_of_best_epochs:max(epochs) + 1, "categorical_accuracy"].to_numpy()
            val_range = local_model.loc[line_of_best_epochs:max(epochs) + 1, "val_categorical_accuracy"].to_numpy()
            
            x_points = np.array([x for x in range(line_of_best_epochs, max(epochs) + 1)])
            # print(x_points)
            # print(acc_range)
            # print(val_range)

            acc_ls, acc_m, acc_b = least_squares(x_points, acc_range)
            val_ls, val_m, val_b = least_squares(x_points, val_range)

            # print(acc_ls)
            
            all_acc.add_trace(go.Scatter(x = tuple(acc_ls.keys()), y = tuple(acc_ls.values()), name = rf"{model} LoBF starting @ {line_of_best_epochs}<br>y = {round(acc_m, 6)}x + {round(acc_b, 6)}", legendgroup = model))
            all_acc.add_trace(go.Scatter(x = tuple(val_ls.keys()), y = tuple(val_ls.values()), name = rf"{model} LoBF starting @ {line_of_best_epochs}<br>y = {round(val_m, 6)}x + {round(val_b, 2)}", legendgroup = model))

            # Putting this here in case Dr. G says anything, then I can show him it doesn't work
            # all_acc.add_annotation(x = 0, y = tuple(acc_ls.values()), text = f"y = {round(acc_m, 6)}x + {round(acc_b, 6)}")
            # all_acc.add_annotation(x = 0, y = tuple(val_ls.values()), text = f"y = {round(val_m, 6)}x + {round(val_b, 6)}")


    all_acc.update_xaxes(title = "Epoch")
    all_los.update_xaxes(title = "Epoch")

    all_acc.update_layout(title = "Accuracy", legend = dict(orientation = "h"))
    all_los.update_layout(title = "Loss")

    all_acc.show()
    # all_los.show()
    # print(str(cwd / "ML_Acc.html"))
    all_acc.write_html(str(cwd / "ML_Acc.html"))
    all_los.write_html(str(cwd / "ML_Losses.html"))



if __name__ in "__main__":
    main()
