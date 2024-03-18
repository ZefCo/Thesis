import os
# os.environ["GET_PYTHON_REFRESH"] = "quiet"
import git
import pandas
import pathlib
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

cwd = pathlib.Path.cwd()


# I fucking hate him. I hate Dr G and don't want to do this anymore.

path_to_model_folder = cwd / "FractalModels" "IvE" / "version_slurm_24" / "Model"
path_to_data_folder = cwd.parent / "Data_Files" / "FuckDrG" / "Data_Set_1_histogram_6mer_KTA"

model.load_weights(path_to_model_folder)