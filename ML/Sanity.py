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

df_train = np.random.normal(size=(17980, 380))
df_validation = np.random.normal(size=(17980, 380))

print("Before Adding Dummy Dimension")
# print(df_train.shape)
print(type(df_train))
# print(df_train)
df_train = df_train[..., None]
print(f"\n\nAfter adding Dummy Dimension")
# print(df_train.shape)
print(type(df_train))
# print(df_train)
df_validation = df_validation[..., None]

# y_train = np.random.normal(size=(17980, 1))
# y_validation = np.random.normal(size=(17980, 1))