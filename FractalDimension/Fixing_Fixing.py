import pathlib
cwd = pathlib.Path.cwd()
import os
import re
import numpy as np
import glob
from plotly import graph_objects as go
import timeit
import DistanceClass as distance
# from scipy.spatial import distance
# from statsmodels.graphics.tsaplots import plot_acf
import pandas
import pickle
import itertools
import matplotlib.pyplot as plt
import GeneClass as Gene
import random
import Heatmaps as hm
from typing import Tuple
from MomentCalculations import moments_v2
from MomentCalculations import moments
from MomentCalculations import _unrenormalize


with open(cwd / "Dicts" / "Exon_6mer.pkl", "rb") as file:
    something = pickle.load(file)

# with open(cwd / "Dicts_Fixed" / "Exon_6mer.pkl", "rb") as file:
#     somethingelse = pickle.load(file)


print(something)
# print(somethingelse)