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


def main():
    '''
    Look i knew I keep doing this but never import the things... the project got rather intense and I never was able to keep up with all
    the demands while keeping everything clean. There should be an import script so I don't recrete the thing over and over again, and probably
    a few others, but oh well.
    '''


def moment(n: int = 1, regions = None):
    '''
    Takes a region and finds the moment of that region.

    m = [sum(region)**n / sum(area)]**(1/n)

    It will have to import the whole data, then isolate a single area and count that. Maybe I shoudl count by region, add them, then calculate each moment.
    This would divide the space into R regions
    '''
    u: float = 1 / n  # think of this as an upside down n


if __name__ in "__main__":
    main()