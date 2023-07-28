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
import pandas
import pickle



def main():
    '''
    '''

    gene = import_data()
    # print(gene)



def import_data():
    '''
    Import data and grab one gene. Don't care which one. Later I'll figure out how to do isoforms
    '''

    data_file = cwd.parent / "Data_Files" / "Gene_Files" / "Hg19" / "Known_Genes_hg19_ensGene.pkl"

    with open(data_file, "rb") as p:
        data: dict = pickle.load(p)

    print(data)

    # gene_keys = data.keys()

    # data = data[gene_keys[1]]

    # return data




if __name__ in "__main__":
    main()