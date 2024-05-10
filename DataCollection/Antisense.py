import re
import ucsc_restapi as upi
import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import pickle
import random
import GeneClass as Gene
from GData import random_keys


def main():
    '''
    Beacuse I don't understand antisense genes.
    '''
    element_drive = pathlib.Path("F:\Gene_Data_Sets")
    gene_drive = pathlib.Path("F:\GeneData")
    file = element_drive / "Data_Set_1_cleaned_dict.pkl"

    minus_dict = {}
    minus_out = gene_drive / "Minus.pkl"
    positive_dict = {}
    positive_out = gene_drive / "Positive.pkl"

    with open(file, "rb") as pfile:
        data_dict: dict = pickle.load(pfile)
    
    selected_data = random_keys(data_dict, 500)
    gene: Gene.Gene

    for name, gene in selected_data.items():
        if gene.strand in "+":
            minus_dict[name] = gene
        elif gene.strand in "-":
            positive_dict[name] = gene

    with open(minus_out, "wb") as file:
        pickle.dump(minus_dict, file)
    with open(positive_out, "wb") as file:
        pickle.dump(positive_dict, file)


if __name__ in "__main__":
    main()