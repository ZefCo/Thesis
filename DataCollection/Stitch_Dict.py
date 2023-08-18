import pandas
import re
import pathlib
cwd = pathlib.Path.cwd()
# import Permutations as perm
import GeneClass as Gene
import pickle
import requests
import os


def main():
    '''
    '''
    # with open(pathlib.Path(f"/media/ethanspeakman/Elements/GeneData/Known_Genes_hg19_NCBIGene_DICT.pkl"), "rb") as p:
    #     data = pickle.load(p)

    # print(data["NM_001351428.2"])
    # print(data["NM_001351428.2"].full_seq)


def select_data(n = 20_000):
    '''
    Randomly selects a certain ammount of the genes from the data. By default the amount if 20,000. This is then split into two even sub groups.
    '''


def stich_dict(directory_path: pathlib.Path):
    '''
    '''

    files = []
    z = {}

    for file in os.listdir(directory_path):
        if file.endswith(".pkl"):
            files.append(file)
    
    for file in files:
        with open(file, "rb") as p:
            x = pickle.load(p)

        z = z | x


    with open(directory_path / "Dict_Of_All_Genes.pkl", 'wb') as p:
        pickle.dump(z, p)



if __name__ in "__main__":
    main()