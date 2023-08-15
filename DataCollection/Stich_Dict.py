import pandas
import re
import pathlib
cwd = pathlib.Path.cwd()
# import Permutations as perm
import GeneClass as Gene
import pickle
import requests


def main():
    '''
    '''
    with open(pathlib.Path(f"/media/ethanspeakman/Elements/GeneData/Known_Genes_hg19_NCBIGene_DICT.pkl"), "rb") as p:
        data = pickle.load(p)

    print(data["NM_001351428.2"])
    print(data["NM_001351428.2"].full_seq)


if __name__ in "__main__":
    main()