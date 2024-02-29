import pandas
import re
import pathlib
cwd = pathlib.Path.cwd()
# import Permutations as perm
import GeneClass as Gene
import pickle
import requests
import random
import os
import GeneCollection



def main():
    '''
    '''
    "D:\Coding\Thesis"
    output_file = cwd.parent / "Data_Files" / "Primates" / "Genetics" / "Homo_sapiens" / "CDS_UTR" / "Test_Data.pkl"
    GeneCollection.hg19_sequences(cwd.parent / "Data_Files\Gene_Files\HG19\Known_Genes_hg19.csv",
                                  output_file, 
                                  ref_track="ncib",
                                  gene_start = 0,
                                  gene_stop = 10,
                                  species = "Homo_sapiens",
                                  genome = "hg19")


if __name__ in "__main__":
    main()