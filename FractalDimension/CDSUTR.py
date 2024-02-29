import pandas
import numpy as np
import os
import pathlib
cwd = pathlib.Path.cwd()
import MomentCalculations as MC
import Heatmaps as heat
import TimeEmbedding as TE
import pickle
import GeneClass as Gene
from typing import Tuple


def main():
    '''
    '''
    gene: Gene.Gene
    with open(cwd.parent / "Data_Files\Primates\Genetics\Homo_sapiens\CDS_UTR\Test_Data.pkl", "rb") as file:
        data: dict = pickle.load(file)

    for name, gene in data.items():
        print(gene.ncibname)
        print(gene.strand)
        print(gene.chrm)
        # print(gene.urls)
        for url in gene.urls:
            print(url)
        # try:
        #     print(gene.utr3_seq[0])
        #     print("----")
        # except Exception as e:
        #     print("cannot get utr3")
        # try:
        #     print(gene.cds_seq[0])
        #     print("----")
        # except Exception as e:
        #     print("Cannot get cds")
        # try:
        #     print(gene.utr5_seq[0])
        #     print("----")
        # except Exception as e:
        #     print("cannot get utr5")
        # try:
        #     print(gene.exon_seq[0])
        #     print("----")
        # except Exception as e:
        #     print("Cannot get exon")
        # try:
        #     print(gene.intron_seq[0])
        #     print(gene.intron_cords[0])
        #     print("----")
        # except Exception as e:
        #     print("Cannot get intron")

        print("\n\n________")
        # print(gene.cds_seq)
        # print(gene.utr3_seq)
        # print(gene.utr5_seq)


if __name__ in "__main__":
    main()