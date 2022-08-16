import pandas
import requests
import RQuery
import CommonMethods as CM
import blat_api as bpi
import ucsc_restapi as upi
from FusionClass import FusionGene as Fusions
from GeneClass import Gene
import pathlib



def main(infile):
    '''
    '''
    with open(infile) as datafile:
        fusion_data = pandas.read_csv(datafile, header=0)

    rows, _ = fusion_data.shape

    for row in range(310, rows):
        row_of_interest = fusion_data.iloc[row, :].copy()

        hgene, tgene = row_of_interest["Hgene"], row_of_interest["Tgene"]
        henst, tenst = row_of_interest["Henst"], row_of_interest["Tenst"]
        hchrm, tchrm = row_of_interest["Hchr"], row_of_interest["Tchr"]
        hstrand, tstrand = row_of_interest["Hstrand"], row_of_interest["Tstrand"]
        seq = row_of_interest["Seq"]

        fusion = Fusions(hgene = hgene, tgene = tgene, henst = henst, tenst = tenst, hchrm = hchrm, tchrm = tchrm, hstrand = hstrand, tstrand = tstrand, seq = seq)
        print(f"####\n{fusion.gnames}\n{fusion.enames}")
        fusion.blat()

        # How to handel JSON Decoder errors: skip? Dump?

        print(f"\n~~Finished Row {row} of {rows}~~\n####\n")

        # exit()


if __name__ in '__main__':
    fusion_file = pathlib.Path.cwd().parent / "Data_Files" / "UTData_cds.csv"
    main(fusion_file)