from importlib.resources import path
import pandas
import requests
import RQuery
import CommonMethods as CM
import blat_api as bpi
import ucsc_restapi as upi
from FusionClass import FusionGene as Fusions
from GeneClass import Gene
import pathlib



def running_totals(add_var, var_str):
    '''
    '''
    add_var += 1
    print(f"Current Running Total of {var_str}: {add_var}")

    return add_var



def main(infile: str or pathlib.Path, min_length: int = 1000):
    '''
    '''
    with open(infile) as datafile:
        fusion_data = pandas.read_csv(datafile, header=0)

    fusion_data = fusion_data[fusion_data["SeqLen"] >= min_length]

    rows, _ = fusion_data.shape

    unknowns: int = 0
    cis: int = 0
    ta: int = 0
    te: int = 0

    for row in range(0, rows):
        row_of_interest = fusion_data.iloc[row, :].copy()

        hgene, tgene = row_of_interest["Hgene"], row_of_interest["Tgene"]
        henst, tenst = row_of_interest["Henst"], row_of_interest["Tenst"]
        hchrm, tchrm = row_of_interest["Hchr"], row_of_interest["Tchr"]
        hstrand, tstrand = row_of_interest["Hstrand"], row_of_interest["Tstrand"]
        seq = row_of_interest["Seq"]

        fusion = Fusions(hgene = hgene, tgene = tgene, henst = henst, tenst = tenst, hchrm = hchrm, tchrm = tchrm, hstrand = hstrand, tstrand = tstrand, seq = seq)
        print(f"####\n{fusion.hgene}_{fusion.tgene}\n{fusion.henst}_{fusion.tenst}")
        fusion.blat()

        fusion.classify()

        fusion.distance_measure()

        if fusion.classification in "Unknown":
            unknowns = running_totals(unknowns, "Unknown")
        elif fusion.classification in "C-SAG":
            cis = running_totals(cis, "Cis-SAG")
        elif fusion.classification in "T-E":
            te = running_totals(te, "Trans-E")
        elif fusion.classification in "T-A":
            ta = running_totals(ta, "Trans-A")

        if fusion.clean_blat and fusion.clean_enst:
            outfile = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / f"UT_BE_min{min_length}.csv"
        elif fusion.clean_blat and not fusion.clean_enst:
            outfile = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / f"EnstFailed.csv"
        elif not fusion.clean_blat and not fusion.clean_enst:
            outfile = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "BlatFailed.csv"

        # fusion.write_to_csv(outfile = outfile)

        print(f"\n~~Finished Row {row} of {rows}~~\n####\n")

        # exit()




if __name__ in '__main__':
    fusion_file = pathlib.Path.cwd().parent / "Data_Files" / "UTData_cds.csv"
    main(fusion_file)