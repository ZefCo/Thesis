from importlib.resources import path
import pandas
import requests
import RQuery
import CommonMethods as CM
import blat_api as bpi
import ucsc_restapi as upi
from FusionClass import FusionGene as Fusions
from GeneClass import Gene
from BlatClass import Blat
from ScoreingClasses import DissSimilarityScore
import pathlib



def running_totals(add_var, var_str):
    '''
    '''
    add_var += 1
    print(f"Current Running Total of {var_str}: {add_var}")

    return add_var



def gapSearch(infile: str or pathlib.Path):
    '''
    '''
    with open(infile) as datafile:
        fusion_data = pandas.read_csv(datafile, header=0)

    rows, _ = fusion_data.shape

    for row in range(rows):
        row_of_interest = fusion_data.iloc[row, :].copy()

        hgene, tgene = row_of_interest["hgene"], row_of_interest["tgene"]
        henst, tenst = row_of_interest["henst"], row_of_interest["tenst"]
        hchrm, tchrm = row_of_interest["hchrm"], row_of_interest["tchrm"]
        hstrand, tstrand = row_of_interest["hstrand"], row_of_interest["tstrand"]
        seq, ctype, source = row_of_interest["seq"], row_of_interest["Ctype"], row_of_interest["Source"]

        classification, short_distance, head2tailDistance = row_of_interest["classification"], row_of_interest["shortDistance"], row_of_interest["head2tailDistance"]

        fusion = Fusions(hgene = hgene, tgene = tgene, henst = henst, tenst = tenst, hchrm = hchrm, tchrm = tchrm, hstrand = hstrand, tstrand = tstrand, seq = seq, ctype = ctype, source = source)

        # Fuck I forgot that these things go out as lists but come back in as strings
        convert_list = ["hexonStarts", "texonStarts", "hexonEnds", "texonEnds", "hexonFrames", "texonFrames", "hblockSizes", "tblockSizes", "htStarts", "ttStarts"]
        for strtup in convert_list:
            row_of_interest[strtup] = CM.convert2list(row_of_interest[strtup])

        head_gene = Gene(row_of_interest["hgene"],
                        gname = row_of_interest["hgname"],
                        txStart = row_of_interest["htxStart"],
                        txEnd = row_of_interest["htxEnd"],
                        cdsStart =row_of_interest["hcdsStart"],
                        cdsEnd = row_of_interest["hcdsEnd"],
                        exonCount = row_of_interest["hexonCount"],
                        exonStarts = row_of_interest["hexonStarts"],
                        exonEnds = row_of_interest["hexonEnds"],
                        exonFrames = row_of_interest["hexonFrames"])

        tail_gene = Gene(row_of_interest,
                        gname = row_of_interest["tgname"],
                        txStart = row_of_interest["ttxStart"],
                        txEnd = row_of_interest["ttxEnd"],
                        cdsStart = row_of_interest["tcdsStart"],
                        cdsEnd = row_of_interest["tcdsEnd"],
                        exonCount = row_of_interest["texonCount"],
                        exonStarts = row_of_interest["texonStarts"],
                        exonEnds = row_of_interest["texonEnds"],
                        exonFrames = row_of_interest["texonFrames"])

        head_blat = Blat(qName = row_of_interest["hgene"],
                        tStart = row_of_interest["htStart"],
                        tEnd = row_of_interest["htEnd"],
                        blockCount = row_of_interest["hblockCount"],
                        blockSizes = row_of_interest["hblockSizes"],
                        tStarts = row_of_interest["htStarts"])

        tail_blat = Blat(qName = row_of_interest["tgene"],
                        tStart = row_of_interest["ttStart"],
                        tEnd = row_of_interest["ttEnd"],
                        blockCount = row_of_interest["tblockCount"],
                        blockSizes = row_of_interest["tblockSizes"],
                        tStarts = row_of_interest["ttStarts"])



        fusion.second_import(classification, short_distance, head2tailDistance, head_gene, tail_gene, head_blat, tail_blat)
        print(f"####\n{fusion.hgene}_{fusion.tgene}\n{fusion.henst}_{fusion.tenst}")

        score: DissSimilarityScore = DissSimilarityScore(fusion)
        score.slip_junction()
        print(f"{score.Hslip}{score.Tslip}")
        score.write_slip(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "SlipJunctionWSeq_fs.csv")

        print(f"\n~~Finished Row {row} of {rows}~~\n####\n")

        # if row == 10:
        #     exit()



def databaseBuild(infile: str or pathlib.Path, min_length: int = 100, start_index: int = 0):
    '''
    Row 101 of the Data hit a NoneType attribute error at scoring classes line 125.
    '''
    with open(infile) as datafile:
        fusion_data = pandas.read_csv(datafile, header=0)

    fusion_data = fusion_data[fusion_data["SeqLen"] >= min_length]

    rows, _ = fusion_data.shape

    unknowns: int = 0
    cis: int = 0
    ta: int = 0
    te: int = 0

    for row in range(start_index, rows + 1):
        row_of_interest = fusion_data.iloc[row, :].copy()

        hgene, tgene = row_of_interest["Hgene"], row_of_interest["Tgene"]
        henst, tenst = row_of_interest["Henst"], row_of_interest["Tenst"]
        hchrm, tchrm = row_of_interest["Hchr"], row_of_interest["Tchr"]
        hstrand, tstrand = row_of_interest["Hstrand"], row_of_interest["Tstrand"]
        seq, ctype, source = row_of_interest["Seq"], row_of_interest["Ctype"], row_of_interest["Source"]

        fusion = Fusions(hgene = hgene, tgene = tgene, henst = henst, tenst = tenst, hchrm = hchrm, tchrm = tchrm, hstrand = hstrand, tstrand = tstrand, seq = seq, ctype = ctype, source = source)
        print(f"####\n{fusion.hgene}_{fusion.tgene}\n{fusion.henst}_{fusion.tenst}")

        fusion.blat()
        fusion.classify()
        fusion.distance_measure()


        if (fusion.classification is None) or (fusion.classification in "Unknown"):
            unknowns = running_totals(unknowns, "Unknown")
        elif fusion.classification in "C-SAG":
            cis = running_totals(cis, "Cis-SAG")
        elif fusion.classification in "T-E":
            te = running_totals(te, "Trans-E")
        elif fusion.classification in "T-A":
            ta = running_totals(ta, "Trans-A")


        if fusion._clean_blat and fusion._clean_enst:
            outfusion: pathlib.Path = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / f"UT_BE_min{min_length}_fs.csv"
            outscore: pathlib.Path = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / f"Fusion_Scores_min{min_length}_fs.csv"
            
            score: DissSimilarityScore = DissSimilarityScore(fusion)
            score.score()

            score.write_score(outfile = outscore)
        
        elif fusion._clean_blat and not fusion._clean_enst:
            outfusion: pathlib.Path = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / f"EnstFailed_fs.csv"
        
        elif not fusion._clean_blat and not fusion._clean_enst:
            outfusion: pathlib.Path = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "BlatFailed_fs.csv"

        fusion.write_to_database(outfile = outfusion)
        # exit()


        print(f"\n~~Finished Row {row} of {rows}~~\n####\n")

        if row > 100:
            exit()
        # exit()

        # if row == 1:
        #     exit()


def scoreing(infile: str or pathlib.Path, outfile: str or pathlib.Path):
    '''
    '''

    if not isinstance(outfile, pathlib.Path):
        outfile: pathlib.Path = pathlib.Path(outfile)

    with open(infile) as datafile:
        fusion_data = pandas.read_csv(datafile, header=0)

    rows, _ = fusion_data.shape

    for row in range(rows):
        row_of_interest = fusion_data.iloc[row, :].copy()

        hgene, tgene = row_of_interest["hgene"], row_of_interest["tgene"]
        henst, tenst = row_of_interest["henst"], row_of_interest["tenst"]
        hchrm, tchrm = row_of_interest["hchrm"], row_of_interest["tchrm"]
        hstrand, tstrand = row_of_interest["hstrand"], row_of_interest["tstrand"]
        seq, ctype, source = row_of_interest["seq"], row_of_interest["Ctype"], row_of_interest["Source"]

        classification, short_distance, head2tailDistance = row_of_interest["classification"], row_of_interest["shortDistance"], row_of_interest["head2tailDistance"]

        fusion = Fusions(hgene = hgene, tgene = tgene, henst = henst, tenst = tenst, hchrm = hchrm, tchrm = tchrm, hstrand = hstrand, tstrand = tstrand, seq = seq, ctype = ctype, source = source)

        # Fuck I forgot that these things go out as lists but come back in as strings
        convert_list = ["hexonStarts", "texonStarts", "hexonEnds", "texonEnds", "hexonFrames", "texonFrames", "hblockSizes", "tblockSizes", "htStarts", "ttStarts"]
        for strtup in convert_list:
            row_of_interest[strtup] = CM.convert2list(row_of_interest[strtup])

        head_gene = Gene(row_of_interest["hgene"],
                        gname = row_of_interest["hgname"],
                        txStart = row_of_interest["htxStart"],
                        txEnd = row_of_interest["htxEnd"],
                        cdsStart =row_of_interest["hcdsStart"],
                        cdsEnd = row_of_interest["hcdsEnd"],
                        exonCount = row_of_interest["hexonCount"],
                        exonStarts = row_of_interest["hexonStarts"],
                        exonEnds = row_of_interest["hexonEnds"],
                        exonFrames = row_of_interest["hexonFrames"])

        tail_gene = Gene(row_of_interest,
                        gname = row_of_interest["tgname"],
                        txStart = row_of_interest["ttxStart"],
                        txEnd = row_of_interest["ttxEnd"],
                        cdsStart = row_of_interest["tcdsStart"],
                        cdsEnd = row_of_interest["tcdsEnd"],
                        exonCount = row_of_interest["texonCount"],
                        exonStarts = row_of_interest["texonStarts"],
                        exonEnds = row_of_interest["texonEnds"],
                        exonFrames = row_of_interest["texonFrames"])

        head_blat = Blat(qName = row_of_interest["hgene"],
                        tStart = row_of_interest["htStart"],
                        tEnd = row_of_interest["htEnd"],
                        blockCount = row_of_interest["hblockCount"],
                        blockSizes = row_of_interest["hblockSizes"],
                        tStarts = row_of_interest["htStarts"])

        tail_blat = Blat(qName = row_of_interest["tgene"],
                        tStart = row_of_interest["ttStart"],
                        tEnd = row_of_interest["ttEnd"],
                        blockCount = row_of_interest["tblockCount"],
                        blockSizes = row_of_interest["tblockSizes"],
                        tStarts = row_of_interest["ttStarts"])



        fusion.second_import(classification, short_distance, head2tailDistance, head_gene, tail_gene, head_blat, tail_blat)
        print(f"####\n{fusion.hgene}_{fusion.tgene}\n{fusion.henst}_{fusion.tenst}")

        score: DissSimilarityScore = DissSimilarityScore(fusion)
        score.score()
        score.write_score(outfile = outfile)


        # print(vars(fusion))
        # print(fusion.hgene)

        # if row == 0:
        #     exit()


def addCancerTypes():
    '''
    Because now suddenly the cancer types are important!

    This is going to be ad hoc, because what I should do is add things like the cancer and seq in the future databases when I make the databases, rather then add them in at a later date.
    So this isn't going to be a very modular code, it's going to work for what I have in front of me.
    '''

    original_file: pathlib.Path = pathlib.Path.cwd().parent / "Data_Files" / "UTData_cds.csv"
    updated_file: pathlib.Path = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_BE_min100.csv"

    with open(original_file) as og:
        original_data: pandas.DataFrame = pandas.read_csv(og, header = 0)

    if updated_file.suffix in ".xlsx":
        with open(updated_file) as ex:
            updated_data: pandas.DataFrame = pandas.read_excel(ex, sheet_name = "SlipJunctionWSeq", header = 0)
    elif updated_file.suffix in ".csv":
        with open(updated_file) as ex:
            updated_data: pandas.DataFrame = pandas.read_csv(ex, header = 0)

    outpath = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_BE_min100_Ctype.csv"
    # pandas.DataFrame(columns=tuple(slippage_data.columns)).to_csv(outpath, header = True, index = False)

    # print(slippage_data)

    # original_data["CName"] = original_data["Hgene"] + "_" + original_data["Tgene"]
    # original_data["Censt"] = original_data["Henst"] + "_" + original_data["Tenst"]

    # print(original_data)

    rows, _ = updated_data.shape

    score_2oh = pandas.DataFrame()
    for row in range(rows):
        row_of_interest = updated_data.iloc[row, :].copy()
        matching_row = original_data.index[(original_data["Seq"] == row_of_interest["seq"]) & \
                                            (original_data["Hgene"] == row_of_interest["hgene"]) & \
                                            (original_data["Tgene"] == row_of_interest["tgene"]) & \
                                            (original_data["Henst"] == row_of_interest["henst"]) & \
                                            (original_data["Tenst"] == row_of_interest["tenst"]) & \
                                            (original_data["Hstrand"] == row_of_interest["hstrand"]) & \
                                            (original_data["Tstrand"] == row_of_interest["tstrand"]) & \
                                            (original_data["Hchr"] == row_of_interest["hchrm"]) & \
                                            (original_data["Tchr"] == row_of_interest["tchrm"])].to_list()
        # matching_row = original_data[original_data["Seq"] == row_of_interest["Seq"]]
        # match_index = matching_row.i
        # print(original_data.loc[matching_row, "Ctype"].iat[0])
        # print(type(original_data.loc[matching_row, "Ctype"].iat[0]))

        row_of_interest["Ctype"] = original_data.loc[matching_row, "Ctype"].iat[0]
        row_of_interest["Source"] = original_data.loc[matching_row, "Source"].iat[0]
        # print(matching_row)
        # print(row_of_interest)

        score_2oh = pandas.concat([score_2oh, row_of_interest.to_frame().T], axis = 0)

        # if row > 5:
        #     print(slippage_2oh)
        #     slippage_2oh.to_csv(outpath, index = False)
        #     exit()

    score_2oh.to_csv(outpath, index = False)






if __name__ in '__main__':
    # fusion_file = pathlib.Path.cwd().parent / "Data_Files" / "UTData_fs.csv"
    # databaseBuild(fusion_file, start_index=0)
    # refusion_file = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_BE_min100.csv"
    # gapSearch(refusion_file)
    # addCancerTypes()
    scoreing(infile = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_BE_min100_Ctype.csv", outfile = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Scoring_min100_83022_15k.csv")
