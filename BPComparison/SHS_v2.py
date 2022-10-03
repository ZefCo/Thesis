import enum
from importlib.resources import path
from re import T
import pandas
from pyparsing import col
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


length = 10
random_sample = 100

with open(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_BE_min100_Ctype.csv") as csvfile:
    collected_data = pandas.read_csv(csvfile, header = 0)

unique_data = collected_data.groupby(["tgene", "hgene"]).first()
# print(unique_data["hgene"])
unique_data = unique_data.reset_index()
rample_data = unique_data.sample(n = random_sample)
rows, cols = rample_data.shape

for row in range(rows):
    row_of_interest = rample_data.iloc[row, :].copy()
    fusion = Fusions(row_of_interest["hgene"], row_of_interest["tgene"], row_of_interest["seq"], row_of_interest["henst"], row_of_interest["tenst"], row_of_interest["hstrand"], row_of_interest["tstrand"], row_of_interest["hchrm"], row_of_interest["tchrm"], row_of_interest["ctype"], row_of_interest["source"])

    classification, short_distance, head2tailDistance = row_of_interest["classification"], row_of_interest["shortDistance"], row_of_interest["head2tailDistance"]

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

    print(f"Head Seq: {fusion.hgene}      {fusion.henst}      {fusion.hstrand}      {fusion.hchrm}")
    print(f"Tail Seq: {fusion.tgene}      {fusion.tenst}      {fusion.tstrand}      {fusion.tchrm}")

    # print(fusion.head_blat.tStarts)
    hindex, tindex = fusion._align_blat()

    try:
        if (((fusion.hstrand in "+") and ((hindex + 1) <= fusion.head_gene.exonCount)) or ((fusion.hstrand in "-") and ((hindex - 1) >= 0))):
            if fusion.hstrand in "+":
                hSeqStart = fusion.head_gene.exonEnds[hindex] - length
                hSeqEnd = fusion.head_gene.exonEnds[hindex]

                hIntStart = fusion.head_gene.exonEnds[hindex]
                hIntEnd = fusion.head_gene.exonEnds[hindex] + length

                hExoStart = fusion.head_gene.exonStarts[hindex + 1]
                hExoEnd = fusion.head_gene.exonStarts[hindex + 1] + length

                hSeq_rev_dir: bool = True
                hInt_rev_dir: bool = False
                hExo_rev_dir: bool = False

            elif fusion.hstrand in "-":
                hSeqStart = fusion.head_gene.exonStarts[hindex]
                hSeqEnd = fusion.head_gene.exonStarts[hindex] + length

                hIntStart = fusion.head_gene.exonStarts[hindex] - length
                hIntEnd = fusion.head_gene.exonStarts[hindex]

                hExoStart = fusion.head_gene.exonEnds[hindex - 1] - length
                hExoEnd = fusion.head_gene.exonEnds[hindex - 1]

                hSeq_rev_dir: bool = False
                hInt_rev_dir: bool = True
                hExo_rev_dir: bool = True

        H3_Seq, H3_Seq_url = upi.sequence(chrom = fusion.hchrm, start = hSeqStart, end = hSeqEnd, strand = fusion.hstrand, reverse_dir = hSeq_rev_dir)
        H5_Int, H5_Int_url = upi.sequence(chrom = fusion.hchrm, start = hIntStart, end = hIntEnd, strand = fusion.hstrand, reverse_dir = hInt_rev_dir)
        H5_Exo, H5_Exo_url = upi.sequence(chrom = fusion.hchrm, start = hExoStart, end = hExoEnd, strand = fusion.hstrand, reverse_dir = hExo_rev_dir)

    except Exception as e:
        print(f"Error happened at Head, so I'm skipping it\nError type {type(e)}")
        H3_Seq, H5_Int, H5_Exo = "", "", ""
        for _ in range(length):
            H3_Seq, H5_Int, H5_Exo = f"{H3_Seq}X", f"{H5_Int}X", f"{H5_Exo}X"
        # break


            # hExoStart: int = fusion.head_gene.exonStarts[head5primeIndex]
            # hExoEnd: int = fusion.head_gene.exonEnds[head5primeIndex]

    try:
        if (((fusion.tstrand in "+") and ((tindex - 1) >= 0)) or ((fusion.tstrand in "-") and ((tindex + 1) <= fusion.tail_gene.exonCount))):
            if fusion.tstrand in "+":
                tSeqStart = fusion.tail_gene.exonStarts[tindex]
                tSeqEnd = fusion.tail_gene.exonStarts[tindex] + length

                tIntStart = fusion.tail_gene.exonStarts[tindex] - length
                tIntEnd = fusion.tail_gene.exonStarts[tindex]

                tExoStart = fusion.tail_gene.exonEnds[tindex - 1] - length
                tExoEnd = fusion.tail_gene.exonEnds[tindex - 1]

                tSeq_rev_dir: bool = False
                tInt_rev_dir: bool = True
                tExo_rev_dir: bool = True

            elif fusion.tstrand in "-":
                tSeqStart = fusion.tail_gene.exonEnds[tindex] - length
                tSeqEnd = fusion.tail_gene.exonEnds[tindex]

                tIntStart = fusion.tail_gene.exonEnds[tindex]
                tIntEnd = fusion.tail_gene.exonEnds[tindex] + length

                tExoStart = fusion.tail_gene.exonStarts[tindex + 1]
                tExoEnd = fusion.tail_gene.exonStarts[tindex + 1] + length

                tSeq_rev_dir: bool = True
                tInt_rev_dir: bool = False
                tExo_rev_dir: bool = False

        T5_Seq, T5_Seq_url = upi.sequence(chrom = fusion.tchrm, start = tSeqStart, end = tSeqEnd, strand = fusion.tstrand, reverse_dir = tSeq_rev_dir)
        T3_Int, T3_Int_url = upi.sequence(chrom = fusion.tchrm, start = tIntStart, end = tIntEnd, strand = fusion.tstrand, reverse_dir = tInt_rev_dir)
        T3_Exo, T3_Exo_url = upi.sequence(chrom = fusion.tchrm, start = tExoStart, end = tExoEnd, strand = fusion.tstrand, reverse_dir = tExo_rev_dir)

    except Exception as e:
        print(f"Error happened at Tail, so I'm skipping it\nError type {type(e)}")
        T5_Seq, T3_Int, T3_Exo = "", "", ""
        for _ in range(length):
            T5_Seq, T3_Int, T3_Exo = f"{T5_Seq}X", f"{T3_Int}X", f"{T3_Exo}X"


    # print(f"Head Seq: {fusion.hgene}      {fusion.henst}      {fusion.tstrand}")
    # print("Head Exon ~ Intron ~ Exon + 1")
    # print(hSeqStart, hSeqEnd)
    # print(hIntStart, hIntEnd)
    # print(hExoStart, hExoEnd)
    # print(f"Tail Seq: {fusion.tgene}      {fusion.tenst}      {fusion.tstrand}")
    # print("Tail Exon ~ Intron ~ Exon - 1")
    # print(tSeqStart, tSeqEnd)
    # print(tIntStart, tIntEnd)
    # print(tExoStart, tExoEnd)



    # print(f"H3_Seq_url:\n{H3_Seq_url}")
    # print(f"H5_Int_url:\n{H5_Int_url}")
    # print(f"H5_Exo_url:\n{H5_Exo_url}")
    # print(f"H3_Seq Seq:   {H3_Seq}")
    # print(f"H5_Int Seq:   {H5_Int}")
    # print(f"H5_Exo Seq:   {H5_Exo}")

    # print(f"T5_Seq_url:\n{T5_Seq_url}")
    # print(f"T3_Int_url:\n{T3_Int_url}")
    # print(f"T3_Exo_url:\n{T3_Exo_url}")
    # print(f"T5_Seq Seq:   {T5_Seq}")
    # print(f"T3_Int Seq:   {T3_Int}")
    # print(f"T3_Exo Seq:   {T3_Exo}")

    score = DissSimilarityScore(fusion)

    score.sequences = {"H3_Seq": H3_Seq,
                    "H5_Intron": H5_Int,
                    "H5_Exon": H5_Exo,
                    "T5_Seq": T5_Seq,
                    "T3_Intron": T3_Int,
                    "T3_Exon": T3_Exo}

    score.score()
    score.write_score(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Exon2ExonJunctionScores_v2.csv")

    print(f"\n~~Finished Row {row} of {rows}~~\n####\n")


    # print(score.scores)

