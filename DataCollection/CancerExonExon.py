import re
import ucsc_restapi as upi
import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import pickle


def main():
    '''
    Not sure how best to tackle and save this information.

    Import the data

    Take only those with high confidence

    With the sequence: it consists of NNN___NNNN|NNNNN
    Take only 6+ and 6- from the |
    '''
    cancer_file = cwd.parent / "Data_Files" / "FullfusionlistXS_Ship1_2_3_4_5_6_02272024_ShipmentInfoxlsx.xlsx"


    cancer_data = pandas.read_excel(cancer_file, sheet_name = "Sheet 1")
    cancer_data["CANCER TYPE"] = cancer_data["CANCER TYPE"].astype("category")
    cancer_data["confidence"] = cancer_data["confidence"].astype("category")

    # rows, cols = cancer_data.shape
    # bad = 0
    # for row in range(rows):
    #     fusion_transcript = cancer_data.loc[row, "fusion_transcript"]
    #     if re.search(r"\|", fusion_transcript):
    #         pass
    #     else:
    #         bad += 1

    # print(bad)

    cancer_data["Tail"] = cancer_data["fusion_transcript"].apply(lambda x: re.split(r"\|", x)[1] if re.search(r"\|", x) else pandas.NA)
    cancer_data["Head"] = cancer_data["fusion_transcript"].apply(lambda x: re.split(r"\|", x)[0] if re.search(r"\|", x) else pandas.NA)

    cancer_data = cancer_data.dropna()
    # print(cancer_data)

    cancer_pickle = cwd.parent / "Data_Files" / "Fusion_AllConfidence.pkl"
    cancer_data.to_pickle(cancer_pickle)

    confidence = tuple(cancer_data["confidence"].unique())
    for con in confidence:
        con_data: pandas.DataFrame = cancer_data[cancer_data["confidence"] == con]
        con_data = con_data.reset_index()
        conf_pickle = cwd.parent / "Data_Files" / f"Fusion_{con}.pkl"
        con_data.to_pickle(conf_pickle)




if __name__ in "__main__":
    main()