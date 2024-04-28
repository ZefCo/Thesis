import ExonExonData as eed
import pathlib
import re
cwd = pathlib.Path.cwd()
import pandas
import pickle
from Heatmaps import heatmapv2 as heatmap
from Heatmaps import _undigitize_seq as undigit
from TimeEmbedding import time_embedding as te
from ExonExonData import gen_plot


def main():
    '''
    '''
    all_levels_file = cwd.parent / "Data_Files"/ "UTData_fs.pkl"
    all_levels_plot = cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "UT_Data.png"
    create_plot(all_levels_file, all_levels_plot, title = "Exon to Exon Comparison\nUT Database")

    # Lab data from the University
    # print("All Confidence Levels")
    # all_levels_file = cwd.parent / "Data_Files"/ "Fusion_AllConfidence.pkl"
    # all_levels_plot = cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "AllConfidence.png"
    # create_plot(all_levels_file, all_levels_plot)
    
    # print("High Confidence Levels")
    # high_levels_file = cwd.parent / "Data_Files"/ "Fusion_high.pkl"
    # high_levels_plot = cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "HighConfidence.png"
    # create_plot(high_levels_file, high_levels_plot, title = f"Cancer: High Confidence")
    
    # print("Low Confidence Levels")
    # low_levels_file = cwd.parent / "Data_Files"/ "Fusion_low.pkl"
    # low_levels_plot = cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "LowConfidence.png"
    # create_plot(low_levels_file, low_levels_plot, title = f"Cancer: Low Confidence")
    
    # print("Medium Confidence Levels")
    # medium_levels_file = cwd.parent / "Data_Files"/ "Fusion_medium.pkl"
    # medium_levels_plot = cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "MediumConfidence.png"
    # create_plot(medium_levels_file, medium_levels_plot, title = f"Cancer: Medium Confidence")



def clean_frame(frame: pandas.DataFrame):
    '''
    Cleans the frame, in case something like a float slips in
    '''
    frame = frame.dropna()
    frame = frame.reset_index()

    return frame


def create_plot(cancer_data: pathlib.Path or pandas.DataFrame, plot_file: pathlib.Path, title: str = None, *args, **kwargs):
    '''
    '''

    if isinstance(cancer_data, pathlib.Path):
        cancer_data: pandas.DataFrame = pandas.read_pickle(cancer_data)

    cancer_data = clean_frame(cancer_data)

    cancer_data["HL"] = cancer_data["Head"].apply(lambda x: len(x))
    cancer_data["TL"] = cancer_data["Tail"].apply(lambda x: len(x))
    cancer_data = cancer_data[(cancer_data["HL"] > 5) & (cancer_data["TL"] > 5)]
    cancer_data = cancer_data.reset_index()
    cancer_data = cancer_data.drop("index", axis = 1)

    cancer_data["Head_Seq"] = cancer_data["Head"].apply(lambda x: x[len(x) - 6::])
    cancer_data["Tail_Seq"] = cancer_data["Tail"].apply(lambda x: x[0:6])

    cancer_data["Fusion_Seq"] = cancer_data["Head_Seq"] + cancer_data["Tail_Seq"]
    cancer_data["Fusion_Seq"] = cancer_data["Fusion_Seq"].apply(lambda x: pandas.NA if re.search(r"\_|\.", x) else x)
    cancer_data = cancer_data.dropna(how = "any")

    try:
        cancer_data = cancer_data.reset_index()
        cancer_data = cancer_data.drop("index", axis = 1)
    except Exception as e:
        print("First error")
        print(type(e))

    cancer_data["XY"] = cancer_data["Fusion_Seq"].apply(lambda x: te(x)[0])
    cancer_data["X"] = cancer_data["XY"].apply(lambda x: x[1]) # remember, things got switched around. I forgot this the first time I ran it through
    cancer_data["Y"] = cancer_data["XY"].apply(lambda x: x[0])

    cancer_data["X"] = cancer_data["X"].apply(lambda x: pandas.NA if x > 1 else x)
    cancer_data["Y"] = cancer_data["Y"].apply(lambda x: pandas.NA if x > 1 else x)
    cancer_data = cancer_data.dropna(how = "any")
    try:
        cancer_data = cancer_data.reset_index()
        cancer_data = cancer_data.drop("index", axis = 1)
    except Exception as e:
        print("Second Error")
        print(type(e))

    # print_data = cancer_data[["Fusion_Seq", "X", "Y"]]
    # print(print_data.loc[16, "Fusion_Seq"])

    # # print_data.to_csv("Sanity.csv")
    # # print(cancer_data)
    # exit()

    # # print(cancer_data)
    # print(cancer_data.loc[0, "Head_Seq"], cancer_data.loc[0, "Tail_Seq"], cancer_data.loc[0, "Fusion_Seq"])
    # print(cancer_data.loc[0, "X"], cancer_data.loc[0, "Y"])

    if title is None:
        title = "Exon-Exon: Cancer"

    gen_plot(cancer_data, file_name = plot_file, title = title, inches = 10)

if __name__ in "__main__":
    main()