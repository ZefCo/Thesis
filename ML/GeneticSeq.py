import pandas
import pathlib
import re

def main():
    '''
    '''

    cwd = pathlib.Path.cwd()

    train_data: pandas.DataFrame = pandas.read_pickle(cwd.parent / "Data_Files" / "TrainingData" / "TrainingGeneDataSample.pkl")
    train_data.pop("ename"), train_data.pop("gname")

    train_data_cols = tuple(train_data.columns)

    utr3_list, cds_list, intron_list, utr5_list = list(), list(), list(), list()

    for col in train_data_cols:
        if re.search("Cord", col):
            train_data.pop(col)
        elif re.search("Seq", col):
            if re.search("CDS", col):
                cds_list.append(col)
            elif re.search("Intron", col):
                intron_list.append(col)
            elif re.search("UTR3", col):
                utr3_list.append(col)
            elif re.search("UTR5", col):
                utr5_list.append(col)

    cds_list.sort(key = len)
    intron_list.sort(key = len)
    utr5_list.sort(key = len)
    utr3_list.sort(key = len)

    full_seq = f""

    

    # print(train_data)



if __name__ in "__main__":
    main()