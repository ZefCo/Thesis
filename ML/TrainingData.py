import pathlib
import pandas
import re


def import_file(file_path: str or pathlib.Path, seq_name: str) -> pandas.DataFrame:
    '''
    '''

    return_data = pandas.read_csv(str(file_path), header = 0)

    column_names = tuple(return_data.columns)
    new_col_names = dict()

    for name in column_names:
        if re.match("Cord", name):
            return_data.pop(name)
        elif re.match("Seq", name):
            new_col_names[name] = f"{seq_name}_{name}"
        # else:
        #     new_col_names[name] = name

    return_data.pop("GName")
    return_data.pop("EName")

    return_data = return_data.rename(columns=new_col_names)

    return return_data


def main():
    '''
    '''

    cwd = pathlib.Path.cwd()
    train_files = cwd.parent / "Data_Files" / "TrainingData"

    utr5_data = import_file(str(train_files / "UTR5_seq.csv"), "UTR5")
    utr3_data = import_file(str(train_files / "UTR3_seq.csv"), "UTR3")
    cds_data = import_file(str(train_files / "CDS_seq.csv"), "CDS")
    intron_data = import_file(str(train_files / "Intron_seq.csv"), "INT")

    # ran_sample = utr5_data.sample(n = 1000)

    # testing_data = pandas.DataFrame()

    rows, _ = utr5_data.shape

    for row in range(rows):
        utr5_sample = utr5_data.iloc[row, :]
        name, nname, strand, chr = utr5_sample["Name"], utr5_sample["NName"], utr5_sample["Strand"], utr5_sample["Chr"]
        print(utr5_sample)


        cds_sample = cds_data[(cds_data["Name"] == name) & (cds_data["NName"] == nname) & (cds_data["Strand"] == strand) & (cds_data["Chr"] == chr)].squeeze()
        if cds_sample.empty:
                cds_sample["Name"] = name
                cds_sample["NName"] = nname
                cds_sample["Strand"] = strand
                cds_sample["Chr"] = chr
        print(f"\n")
        print(cds_sample)

        intron_sample = intron_data[(intron_data["Name"] == name) & (intron_data["NName"] == nname) & (intron_data["Strand"] == strand) & (intron_data["Chr"] == chr)].squeeze()
        if intron_sample.empty:
                intron_sample["Name"] = name
                intron_sample["NName"] = nname
                intron_sample["Strand"] = strand
                intron_sample["Chr"] = chr
        print(f"\n")
        print(intron_sample)

        utr3_sample = utr3_data[(utr3_data["Name"] == name) & (utr3_data["NName"] == nname) & (utr3_data["Strand"] == strand) & (utr3_data["Chr"] == chr)].squeeze()
        if utr3_sample.empty:
                utr3_sample["Name"] = name
                utr3_sample["NName"] = nname
                utr3_sample["Strand"] = strand
                utr3_sample["Chr"] = chr
        print(f"\n")
        print(utr3_sample)
        print(f"\n")



        if row == 3:
            break
        #     print(f"\n")
        #     print(utr5_sample)
        #     print(f"\n")
        #     print(cds_sample)
        #     print(f"\n")
        #     print(intron_sample)
        #     print(f"\n")
        #     print(utr3_sample)
        #     print(f"\n")



        # row_to_merge = pandas.concat([utr5_sample, cds_sample, intron_sample, utr3_sample])

        # print(row_to_merge)
        # print(type(row_to_merge))

        # # utr5_file, index = False, header = False, mode = 'a'

        # if row == 10:
        #     break



        # try:
        #     testing_data = pandas.concat([testing_data, row_to_merge.to_frame().T])
        # except Exception as e:
        #     print(type(e))
        #     print(e)
        #     print(utr5_sample)
        #     print(cds_sample)
        #     print(intron_sample)
        #     print(utr3_sample)
        #     break

    # testing_data.to_csv(str(train_files / "Training_Data.csv"), index = False, header = True)


def sanity():
     '''
     '''





if __name__ in "__main__":
    main()