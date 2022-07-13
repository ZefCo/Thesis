import pandas
import pathlib



genome = 'hg19'
csv_dir = pathlib.Path.cwd() / "Data_Files" / "Sequence_Files" / f"{genome.upper()}" / "SubFiles"

subfiles = csv_dir.iterdir()

# print(subfiles)
master_frame = pandas.DataFrame()

for subfile in subfiles:
    if str(subfile.suffix) in '.csv':
        # print(subfile)
        with open(subfile, 'r') as csv_subfile:
            subdata = pandas.read_csv(csv_subfile)

    master_frame = pandas.concat([master_frame, subdata], axis = 0)

# print(master_frame["Exon_1"])

master_frame.to_csv(csv_dir.parent / f"Known_Genes_{genome}.csv")