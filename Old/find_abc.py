import pandas
import pathlib

with open(pathlib.Path.cwd() / "Data_Files" / "UTData_cds.csv") as csvfile:
    ut_cds = pandas.read_csv(csvfile, header = 0)

ut_cds = ut_cds.drop_duplicates(["Hgene", "Tgene"])

head_cds = ut_cds[ut_cds.duplicated('Tgene')]
print(head_cds.shape)
head_cds.to_csv(pathlib.Path.cwd() / "Data_Files" / "Tail_Duplicates.csv", index = True)