# Because I should have done this to begin with

import pandas
import pathlib

cwd = pathlib.Path.cwd()
root = cwd.parent

nucFiles = root / "Data_Files" / "NucComp"

csv_files = nucFiles.rglob("*.csv")

for j, file in enumerate(csv_files):
    print(file)
    with open(file) as csv:
        # pass
        try:
            local_frame = pandas.read_csv(csv)
        except Exception as e:
            local_frame = None
            print(e)
            print(type(e))

    if local_frame is not None:

        column_names = list(local_frame.columns)

        new_cols = {"Column1": "Name",	"Column2": "GName",	"Column3": "EName",	"Column4": "NName",	"Column5": "Chr",	"Column6": "Strand"}
        new_col_range = int((len(column_names) - 6) / 2)

        new_colsS = {f"Column{(2*i + 1) + 6}": f"Seq_{i}" for i in range(new_col_range)}
        new_colsC = {f"Column{(2*i) + 8}": f"Cord_{i}" for i in range(new_col_range)}

        # print(new_colsS)

        local_frame.rename(columns = new_cols, inplace = True)
        local_frame.rename(columns = new_colsS, inplace = True)
        local_frame.rename(columns = new_colsC, inplace = True)

        # print(local_frame)
        # print(local_frame.columns)

        rows, cols = local_frame.shape
    
        for row in range(rows):
            row_of_interest = local_frame.iloc[row, :].copy()

            if row_of_interest["Strand"] in "-":
                for col in new_colsS.values():
                    if isinstance(row_of_interest[col], str):
                        seq = row_of_interest[col]
                        seq = seq[::-1]
                        row_of_interest[col] = seq

            local_frame.loc[row, :] = row_of_interest

        # new_file = file.basename

        file_name = file.stem
        file_name = f"{file_name}_rev.csv"
        # print(file_name)

        local_frame.to_csv(nucFiles / file_name, index = False, header = True)

