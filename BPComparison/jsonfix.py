import pandas
import pathlib
import blat_api as ba

with open(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_Errors_8322.csv") as csvfile:
    utdata = pandas.read_csv(csvfile, header = 0)

second_error_file = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_Error_8422.csv"
if not second_error_file.is_file():
    pandas.DataFrame(columns=["Hgene", "Henst", "Hchr", "Hstrand", "Tgene", "Tenst", "Tchr", "Tstrand", "Seq", "URL", "Error", "Type"]).to_csv(second_error_file, index = None, header=True)

utdata = utdata[utdata["Type"] == "<class 'json.decoder.JSONDecodeError'>"]

start_row = 3341
rows, _ = utdata.shape

for row in range(start_row, rows):
    row_of_interest = utdata.iloc[row, :]

    url = row_of_interest["URL"]
    print(f"Re-trying row {row}")
    try:
        blat: pandas.DataFrame = ba.blat_query(qurl=url)
    except Exception as e:
        print(type(e))
        print(e)
        blat = None

        row_of_interest["Error"] = e
        row_of_interest["Type"] = type(e)

        row_of_interest.to_frame().T.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_Error_8422.csv", header = None, index = None, mode = 'a')

    if isinstance(blat, pandas.DataFrame):
        row_of_interest = row_of_interest[["Hgene", "Henst", "Hchr", "Hstrand", "Tgene", "Tenst", "Tchr", "Tstrand", "Seq", "URL"]]
        row_of_interest.to_frame().T.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_Blat_1000_8322.csv", header = None, index = None, mode = 'a')

 