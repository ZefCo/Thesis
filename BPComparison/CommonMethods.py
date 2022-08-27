import pathlib
import pandas
import re

def main():
    '''
    '''



def convert2list(sequence: str, ) -> tuple:
    '''
    Because on some of these responses come back as lists, but their not type(list), their type(str). This converts
    them to a list. But it converts EVERYTHING in that column to a list for simplicity, so an entry of 1 is not in a
    list of 1. Trust me, in the long run this works out for the best.

    This was originally in the Blat API, but honestly I need to use it in several places so I moved it here. Easier access
    for everyone.

    OK technically this is a tuple, but list, tuple, set, tomato, potatoe
    
    Yes I know what the actually phrase is!
    '''
    # print(sequence, type(sequence))
    sequence = re.sub(r"\(|\)", "", sequence)
    # print(sequence, type(sequence))
    # exit()

    sequence = re.split(',', sequence)

    seqlist = []
    for strint in sequence:
        
        try:
            strint = int(strint)

        except ValueError as e:
            # print("!!! Value Error !!!")
            # print(f"Error: {e}\tType: {type(e)}\nInput: {strint}")
            strint = None

        except Exception as e:
            # print("!!! Other Error !!!")
            # print(f"Error: {e}\tType: {type(e)}\nInput: {strint}")
            strint = None

        seqlist.append(strint)

    seqlist = tuple(seqlist)

    if seqlist[-1] is None:
        seqlist = list(seqlist)
        seqlist = seqlist[:-1]
        seqlist = tuple(seqlist)

    # print(seqlist)

    return seqlist



def import_csv(filepath: str or pathlib.Path, header = None, index_col = None, col_sets: list or tuple = None) -> pandas.DataFrame:
    '''
    '''
    type_check = lambda variable, intended_type: variable if isinstance(variable, intended_type) else None

    header = type_check(header, int)
    index_col = type_check(index_col, str)

    with open(filepath) as csvfile:
        csvdata: pandas.DataFrame = pandas.read_csv(csvfile, header = header, index_col = index_col)

    if isinstance(col_sets, (list, tuple)):
        rows, _ = csvdata.shape
        print(f"Converting Cols to Tuples\nThis may take a while depending on number of rows: {rows}")
        for row in range(rows):
            # csvdata[col_sets] = convert2list()
            row_of_interest = csvdata.iloc[row, :].copy().squeeze()
            for index_sex in col_sets:
                row_of_interest[index_sex] = convert2list(row_of_interest[index_sex])

            csvdata.iloc[row, :] = row_of_interest

    return csvdata


def compliment_strand(sequence: str) -> str:
    '''
    '''

    conversion_dictionary = {"A": "T",
                             "C": "G",
                             "G": "C",
                             "T": "A"}

    return_seq = f""

    for n in sequence:
        if n in conversion_dictionary:
            return_seq = f"{return_seq}{conversion_dictionary[n]}"

    return return_seq



if __name__ in '__main__':
    main()