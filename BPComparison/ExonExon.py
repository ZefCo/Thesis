from requests import Response
import blat_api as ba
import ucsc_restapi as api
import pathlib
import pandas
import RQuery
import json
import re

import logging
import traceback
from inspect import currentframe, getframeinfo

# https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/

# Establish logging: because it's better then print
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(name)s:%(message)s')

file_handler = logging.FileHandler(f'Exon2ExonComparison.log')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)

error_attempts = 5

def pregen(filename: str, headers: list):
    '''
    I'm tired to typing the same thing over and over again
    '''

    out_file = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / filename

    if not out_file.is_file():
        pandas.DataFrame(columns = headers).to_csv(out_file, header = True, index = False)

    return out_file



def blating():
    '''
    We're going to take this a littel slower. Instead of just doing everything at once:

    Fuck: I forgot that you can't do this one step at a time. I have 40,000 fusions, which means I potentially have 40,000 blats, and blats do NOT come back clean. They could line up
    to several places. In theory I'll have a head and tail blat, but I can have a blat that matches to the same chr multiple times.

    Let's dump the multiple blats to another file. I'll sort through them later.
    
    Idea on how to handle to much blat: do it head/tail. I know I didn't want to do it like this but write another method to break the frame up by head and tail, then
    look at the head version then the tail version.

    '''
    # Where to start in the UT Database, so I can stop the code and start it right where I left off
    # start_index = 5552
    start_index = 0
    # Min length for the sequences of interest
    min_length = 1000
    # Collection of headers for output files
    utheaders = ["Hgene", "Henst", "Hchr", "Hstrand", "Tgene", "Tenst", "Tchr", "Tstrand", "Seq"]
    ht_blat_headers = ["HblockCount", "HblockSizes", "HtStarts", "TblockCount", "TblockSizes", "TtStarts"]
    ht_enst_headers = ["HenstURL" + 'HexonCount', 'HexonStarts', 'HexonEnds', 'HexonFrames', "TenseURL", 'TexonCount', 'TexonStarts', 'TexonEnds', 'TexonFrames']

    utheadurl = utheaders + ["BlatURL"]
    utheaderrors = utheadurl + ["Error", "Type"]

    # Generate output files. Append to these later. It's easier to append if something goes wrong, that way data lost isn't lost forever
    out_zlat = pregen(filename = "NoBlat.csv", headers = utheadurl)
    out_error = pregen(filename = "UTErrors.csv", headers = utheaderrors)
    out_blat = pregen(filename = f"UTBlat_{min_length}.csv", headers = utheadurl + ht_blat_headers + ht_enst_headers)
    out_three = pregen(filename = "ThreePlusBlat.csv", headers = utheadurl)
    
    with open(pathlib.Path.cwd().parent / "Data_Files" / "UTData_cds.csv") as csvfile:
        utdata = pandas.read_csv(csvfile, header = 0)

    utdata = utdata[utdata["SeqLen"] >= min_length]
    utdata = utdata[utheaders]

    rows, _ = utdata.shape

    for row in range(start_index, rows):
        row_of_interest = utdata.iloc[row, :]

        sequence: str = row_of_interest["Seq"]
        qurl = ba.gen_url(qseq = sequence)

        row_of_interest: pandas.Series = row_of_interest.squeeze()
        row_of_interest["BlatURL"] = qurl
        hgene, henst, hstrand, hchr, tgene, tenst, tstrand, tchr = row_of_interest["Hgene"], row_of_interest["Henst"], row_of_interest["Hstrand"], row_of_interest["Hchr"], row_of_interest["Tgene"], row_of_interest["Tenst"], row_of_interest["Tstrand"], row_of_interest["Tchr"]

        print(f"\n####\n{hgene}_{tgene}\n{henst}_{tenst}\n####")

        blat = blat_attempt(qurl = qurl, hgene=hgene, henst=henst, tgene=tgene, tenst=tenst)

        if not isinstance(blat, pandas.DataFrame):
            row_of_interest["Error"] = blat
            row_of_interest["Type"] = type(blat)

            row_of_interest.to_frame().T.to_csv(out_error, header = None, index = None, mode = 'a')

            continue

        # This part of the code is only reaches if there were no blat errors
        # Find the blat we actually want: because blat doesn't include any Gene names or ENST we need to isolate this to only some genes
        blat = blat[((blat["strand"] == hstrand) & (blat["tName"] == hchr)) | ((blat["strand"] == tstrand) & (blat["tName"] == tchr))]

        if blat.shape[0] < 2:
            print("$$$ Not enough blat $$$")
            row_of_interest.to_frame().T.to_csv(out_zlat, header = None, index = None, mode = 'a')

        elif blat.shape[0] > 2:
            print("$$$ To much blat $$$")
            row_of_interest.to_frame().T.to_csv(out_three, header = None, index = None, mode = 'a')

        elif blat.shape[0] == 2:
            print("~~~ Clean Blat ~~~")

            enst_series = pandas.Series = ensting(blat, henst, tenst)
            enst_index = tuple(enst_series.index)
            # print(enst_series)


            # row_of_interest: pandas.Series = pandas.concat([row_of_interest, enst_series])
            # print(row_of_interest)

            # row_of_interest.to_frame().T.to_csv(out_blat, header = None, index = None, mode = 'a')
            exit()


        print(f"Finished UT Database row {row}")



def blat_attempt(qurl: str, hgene: str = None, henst: str = None, tgene: str = None, tenst: str = None, error_attempts: int = 5) -> pandas.DataFrame or None:
    '''
    this block attempts to blat the sequence. If for some reason it doesn't blat it sends the url to an error file
    '''
    for attempt in range(error_attempts):
        try:
            blat: pandas.DataFrame = ba.blat_query(qurl = qurl)
            if attempt > 0:
                logger_output(message_title = "Errors at blat", data = f"{attempt} number of attempts to blat the following CmRNA:\n{hgene}_{tgene}\n{henst}_{tenst}")
            break

        except json.decoder.JSONDecodeError as e:
            print("$$$$ JSON Decoder Error $$$$")
            print("          retrying          ")
            blat = e

        except Exception as e:
            print("$$$$      New Error     $$$$")
            logger_output(message_title="New Error at BLAT", data=f"Fusion: {hgene}_{tgene}\t{henst}_{tenst}\nError: {e}\n\tType: {type(e)}\n\nTraceback: {traceback.format_exc()}")

            blat = e
            
            break
    
    else:
        print("$$$ Failed to blat, sending to log $$$")
        logger_output(message_title=f"Failed to Blat of the following CmRNA after {attempt + 1} number of attempts", data=f"Fusion: {hgene}_{tgene}\t{henst}_{tenst}\nError: {blat}\n\tType: {type(blat)}")

    return blat



def ensting(blat: pandas.DataFrame, henst: str, tenst: str) -> pandas.Series or None:
    '''
    This is going to return the ENST frame

    This needs to identify which is which and change the INDEX to be the ENST, then add to the dataframe which is which based on the Index
    '''

    # I'm doing this as a tuple, because a tuple preserves the order. The order is important. I can switch to a set at a later moment, but I'll keep it as a tuple most often
    string2tuple = lambda input_str: tuple(re.split(',', input_str))

    flat = blat.unstack()
    # print(flat)
    # exit()

    # print(blat.T)
    enst_index = []
    comparison_start, comparison_end = [], []
    rows, _ = blat.shape
    index_of_sets = ['exonStarts', 'exonEnds', 'exonFrames']
    index_of_interest = ['exonCount'] + index_of_sets

    head_names, tail_names = [f"H{index}" for index in index_of_interest], [f"T{index}" for index in index_of_interest]

    order_index = ["HenstURL"] + head_names + ["TenstURL"] + tail_names

# , dtypes=dict(zip(index_of_interest, ["str", "set", "set", "set"]))
    return_frame = pandas.DataFrame(dtype='str')
    enst_index = []
    for row in range(rows):
        row_of_interest: pandas.Series = blat.iloc[row, :].squeeze()

        enstURL = api.ens_tracks(chrom=row_of_interest["tName"], start=row_of_interest["tStart"], end=row_of_interest["tStart"] + row_of_interest["blockSizes"][0])

        enst_frame: pandas.DataFrame = api.convert2frame(RQuery.query(enstURL))
        enst_set: set = set(enst_frame["name"])

        enst_intersection = enst_set.intersection([henst, tenst])

        if len(enst_intersection) == 1:
            if henst in enst_intersection:
                enst_index.append(henst)
                ienst = henst
                # iindex = dict(zip(index_of_interest, head_names))
                # iurl = "HenstURL"
                
            elif tenst in enst_intersection:
                enst_index.append(tenst)
                # return_frame.iloc[row, :] = enst_frame[enst_frame["name"] == tenst]
                ienst = tenst
                # iindex = dict(zip(index_of_interest, tail_names))
                # iurl = "TenstURL"

            # print(enst_frame[enst_frame["name"] == ienst])
            return_frame = pandas.concat([return_frame, enst_frame[enst_frame["name"] == ienst]])
            

            # enst_row: pandas.Series = enst_frame[enst_frame["name"] == ienst].squeeze()
            # # print(enst_row)
            # enst_row = enst_row[index_of_interest]
            
            # for iset in index_of_sets:
            #     # print(f"{iset}\n{enst_row[iset]}")
            #     enst_row[iset] = string2tuple(enst_row[iset])
            #     # print(enst_row[iset])

            # enst_row = enst_row.rename(iindex)
            # enst_row[iurl] = enstURL

            # return_frame = pandas.concat([return_frame, enst_row])

        else:
            print(f"I dont know what happened... length = {len(enst_intersection)}")

            return_frame = None

            break

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # OK I've lost track of what I'm trying to do... I think I'm trying to take the ENST frame and the Blat frame and merge them. For reference use this URL
    # to see what the ENST frame looks like.
    # https://api.genome.ucsc.edu/getData/track?track=ensGene;genome=hg19;chrom=chr1;start=94461663;end=94461751
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    
    blat["Ident"] = enst_index
    blat = blat.set_index(blat["Ident"])

    return_frame = return_frame[index_of_interest]

    print(return_frame.T)
    exit()

    # enst_frame = enst_frame[enst_frame["name"].isin(enst_index)]
    # print(enst_frame)

    # print(blat)
    # exit()

    # if isinstance(return_frame, pandas.Series):
    #     return_frame = return_frame[order_index]


    return return_frame
        

def logger_output(message_title=None, data=None):
    '''
    For outputing log messages
    '''
    if (message_title is None) & (data is not None):
        logger.info(f'\n{data}')
    elif (message_title is not None) & (data is None):
        logger.info(f'\n{message_title}')
    elif (message_title is not None) & (data is not None):
        logger.info(f'\n{message_title}\n{data}')
    else:
        logger.info(f'Something happened that was log worthy, but no message')



# def string_2_set(string_input: str) -> set:
#     '''
#     '''
#     out_set = re.split(',', string_input)




def main():
    '''
    '''
    blating()
    # ensting()




if __name__ in '__main__':
    main()