from ast import dump
from xmlrpc.client import Server
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

# It's unavoidable: I must use classes. I have to much I'm passing around and I need to pull it out at different times.
# I need a class, I need access to self, I need attributes.

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
    ht_blat_headers = ["HtStart", "HblockCount", "HblockSizes", "HtStarts", "TtStart", "TblockCount", "TblockSizes", "TtStarts"]
    ht_enst_headers = ["HenstURL", 'HexonCount', 'HexonStarts', 'HexonEnds', 'HexonFrames', "TenstURL", 'TexonCount', 'TexonStarts', 'TexonEnds', 'TexonFrames']

    utheadurl = utheaders + ["BlatURL"]
    utheaderrors = utheaders + ["LastURL", "Error", "Type"]

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

            enst_series: pandas.DataFrame = ensting(blat, henst, tenst)
            # enst_index = tuple(enst_series.index)
            # print(f"Blat:\n{blat.T}")
            # print(f"Enst:\n{enst_series.T}")

            if isinstance(enst_series, pandas.Series):
                print("~~~ Clean ENST ~~~")

                row_of_interest: pandas.Series = pandas.concat([row_of_interest, enst_series], axis = 0)
                row_of_interest = row_of_interest[utheadurl + ht_blat_headers + ht_enst_headers]
                # print(row_of_interest.T)
                # print(list(row_of_interest.index))
                # exit()
                # print(list(row_of_interest.index))
                # print(enst_series)

                row_of_interest.to_frame().T.to_csv(out_blat, header = None, index = None, mode = 'a')
                # exit()

            elif isinstance(enst_series, tuple):

                row_of_interest["LastURL"], row_of_interest["Error"], row_of_interest["Type"] = enst_series[0], enst_series[1], enst_series[2]

                row_of_interest.to_frame().T.to_csv(out_error, header = None, index = None, mode = 'a')


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



def enst_attempt(chrom: str, start: int, end: int, error_attempts: int = 5):
    '''
    '''
    try:
        enstURL = api.ens_tracks(chrom=chrom, start=start, end=end)
    
    except Exception as e:
        print("$$$$      New Error     $$$$")
        logger_output(message_title="New Error at ENST", data=f"URL:\nChr: {chrom}\tStart: {start}\tEnd: {end}\nError: {e}\n\tType: {type(e)}\n\nTraceback: {traceback.format_exc()}")

        enstURL = None
    
    return enstURL



def enst_convert_attemp(enstURL, error_attempts: int = 5):
    '''
    '''

    for attempt in range(error_attempts):

        try:
            enst_frame: pandas.DataFrame = api.convert2frame(RQuery.query(enstURL))
            break

        except AttributeError as e:
            print("$$$$ Attribute Error $$$$")
            # print("trying again")
            # enst_frame = e
            print(traceback.format_exc())
            print(enstURL)
            exit()
        
        except Exception as e:
            print("$$$$   New ENST Error   $$$$")

            print(f"Error: {e}\nType: {type(e)}")

            print(traceback.format_exc())

            # logger_output(message_title=f"Error at trying to convert ENST Response to Dataframe", data=f"Fusion: {henst}_{tenst}\nError: {e}\n\tType: {type(e)}\nURL: {enstURL}\n")
            # return_series: tuple = (enstURL, e, type(e))
            # dump_frame = None
            # enst_frame = None
            exit()

    else:
        print("$$$$ Failed to convert to frame $$$$")
        exit()

    return enst_frame




def ensting(oblat: pandas.DataFrame, henst: str, tenst: str) -> pandas.Series or None or tuple:
    '''
    This does quit a bit: it does two things.

    First it pulls up the enst track and identifies which gene we're looking at. BLAT does not do this, it only aligns it to the genome. Looking at the ENST track from the BLAT
    allows us to identify which specific gene we're looking at.

    Next it combines this data with the BLAT data to make a smaller, useful subset. Why do it like this? Because again, BLAT only aligns to the human genome, it does not tell you which
    gene you're actually looking at. The BLAT data needs to be preserved, but it needs to be done via the identify. So at some point I was going to have to combine the ENST and the 
    BLAT data. Might as well do it here.

    '''

    # I'm doing this as a tuple, because a tuple preserves the order. The order is important. I can switch to a set at a later moment, but I'll keep it as a tuple most often
    string2tuple = lambda input_str: tuple(re.split(',', input_str))

    blat = oblat.copy()
    # print(blat.T)
    # exit()

    enst_index = []
    rows, _ = blat.shape
    blat_of_interest = ["tStart", "blockCount", "blockSizes", "tStarts"]
    enst_of_sets = ['exonStarts', 'exonEnds', 'exonFrames']
    enst_of_interest = ["enstURL", 'exonCount'] + enst_of_sets

    head_names, tail_names = [f"H{index}" for index in blat_of_interest + enst_of_interest], [f"T{index}" for index in blat_of_interest + enst_of_interest]

    order_index = head_names + tail_names

    dump_frame = pandas.DataFrame(dtype='str')
    return_series = pandas.Series(dtype='str')
    enst_index, url_index = [], []
    for row in range(rows):
        row_of_interest: pandas.Series = blat.iloc[row, :].squeeze()

        enstURL: str = enst_attempt(chrom=row_of_interest["tName"], start=row_of_interest["tStart"], end=row_of_interest["tStart"] + row_of_interest["blockSizes"][0])

        # print(type(enstURL))
        # exit()

        if isinstance(enstURL, str):
            enst_frame = enst_convert_attemp(enstURL)
        else:
            break
        
        enst_set: set = set(enst_frame["name"])

        enst_intersection = enst_set.intersection([henst, tenst])

        if len(enst_intersection) == 1:
            if henst in enst_intersection:
                enst_index.append(henst), url_index.append(enstURL)
                ienst = henst
                
            elif tenst in enst_intersection:
                enst_index.append(tenst), url_index.append(enstURL)
                ienst = tenst

            enst_frame = enst_frame[enst_frame["name"] == ienst]

            for index in enst_of_sets:
                enst_frame.loc[enst_frame.index[0], index] = string2tuple(enst_frame.loc[enst_frame.index[0], index])
            
            dump_frame = pandas.concat([dump_frame, enst_frame[enst_frame["name"] == ienst]])

        else:
            print(f"$$$$ ENST Issue $$$$\nlength = {len(enst_intersection)}\nPrinting to Log")
            logger_output(message_title=f"No length to intersection at the following URL", data=f"Fusion: {henst}_{tenst}\nURL:\n{enstURL}")

            dump_frame = None
            enst_frame = None

            break

    if isinstance(dump_frame, pandas.DataFrame): # and isinstance(enst_frame, list):

        # for dow in range(dump_frame.shape[0]):
        #     for index in enst_of_sets:
        #         print(dump_frame.loc[dow, index])
        #         dump_frame.loc[dow, index] = string2tuple(dump_frame.loc[dow, index])

        blat = blat[blat_of_interest]
        blat["Ident"] = enst_index
        blat = blat.set_index(blat["Ident"])
        del blat["Ident"]

        dump_frame["enstURL"] = url_index
        dump_frame = dump_frame[enst_of_interest]
        dump_frame["Ident"] = enst_index
        dump_frame = dump_frame.set_index(dump_frame["Ident"])
        del dump_frame["Ident"]

        dump_frame = pandas.concat([blat, dump_frame], axis = 1)

        rows, _ = dump_frame.shape
        for row in range(rows):
            future_series: pandas.Series = dump_frame.iloc[row, :].squeeze()
            new_index = head_names if future_series.name in henst else tail_names

            future_series.index = new_index

            return_series = pandas.concat([return_series, future_series], axis = 0)

        return_series = return_series[order_index]

    else:
        return_series = None

    return return_series
        

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



def main():
    '''
    '''
    blating()
    # ensting()




if __name__ in '__main__':
    main()