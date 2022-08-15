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
    start_index = 3530
    # start_index = 0
    # Min length for the sequences of interest
    min_length = 1000
    # Collection of headers for output files
    utheaders = ["Hgene", "Henst", "Hchr", "Hstrand", "Tgene", "Tenst", "Tchr", "Tstrand", "Seq"]

    blat_of_interest = ["qStart", "qEnd", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts"]
    ht_blat_headers = [f"H{bOi}" for bOi in blat_of_interest] + [f"T{bOi}" for bOi in blat_of_interest]
    enst_of_interest = ["enstURL", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "exonFrames", "cdsStartStat", "cdsEndStat", "name2"]
    ht_enst_headers = [f"H{eOi}" for eOi in enst_of_interest] + [f"T{eOi}" for eOi in enst_of_interest]

    utheadurl = utheaders + ["BlatURL", "HenstURL", "TenstURL", ]
    utheaderrors = utheaders + ["Error", "Type"]

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

            enst_return = ensting(blat, henst, tenst)

            enst: pandas.DataFrame = enst_return[0]

            if isinstance(enst, pandas.DataFrame):
                if enst.shape[0] < 2:
                    print("$$$ Not enough enst $$$")
                    for i, url in enumerate(enst_return[1]):
                        row_of_interest[f"URL_{i}"] = url

                    # OK so ZLAT, which stands for Zero BLAT, might not be the best location for this, but IDC.
                    row_of_interest.to_frame().T.to_csv(out_zlat, header = None, index = None, mode = 'a')

                elif enst.shape[0] > 2:
                    print("$$$ To much enst $$$")
                    for i, url in enumerate(enst_return[1]):
                        row_of_interest[f"URL_{i}"] = url

                    # Again maybe not the best place for this but why not?
                    row_of_interest.to_frame().T.to_csv(out_three, header = None, index = None, mode = 'a')


                elif enst.shape[0] == 2:
                    print("~~~ Clean ENST ~~~")
                    blat = blat[blat_of_interest]
                    blat.index = list(enst.index)

                    dump_frame = pandas.concat([blat, enst], axis = 1)
                    # print(enst)

                    # Take the dump frame and make it a series, then add this to the row of interest
                    # Write the row of interest to the out file. How to handle the headers? If it's the first one, add the headers, else don't

                    dows, _ = dump_frame.shape
                    for dow in range(dows):
                        dow_of_interest: pandas.Series = dump_frame.iloc[dow, :].squeeze()
                        hORt = dow_of_interest.name
                        new_index = {old_index:f"{hORt}{old_index}" for old_index in dow_of_interest.index}
                        dow_of_interest = dow_of_interest.rename(index = new_index)

                        row_of_interest = pandas.concat([row_of_interest, dow_of_interest])

                    # print(row_of_interest.index)
                    # print(row_of_interest.T)
                    row_of_interest = row_of_interest[utheadurl + ht_blat_headers + ht_enst_headers]
                    # print(row_of_interest)
                    row_of_interest.to_frame().T.to_csv(out_blat, header = None, index = None, mode = 'a')
            
            else:
                # this should dump the error
                print("$$$ Printing Error to Error files $$$")
                row_of_interest["Error"] = enst
                row_of_interest["Type"] = type(enst)
                for i, url in enumerate(enst_return[1]):
                    row_of_interest[f"URL_{i}"] = url

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



def enst_attemp(enstURL, error_attempts: int = 5):
    '''
    '''

    for attempt in range(error_attempts):

        try:
            enst_frame: pandas.DataFrame = api.convert2frame(RQuery.query(enstURL))
            break

        except AttributeError as e:
            print("$$$$ Attribute Error $$$$")
            print("Possibly wrong blat information")
            enst_frame = e
            break
        
        except Exception as e:
            print("$$$$   New ENST Error   $$$$")

            print(f"Error: {e}\nType: {type(e)}")

            print(traceback.format_exc())

            exit()

    else:
        print("$$$$ Failed to convert to frame $$$$")
        exit()

    return enst_frame




def ensting(blat_input: pandas.DataFrame, henst: str, tenst: str) -> pandas.Series or None or tuple:
    '''
    This does quit a bit: it does two things.

    First it pulls up the enst track and identifies which gene we're looking at. BLAT does not do this, it only aligns it to the genome. Looking at the ENST track from the BLAT
    allows us to identify which specific gene we're looking at.

    Next it combines this data with the BLAT data to make a smaller, useful subset. Why do it like this? Because again, BLAT only aligns to the human genome, it does not tell you which
    gene you're actually looking at. The BLAT data needs to be preserved, but it needs to be done via the identify. So at some point I was going to have to combine the ENST and the 
    BLAT data. Might as well do it here.

    '''

    # I'm doing this as a tuple, because a tuple preserves the order. The order is important. I can switch to a set at a later moment, but I'll keep it as a tuple most often
    # string2tuple = lambda input_str: tuple(re.split(',', input_str))

    blat = blat_input.copy()
    # print(blat.T)
    # exit()

    enst_index = []
    rows, _ = blat.shape
    # blat_of_interest = ["tStart", "blockCount", "blockSizes", "tStarts"]
    enst_of_sets = ['exonStarts', 'exonEnds', 'exonFrames']
    enst_of_interest = ["enstURL", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount"] + enst_of_sets + ["name2", "cdsStartStat", "cdsEndStat"]

    # head_names, tail_names = [f"H{index}" for index in blat_of_interest + enst_of_interest], [f"T{index}" for index in blat_of_interest + enst_of_interest]

    # order_index = head_names + tail_names

    # dump_frame = pandas.DataFrame(dtype='str')
    # return_series = pandas.Series(dtype='str')
    enst_index, url_index = [], []

    enst_frame = pandas.DataFrame()
    continuing = True


    for row in range(rows):
        row_of_interest: pandas.Series = blat.iloc[row, :].squeeze()

        enstURL = api.ens_tracks(chrom=row_of_interest["tName"], start=row_of_interest["tStart"], end=row_of_interest["tStart"] + row_of_interest["blockSizes"][0])
        url_index.append(enstURL)

        local_enst_frame: pandas.DataFrame = enst_attemp(enstURL)

        if isinstance(local_enst_frame, pandas.DataFrame):
            local_enst_frame["enstURL"] = enstURL
            enst_frame = pandas.concat([enst_frame, local_enst_frame])
        else:
            # Doing it like this because I want both URLs
            continuing = False

    if continuing:

        # print("!!! Continuing the ENSTing !!!!")

        enst_frame = enst_frame[(enst_frame["name"] == henst) | (enst_frame["name"] == tenst)]
        # Updates index: it's very common to have the head and tail have the same index
        enst_frame.reset_index(drop = True, inplace = True)

        enst_index = enst_frame.index

        # changing index to head and tail. Because we're preserving their positions, we can later pull this index and apply it to the blat to identify which is the head
        # and which is the tail
        new_index = {old_index: ("H" if enst_frame.loc[old_index, "name"] in henst else "T" if enst_frame.loc[old_index, "name"] in tenst else f"IDK{old_index}") for old_index in enst_index}

        enst_frame = enst_frame.rename(index = new_index)

        # print(enst_frame)
        # exit()

        for index in list(enst_frame.index):
            for sets in enst_of_sets:
                # print(enst_frame.loc[index, sets])
                enst_frame.loc[index, sets] = RQuery.convert2list(enst_frame.loc[index, sets])

        # enst_frame.rename(index={old_index: new_index}, inplace=True)
        enst_frame = enst_frame[enst_of_interest]
        # print(enst_frame.T)
        # print(blat.T)
        # exit()
    
    else:
        # This should now have an error in it
        enst_frame = local_enst_frame

    return enst_frame, url_index

        

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