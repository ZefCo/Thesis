from requests import Response
import blat_api as ba
import ucsc_restapi as api
import pathlib
import pandas
import RQuery
import json

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
    start_index = 5552
    min_length = 1000
    output_headers = []
    utheaders = ["Hgene", "Henst", "Hchr", "Hstrand", "Tgene", "Tenst", "Tchr", "Tstrand", "Seq"]
    utheadurl = utheaders + ["URL"]
    utheaderrors = utheadurl + ["Erro", ["Type"]]
    # utheaderr = utheaders + ["Error", "Type"]
    # zeroheader = utheaders + ["URL"]

    out_zlat = pregen(filename = "ZeroBlat_2.csv", headers = utheadurl)
    out_error = pregen(filename = "UT_Errors_2.csv", headers = utheaderrors)
    out_blat = pregen(filename = f"UT_Blat_{min_length}_2.csv", headers = utheadurl)
    out_three = pregen(filename = "ThreePlusBlat_2.csv", headers = utheadurl)
    
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
        row_of_interest["URL"] = qurl
        hgene, henst, hstrand, hchr, tgene, tenst, tstrand, tchr = row_of_interest["Hgene"], row_of_interest["Henst"], row_of_interest["Hstrand"], row_of_interest["Hchr"], row_of_interest["Tgene"], row_of_interest["Tenst"], row_of_interest["Tstrand"], row_of_interest["Tchr"]

        print(f"\n####\n{hgene}_{tgene}\n{henst}_{tenst}\n####")

        for attempt in range(error_attempts):
            try:
                blat: pandas.DataFrame = ba.blat_query(qurl = qurl)
                if attempt > 0:
                    logger_output(message_title = "Errors at blat", data = f"{attempt} number of attempts to blat the following url:\n{qurl}")
                break

            except json.decoder.JSONDecodeError as e:
                print("$$$$ JSON Decoder Error $$$$")
                print("          retrying          ")
                last_error_message = e
                # logger_output(message_title="JSON Decoder Error at creating BLAT", data=f"Check other logs for additional details?")

                # row_of_interest["Error"] = e
                # row_of_interest["Type"] = type(e)

                # row_of_interest.to_frame().T.to_csv(out_error, header = None, index = None, mode = 'a')
                
                # continue

            except Exception as e:
                print("$$$$      New Error     $$$$")
                logger_output(message_title="New Error at BLAT", data=f"Fusion: {hgene}_{tgene}\t{henst}_{tenst}\nError: {e}\n\tType: {type(e)}\n\nTraceback: {traceback.format_exc()}")
                
                row_of_interest["Error"] = e
                row_of_interest["Type"] = type(e)

                row_of_interest.to_frame().T.to_csv(out_error, header = None, index = None, mode = 'a')

                break
        
        else:
            print("$$$ Failed to blat, sending to log $$$")
            logger_output(message_title=f"Failed to Blat of the following URL after {attempt + 1} number of attempts", data=f"Fusion: {hgene}_{tgene}\t{henst}_{tenst}\nError: {last_error_message}\n\tType: {type(last_error_message)}")

            row_of_interest["Error"] = last_error_message
            row_of_interest["Type"] = type(last_error_message)

            row_of_interest.to_frame().T.to_csv(out_error, header = None, index = None, mode = 'a')

            continue


        blat = blat[((blat["strand"] == hstrand) & (blat["tName"] == hchr)) | ((blat["strand"] == tstrand) & (blat["tName"] == tchr))]

        if blat.shape[0] == 2:
            print("~~~ Clean Blat ~~~")
            row_of_interest.to_frame().T.to_csv(out_blat, header = None, index = None, mode = 'a')

        elif blat.shape[0] < 2:
            print("$$$ Not enough blat $$$")
            row_of_interest.to_frame().T.to_csv(out_zlat, header = None, index = None, mode = 'a')

        elif blat.shape[0] > 2:
            print("$$$ To much blat $$$")
            row_of_interest.to_frame().T.to_csv(out_three, header = None, index = None, mode = 'a')

        print(f"Finished UT Database row {row}")



def ensting():
    '''
    Going to try this with the clean blat ones first. Ambiguous ones can be done later.

    Also want to classify them: Cis-SAG, EGenic (intergenic), AGenic (intragenic)
    '''

    with open(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UT_Blat_1000_2.csv") as csvfile:
        utdata = pandas.read_csv(csvfile, header = 0)

    rows, _ = utdata.shape

    for row in range(rows):
        row_of_interest = utdata.iloc[row, :]
        hgene, henst, hstrand, hchr, tgene, tenst, tstrand, tchr = row_of_interest["Hgene"], row_of_interest["Henst"], row_of_interest["Hstrand"], row_of_interest["Hchr"], row_of_interest["Tgene"], row_of_interest["Tenst"], row_of_interest["Tstrand"], row_of_interest["Tchr"]
        blatURL = row_of_interest["URL"]
        print(f"\n####\n{hgene}_{tgene}\n{henst}_{tenst}\n####")


        for attempt in range(error_attempts):
            try:
                blat: pandas.DataFrame = ba.blat_query(qurl = blatURL)
                if attempt > 0:
                    print(f"Took {attempt} # of attempts to blat correctly")

                break

            except json.decoder.JSONDecodeError as e:
                print("$$$$ JSON Decoder Error $$$$")
                print("          retrying          ")
                last_error_message = e

            except Exception as e:
                print("$$$$      New Error     $$$$")
                print(f"Error: {e}\n\tType: {type(e)}")

                break
        
        else:
            print("$$$ Failed to blat, skipping $$$")

            continue

        blat = blat[((blat["strand"] == hstrand) & (blat["tName"] == hchr)) | ((blat["strand"] == tstrand) & (blat["tName"] == tchr))]

        if blat.shape[0] == 2:
            # print("~~~ Clean Blat ~~~")
            print(blat.T)
            bows, _ = blat.shape
            for bow in range(bows):
                bow_of_interest: pandas.Series = blat.iloc[bow, :].squeeze()
                # print(type(bow_of_interest["blockSizes"]))
                enstURL = api.ens_tracks(chrom=bow_of_interest["tName"], start=bow_of_interest["tStart"], end=bow_of_interest["tStart"] + bow_of_interest["blockSizes"][0])

                # get ENST URL, filter to only the enst being used for the fusion row.
                print(enstURL)



            # row_of_interest.to_frame().T.to_csv(out_blat, header = None, index = None, mode = 'a')
            # Need to reindex the blat frame with the ENST
            # I don't know which end I have

        elif blat.shape[0] < 2:
            print("$$$ Not enough blat $$$")
            # row_of_interest.to_frame().T.to_csv(out_zlat, header = None, index = None, mode = 'a')

        elif blat.shape[0] > 2:
            print("$$$ To much blat $$$")
            # row_of_interest.to_frame().T.to_csv(out_three, header = None, index = None, mode = 'a')

        print(f"Finished UT Database row {row}")

        if row >= 0:
            exit()

        

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
    # blating()
    ensting()




if __name__ in '__main__':
    main()