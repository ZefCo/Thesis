from datetime import datetime
from doctest import master
from threading import local
from requests import Response
import blat_api as ba
import ucsc_restapi as api
import pathlib
import pandas
import RQuery
import json
import re
import CommonMethods as CM
# import rpy2
# from rpy2 import robjects
# from rpy2.robjects.packages import STAP
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.conversion import localconverter


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

start_time = datetime.now()
start_time = f"{start_time.strftime(f'%H%M%S')}_{start_time.strftime(f'%d%m%Y')}"
log_folder = pathlib.Path.cwd() / "LogFiles"

if not log_folder.is_dir():
    log_folder.mkdir(parents=True, exist_ok=True)

file_handler = logging.FileHandler(log_folder / f'EEBP_{start_time}.log')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)



class CmRNA():
    def __init__(self, file_path: str or pathlib.Path, error_attemps = 5) -> None:
        self.error_attempts = error_attemps
        
        self.master_headers = ["Hgene", "Henst", "Hchr", "Hstrand", "Tgene", "Tenst", "Tchr", "Tstrand", "Seq", "SeqLen"]
        self.blat_of_interest = ["qStart", "qEnd", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts"]
        self.hlat_headers, self.tlat_headers = [f"H{bOi}" for bOi in self.blat_of_interest], [f"T{bOi}" for bOi in self.blat_of_interest]

        self.enst_of_sets = ['exonStarts', 'exonEnds', 'exonFrames']
        self.enst_of_interest = ["enstURL", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount"] + self.enst_of_sets + ["name2", "cdsStartStat", "cdsEndStat"]


        self.data = self.import_data(file_path)
        # self.data: robjects.DataFrame = self.rFun.importData(str(file_path))
        # pass



    def blat(self, fusionData: pandas.DataFrame, min_length = 20, start_row: int = ..., end_row: int = ...):
        '''
        
        '''

        if isinstance(start_row, int):
            start_row = start_row
        else:
            start_row = 0

        # master_headers = ["Hgene", "Henst", "Hchr", "Hstrand", "Tgene", "Tenst", "Tchr", "Tstrand", "Seq", "SeqLen"]

        try:
            fusionData = fusionData[self.master_headers]
        except Exception as e:
            print(f"!!!!\tNew Error at line {getframeinfo(currentframe()).lineno - 2}\t!!!!")
            self.logger_output(message_title="Error when trying to filter data", data=f"Requested Master Headers (note: Python list type, not R Str Vector):\n{self.master_headers}\nDataFrame Headers:\n{tuple(fusionData.columns)}Error: {e}\tType: {type(e)}\nTraceback:\n{traceback.format_exc()}")

        fusionData = fusionData[fusionData["SeqLen"] >= min_length]

        rows, _ = fusionData.shape
        if (isinstance(end_row, int)) and end_row > start_row:
            end_row = end_row
        else:
            end_row = rows

        for row in range(start_row, end_row):
            row_of_interest: pandas.Series = fusionData.iloc[row, :].copy().squeeze()
            local_blat: pandas.DataFrame = self.blat_data(data_to_blat=row_of_interest)

            if local_blat.shape[0] == 2:
                print("~~~~\tClean BLAT\t~~~~")

                row_of_interest["BlatURL"] = local_blat.loc[0, "BlatURL"]

                local_enst: pandas.DataFrame = self.enst_indentification(local_blat)

            elif local_blat.shape[0] > 2:
                print("~~~~\tToo much BLAT\t~~~~")
                
                row_of_interest["BlatURL"] = local_blat.loc[0, "BlatURL"]
                self.error_outputs(row_of_interest, error_file="TooMuchBlat.csv")

                continue

            elif local_blat.shape[0] < 2:
                print("~~~~\tToo little BLAT\t~~~~")

                row_of_interest["BlatURL"] = local_blat.loc[0, "BlatURL"]
                self.error_outputs(row_of_interest, error_file="ZeroBlat.csv")

                continue


            if isinstance(local_enst, tuple):
                print("!!!!\tError during ENST Identification\t!!!!")
                print("    \tSending to logs")
                self.error_outputs(row_of_interest, error_file="ENST_Error.csv", return_error=local_enst)

            elif isinstance(local_enst, pandas.DataFrame):

                local_enst = local_enst[(local_enst["name"] == row_of_interest["Henst"]) | (local_enst["name"] == row_of_interest["Tenst"])]

                print(f"~~~~\tRow of Interest\t~~~~\n{row_of_interest}")
                print(f"~~~~\tLocal BLAT\t~~~~\n{local_blat.T}")
                print(f"~~~~\tLocal Enst\t~~~~\n{local_enst.T}")


                

            print(f"Finished row {row} of {end_row}")




    def blat_attempt(self, blat_target: str, hgene: str = None, henst: str = None, tgene: str = None, tenst: str = None) -> pandas.DataFrame or None:
        '''
        this block attempts to blat the sequence. If for some reason it doesn't blat it sends the url to an error file
        '''
        for attempt in range(self.error_attempts):
            try:
                blat: pandas.DataFrame = ba.blat_query(qseq = blat_target)
                if attempt > 0:
                    self.logger_output(message_title = "Errors at blat", data = f"{attempt} number of attempts to blat the following CmRNA:\n{hgene}_{tgene}\n{henst}_{tenst}")
                break

            except json.decoder.JSONDecodeError as e:
                print("!!!!\tJSON Decoder Error During Blat Attempt\t!!!!")
                print("    \tretrying")
                blat = e

            except Exception as e:
                print(f"!!!!\tNew Error at line {getframeinfo(currentframe()).lineno - 11} During Blat Attempt\t!!!")
                self.logger_output(message_title="New Error at BLAT", data=f"Fusion: {hgene}_{tgene}\t{henst}_{tenst}\nError: {e}\n\tType: {type(e)}\n\nTraceback: {traceback.format_exc()}")
                blat = e
                
                break
        
        else:
            print("!!!!\tFailed to blat, sending to log\t!!!!")
            self.logger_output(message_title=f"Failed to Blat of the following CmRNA after {attempt + 1} number of attempts", data=f"Fusion: {hgene}_{tgene}\t{henst}_{tenst}\nError: {blat}\n\tType: {type(blat)}")

        return blat



    def blat_data(self, data_to_blat: pandas.Series) -> pandas.DataFrame:
        '''
        '''

        hgene, henst, hstrand, hchr, tgene, tenst, tstrand, tchr = data_to_blat["Hgene"], data_to_blat["Henst"], data_to_blat["Hstrand"], data_to_blat["Hchr"], data_to_blat["Tgene"], data_to_blat["Tenst"], data_to_blat["Tstrand"], data_to_blat["Tchr"]

        print(f"\n####\n{hgene}_{tgene}\n{henst}_{tenst}\n####")

        seq: str = data_to_blat["Seq"]
        burl: str = ba.gen_url(seq)
        data_to_blat["BlatURL"] = burl

        blat: pandas.DataFrame = self.blat_attempt(seq, hgene=hgene, henst=henst, tgene=tgene, tenst=tenst)

        if isinstance(blat, pandas.DataFrame):
            blat = blat[((blat["strand"] == hstrand) & (blat["tName"] == hchr)) | ((blat["strand"] == tstrand) & (blat["tName"] == tchr))]
            blat.reset_index(drop = True, inplace = True)
            blat["BlatURL"] = burl
        
        else:
            self.error_outputs(data_to_blat, error_file_path = "BlatError.csv", return_error= blat)

        return blat



    def check_file(self, filename: str, headers: list):
        '''
        Check to see if the file exists: if it doesn't then creates it.
        '''

        out_file = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / filename

        if not out_file.is_file():
            pandas.DataFrame(columns = headers).to_csv(out_file, header = True, index = False)

        return out_file



    def enst_indentification(self, fusionData: pandas.DataFrame, start_row: int = ..., end_row: int = ..., min_length: int = 20) -> pandas.DataFrame or list:
        '''
        Need to re-write this. The BLAT doesn't tell me what gene I am looking at, only the corrdinates. I need to convert those corredinates to a ENST.

        '''

        # print(fusionData)

        # fusionData["enstURL"] = None
        enst_data: pandas.DataFrame = pandas.DataFrame()
        enst_urls: list = list()
        enst_errors: list = list()

        if isinstance(start_row, int):
            start_row = start_row
        else:
            start_row = 0

        # fusionData = fusionData[fusionData["SeqLen"] >= min_length]

        rows, _ = fusionData.shape
        if (isinstance(end_row, int)) and end_row > start_row:
            end_row = end_row
        else:
            end_row = rows

        for row in range(start_row, end_row):
            row_of_interest = fusionData.iloc[row, :].copy().squeeze()

            chr, start, end = row_of_interest["tName"], row_of_interest["tStart"], row_of_interest["tStart"] + row_of_interest["blockSizes"][0]

            enstURL = api.ens_tracks(chrom=chr, start=start, end=end)
            enst_urls.append(enstURL)

            local_enst: pandas.DataFrame = self.enst_attempt(enstURL)

            if isinstance(local_enst, pandas.DataFrame) and isinstance(enst_data, pandas.DataFrame):
                # everything is fine: append the data and go on
                local_enst["enstURL"] = enstURL
                enst_data = pandas.concat([enst_data, local_enst])
                enst_data.reset_index(drop = True, inplace = True)

            elif (not isinstance(local_enst, pandas.DataFrame)) and (isinstance(enst_data, pandas.DataFrame)):
                # everything is NOT fine: error found, turns into an error report
                enst_errors.append(local_enst)

                enst_data = None

            elif (not isinstance(local_enst, pandas.DataFrame)) and (not isinstance(enst_data, pandas.DataFrame)):
                # already have an error report, adding more to it
                enst_errors.append(local_enst)

            elif isinstance(local_enst, pandas.DataFrame) and (not isinstance(enst_data, pandas.DataFrame)):
                # already have an error report, but nothing to add to
                pass

        if not isinstance(enst_data, pandas.DataFrame):
            enst_data = (enst_urls, enst_errors)


        return enst_data


        # print(enst_data.T)
        # exit()


            # fusionData.loc[row, "enstURL"] = enstURL

        # print(fusionData.T)


    def enst_attempt(self, enstURL):
        '''
        '''

        for attempt in range(self.error_attempts):

            try:
                # something = RQuery.query(enstURL)
                # print(something.text)
                # print(type(something))
                # exit()
                enst_frame: pandas.DataFrame = api.convertRequest(RQuery.query(enstURL))
                if attempt > 0:
                    self.logger_output(message_title = "Errors at ENST identification", data = f"{attempt} number of attempts to ENST Identify the following CmRNA URL:\n{enstURL}")
                break

            except AttributeError as e:
                print(f"!!!!\tAttribute Error\t!!!!")
                print(f"    \tPossibly wrong blat information")
                enst_frame = e

                break
            
            except Exception as e:
                print(f"!!!!\tNew ENST Error\t!!!!")
                print(f"    \tPrinting to log")
                enst_frame = e

                break

        else:
            print("!!!!\tFailed to convert to frame\t!!!!")

            # exit()

        return enst_frame




    def enst_data(self, fusionData: pandas.DataFrame):
        '''

        '''




    def error_outputs(self, data: pandas.Series = ..., error_file: str = ..., return_error: Exception or list or tuple = ..., *args, **kwargs):
        '''
        For outputing to the error csv. By error I mean and any unexpected result
        '''

        if return_error is not ...:

            if isinstance(return_error, list):
                return_error = tuple(return_error)
                return_type = tuple([type(rerror) for rerror in return_error])
            else:
                return_type = type(return_type)

            data["Error"] = return_error
            data["Error Type"] = return_type

        error_file = self.check_file(error_file, list(data.index))

        data.to_frame().T.to_csv(error_file, header = None, index = None, mode = 'a')



    def import_data(self, file_path: str or pathlib.Path) -> pandas.DataFrame:
        '''
        '''
        with open(file_path) as file:
            data = pandas.read_csv(file, header = 0)

        return data



    def logger_output(self, message_title=None, data=None):
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
    cmrna = CmRNA(pathlib.Path.cwd().parent / "Data_Files" / "UTData_cds.csv")
    # print(type(cmrna.data))
    cmrna.data = cmrna.blat(cmrna.data, min_length=1000, end_row=10)


if __name__ in '__main__':
    main()