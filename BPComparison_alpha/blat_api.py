from importlib.resources import path
import pandas
import pathlib
import requests
import json
import re
import RQuery
import logging
import traceback
import CommonMethods as CM
from inspect import currentframe, getframeinfo
from datetime import datetime

# Establish logging: because it's better then print
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(name)s:%(message)s')

start_time = datetime.now()
start_time = f"{start_time.strftime(f'%H%M%S')}_{start_time.strftime(f'%d%m%Y')}"
log_folder = pathlib.Path.cwd() / "LogFiles"

if not log_folder.is_dir():
    log_folder.mkdir(parents=True, exist_ok=True)


file_handler = logging.FileHandler(log_folder /f'BlatAPI_{start_time}.log')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)




def main():
    '''
    '''

    # with open(pathlib.Path.cwd().parent / "Data_Files" / "UTData_cds_test.csv") as csvfile:
    #     utdata = pandas.read_csv(csvfile, header = 0)

    # utdata = utdata[utdata["SeqLen"] > 100]

    # blat_queries(utdata)

    bcar_acab_seq = 'ATGGCTGCAGGAAAATTTGCAAGCCTTCCCAGAAACATGCCGGTGAATCACCAGTTCCCCCTGGCCTCATCCATGGACCTTCTGAGCAGCAGGTCCCCTCTCGCTGAGCATCGCCCAGATGCCTATCAAGATGTGTCTATACATGGCACCCTTCCACGGAAGAAAAAAGGTCCTCCTCCCATAAGGTCCTGTGATGACTTCAGTCACATGGGCACCCTCCCCCACTCCAAATCCCCACGGCAGAACTCGCCTGTGACCCAGGATGGCATCCAGGAGAGCCCATGGCAGGACCGGCACGGCGAAACCTTCACCTTCAGGGATCCACATCTTCTGGACCCAACTGTGGAATATGTGAAGTTCTCCAAGGAGAGGCACATCATGGACAGGACCCCCGAGAAACTGAAGAAGGAGCTGGAGGAGGAGCTGCTCCTGAGCAGCGAGGACCTGCGCAGCCATGCCTGGTACCACGGCCGCATCCCCCGACAGGTGTCTGAAAACCTTGTGCAGCGAGATGGTGACTTCCTAGTTCGTGACTCTCTGTCCAGCCCTGGGAACTTTGTCCTGACCTGTCAGTGGAAGAACCTCGCTCAGCACTTCAAAATCAACCGGACAGTTCTGCGACTCAGCGAGGCCTACAGCCGCGTGCAGTACCAGTTCGAGATGGAGAGCTTCGACTCCATCCCCGGCCTGGTGCGCTGCTACGTGGGCAACCGCCGGCCCATCTCCCAGCAGAGTGGCGCCATCATCTTCCAGCCCATCAACAGGACGGTGCCTCTGCGGTGCCTGGAGGAGCATTATGGCACCTCCCCAGGCCAGGCCCGGGAGGGCAGCCTCACCAAGGGAAGGCCGGATGTGGCCAAGAGGCTGAGCCTCACCATGGGTGGCGTCCAGGCCCGAGAGCAGAATTTGCCCAGGGGAAACCTCCTCAGGGAGTACCCCTGTGGCAACTCAACACCCTGGAAGACTCCTTCTGTGTCCCCAAACATCACCCAGCTGTTCCAGAAGCAGAAATGGACACAGGTCAACCCTTCACCATCCTGCAGGTGCAGCACCAGGGAGAAGCTCACCATGCTGCCAGAGTGCCCCGAGGGTGCCGGGGGCCTCCCGCCCCCCCAGAGAACACAGCGCAGCACGGAAATTCTACAAGACCTGACGGACAGGAACATCTCCGACTTCTTGGTAAAAACGTATCCTGCTCTTATAAGAAGCAGCTTAAAGAGCAAATTCTGGGTCAATGAACAGAGGTATGGAGGAATTTCCATTGGAGGAAAGCTCCCAGTCGTCCCCATCACGGGGGAAGCACTTGTTGGGTTTTTAAGCGACCTTGGCCGGATCATGAATGTGAGCGGGGGCCCTATCACTAGAGAGGCCTCTAAAGAAATACCTGATTTCCTTAAACATCTAGAAACTGAAGACAACATTAAGGTGTGGTTTAATAACAAAGGCTGGCATGCCCTGGTCAGCTTTCTCAATGTGGCCCACAACGCCATCTTACGGGCCAGCCTGCCTAAGGACAGGAGCCCCGAGGAGTATGGAATCACCGTCATTAGCCAACCCCTGAACCTGACCAAGGAGCAGCTCTCAGAGATTACAGTGCTGACCACTTCAGTGGATGCTGTGGTTGCCATCTGCGTGATTTTCTCCATGTCCTTCGTCCCAGCCAGCTTTGTCCTTTATTTGATCCAGGAGCGGGTGAACAAATCCAAGCACCTCCAGTTTATCAGTGGAGTGAGCCCCACCACCTACTGGGTGACCAACTTCCTCTGGGACATCATGAATTATTCCGTGAGTGCTGGGCTGGTGGTGGGCATCTTCATCGGGTTTCAGAAGAAAGCCTACACTTCTCCAGAAAACCTTCCTGCCCTTGTGGCACTGCTCCTGCTGTATGGATGGGCGGTCATTCCCATGATGTACCCAGCATCCTTCCTGTTTGATGTCCCCAGCACAGCCTATGTGGCTTTATCTTGTGCTAATCTGTTCATCGGCATCAACAGCAGTGCTATTACCTTCATCTTGGAATTATTTGAGAATAACCGGACGCTGCTCAGGTTCAACGCCGTGCTGAGGAAGCTGCTCATTGTCTTCCCCCACTTCTGCCTGGGCCGGGGCCTCATTGACCTTGCACTGAGCCAGGCTGTGACAGATGTCTATGCCCGGTTTGGTGAGGAGCACTCTGCAAATCCGTTCCACTGGGACCTGATTGGGAAGAACCTGTTTGCCATGGTGGTGGAAGGGGTGGTGTACTTCCTCCTGACCCTGCTGGTCCAGCGCCACTTCTTCCTCTCCCAATGGATTGCCGAGCCCACTAAGGAGCCCATTGTTGATGAAGATGATGATGTGGCTGAAGAAAGACAAAGAATTATTACTGGTGGAAATAAAACTGACATCTTAAGGCTACATGAACTAACCAAGATTTATCCAGGCACCTCCAGCCCAGCAGTGGACAGGCTGTGTGTCGGAGTTCGCCCTGGAGAGTGCTTTGGCCTCCTGGGAGTGAATGGTGCCGGCAAAACAACCACATTCAAGATGCTCACTGGGGACACCACAGTGACCTCAGGGGATGCCACCGTAGCAGGCAAGAGTATTTTAACCAATATTTCTGAAGTCCATCAAAATATGGGCTACTGTCCTCAGTTTGATGCAATTGATGAGCTGCTCACAGGACGAGAACATCTTTACCTTTATGCCCGGCTTCGAGGTGTACCAGCAGAAGAAATCGAAAAGGTTGCAAACTGGAGTATTAAGAGCCTGGGCCTGACTGTCTACGCCGACTGCCTGGCTGGCACGTACAGTGGGGGCAACAAGCGGAAACTCTCCACAGCCATCGCACTCATTGGCTGCCCACCGCTGGTGCTGCTGGATGAGCCCACCACAGGGATGGACCCCCAGGCACGCCGCATGCTGTGGAACGTCATCGTGAGCATCATCAGAGAAGGGAGGGCTGTGGTCCTCACATCCCACAGCATGGAAGAATGTGAGGCACTGTGTACCCGGCTGGCCATCATGGTAAAGGGCGCCTTTCGATGTATGGGCACCATTCAGCATCTCAAGTCCAAATTTGGAGATGGCTATATCGTCACAATGAAGATCAAATCCCCGAAGGACGACCTGCTTCCTGACCTGAACCCTGTGGAGCAGTTCTTCCAGGGGAACTTCCCAGGCAGTGTGCAGAGGGAGAGGCACTACAACATGCTCCAGTTCCAGGTCTCCTCCTCCTCCCTGGCGAGGATCTTCCAGCTCCTCCTCTCCCACAAGGACAGCCTGCTCATCGAGGAGTACTCAGTCACACAGACCACACTGGACCAGGTGTTTGTAAATTTTGCTAAACAGCAGACTGAAAGTCATGACCTCCCTCTGCACCCTCGAGCTGCTGGAGCCAGTCGACAAGCCCAGGACTGA'
    bcar_acab_url = gen_url(bcar_acab_seq)
    print(bcar_acab_url)

    # test = blat_query(bcar_acab_url)

    # if isinstance(test, requests.models.Response):
    #     blat: dict = json.loads(test.text)
    #     britta = blat['blat']

    #     att_people = {}
    #     for i, b in enumerate(britta):
    #         gdb = dict(zip(blat['fields'], b))
    #         att_people[f"{i}"] = gdb

    #     for qindex, qresponse in att_people.items():
    #         print(f"Blat index {qindex}")
    #         for fieldname, fieldvalue in qresponse.items():
    #             print(f"\t{fieldname}: {fieldvalue}")




def blat_queries(dataframe: pandas.DataFrame, qtype = "DNA", qdb = "hg19") -> pandas.DataFrame:
    '''
    For sending in a ton of different sequences via a dataframe

        hg19 comes back with:
    track -> blat
    genome -> set by user
    fields: list
    blat: list of list
        matches
        misMatches
        repMatches
        nCount
        qNumInsert
        qBaseInsert
        tNumInsert
        tBaseInsert
        -> strand <-        which strand of DNA it's on, will be used for reference
        qName
        qSize
        qStart
        qEnd
        -> tName <-         which chrm it is on, will be used for reference
        tSize
        tStart
        tEnd
        -> blockCount <-    how many blocks it aligns to, mostly exons, useful for finding the cds
        -> blockSizes <-    how big each block is, useful for finding end. By default this is a string, should be a list/set/tuple
        qStarts
        -> tStarts <-       where each block starts, useful for finding the sequence. By default this is a string, should be a list/set/tuple

        note these are ALL read left to right, but the - strand is actually RIGHT to LEFT, so the 
        chrm region is REVERSED

        They seem to be off by 1 nucleotide

        it's all -1 in the start position: does the web api index start at 0?
        Programically it appears to start at 0.

        Yes, it most certinaly appears to start index at 0, and it reports the correct sequence when this is put
        into the getData/sequence url. Now to find a way to grab the correct ENST

    '''

    rows, _ = dataframe.shape

    for row in range(rows):
        row_of_interest = dataframe.iloc[row, :]

        qseq = row_of_interest["Seq"]

        query_frame = blat_query(qseq)
        
        print(type(query_frame))
        print(query_frame)
        if row == 0:
            break



def blat_query(qseq: str = None, qurl: str = None, qdb = 'hg19', qtype = 'DNA') -> pandas.DataFrame:
    '''
    If this receives a sequence, then it generates a URL
    Else (if the sequence is none and) IF this receives a url, then it will use the url
    So if a sequence is sent in this will generate a new URL EVEN IF one is already provided
    '''

    if isinstance(qseq, str):
        query_url = gen_url(qseq, qdb=qdb, qtype=qtype)
    
    elif isinstance(qurl, str):
        query_url = qurl        

    try:
        query = RQuery.query(query_url)
    except UnboundLocalError as e:
        print("!!!!\tUnbound Local Error\t!!!!")
        logger_output(message_title="Unbound Local Error when trying to query UCSC Database", data=f"Query URL:\n{query_url}")

        query = None

    except Exception as e:
        print("!!!!\tNew Error in Blat API\t!!!!")
        logger_output(message_title="New Error when trying to query UCSC Database", data=f"Error: {e}\n\tType: {type(e)}\nQuery URL:\n{query_url}")

        query = None


    if isinstance(query, requests.models.Response):
        return_data = convert2frame(query)

    else:
        return_data = None

    return return_data


def convert2frame(query: requests.Response) -> pandas.DataFrame:
    '''
    Converts the blat query into a Dataframe. This is somewhat different from the rest api
    because of layout. The blat API returns

    track -> blat
    genome -> set by user
    fields: list
    blat: list of list
        matches
        misMatches
        repMatches
        nCount
        qNumInsert
        qBaseInsert
        tNumInsert
        tBaseInsert
        -> strand <-        which strand of DNA it's on, will be used for reference
        qName
        qSize
        qStart
        qEnd
        -> tName <-         which chrm it is on, will be used for reference
        tSize
        tStart
        tEnd
        -> blockCount <-    how many blocks it aligns to, mostly exons, useful for finding the cds
        -> blockSizes <-    how big each block is, useful for finding end. By default this is a string, should be a list/set/tuple
        qStarts
        -> tStarts <-       where each block starts, useful for finding the sequence. By default this is a string, should be a list/set/tuple

    and this is converted to a dataframe. Also note: this is for hg19, if the genome is not hg19 I'm not sure what
    the response will be.
    '''
    blat: dict = json.loads(query.text)
    britta = blat['blat']

    att_people = {}
    for i, b in enumerate(britta):
        gdb = dict(zip(blat['fields'], b))
        gdb["blockSizes"] = CM.convert2list(gdb["blockSizes"])
        gdb["tStarts"] = CM.convert2list(gdb["tStarts"])
        gdb["qStarts"] = CM.convert2list(gdb["qStarts"])

        att_people[f"{i}"] = gdb

    blat_frame = pandas.DataFrame()

    for qindex, qresponse in att_people.items():
        qseries = pandas.Series(qresponse, name = qindex)
        # print(qseries)
        blat_frame = pandas.concat([blat_frame, qseries.to_frame()], axis = 1)
    
    blat_frame = blat_frame.T

    return blat_frame



def gen_url(qseq, qdb = 'hg19', qtype = 'DNA') -> str:
    '''
    Just because I want to see the URL, and personally I don't like returning variable number of tuples
    '''

    main_url = f"https://genome.ucsc.edu/cgi-bin/hgBlat?"

    qtype, qdb, qoutput = f"type={qtype}", f"db={qdb}", f"output=json"

    qseq = f"userSeq={qseq}"

    query_url = f"{main_url}{qseq}&{qtype}&{qdb}&{qoutput}"

    return query_url


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


        


if __name__ in '__main__':
    main()