from cmath import exp
import requests
import logging
import traceback
from inspect import currentframe, getframeinfo
import re

# Establish logging: because it's better then print
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(name)s:%(message)s')

file_handler = logging.FileHandler(f'RQuery.log')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)



def query(url: str) -> requests.models.Response:
    '''
    A try/except to handle the querying of UCSC
    '''
    try:
        query: requests.models.Response = requests.get(url)
    except Exception as e:
        query = None
        logger_output(message_title="New Exception", data=f"Line 24\nError: {e}\n\tType: {type(e)}\n\nTraceback:\n{traceback.format_exc()}")

    return query


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

    sequence = re.split(',', sequence)

    seqlist = []
    for strint in sequence:
        
        try:
            strint = int(strint)

        except ValueError as e:
            strint = None

        except Exception as e:
            # print(f"Error: {e}\tType: {type(e)}")
            strint = None

        seqlist.append(strint)

    seqlist = tuple(seqlist)

    return seqlist



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
    pass