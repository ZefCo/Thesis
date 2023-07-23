from cmath import exp
from textwrap import wrap
import requests
import logging
import traceback
from inspect import currentframe, getframeinfo
import re
import time

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

    url_query: requests.models.Response = query_response(url)

    if url_query.status_code == 429:
        # Apparently UCSC thought it would be hilarious to report their wait time as a string...
        wait_time = int(url_query.headers["Retry-After"])

        print(f"Too many requests - request made to wait {wait_time} seconds - waiting {wait_time * 2} seconds")

        time.sleep(wait_time * 2)

        print(f"Trying once again - if this fails then look into fixing")

        url_query: requests.models.Response = query_response(url)

    return url_query


def query_response(url: str):
    '''
    '''
    try:
        query: requests.models.Response = requests.get(url)
    except Exception as e:
        query = None
        logger_output(message_title="New Exception", data=f"Line {getframeinfo(currentframe()).lineno -3}\nError: {e}\n\tType: {type(e)}\n\nTraceback:\n{traceback.format_exc()}")

    return query




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