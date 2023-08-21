from cmath import exp
from textwrap import wrap
import requests
from requests.adapters import HTTPAdapter
import logging
import traceback
from inspect import currentframe, getframeinfo
import re
import time
import urllib3

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


def query_response(url: str, max_connctions: int = 2, timeout: float = 30):
    '''
    Gets the actual querry.

    The requests exceptions Connect Timeout error should pause the script for 10 minutes and then try again to get the request. I think that's what is happeneing: too many requests from the same IP.
    '''
    session = requests.Session()
    retry = urllib3.util.retry.Retry(connect = max_connctions, backoff_factor = timeout)
    adapter = HTTPAdapter(max_retries = retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)

    query = session.get(url)

    # query: requests.models.Response = requests.get(url, timeout = timeout)

    # for _ in range(max_iterations):
    #     try:
    #         query: requests.models.Response = requests.get(url, timeout = timeout)
    #         break
        
    #     except requests.exceptions.ConnectTimeout as e:
    #         print(f"\tPausing for 10 minutes due to HTTP Error")
    #         logger_output(message_title="#### HTTP Error", data=f"--- URL: {url}\nLine {getframeinfo(currentframe()).lineno -5}?\nError: {e}\n\tType: {type(e)}\n\nTraceback:\n{traceback.format_exc()}\n### End Log ###\n\n")
    #         time.sleep(600)
    #     except requests.exceptions.ReadTimeout as e:
    #         print("\tRead Timeout: Pausing for 10 minutes")
    #         logger_output(message_title="#### Read Timeout", data=f"--- URL: {url}\nLine {getframeinfo(currentframe()).lineno -5}?\nError: {e}\n\tType: {type(e)}\n\nTraceback:\n{traceback.format_exc()}\n### End Log ###\n\n")
    #         time.sleep(60)
    #     except Exception as e:
    #         query = None
    #         logger_output(message_title="#### New Exception ####", data=f"--- URL: {url}\nLine {getframeinfo(currentframe()).lineno -9}?\nError: {e}\n\tType: {type(e)}\n\nTraceback:\n{traceback.format_exc()}\n### End Log ###\n\n")
    #         break
    # else:
    #     logger_output(message_title="Failed to get query", data=f"--- URL: {url}\nLine {getframeinfo(currentframe()).lineno -13}?\nError: {e}\n\tType: {type(e)}\n\nTraceback:\n{traceback.format_exc()}\n### End Log ###\n\n")
    #     query = None

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