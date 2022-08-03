import requests


def query(url: str) -> requests.models.Response:
    '''
    A try/except to handle the querying of UCSC
    '''
    try:
        query: requests.models.Response = requests.get(url)
    except Exception as e:
        print("####\nException found")
        print(e)
        print(type(e))
        query = None

    return query


if __name__ in '__main__':
    pass