import json
from re import T
from urllib.request import Request
import requests
import pandas
import RQuery
from typing import Tuple


def base_urls():
    '''
    '''

    base_url = f'https://api.genome.ucsc.edu/'
    list_tracks = "list/tracks"
    get_track = "getData/track"
    get_sequence = "getData/sequence"

    return base_url, list_tracks, get_track, get_sequence


def compliment_strand(sequence: str) -> str:
    '''
    '''
    compliment_sequence = ''

    for n in sequence:
        if n in "A":
            compliment_sequence = f"{compliment_sequence}T"
        if n in "C":
            compliment_sequence = f"{compliment_sequence}G"
        if n in "G":
            compliment_sequence = f"{compliment_sequence}C"
        if n in "T":
            compliment_sequence = f"{compliment_sequence}A"

    return compliment_sequence


def ens_tracks(genome: str = "hg19", chrom: str = None, start: int = None, end: int = None) -> Tuple[pandas.DataFrame, str]:
    '''
    For creating a ENS Track URL. Returns a str of the url.

    chrom is the Chromosome in 3 character & 1-2 digit name for the chromosome, i.e. chr1 or chrX or chr15. 
    This should be accompianed by the start and end location on the Chromosome for the desiered track. Yes
    you this will return a url without that info (and it probably would work), but it would be a cluttered 
    mess of a response.
    '''

    ens_url = enst_tracks_url(genome = genome, chrom = chrom, start = start, end = end)

    try:
        query = RQuery.query(ens_url)
    except UnboundLocalError as e:
        print("!!!!\tUnbound Local Error\t!!!!")
        # logger_output(message_title="Unbound Local Error when trying to query UCSC Database", data=f"Query URL:\n{query_url}")

        query = e

    except Exception as e:
        print("!!!!\tNew Error in Blat API\t!!!!")
        # logger_output(message_title="New Error when trying to query UCSC Database", data=f"Error: {e}\n\tType: {type(e)}\nQuery URL:\n{query_url}")

        query = e

    if isinstance(query, requests.models.Response):
        return_data: pandas.DataFrame = convertRequest(query)

    else:
        return_data = query

    return return_data, ens_url



def enst_tracks_url(genome: str = "hg19", chrom: str = None, start: int = None, end: int = None) -> str:
    '''
    '''
    # can't really partial unpack here because of it's position
    base_url, _, get_track, _ = base_urls()
    ens_url = "track=ensGene"

    genome = f"genome={genome}"

    if (chrom is not None) and (start is not None) and (end is not None):
        chrom, start, end = f"chrom={chrom}", f"start={start}", f"end={end}"
         
        ens_url = f"{base_url}{get_track}?{ens_url};{genome};{chrom};{start};{end}"

    else:
        ens_url = f"{base_url}{get_track}?{ens_url};{genome}"

    return ens_url



def geneid_track(genome: str = "hg19", chrom: str = None, start: int = None, end: int = None) -> str:
    '''
    '''

    base_url, _, get_track, _ = base_urls()
    geneid_url = "track=geneid"

    genome = f"genome={genome}"

    if (chrom is not None) and (start is not None) and (end is not None):
        chrom, start, end = f"chrom={chrom}", f"start={start}", f"end={end}"
         
        geneid_url = f"{base_url}{get_track}?{geneid_url};{genome};{chrom};{start};{end}"

    else:
        geneid_url = f"{base_url}{get_track}?{geneid_url};{genome}"

    return geneid_url



def sequence(genome: str = "hg19", chrom: str = None, start: int = None, end: int = None, strand: str = None, reverse_dir = False) -> Tuple[str, str]:
    '''
    '''

    seqURL = sequence_url(genome = genome, chrom = chrom, start = start, end = end)

    try:
        query = RQuery.query(seqURL)
    except UnboundLocalError as e:
        print("!!!!\tUnbound Local Error\t!!!!")
        # logger_output(message_title="Unbound Local Error when trying to query UCSC Database", data=f"Query URL:\n{query_url}")

        query = e

    except Exception as e:
        print("!!!!\tNew Error in Blat API\t!!!!")
        # logger_output(message_title="New Error when trying to query UCSC Database", data=f"Error: {e}\n\tType: {type(e)}\nQuery URL:\n{query_url}")

        query = e

    if isinstance(query, requests.models.Response):
        query = convertRequest(query)

    if strand in "-":
        query = compliment_strand(query)
        # query = query[::-1]

    if reverse_dir:
        query = query[::-1]

    return query, seqURL




def sequence_url(genome: str = "hg19", chrom: str = None, start: int = None, end: int = None) -> str:
    '''
    For getting a specific sequence.
    '''
    base_url, *_, get_sequence = base_urls()

    genome = f"genome={genome}"

    if (chrom is not None) and (start is not None) and (end is not None):
        chrom, start, end = f"chrom={chrom}", f"start={start}", f"end={end}"

        seq_url = f"{base_url}{get_sequence}?{genome};{chrom};{start};{end}"

    else:
        seq_url = f"{base_url}{get_sequence}?{genome}"

    return seq_url


def convert2frame(track: list) -> pandas.DataFrame:
    '''
    '''
    if isinstance(track, list):
        # Because the UCSC people thought it would be HILarious if they made this a list...
        track: dict = dict(zip(range(len(track)), track))

    track_frame = pandas.DataFrame()
    for key, t in track.items():
        track_data = pandas.Series(t, name = key)

        track_frame = pandas.concat([track_frame, track_data.to_frame()], axis = 1)

    # Yes I know there is a way I can transpose the stuff as I work on it, but this is just easier.
    track_frame = track_frame.T

    return track_frame


def convertRequest(query: requests.Response) -> pandas.DataFrame or str:
    '''
    '''
    track: dict = json.loads(query.text)

    try:
        # Litterally don't care how many items are returned, we'll grab that later. Only need to know what index to grab
        # the data from
        _ = track["itemsReturned"]
        index_key = -2
    except KeyError as e:
        index_key = -1
    except Exception as e:
        index_key = -1
        print(f"Exception Raised: {e}\tType: {type(e)}\nSetting index key to -1 and trying")

    track = track[tuple(track.keys())[index_key]]

    if isinstance(track, list):
        # If it's a list: we need to convert it to something more usable
        track: pandas.DataFrame = convert2frame(track)
    elif isinstance(track, str):
        # If it's a string, it probably already is usable
        track: str = track

    return track




if __name__ in '__main__':
    # url_gen = UCSCURLGenerator()

    # ensURL = url_gen.ens_tracks(chrom="chr15", start=65849223, end = 65849223 + 18)
    # print(ensURL)

    # _, seqURL = ens_tracks(chrom = "chr1", start = 94140169, end = 94140169 + 317)
    # print(seqURL)

    # baseURL, trackURL, *_ = base_urls()
    # trackURL = f"{baseURL}{trackURL}?genome=hg19"
    # print(trackURL)

    geneid_url = geneid_track(chrom = "chr1", start = 94140169, end = 94140169 + 317)
    # print(geneid_url)

