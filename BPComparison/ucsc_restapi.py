import json
import requests
import pandas


def base_urls():
    '''
    '''

    base_url = f'https://api.genome.ucsc.edu/'
    list_tracks = "/list/tracks"
    get_track = "getData/track"
    get_sequence = "getData/sequence"

    return base_url, list_tracks, get_track, get_sequence


def ens_tracks(genome: str = "hg19", chrom: str = None, start: int = None, end: int = None) -> str:
    '''
    For creating a ENS Track URL. Returns a str of the url.

    chrom is the Chromosome in 3 character & 1-2 digit name for the chromosome, i.e. chr1 or chrX or chr15. 
    This should be accompianed by the start and end location on the Chromosome for the desiered track. Yes
    you this will return a url without that info (and it probably would work), but it would be a cluttered 
    mess of a response.
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


def sequence(genome: str = "hg19", chrom: str = None, start: int = None, end: int = None) -> str:
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

    seqURL = ens_tracks(chrom = "chr1", start = 94140169, end = 94140169 + 317)
    print(seqURL)

