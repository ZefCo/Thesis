import random
import numpy
import pathlib
cwd = pathlib.Path.cwd()
import pandas
from MomentCalculations import _unrenormalize
from MomentCalculations import moment
from Heatmaps import _import_data
# A collection of functions that are useful. I made these in a jupyter notebook and need them in other areas.
# Probably should just make a grab-bag of functions. That's probalby what I should call this, GrabBag.py


def psuedo_exon(intron: str, kmer: int = 120, get_length: bool = False, *args, **kwargs):
    '''
    This will grab a random length from the intron (the length will have an average of 120 bps based off a normal distribution) and from that data moments will be calculated.

    so: short_intron = intron[random_start: random_end]

    length = normal(average = 120, sigma ~ 25)
    random_start = randint(0, len(intron) - length)
    random_end = random_start + length

    See if these shortened ones are any different then the exons.

    if get_length is true, then kmer is used as the average for determining the length

    If it can't get a psuedo exon string, None is returned
    '''

    if get_length:
        kmer = gen_length(average = kmer, *args, **kwargs)

    if kmer > len(intron):
        kmer = len(intron)

    start = random.randint(0, kmer)
    end = start + kmer

    try:
        fake_exon = intron[start: end]
    except Exception as e:
        fake_exon = None

    return fake_exon



def gen_seq(kmer: int = 6):
    '''
    '''
    if kmer < 0:
        kmer = abs(kmer)
    if kmer == 0:
        kmer = 6

    seq = ""

    for _ in range(kmer):
        n = gen_nuc()
        seq = f"{n}{seq}"

    return seq



def gen_nuc():
    '''
    Generates a nucleotide. The reason it generates the nucleotide rather then the number (which would appear to be more convienent) is because
    this allows me to save the sequences, look at them later, and use the existing scripts easily.
    '''
    n = random.randint(1, 4)

    if n == 1:
        return "A"
    elif n == 2:
        return "G"
    elif n == 3:
        return "T"
    else:
        return "C"
    

def gen_length(average: int = 120, sigma: int = 25):
    '''
    Generates the length, using a gaussian distribution.

    By trail and error: exons -> 120 and 25, introns -> 1000, 200
    '''

    return int(numpy.random.normal(average, scale = sigma))


def gen_region():
    '''
    randomly chooses exon or intron
    '''
    r = random.randint(1, 2)

    return "exon" if r == 1 else "intron"


def heat_dataset(source_data: pathlib.Path or pandas.DataFrame, length: int = 12, *args, **kwargs) -> pandas.DataFrame:
    '''
    I have lost how I originally created the dataset for the heatmaps, so here's a function for doing that.

    This assumes the data is a pickle or a dataframe.

    Except never mind, the frame data will work nicely. I'm just going to use this to filter down to a nice length for the sequences
    '''
    if isinstance(source_data, pathlib.Path):
        source_data: pandas.DataFrame = pandas.read_pickle(source_data)


    source_data["Length"] = source_data["Seq"].apply(lambda x: len(x))
    source_data = source_data[source_data["Length"] >= length]
    source_data = source_data.reset_index() 

    return source_data


def moments_generic(file: pathlib.Path, 
                    ms: list, 
                    k: int = 1, 
                    unlog: bool = False, logy: bool = False, N_value: bool = False, 
                    *args, **kwargs) -> list:
    '''
    Just does one moment calulation at a time
    '''
    return_list = list()

    if isinstance(file, pathlib.Path):
        data = _import_data(file, just_import = True, *args, **kwargs)
        data = pandas.DataFrame(data)
    elif isinstance(file, pandas.DataFrame):
        data = file
    
    if unlog:
        data = _unrenormalize(data, 2*k, *args, **kwargs)
    
    for m in ms:

        if N_value:
            N = 4**(2*k)
            N = N**((1/m) - 1)
            print(f"\t\tm = {m}\tN = {N}")
        else:
            N = 1
            print(f"\t\tm = {m}")

        data_v = moment(data, m = m, unlog = False, N = N, *args, **kwargs)

        if logy:
            data_v = numpy.log2(data_v)

        return_list.append(data_v)

    return return_list

    # return me, mi, mn, uni, pd