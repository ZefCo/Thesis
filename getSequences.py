# # example
# # http://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrM;start=4321;end=5678
# genome = 'hg38'
# chr = '2'
# start_point = '1'
# end_point = '1'


# # First I need the start and end position of the coding sequences
# # Also probably want the full sequence at some point
# ucsc_root = 'https://api.genome.ucsc.edu/'
# end_point_track = '/getData/track'
# end_point_seq = '/getData/sequence'

# # https://api.genome.ucsc.edu/getData/tracks?genome=hg38;track=geneid

# Need the gene name, the chromosome, the start, and the end position

# Import data to dataframe
# For each row in frame
    # make link to get sequences
    # get tx
    # get cds
    # get exons (seqeratly)
# Save to csv files

import requests
import pandas
import pathlib
import re


genome = 'hg19'

with open(pathlib.Path.cwd() / 'Data_Files' / 'Gene_Files' / 'HG19' / 'Known_Genes_hg19.csv') as kghg19:
    known_genes_locations = pandas.read_csv(kghg19, sep=',', header=0)

rows, columns = known_genes_locations.shape

# link = f'http://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chr};start={start_point};end={end_point}'

output_folder = pathlib.Path.cwd() / 'Data_Files' / 'Sequence_Files' / 'HG19'

txseq, cdsseq = 'TXSeq', 'CDSSeq'

column_names = list(known_genes_locations.columns) + [txseq, cdsseq]
return_frame = pandas.DataFrame(columns=column_names)
# print(return_frame)

for row in range(rows):
    # need: ['chrom']   ['txStart'] & ['txEnd']     ['cdsStart'] & ['cdsEnd']   ['exonStarts'] & ['exonEnds']
    row_of_interest = known_genes_locations.iloc[row, :]
    
    chrom = row_of_interest['chrom']

    link = f'http://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chrom};'

    txStart, txEnd = row_of_interest['txStart'], row_of_interest['txEnd']
    cdsStart, cdsEnd = row_of_interest['cdsStart'], row_of_interest['cdsEnd']
    exon_starts, exon_ends = re.split(',', row_of_interest['exonStarts'])[: len(re.split(',', row_of_interest['exonStarts'])) - 1], re.split(',', row_of_interest['exonEnds'])[: len(re.split(',', row_of_interest['exonEnds'])) - 1]

    if txEnd > txStart:
        txlink = f'{link}start={txStart};end={txEnd}'
        txresponse: str = requests.get(txlink).json()['dna']
        print(f'Length of TXResponse: {len(txresponse)}')

        row_of_interest.loc[txseq] = txresponse
    else:
        row_of_interest.loc[txseq] = 0

    if cdsEnd > cdsStart:
        cdslink = f'{link}start={cdsStart};end={cdsEnd}'
        cdsresponse: str = requests.get(cdslink).json()['dna']
        print(f'Length of TXResponse: {len(txresponse)}')

        row_of_interest[cdsseq] = cdsresponse
    else:
        row_of_interest[cdsseq] = 0

    
    exonresponse = ''
    for i, exon in enumerate(exon_starts):
        colname = f'exon{i + 1}'

        if exon_ends[i] > exon:
            exonlink = f'{link}start={exon};end={exon_ends[i]}'
            localexon: str = requests.get(exonlink).json()['dna']
            row_of_interest[colname] = localexon
            # print(exonresponse)
            # exonresponse = f'{exonresponse}{localexon}'

    # if len(exonresponse) > 0:
    #     row_of_interest[exonseq] = exonresponse
    #     print(f'Length of TXResponse: {len(txresponse)}')

        else:
            row_of_interest[colname] = 0

    print(row_of_interest)


    # return_frame.loc[len(return_frame.index)] = row_of_interest
    return_frame = pandas.concat([return_frame, row_of_interest.to_frame().T], axis = 0, ignore_index = True)

# print(known_genes_locations)
print(return_frame)

return_frame.to_csv(output_folder / 'KnowGeneSeq_hg19.csv', index = False)
