import requests
import pandas
import pathlib
import re


genome = 'hg19'
length = 12

output_folder = pathlib.Path.cwd() / 'Data_Files' / 'Sequence_Files' / 'HG19'

files = ['GeneTest']
# files = [f'chr{chrm}_hg19_genes' for chrm in range(1, 23)]
# files = files + ['chrX_hg19_genes', 'chrY_hg19_genes']
start_index = {}
# start_index = {'chr1_hg19_genes': 7930, 
#                'chrX_hg19_genes': 2803, 
#                'chrY_hg19_genes': 382, 
#                'chr2_hg19_genes': 5826,
#                'chr3_hg19_genes': 5343,
#                'chr4_hg19_genes': 3176,
#                'chr5_hg19_genes': 3812,
#                'chr6_hg19_genes': 4189,
#                'chr7_hg19_genes': 4126,
#                'chr8_hg19_genes': 3036,
#                'chr9_hg19_genes': 3279,
#                'chr10_hg19_genes': 3300} 
# files.reverse()


def print_frame(dataframe: pandas.DataFrame, output_path: pathlib.Path):
    '''
    This way if I don't want to create the csv file I can just adjust one line as opposed to two
    '''
    dataframe.to_csv(output_path, index = False)
    pass



for csvfile in files:
    file_path = pathlib.Path.cwd() / 'Data_Files' / 'Gene_Files' / 'HG19' / f'{csvfile}.csv'

    with open(file_path) as kghg19:
        known_genes_locations = pandas.read_csv(kghg19, sep=',', header=0)

    rows, _ = known_genes_locations.shape

    # link = f'http://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chr};start={start_point};end={end_point}'


    column_names = list(known_genes_locations.columns)
    gene_frame = pandas.DataFrame(columns=column_names)
    print_row_start = 0
    newindex = []
    # print(gene_frame)

    if csvfile in start_index.keys():
        print_row_start = start_index[csvfile]

    for row in range(print_row_start, rows):

        newindex.append(row)

        # need: ['chrom']   ['txStart'] & ['txEnd']     ['cdsStart'] & ['cdsEnd']   ['exonStarts'] & ['exonEnds']
        row_of_interest = known_genes_locations.iloc[row, :].copy()
        name = row_of_interest['name2']
        
        chrom = row_of_interest['chrom']

        link = f'http://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chrom};'

        exon_starts, exon_ends = re.split(',', row_of_interest['exonStarts'])[: len(re.split(',', row_of_interest['exonStarts'])) - 1], re.split(',', row_of_interest['exonEnds'])[: len(re.split(',', row_of_interest['exonEnds'])) - 1]
        
        # exonresponse = ''
        for i, exon_start in enumerate(exon_starts):
            fivePrimeTitle = f'exon5Prime{i + 1}'
            threePrimeTitle = f'exon3Prime{i + 1}_R'

            exon_start, exon_end = int(exon_start), int(exon_ends[i])

            exon_length = exon_end - exon_start

            if exon_length > 0:

                if (exon_length >= length):
                    fivePrime_hyperlink = f'{link}start={exon_start};end={exon_start + length}'
                    threePrime_hyperlink = f'{link}start={exon_end - length};end={exon_end}'
                    
                else:
                    fivePrime_hyperlink = f'{link}start={exon_start};end={exon_end}'
                    threePrime_hyperlink = f'{link}start={exon_start};end={exon_end}'
                    
                fivePrime_seq: str = requests.get(fivePrime_hyperlink).json()['dna']
                # print(f"{name}\n5' link = {fivePrime_hyperlink}")
                # print(f"{name}\n5' start = {fivePrime_seq}")
                threePrime_seq: str = requests.get(threePrime_hyperlink).json()['dna'][::-1]
                # print(f"3' link = {threePrime_hyperlink}")
                # print(f"3' end REVERESED = {threePrime_seq}")
                # print(f"3' end CORRECT = {threePrime_seq[::-1]}")

                row_of_interest[fivePrimeTitle] = fivePrime_seq
                row_of_interest[threePrimeTitle] = threePrime_seq

                if (len(fivePrime_seq) != length) or (len(threePrime_seq) != length):
                    print(f"Something weird at index {i + 1}")
                    print(f"Exon 5' Start Length = {len(fivePrime_seq)}\tExon 3' End Length = {len(threePrime_seq)}\tExon Length = {exon_length}")
                    print(f"Links used:\n\t{fivePrime_hyperlink}\n\t{threePrime_hyperlink}")

            else:
                row_of_interest[fivePrimeTitle] = 0
                row_of_interest[threePrimeTitle] = 0
        
        gene_frame = pandas.concat([gene_frame, row_of_interest.to_frame().T], axis = 0, ignore_index = True)

        if (((row % 100) == 0) and (row != 0) and (print_row_start != rows)):

            print_row_end = row
            if print_row_end > print_row_start:
                # print(print_row_start, print_row_end)
                # print(newindex)

                print(f"Completed {row} of {rows} for file {csvfile}: {round((row / rows) * 100, 2)} % complete")
                print(f"Printing Smaller frame of rows {print_row_start}:{print_row_end}")
                # print(gene_frame)

                # sub_gene_frame = gene_frame.iloc[row - 100:row, :].copy()
                # print(sub_gene_frame)

                print_frame(gene_frame, output_folder / 'SubFiles' / f'{csvfile}_{print_row_start}_{print_row_end}.csv')

                print_row_start = print_row_end

                gene_frame = pandas.DataFrame(columns=column_names)

        elif (row + 1 == rows):
            print(f"Ending file {csvfile}")
            try:
                print(f'starting position for final group: {print_row_start}\nending position for final group: {print_row_end}')
            except NameError:
                # print("print row end is not defined, need to set it here")
                # print(print_row_start)
                print_row_end = rows
                print(f'starting position for final group: {print_row_start}\tending position for final group: {print_row_end}')
            # print(gene_frame)

            # gene_frame.to_csv(output_folder / 'SubFiles' / f'{csvfile}_{print_row_start}_{print_row_end}_len{length}.csv', index = False)
            print_frame(gene_frame, output_folder / 'SubFiles' / f'{csvfile}_{print_row_start}_{print_row_end}_len{length}.csv')

            gene_frame = pandas.DataFrame(columns=column_names)

        # else:
        #     print("Still working")


        # if (row == rows) and (print_row_start > rows):


    # print(known_genes_locations)
    # print(gene_frame)

    # gene_frame.to_csv(output_folder / "SubFiles" / f'{csvfile}.csv', index = False)
    # print(f"Created file {csvfile} in folder {output_folder}")
    # # print(output_folder / f'{csvfile}.csv')