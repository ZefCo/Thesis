# So new stragety: just get everything and keep it locally.

# For each chromosome, get the known genes
# Iterate through chr1 to chr22 + chrX + chrY

# https://api.genome.ucsc.edu/getData/track?genome=hg38;track=knownGene;chrom={chr#}

# This will yeild a large json: will want knownGene['knownGene'][iteration]

# Will keep
# chrom
# chromStart
# chromEnd
# name <- probably use this?
# strand
# 

import requests
import pandas
import pathlib

chroms = [f'chr{i}' for i in range(1, 23)]
chroms.append('chrX')
chroms.append('chrY')

# hg38 = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes', 'chromStarts', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames', 'type', 'geneName', 'geneName2', 'geneType', 'transcriptClass', 'source', 'transcriptType', 'tag', 'level', 'tier']
# hg19 = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'proteinID', 'alignID']
hg19 = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']

genome = 'hg19'
columns = hg19
folder = genome.capitalize()

known_genes = pandas.DataFrame(columns=columns) 

for chrom in chroms:

    chr_genes = pandas.DataFrame(columns=columns) 
    
    ucsc_url = f'https://api.genome.ucsc.edu/getData/track?genome={genome};track=ncbiRefSeqCurated;chrom={chrom}'

    response: dict = requests.get(ucsc_url).json()['ncbiRefSeqCurated']

    print(f'{chrom}: {len(response)}')

    for i, geneData in enumerate(response):
        chr_genes.loc[len(chr_genes.index) + 1] = geneData.values()

    chr_genes.to_csv(pathlib.Path.cwd() / "Data_Files" / "Gene_Files" / folder / f"{chrom}_{genome}_genes.csv")

    print(f"Finished {chrom}")

    known_genes = pandas.concat([known_genes, chr_genes], axis=0)


known_genes.to_csv(pathlib.Path.cwd() / "Data_Files" / "Gene_Files" / folder / "Known_Genes_hg19.csv")

# My parents are pathetic peices of shit and I wish I was aborted.
