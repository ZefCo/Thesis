
import requests
import pandas
import pathlib
cwd = pathlib.Path.cwd()

def main():
    '''
    '''
    # chroms = [f'chr{i}' for i in range(1, 23)]
    # chroms.append('chrX')
    # chroms.append('chrY')
    fly_chroms = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"]
    mouse_chroms = [f'chr{i}' for i in range(1, 20)]
    mouse_chroms.append('chrX')
    mouse_chroms.append('chrY')

    # # hg38 = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes', 'chromStarts', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames', 'type', 'geneName', 'geneName2', 'geneType', 'transcriptClass', 'source', 'transcriptType', 'tag', 'level', 'tier']
    # # hg19 = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'proteinID', 'alignID']
    hg19 = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']  # can just reuse this one for the dm6

    # genome = "dm6" #flies
    genome = "mm39" #mouse
    # genome = 'hg19' #humans
    # columns = hg19
    # folder = genome.capitalize()

    fly_folder = cwd.parent / "Data_Files" / "Fly"
    fly_folder.mkdir(parents = True, exist_ok = True)
    mouse_folder = cwd.parent / "Data_Files" / "Mouse"
    mouse_folder.mkdir(parents = True, exist_ok = True)

    track = "ncbiRefSeqCurated"
    get_data(genome, mouse_chroms, hg19, mouse_folder, track)


def get_data(genome: str, chroms: list, columns: list, folder_path: pathlib.Path, track: str):
    '''
    '''

    known_genes = pandas.DataFrame(columns=columns) 

    for chrom in chroms:

        chr_genes = pandas.DataFrame(columns=columns) 
        
        ucsc_url = f'https://api.genome.ucsc.edu/getData/track?genome={genome};track={track};chrom={chrom}'

        response: dict = requests.get(ucsc_url).json()['ncbiRefSeqCurated']

        print(f'{chrom}: {len(response)}')

        for i, geneData in enumerate(response):
            chr_genes.loc[len(chr_genes.index) + 1] = geneData.values()

        chr_genes.to_csv(folder_path / f"{chrom}_{genome}_genes.csv")

        print(f"Finished {chrom}")

        known_genes = pandas.concat([known_genes, chr_genes], axis=0)


    known_genes.to_csv(folder_path / f"Known_Genes_{genome}.csv")


if __name__ in "__main__":
    main()
