
import requests
import pandas
import pathlib
cwd = pathlib.Path.cwd()

def main():
    '''
    '''
    # # hg38 = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes', 'chromStarts', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames', 'type', 'geneName', 'geneName2', 'geneType', 'transcriptClass', 'source', 'transcriptType', 'tag', 'level', 'tier']
    # # hg19 = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'proteinID', 'alignID']
    hg19 = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']  # can just reuse this one for the dm6

    # genome = 'hg19' #humans
    # chroms = [f'chr{i}' for i in range(1, 23)]
    # chroms.append('chrX')
    # chroms.append('chrY')
    # folder = genome.capitalize()

    # genome = "dm6" #flies
    # fly_chroms = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"]
    # fly_folder = cwd.parent / "Data_Files" / "Fly"
    # fly_folder.mkdir(parents = True, exist_ok = True)

    # genome = "mm39" #mouse
    # mouse_chroms = [f'chr{i}' for i in range(1, 20)]
    # mouse_chroms.append('chrX')
    # mouse_chroms.append('chrY')
    # mouse_folder = cwd.parent / "Data_Files" / "Mouse"
    # mouse_folder.mkdir(parents = True, exist_ok = True)

    # genome = "xenTro9" # Western frog
    # frog_chroms = [f"chr{i}" for i in range(1, 11)]
    # frog_chroms.append("chrM")
    # frog_folder = cwd.parent / "Data_Files" / "Frog"
    # frog_folder.mkdir(parents=True, exist_ok=True)

    # genome = "rn7"
    # chroms = [f"chr{i}" for i in range(1, 21)]
    # chroms.append("chrX")
    # chroms.append("chrY")
    # chroms.append("chrM")
    # folder = cwd.parent / "Data_Files" / "Rat"
    # folder.mkdir(parents=True, exist_ok=True)

    # genome = "panTro6"
    # chroms = ["chr1", "chr2A", "chr2B", "chr3", "chr4", "chr5", "chrX", "chrY", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    # folder = cwd.parent / "Data_Files" / "Chimp"
    # folder.mkdir(parents=True, exist_ok=True)

    # genome = "danRer7"
    # chroms = [f"chr{i}" for i in range(1, 26)]
    # folder = cwd.parent / "Data_Files" / "Zebrafish"
    # folder.mkdir(parents=True, exist_ok=True)  #doesn't have a ref seq curated!

    genome = "sacCer3"
    chroms = ["chrIV", "chrM", "chrI", "chrVI", "chrIII", "chrIX", "chrVIII", "chrV", "chrXI", "chrX", "chrXIV", "chrII", "chrXIII", "chrXVI", "chrXII", "chrVII", "chrXV"]
    folder = cwd.parent / "Data_Files" / "Yeast"
    folder.mkdir(parents=True, exist_ok=True)

    # columns = hg19

    track = "ncbiRefSeqCurated"
    get_data(genome, chroms, hg19, folder, track)


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
