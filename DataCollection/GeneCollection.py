import pandas
import re
import pathlib
cwd = pathlib.Path.cwd()
# import Permutations as perm
import GeneClass as Gene
import pickle
import requests


def main():
    '''
    '''
    start_gene = 7828
    # dict_screwup()
    # pickle_file, csv_file = getKnownGene()
    hg19_sequences(cwd.parent / "Data_Files" / "Gene_Files" / "Hg19" / "Known_Genes_hg19_ncbiRefSeqCurated.pkl",
                   pathlib.Path(f"/media/ethanspeakman/Elements/GeneData/Known_Genes_hg19_NCBIGene_DICT_{start_gene}.pkl"), 
                   ref_track="ncib",
                   gene_start = start_gene)
    # hg19_sequences(cwd.parent / "Data_Files" / "Gene_Files" / "Hg19" / "Known_Genes_hg19_ncbiRefSeqCurated.pkl")


def dict_screwup():
    '''
    '''
    with open(pathlib.Path("/media/ethanspeakman/Elements/GeneData/Known_Genes_hg19_ensGene_DICT.pkl"), "rb") as p:
        something = pickle.load(p)

    print(something)
    exit()
    # external_drive = pathlib.Path("/media/ethanspeakman/Elements/GeneData")
    # print(external_drive)

    # test_dir = external_drive / "Test"

    # test_dir.mkdir(parents=True, exist_ok=True)
    # # external_drive.is_mount()

    # # d = dict()
    # # key = "key"
    # # value = 0
    # # version = 0
    # # for _ in range(10):
    # #     if key in d.keys():
    # #         key = re.split("_", key)[0]
    # #         key = f"{key}_{version}"
    # #         version += 1
    # #     d[key] = value
    # #     value += 1

    # # print(d)




def getKnownGene(genome = "hg19", track = "ncbiRefSeqCurated"):
    '''
    This will get all the known gene data and load it to a pkl file for later use.

    It will return the names of the files that it outputs, so I can be lazy and use one script to do two things. THAT'S HOW YOU DO IT BENDER!
    '''
    chroms = [f'chr{i}' for i in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrY')

    if genome in "hg38":
        colnames = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes', 'chromStarts', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames', 'type', 'geneName', 'geneName2', 'geneType', 'transcriptClass', 'source', 'transcriptType', 'tag', 'level', 'tier']

    else:
        colnames = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']

    columns = colnames
    folder = genome.capitalize()
    data_folder = cwd.parent / "Data_Files" / "Gene_Files" / folder
    data_folder.mkdir(parents = True, exist_ok = True)

    known_genes = pandas.DataFrame(columns=columns) 

    for chrom in chroms:

        chr_genes = pandas.DataFrame(columns=columns) 
        
        ucsc_url = f'https://api.genome.ucsc.edu/getData/track?genome={genome};track={track};chrom={chrom}'

        response: dict = requests.get(ucsc_url).json()[track]

        print(f'{chrom}: {len(response)}')

        for i, geneData in enumerate(response):
            chr_genes.loc[len(chr_genes.index) + 1] = geneData.values()

        chr_genes.to_csv(data_folder / f"{chrom}_{genome}_genes_{track}.csv")

        print(f"Finished {chrom}")

        known_genes = pandas.concat([known_genes, chr_genes], axis=0)

    known_genes.to_csv(data_folder / f"Known_Genes_{genome}_{track}.csv")
    known_genes.to_pickle(data_folder / f"Known_Genes_{genome}_{track}.pkl")

    return data_folder / f"Known_Genes_{genome}_{track}.pkl", data_folder / f"Known_Genes_{genome}_{track}.csv"




def hg19_sequences(gene_file: pathlib.Path, output_file: pathlib.Path, ref_track = "ncib", gene_start = 0):
    '''
    This will go to UCSC Genome Browser and grab all HG19 genes, then save them in a pickle file. They will not be chopped up, but instead will
    be preserved in all their glory for chopping up later.
    '''
    # percent_diff = lambda x1, x2: (abs(x1 - x2)) / (0.5*(x1 + x2))

    # random_samples = 200

    # print(gene_file.parent)
    # exit()
    # pickle_dict = {"1": 1, "2": 2}
    # with open(output_file, 'wb') as p:
    #     pickle.dump(pickle_dict, p)

    # exit()

    print(f"writing file to\n\t{output_file}\nAfter ever iteration")


    if gene_file.suffix in ".pkl":
        with open(gene_file, "rb") as p:
            known_genes: pandas.DataFrame = pickle.load(p)
    elif gene_file.suffix in ".csv":
        with open(gene_file) as known_file:
        # with open(root / "Data_Files" / "Sequence_Files" / "HG19" / "hg_test.csv") as known_file:
            known_genes: pandas.DataFrame = pandas.read_csv(known_file, sep = ',', header=0, index_col=0)

    headers = ['name','name2', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    known_genes = known_genes[headers]

    if ref_track in "enst":
        known_genes = known_genes.rename(columns={"name": "ename", "name2": "gname"})
        # headers = ["enstname", "engname", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "cdsStartStat", "cdsEndStat", "exonFrames"]
        known_genes["name"] = None
        known_genes["ncibname"] = None
    else:
        known_genes = known_genes.rename(columns={"name": "ncibname", "name2": "name"})
        # headers = ["name", "ncibname", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "cdsStartStat", "cdsEndStat", "exonFrames"]
        known_genes["ename"] = None
        known_genes["gname"] = None

    # print(known_genes)
    # print(known_genes.columns)
    # exit()

    known_genes["chrom"] = known_genes["chrom"].astype("category")

    known_genes["exonCount"] = known_genes["exonCount"].fillna(0)

    known_genes = known_genes[known_genes["cdsStartStat"] == "cmpl"]
    known_genes = known_genes[known_genes["exonCount"] >= 1]

    # print(known_genes)
    # exit()

    rows, _ = known_genes.shape

    pickle_dict = dict()
    # gene_rows = list()
    # type_sample = 20
    # running_index = 0

    unique_index = 0
    for row in range(gene_start, rows):
        print(f"Working on row {row} / {rows}")
        row_of_interest = known_genes.iloc[row, :]

        gene_of_interest: Gene = Gene.Gene(name = row_of_interest["name"],
                                           ncibname = row_of_interest["ncibname"],
                                           ename = row_of_interest["ename"],
                                           gname = row_of_interest["gname"], 
                                           chrm = row_of_interest["chrom"], 
                                           strand = row_of_interest["strand"], 
                                           txStart = row_of_interest["txStart"],
                                           txEnd = row_of_interest["txEnd"],
                                           cdsStart = row_of_interest["cdsStart"],
                                           cdsEnd = row_of_interest["cdsEnd"],
                                           exonCount = row_of_interest["exonCount"],
                                           exonStarts = row_of_interest["exonStarts"],
                                           exonEnds = row_of_interest["exonEnds"],
                                           exonFrames = row_of_interest["exonFrames"])
        
        gene_of_interest.sequence_breakdown()

        if gene_of_interest.ncibname in pickle_dict.keys():
            key = re.split("_", gene_of_interest.ncibname)[0]
            key = f"{key}_{unique_index}"
            unique_index += 1
        else:
            key = gene_of_interest.ncibname

        pickle_dict[key] = gene_of_interest

        # if row > 10:
        #     break
        try:
            with open(output_file, 'wb') as p:
                pickle.dump(pickle_dict, p)
        except Exception as e:
            print(type(e))
            print(f"Unable to write to file at this time: please check location can be written to")
            
    print(f"Wrote file to\n\t{output_file}")


if __name__ in '__main__':
    main()
