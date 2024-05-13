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
    Might need to look at this a bit more: at gene 31298 I start getting Type Errors and it never clears up. 31297 is plent of genes to start
    with, but why is it coming up with an error now?

    But as of 8/19/23 31298 seems to be fine when doing this at home. Now I get an error at 37188...

    I keep getting various HTTP errors... have I hit a maximum number of times I can access this data? There can't be a hard limit on it, so what is going on?

    OK So I'm getting connection erros and I'm not totally sure how to handle all of them. Basically I'm getting kicked off the UCSC Genome Browser, and I think it's either because
    I'm pinging them too much so they kick me off, or I have weak internet and I get kicked off. Either way, I keep getting these inturruptions. One thing I could do: if I loose connection,
    just skip that intron/exon.

    Stopped at 40587 at school. Stopped at a different number at home. Both from a json decode error.
    '''
    # start_gene = 10872  # left off at 10593 for the sense strands
    # genome = "hg19"
    # species = "Homo_sapiens"

    # start_gene = 0
    # genome = "mm39"
    # species = "Mouse"

    # start_gene = 0
    # genome = "xenTro9"
    # species = "Frog"

    start_gene = 0
    genome = "dm6"
    species = "Fly"

    # start_gene = 0
    # genome = "rn7"
    # species = "Rat"

    # start_gene = 0
    # genome = "sacCer3"
    # species = "Yeast"

    # start_gene = 0
    # genome = "calJac4"
    # species = "Callithrix_jacchus"

    # start_gene = 0
    # genome = "gorGor6"
    # species = "Gorilla_gorilla_gorilla"

    # start_gene = 0
    # genome = "macFas5"
    # species = "Macaca_fascicularis"

    # start_gene = 0
    # genome = "rheMac8"
    # species = "Macaca_mulatta"

    # start_gene = 0
    # genome = "nomLeu3"
    # species = "Nomascus_leucogenys"

    # start_gene = 0
    # genome = "panPan3"
    # species = "Pan_paniscus"

    # start_gene = 0
    # genome = "panTro6"
    # species = "Pan_troglodytes"

    # start_gene = 0
    # genome = "papAnu4"
    # species = "Papio_anubis"

    # start_gene = 0
    # genome = "ponAbe3"
    # species = "Pongo_pygmaeus_abelii"

    # output_file = cwd.parent / "Data_Files" / "Primates" / "Genetics" / f"{species}" / f"Antisense_Investigation.pkl"
    # # dict_screwup()
    # # pickle_file, csv_file = getKnownGene()
    # # hg19_sequences(cwd.parent / "Data_Files" / "Gene_Files" / "Hg19" / "Known_Genes_hg19.csv",
    # #                pathlib.Path(f"D:/Downloads/GeneData/Known_Genes_hg19_NCBIGene_DICT_{start_gene}.pkl"), 
    # #                ref_track="ncib",
    # #                gene_start = start_gene)
    # # hg19_sequences(cwd.parent / "Data_Files" / "Gene_Files" / "Hg19" / "Known_Genes_hg19_ncbiRefSeqCurated.pkl")

    # hg19_sequences(cwd.parent / "Data_Files" / "Primates" / "Known_Genes" / species / f"Known_Genes_{genome}.csv",
    #                output_file, 
    #                ref_track="ncib",
    #                gene_start = start_gene,
    #                species = species,
    #                genome = genome)
 
    # # For primates
    # output_file = cwd.parent / "Data_Files" / "Primates" / "Genetics" / f"{species}" / f"Sense_Investigation.pkl"  #this is actually a spelling error but I'm leaving it for consistency
    # # For other species
    intput_file = cwd.parent / "Data_Files" / species / f"Known_Genes_{genome}.csv"
    output_file = cwd.parent / "Data_Files" / species / f"Sense_Investigation.pkl"

    sense_sequences(intput_file,
                   output_file, 
                   ref_track="ncib",
                   gene_start = start_gene,
                   species = species,
                   genome = genome)


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




def getKnownGene(genome = "hg19", track = "ncbiRefSeqCurated", chroms: list = None):
    '''
    This will get all the known gene data and load it to a pkl file for later use.

    It will return the names of the files that it outputs, so I can be lazy and use one script to do two things. THAT'S HOW YOU DO IT BENDER!
    '''
    if chroms is None:
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




def hg19_sequences(gene_file: pathlib.Path, output_file: pathlib.Path, ref_track = "ncib", gene_start = 0, gene_stop: int = None, species: str = "human", genome: str = "hg19"):
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
    folder_target = output_file.parent
    folder_target.mkdir(parents = True, exist_ok = True)
    
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
    known_genes = known_genes[known_genes["cdsEndStat"] == "cmpl"]
    known_genes = known_genes[known_genes["exonCount"] >= 1]

    # print(known_genes)
    # exit()

    if isinstance(gene_stop, int):
        rows = gene_stop
    else:
        rows, _ = known_genes.shape


    pickle_dict = dict()
    # gene_rows = list()
    # type_sample = 20
    # running_index = 0

    unique_index = 0
    for row in range(gene_start, rows):
        print(f"Working on row {row} / {rows} - {species}")
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
                                        exonFrames = row_of_interest["exonFrames"],
                                        genome = genome)
        try:
            gene_of_interest.sequence_breakdown()
        except Exception as e:
            print(type(e))
            print(e)
            # exit()
            continue

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
            exit()
                
        print(f"Wrote file to\n\t{output_file}")


def sense_sequences(gene_file: pathlib.Path, output_file: pathlib.Path, ref_track = "ncib", gene_start = 0, gene_stop: int = None, species: str = "human", genome: str = "hg19"):
    '''
    Just does the + sequences
    This will go to UCSC Genome Browser and grab all HG19 genes, then save them in a pickle file. They will not be chopped up, but instead will
    be preserved in all their glory for chopping up later.
    '''
    print(f"writing file to\n\t{output_file}\nAfter ever iteration")
    folder_target = output_file.parent
    folder_target.mkdir(parents = True, exist_ok = True)
    
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

    known_genes["chrom"] = known_genes["chrom"].astype("category")

    known_genes["exonCount"] = known_genes["exonCount"].fillna(0)

    known_genes = known_genes[known_genes["cdsStartStat"] == "cmpl"]
    known_genes = known_genes[known_genes["cdsEndStat"] == "cmpl"]
    known_genes = known_genes[known_genes["exonCount"] >= 1]

    if isinstance(gene_stop, int):
        rows = gene_stop
    else:
        rows, _ = known_genes.shape

    pickle_dict = dict()

    unique_index = 0
    for row in range(gene_start, rows):
        print(f"Working on row {row} / {rows} - {species}")
        row_of_interest = known_genes.iloc[row, :]
        if row_of_interest["strand"] in "+":

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
                                            exonFrames = row_of_interest["exonFrames"],
                                            genome = genome)
            try:
                gene_of_interest.sequence_breakdown()
            except Exception as e:
                print(type(e))
                print(e)
                continue

            if gene_of_interest.ncibname in pickle_dict.keys():
                key = re.split("_", gene_of_interest.ncibname)[0]
                key = f"{key}_{unique_index}"
                unique_index += 1
            else:
                key = gene_of_interest.ncibname

            pickle_dict[key] = gene_of_interest

            try:
                with open(output_file, 'wb') as p:
                    pickle.dump(pickle_dict, p)
            except Exception as e:
                print(type(e))
                print(f"Unable to write to file at this time: please check location can be written to")
                exit()
                
        print(f"Wrote file to\n\t{output_file}")



def anti_sequences(gene_file: pathlib.Path, output_file: pathlib.Path, ref_track = "ncib", gene_start = 0, gene_stop: int = None, species: str = "human", genome: str = "hg19"):
    '''
    Only does the Anit-sense sequences
    This will go to UCSC Genome Browser and grab all HG19 genes, then save them in a pickle file. They will not be chopped up, but instead will
    be preserved in all their glory for chopping up later.
    '''
    print(f"writing file to\n\t{output_file}\nAfter ever iteration")
    folder_target = output_file.parent
    folder_target.mkdir(parents = True, exist_ok = True)
    
    if gene_file.suffix in ".pkl":
        with open(gene_file, "rb") as p:
            known_genes: pandas.DataFrame = pickle.load(p)
    elif gene_file.suffix in ".csv":
        with open(gene_file) as known_file:
            known_genes: pandas.DataFrame = pandas.read_csv(known_file, sep = ',', header=0, index_col=0)

    headers = ['name','name2', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    known_genes = known_genes[headers]

    if ref_track in "enst":
        known_genes = known_genes.rename(columns={"name": "ename", "name2": "gname"})
        known_genes["name"] = None
        known_genes["ncibname"] = None
    else:
        known_genes = known_genes.rename(columns={"name": "ncibname", "name2": "name"})
        known_genes["ename"] = None
        known_genes["gname"] = None

    known_genes["chrom"] = known_genes["chrom"].astype("category")

    known_genes["exonCount"] = known_genes["exonCount"].fillna(0)

    known_genes = known_genes[known_genes["cdsStartStat"] == "cmpl"]
    known_genes = known_genes[known_genes["cdsEndStat"] == "cmpl"]
    known_genes = known_genes[known_genes["exonCount"] >= 1]

    if isinstance(gene_stop, int):
        rows = gene_stop
    else:
        rows, _ = known_genes.shape

    pickle_dict = dict()

    unique_index = 0
    for row in range(gene_start, rows):
        print(f"Working on row {row} / {rows} - {species}")
        row_of_interest = known_genes.iloc[row, :]
        if row_of_interest["strand"] in "-":

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
                                            exonFrames = row_of_interest["exonFrames"],
                                            genome = genome)
            try:
                gene_of_interest.sequence_breakdown()
            except Exception as e:
                print(type(e))
                print(e)
                continue

            if gene_of_interest.ncibname in pickle_dict.keys():
                key = re.split("_", gene_of_interest.ncibname)[0]
                key = f"{key}_{unique_index}"
                unique_index += 1
            else:
                key = gene_of_interest.ncibname

            pickle_dict[key] = gene_of_interest

            try:
                with open(output_file, 'wb') as p:
                    pickle.dump(pickle_dict, p)
            except Exception as e:
                print(type(e))
                print(f"Unable to write to file at this time: please check location can be written to")
                exit()
                
        print(f"Wrote file to\n\t{output_file}")



if __name__ in '__main__':
    main()
