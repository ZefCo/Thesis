import pathlib
cwd = pathlib.Path.cwd()
import pickle
import pandas
import ucsc_restapi as API

def main():
    '''
    This looks at and grabs the intergenomic data.

    This is going to be a bit more complicated then the genetic data.

    It will grab about 50 genes from each chromosome. Then it looks at the end of that gene and finds the next closest gene. Hopefully the genes are in positional order, but it does it the long
    way because I'm afriad they might not be.
    '''
    output_file = cwd.parent / "Data_Files" / "Primates" / "Homo_Sapiens" / "Inter" / "InterData_NCIB.pkl"
    # known_gene = cwd.parent.parent / "Thesis_Data" / "Known_Genes_hg19.csv"
    # getIntergenomic(known_gene, output_file)
    # with open(output_file, "rb") as file:
    #     data: pandas.DataFrame = pickle.load(file)

    # data = data.sample(n = 250)
    # data = data.reset_index()
    # rows, cols = data.shape
    # print(data)

    # for row in range(rows):
    #     row_of_interest = data.iloc[row, :]

    select_data(output_file)



def select_data(input_file: pathlib.Path, n: int = 250):
    '''
    This part selects the data
    '''
    with open(input_file, "rb") as file:
        data: pandas.DataFrame = pickle.load(file)
        # print(data)
        # exit()

    subdata = data.sample(n = n)
    subindex = tuple(subdata.index)

    file_path = cwd.parent / "Data_Files" / "Intergenomic" / "Seq2Check"
    right_file = file_path / "Right"
    right_file.mkdir(parents = True, exist_ok = True)
    left_file = file_path / "Left"
    left_file.mkdir(parents = True, exist_ok = True)
    

    # Goes through the lines and grabs the minimum positive distance between txS and txE: txS - txE
    # That becomes the new start and end:
    # txE = txS_inter
    # txS = txE_inter
    # and it will be named the genes that it sits between
    for index in subindex:
        txStart2 = subdata.loc[index, "txStart"]
        txEnd2 = subdata.loc[index, "txEnd"]
        chrom = subdata.loc[index, "chrom"]

        dataN: pandas.DataFrame = data.loc[data.index != index, :]
        dataN = dataN[dataN["chrom"] == chrom]
        center_name = subdata.loc[index, "ncibname"]

        dataN["delta_Left"] = txStart2 - dataN["txEnd"]  # looking for the left sided intergenomic data
        dataN["delta_Right"] = dataN["txStart"] - txEnd2  # looking for the right sided data

        # print(dataN)
        try:
            closest_start = min(filter(lambda x: x > 0, dataN["delta_Left"]))  # this looks for the start position of the left intergenomic region
            closest_end = min(filter(lambda x: x > 0, dataN["delta_Right"]))  # this look for the end position of the right intergenomic region
        except Exception as e:
            continue

        left_start_index = dataN.isin([closest_start]).any(axis=1).idxmax()
        right_end_index = dataN.isin([closest_end]).any(axis=1).idxmax()

        left_start = dataN.loc[left_start_index, "txEnd"]
        left_name = dataN.loc[left_start_index, "ncibname"]
        right_end = dataN.loc[right_end_index, "txStart"]
        right_name = dataN.loc[right_end_index, "ncibname"]

        print("\tLeft Query")
        left_query, left_url = API.sequence(chrom = chrom, start = left_start, end = txStart2, strand = "+")
        print("\tWriting Left Query")
        with open(left_file / f"{left_name}_{center_name}", "w+") as file:
            file.write(f"{left_name}_{center_name}\n{chrom}\n\n{left_url}\n\n\n{left_query}")

        print("\t\tRight Query")
        right_query, right_url = API.sequence(chrom = chrom, start = txEnd2, end = right_end, strand = "+")
        print("\t\tWriting Right Query")
        with open(right_file / f"{center_name}_{right_name}", "w+") as file:
            file.write(f"{center_name}_{right_name}\n{chrom}\n\n{right_url}\n\n\n{right_query}")

        # exit()

        # need to get the index of these values.

        


def write2file(input_file: pathlib.Path):
    '''
    Writes each intergenomic sequence to a file, because it will be easier for me to check that way.
    '''
    with open(input_file, "rb") as file:
        data: pandas.DataFrame = pickle.load(file)

    data = data.sample(n = 250)
    data = data.reset_index()
    rows, cols = data.shape

    for row in range(rows):
        chrom = data.loc[row, "chrom"]
        txStart = data.loc[row, "txStart"]
        txEnd = data.loc[row, "txEnd"]
        seq = API.sequence()


def getIntergenomic(gene_file: pathlib.Path, output_file: pathlib.Path, ref_track = "ncib", gene_start: int = 0, species: str = "Homo Sapiens",
                    *args, **kwargs):
    '''
    I'm just going to hard code a lot of the information here for now. Getting intergenomic data from other species for the most part doesn't make sense.
    '''

    inter_positions = pandas.DataFrame(columns = ["name", "txStart", "txEnd", "chrom", "strand"])


    print(f"writing file to\n\t{output_file}")
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
    chromes = tuple(known_genes["chrom"].unique())

    known_genes["exonCount"] = known_genes["exonCount"].fillna(0)

    # I'm going to be a bit more agressive about what to include: I want to make sure I'm not including *any* genes, so I'm going to widen what my definition of a gene is
    # known_genes = known_genes[known_genes["cdsStartStat"] == "cmpl"]
    # known_genes = known_genes[known_genes["cdsEndStat"] == "cmpl"]
    # known_genes = known_genes[known_genes["exonCount"] >= 1]

    new_row = 0
    errors = 0
    for c, chrom in enumerate(chromes):
        chrom_genes = known_genes[known_genes["chrom"] == chrom]
        unique_names = tuple(chrom_genes["ncibname"].unique())

        for n, name in enumerate(unique_names):
            name_genes: pandas.DataFrame = chrom_genes[chrom_genes["ncibname"] == name]
            try:
                txStart = min(tuple(name_genes["txStart"].to_list()))
            except Exception as e:
                print(name)
                print(name_genes["txStart"])
                print("\n\n")
                print(name_genes)
                print("\n\n")
                print(chrom_genes)
                exit()
                # print(type(e))
                # print(e)
                errors += 1
                continue

            try:
                txEnd = max(tuple(name_genes["txEnd"].to_list()))
            except Exception as e:
                print(name)
                print(name_genes["txEnd"])
                print("\n\n")
                print(name_genes)
                print("\n\n")
                print(chrom_genes)
                exit()
                # print(type(e))
                # print(e)
                errors += 1
                continue

            # if isinstance(name_genes, pandas.Series):
            #     print("Series")
            #     print(name_genes)
            # else:
            #     print("Dataframe")
            #     print(name_genes)
            #     print(name_genes["strand"])
            #     print(type(name_genes["strand"]))
            #     exit()
                
            strand = tuple(name_genes["strand"].to_list())
            if len(strand) > 0:
                for s in strand:
                    if strand[0] not in s:
                        continue
                else:
                    strand = strand[0]
            else:
                continue

            inter_positions.loc[new_row, "ncibname"] = name
            inter_positions.loc[new_row, "txStart"] = txStart
            inter_positions.loc[new_row, "txEnd"] = txEnd
            inter_positions.loc[new_row, "chrom"] = chrom
            inter_positions.loc[new_row, "strand"] = strand
            new_row += 1
            
    with open(output_file, 'wb') as p:
        pickle.dump(inter_positions, p)
            
    print(f"Wrote file to\n\t{output_file}")
    print(f"There were {errors} when getting this data, and I'm to lazy to figure out specifically why")




if __name__ in "__main__":
    main()