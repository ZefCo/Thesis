import pandas
import re
import pathlib
import Permutations as perm
import GeneClass as Gene



def main():
    '''
    '''
    cwd = pathlib.Path.cwd()
    root = cwd.parent
    # nuc_perm = perm.nucleotide_permutations()
    # print(nuc_perm)
    with open(root / "Data_Files" / "Sequence_Files" / "HG19" / "Known_Genes_hg19_noSeq.csv") as known_file:
    # with open(root / "Data_Files" / "Sequence_Files" / "HG19" / "hg_test.csv") as known_file:
        known_genes: pandas.DataFrame = pandas.read_csv(known_file, sep = ',', header=0, index_col=0)

    headers = ["name2", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "cdsStartStat", "cdsEndStat", "exonFrames"]
    known_genes = known_genes[headers]

    known_genes["chrom"] = known_genes["chrom"].astype("category")

    training_genes = pandas.DataFrame()
    for c in known_genes["chrom"].cat.categories:
        subframe = known_genes[known_genes["chrom"] == c]
        subframe = subframe.sample(n = 200)

        training_genes = pandas.concat([training_genes, subframe])

    rows, cols = training_genes.shape

    pickle_frame = pandas.DataFrame()
    # gene_rows = list()

    for row in range(rows):
        row_of_interest = training_genes.iloc[row, :]

        gene_of_interest: Gene = Gene.Gene(row_of_interest["name2"],
                                           ncibname=row_of_interest["name"], 
                                           chrm=row_of_interest["chrom"], 
                                           strand=row_of_interest["strand"], 
                                           txStart=row_of_interest["txStart"],
                                           txEnd=row_of_interest["txEnd"],
                                           cdsStart=row_of_interest["cdsStart"],
                                           cdsEnd=row_of_interest["cdsEnd"],
                                           exonCount=row_of_interest["exonCount"],
                                           exonStarts=row_of_interest["exonStarts"],
                                           exonEnds=row_of_interest["exonEnds"],
                                           exonFrames=row_of_interest["exonFrames"])


        gene_of_interest.sequence_breakdown()
        # gene_of_interest.write_sequences(root / "Data_Files" / "NucComp")

        gene_row: pandas.Series = gene_of_interest.new_row()
        # print(type(gene_row))

        pickle_frame = pandas.concat([pickle_frame, gene_row.to_frame().T])
        # gene_rows.append(gene_row)

        # if row == 3:
        #     break

        # pickle_frame = pickle_frame.reset_index()
        print(pickle_frame.shape)
        # pickle_set = set(pickle_frame.columns)
        pickle_frame.to_pickle("TrainingGeneData_v2.pkl")  # Updating Pickle file every iteration.

    # test_frame: pandas.DataFrame = pandas.read_pickle(cwd / "FUH.pkl")
    # print(test_frame.shape)
    # test_set = set(test_frame.columns)

    # print(test_set - pickle_set)



if __name__ in '__main__':
    main()