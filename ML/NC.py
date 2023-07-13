import pandas
import re
import pathlib
cwd = pathlib.Path.cwd()
import Permutations as perm
import GeneClass as Gene
import pickle



def single_image():
    '''
    '''
    random_samples = 1


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
        subframe = subframe.sample(n = random_samples)

        training_genes = pandas.concat([training_genes, subframe])
        break

    # print(training_genes)
    rows, _ = training_genes.shape

    for row in range(rows):
        print(f"Working on row {row}")
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

        # gene_row: pandas.Series = gene_of_interest.new_IE_row()

    # print(gene_of_interest.full_seq)
    # print(gene_of_interest.exon_seq[0])
    # print(gene_of_interest.intron_seq[0])
    # print(gene_row)

    # seq_frame = pandas.DataFrame([gene_of_interest.full_seq, gene_of_interest.exon_seq[0], gene_of_interest.intron_seq[0]], columns=["Full", "Exon", "Intron"])
    # seq_frame.to_csv()

    print(len(gene_of_interest.full_seq[0]))
    print(len(gene_of_interest.exon_seq[0]))
    print(len(gene_of_interest.intron_seq[0]))

    with open(cwd / "Seqs.txt", "w+") as textfile:
        textfile.write(gene_of_interest.full_seq[0].upper())
        textfile.write("\n\n\n")
        textfile.write(gene_of_interest.exon_seq[0].upper())
        textfile.write("\n\n\n")
        textfile.write(gene_of_interest.intron_seq[0].upper())


def main():
    '''
    '''
    random_samples = 200


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
        subframe = subframe.sample(n = random_samples)

        training_genes = pandas.concat([training_genes, subframe])

    rows, _ = training_genes.shape

    print(training_genes.shape)
    # print(pandas.Categorical(training_genes["chrom"]).unique())
    # exit()

    pickle_frame = pandas.DataFrame(columns=["Seq", "Type"])
    # gene_rows = list()
    type_sample = 20
    running_index = 0

    for row in range(rows):
        print(f"Working on row {row}")
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

        gene_row: pandas.Series = gene_of_interest.new_IE_row()

        row_index = gene_row.index

        for index in row_index:
            # print(index)
            if re.search("Seq", index):
                pass
            else:
                gene_row = gene_row.drop(index)

        # print(gene_row)
        # print(gene_row.shape)
        selected_samples = gene_row.shape[0]
        if selected_samples > type_sample:
            selected_samples = type_sample

        gene_row = gene_row.sample(n = selected_samples)
        # print(gene_row)

        row_index = gene_row.index
        for index in row_index:
            if re.search(r"Seq\d+", index):
                seq_type = re.sub(r"Seq\d+", "", index)

            # pickle_frame = pickle_frame.append({})
            pickle_frame.loc[running_index] = [gene_row[index], seq_type]
            running_index += 1

            # pickle_frame.append() = gene_row[index]
            # pickle_frame[running_index, "Type"] = seq_type

        # print(gene_row)
        # print(pickle_frame)

            

        # exit()
        # print(type(gene_row))

        # pickle_frame = pandas.concat([pickle_frame, gene_row.to_frame().T])
        # gene_rows.append(gene_row)

        # if row == 3:
        #     break

        # pickle_frame = pickle_frame.reset_index()
        # print(pickle_frame.shape)
        # pickle_set = set(pickle_frame.columns)
    pickle_frame.to_pickle("TrainingGeneData_v5.pkl")

    # test_frame: pandas.DataFrame = pandas.read_pickle(cwd / "FUH.pkl")
    # print(test_frame.shape)
    # test_set = set(test_frame.columns)

    # print(test_set - pickle_set)



def same_length_same_gene(threshold: float = 0.15):
    '''
    '''

    percent_diff = lambda x1, x2: (abs(x1 - x2)) / (0.5*(x1 + x2))

    random_samples = 200


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
        subframe = subframe.sample(n = random_samples)

        training_genes = pandas.concat([training_genes, subframe])

    rows, _ = training_genes.shape

    print(training_genes.shape)
    # print(pandas.Categorical(training_genes["chrom"]).unique())
    # exit()

    pickle_frame = pandas.DataFrame(columns=["Seq", "Type"])
    # gene_rows = list()
    type_sample = 20
    running_index = 0

    for row in range(rows):
        print(f"Working on row {row}")
        new_gene_row = pandas.Series()
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

        gene_row: pandas.Series = gene_of_interest.new_IE_row()


        row_index = gene_row.index
        exon_len_index = []
        intron_len_index = []

        for index in row_index:
    #         # print(index)
            if re.search("Seq", index):
                if re.search("Exon", index):
                    len_header = re.sub("ExonSeq", "ExonLen", index)
                    exon_len_index.append(len_header)

                elif re.search("Intron", index):
                    len_header = re.sub("IntronSeq", "IntronLen", index)
                    intron_len_index.append(len_header)

                gene_row[len_header] = len(gene_row[index])

            else:
                gene_row = gene_row.drop(index)

        # print(new_gene_row)
        # print(gene_row)

        for eength_index in exon_len_index:
            eength = gene_row[eength_index]

            if eength >= 100:

                local_intron_len_index = intron_len_index
                
                for i, ienght_index in enumerate(local_intron_len_index):
                    iength = gene_row[ienght_index]

                    if percent_diff(eength, iength) <= threshold:
                        exon_header = re.sub("ExonLen", "ExonSeq", eength_index)
                        intron_header = re.sub("IntronLen", "IntronSeq", ienght_index)
                        pickle_frame.loc[running_index] = [gene_row[exon_header], "Exon"]
                        running_index += 1

                        pickle_frame.loc[running_index] = [gene_row[intron_header], "Intron"]
                        running_index += 1

                        intron_len_index.pop(i)

                        break

    pickle_frame.to_pickle("TrainingGeneData_vSLSG.pkl")

def same_length_same_gene_his(threshold: float = 0.15):
    '''
    '''

    percent_diff = lambda x1, x2: (abs(x1 - x2)) / (0.5*(x1 + x2))

    random_samples = 200


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
        subframe = subframe.sample(n = random_samples)

        training_genes = pandas.concat([training_genes, subframe])

    rows, _ = training_genes.shape

    print(training_genes.shape)
    # print(pandas.Categorical(training_genes["chrom"]).unique())
    # exit()

    pickle_frame = pandas.DataFrame(columns=["Seq", "Type"])
    # gene_rows = list()
    type_sample = 20
    running_index = 0

    for row in range(rows):
        print(f"Working on row {row}")
        new_gene_row = pandas.Series()
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

        gene_row: pandas.Series = gene_of_interest.new_IE_row()


        row_index = gene_row.index

        min_length = 100
        max_length = min_length

        for index in row_index:
            if re.search("Seq", index):
                if re.search("Exon", index):
                    gene_len = len(gene_row[index])
                    if max_length < gene_len:
                        max_length = gene_len
            else:
                gene_row = gene_row.drop(index)


        for index in row_index:
            new_row = {"Seq": [], "Type": []}
    #         # print(index)
            if re.search("Seq", index):
                gene_length = len(gene_row[index])

                if (gene_length >= min_length) and (gene_length <= max_length):
                    if re.search(r"Seq\d+", index):
                        seq_type = re.sub(r"Seq\d+", "", index)

                        new_row["Seq"] = gene_row[index]
                        new_row["Type"] = seq_type

                        pickle_frame = pandas.concat([pickle_frame, pandas.DataFrame([new_row])], ignore_index = True)
                        # pickle_frame = pickle_frame.append()

                        # pickle_frame[running_index] = [gene_row[index], seq_type]
                        # running_index += 1
                        # print(pickle_frame)

        # if row > 1:
            # exit()

    pickle_frame.to_pickle("TrainingGeneData_SLSGHis.pkl")


def generate_gene_data(random_samples = 200):
    '''
    '''


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
        subframe = subframe.sample(n = random_samples)

        training_genes = pandas.concat([training_genes, subframe])

    rows, _ = training_genes.shape

    print(training_genes.shape)
    # print(pandas.Categorical(training_genes["chrom"]).unique())
    # exit()

    pickle_dict = dict()
    # gene_rows = list()
    # type_sample = 20
    running_index = 0

    for row in range(rows):
        print(f"Working on row {row}")
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

        pickle_dict[running_index] = gene_of_interest
        running_index += 1

    with open("GenePickle.pkl", "wb") as handle:
        pickle.dump(pickle_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)








if __name__ in '__main__':
    # main()
    # single_image()
    # same_length_same_gene()
    # same_length_same_gene_his()
    generate_gene_data()
