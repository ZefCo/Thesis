import pandas
import re
import pathlib
cwd = pathlib.Path.cwd()
# import Permutations as perm
import GeneClass as Gene
import pickle
import requests
import random
import os


def main():
    '''
    '''

    mouse = "Mouse"
    fly = "Fly"
    n = 10_000

    # species = "Mouse"
    # genome = "mm39"

    species = "Fly"
    genome = "dm6"

    # redo(f"D:\Downloads\GeneData\{species}\Master{species}Dict.pkl", genome, f"D:\Downloads\GeneData\{species}\Fixed_Master{species}Dict.pkl")
    # with open(pathlib.Path(f"D:\Downloads\GeneData\Fly/MasterflyDict.pkl"), "rb") as p:
    #     data = pickle.load(p)

    # print(type(data))
    # print(len(data))
    # first_key = tuple(data.keys())[0]
    # print(first_key)
    # print(data[first_key].exon_seq)
    # print(type(data[first_key]))

    # print(data["NM_001351428.2"])
    # print(data["NM_001351428.2"].full_seq)
    # stitch_frame("/media/ethanspeakman/Elements/Gene_Data_Sets/Combined/", "/media/ethanspeakman/Elements/Gene_Data_Sets/Combined/Combined_Hist.pkl")
    # select_data("G:/Known_Genes_Master.pkl", "G:/Gene_Data_Sets", n = 40_000)
    # test_dict(pathlib.Path("D:\Downloads\GeneData\Fly\SelectedFlyDict\Data_set_1_cleaned_dict.pkl"))
    # test_dict(pathlib.Path("G:\Gene_Data_Sets\Data_set_2_cleaned_dict.pkl"))

    # stitch_dict(f"D:\Downloads\GeneData\{fly}", f"D:\Downloads\GeneData\{fly}\Master{fly}Dict.pkl")
    # stitch_dict(f"D:\Downloads\GeneData\{mouse}", f"D:\Downloads\GeneData\{mouse}\Master{mouse}Dict.pkl")
    select_data(f"D:\Downloads\GeneData\{fly}\Fixed_Master{fly}Dict.pkl", f"D:\Downloads\GeneData\{fly}\Selected{fly}Dict", n  = n)
    select_data(f"D:\Downloads\GeneData\{mouse}\Fixed_Master{mouse}Dict.pkl", f"D:\Downloads\GeneData\{mouse}\Selected{mouse}Dict", n  = n)


def redo(file_path: pathlib.Path, genome: str, new_file_path: pathlib.Path):
    '''
    So when I first did this I didn't include the genome. The problem with that is the EVERYTHING was done via hg19. hg19 has no mouse or fly data, meaning everything was blank.

    So this just takes the dictionary that was already done and redoes all the stuff with the proper genome.
    '''
    name: str
    gene_data: Gene.Gene
    with open(file_path, "rb") as pkl_file:
        gene_dict: dict = pickle.load(pkl_file)

    new_gene_data = dict()

    for name, gene_data in gene_dict.items():
        gene_of_interest: Gene.Gene = Gene.Gene(name = gene_data.name,
                                           ncibname = gene_data.ncibname,
                                           ename = gene_data.ename,
                                           gname = gene_data.gname, 
                                           chrm = gene_data.chrm, 
                                           strand = gene_data.strand, 
                                           txStart = gene_data.txStart,
                                           txEnd = gene_data.txEnd,
                                           cdsStart = gene_data.cdsStart,
                                           cdsEnd = gene_data.cdsEnd,
                                           exonCount = gene_data.exonCount,
                                           exonStarts = gene_data.exonStarts,
                                           exonEnds = gene_data.exonEnds,
                                           exonFrames = gene_data.exonFrames,
                                           genome = genome)
        
        try:
            gene_of_interest.sequence_breakdown()
        except Exception as e:
            print(type(e))
            print(e)
            # exit()
            continue

        new_gene_data[name] = gene_of_interest

        try:
            with open(new_file_path, 'wb') as p:
                pickle.dump(new_gene_data, p)
        except Exception as e:
            print(type(e))
            print(f"Unable to write to file at this time: please check location can be written to")
            exit()



def test_dict(file_path: pathlib.Path, test_key = None):
    '''
    Opens and prints the contents of the dictioanry. Makes sure things are where they should be.
    '''
    def check_key(dic: dict, key: str):
        if key in dic.keys():
            print(f"{dic[key].full_seq[0]}")
        else:
            print(f"Key {key} not present")

    with open(file_path, "rb") as p:
        x: dict = pickle.load(p)

    keys = tuple(x.keys())

    print(f"Full sequence and object:")
    print(x[keys[0]].full_seq[0])
    print(x[keys[0]])
    # print(x)
    print(f"Total len of {file_path}: {len(keys)}")

    if test_key is not None:
        check_key(x, test_key)



def convert2dataframe(dictionary: dict, report_file: pathlib.Path, min_length: int = 10, columns = ["NCIBName", "Classificaion", "Seq"], histogram: bool = True):
    '''
    I was going to embed this in the select data method, but a) I don't like doing that and b) this is going to be too big to do that.

    This will output a pandas Dataframe that has the sequence (str), exon/intron classification (category/str), and the source gene (category/str)

    This also generates a report of the data that was put into it.

    You can pre-apply the histogram method if you want.
    '''
    gene: Gene.Gene

    genetic: pandas.DataFrame = pandas.DataFrame(columns = columns)
    with open(report_file, "w+") as txtf:  # I want a report of what I am generating.
        for ncib, gene in dictionary.items():
            max_length = 0

            if isinstance(gene.exon_seq, list) and isinstance(gene.intron_seq, list):

                gene_data_line = f"{gene.name}\t{gene.ename}\t{gene.gname}\t{gene.ncibname}\nChrome: {gene.chrm}\tStrand: {gene.strand}\nReminder: Index here starts at 1\n\n"
                txtf.write(gene_data_line)
                exon_count, intron_count = 0, 0

                for e in gene.exon_seq:
                    if isinstance(e, str):
                        exon_count += 1

                        exon_len = len(e)
                        if exon_len >= min_length:
                            exon_data_line = f"Exon {exon_count}: {exon_len}\n"
                            txtf.write(exon_data_line)

                            if exon_len >= max_length:
                                max_length = exon_len

                            new_row = {f"{columns[0]}": [ncib], f"{columns[1]}": ["exon"], f"{columns[2]}": [e]}
                            new_row = pandas.DataFrame(new_row)
                            genetic = pandas.concat([genetic, new_row], axis = 0, ignore_index = True)

                
                if histogram:
                    for i in gene.intron_seq:
                        if isinstance(i, str):
                            intron_count += 1
                            intron_len = len(i)

                            if (max_length >= intron_len) and (intron_len >= min_length):
                                intron_data_line = f"Intron {intron_count}: {intron_len}\n"
                                txtf.write(intron_data_line)

                                new_row = {f"{columns[0]}": [ncib], f"{columns[1]}": ["intron"], f"{columns[2]}": [e]}
                                new_row = pandas.DataFrame(new_row)
                                genetic = pandas.concat([genetic, new_row], axis = 0, ignore_index = True)
                else:
                    for i in gene.intron_seq:
                        if isinstance(i, str):
                            intron_count += 1
                            intron_len = len(i)

                            if intron_len >= min_length:
                                intron_data_line = f"Intron {intron_count}: {intron_len}\n"
                                txtf.write(intron_data_line)

                                new_row = {f"{columns[0]}": [ncib], f"{columns[1]}": ["intron"], f"{columns[2]}": [e]}
                                new_row = pandas.DataFrame(new_row)
                                genetic = pandas.concat([genetic, new_row], axis = 0, ignore_index = True)


    return genetic



def clean_dict(dictionary: dict, report_file: pathlib.Path, min_length: int = 10):
    '''
    Just in case something weird slips thorugh, like a None or something.
    '''
    gene: Gene.Gene
    clean_dictionary: dict = {}

    with open(report_file, "w+") as txtf:  # I want a report of what I am generating.
        for ncib, gene in dictionary.items():

            if isinstance(gene.exon_seq, list) and isinstance(gene.intron_seq, list):

                gene_data_line = f"{gene.name}\t{gene.ename}\t{gene.gname}\t{gene.ncibname}\nChrome: {gene.chrm}\tStrand: {gene.strand}\nReminder: Index here starts at 1\n\n"
                txtf.write(gene_data_line)
                exon_count, intron_count = 0, 0

                for e in gene.exon_seq:
                    if isinstance(e, str):
                        exon_count += 1

                        exon_len = len(e)
                        if exon_len >= min_length:
                            exon_data_line = f"Exon {exon_count}: {exon_len}\n"
                            txtf.write(exon_data_line)

                for i in gene.intron_seq:
                    if isinstance(i, str):
                        intron_count += 1
                        intron_len = len(i)

                        if intron_len >= min_length:
                            intron_data_line = f"Intron {intron_count}: {intron_len}\n"
                            txtf.write(intron_data_line)

                clean_dictionary[ncib] = gene
    
    return clean_dictionary





def select_data(file_path: pathlib.Path, output_path: pathlib.Path, n = 25_000, histogram_method: bool = True):
    '''
    Randomly selects a certain ammount of the genes from the data. By default the amount if 25,000. This is then split into two even sub groups. This is then output to 4 files: the two datasets, and two additional datasets
    filtered by the histogram method. The first two cleaned datasets are stored in dictioanry form - you'll have to go back and turn those into dataframes or some other type to be useful. The histogramed datasets are output
    into a dataframe. I did this for a reason: the histogram method was successfully used in the machine learning algorithms, and so became the "standard" we looked too. However, I began to realize using the non-histogram data
    might be usful, so that's why you get 4 files and not all of them useful right away.
    '''

    def random_keys(dictionary: dict, n: int, half = False) -> dict:
        keys = tuple(dictionary.keys())
        if n > len(keys):
            n = len(keys)
        random_keys = random.choices(keys, k = n)

        if half:
            random_keys = random.choices(random_keys, k = n)  # mixing up the order

            random_keys_1 = random_keys[0: int(n / 2)]
            random_keys_2 = random_keys[int(n / 2): len(random_keys) + 1]

            return {key: dictionary[key] for key in random_keys_1}, {key: dictionary[key] for key in random_keys_2}
        else:

            return {key: dictionary[key] for key in random_keys}
        

    if not isinstance(file_path, pathlib.Path):
        file_path = pathlib.Path(file_path)
    if not isinstance(output_path, pathlib.Path):
        output_path = pathlib.Path(output_path)

    print(f"Checking if {output_path} exists, and if not creating it.")
    output_path.mkdir(parents = True, exist_ok = True)

    print(f"Loading data from {file_path}")
    with open(file_path, "rb") as p:
        data: dict = pickle.load(p)

    print(f"Selecting {n} random keys")
    data = random_keys(data, n, half = False)

    print(f"Creating two seperate datasets from the previously selected data")
    data_1, data_2 = random_keys(data, n, half = True)

    print(f"Cleaning dataset 1")
    data_1 = clean_dict(data_1, output_path / "Report_1.txt")
    data_1_file = "Data_Set_1_cleaned_dict.pkl"
    print(f"Writing Dataset 1 to {output_path}/{data_1_file}")
    with open(output_path / data_1_file, "wb") as p:
        pickle.dump(data_1, p)

    print(f"Converting Dataset 1 to a dataframe")
    data_1_prime = convert2dataframe(data_1, output_path / "Report_1_cleaned.txt", histogram = False)
    data_1_file = "Data_Set_1_frame.pkl"
    print(f"Writing converted Dataset 1 to {output_path}/{data_1_file}")
    with open(output_path / data_1_file, "wb") as p:
        data_1_prime.to_pickle(p)

    print(f"Cleaning dataset 2")
    data_2 = clean_dict(data_2, output_path / "Report_2.txt")
    data_2_file = "Data_Set_2_cleaned_dict.pkl"
    print(f"Writing Dataset 2 to {output_path}/{data_2_file}")
    with open(output_path / data_2_file, "wb") as p:
        pickle.dump(data_2, p)

    print(f"Converting Dataset 2 to a dataframe")
    data_2_prime = convert2dataframe(data_2, output_path / "Report_2_cleaned.txt", histogram = False)
    data_2_file = "Data_Set_2_frame.pkl"
    print(f"Writing converted Dataset 2 to {output_path}/{data_2_file}")
    with open(output_path / data_2_file, "wb") as p:
        data_2_prime.to_pickle(p)

    if histogram_method:
        print(f"Selecting genetic data based on histogram method")
        print(f"Reports of selected data will be written to {output_path}/Report_1_histogram.txt and Report_2_histogram.txt")

        data_1 = convert2dataframe(data_1, output_path / "Report_1_histogram.txt", histogram = True)
        data_1_file = "Data_Set_1_histogram.pkl"
        print(f"Writing Dataset 1 to {output_path / data_1_file}")
        with open(output_path / data_1_file, "wb") as p:
            pickle.dump(data_1, p)

        data_2 = convert2dataframe(data_2, output_path / "Report_2_histogram.txt", histogram = True)
        data_2_file = "Data_Set_2_histogram.pkl"
        print(f"Writing Dataset 2 to {output_path / data_2_file}")
        with open(output_path / data_2_file, "wb") as p:
            pickle.dump(data_2, p)



def stitch_dict(directory_path: pathlib.Path, output_file: pathlib.Path):
    '''
    Combines multiple dictionaries into one.
    '''
    if not isinstance(directory_path, pathlib.Path):
        directory_path = pathlib.Path(directory_path)
    if not isinstance(output_file, pathlib.Path):
        output_file = pathlib.Path(output_file)
    
    print(f"Outputing to {output_file}: first checking if {output_file.parent} exists")
    output_file.parent.mkdir(parents = True, exist_ok = True)  # makes sure the directory is there

    files = []
    z = {}

    print(f"Finding .pkl files in {directory_path}")
    for file in os.listdir(directory_path):
        if file.endswith(".pkl"):
            files.append(directory_path / file)
    
    print(f"Found {len(files)} pikle files, now combining them")
    for file in files:
        with open(file, "rb") as p:
            x = pickle.load(p)

        z = z | x

    print(f"Writing master pikle file to {output_file}")
    with open(output_file, 'wb') as p:
        pickle.dump(z, p)




def stitch_frame(directory_path: pathlib.Path, output_file: pathlib.Path):
    '''
    '''
    if not isinstance(directory_path, pathlib.Path):
        directory_path = pathlib.Path(directory_path)
    if not isinstance(output_file, pathlib.Path):
        output_file = pathlib.Path(output_file)
    
    print(f"Outputing to {output_file}: first checking if {output_file.parent} exists")
    output_file.parent.mkdir(parents = True, exist_ok = True)  # makes sure the directory is there

    files = []
    z = pandas.DataFrame()

    print(f"Finding .pkl files in {directory_path}")
    for file in os.listdir(directory_path):
        if file.endswith(".pkl"):
            files.append(directory_path / file)
    
    print(f"Found {len(files)} pikle files, now combining them")
    for file in files:
        with open(file, "rb") as p:
            x = pickle.load(p)

        z = pandas.concat([z, x], ignore_index = True)

    z = z.reset_index()

    print(f"Writing master pikle file to {output_file}")
    with open(output_file, 'wb') as p:
        pickle.dump(z, p)



if __name__ in "__main__":
    main()