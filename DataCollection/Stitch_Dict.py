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
    # with open(pathlib.Path(f"/media/ethanspeakman/Elements/GeneData/Known_Genes_hg19_NCBIGene_DICT.pkl"), "rb") as p:
    #     data = pickle.load(p)

    # print(data["NM_001351428.2"])
    # print(data["NM_001351428.2"].full_seq)
    # test_dict(pathlib.Path("/media/ethanspeakman/Elements/GeneData/Known_Genes_hg19_NCBIGene_DICT_31298.pkl"), "NM_001329675.2")
    # stitch_dict("/media/ethanspeakman/Elements/GeneData", "/media/ethanspeakman/Elements/GeneData/Known_Genes_Master.pkl")
    select_data("/media/ethanspeakman/Elements/GeneData/Known_Genes_hg19_NCBIGene_DICT.pkl", "/media/ethanspeakman/Elements/GeneData/", n = 20)


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

    print(x)
    print(x[keys[0]])
    print(x[keys[0]].full_seq[0])

    if test_key is not None:
        check_key(x, test_key)


def histogram(dictionary: dict, report_file: pathlib.Path, min_length: int = 10):
    '''
    I was going to embed this in the select data method, but a) I don't like doing that and b) this is going to be too big to do that.

    This will output a pandas Dataframe that has the sequence (str), exon/intron classification (category/str), and the source gene (category/str)

    This also generates a report of the data that was put into it
    '''
    gene: Gene.Gene
    columns = ["NCIBName", "Classificaion", "Seq"]  #instead of manually typing these out for multiple lines, I'm making one thing here because it's going to be easier if I have to adjust later

    genetic: pandas.DataFrame = pandas.DataFrame(columns = columns)
    with open(report_file, "w+") as txtf:  # I want a report of what I am generating.
        for ncib, gene in dictionary.items():
            max_length = 0

            gene_data_line = f"{gene.name}\t{gene.ename}\t{gene.gname}\t{gene.ncibname}\nChrome: {gene.chrm}\tStrand: {gene.strand}\nReminder: Index here starts at 1\n\n"
            txtf.write(gene_data_line)
            exon_count, intron_count = 0, 0

            for e in gene.exon_seq:
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


            for i in gene.intron_seq:
                intron_count += 1
                intron_len = len(i)

                if (max_length >= intron_len) and (intron_len >= min_length):
                    intron_data_line = f"Intron {intron_count}: {intron_len}\n"
                    txtf.write(intron_data_line)

                    new_row = {f"{columns[0]}": [ncib], f"{columns[1]}": ["intron"], f"{columns[2]}": [e]}
                    new_row = pandas.DataFrame(new_row)
                    genetic = pandas.concat([genetic, new_row], axis = 0, ignore_index = True)
    
    
    return genetic







def select_data(file_path: pathlib.Path, output_path: pathlib.Path, n = 20_000):
    '''
    Randomly selects a certain ammount of the genes from the data. By default the amount if 20,000. This is then split into two even sub groups.
    '''

    def random_keys(dictionary: dict, half = False) -> dict:
        keys = tuple(dictionary.keys())
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

    with open(file_path, "rb") as p:
        data: dict = pickle.load(p)

    data = random_keys(data, half = False)

    data_1, data_2 = random_keys(data, half = True)

    data_1 = histogram(data_1, cwd / "Report_1.txt")
    data_2 = histogram(data_2, cwd / "Report_2.txt")

    with open(output_path / "Data_Set_1.pkl", "wb") as p:
        pickle.dump(data_1, p)

    with open(output_path / "Data_Set_2.pkl", "wb") as p:
        pickle.dump(data_2, p)


def stitch_dict(directory_path: pathlib.Path, output_file: pathlib.Path):
    '''
    '''
    if not isinstance(directory_path, pathlib.Path):
        directory_path = pathlib.Path(directory_path)
    if not isinstance(output_file, pathlib.Path):
        output_file = pathlib.Path(output_file)

    files = []
    z = {}

    for file in os.listdir(directory_path):
        if file.endswith(".pkl"):
            files.append(directory_path / file)
    
    for file in files:
        with open(file, "rb") as p:
            x = pickle.load(p)

        z = z | x


    with open(output_file, 'wb') as p:
        pickle.dump(z, p)



if __name__ in "__main__":
    main()