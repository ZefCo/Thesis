import re
import ucsc_restapi as upi
import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import pickle
import random
import GeneClass as Gene



def main():
    '''
    '''
    # data_files = cwd.parent / "Data_Files"
    # ut_file = data_files / "UTData_fs.pkl"
    # ut_data: pandas.DataFrame = pandas.read_pickle(ut_file)
    # ut_data = ut_data.rename(columns={"Unnamed: 18": "Posterior_20", "Unnamed: 19": "Anterior_20"})
    # ut_data.to_pickle(ut_file)

    # folder_path = pathlib.Path("F:\Gene_Data_Sets")
    # data: pathlib.Path = folder_path / "Exon2Exon.pkl"
    # data = pandas.read_pickle(data)
    # # save_file_path = folder_path / "E2E_DS1_k6.pkl"
    # print(data)
    # master_data: pandas.DataFrame = gen_te_points(data = data, n = 1_000)
    # folder_path = pathlib.Path("F:\Gene_Data_Sets")
    # e2e_file_path = folder_path / "E2E_DS1_k6.pkl"

    # data = pandas.read_pickle(e2e_file_path)
    # print(data)

    # CancerData()
    min_length: int = 20
    NormalData(min_length = min_length)


def CancerData():
    '''
    '''
    data_files = cwd.parent / "Data_Files"

    lab_file = data_files / "Fusion_AllConfidence.pkl"
    lab_file_out = data_files / "Lab_Cancer.xlsx"

    ut_file = data_files / "UTData_fs.pkl"
    ut_file_out = data_files / "UT_Cancer.xlsx"

    lab_cancer = pandas.DataFrame(columns = ["Name_A", "Chrm_A", "Name_P", "Chrm_P", "Posterior_10", "Anterior_10"])

    lab_data: pandas.DataFrame = pandas.read_pickle(lab_file)
    lab_data["Posterior_20"] = lab_data["Head"].apply(lambda x: x[len(x) - 20:len(x)])
    lab_data["Posterior_20"] = lab_data["Posterior_20"].apply(lambda x: pandas.NA if re.search(r"\_|\.", x) else x)
    lab_data["Anterior_20"] = lab_data["Tail"].apply(lambda x: x[0: 20])
    lab_data["Anterior_20"] = lab_data["Anterior_20"].apply(lambda x: pandas.NA if re.search(r"\_|\.", x) else x)

    lab_data = lab_data.dropna()
    lab_data = lab_data.reset_index()

    # print(lab_data.head())
    # print(lab_data.columns)

    lab_data = lab_data[["gene_id1", "gene_id2", "Posterior_20", "Anterior_20"]]
    lab_data = lab_data.rename({"gene_id1": "Post_name", "gene_id2": "Ant_Name"})
    
    ut_data: pandas.DataFrame = pandas.read_pickle(ut_file)
    ut_data = ut_data[["Henst", "Tenst", "Hchr", "Tchr", 'Posterior_20', 'Anterior_20']]
    ut_data = ut_data.rename({"Henst": "Post_Name", "Tenst": "Ant_Name", "Hchr": "Post_Chr", "Tchr": "Ant_Chr"})
    
    lab_data.to_excel(data_files / "Cancer_Data_Lab.xlsx")
    ut_data.to_excel(data_files / "Cancer_Data_UT.xlsx")




def NormalData(source_data: pathlib.Path, exon_out_path: pathlib.Path, intron_out_path: pathlib.Path, min_length: int = 6, *args, **kwargs):
    '''
    '''
    def ant_post(ant: str, post: str, *args, **kwargs) -> tuple:
        '''
        grabs the posterior and anterior

        Returns the Anterior and Posterior sequeces in that order
        '''
        anterior = ant[0: min_length].upper()
        posterior = post[len(post) - min_length: len(post)].upper()

        return anterior, posterior

    with open(source_data, "rb") as pfile:
        data_dict: dict = pickle.load(pfile)
    
    selected_data = random_keys(data_dict, 500)

    drg_introns = pandas.DataFrame(columns = ["Name", "Chrm", "Strand", f"Posterior_{min_length}", "Intron", f"Anterior_{min_length}"])
    drg_exons = pandas.DataFrame(columns = ["Name", "Chrm", "Strand", f"Posterior_{min_length}", "Exon", f"Anterior_{min_length}"])
    I, E = 0, 0
    gene: Gene.Gene
    for name, gene in selected_data.items():
        introns = gene.intron_seq
        exons = gene.exon_seq
        chrome = gene.chrm
        strand = gene.strand

        for i, intron in enumerate(introns):
            try:
                len_introns = len(intron)
            except Exception as e:
                len_introns = 0

            try:
                len_ant_in = len(introns[i + 1])
            except Exception as e:
                len_ant_in = 0

            try:
                len_pos_in = len(introns[i - 1])
            except Exception as e:
                len_pos_in = 0

            try:
                len_ant_ex = len(exons[i + 1])
            except Exception as e:
                len_ant_ex = 0
            
            try:
                len_pos_ex = len(exons[i])
            except Exception as e:
                len_pos_ex = 0

            ## Introns
            if ((len_pos_ex >= min_length) and (len_introns >= min_length) and (len_ant_ex >= min_length)):
                if strand in "-":
                    try:
                        # anterior = exons[i][0: min_length].upper() #[len(exons[i]) - 10: len(exons[i])][::-1]  # last 10
                        # posterior = exons[i + 1][len(exons[i + 1]) - min_length: len(exons[i + 1])].upper() #[0:10][::-1]  # first 10
                        anterior, posterior = ant_post(exons[i], exons[i + 1])
                    except Exception as e:
                        continue

                else:
                    try:
                        anterior, posterior = ant_post(exons[i + 1], exons[i])
                    except Exception as e:
                        continue

                things_to_add = [name, chrome, strand, posterior, intron, anterior]
                drg_introns.loc[len(drg_introns.index)] = things_to_add
    
            if i > 0:  ## because if we just did this we would get the 0th position exon which would give an index error
                       ## Note the len of anterior exon being used: that's beause 
            ## Exons
                if ((len_pos_in >= min_length) and (len_ant_ex >= min_length) and (len_ant_in >= min_length)):
                    try:
                        exon = exons[i + 1]
                    except Exception as e:
                        continue
                    if strand in "-":
                        try:
                            anterior, posterior = ant_post(introns[i], introns[i + 1])
                        except Exception as e:
                            continue
                    else:
                        try:
                            anterior, posterior = ant_post(introns[i + 1], introns[i])
                        except Exception as e:
                            continue

                things_to_add = [name, chrome, strand, posterior, exon, anterior]
                drg_exons.loc[len(drg_introns.index)] = things_to_add

    exon_out_path.parent.mkdir(parents = True, exist_ok = True)
    intron_out_path.parent.mkdir(parents = True, exist_ok = True)

    drg_introns.to_excel(intron_out_path, sheet_name = "sheet1", index = False)
    drg_exons.to_excel(exon_out_path, sheet_name = "sheet1", index = False)



def compliment(sequence: str):
    '''
    This is antiquated and should not be used.

    Somewhere I screwup the data with the antisense genes.

    This fixes it, somehow? I don't understand biology, nor do I want to work with biologist at all. I don't find it interesting and never should have picked up this project.

    I might not need this, I might be able to just do the list[::-1][slice], but I'm over this project so deal with it.
    '''
    sequence = sequence.upper()
    ecneuqes = ""
    for n in sequence:
        if n in "A":
            ecneuqes = f"T{ecneuqes}"
        if n in "T":
            ecneuqes = f"A{ecneuqes}"
        if n in "C":
            ecneuqes = f"G{ecneuqes}"
        if n in "G":
            ecneuqes = f"C{ecneuqes}"

    return ecneuqes



def random_keys(dictionary: dict, n: int, half = False) -> dict:
    '''
    '''
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


if __name__ in "__main__":
    main()