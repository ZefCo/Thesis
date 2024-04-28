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
    element_drive = pathlib.Path("F:\Gene_Data_Sets")
    gene_drive = pathlib.Path("F:\GeneData")
    file = element_drive / "Data_Set_1_cleaned_dict.pkl"
    exon_out = gene_drive / "Exon.xlsx"
    intron_out = gene_drive / "Intron.xlsx"

    with open(file, "rb") as pfile:
        data_dict: dict = pickle.load(pfile)
    
    selected_data = random_keys(data_dict, 5000)

    drg_introns = pandas.DataFrame(columns = ["Name", "Chrm", "Posterior_10", "Intron", "Anterior_10"])
    drg_exons = pandas.DataFrame(columns = ["Name", "Chrm", "Posterior_10", "Exon", "Anterior_10"])
    I, E = 0, 0
    gene: Gene.Gene
    for name, gene in selected_data.items():
        introns = gene.intron_seq
        exons = gene.exon_seq
        chrome = gene.chrm

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
                len_ant_ex = len(exons[i])
            except Exception as e:
                len_ant_ex = 0
            
            try:
                len_pos_ex = len(exons[i + 1])
            except Exception as e:
                len_pos_ex = 0

            if ((len_introns > 0) and (len_pos_ex >= 10) and (len_ant_ex >= 10)):
                posterior = exons[i][len(exons[i]) - 10: len(exons[i])]  # last 10
                anterior = exons[i + 1][0:10]  # first 10
                things_to_add = [name, chrome, posterior, intron, anterior]
                drg_introns.loc[len(drg_introns.index)] = things_to_add
    
                # drg_introns.loc[I, "Name"] = name
                # drg_introns.loc[I, "Chrm"] = chrome
                # drg_introns.loc[I, "Posterior_10"] = posterior
                # drg_introns.loc[I, "Intron"] = intron
                # drg_introns.loc[I, "Anterior_10"] = anterior
                # I += 1

            if i < len(introns):

                if ((len_pos_ex > 0) and (len_introns >= 10) and (len_ant_in >= 10)):
                    posterior = introns[i][len(introns[i]) - 10: len(introns[i])]  # last 10
                    anterior = introns[i + 1][0:10]  # first 10
                    things_to_add = [name, chrome, posterior, exons[i + 1], anterior]
                    drg_exons.loc[len(drg_exons.index)] = things_to_add

                    # drg_exons.loc[E, "Name"] = name
                    # drg_exons.loc[E, "Chrm"] = chrome
                    # drg_exons.loc[E, "Posterior_10"] = posterior
                    # drg_exons.loc[E, "Exon"] = exons[i + 1]
                    # drg_exons.loc[E, "Anterior_10"] = anterior
                    # E += 1

    drg_introns.to_excel(intron_out, sheet_name = "sheet1", index = False)
    drg_exons.to_excel(exon_out, sheet_name = "sheet1", index = False)


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