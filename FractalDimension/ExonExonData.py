import pandas
import pathlib
cwd = pathlib.Path.cwd()
import GeneClass as Gene
import pickle
import random
import numpy as np
from matplotlib import pyplot as plt
from TimeEmbedding import time_embedding as te
from Heatmaps import _heat_data as he
from Heatmaps import _heat_data_v2 as he2
from Heatmaps import _dict_deep_merge as merge
from Heatmaps import _reorder_frame as reorder
from Heatmaps import _init_dict as initD
from Heatmaps import _undigitize_seq as undigit


def main():
    '''
    Look at the head of one exon, then the tail of the previous exon: make those x and y respectivly and plot.

    Randomly select many genes.

    Plot these all together and by chromosome.

    Keep track of how many genes you take and for which chromosome they come from. 
    '''
    folder_path = pathlib.Path("F:\Gene_Data_Sets")
    data: pathlib.Path = folder_path / "Data_Set_1_cleaned_dict.pkl"
    save_file_path = folder_path / "E2E_DS1_k6.pkl"
    master_data: pandas.DataFrame = gen_te_points(data = data, n = 1_000)
    # print(master_data)
    master_data.to_pickle(save_file_path)  # saves the data

    # file_plot_names = []

    plots(save_file_path, file_path = cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "E2ENormal.png")




def plots(data: pathlib.Path or pandas.DataFrame, *args, **kwargs):
    '''
    Generates the plots from an input file
    
    Actually it feeds it to the gen_plots function

    This is not a thing you can just feed data into: you're going to have to make adjustments to it. Sorry. Probably should be it's own script really.
    '''

    if isinstance(data, pathlib.Path):
        with open(data, "rb") as file:
            data: pandas.DataFrame = pickle.load(file)


    data["Chrm"] = data["Chrm"].astype("category")
    chromosome = tuple(data["Chrm"].unique())
    gen_plot(data, None)

    for c, chrm in enumerate(chromosome):
        subdata = data[data["Chrm"] == chrm]
        gen_plot(subdata, None)




def gen_plot(data: pandas.DataFrame, file_name: pathlib.Path = None, inches: int = 5, linestyle = "", marker = ".", title: str = None, *args, **kwargs):
    '''
    This does the actual plots
    '''
    fig, ax = plt.subplots()
    fig.set_size_inches(inches, inches)


    x: tuple = tuple(data["X"].to_list())
    y: tuple = tuple(data["Y"].to_list())

    plt.plot(x, y, linestyle = linestyle, marker = marker, color = "k", markersize = 0.5)

    plt.xlabel("X")
    plt.ylabel("Y")

    if isinstance(title, str):
        plt.title(title)

    if isinstance(file_name, pathlib.Path):
        file_name.parent.mkdir(exist_ok = True, parents = True)
        plt.savefig(file_name)
        plt.close()
    else:
        plt.show()


def gen_he_points(genes, str_names: bool = False, *args, **kwargs) -> pandas.DataFrame:
    '''
    Generates a heatmap

    Takes 
    '''
    _, _, master_dict = initD(*args, **kwargs)  # really just need the permutations
    heat_data = pandas.DataFrame(0, columns = tuple(master_dict.keys()), index = tuple(master_dict.keys()))

    total_genes = 0
    for g, gene in enumerate(genes):
        heat_data = heat2heat_v2(gene, heat_data, *args, **kwargs)
        total_genes += 1
        # if isinstance(heat_data, pandas.DataFrame):
            #     master_dict = merge(master_dict, gene_he)

    # master_data: pandas.DataFrame = pandas.DataFrame(master_dict)
    heat_data = reorder(heat_data, *args, **kwargs)

    # This part converts the sequences, which are in a numberical form - I call it digital - into a string form.

    if str_names:
        digital_seqs = tuple(heat_data.index)
        new_seq = dict()
        for dig_seq in digital_seqs:
            str_seq = undigit(dig_seq)
            new_seq[dig_seq] = str_seq
        
        heat_data = heat_data.rename(columns=new_seq, index=new_seq)

    print(f"Total genes used: {total_genes}")

    return heat_data
    


def gen_te_points(genes: list, *args, **kwargs) -> pandas.DataFrame:
    '''
    Generates a series of points.
    '''
    master_data = pandas.DataFrame()

    total_genes = 0
    for g, gene in enumerate(genes):
        local_xy = exon2exon(gene, *args, **kwargs)

        if isinstance(local_xy, pandas.DataFrame):
            local_xy["Gndex"] = g
            master_data = pandas.concat([master_data, local_xy])  # joins them by rows
            total_genes += 1


    print(f"Total genes used: {total_genes}")

    return master_data



def heat2heat_v2(gene: Gene.Gene, dataframe: pandas.DataFrame, k: int = 6, *args, **kwargs):
    '''
    this works soley with dataframse.
    '''
    exons: list = gene.exon_seq
    index: int = len(exons)
    
    for i in range(1, index):
        # we're just going to look at exons and introns that are at minimum len >= 6
        try:
            len_x = len(exons[i])
            len_y = len(exons[i - 1])
        except Exception as e:
            # return dataframe
            continue

        if (len_x >= k) and (len_y >= k):

            exon_x = exons[i][:k]  # each entry in the list is a string, this grabs the first k nucleotides
            exon_y = exons[i - 1][-k:]  # this grabs the last k nucleotides

            dataframe = he2(f"{exon_y}{exon_x}", dataframe, k_p = k, k_m = k)  # I should have stayed in Burbank.

        else:
            # return None
            continue
        
    return dataframe




def heat2heat(gene: Gene.Gene, local_xy: dict, k: int = 6, *args, **kwargs) -> dict:
    '''
    Takes in a gene and looks at the exons, then gets a heatmap dataset for beginning and end of the exons.

    Returns this dictionary. If any of the exons are not long enough, it discards the whole gene.

    Local_xy is because I need to init the dictionary with all possible permuations. Instead of doing that over and over again, you need
    to init that dictionary twice before calling this function: one is to be used as the blank local_xy, while the other is the master_xy you add too.
    '''
    exons: list = gene.exon_seq
    index: int = len(exons)
    # xy_points: pandas.DataFrame = pandas.DataFrame(None, index = [i for i in range(index - 1)], columns = ["X", "Y", "Chrm"])

    for i in range(1, index):
        # we're just going to look at exons and introns that are at minimum len >= 6
        try:
            len_x = len(exons[i])
            len_y = len(exons[i - 1])
        except Exception as e:
            return None

        if (len_x >= k) and (len_y >= k):

            exon_x = exons[i][:k]  # each entry in the list is a string, this grabs the first k nucleotides
            exon_y = exons[i - 1][-k:]  # this grabs the last k nucleotides

            small_xy: dict = he(f"{exon_y}{exon_x}", k_p = k, k_m = k)  # I should have stayed in Burbank.

            local_xy = merge(local_xy, small_xy)        
        else:
            return None
        
    return local_xy



def exon2exonSeq(gene: Gene.Gene, k: int = 10, *args, **kwargs) -> pandas.DataFrame:
    '''
    Just returns the sequence
    '''
    exons: list = gene.exon_seq
    index: int = len(exons)
    seq_points: pandas.DataFrame = pandas.DataFrame(None, index = [i for i in range(index - 1)], columns = ["Name", "Chrm", f"Posterior_{k}", f"Anterior_{k}"])

    for i in range(1, index):
            # we're just going to look at exons and introns that are at minimum len >= k
        try:
            len_x = len(exons[i])
            len_y = len(exons[i - 1])
        except Exception as e:
            return None

        if (len_x >= k) and (len_y >= k):

            exon_x = exons[i][:k]  # each entry in the list is a string, this grabs the first k nucleotides
            exon_y = exons[i - 1][-k:]  # this grabs the last k nucleotides

            seq_points.loc[i, "Name"] = gene.name
            seq_points.loc[i, "Chrm"] = gene.chrm
            seq_points.loc[i, f"Posterior_{k}"] = exon_y.upper()
            seq_points.loc[i, f"Anterior_{k}"] = exon_x.upper()

    return seq_points





def exon2exon(gene: Gene.Gene, k: int = 6, *args, **kwargs) -> pandas.DataFrame:
    '''
    Takes a gene class. Iterates over the exon list and grabs the last k-mer of the previous exon and the first k-mer of the next exon

    Returns an dataframe with "X", "Y", "Chrm"

    Uses the time embedding method. Can't really think of a way to put the heat embedding in this one since they return fundamentally different things.
    '''
    exons: list = gene.exon_seq
    index: int = len(exons)
    xy_points: pandas.DataFrame = pandas.DataFrame(None, index = [i for i in range(index - 1)], columns = ["X", "Y", "Chrm", "Seq"])

    for i in range(1, index):
        # we're just going to look at exons and introns that are at minimum len >= 6
        try:
            len_x = len(exons[i])
            len_y = len(exons[i - 1])
        except Exception as e:
            return None

        if (len_x >= k) and (len_y >= k):

            exon_x = exons[i][:k]  # each entry in the list is a string, this grabs the first k nucleotides
            exon_y = exons[i - 1][-k:]  # this grabs the last k nucleotides

            # try:
            #     exon_x = exons[i][:k]  # each entry in the list is a string, this grabs the first k nucleotides
            # except Exception as e:
            #     continue
            # try:
            #     exon_y = exons[i - 1][-k:]  # this grabs the last k nucleotides
            # except Exception as e:
            #     continue  # so if there's any error it just skips to the next one

            # print(f"{exon_y}{exon_x}")  # remember: half way through I had to switch the x and y coords... fuck I hate this project
            xy = te(f"{exon_y}{exon_x}", k_p = k, k_m = k)  # I should have stayed in Burbank.

            if isinstance(xy, np.ndarray):
                xy_points.loc[i - 1, "X"], xy_points.loc[i - 1, "Y"], xy_points.loc[i - 1, "Chrm"] = xy[0][1], xy[0][0], gene.chrm
                xy_points.loc[i - 1, "Seq"] = f"{exon_y}{exon_x}"
        
        else:
            return None

    # print(xy_points)

    return xy_points




def collect_genes(data: pathlib.Path or dict, n: int = 1_000, *args, **kwargs) -> tuple:
    '''
    n is the nubmer of genes to import
    '''
    if isinstance(data, pathlib.Path):
        with open(data, "rb") as file:
            data: dict = pickle.load(file)

    gene = set()
    value: Gene.Gene
    for key, value in data.items():
        gene.add(value.chrm)

    # print(gene)
    # exit()

    data = tuple(random.choices(tuple(data.values()), k = n))  # I think this will work.

    return data


if __name__ in "__main__":
    main()