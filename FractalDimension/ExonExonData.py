import pandas
import pathlib
cwd = pathlib.Path.cwd()
import GeneClass as Gene
import pickle
import random
import numpy as np
from matplotlib import pyplot as plt
from TimeEmbedding import time_embedding as te


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
    master_data: pandas.DataFrame = gen_points(data = data, n = 1_000)
    # print(master_data)
    master_data.to_pickle(save_file_path)  # saves the data

    # file_plot_names = []

    # plots(save_file_path)




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



def gen_points(*args, **kwargs) -> pandas.DataFrame:
    '''
    '''
    genes = collect_genes(*args, **kwargs)  # I'm trying something
    master_data = pandas.DataFrame()

    for g, gene in enumerate(genes):
        local_xy = exon2exon(gene, *args, **kwargs)

        if isinstance(local_xy, pandas.DataFrame):
            local_xy["Gndex"] = g
            master_data = pandas.concat([master_data, local_xy])  # joins them by rows

    return master_data



def exon2exon(gene: Gene.Gene, k: int = 6, *args, **kwargs) -> pandas.DataFrame:
    '''
    Takes a gene class. Iterates over the exon list and grabs the last k-mer of the previous exon and the first k-mer of the next exon

    Returns an dataframe with "X", "Y", "Chrm"
    '''
    exons: list = gene.exon_seq
    index: int = len(exons)
    xy_points: pandas.DataFrame = pandas.DataFrame(None, index = [i for i in range(index - 1)], columns = ["X", "Y", "Chrm"])

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