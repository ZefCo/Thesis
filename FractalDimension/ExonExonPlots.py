import ExonExonData as eed
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import pickle
from Heatmaps import heatmapv2 as heatmap
from Heatmaps import _undigitize_seq as undigit
from TimeEmbedding import time_embedding_1D as te_1d


def main():
    '''
    These are not going to be one size fits all scripts. I want to start doing a file with all my functions and then other files that call those functions. Those other files will be easier to manipulate.
    '''
    # # This is becase I screwed up and forgot to rename the columns with the sequence
    # data = pandas.read_csv(cwd / "Raw_Exon2Exon.csv", dtype={"Unnamed: 0": str})
    # data = data.set_index("Unnamed: 0")
    # digital_seqs = tuple(data.index)
    # # # print(digital_seqs)
    # # # print(type(digital_seqs[0]))
    # new_seq = dict()
    # for dig_seq in digital_seqs:
    #     str_seq = undigit(dig_seq)
    #     new_seq[dig_seq] = str_seq
    # data = data.rename(columns = new_seq, index = new_seq)
    # # print(data)
    # # data.to_csv(cwd / "Raw_Exon2Exon.csv")
    data_set = 2
    k = 6
    folder_path = pathlib.Path("F:\Gene_Data_Sets")  # windows path
    # folder_path = cwd.parent.parent / "Thesis_Data" / "Gene_Data_Sets"  #Linux path
    # save_file_path = folder_path / f"E2E_DS{data_set}_k{k}_Seq.pkl"


    # heatmaps(genes, cwd / "Sanity_HE.png", k = k)
    # plots(genes, cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "E2ENormal.pkl", cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "E2ENormal.png")
    # heatmap(data, colors=["white", "black"], bounds = [0, 0.1, 0.2], fileoutput=cwd / "Sanity.png")


    # # Creating Exon2Exon Normal data of just sequences
    # save_file_path = folder_path / f"Normal_Exon2Exon.xlsx"
    # data: pathlib.Path = folder_path / f"Data_Set_{data_set}_cleaned_dict.pkl"
    # genes = eed.collect_genes(data, n = 10_000)  # I'm trying something
    # master_data = pandas.DataFrame()
    # for i, g in enumerate(genes):
    #     something = eed.exon2exonSeq(g, k = 20)
    #     master_data = pandas.concat([master_data, something], ignore_index = True)

    #     master_data = master_data.dropna()

    # master_data.to_excel(save_file_path, index = False)


    drg_data = pandas.read_excel(cwd.parent / "Data_Files" / "Normal_Exon2Exon_v2.xlsx", sheet_name = "Sheet1")
    anterior = drg_data["Anterior_10"]
    posterior = drg_data["Posterior_10"]

    anterior = anterior.apply(lambda x : x[0:6])
    posterior = posterior.apply(lambda x : x[len(x) - 6: len(x)])

    master_data = pandas.DataFrame()
    master_data["Anterior"] = anterior
    master_data["Posterior"] = posterior

    plots_v2(master_data, cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "Dump.pkl", cwd / "TE_Images_ForPaper" / "Exon2Exon" / "Cancer" / "E2ETest.png")


    # print(type(anterior))


def plots_v2(genes: pandas.DataFrame, save_file_path: pathlib.Path, plot_path: pathlib.Path, *args, **kwargs):
    '''
    assumes the data is already formated as a dataframe with a posterior and anterior sequence
    '''

    save_file_path.parent.mkdir(parents = True, exist_ok = True)

    marker = ","
    title = "Next Exon (X) to Previous Exon (Y)"

    inches = 10
    
    print("Generating points")
    # master_data: pandas.DataFrame = eed.gen_te_points(genes, *args, **kwargs)
    genes["Fusion_Seq"] = genes["Posterior"] + genes["Anterior"]
    genes["X"] = genes["Anterior"].apply(lambda x: te_1d(x)[0])
    genes["Y"] = genes["Posterior"].apply(lambda x: te_1d(x, w_backwards = True)[0])

    print("Saving to file")
    genes.to_pickle(save_file_path)  # saves the data
    genes.to_csv(cwd / "Sanity_TE.csv")

    print("Reopening file")
    with open(save_file_path, "rb") as file:
        master_data: pandas.DataFrame = pickle.load(file)

    # try:
    #     genes["Chrm"] = genes["Chrm"].astype("category")
    # except Exception as e:
    #     pass
    # chromosome = tuple(genes["Chrm"].unique())
    # gene_index = tuple(genes["Gndex"].unique())
    print("Generating Plots")
    eed.gen_plot(genes, 
                 plot_path, 
                 inches = inches, 
                 title = title, 
                 marker = marker)
    



def heatmaps(genes, plot_path: pathlib.Path, *args, **kwargs):
    '''
    Generates a 'raw' heatmap, just the data, in a csv/xlsx file for Dr. G.
    '''

    print("Generating Heat Points")
    data: pandas.DataFrame = eed.gen_he_points(genes, str_names = True, *args, **kwargs)

    # data = pandas.read_csv(cwd / "RawExon2Exon.csv", index_col="Unnamed: 0")
    # print(data)
    
    print("Generating Heat Maps")
    heatmap(data, colors=["white", "black"], bounds = [0, 0.01, 0.1], fileoutput = plot_path)
    print("Saving to file")
    data.to_csv(cwd / "Sanity_HE.csv")



    # with open(data, "rb") as file:
    #     data: pandas.DataFrame = pickle.load(file)

    # data = eed.collect_genes(data)

    # print(data)



def plots(genes: list, save_file_path: pathlib.Path, plot_path: pathlib.Path, *args, **kwargs):
    '''
    genes is a list of Gene.Genes

    redo this with 4-mers and 5-mers. Adjust box size as appropriate.
    '''
    marker = ","
    title = "Next Exon (X) to Previous Exon (Y)"

    inches = 10
    
    print("Generating points")
    master_data: pandas.DataFrame = eed.gen_te_points(genes, *args, **kwargs)
    print("Saving to file")
    master_data.to_pickle(save_file_path)  # saves the data
    master_data.to_csv(cwd / "Sanity_TE.csv")

    print("Reopening file")
    with open(save_file_path, "rb") as file:
        master_data: pandas.DataFrame = pickle.load(file)

    master_data["Chrm"] = master_data["Chrm"].astype("category")
    chromosome = tuple(master_data["Chrm"].unique())
    gene_index = tuple(master_data["Gndex"].unique())
    print("Generating Plots")
    eed.gen_plot(master_data, 
                 plot_path, 
                 inches = inches, 
                 title = title, 
                 marker = marker)

    # for c, chrm in enumerate(chromosome):
    #     subdata = data[data["Chrm"] == chrm]
    #     eed.gen_plot(subdata, cwd / "TE_Images_ForPaper" / "Exon2Exon"/ f"K{k}" / f"DS{data_set}_{chrm}_k{k}.png", inches = inches)



if __name__ in "__main__":
    main()