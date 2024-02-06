import ExonExonData as eed
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import pickle
from Heatmaps import heatmapv2 as heatmap
from Heatmaps import _undigitize_seq as undigit


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
    # folder_path = pathlib.Path("F:\Gene_Data_Sets")  # windows path
    folder_path = cwd.parent.parent / "Thesis_Data" / "Gene_Data_Sets"  #Linux path

    save_file_path = folder_path / f"E2E_DS{data_set}_k{k}.pkl"
    data: pathlib.Path = folder_path / f"Data_Set_{data_set}_cleaned_dict.pkl"
    genes = eed.collect_genes(data, n = 10_000)  # I'm trying something

    
    heatmaps(genes, cwd / "Sanity_HE.png", k = k)
    plots(genes, cwd / "Dumping.pkl", cwd / "Sanity_TE_v2.png")
    # heatmap(data, colors=["white", "black"], bounds = [0, 0.1, 0.2], fileoutput=cwd / "Sanity.png")



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



def plots(genes, save_file_path: pathlib.Path, plot_path: pathlib.Path, *args, **kwargs):
    '''
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