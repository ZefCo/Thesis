import ExonExonData as eed
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import pickle
from Heatmaps import heatmapv2 as heatmap


def main():
    '''
    These are not going to be one size fits all scripts. I want to start doing a file with all my functions and then other files that call those functions. Those other files will be easier to manipulate.
    '''
    heatmaps()


def heatmaps():
    '''
    Generates a 'raw' heatmap, just the data, in a csv/xlsx file for Dr. G.
    '''
    data_set = 2
    k = 6
    # folder_path = pathlib.Path("F:\Gene_Data_Sets")  # windows path
    folder_path = cwd.parent.parent / "Thesis_Data" / "Gene_Data_Sets"  #Linux path

    save_file_path = folder_path / f"E2E_DS{data_set}_k{k}.pkl"
    data: pathlib.Path = folder_path / f"Data_Set_{data_set}_cleaned_dict.pkl"
    data: pandas.DataFrame = eed.gen_he_points(data = data, k = k, n = 10_000)

    # data = pandas.read_csv(cwd / "RawExon2Exon.csv", index_col="Unnamed: 0")
    # print(data)

    heatmap(data, colors=["green", "darkred"], bounds = [0, 0.1, 0.2], fileoutput=cwd / "Sanity.png")

    data.to_csv(cwd / "Raw_Exon2Exon.csv")



    # with open(data, "rb") as file:
    #     data: pandas.DataFrame = pickle.load(file)

    # data = eed.collect_genes(data)

    # print(data)



def plots():
    '''
    redo this with 4-mers and 5-mers. Adjust box size as appropriate.
    '''
    data_set = 2
    k = 6
    marker = "s"
    title = "Next Exon (X) to Previous Exon (Y)"

    folder_path = pathlib.Path("F:\Gene_Data_Sets")
    data: pathlib.Path = folder_path / f"Data_Set_{data_set}_cleaned_dict.pkl"
    save_file_path = folder_path / f"E2E_DS{data_set}_k{k}.pkl"
    inches = 10
    n = 10_000
    
    print("Generating points")
    master_data: pandas.DataFrame = eed.gen_points(data = data, n = n, k = k)
    print("Saving to file")
    master_data.to_pickle(save_file_path)  # saves the data

    print("Reopening file")
    with open(save_file_path, "rb") as file:
        data: pandas.DataFrame = pickle.load(file)

    data["Chrm"] = data["Chrm"].astype("category")
    chromosome = tuple(data["Chrm"].unique())
    gene_index = tuple(data["Gndex"].unique())
    print("Generating Plots")
    eed.gen_plot(data, 
                 cwd / "TE_Images_ForPaper" / "Exon2Exon"/ f"K{k}" / f"E2E_DS{data_set}_k{k}.png", 
                 inches = inches, 
                 title = title, 
                 marker = marker)

    for c, chrm in enumerate(chromosome):
        subdata = data[data["Chrm"] == chrm]
        eed.gen_plot(subdata, cwd / "TE_Images_ForPaper" / "Exon2Exon"/ f"K{k}" / f"DS{data_set}_{chrm}_k{k}.png", inches = inches)



if __name__ in "__main__":
    main()