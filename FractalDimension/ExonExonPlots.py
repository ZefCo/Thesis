import ExonExonData as eed
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import pickle

def main():
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