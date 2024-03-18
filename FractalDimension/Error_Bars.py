import pandas
import pathlib
cwd = pathlib.Path.cwd()
import pickle
import Heatmaps as heat
import MomentCalculations as MC
import sys
import os
import re


def main():
    '''
    '''
    # collect_data(cwd / "TE_Images_ForPaper" / "Error_Data", cwd / "TE_Images_ForPaper" / "Dict" / "Seq_For_Images_n100000_minLength12.pkl", n = 50_000)
    
    sub_moments = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0]
    sub_data_dir = cwd / "TE_Images_ForPaper" / "Error_Data"
    sub_stats = error_data(sub_data_dir, sub_moments, pickle_save = cwd / "TE_Images_ForPaper" / "Dict" / "SubStatsFrame.pkl")
    plot_with_error([sub_stats])



def plot_with_error(errors: list):
    '''
    Adds the error bars to the plot.
    I'm going to be lazy right now (1/9/24) and just hard code this.
    '''
    ms = [n / 100 for n in range(int(100*(0.5)), int(100*2) + 1)]

    MC.moments_v2(cwd / "Dicts", ms, min_k = 6, max_k = 6, N_value = True, x_ticks = {0.5: 0.5, 1.0: 1.0, 1.5: 1.5, 2.0: 2.0}, error_bars = errors, y_ticks={0: 0, 5: 5})


def error_data(sub_data_dir: pathlib.Path, sub_moments: list, k: int = 6, pickle_save: pathlib.Path = None) -> pandas.DataFrame:
    '''
    Figures out the standard deviation of the ranges of moments. Calculates the moment for each sampled data and then puts it into a larger
    dataframe to figure out the standard deviation (by row).

    The dataframe that stores all the moments and stds uses the index as the moment calculation. Doesn't find the moments for every single point,
    just some of them.

    sub_moments is a list of moments to calculate. I don't want to deal with some weird for loop and ranges and whatnot.

    Pickle can be used to save the data. Returns a dataframe
    '''
    sub_name: pathlib.Path

    sub_samples = list(sub_data_dir.iterdir())
    sub_names = [sub_name.name for sub_name in sub_samples]

    sub_ranges = pandas.DataFrame(0, columns = sub_moments, index = sub_names)
    sub_ranges["Class"] = "None"

    for sub_name in sub_names:
        classification = re.split("_", re.split(r"\.", sub_name)[0])[0].lower()  # spliting the name on the . and the first underscore to get Exon or Intron
        sub_ranges.loc[sub_name, "Class"] = classification

    for i, sub in enumerate(sub_samples):
        for m in sub_moments:
            N = 4**(2*k)
            N = N**((1/m) - 1)

            local_m = MC.moment(sub, m, N)

            sub_ranges.loc[sub_names[i], m] = local_m

    sub_exons = sub_ranges[sub_ranges["Class"] == "exon"]
    sub_exons_std: pandas.Series = sub_exons.std(axis = 0, numeric_only = True)
    sub_exons_mean: pandas.Series = sub_exons.mean(axis = 0, numeric_only = True)
    sub_introns = sub_ranges[sub_ranges["Class"] == "intron"]
    sub_introns_std: pandas.Series = sub_introns.std(axis = 0, numeric_only = True)
    sub_introns_mean: pandas.Series = sub_introns.mean(axis = 0, numeric_only = True)

    # print(sub_exons_std)
    # print(sub_introns_std)
    # print(type(sub_introns_std))

    sub_stats = {"Exon_Mean": sub_exons_mean, "Exon_STD": sub_exons_std, "Intron_Mean": sub_introns_mean, "Intron_STD": sub_introns_std}
    sub_stats = pandas.concat(sub_stats, axis = 1)

    if isinstance(pickle_save, pathlib.Path):
        sub_stats.to_pickle(pickle_save)

    return sub_stats


def collect_data(file_dir: pathlib.Path, source_data: pathlib.Path, n: int = 10_000, k: int = 6):
    '''
    Takes a file path which will be used to output the file, and a source path where the data will be taken from.

    A temp file is saved so the script can be run parallel to others, saving time instead of doing it one at a time.
    '''
    # file_dir.mkdir(parents = True, exist_ok = True)
    files = int(len(list(file_dir.iterdir())) / 2)

    MC.write_pickel(pandas.DataFrame({"Nothing": [0]}), file_dir / f"Exon_EB_{files}.pkl")
    MC.write_pickel(pandas.DataFrame({"Nothing": [0]}), file_dir / f"Intron_EB_{files}.pkl")

    _, exon, intron, _, _, _, _, _, _ = heat.heat_embedding(source_data,
                                                            n = n,
                                                            k_m = k, k_p = k,
                                                            just_import = False)
    
    exon = pandas.DataFrame(exon)
    exon = heat._reorder_frame(exon)
    exon = MC._unrenormalize(exon, 2*k, log2=False)
    MC.write_pickel(exon, file_dir / f"Exon_EB_{files}.pkl")

    intron = pandas.DataFrame(intron)
    intron = heat._reorder_frame(intron)
    intron = MC._unrenormalize(intron, 2*k, log2=False)
    MC.write_pickel(intron, file_dir / f"Intron_EB_{files}.pkl")


if __name__ in "__main__":
    main()