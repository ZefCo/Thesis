import pandas
import numpy as np
import os
import pathlib
cwd = pathlib.Path.cwd()
import MomentCalculations as MC
import Heatmaps as heat
from typing import Tuple


def main():
    '''
    Does the Moment calculations for the fly and Mouse.
    '''
    # _, _ = gen_heat_data(cwd / "TE_Images_ForPaper" / "Dict" / "Mouse_Seq_For_Images.pkl", cwd / "Mouse_Dicts" / "Exon_6mer.pkl", cwd / "Mouse_Dicts" / "Intron_6mer.pkl")
    # _, _ = gen_heat_data(cwd / "TE_Images_ForPaper" / "Dict" / "Fly_Seq_For_Images.pkl", cwd / "Fly_Dicts" / "Exon_6mer.pkl", cwd / "Fly_Dicts" / "Intron_6mer.pkl")

    # gen_moments(cwd / "Mouse_Dicts", mer_out = cwd / "TE_Images_ForPaper" / "Moments" / "Mouse_Moments.png", pd_out = None)
    # gen_moments(cwd / "Fly_Dicts", mer_out = cwd / "TE_Images_ForPaper" / "Moments" / "Fly_Moments.png", pd_out = None)
    all_species_plots()


def all_species_plots():
    '''
    Plots all 3 species on a single plot.

    Probably should have broken the methods up earlier so that I can grab the actual plot data, then plot it.

    I'm probably just going to hard code this
    '''
    step = 100
    max_n = 2
    min_n = 0.5
    ms = [n / step for n in range(int(step*(min_n)), int(step*max_n) + 1)]


    fly_dir = cwd / "Fly_Dicts"
    mou_dir = cwd / "Mouse_Dicts"
    hum_dir = cwd / "Dicts"

    hum_me, hum_mi, uni, _ = MC.moments_v3(hum_dir, ms, 6, 6, N_value = True)
    mou_me, mou_mi, _, _ = MC.moments_v3(mou_dir, ms, 6, 6, N_value = True)
    fly_me, fly_mi, _, _ = MC.moments_v3(fly_dir, ms, 6, 6, N_value = True)

    me = {"Human": hum_me, "Mouse": mou_me, "Fly": fly_me}
    mi = {"Human": hum_mi, "Mouse": mou_mi, "Fly": fly_mi}

    MC.multiple_species_plots(ms, me, mi, uni, cwd / "TE_Images_ForPaper" / "AllSpecies", x_ticks={0.5: 0.5, 1.0: 1.0, 1.5: 1.5, 2.0: 2.0})





def gen_moments(filepath: pathlib.Path, *args, **kwargs):
    '''
    Takes the heat data and generates the moments
    I really don't need a method for this... I'm not sure why I am going to stick with it hard coded like this... but whatever
    '''
    ms = [n / 100 for n in range(int(100*(0.5)), int(100*2) + 1)]

    MC.moments_v2(filepath, ms, min_k = 6, max_k = 6, N_value = True, x_ticks = {0.5: 0.5, 1.0: 1.0, 1.5: 1.5, 2.0: 2.0}, *args, **kwargs)


def gen_heat_data(filepath: pathlib.Path, exon_output: pathlib.Path = None, intron_output: pathlib.Path = None) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
    '''
    Generates the heat data
    '''
    if isinstance(exon_output, pathlib.Path):
        exon_dir: pathlib.Path = exon_output.parent
        exon_dir.mkdir(exist_ok = True, parents = True)
    if isinstance(intron_output, pathlib.Path):
        intron_dir: pathlib.Path = intron_output.parent
        intron_dir.mkdir(exist_ok = True, parents = True)

    _, exon, intron, _, _, _, _, _, _ = heat.heat_embedding(filepath, just_import = True)

    exon = pandas.DataFrame(exon)
    exon = heat._reorder_frame(exon)
    exon = MC._unrenormalize(exon, 12, log2 = False)

    intron = pandas.DataFrame(intron)
    intron = heat._reorder_frame(intron)
    intron = MC._unrenormalize(intron, 12, log2 = False)

    if isinstance(exon_output, pathlib.Path):
        exon.to_pickle(exon_output)
    if isinstance(intron_output, pathlib.Path):
        intron.to_pickle(intron_output)

    return exon, intron


if __name__ in "__main__":
    main()