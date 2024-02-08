import pandas
import numpy as np
import os
import pathlib
cwd = pathlib.Path.cwd()
import MomentCalculations as MC
import Heatmaps as heat
import TimeEmbedding as TE
import pickle
import GeneClass as Gene
from typing import Tuple


def main():
    '''
    This does more then just mouse and fly data. Now it looks at all different speiceis. I wont change the name, but it is supposed to be for the phylogenetic trees.
    
    Does the Moment calculations for the fly and Mouse.

    '''
    # primate_filepath = cwd.parent / "Data_Files" / "Primates" / "Genetics" # species / f"{species}_{genome}_Frame.pkl"
    # species = ["Callithrix_jacchus", "Gorilla_gorilla_gorilla" ,"Macaca_fascicularis" ,"Macaca_mulatta" ,"Nomascus_leucogenys" ,"Pan_paniscus" ,"Pan_troglodytes" ,"Papio_anubis" ,"Pongo_pygmaeus_abelii"]
    # genome = ["calJac4", "gorGor6", "macFas5", "rheMac8", "nomLeu3", "panPan3", "panTro6", "papAnu4", "ponAbe3"]

    # primate_sg = {s:genome[i] for i, s in enumerate(species)}
    # for prim, gen in primate_sg.items():
    #     local_primate_filepath = primate_filepath / prim
    #     local_primate_he: pathlib.Path = local_primate_filepath / "HEData"

    # # with open(local_primate_he / "Exon_6mer.pkl", "rb") as file:
    # #     data = pickle.load(file)

    # # print(data)
    # # print("\n\n")
    # print(prim, gen)
    # with open(local_primate_he / "Intron_6mer.pkl", "rb") as file:
    #     data = pickle.load(file)

    # print(type(data))

    # print("\n\n\nHumans\n")
    # with open(cwd / "HS_Dicts_EF" / "Exon_6mer.pkl", "rb") as file:
    #     data = pickle.load(file)

    # print(type(data))


    # species_TE(cwd.parent / "Data_Files" / "Primates" /"Genetics" / species / f"{species}_{genome}_Frame.pkl",
    #            exon_outfile = cwd.parent / "Data_Files" / "Primates" /"Genetics" / species / "Exons.pkl",
    #            intron_outfile = cwd.parent / "Data_Files" / "Primates" /"Genetics" / species / "Introns.pkl")
    

    # TE.matplotfigure(frame = cwd.parent / "Data_Files" / "Primates" /"Genetics" / species / "Exons.pkl",
    #                  dir = cwd.parent / "Data_Files" / "Primates" /"Genetics" / species,
    #                  file_name = f"{species}_Exons.png")

    # TE.matplotfigure(frame = cwd.parent / "Data_Files" / "Primates" /"Genetics" / species / "Introns.pkl",
    #                  dir = cwd.parent / "Data_Files" / "Primates" /"Genetics" / species,
    #                  file_name = f"{species}_Introns.png")

    # # _, _ = gen_heat_data(cwd / "TE_Images_ForPaper" / "Dict" / "Mouse_Seq_For_Images.pkl", 
    # #                      cwd / "Mouse_Dicts_EF" / "Exon_6mer.pkl", cwd / "Mouse_Dicts_EF" / "Intron_6mer.pkl", just_import = False, n = 20_000)
    # # _, _ = gen_heat_data(cwd / "TE_Images_ForPaper" / "Dict" / "Fly_Seq_For_Images.pkl", 
    # #                      cwd / "Fly_Dicts_EF" / "Exon_6mer.pkl", cwd / "Fly_Dicts_EF" / "Intron_6mer.pkl", just_import = False, n = 20_000)
    # # _, _ = gen_heat_data(cwd / "TE_Images_ForPaper" / "Dict" / "Seq_For_Images_n100000_minLength12.pkl", 
    # #                      cwd / "HS_Dicts_EF" / "Exon_6mer.pkl", cwd / "HS_Dicts_EF" / "Intron_6mer.pkl", just_import = False, n = 20_000)


    # # gen_moments(cwd / "Mouse_Dicts", mer_out = cwd / "TE_Images_ForPaper" / "Moments" / "Mouse_Moments.png", pd_out = None)
    # # gen_moments(cwd / "Fly_Dicts", mer_out = cwd / "TE_Images_ForPaper" / "Moments" / "Fly_Moments.png", pd_out = None)
    # all_species_plots()

    primates(just_import = False, n = 10_000, reload = True, N_value = True, great_apes = True, filename = "GreatApesMoments.png")
    # heat.heatmapv2(cwd.parent / "Data_Files" / "Primates" / "Genetics" / "Homo_sapiens" / "HEData" / "Exon_6mer.pkl", colors = ["white", "black"], bounds = [0, 1e-11, 1e-10], fileoutput=cwd / "TE_Images" / "Sanityv1.png")
    # heat.heatmapv2(cwd / "HS_Dicts_EF" / "Exon_6mer.pkl", colors = ["white", "black"], bounds = [0, 1e-11, 1e-10], fileoutput=cwd / "TE_Images" / "Sanityv2.png")
    # with open(cwd / "TE_Images_ForPaper" / "Dict" / "Seq_For_Images_n100000_minLength12.pkl", "rb") as file:
    #     data = pickle.load(file)
    #     print(data)


def primates(reload: bool = False, great_apes: bool = False, *args, **kwargs):
    '''
    This looks at the primates and does a lot.
    '''
    check_ymax = lambda x, y: x if x > y else y
    ymax = 1

    primate_filepath = cwd.parent / "Data_Files" / "Primates" / "Genetics" # species / f"{species}_{genome}_Frame.pkl"
    if great_apes:
        species = ["Gorilla_gorilla_gorilla","Pan_troglodytes"]
        genome = ["gorGor6", "panTro6"]
        ga_markers = [".", "o", "+", "x", "D", "s"]
    else:
        species = ["Callithrix_jacchus", "Gorilla_gorilla_gorilla" ,"Macaca_fascicularis" ,"Macaca_mulatta" ,"Nomascus_leucogenys" ,"Pan_paniscus" ,"Pan_troglodytes" ,"Papio_anubis" ,"Pongo_pygmaeus_abelii"]
        genome = ["calJac4", "gorGor6", "macFas5", "rheMac8", "nomLeu3", "panPan3", "panTro6", "papAnu4", "ponAbe3"]

    primate_sg = {s:genome[i] for i, s in enumerate(species)}

    me = dict()
    mi = dict()

    step = 100
    max_n = 2
    min_n = 0.5
    ms = [n / step for n in range(int(step*(min_n)), int(step*max_n) + 1)]

    if reload:
            pass
    else:
        print("Data on Homo Sapiens")
        exon, intron = gen_heat_data(cwd / "TE_Images_ForPaper" / "Dict" / "Seq_For_Images_n100000_minLength12.pkl", *args, **kwargs)
        exon.to_pickle(primate_filepath / "Homo_sapiens" / "HEData" / "Exon_6mer.pkl")
        intron.to_pickle(primate_filepath / "Homo_sapiens" / "HEData" / "Intron_6mer.pkl")

    local_primate_he = primate_filepath / "Homo_sapiens" / "HEData"
    local_me, local_mi, uni, _ = MC.moments_v3(local_primate_he, ms, 6, 6, *args, **kwargs)

    if great_apes:
        me["H.S."] = {"data": local_me, "marker": ga_markers[0]}
        mi["H.S."] = {"data": local_mi, "marker": ga_markers[1]}
    else:
        me["H.S."] = {"data": local_me, "marker": None}
        mi["H.S."] = {"data": local_mi, "marker": None}

    ymax = check_ymax(ymax, max(max(tuple(local_me.values()))))
    ymax = check_ymax(ymax, max(max(tuple(local_mi.values()))))

    print("Finished Homo Sapiens")

    for prim, gen in primate_sg.items():
        local_primate_filepath = primate_filepath / prim
        local_primate_he: pathlib.Path = local_primate_filepath / "HEData"

        if reload:
            pass
        else:
            local_primate_he.mkdir(parents = True, exist_ok = True)

            print(f"Generating Heat Data for {prim}")
            exon, intron = gen_heat_data(local_primate_filepath / f"{prim}_{gen}_Frame.pkl", *args, **kwargs)
            print("Saving Heat Data")
            exon.to_pickle(local_primate_he / "Exon_6mer.pkl")
            intron.to_pickle(local_primate_he / "Intron_6mer.pkl")

        local_me, local_mi, _, _ = MC.moments_v3(local_primate_he, ms, 6, 6, *args, **kwargs)
    
        ymax = check_ymax(ymax, max(max(tuple(local_me.values()))))
        ymax = check_ymax(ymax, max(max(tuple(local_mi.values()))))

        if great_apes:
            if prim in "Pan_troglodytes":
                me[prim] = {"data": local_me, "marker": ga_markers[2]}
                mi[prim] = {"data": local_mi, "marker": ga_markers[3]}
            else:
                me[prim] = {"data": local_me, "marker": ga_markers[4]}
                mi[prim] = {"data": local_mi, "marker": ga_markers[5]}
        else:
            me[prim] = {"data": local_me, "marker": None}
            mi[prim] = {"data": local_mi, "marker": None}



    # hs_me, hs_mi, uni, _ = MC.moments_v3(cwd / "HS_Dicts_EF", ms, 6, 6, N_value = True)
    # me["H.S."] = {"data": hs_me, "marker": None}
    # mi["H.S."] = {"data": hs_mi, "marker": None}

    MC.multiple_species_plots(ms, me, mi, uni, cwd / "TE_Images_ForPaper" / "AllSpecies", 
                              x_ticks={0.5: 0.5, 1.0: 1.0, 1.5: 1.5, 2.0: 2.0}, 
                              y_ticks={0:0, int(ymax): int(ymax)},
                              legend = True, inches = 8,
                              *args, **kwargs)


    # with open(cwd.parent / "Data_Files" / "Primates" / "Genetics" / species / f"{species}_{genome}_Frame.pkl", "rb") as file:
    #     data = pickle.load(file)
    # print(data)

def species_TE(data: pathlib.Path, *args, **kwargs):
    '''
    Since we're doing Phylogenetic trees, this is a more generalized version of the all speices plots.
    '''
    if isinstance(data, pathlib.Path):
        with open(data, "rb") as file:
            data = pickle.load(file)


    TE.time_embedding_v2(data, n = 10_000, *args, **kwargs)


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


    fly_dir = cwd / "Fly_Dicts_EF"
    mou_dir = cwd / "Mouse_Dicts_EF"
    hum_dir = cwd / "HS_Dicts_EF"

    hum_me, hum_mi, uni, _ = MC.moments_v3(hum_dir, ms, 6, 6, N_value = True)
    mou_me, mou_mi, _, _ = MC.moments_v3(mou_dir, ms, 6, 6, N_value = True)
    fly_me, fly_mi, _, _ = MC.moments_v3(fly_dir, ms, 6, 6, N_value = True)

    # me = {"Human": {"data": hum_me, "marker": "."}, "Mouse": {"data": mou_me, "marker": "+"}, "Fly": {"data": fly_me, "marker": "d"}}
    # mi = {"Human": {"data": hum_mi, "marker": "o"}, "Mouse": {"data": mou_mi, "marker": "x"}, "Fly": {"data": fly_mi, "marker": "s"}}

    me = {"Human": {"data": hum_me, "marker": None}, "Mouse": {"data": mou_me, "marker": None}, "Fly": {"data": fly_me, "marker": None}}
    mi = {"Human": {"data": hum_mi, "marker": None}, "Mouse": {"data": mou_mi, "marker": None}, "Fly": {"data": fly_mi, "marker": None}}

    MC.multiple_species_plots(ms, me, mi, uni, cwd / "TE_Images_ForPaper" / "AllSpecies", x_ticks={0.5: 0.5, 1.0: 1.0, 1.5: 1.5, 2.0: 2.0},  y_ticks = {0: 0, 7: 7})





def gen_moments(filepath: pathlib.Path, *args, **kwargs):
    '''
    Takes the heat data and generates the moments
    I really don't need a method for this... I'm not sure why I am going to stick with it hard coded like this... but whatever
    '''
    ms = [n / 100 for n in range(int(100*(0.5)), int(100*2) + 1)]

    MC.moments_v2(filepath, ms, min_k = 6, max_k = 6, N_value = True, x_ticks = {0.5: 0.5, 1.0: 1.0, 1.5: 1.5, 2.0: 2.0}, *args, **kwargs)


def gen_heat_data(filepath: pathlib.Path, exon_output: pathlib.Path = None, intron_output: pathlib.Path = None, *args, **kwargs) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
    '''
    Generates the heat data
    '''
    if isinstance(exon_output, pathlib.Path):
        exon_dir: pathlib.Path = exon_output.parent
        exon_dir.mkdir(exist_ok = True, parents = True)
    if isinstance(intron_output, pathlib.Path):
        intron_dir: pathlib.Path = intron_output.parent
        intron_dir.mkdir(exist_ok = True, parents = True)

    _, exon, intron, _, _, _, _, _, _ = heat.heat_embedding_v2(filepath, *args, **kwargs)

    # exon = pandas.DataFrame(exon)
    # exon = heat._reorder_frame(exon)
    # print("un renormalizing Exon")
    # exon = MC._unrenormalize(exon, 12, log2 = False)

    # intron = pandas.DataFrame(intron)
    # # intron = heat._reorder_frame(intron)
    # print("Un renomalizing intron")
    # intron = MC._unrenormalize(intron, 12, log2 = False)

    # exit()

    if isinstance(exon_output, pathlib.Path):
        exon.to_pickle(exon_output)
    if isinstance(intron_output, pathlib.Path):
        intron.to_pickle(intron_output)

    return exon, intron


if __name__ in "__main__":
    main()