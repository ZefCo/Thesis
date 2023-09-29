import pandas
import numpy as np
import pathlib
cwd = pathlib.Path.cwd()
import re
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import itertools


def main():
    '''
    '''
    kmer = 1

    # # create_data()
    # print_data(kmer)
    # perms = nucleotide_permutations("ACGT", 2)
    # print(perms)
    window_distribution(kmer, shaded = False)
    # recreate(str(cwd / "Window_Plots" / f"{kmer}mer_Frame.csv"), kmer, labels = True)


def nucleotide_counter(sequence: str, window_size: int):
    '''
    '''
    keys: set = set()
    counter = dict()
    master_count = 0


    for i in range(len(sequence) - window_size):
        seq = sequence[i: i + window_size].upper()

        if seq not in keys:
            keys.add(seq)
            counter[seq] = 1
            master_count += 1

        else:
            counter[seq] += 1
            master_count += 1

    # for key, value in counter.items():
    #     counter[key] = value / master_count

    return counter


def create_data(kmer):
    '''
    '''
    global_data: pandas.DataFrame = pandas.read_pickle(cwd / "TrainingGeneData_v5.pkl")
    # print(global_data.shape)

    keep = np.where(global_data["Seq"].str.len() >= 100)[0]
    # print(keep)
    global_data = global_data.iloc[keep, :]

    # rows, cols = global_data.shape

    global_data = global_data.reset_index()

    # print(global_data["Type"].unique())

    # print(global_data.shape)
    # print(global_data)

    # print(global_data)

    intron, intname = 0, "Intron"
    exon, exname = 0, "Exon"

    global_mer = dict()
    exon_mer = dict()
    intron_mer = dict()

    exon_data = global_data[global_data["Type"] == exname]
    exon_data = exon_data.reset_index()
    eows, cols = exon_data.shape

    print(f"Starting Exon run\nRows = {eows}")
    for row in range(eows):
        seq_of_interest = exon_data.loc[row, "Seq"]
        nuc_count = nucleotide_counter(seq_of_interest, kmer)

        for key, value in nuc_count.items():
            if key in exon_mer.keys():
                exon_mer[key] += value
                global_mer[key] += value
            else:
                exon_mer[key] = value
                global_mer[key] = value

        # if row < 10:
            # break


    intron_data = global_data[global_data["Type"] == intname]
    intron_data = intron_data.reset_index()
    iows, cols = intron_data.shape

    print(f"Starting Intron run\nRows = {iows}")
    for row in range(iows):
        seq_of_interest = intron_data.loc[row, "Seq"]
        nuc_count = nucleotide_counter(seq_of_interest, kmer)

        for key, value in nuc_count.items():
            if key in intron_mer.keys():
                intron_mer[key] += value
            else:
                intron_mer[key] = value

            if key in global_mer.keys():
                global_mer[key] += value
            else:
                global_mer[key] = value

        # if row < 10:
            # break

    counts_frame = pandas.DataFrame([global_mer, exon_mer, intron_mer]).T
    counts_frame.columns = ["Global", "Exon", "Intron"]
    counts_frame = counts_frame.fillna(0)
    # print(counts_frame)
    # counts_frame.to_pickle(cwd / "KMerCounts.pkl")

    # print(exon_mer)
    # print(intron_mer)
    # print(global_mer)

    return counts_frame


def print_data(kmer):
    '''
    '''
    kdiff = kmer - 2
    # print(kdiff)
    count_frame = pandas.read_pickle(str(cwd / "KMerCounts.pkl"))
    count_frame = count_frame[count_frame["Global"] > 2]
    count_frame = count_frame.drop(index = ("NNNNNN"))

    count_global_sum = count_frame["Global"].sum()
    count_exon_sum = count_frame["Exon"].sum()
    count_intron_sum = count_frame["Intron"].sum()

    count_frame["EG%"] = 100 * count_frame["Exon"] / count_global_sum
    count_frame["EE%"] = 100 * count_frame["Exon"] / count_exon_sum
    count_frame["IG%"] = 100 * count_frame["Intron"] / count_global_sum
    count_frame["II%"] = 100 * count_frame["Intron"] / count_intron_sum

    # print(count_frame)
    windows = list(count_frame.index)

    perms5, perms3 = nucleotide_permutations("ACGT", 2), nucleotide_permutations("ACGT", 2)
    eerms5, eerms3 = nucleotide_permutations("ACGT", 2), nucleotide_permutations("ACGT", 2)
    ierms5, ierms3 = nucleotide_permutations("ACGT", 2), nucleotide_permutations("ACGT", 2)

    for window in windows:
        for key in perms5.keys():
            if re.match(fr"{key}\w{{{kdiff}}}", window):
                perms5[key] += count_frame.loc[window, "Global"]
                eerms5[key] += count_frame.loc[window, "Exon"]
                ierms5[key] += count_frame.loc[window, "Intron"]

            if re.match(fr"\w{{{kdiff}}}{key}", window):
                perms3[key] += count_frame.loc[window, "Global"]
                eerms3[key] += count_frame.loc[window, "Exon"]
                ierms3[key] += count_frame.loc[window, "Intron"]


    for key in perms5.keys():
        perms5[key] = perms5[key] / count_global_sum
        eerms5[key] = eerms5[key] / count_exon_sum
        ierms5[key] = ierms5[key] / count_intron_sum
        
        perms3[key] = perms3[key] / count_global_sum
        eerms3[key] = eerms3[key] / count_exon_sum
        ierms3[key] = ierms3[key] / count_intron_sum


    # cannonical_5g, cannonical_5e, cannonical_5i = 0, 0, 0
    # cannonical_3g, cannonical_3e, cannonical_3i = 0, 0, 0

    # for window in windows:
    #     if re.match(fr"GT\w{{{kdiff}}}", window):
    #         cannonical_5g += count_frame.loc[window, "Global"]
    #         cannonical_5e += count_frame.loc[window, "Exon"]
    #         cannonical_5i += count_frame.loc[window, "Intron"]

    #     if re.match(fr"\w{{{kdiff}}}AG", window):
    #         cannonical_3g += count_frame.loc[window, "Global"]
    #         cannonical_3e += count_frame.loc[window, "Exon"]
    #         cannonical_3i += count_frame.loc[window, "Intron"]

    # c5g, c5e, c5i = cannonical_5g / count_global_sum, cannonical_5e / count_exon_sum, cannonical_5i / count_intron_sum
    # c3g, c3e, c3i = cannonical_3g / count_global_sum, cannonical_3e / count_exon_sum, cannonical_3i / count_intron_sum

    # print(f"Global occurance of Cannonical Splicesite:\n\t5' = {cannonical_5g} - {c5g}\t3' = {cannonical_3g} - {c3g}")
    # print(f"Global occurance of Cannonical Splicesite in Exons:\n\t5' = {cannonical_5e} - {c5e}\t3' = {cannonical_3e} - {c3e}")
    # print(f"Global occurance of Cannonical Splicesite in Introns:\n\t5' = {cannonical_5i} - {c5i}\t3' = {cannonical_3i} - {c3i}")

    aveglo5 = sum(perms5.values()) / len(perms5.keys())
    maxglo5 = max(perms5.values())
    minglo5 = min(perms5.values())

    aveexo5 = sum(eerms5.values()) / len(eerms5.keys())
    maxexo5 = max(eerms5.values())
    minexo5 = min(eerms5.values())

    aveint5 = sum(ierms5.values()) / len(ierms5.keys())
    maxint5 = max(ierms5.values())
    minint5 = min(ierms5.values())

    aveglo3 = sum(perms3.values()) / len(perms3.keys())
    maxglo3 = max(perms3.values())
    minglo3 = min(perms3.values())

    aveexo3 = sum(eerms3.values()) / len(eerms3.keys())
    maxexo3 = max(eerms3.values())
    minexo3 = min(eerms3.values())

    aveint3 = sum(ierms3.values()) / len(ierms3.keys())
    maxint3 = max(ierms3.values())
    minint3 = min(ierms3.values())

    x = ["GT-NNNN", "Max XX-NNNN 5'", "Ave XX-NNNN 5'", "Min XX-NNNN 5'", "NNNN-AG", "Max NNNN-XX 3'", "Ave NNNN-XX 3'", "Min NNNN-XX 3'"]
    y1 = [perms5["GT"], maxglo5, aveglo5, minglo5, perms3["AG"], maxglo3, aveglo3, minglo3]
    y2 = [eerms5["GT"], maxexo5, aveexo5, minexo5, eerms3["AG"], maxexo3, aveexo3, minexo3]
    y3 = [ierms5["GT"], maxint5, aveint5, minint5, ierms3["AG"], maxint3, aveint3, minint3]

    fig = go.Figure()
    fig.add_trace(go.Bar(x = x, y = y1, name = "Global", text = y1))
    fig.add_trace(go.Bar(x = x, y = y2, name = "Exon", text = y2))
    fig.add_trace(go.Bar(x = x, y = y3, name = "Intron", text = y3))
    fig.add_hline(y = 0.25*0.25, line_dash = "dot", annotation_text = "Random Occurance", annotation_position = "bottom right")  #, 
    fig.add_vrect(x0 = -0.5, x1 = 0.5, col = "all", fillcolor = "red", opacity = 0.25, annotation_text = "Cannonical 5'")
    fig.add_vrect(x0 = 3.5, x1 = 4.5, col = "all", fillcolor = "green", opacity = 0.25, annotation_text = "Cannonical 3'")
    fig.update_layout(barmode = "group", title = "Occurences of 6-Mer window Sequences")
    # fig.update_traces(opacity = 0.75)
    fig.show()

    # should eventually compare this to all other combinations and display their average, min, and max

    # print(perms5)
    # print(perms3)
    # print(eerms5)
    # print(eerms3)
    # print(ierms5)
    # print(ierms3)



def nucleotide_permutations(sequence: str = "ACGT", length: int = 3) -> dict:
    nuc_perm = dict()

    if len(sequence) < length:
        return None

    perms = itertools.product(sequence, repeat=length)
    for p in perms:

        key = ""
        for n in p:
            key = f"{key}{n}"

        nuc_perm[key] = 0

    return nuc_perm


def window_distribution(kmer, shaded = True):
    '''
    '''

    average = 0.25**kmer

    window_plots_dir = cwd / "Window_Plots"
    window_plots_dir.mkdir(parents = True, exist_ok = True)

    data = create_data(kmer)

    motiffs = list(data.index)
    for motiff in motiffs:
        if re.search("N", motiff):
            data = data.drop(index = (motiff))

    count_global_sum = data["Global"].sum()
    count_exon_sum = data["Exon"].sum()
    count_intron_sum = data["Intron"].sum()

    data["G%"] = (data["Global"] / count_global_sum) - average
    data["E%"] = (data["Exon"] / count_exon_sum) - average
    data["I%"] = (data["Intron"] / count_intron_sum) - average

    data = data.sort_index()

    data.to_csv(str(window_plots_dir / f"{kmer}mer_Frame.csv"), header = True)

    fig = go.Figure()
    fig.add_trace(go.Bar(x = list(data.index), y = data["G%"], name = "Full Seq"))
    fig.add_trace(go.Bar(x = list(data.index), y = data["E%"], name = "Exon"))
    fig.add_trace(go.Bar(x = list(data.index), y = data["I%"], name = "Intron"))

    # fig = make_subplots(rows = 3, cols = 1)
    # fig.add_trace(go.Bar(x = list(data.index), y = data["G%"], name = "Global"), row = 1, col = 1)
    # fig.add_trace(go.Bar(x = list(data.index), y = data["E%"], name = "Exon"), row = 2, col = 1)
    # fig.add_trace(go.Bar(x = list(data.index), y = data["I%"], name = "Intron"), row = 3, col = 1)

    if shaded:
        if kmer > 1:
            sections = (4**kmer)/4
            fig.add_vrect(x0 =            - 0.5, x1 =   sections - 0.5, col = "all", fillcolor = "red",    opacity = 0.25, annotation_text = "A")
            fig.add_vrect(x0 =   sections - 0.5, x1 = 2*sections - 0.5, col = "all", fillcolor = "green",  opacity = 0.25, annotation_text = "C")
            fig.add_vrect(x0 = 2*sections - 0.5, x1 = 3*sections - 0.5, col = "all", fillcolor = "blue",   opacity = 0.25, annotation_text = "G")
            fig.add_vrect(x0 = 3*sections - 0.5, x1 = 4*sections - 0.5, col = "all", fillcolor = "yellow", opacity = 0.25, annotation_text = "T")


    fig.update_layout(barmode = "overlay", title = f"Occurences of {kmer}-Mer window Sequences <br> Subtracted by a mean of {average}")
    fig.update_traces(opacity = 0.9)
    fig.update_xaxes(showticklabels = False)
    fig.show()
    fig.write_html(str(window_plots_dir / f"{kmer}mer_Bars.html"))


def recreate(filepath: pathlib.Path, kmer, shaded = False, labels = False):
    '''
    '''
    # data.to_csv(str(window_plots_dir / f"{kmer}mer_Frame.csv"), header = True)
    data = pandas.read_csv(filepath, header = 0)
    # print(data)

    average = 0.25**kmer


    fig = go.Figure()
    fig.add_trace(go.Bar(x = data["Unnamed: 0"], y = data["G%"], name = "Full Seq"))
    fig.add_trace(go.Bar(x = data["Unnamed: 0"], y = data["E%"], name = "Exon"))
    fig.add_trace(go.Bar(x = data["Unnamed: 0"], y = data["I%"], name = "Intron"))

    # fig = make_subplots(rows = 3, cols = 1)
    # fig.add_trace(go.Bar(x = list(data.index), y = data["G%"], name = "Global"), row = 1, col = 1)
    # fig.add_trace(go.Bar(x = list(data.index), y = data["E%"], name = "Exon"), row = 2, col = 1)
    # fig.add_trace(go.Bar(x = list(data.index), y = data["I%"], name = "Intron"), row = 3, col = 1)

    if shaded:
        if kmer > 1:
            sections = (4**kmer)/4
            fig.add_vrect(x0 =            - 0.5, x1 =   sections - 0.5, col = "all", fillcolor = "red",    opacity = 0.25, annotation_text = "A")
            fig.add_vrect(x0 =   sections - 0.5, x1 = 2*sections - 0.5, col = "all", fillcolor = "green",  opacity = 0.25, annotation_text = "C")
            fig.add_vrect(x0 = 2*sections - 0.5, x1 = 3*sections - 0.5, col = "all", fillcolor = "blue",   opacity = 0.25, annotation_text = "G")
            fig.add_vrect(x0 = 3*sections - 0.5, x1 = 4*sections - 0.5, col = "all", fillcolor = "yellow", opacity = 0.25, annotation_text = "T")


    fig.update_layout(barmode = "overlay", title = f"Occurences of {kmer}-Mer window Sequences <br> Subtracted by a mean of {average}")
    fig.update_traces(opacity = 0.9)
    fig.update_xaxes(showticklabels = labels)
    fig.show()
    # fig.write_html(str(window_plots_dir / f"{kmer}mer_Bars.html"))



if __name__ in "__main__":
    main()



