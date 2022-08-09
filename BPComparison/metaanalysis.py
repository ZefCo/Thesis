# This will look at the UT database and create some heatmaps such as:


# Chr A vs Chr B
# direction A vs direction B (non chr)
# direction of Chr A vs direction of Chr B

# Also by Cancer type

import pandas
import plotly.express as px
import pathlib
import numpy
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def main():
    '''
    '''
    show = True
    min_length = 1000
    with open(pathlib.Path.cwd().parent / "Data_Files" / "UTData_cds.csv") as utdatafile:
        ut_data = pandas.read_csv(utdatafile, header=0)

    ut_data = ut_data[ut_data["SeqLen"] >= min_length]

    # chrmlist = [f'chr{i + 1}' for i in range(21)]
    ut_data = create_catagory(ut_data, "Hchr")
    ut_data = create_catagory(ut_data, "Tchr")
    ut_data = create_catagory(ut_data, "Hstrand")
    ut_data = create_catagory(ut_data, "Tstrand")

    ut_c_heatmap = heatMapGen(ut_data, "Hchr", "Tchr", diffzero=False)
    ut_c_heatmap_nosame = heatMapGen(ut_data, "Hchr", "Tchr", diffzero=True)

    chrmlist = [f"chr{i + 1}" for i in range(21)]
    chrmlist.append("chrX"), chrmlist.append("chrY")

    ut_c_heatmap = ut_c_heatmap.loc[chrmlist, chrmlist]
    ut_c_heatmap_nosame = ut_c_heatmap_nosame.loc[chrmlist, chrmlist]

    chrm_heatmap = px.imshow(ut_c_heatmap, text_auto=True, labels=dict(x="Tail Chrm", y = "Head Chrm", color = "No Normalization"))
    chrm_heatmap.update_xaxes(side = "top")
    chrm_heatmap.update_layout(title=f"Fusion Frequency across chromosomes ~ Min Seq Len = {min_length}")

    chrm_heatmap_nosame = px.imshow(ut_c_heatmap_nosame, text_auto=True, labels=dict(x="Tail Chrm", y = "Head Chrm", color = f"No Normalization"))
    chrm_heatmap_nosame.update_xaxes(side = "top")
    chrm_heatmap_nosame.update_layout(title=f"Fusion Frequency across chromosomes ~ Min Seq Len = {min_length} ~ Same Chr value = 0")

    ut_logmap = normalization(ut_c_heatmap, method="log")
    chrm_logmap = px.imshow(ut_logmap, text_auto=False, labels=dict(x="Tail Chrm", y = "Head Chrm", color = "Log 2"))
    chrm_logmap.update_xaxes(side = "top")
    chrm_logmap.update_layout(title=f"Fusion Frequency across chromosomes ~ Min Seq Len = {min_length}")

    ut_logmap_nosame = normalization(ut_c_heatmap_nosame, method="log")
    chrm_logmap_nosame = px.imshow(ut_logmap_nosame, text_auto=False, labels=dict(x="Tail Chrm", y = "Head Chrm", color = f"Log 2"))
    chrm_logmap_nosame.update_xaxes(side = "top")
    chrm_logmap_nosame.update_layout(title=f"Fusion Frequency across chromosomes ~ Min Seq Len = {min_length} ~ Same Chr value = 0")

    # ut_rangemap = normalization(ut_c_heatmap, method = "range")
    # chrm_rangemap = px.imshow(ut_rangemap, text_auto=False, labels=dict(x="Tail Chrm", y = "Head Chrm", color = "Normal to Range"))
    # chrm_rangemap.update_xaxes(side = "top")
    # chrm_rangemap.show()

    # ut_c_zmap = normalization(ut_c_heatmap, method = "zcore")
    # chrm_zmap = px.imshow(ut_c_zmap, text_auto=False, labels=dict(x="Tail Chrm", y = "Head Chrm", color = "Normal to Z-Score"))
    # chrm_zmap.update_xaxes(side = "top")
    # chrm_zmap.show()

    ut_s_heatmap = heatMapGen(ut_data, "Hstrand", "Tstrand")
    # # ut_s_heatmap = strand_count(ut_data)

    strand_heatmap = px.imshow(ut_s_heatmap, text_auto=True, labels=dict(x="Tail strand", y = "Head strand", color = "No Normalization"))
    strand_heatmap.update_xaxes(side = "top")

    # ut_s_zmap = normalization(ut_s_heatmap, method = "zcore")
    # strand_zmap = px.imshow(ut_s_zmap, text_auto=False, labels=dict(x="Tail strand", y = "Head strand", color = "Normal to Z-Score"))
    # strand_zmap.update_xaxes(side = "top")
    # strand_zmap.show()


    if show:
        strand_heatmap.show()
        strand_heatmap.write_html(pathlib.Path.cwd() / "AllStrandFrequency.html")

        chrm_logmap.show()
        chrm_logmap_nosame.show()

        chrm_logmap.write_html(pathlib.Path.cwd() / "FFxChrom_LogTransform.html")
        chrm_logmap_nosame.write_html(pathlib.Path.cwd() / "FFxChrom_LogTransform_noSame.html")

        chrm_heatmap.show()
        chrm_heatmap_nosame.show()

        chrm_heatmap.write_html(pathlib.Path.cwd() / "FFxChrom_Raw.html")
        chrm_heatmap_nosame.write_html(pathlib.Path.cwd() / "FFxChrom_Raw_noSame.html")

        strand_by_chrm(ut_data)



def create_catagory(dataframe: pandas.DataFrame, new_cat: str):
    '''
    '''

    dataframe[new_cat] = dataframe[new_cat].astype("category")

    return dataframe



def heatMapGen(dataframe: pandas.DataFrame, hcol: str, tcol: str, diffzero: bool = False):
    '''
    '''
    heads, tails = list(dataframe[hcol].value_counts().index), list(dataframe[tcol].value_counts().index)

    heatmap = pandas.DataFrame(0, columns=heads, index=tails)

    for head in heads:
        subframe: pandas.DataFrame = dataframe[dataframe[hcol] == head]
        for tail in tails:
            heatvalues = subframe[tcol].value_counts()[tail]
            if (tail == head) and diffzero:
                heatvalues = 0


            heatmap.loc[head, tail] = heatvalues

    return heatmap



def normalization(dataframe: pandas.DataFrame, method = "range", stack = True):
    '''
    For now the only method that works is setting the normal to the range.
    '''
    methods = ["range", "clipping", "log", "zcore", "sum"]
    if method not in methods:
        method = "range"
    if method in ["clipping"]:
        method = "range"
    # else:
    #     method = "range"

    if stack:
        mathframe: pandas.Series = dataframe.stack()
    else:
        mathframe: pandas.Series = dataframe.squeeze()

    if method in "range":
        min_value = mathframe.min()
        max_value = mathframe.max()
        minmax = max_value - min_value

        dataframe = (dataframe - min_value) / minmax

    elif method in "zcore":
        sd_value = mathframe.std()
        # sd_value = sd_value.std()
        mean_value = mathframe.mean()

        dataframe = (dataframe - mean_value) / sd_value

    elif method in "sum":
        mean_value = mathframe.sum()

        dataframe = dataframe / mean_value

    elif method in "log":
        dataframe = numpy.log2(dataframe + 1)

    return dataframe



def strand_by_chrm(dataframe: pandas.DataFrame):
    '''
    '''

    hchrms, tchrms = list(dataframe["Hchr"].value_counts().index), list(dataframe["Tchr"].value_counts().index)
    # hstrands, tstrands = list(dataframe["Hstrand"].value_counts().index), list(dataframe["Tstrand"].value_counts().index)
    rows, cols = 4, 6
    hchrms.sort(), tchrms.sort()



    hfig = make_subplots(rows = rows, cols = cols, subplot_titles=hchrms)

    hheatmaps = {}
    for head in hchrms:
        subframe: pandas.DataFrame = dataframe[dataframe["Hchr"] == head]
        subframe = heatMapGen(subframe, "Hstrand", "Tstrand")
        subframe = subframe.loc[["+", "-"], ["+", "-"]]

        hheatmaps[head] = subframe
        # heatmaps.append(heatMapGen(subframe, "Hstrand", "Tstrand"))

    # hindex = list(heatmaps.keys())
    # hindex.sort()

    hndex = 0
    for row in range(rows):
        for col in range(cols):
            if hndex > 23:
                break
            else:
                hfig.add_trace(go.Heatmap(x = list(hheatmaps[hchrms[hndex]].index), y = list(hheatmaps[hchrms[hndex]].columns), z=hheatmaps[hchrms[hndex]], text=hheatmaps[hchrms[hndex]], texttemplate="%{text}"), row = row + 1, col = col + 1)
            hndex += 1



    tfig = make_subplots(rows = rows, cols = cols, subplot_titles=tchrms)

    theatmaps = {}
    for tail in tchrms:
        subframe: pandas.DataFrame = dataframe[dataframe["Tchr"] == tail]
        subframe = heatMapGen(subframe, "Tstrand", "Hstrand")
        subframe = subframe.loc[["+", "-"], ["+", "-"]]

        theatmaps[tail] = subframe
        # heatmaps.append(heatMapGen(subframe, "Hstrand", "Tstrand"))

    # hindex = list(heatmaps.keys())
    # hindex.sort()

    tndex = 0
    for row in range(rows):
        for col in range(cols):
            if tndex > 23:
                break
            else:
                tfig.add_trace(go.Heatmap(x = list(theatmaps[tchrms[tndex]].index), y = list(theatmaps[tchrms[tndex]].columns), z=theatmaps[tchrms[tndex]], text=theatmaps[tchrms[tndex]], texttemplate="%{text}"), row = row + 1, col = col + 1)
            tndex += 1

    hfig.update_layout(title = "Strand Frequency by Head Chromosome")
    tfig.update_layout(title = "Strand Frequency by Tail Chromosome")
    hfig.show()
    tfig.show()

    hfig.write_html(pathlib.Path.cwd() / "SFxHChrom.html")
    hfig.write_html(pathlib.Path.cwd() / "SFxTChrom.html")



if __name__ in '__main__':
    main()