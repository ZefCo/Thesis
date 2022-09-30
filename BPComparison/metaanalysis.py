# import rpy2
import pandas
import pathlib
# from RFunClass import RFun
import plotly.express as px


# want min
# want max
# qaurtile of data
# Done by each column
# A lot of this can be done in ggplot

# import data to R
# apply catagories to Classification
# For each classification + global
    # y = column

# For Global Seq < 10
# Global count of how many have seq from head
# Global count of how many have seq from tail
# Global count of how many have seq from both
    # Do this by cancer too


def score_histogram(file_path = pathlib.Path or str):
    '''
    '''
    if isinstance(file_path, str):
        file_path = pathlib.Path(file_path)
    
    # import file in different ways, depending on if its csv or xlx
    if file_path.suffix in ".csv":
        with open(file_path) as indata:
            data = pandas.read_csv(indata, header=0)
    elif file_path.suffix in ".xlsx":
        with open(file_path) as indata:
            data = pandas.read_excel(indata, sheet_name = "idk", header = 0)

    # grab only the stuff between [1, 10] length. This data is copy and pasted from the Slip Seq File (same order)
    data = data[(data["CSlipLen"] < 11) & (data["CSlipLen"] > 0)]

    # empty dataframes to hold the counts later
    ht3Intron_freq_data = pandas.DataFrame()
    ht3Exon_freq_data = pandas.DataFrame()
    th5Intron_freq_data = pandas.DataFrame()
    th5Exon_freq_data = pandas.DataFrame()

    # Being a little lazy here: I'm going to create them as catagories, that way I can just easily count the number of catagories
    data["H3_Seq_T3_Intron"], data["H3_Seq_T3_Exon"], data["H5_Intron_T5_Seq"], data["H5_Exon_T5_Seq"], data["Ctype"] = pandas.Categorical(data["H3_Seq_T3_Intron"]), pandas.Categorical(data["H3_Seq_T3_Exon"]), pandas.Categorical(data["H5_Intron_T5_Seq"]), pandas.Categorical(data["H5_Exon_T5_Seq"]), pandas.Categorical(data["Ctype"])

    # grab the cancer types, include a global option so we can look at everything
    c_types = tuple(["Global"] + list(data["Ctype"].cat.categories))

    # Iterate through the types of cancer, counting the frequency of the scores
    for c_type in c_types:
        if c_type in "Global":
            ht3Intron_score_freq = data["H3_Seq_T3_Intron"].value_counts()
            ht3Intron_score_freq.name = c_type
            ht3Intron_freq_data = pandas.concat([ht3Intron_freq_data, ht3Intron_score_freq.to_frame().T])

            ht3Exon_score_freq = data["H3_Seq_T3_Exon"].value_counts()
            ht3Exon_score_freq.name = c_type
            ht3Exon_freq_data = pandas.concat([ht3Exon_freq_data, ht3Exon_score_freq.to_frame().T])

            th5Intron_score_freq = data["H5_Intron_T5_Seq"].value_counts()
            th5Intron_score_freq.name = c_type
            th5Intron_freq_data = pandas.concat([th5Intron_freq_data, th5Intron_score_freq.to_frame().T])

            th5Exon_score_freq = data["H5_Exon_T5_Seq"].value_counts()
            th5Exon_score_freq.name = c_type
            th5Exon_freq_data = pandas.concat([th5Exon_freq_data, th5Exon_score_freq.to_frame().T])

        else:
            subdata = data[data["Ctype"] == c_type]
            
            ht3Intron_score_freq = subdata["H3_Seq_T3_Intron"].value_counts()
            ht3Intron_score_freq.name = c_type
            ht3Intron_freq_data = pandas.concat([ht3Intron_freq_data, ht3Intron_score_freq.to_frame().T])

            ht3Exon_score_freq = subdata["H3_Seq_T3_Exon"].value_counts()
            ht3Exon_score_freq.name = c_type
            ht3Exon_freq_data = pandas.concat([ht3Exon_freq_data, ht3Exon_score_freq.to_frame().T])

            th5Intron_score_freq = subdata["H5_Intron_T5_Seq"].value_counts()
            th5Intron_score_freq.name = c_type
            th5Intron_freq_data = pandas.concat([th5Intron_freq_data, th5Intron_score_freq.to_frame().T])

            th5Exon_score_freq = subdata["H5_Exon_T5_Seq"].value_counts()
            th5Exon_score_freq.name = c_type
            th5Exon_freq_data = pandas.concat([th5Exon_freq_data, th5Exon_score_freq.to_frame().T])


    ht3Intron_freq_data = ht3Intron_freq_data.T
    ht3Intron_norm_data = ht3Intron_freq_data / ht3Intron_freq_data.sum(axis = 0)
    ht3Intron_freq_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "HT3_Intron_Data.csv")
    ht3Intron_norm_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "HT3_Intron_Norm.csv")

    ht3Exon_freq_data = ht3Exon_freq_data.T
    ht3Exon_norm_data = ht3Exon_freq_data / ht3Exon_freq_data.sum(axis = 0)
    ht3Exon_freq_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "HT3_Exon_Data.csv")
    ht3Exon_norm_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "HT3_Exon_Norm.csv")
    
    th5Intron_freq_data = th5Intron_freq_data.T
    th5Intron_norm_data = th5Intron_freq_data / th5Intron_freq_data.sum(axis = 0)
    th5Intron_freq_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "TH5_Intron_Data.csv")
    th5Intron_norm_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "TH5_Intron_Norm.csv")
    
    th5Exon_freq_data = th5Exon_freq_data.T
    th5Exon_norm_data = th5Exon_freq_data / th5Exon_freq_data.sum(axis = 0)
    th5Exon_freq_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "TH5_Exon_Data.csv")
    th5Exon_norm_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "TH5_Exon_Norm.csv")

    for c_type in c_types:
        # subdata = pandas.DataFrame(g_freq_data[c_type])
        subhidata, subhedata, subtidata, subtedata = ht3Intron_freq_data[c_type], ht3Exon_freq_data[c_type], th5Intron_freq_data[c_type], th5Exon_freq_data[c_type]
        subhidano, subhedano, subtidano, subtedano = ht3Intron_norm_data[c_type], ht3Exon_norm_data[c_type], th5Intron_norm_data[c_type], th5Exon_norm_data[c_type]

        # print(subdata)

        hig = px.histogram(subhidata, x = tuple(subhidata.index), y = c_type, nbins = len(tuple(subhidata.index)), title = f"{c_type}: Head 3 Prime to Tail 3 Prime Intron Scores")
        hig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Scores" / f"HT3I_{c_type}.html")
        hin = px.histogram(subhidano, x = tuple(subhidano.index), y = c_type, nbins = len(tuple(subhidata.index)), title = f"{c_type}: Normalized Head 3 Prime to Tail 3 Prime Intron Scores")
        hin.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Scores" / f"HT3I_{c_type}_norm.html")

        heg = px.histogram(subhedata, x = tuple(subhedata.index), y = c_type, nbins = len(tuple(subhedata.index)), title = f"{c_type}: Head 3 Prime to Tail 3 Prime Exon Scores")
        heg.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Scores" / f"HT3E_{c_type}.html")
        hen = px.histogram(subhedano, x = tuple(subhedano.index), y = c_type, nbins = len(tuple(subhedano.index)), title = f"{c_type}: Normalized Head 3 Prime to Tail 3 Prime Exon Scores")
        hen.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Scores" / f"HT3E_{c_type}_norm.html")

        tig = px.histogram(subtidata, x = tuple(subtidata.index), y = c_type, nbins = len(tuple(subtidata.index)), title = f"{c_type}: Tail 5 Prime to Head 5 Prime Intron Scores")
        tig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Scores" / f"TH5I_{c_type}.html")
        tin = px.histogram(subtidano, x = tuple(subtidano.index), y = c_type, nbins = len(tuple(subtidano.index)), title = f"{c_type}: Normalized Tail 5 Prime to Head 5 Prime Intron Scores")
        tin.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Scores" / f"TH5I_{c_type}_norm.html")

        teg = px.histogram(subtedata, x = tuple(subtedata.index), y = c_type, nbins = len(tuple(subtedata.index)), title = f"{c_type}: Tail 5 Prime to Head 5 Prime Exon Scores")
        teg.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Scores" / f"TH5E_{c_type}.html")
        ten = px.histogram(subtedano, x = tuple(subtedano.index), y = c_type, nbins = len(tuple(subtedano.index)), title = f"{c_type}: Normalized Tail 5 Prime to Head 5 Prime Exon Scores")
        ten.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Scores" / f"TH5E_{c_type}_norm.html")



def seq_histogram(file_path = pathlib.Path or str):
    '''
    '''
    if isinstance(file_path, str):
        file_path = pathlib.Path(file_path)

    if file_path.suffix in ".csv":
        with open(file_path) as indata:
            data = pandas.read_csv(indata, header=0)
    elif file_path.suffix in ".xlsx":
        with open(file_path) as indata:
            data = pandas.read_excel(indata, sheet_name = "SlippageGenes_2", header = 0)

    # data = data[data["CSlipLen"] >= 0]
    # print(total_count)
    data = data[(data["CSlipLen"] < 11) & (data["CSlipLen"] > 0)]

    data["HSlip"], data["TSlip"], data["CSlip"], data["Ctype"] = pandas.Categorical(data["HSlip"]), pandas.Categorical(data["TSlip"]), pandas.Categorical(data["CSlip"]), pandas.Categorical(data["Ctype"])

    # hSeqs_Gs = tuple(data["HSlip"].cat.categories)
    # tSeqs_Gs = tuple(data["TSlip"].cat.categories)
    # cSeqs_Gs = tuple(data["CSlip"].cat.categories)

    # hSlipCount_Gs, _ = data[data["HSlipLen"] > 0].shape
    # tSlipCount_Gs, _ = data[data["TSlipLen"] > 0].shape
    # cSlipCount_Gs, _ = data[data["CSlipLen"] > 0].shape

    c_types = tuple(["Global"] + list(data["Ctype"].cat.categories))

    g_freq_data = pandas.DataFrame()
    h_freq_data = pandas.DataFrame()
    t_freq_data = pandas.DataFrame()

    for c_type in c_types:
        if c_type in "Global":
            g_seq_freq = data["CSlip"].value_counts()
            g_seq_freq.name = "Global"
            g_freq_data = pandas.concat([g_freq_data, g_seq_freq.to_frame().T])

            h_seq_freq = data["HSlip"].value_counts()
            h_seq_freq.name = "Global"
            h_freq_data = pandas.concat([h_freq_data, h_seq_freq.to_frame().T])
            
            t_seq_freq = data["TSlip"].value_counts()
            t_seq_freq.name = "Global"
            t_freq_data = pandas.concat([t_freq_data, t_seq_freq.to_frame().T])

        else:
            subdata = data[data["Ctype"] == c_type]

            g_seq_freq = subdata["CSlip"].value_counts()
            g_seq_freq.name = c_type
            g_freq_data = pandas.concat([g_freq_data, g_seq_freq.to_frame().T])
            
            h_seq_freq = subdata["HSlip"].value_counts()
            h_seq_freq.name = c_type
            h_freq_data = pandas.concat([h_freq_data, h_seq_freq.to_frame().T])
            
            t_seq_freq = subdata["TSlip"].value_counts()
            t_seq_freq.name = c_type
            t_freq_data = pandas.concat([t_freq_data, t_seq_freq.to_frame().T])

    g_freq_data = g_freq_data.T
    g_freq_norm = g_freq_data / g_freq_data.sum(axis = 0)
    g_freq_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "GFreqData.csv")
    g_freq_norm.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "GFreqNorm.csv")

    h_freq_data = h_freq_data.T
    h_freq_norm = h_freq_data / h_freq_data.sum(axis = 0)
    h_freq_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "HFreqData.csv")
    h_freq_norm.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "HFreqNorm.csv")


    t_freq_data = t_freq_data.T
    t_freq_norm = t_freq_data / t_freq_data.sum(axis = 0)
    t_freq_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "TFreqData.csv")
    t_freq_norm.to_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "HistogramData" / "TFreqNorm.csv")
    
    # OK Yes I am doing the same thing a second time
    c_types = tuple(g_freq_data.columns)

    for c_type in c_types:
        # subdata = pandas.DataFrame(g_freq_data[c_type])
        # subgdata, subhdata, subtdata = g_freq_data[c_type].nlargest(n = 10), h_freq_data[c_type].nlargest(n = 10), t_freq_data[c_type].nlargest(n = 10)
        subgdata, subhdata, subtdata = g_freq_data[c_type], h_freq_data[c_type], t_freq_data[c_type]
        subgdata, subhdata, subtdata = subgdata.loc[lambda x: x / subgdata.sum() >= 0.01], subhdata.loc[lambda x: x / subhdata.sum() >= 0.01], subtdata.loc[lambda x: x / subhdata.sum() >= 0.01]

        # subgdata, subhdata, subtdata = subgdata.loc[lambda x: x >=5], subhdata.loc[lambda x: x >=5], subtdata.loc[lambda x: x >=5]
        subgdano, subhdano, subtdano = g_freq_norm[c_type].nlargest(n = 10), h_freq_norm[c_type].nlargest(n = 10), t_freq_norm[c_type].nlargest(n = 10)

        # print(type(subgdata))

        # print(subdata)

        fig = px.histogram(subgdata, x = tuple(subgdata.index), y = c_type, title = f"{c_type}: Fusion Junction Slippage", color = tuple(subgdata.index))
        fig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Seq" / f"FusJun_{c_type}.html")
        fin = px.histogram(subgdano, x = tuple(subgdano.index), y = c_type, title = f"{c_type}: Normalized Fusion Junction Slippage", color = tuple(subgdano.index))
        fin.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Seq" / f"FusNorm_{c_type}.html")

        hig = px.histogram(subhdata, x = tuple(subhdata.index), y = c_type, title = f"{c_type}: Head Slippage", color = tuple(subhdata.index))
        hig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Seq" / f"HeadJun_{c_type}.html")
        hin = px.histogram(subhdano, x = tuple(subhdano.index), y = c_type, title = f"{c_type}: Normalized Head Slippage", color = tuple(subhdano.index))
        hin.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Seq" / f"HeadNorm_{c_type}.html")

        tig = px.histogram(subtdata, x = tuple(subtdata.index), y = c_type, title = f"{c_type}: Tail Slippage", color = tuple(subtdata.index))
        tig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Seq" / f"TailJun_{c_type}.html")
        tin = px.histogram(subtdano, x = tuple(subtdano.index), y = c_type, title = f"{c_type}: Normalized Tail Slippage", color = tuple(subtdano.index))
        tin.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / "Seq" / f"TailNorm_{c_type}.html")

        fig.update_layout(title_text=c_type)
        # fig.show()
        # hig.show()
        # tig.show()



    # print(f"Head Slip = {hSlipCount_Gs} / {total_counts}\tTail Slip = {tSlipCount_Gs} / {total_counts}\tTotal Slip = {cSlipCount_Gs} / {total_counts}")

    # ctypes = tuple(data["Ctype"].cat.categories)



if __name__ in '__main__':
    seq_histogram(file_path = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "SlippageGenes_2.csv")
    # score_histogram(file_path = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Scoring_min100_83022_15k.csv")
