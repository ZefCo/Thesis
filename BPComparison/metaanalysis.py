# import rpy2
from re import sub
from sys import path_hooks
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


def main(file_path = pathlib.Path or str):
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

    data = data[data["CSlipLen"] >= 0]
    total_counts, _ = data.shape
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
        subgdata, subhdata, subtdata = g_freq_data[c_type], h_freq_data[c_type], t_freq_data[c_type]
        subgdano, subhdano, subtdano = g_freq_norm[c_type], h_freq_norm[c_type], t_freq_norm[c_type]

        # print(subdata)

        fig = px.histogram(subgdata, x = tuple(subgdata.index), y = c_type, title = f"{c_type}: Fusion Junction Slippage")
        fig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / f"FusJun_{c_type}.html")
        fin = px.histogram(subgdano, x = tuple(subgdano.index), y = c_type, title = f"{c_type}: Normalized Fusion Junction Slippage")
        fin.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / f"FusNorm_{c_type}.html")

        hig = px.histogram(subhdata, x = tuple(subhdata.index), y = c_type, title = f"{c_type}: Head Slippage")
        hig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / f"HeadJun_{c_type}.html")
        fig = px.histogram(subhdano, x = tuple(subhdano.index), y = c_type, title = f"{c_type}: Normalized Head Slippage")
        fig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / f"HeadNorm_{c_type}.html")

        tig = px.histogram(subtdata, x = tuple(subtdata.index), y = c_type, title = f"{c_type}: Tail Slippage")
        tig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / f"TailJun_{c_type}.html")
        fig = px.histogram(subtdano, x = tuple(subtdano.index), y = c_type, title = f"{c_type}: Normalized Tail Slippage")
        fig.write_html(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Histograms" / f"TailNorm_{c_type}.html")

        # fig.update_layout(title_text=c_type)
        # fig.show()



    # print(f"Head Slip = {hSlipCount_Gs} / {total_counts}\tTail Slip = {tSlipCount_Gs} / {total_counts}\tTotal Slip = {cSlipCount_Gs} / {total_counts}")

    # ctypes = tuple(data["Ctype"].cat.categories)



if __name__ in '__main__':
    main(file_path = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "SlippageGenes_2.csv")
