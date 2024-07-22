import ExonExonData as eed
import pathlib
cwd = pathlib.Path.cwd()
import pandas
import pickle
from Heatmaps import heatmapv2 as heatmap
from Heatmaps import _undigitize_seq as undigit
from TimeEmbedding import time_embedding_1D as te_1d
from ExonExonPlots import plots_v2
from TimeEmbedding import time_embedding_v2
from TimeEmbedding import matplotfigure
from RecoverSeq import recover_seq
import gc
from memory_profiler import profile

def limit(data: dict, x_lim: list, y_lim: list):
    '''
    Takes in the dict. This will limit the xy cords.
    The reason for doing this: instead of doing a limit on the same data set 16 times (for a total of 96 minutes)
    I can do it once then filter the data as needed each time.
    '''
    filtered_data = {}
    for key, xy in data.items():
        xy = xy[(x_lim[1] >= xy[:, 0]) & (xy[:, 0] >= x_lim[0])]
        xy = xy[(y_lim[1] >= xy[:, 1]) & (xy[:, 1] >= y_lim[0])]

        filtered_data[key] = xy

    return filtered_data

def main():
    '''
    '''
    e_data_file = cwd.parent / "Data_Files" / "Primates" / "Homo_sapiens" / "EI_Data" / "E_Data_dict.pkl"
    i_data_file = cwd.parent / "Data_Files" / "Primates" / "Homo_sapiens" / "EI_Data" / "I_Data_dict.pkl"

    k_label = 2
    steps = 4**k_label
    delta = 1 / steps
    inches = 5

    # # In case we don't want to start at the center.
    # x_start = 4*delta
    # x_stop = 5*delta
    # x_steps = steps - 4
    # y_start = 12*delta
    # y_stop = 13*delta
    # y_steps = steps - 12

    x_start = 0*delta
    x_stop = 1*delta
    x_steps = steps
    y_start = 0*delta
    y_stop = 1*delta
    y_steps = steps

    with open(e_data_file, "rb") as file:
        e_frame = pickle.load(file)

    with open(i_data_file, "rb") as file:
        i_frame = pickle.load(file)


    for i in range(0, x_steps):
        x_box = [x_start, x_stop]
        for j in range(0, y_steps):
            y_box = [y_start, y_stop]

            x_seq = recover_seq(x_box[1], k = k_label)
            y_seq = recover_seq(y_box[1], k = k_label)

            # print([x_start, x_stop], x_seq[0], [y_start, y_stop], y_seq[0])

            temp_e_frame = limit(e_frame, y_box, x_box)
            temp_i_frame = limit(i_frame, y_box, x_box)

            exon_file_path = cwd / "TE_Images_ForPaper" / "TE_Zoom" / "Exon"
            for n in x_seq[0]:
                exon_file_path = exon_file_path / f"{n}"
            exon_file_path.mkdir(parents = True, exist_ok = True)

            intron_file_path = cwd / "TE_Images_ForPaper" / "TE_Zoom" / "Intron"
            for n in x_seq[0]:
                intron_file_path = intron_file_path / f"{n}"
            intron_file_path.mkdir(parents = True, exist_ok = True)

            # print(file_path / "Exon" / f"{x_seq[0]}_{y_seq[0]}_exon.png")
            # print(file_path / "Exon" / f"{x_seq[0]}_{y_seq[0]}_intron.png")

            matplotfigure(temp_e_frame, dir = exon_file_path, file_name = f"{y_seq[0][::-1]}_{x_seq[0]}_exon.png", 
                          inches = inches, default_title = False, title = "", x_title = f"{x_seq[0]}", y_title = f"{y_seq[0][::-1]}",
                          x_tick_marks = {x: x for x in x_box}, y_tick_marks = {y: y for y in y_box}, 
                          bottom = False, left = False, labelbottom = False, labelleft = False)
            matplotfigure(temp_i_frame, dir = intron_file_path, file_name = f"{y_seq[0][::-1]}_{x_seq[0]}_intron.png", 
                          inches = inches, default_title = False, title = "", x_title = f"{x_seq[0]}", y_title = f"{y_seq[0][::-1]}",
                          x_tick_marks = {x: x for x in x_box}, y_tick_marks = {y: y for y in y_box}, 
                          bottom = False, left = False, labelbottom = False, labelleft = False)
            
            y_start = y_stop
            y_stop += delta
            del temp_e_frame
            del temp_i_frame
            gc.collect()
        else:
            y_start = 0
            y_stop = delta
            y_steps = steps
            gc.collect()

        x_start = x_stop
        x_stop += delta
        x_steps = steps
        gc.collect()


@profile
def my_func():
    a = [1] * (10 ** 6)
    b = [2] * (2 * 10 ** 7)
    del b
    return a


if __name__ in "__main__":
    main()