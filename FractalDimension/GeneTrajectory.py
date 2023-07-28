import pathlib
cwd = pathlib.Path.cwd()
import os
import re
import numpy as np
import glob
from plotly import graph_objects as go
import timeit
import DistanceClass as distance
# from scipy.spatial import distance
import pandas
import pickle
import GeneClass as Gene
import TimeEmbedding as TE



def main():
    '''
    '''

    gene: pandas.Series = import_data()
    gene = Gene.Gene(name = None,
                     ename = gene["name"],
                     gname = gene["name2"],
                     ncibname = None, 
                     chrm = gene["chrom"],
                     strand = gene["strand"],
                     txStart = gene["txStart"],
                     txEnd = gene["txEnd"],
                     cdsStart = gene["cdsStart"],
                     cdsEnd = gene["cdsEnd"],
                     exonCount = gene["exonCount"],
                     exonStarts = gene["exonStarts"],
                     exonEnds = gene["exonEnds"],
                     exonFrames = gene["exonFrames"])
    
    gene.sequence_breakdown()

    gene_plots(gene, mode = "markers+lines", kp = 6, km = 6)



def gene_plots(gene: Gene.Gene, kp = 6, km = 9, mode = "markers+lines", marker_size = 1, line_width = 0.25):
    '''
    '''

    e, i = 0, 0

    fig = go.Figure()
    # print(gene.exon_seq)
    
    for exon in gene.exon_seq:
        xy = TE.time_embedding_v3(exon, 6, 9)
        if xy is not None:
            fig.add_trace(go.Scatter(x = xy[:, 0], y = xy[:, 1], name = f"Exon {e}: Len = {len(exon)}", legendgrouptitle_text = "Exons", legendgroup = "Exons", mode = mode, marker = dict(size = marker_size), line = dict(width = line_width)))
            e += 1

    for intron in gene.intron_seq:
        xy = TE.time_embedding_v3(intron, 6, 9)
        if xy is not None:
            fig.add_trace(go.Scatter(x = xy[:, 0], y = xy[:, 1], name = f"Intron {i}: Len = {len(intron)}", legendgrouptitle_text = "Introns", legendgroup = "Introns", mode = mode, marker = dict(size = marker_size), line = dict(width = line_width)))
            i += 1


    xy = TE.time_embedding_v3(gene.full_seq[0], 6, 9)
    fig.add_trace(go.Scatter(x = xy[:, 0], y = xy[:, 1], name = f"Full_Seq", legendgrouptitle_text = "Fulle Sequence", legendgroup = "Sequence Trajectory", mode = mode, marker = dict(size = marker_size), line = dict(width = line_width)))
    fig.update_layout(legend = dict(groupclick = "toggleitem"), title = f"Gene Trajectory<br>K+ = {kp} vs K- = {km}")
    fig.show()
    fig.write_html(cwd / f"Gene_trajectory_{mode}.html")




def import_data():
    '''
    Import data and grab one gene. Don't care which one. Later I'll figure out how to do isoforms
    '''

    data_file = cwd.parent / "Data_Files" / "Gene_Files" / "Hg19" / "Known_Genes_hg19_ensGene.pkl"

    with open(data_file, "rb") as p:
        data: pandas.DataFrame = pickle.load(p)
    data = data.reset_index()

    single_gene = data.loc[200, :]

    return single_gene




if __name__ in "__main__":
    main()