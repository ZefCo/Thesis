from importlib.resources import path
import pandas
import pathlib
import numpy
# import re
import random
# import openpyxl
import itertools

# Instead of doing random once: do it ~ 100 times, then take the average of that number. Compare that number to the original g prime.
# Grab only the ones with a GNum < 0.5

length = 10
output_on = True

weights = numpy.array([(1 / (2**(i))) for i in range(1, length + 1)])
master_freq = ["33P", "33Ran", "55P", "55Ran"]

outputpath = pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len_{length}"



def frequency_results(input_dataframe: pandas.DataFrame, name: str):
    unique_gnum = set()
    colnames = input_dataframe.columns

    for col in colnames:
        unique_columns = set(input_dataframe[col])
        unique_gnum = unique_gnum.union(unique_columns)
    unique_gnum = sorted(unique_gnum)

    unique_gnum = pandas.Series(0, index = unique_gnum)
    # unique_stat = pandas.Series(0, index = ["Mean", "Median", "Mode", "Range", "SD", "Var"])

    for col in colnames:
        unique_columns = input_dataframe[col].value_counts()
        unique_gnum = unique_gnum.add(unique_columns, fill_value = 0)

    unique_gnum.name = name
    
    # print(unique_gnum)

    return unique_gnum



def g_number(agene, bgene):

    # weights = numpy.array([(1 / (2**(i))) for i in range(1, len(agene) + 1)])
    # print(weights)
    gdelta = numpy.array([0 if na == bgene[i] else 1 for i, na in enumerate(agene)])
    # print(gdelta)

    g = numpy.sum(weights * gdelta)

    return g


def meta_stats(dataframe: pandas.DataFrame, name: str):

    dataframe = dataframe.unstack()
    # print(dataframe)
    # print(dataframe.mean())
    meta_data = pandas.Series(data = {"Mean": dataframe.mean(), "Median": dataframe.median(), "Mode": [entry for entry in dataframe.mode()], "Max": dataframe.max(), "Min": dataframe.min(), "SD": dataframe.std(), "Var": dataframe.var()}, name = name)
    
    return meta_data






def random_comparison(agene, bgene, random_tries = 1):
    g_num = g_number(agene, bgene)
    rgs = []
    rg_primes = []
    for i in range(random_tries):
        rgene = random_seq(len(agene))
        rgs.append(g_number(agene, rgene))
        with numpy.errstate(divide='ignore'):
            rgp = numpy.where(rgs[i] != 0., g_num / rgs[i], 100)
        rg_primes.append(rgp)
        

    rgs_bar = numpy.array(rgs).mean()
    rg_primes_bar = numpy.array(rg_primes).mean()

    return g_num, rgs_bar, rg_primes_bar



def random_seq(length):
    seq = f''
    for _ in range(length):
        n = random.randint(1, 4)
        if n == 1:
            n = 'A'
        elif n == 2:
            n = 'C'
        elif n == 3:
            n = 'G'
        elif n == 4:
            n = 'T'
        seq = f'{seq}{n}'

    return seq



with open(pathlib.Path.cwd() / "Data_Files" / "UTData_test.csv") as utdata_file:
    utdata = pandas.read_csv(utdata_file, sep = ',', header = 0)

with open(pathlib.Path.cwd() / "Data_Files" / "Sequence_Files" / "HG19" / f"GeneTestSequences_len{length}.csv") as hg19seq_file:
    exon_refseq = pandas.read_csv(hg19seq_file, sep = ',', header=0)


utdata = utdata.drop_duplicates(subset=['Hgene', 'Tgene'])

# fg_rows -> the rows of the fusion genes
fg_rows, _ = utdata.shape
# hg19_rows -> the rows of the hg19 genome genes
hg19_rows, _ = exon_refseq.shape


# return_frame = pandas.DataFrame()
# return_random = pandas.DataFrame()

# FFP -> Five Prime to Five Prime
# TTP -> Three Prime to Three Prime
# R -> random

# Dictionarys for holding frames to be output to excel sheets later.
FFP_xlsx = {}
FFR_xlsx = {}
TTP_xlsx = {}
TTR_xlsx = {}
Fre_xlsx = {}
Stat_xlsx = {}

# Generates a list of all possible Dissimilarity Numbers for this length
gnums = [numpy.reshape(numpy.array(i), (1, length)) for i in itertools.product([0, 1], repeat = 1*length)]
# Creates an emtpy dataframe for the running total of frequencies for each Dissimilarity Number
zero_data = [0 for _ in range(0, 2**length)]
Frequency = pandas.DataFrame(data = {colname: zero_data for colname in master_freq}, index = [numpy.sum(g * weights) for g in gnums])


for fgrow in range(fg_rows):
    # Grabing a fusion from each row of the UT Database
    fusion_of_interest = utdata.iloc[fgrow, :].copy()

    # Grabbing the H and T gene info
    head_gene, tail_gene = fusion_of_interest['Hgene'], fusion_of_interest['Tgene']
    head_chr, tail_chr = fusion_of_interest['Hchr'], fusion_of_interest['Tchr']
    head_strand, tail_strand = fusion_of_interest['Hstrand'], fusion_of_interest['Tstrand']
    head_enst, tail_enst = fusion_of_interest['Henst'], fusion_of_interest['Tenst']

    # Creating a smaller frame from the Known Genes file based off of the data from the UT Database
    head_subframe = exon_refseq[(exon_refseq['name2'] == head_gene) & (exon_refseq['chrom'] == head_chr) & (exon_refseq['strand'] == head_strand)]
    head_subrows, _ = head_subframe.shape
    tail_subframe = exon_refseq[(exon_refseq['name2'] == tail_gene) & (exon_refseq['chrom'] == tail_chr) & (exon_refseq['strand'] == tail_strand)]
    tail_subrows, _ = tail_subframe.shape

    # Creating a Unique Identifier for the Excel Sheets to be output later
    # Possibly not going to be used: this is an Ensemble Name, which is not in the HG19 data I have
    unique_identifier = f"{head_gene}_{tail_gene}"

    if (head_subrows > 0) and (tail_subrows > 0):
        print("####")
        print(f"Unique Identifier\t{unique_identifier}")

        FFP_frame = pandas.DataFrame()
        FFP_rando = pandas.DataFrame()
        # FFP_prime = pandas.DataFrame()
        TTP_frame = pandas.DataFrame()
        TTP_rando = pandas.DataFrame()
        # TTP_prime = pandas.DataFrame()
        # print(f"Head Info\t{head_gene}\t{head_chr}\t{head_strand}")
        # print(F"Tail Info\t{tail_gene}\t{tail_chr}\t{tail_strand}")
        # print("####")
        
        for how in range(head_subrows):

            # Get the head exon count for each isoform
            hexon_count = head_subframe.iloc[how, :]['exonCount']
            
            # Generate a name list from this row
            head_gene_names = [f'{head_gene}_IF{how}_exon{e}' for e in range(1, hexon_count + 1)]
            # print(f"\t{len(head_gene_names)}")

            # This seems backwards but draw it out: the 5' is the last few nucelotides of the sequence, while the 3' is the first few.
            # Trust me: it makes sense when you draw it out.
            hexon_5prime = [head_subframe.iloc[how, :][f'exon5Prime{e}'] for e in range(1, hexon_count + 1)]
            hexon_3prime = [head_subframe.iloc[how, :][f'exon3Prime{e}_R'] for e in range(1, hexon_count + 1)]
            # print(f"\t{len(hexon_5prime)}, {len(hexon_3prime)}")

            for tow in range(tail_subrows):
                FFP_subframe = pandas.DataFrame()
                FFP_subrando = pandas.DataFrame()
                TTP_subframe = pandas.DataFrame()
                TTP_subrando = pandas.DataFrame()

                # Get the tail exon count
                texon_count = tail_subframe.iloc[tow, :]['exonCount']
                
                # Generate a name list
                tail_gene_names = [f'{tail_gene}_IF{tow}_exon{e}' for e in range(1, texon_count + 1)]
                # old_tail_gene_names = tail_gene_names
                # print(f"\t\t{len(tail_gene_names)}")

                # Getting all sequences in a list for 5' and 3'
                texon_5prime = [tail_subframe.iloc[tow, :][f'exon5Prime{e}'] for e in range(1, texon_count + 1)]
                texon_3prime = [tail_subframe.iloc[tow, :][f'exon3Prime{e}_R'] for e in range(1, texon_count + 1)]
                # print(f"\t\t{len(texon_5prime)}, {len(texon_3prime)}")


                for i, hexon5p in enumerate(hexon_5prime):
                    new_d5prime = []
                    new_d5rando = []
                    new_d3prime = []
                    new_d3rando = []

                    for j, texon5p in enumerate(texon_5prime):

                        if (len(hexon5p) != length) or (len(texon5p) != length):
                            print(f'hexon length: {len(hexon5p)}\ttexon length: {len(texon5p)}')
                            # print(f"Something weird at index {i} jndex {j}")
                            print("One of the above sequences is not of the appropriate length\nSetting G number to 1 and moving on")

                            g_5prime, g_5rando, _ = 1, 1, 1
                            g_3prime, g_3rando, _ = 1, 1, 1

                        else:
                            g_5prime, g_5rando, _ = random_comparison(texon5p, hexon5p)
                            g_3prime, g_3rando, _ = random_comparison(hexon_3prime[i], texon_3prime[j])

                        new_d5prime.append(g_5prime)
                        new_d5rando.append(g_5rando)
                        new_d3prime.append(g_3prime)
                        new_d3rando.append(g_3rando)
                    
                    new_5prime_row = pandas.DataFrame(data = dict(zip(tail_gene_names, new_d5prime)), index = pandas.Index([head_gene_names[i]]))
                    FFP_subframe = pandas.concat([FFP_subframe, new_5prime_row], axis = 0)

                    new_5rando_row = pandas.DataFrame(data = dict(zip(tail_gene_names, new_d5rando)), index = pandas.Index([head_gene_names[i]]))
                    FFP_subrando = pandas.concat([FFP_subrando, new_5rando_row], axis = 0)

                    new_3prime_row = pandas.DataFrame(data = dict(zip(tail_gene_names, new_d3prime)), index = pandas.Index([head_gene_names[i]]))
                    TTP_subframe = pandas.concat([TTP_subframe, new_3prime_row], axis = 0)

                    new_3rando_row = pandas.DataFrame(data = dict(zip(tail_gene_names, new_d3rando)), index = pandas.Index([head_gene_names[i]]))
                    TTP_subrando = pandas.concat([TTP_subrando, new_3rando_row], axis = 0)

                FFP_frame = pandas.concat([FFP_frame.stack(), FFP_subframe.stack()], axis = 0).unstack()
                TTP_frame = pandas.concat([TTP_frame.stack(), TTP_subframe.stack()], axis = 0).unstack()
                FFP_rando = pandas.concat([FFP_rando.stack(), FFP_subrando.stack()], axis = 0).unstack()
                TTP_rando = pandas.concat([TTP_rando.stack(), TTP_subrando.stack()], axis = 0).unstack()
        
        # Taking the transpose of the FFP porition, because the Tail is the second gene and I want that one as the row labels
        FFP_frame = FFP_frame.T
        FFP_rando = FFP_rando.T

        FFP_xlsx[unique_identifier] = FFP_frame
        TTP_xlsx[unique_identifier] = TTP_frame
        FFR_xlsx[f"{unique_identifier}_Random"] = FFP_rando
        TTR_xlsx[f"{unique_identifier}_Random"] = TTP_rando

        FFP_stat = meta_stats(FFP_frame, f"{unique_identifier}_55P")
        FFR_stat = meta_stats(FFP_rando, f"{unique_identifier}_55Ran")
        TTP_stat = meta_stats(TTP_frame, f"{unique_identifier}_33P")
        TTR_stat = meta_stats(TTP_rando, f"{unique_identifier}_33Ran")

        stat_frame = pandas.concat([FFP_stat.to_frame().stack(),
                                    FFR_stat.to_frame().stack(),
                                    TTP_stat.to_frame().stack(),
                                    TTR_stat.to_frame().stack()], axis = 0).unstack()

        Stat_xlsx[unique_identifier] = stat_frame

        FFP_freq = frequency_results(FFP_frame, f"{unique_identifier}_55P")
        FFR_freq = frequency_results(FFP_rando, f"{unique_identifier}_55Ran")
        TTP_freq = frequency_results(TTP_frame, f"{unique_identifier}_33P")
        TTR_freq = frequency_results(TTP_rando, f"{unique_identifier}_33Ran")

        # So this is important: the concatation order is to output the final columns ALPHABETICALLY, which means it's going to be 
        # [UniqueID]_33P    [UniqueID]_33ran    [UniqueID]_55P  [UniqueID]_33Ran
        # So even though I'm putting them in a different order they will output in this order. That's the order that has to be preserved when doing the
        # running totals.
        frequency_frame = pandas.concat([FFP_freq.to_frame().stack(), 
                                         FFR_freq.to_frame().stack(), 
                                         TTP_freq.to_frame().stack(), 
                                         TTR_freq.to_frame().stack()], axis = 0).unstack()
        # print(frequency_frame)
        # frequency_meta = pandas.concat([FFP_stat.to_frame().stack(),
        #                                 FFR_stat.to_frame().stack(),
        #                                 TTP_stat.to_frame().stack(),
        #                                 TTR_stat.to_frame().stack()], axis = 0).unstack()
        # print(frequency_meta)

        Fre_xlsx[unique_identifier] = frequency_frame

        frequency_add = frequency_frame.fillna(0)
        frequency_add.columns = master_freq
        Frequency = Frequency.add(frequency_add, fill_value = 0)
        # print(Frequency)


if output_on:
    if not outputpath.is_dir():
        outputpath.mkdir()
        # Can give a FileNotFoundError if parents don't exist: maybe add try except to create the path?
        # Can also give FileExistsError if the directory already is there, but I think this should be fine

    with pandas.ExcelWriter(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len_{length}" / f"FFP_Matrix_Length-{length}.xlsx") as FFP_writer, \
        pandas.ExcelWriter(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len_{length}" / f"TTP_Matrxi_Length-{length}.xlsx") as TTP_writer, \
        pandas.ExcelWriter(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len_{length}" / f"FFR_Matrxi_Length-{length}.xlsx") as FFR_writer, \
        pandas.ExcelWriter(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len_{length}" / f"TTR_Matrxi_Length-{length}.xlsx") as TTR_writer, \
        pandas.ExcelWriter(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len_{length}" / f"Similarity_Frequency_Length-{length}.xlsx") as Fre_writer, \
        pandas.ExcelWriter(pathlib.Path.cwd() / "Data_Files" / "GNum"/ f"Len_{length}" / f"Statistics_Length-{length}.xlsx") as Sta_writer:

        for sheetname, sheet in FFP_xlsx.items():
            FFP_xlsx[sheetname].to_excel(FFP_writer, sheetname)
            TTP_xlsx[sheetname].to_excel(TTP_writer, sheetname)
            FFR_xlsx[f"{sheetname}_Random"].to_excel(FFR_writer, f"{sheetname}_Random")
            TTR_xlsx[f"{sheetname}_Random"].to_excel(TTR_writer, f"{sheetname}_Random")
            Fre_xlsx[sheetname].to_excel(Fre_writer, sheetname)
            Stat_xlsx[sheetname].to_excel(Sta_writer, sheetname)

        Frequency.to_excel(Fre_writer, "Totals")


# Get the mean, median, mode, range, SD, and variance of all frames.