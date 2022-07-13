import pandas
import pathlib
import numpy
import re
import random

# Instead of doing random once: do it ~ 100 times, then take the average of that number. Compare that number to the original g prime.
# Grab only the ones with a GNum < 0.5

length = 5


def frequency_results(input_dataframe: pandas.DataFrame):
    unique_gnum = set()
    colnames = input_dataframe.columns

    for col in colnames:
        unique_columns = set(input_dataframe[col])
        unique_gnum = unique_gnum.union(unique_columns)
    unique_gnum = sorted(unique_gnum)

    unique_gnum = pandas.Series(0, index = unique_gnum)

    for col in colnames:
        unique_columns = input_dataframe[col].value_counts()
        unique_gnum = unique_gnum.add(unique_columns, fill_value = 0)

    return(unique_gnum)



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


def g_number(agene, bgene):

    weights = numpy.array([(1 / (2**(i))) for i in range(1, len(agene) + 1)])
    # print(weights)
    gdelta = numpy.array([0 if na == bgene[i] else 1 for i, na in enumerate(agene)])
    # print(gdelta)

    g = numpy.sum(weights * gdelta)

    return g


def random_comparison(agene, bgene, random_tries = 1):
    g_num = g_number(agene, bgene)
    rgs = []
    rg_primes = []
    for i in range(random_tries):
        rgene = random_seq(len(agene))
        rgs.append(g_number(agene, rgene))
        rg_primes.append(g_num / rgs[i])

    rgs_bar = numpy.array(rgs).mean()
    rg_primes_bar = numpy.array(rg_primes).mean()

    return g_num, rgs_bar, rg_primes_bar




with open(pathlib.Path.cwd() / "Data_Files" / "UTData_test.csv") as utdata_file:
    utdata = pandas.read_csv(utdata_file, sep = ',', header = 0)

with open(pathlib.Path.cwd() / "Data_Files" / "Sequence_Files" / "HG19" / f"GeneTest_0_58_len{length}.csv") as hg19seq_file:
    exon_refseq = pandas.read_csv(hg19seq_file, sep = ',', header=0)


fg_rows, _ = utdata.shape
hg19_rows, _ = exon_refseq.shape


# return_frame = pandas.DataFrame()
# return_random = pandas.DataFrame()

HH_frame = pandas.DataFrame()
HH_random = pandas.DataFrame()
HH_prime = pandas.DataFrame()
TT_frame = pandas.DataFrame()
TT_random = pandas.DataFrame()
TT_prime = pandas.DataFrame()

HH_column_names = []
TT_column_names = []

for fgrow in range(fg_rows):
    fusion_of_interest = utdata.iloc[fgrow, :].copy()
    head_gene, tail_gene = fusion_of_interest['Hgene'], fusion_of_interest['Tgene']
    head_chr, tail_chr = fusion_of_interest['Hchr'], fusion_of_interest['Tchr']
    head_strand, tail_strand = fusion_of_interest['Hstrand'], fusion_of_interest['Tstrand']

    head_subframe = exon_refseq[(exon_refseq['name2'] == head_gene) & (exon_refseq['chrom'] == head_chr) & (exon_refseq['strand'] == head_strand)]
    head_subrows, _ = head_subframe.shape
    tail_subframe = exon_refseq[(exon_refseq['name2'] == tail_gene) & (exon_refseq['chrom'] == tail_chr) & (exon_refseq['strand'] == tail_strand)]
    tail_subrows, _ = tail_subframe.shape

    for how in range(head_subrows):
        hexon_count = head_subframe.iloc[how, :]['exonCount']
        
        head_gene_names = [f'{head_gene}_{how}_exon{e}' for e in range(1, hexon_count + 1)]
        TT_column_names = TT_column_names + head_gene_names

        # This seems backwards but draw it out: the 5' is the last few nucelotides of the sequence, while the 3' is the first few.
        # Trust me: it makes sense when you draw it out.
        hexon_5prime = [head_subframe.iloc[how, :][f'exonEnd{e}_R'] for e in range(1, hexon_count + 1)]
        hexon_3prime = [head_subframe.iloc[how, :][f'exonStart{e}'] for e in range(1, hexon_count + 1)]

        for tow in range(tail_subrows):
            texon_count = tail_subframe.iloc[tow, :]['exonCount']
            
            tail_gene_names = [f'{tail_gene}_{tow}_exon{e}' for e in range(1, texon_count + 1)]
            # measured_3prime_column_name = [f'{tail_gene}_{tow}_exon{e}' for e in range(1, texon_count + 1)]
            HH_column_names = HH_column_names + tail_gene_names

            texon_5prime = [tail_subframe.iloc[tow, :][f'exonEnd{e}_R'] for e in range(1, texon_count + 1)]
            texon_3prime = [tail_subframe.iloc[tow, :][f'exonStart{e}'] for e in range(1, texon_count + 1)]

            for i, hexon in enumerate(hexon_5prime):
                for j, texon in enumerate(texon_5prime):

                    if (len(hexon) != length) or (len(texon) != length):
                        print(f'hexon length: {len(hexon)}\ttexon length: {len(texon)}')
                        print(f"Something weird at index {i} jndex {j}")

                    # g_5prime = g_number(hexon, texon)
                    # rexon = random_seq(len(hexon))
                    # g_5rando = g_number(hexon, rexon)
                    g_5prime, g_5rando, g_5prime_bar = random_comparison(hexon, texon)
                    HH_frame.loc[head_gene_names[i], tail_gene_names[j]] = g_5prime
                    HH_random.loc[head_gene_names[i], tail_gene_names[j]] = g_5rando
                    HH_prime.loc[head_gene_names[i], tail_gene_names[j]] = g_5prime_bar

                    # print(f'H5Exon: {hexon} -> T5Exon: {texon} = {g_5prime} & R5Exon: {rexon} = {g_5rando}')

                    # g_5prime = g_number(hexon_3prime[i], texon_3prime[j])
                    # r3exon = random_seq(len(hexon_3prime[i]))
                    # g_5rando = g_number(texon_3prime[i], r3exon)
                    g_3prime, g_3rando, g_3prime_bar = random_comparison(texon_3prime[j], hexon_3prime[i])
                    TT_frame.loc[tail_gene_names[j], head_gene_names[i]] = g_3prime
                    TT_random.loc[tail_gene_names[j], head_gene_names[i]] = g_3rando
                    TT_prime.loc[tail_gene_names[j], head_gene_names[i]] = g_3prime_bar


                    # print(f'T3Exon: {texon_3prime[i]} -> H3Exon: {hexon_3prime[i]} = {g_5prime} & R3Exon: {r3exon} = {g_5rando}')

HH_frame = HH_frame.fillna(1)
HH_Frame_frequency = frequency_results(HH_frame)
HH_random = HH_random.fillna(1)
HH_Random_frequency = frequency_results(HH_random)
TT_frame = TT_frame.fillna(1)
TT_Frame_frequency = frequency_results(TT_frame)
TT_random = TT_random.fillna(1)
TT_Random_frequency = frequency_results(TT_random)


HH_prime = HH_prime.fillna(100)
HH_prime.replace([numpy.inf, -numpy.inf], 100, inplace=True)
HH_prime_frequency = frequency_results(HH_prime)
# print(HH_prime.max())
# print(HH_prime_frequency.max())

TT_prime = TT_prime.fillna(100)
TT_prime.replace([numpy.inf, -numpy.inf], 100, inplace=True)
TT_prime_frequency = frequency_results(TT_prime)


HH_frame.to_csv(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len{length}" / "Head-to-Head_test.csv")
HH_random.to_csv(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len{length}" / "Head-to-Head_random_test.csv")
TT_frame.to_csv(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len{length}" / "Tail-to-Tail_test.csv")
TT_frame.to_csv(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len{length}" / "Tail-to-Tail_random_test.csv")

HH_prime.to_csv(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len{length}" / "Head-to-Head_prime.csv")
TT_prime.to_csv(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len{length}" / "Tail-to-Tail_prime.csv")

data = {"Head2Head": HH_Frame_frequency, "Head2Random": HH_Random_frequency, "Tail2Tail": TT_Frame_frequency, "Tail2Random": TT_Random_frequency}
output_frequency = pandas.concat(data, axis = 1)
output_frequency = output_frequency.sort_index()
output_frequency.to_csv(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len{length}" / "Frequency.csv")

data_prime = {"HeadHeadPrime": HH_prime_frequency, "TailTailPrime": TT_prime_frequency}
output_frequency_prime = pandas.concat(data_prime, axis=1)
output_frequency_prime = output_frequency_prime.sort_index()
output_frequency_prime.to_csv(pathlib.Path.cwd() / "Data_Files" / "GNum" / f"Len{length}" / "FrequencyPrime.csv")