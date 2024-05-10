import pandas
import pathlib
cwd = pathlib.Path.cwd()
import random
import seqlogo
import numpy


def main():
    '''
    '''
    # Load up Data
    # Count each instance of each nucleotide
    # Put into data frame as:
    #   1   2   3   ... N
    # A a1  a2  a3  ... aN
    # G g1  g2  g3  ... gN
    # T t1  t2  t3  ... tN
    # C c1  c2  c3  ... cN
    # fake_5prime = [randomword(6) for _ in range(10)]
    # fake_3prime = [randomword(6) for _ in range(10)]

    # # fake_regions = [randomregion() for _ in range(10)]
    # fake_data = pandas.DataFrame(data = {"Anterior": fake_5prime, "Posterior": fake_3prime})
    # print(fake_data)
    # anterior, posterior = boundary_motiffs(fake_data, normalize=True)

    # print(anterior)
    # print(posterior)
    # # print(randomword(4))
    # # print(randomword(4))
    # # print(randomword(4))
    # # print(randomword(4))


    # # Look I can't duplicate their logo. Not sure why.
    # # For Generating the sequence count frames
    exon2exon_data: pandas.DataFrame = pandas.read_excel(cwd.parent / "Data_Files" / "Normal_Exon2Exon_v2.xlsx", sheet_name = "Sheet1")
    exon2exon_data = exon2exon_data.sample(n = 100)
    exon2exon_data = exon2exon_data.reset_index()
    # print(exon2exon_data)
    # print(exon2exon_data.head())
    # exit()

    exon2exon_data["Anterior"] = exon2exon_data["Anterior_10"].apply(lambda x: x[0: 3])
    exon2exon_data["Posterior"] = exon2exon_data["Posterior_10"].apply(lambda x: x[len(x) - 3: len(x)])

    anterior_counts, posterior_counts = boundary_motiffs(exon2exon_data)

    # anterior_counts.to_pickle(cwd / "Anterior_Counts_dd.pkl")
    # posterior_counts.to_pickle(cwd / "Posterior_Counts_dd.pkl")


    # posterior_counts = pandas.read_pickle(cwd / "Posterior_Counts.pkl")
    posterior_counts = normalize(posterior_counts)
    # posterior_counts.to_pickle(cwd / "Posterior_Norm_dd.pkl")
    print("Posterior")
    print(posterior_counts)
    print("\n")

    # anterior_counts = pandas.read_pickle(cwd / "Anterior_Counts.pkl")
    anterior_counts = normalize(anterior_counts)
    # anterior_counts.to_pickle(cwd / "Anterior_Norm_dd.pkl")
    print("Anterior")
    anterior_counts = anterior_counts.rename(columns={1: 4, 2: 5, 3: 6})
    print(anterior_counts)

    ppm = pandas.concat([posterior_counts, anterior_counts], axis = 1).T
    # ppm = ppm.rename(columns={"A": 0, "C": 1, "G": 2, "T": 3})
    print("\n")
    print(ppm)
    print(ppm.sum(axis = 1))
    ppm = seqlogo.Ppm(ppm)
    seqlogo.seqlogo(ppm, ic_scale = False, format = 'svg', size = 'medium')



def normalize(motifs: pandas.DataFrame):
    '''
    '''
    return motifs / motifs.sum()


def randomword(length):
    '''
    '''
    nucleotides = ["A", "G", "T", "C"]
    return ''.join(random.choice(nucleotides) for i in range(length))


def randomregion():
    '''
    '''
    region = random.randint(1, 2)
    return "exon" if region == 1 else "intron"


def boundary_motiffs(sequences: pandas.DataFrame, k: int = 3, 
                     normalize: bool = False,
                     five_prime: str = "Anterior", three_prime: str = "Posterior"):
    '''
    Mainly for a 3'|5' type boundary, where both are the same length
    '''
    columns = [x + 1 for x in range(k)]
    counts_anterior = pandas.DataFrame(data = 0, index = ["A", "G", "T", "C"], columns = columns)
    counts_posterior = pandas.DataFrame(data = 0, index = ["A", "G", "T", "C"], columns = columns)

    rows, _ = sequences.shape
    for row in range(rows):
        fiveprime_of_interest = sequences.loc[row, five_prime]
        threeprime_of_interest = sequences.loc[row, three_prime]

        for position in range(k):
            anterior_nuc = fiveprime_of_interest[position]
            posterior_nuc = threeprime_of_interest[position]

            counts_anterior.loc[anterior_nuc, position + 1] += 1
            counts_posterior.loc[posterior_nuc, position + 1] += 1

    
    if normalize:
        anterior_norm = counts_anterior.sum()
        posterior_nuc_norm = counts_posterior.sum()

        counts_anterior = counts_anterior / anterior_norm
        counts_posterior = counts_posterior / posterior_nuc_norm

    return counts_anterior, counts_posterior






if __name__ in "__main__":
    main()