import itertools
import pandas
import numpy


def main():
    nuc_perm: dict = nucleotide_permutations()

    print(nuc_perm)


def nucleotide_permutations(length: int = 5) -> dict:
    nuc_perm = dict()
    for len in range(1, length + 1):
        nuc_perm[len] = [thing for thing in itertools.combinations_with_replacement("ACGT", len)]

    return nuc_perm



def similarity_key(length: int = 9, rows: int = 1) -> pandas.DataFrame:
    print_frame = pandas.DataFrame()

    weights = numpy.array([(1 / (2**(i))) for i in range(1, length + 1)])

    insane = [numpy.reshape(numpy.array(i), (rows, length)) for i in itertools.product([0, 1], repeat = rows*length)]
    insane_part2 = [numpy.sum(nuts * weights) for nuts in insane]

    insane = [thing[0] for thing in insane]

    print_frame = pandas.DataFrame(data = {"Similarity": insane, "G Score": insane_part2})
    # print_frame.to_csv(pathlib.Path.cwd() / f"Similarity Vector Key Length - {length}.csv", index = False)

    return print_frame






if __name__ in "__main__":
    main()