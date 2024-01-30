import pandas
import pathlib
cwd = pathlib.Path.cwd()
import getKnownGenev2 as GetGene
import ucsc_restapi as API
import RQuery
import re


def main():
    '''
    Because I need to get a whole lot of species.

    Imports from the file 'Detailed_Genome.csv' with these headers:
    Scientific Name	Domain	Kingdom	Phylum	Class	Order	Family	Genius	Common Name	Latest Genome	Wiki

    And then it will go through them all and collect their known genes.

    Note that file does not include Homo Spaiens. Deal with it.

    '''

    file = cwd / "Detailed_Genome.csv"
    species: str
    species_dir = cwd.parent / "Data_Files" / "SpeciesChroms"
    species_dir.mkdir(parents=True, exist_ok=True)
    with open(file, "r") as file:
        data = pandas.read_csv(file, index_col="Scientific Name")
    
    specieses = tuple(data.index)  # I know it looks stupid but whatever
    for species in specieses:
        file_name = species.replace(" ", "_")
        file_name =  species_dir / f"{file_name}.txt"
        print(f"Species: {species}")
        genome = data.loc[species, "Latest Genome"]
        print_chromes(genome, file_name, species)


def print_chromes(genome: str, file_name: pathlib.Path, species: str):
    '''
    Prints to console the avalible chromosome
    '''
    baseURL, _, _, list_chromes, *_ = API.base_urls()
    chromURL = f"{baseURL}{list_chromes}?genome={genome}"

    chromes = RQuery.query(chromURL)
    chromes = API.convertRequest(chromes, dict_return=True)
    chromes = chromes["chromosomes"]
    chromes_List = list(chromes.keys())

    with open(file_name, "w") as file:
        file.write(f"URL Used:\n{chromURL}\n")
        file.write(f"Species: {species}\n")
        for chrom in chromes_List:
            if re.match(r"chr", chrom):
                if re.match(r"chrUn", chrom) or re.match(r"chrM", chrom) or re.search(r"random", chrom):
                    pass
                else:
                    file.write(f"\t{chrom}\n")
    # print(f"All avalible chromosomes for {genome}:")
    # for chrom in chromes_List:
    #     print(f"\t{chrom}")
    
    # print(f"URL Used:\n{chromURL}")




if __name__ in "__main__":
    main()