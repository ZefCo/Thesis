import pandas
import pathlib
cwd = pathlib.Path.cwd()
import getKnownGenev2 as GetGene
import ucsc_restapi as API
import RQuery
import re

def main():
    '''
    '''
    # target_KPCOFGS(cwd / "Detailed_Genome.csv", "Primates", "Order",  cwd.parent / "Data_Files" / "Primates")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Callithrix_jacchus.txt",      cwd.parent / "Data_Files" / "Primates" / "Callithrix_jacchus")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Chlorocebus_sabaeus.txt",     cwd.parent / "Data_Files" / "Primates" / "Chlorocebus_sabaeus")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Gorilla_gorilla_gorilla.txt", cwd.parent / "Data_Files" / "Primates" / "Gorilla_gorilla_gorilla")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Macaca_fascicularis.txt",     cwd.parent / "Data_Files" / "Primates" / "Macaca_fascicularis")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Macaca_mulatta.txt",          cwd.parent / "Data_Files" / "Primates" / "Macaca_mulatta")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Nasalis_larvatus.txt",        cwd.parent / "Data_Files" / "Primates" / "Nasalis_larvatus")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Nomascus_leucogenys.txt",     cwd.parent / "Data_Files" / "Primates" / "Nomascus_leucogenys")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Pan_paniscus.txt",            cwd.parent / "Data_Files" / "Primates" / "Pan_paniscus")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Pan_troglodytes.txt",         cwd.parent / "Data_Files" / "Primates" / "Pan_troglodytes")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Papio_anubis.txt",            cwd.parent / "Data_Files" / "Primates" / "Papio_anubis")
    read_Species(cwd.parent / "Data_Files" / "Primates" / "Genome" / "Pongo_pygmaeus_abelii.txt",   cwd.parent / "Data_Files" / "Primates" / "Pongo_pygmaeus_abelii")




def read_Species(file_path: pathlib.Path, output_folder: pathlib.Path, *args, **kwargs):
    '''
    This reads the files generated from teh print chromosomes function and uses that info to starts getting genetic data on that species.
    '''
    with open(file_path, "r") as file:
        text = file.readlines()

    chroms = []

    for line in text:
        if re.match("Species", line):
            species = re.sub(r"Species: ", "", line)
            species = re.sub(r"\n", "", species)
        elif re.match("Genome", line):
            genome = re.sub(r"Genome: ", "", line)
            genome = re.sub(r"\n", "", genome)
        elif re.match(r"\t", line):
            local_chrom = re.sub(r"\t", "", line)
            local_chrom = re.sub(r"\n", "", local_chrom)
            chroms.append(local_chrom)
    
    try:
        get_Species_Genes(genome, chroms, output_folder, *args, **kwargs)
    except Exception as e:
        print(type(e))
        print(e)
        print(f"Skipping this species")


def get_Species_Genes(genome: str, 
                      chroms: list,
                      output_folder: pathlib.Path,
                      track: str = "ncbiRefSeqCurated",
                      data_headers: list = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'],
                      *args, **kwargs):
    '''
    Because I'm trying to get things to work cleaner.

    This grabs the known genes of the various species. It needs a list of chromosomes to target, because many species have different chromosomes. Not sure what will happen if the species has no chromosomes...
    Note about the data headers: those are the hg19 headers. They seem to be useful for a variety of species.
    '''
    output_folder.mkdir(parents=True, exist_ok=True)
    GetGene.get_data(genome, chroms, data_headers, output_folder, track)



def target_KPCOFGS(data: pathlib.Path, 
                   target: str,
                   taxonomy: str, 
                   output_path: pathlib.Path, 
                   *args, **kwargs):
    '''
    Kevin, Please Come Over For Gay Sex.

    This targets a specific Taxonomy and pulls all things related to that. Note in the detailed species file Domain and Realm are interchangable. The taxonomy is one of the columns.
    Can only target one taxonomy at a time, sorry.

    Why does it have to be Gay sex? Why not Great? Because it's Kevin. Why not Karen? Well now everyone is just confused.
    '''
    output_path.mkdir(parents=True, exist_ok=True)
    with open(data, "r") as file:
        data = pandas.read_csv(file, index_col="Scientific Name")

    data: pandas.DataFrame = data[data[taxonomy] == target]

    species: str
    specieses = tuple(data.index)  # I know it looks stupid but whatever
    for species in specieses:
        file_name = species.replace(" ", "_")
        file_name =  output_path / f"{file_name}.txt"
        print(f"Species: {species}")
        genome = data.loc[species, "Latest Genome"]
        print_chromes(genome, file_name, species)



def all_species():
    '''
    Because I need to get a whole lot of species.

    Imports from the file 'Detailed_Genome.csv' with these headers:
    Scientific Name	Domain	Kingdom	Phylum	Class	Order	Family	Genius	Common Name	Latest Genome	Wiki

    And then it will go through them all and collect their known genes.

    Note that file does not include Homo Spaiens. Deal with it. Or it might. Either way, deal with.

    This does everything
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
        file.write(f"Genome: {genome}\n")
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