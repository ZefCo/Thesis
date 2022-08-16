import pandas
import Bio.SearchIO
from io import StringIO
import Bio.Blast.Record
from Bio.SearchIO import BlatIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pathlib
import xml
import re
import socket
import urllib.error

with open(pathlib.Path.cwd() / "Data_Files" / "UTData_cds.csv") as utfile:
    ut_data = pandas.read_csv(utfile, header = 0)

with open(pathlib.Path.cwd() / "Data_Files" / "Sequence_Files" / "HG19" / "Known_Genes_hg19.csv") as hg19file:
    hg19_data = pandas.read_csv(hg19file, header = 0)

# Create an empty dataframe that will hold the new data. Not sure if I'm going to keep doing this: I want to append the 
# data to a csv as the wwwblastn is pretty finiky.
printable_frame = pandas.DataFrame(columns = hg19_data.columns)

blast_output_file = pathlib.Path.cwd() / "Data_Files" / "Blast_Genes_HG19.csv"
# Write the file if it does not exist. Include the Headers
if not blast_output_file.is_file():
    printable_frame.to_csv(blast_output_file, header = True, index = False)
    blasted_list = []
else:
    with open(blast_output_file) as initialize_file:
        blasted_list = pandas.read_csv(initialize_file, header = 0)['name'].to_list()


new_ut_file = pathlib.Path.cwd() / "Data_Files" / "New_UTData_cds.csv"
# Write the file if it does not exist. Include the Headers
if not new_ut_file.is_file():
    columns = list(ut_data.columns) + ["HNames", "TNames"]
    pandas.DataFrame(columns = columns).to_csv(new_ut_file, header = True, index = False)

# Take shape of the UT data, going to iterate over this later. Create a start position - so that if the script crashes I can tell it
# where to resume - and use that in the for loop
start_position = 668
ut_rows, _ = ut_data.shape
# Size for the alignment window
window_size = 100
minimum_size = 1000
# Database to be used for BLAST. This database includes human and non human sequences, so we have to do some
# extra work to get what we need. I tired a few other databases, and they did not yeild any results, so this
# is the one we're going to work with
database = "refseq_rna"


def query_result(sequence, genenames, database = "refseq_rna", program = "blastn", format = "blast-xml"):
    '''
    '''
    try:
        blast_result: StringIO = NCBIWWW.qblast(program = program, database = database, sequence = sequence)
    except urllib.error.URLError:
        blast_result = None
        print("URL error avoided")
    except socket.gaierror:
        blast_result = None
        print("Socket Error avoided")
    except Exception as e:
        blast_result = None
        print(f"Add Type Error:\n{type(e)}\nTo exceptions list")
    
    if isinstance(blast_result, StringIO):
        try:
            blast_read: Bio.SearchIO._model.query.QueryResult = Bio.SearchIO.read(blast_result, format = format)
        except Exception as e:
            print(type(e))
        # print(f"Hit Keys: {blast_read.hit_keys}")
    
        return_names = []

        for key in blast_read.hit_keys:
                for name in genenames:
                    if re.search(name, key):
                        # print(key)
                        # print(blast_read[key])
                        # Then use this to figure out which genes I actually need to use? Filter the known gene
                        # and UT database to a smaller set?
                        return_names.append(name)

        if len(return_names) == 0:
            return_names = None

    else:
        return_names = None

    return return_names


for row in range(start_position, ut_rows):

    
    row_of_interest = ut_data.iloc[row, :].copy()

    filtered_data = pandas.DataFrame()
    full_seq = ut_data.loc[row, "Seq"]
    head_name, tail_name = ut_data.loc[row, "Hgene"], ut_data.loc[row, "Tgene"]

    # Grab the relevant data from the known gene layout
    hg19_headdata = hg19_data[hg19_data['name2'] == head_name]
    hg19_taildata = hg19_data[hg19_data['name2'] == tail_name]
    hg19_subdata = pandas.concat([hg19_headdata, hg19_taildata], axis = 0)

    # I like sets, they're fast
    # Creating 3 sets: one for all the unique names, one for all the unique head names, and one for all the unique tail names. The reason for this is that sometimes the
    # BLAST identifies just one set of genes, and that might not include the head or the tail. This way I can cross reference the heads and the tails to make sure that
    # both showed up. Look at the ZNF148_SLC41A3 fusion: the cds sequence will only come back with SLC41A3, which is then put into both the ZNF148 & SLC41A3 data, screwing up
    # any results.
    hg19_unique_names = set(hg19_subdata['name'])
    hg19_unique_heads = set(hg19_headdata['name'])
    hg19_unique_tails = set(hg19_taildata['name'])


    # since we're using the whole sequence: make sure it's bigger then a minimum size
    if len(full_seq) >= minimum_size:
        print(f"### {head_name}_{tail_name} ###")

        head_gene_names: list or None = query_result(full_seq[0:window_size], hg19_unique_heads)
        tail_gene_names: list or None = query_result(full_seq[len(full_seq) - window_size: len(full_seq)], hg19_unique_tails)

        # Make sure there is a head and tail unique name present
        # head_present, tail_present = set(head_gene_names) & hg19_unique_heads, set(tail_gene_names) & hg19_unique_tails

        if isinstance(head_gene_names, list) and isinstance(tail_gene_names, list):
            # if (len(head_gene_names) > 0) and (len(tail_gene_names) > 0):
            print(f"Head Gene List:\n{head_gene_names}")
            print(f"Tail Gene List:\n{tail_gene_names}")
            row_of_interest["HNames"] = head_gene_names
            row_of_interest["TNames"] = tail_gene_names

            # Filter the subset of hg19 data, composed of the fusion genes from the UT Database, to include only the genes that show up in the BLAST.
            filtered_data: pandas.DataFrame = hg19_subdata[hg19_data['name'].isin(head_gene_names) | hg19_data['name'].isin(tail_gene_names)]

            # If those names do not appear in the BLASTed Gene File, add them.
            # I should create a new column in the UT Database with the gene name of everything I've found in BLAST, that way I can easily re run the Exon2Exon script
            # It would be the 4th version but that's OK. It would be the best version so far.

            # So if I did this right:
            # Cast the LIST of names from the BLAST_Genes to a SET -> call that set B
            # Create a SET from the LIST of head and tail genes -> call that set G
            # Find G - B = G'
            blasted_set = set(blasted_list)
            gene_set = set(head_gene_names + tail_gene_names)
            gene_diff = gene_set.difference(blasted_set)
            print(f"Blasted Set:\n{blasted_set}\n\nGene Set:\n{gene_set}\n\nDifference:\n{gene_diff}")

            # Find ever g' in G', add it to the LIST of B for next iteration
            for gene_name in gene_diff:
                blasted_list.append(gene_name)

            # Remove anything that is not in G' from the dataframe
            filtered_data = filtered_data[filtered_data['name'].isin(gene_diff)]

            # If that Dataframe is not empty: write to our output file
            if filtered_data.shape[0] > 0:
                # with open(blast_output_file, "a+") as blast_file:
                filtered_data.to_csv(blast_output_file, header = None, mode = "a", index = False)

                # with open(new_ut_file, "a+") as newer_file:
                row_of_interest.to_frame().T.to_csv(new_ut_file, header = None, mode = "a", index = False)
        else:
            print(f"Length of Head Gene List: {len(head_gene_names) if isinstance(head_gene_names, list) else 'Null Set'}\nLength of Tail Gene List: {len(tail_gene_names) if isinstance(tail_gene_names, list) else 'Null Set'}\nWill be skipping this one due to at least one legnth being 0")

    # Let me know what row we made it to
        print(f"### Finished row {row} of UT Database ####")
        # exit()

        # This code was written to find the actual BP. Still working on that, but right now I want to 
        # Just figure out which specific genes I'm using and rerun it through ExonExonv3.py
        # for window_start in range(int(len(full_seq)/window_size) + 1):
        #     if ((window_start * window_size) + window_size - 1) > len(full_seq):
        #         # print(window_start * window_size, len(full_seq))
        #         sub_seq = full_seq[window_start * window_size: len(full_seq)]
        #     else:
        #         # print((window_start * window_size), (window_start * window_size) + window_size - 1)
        #         sub_seq = full_seq[window_start * window_size: (window_start * window_size) + window_size - 1]
        #     print(sub_seq)

        #     # result = BlatIO(sub_seq)

        #     result: StringIO = NCBIWWW.qblast("blastn", database, sub_seq)
        #     # read_result: Bio.Blast.Record.Blast = NCBIXML.read(result)

        #     # print(f"Length of Description and Alignments with database = {database}")
        #     # print(len(read_result.descriptions))
        #     # print(len(read_result.alignments))

        #     # # So basically all of these are lists, and they look to the be the same size?

        #     # for i, alignment in enumerate(read_result.alignments):
        #     # #     print("Alignment from BLAST")
        #     #     print(read_result.descriptions[i])
        #     #     print(alignment)
        #     # #     print("HSPS and Type (should just be memory location)")
        #     #     print(alignment.hsps[0])
        #     # #     print(type(alignment.hsps))
        #     #     # print(len(alignment.hsps))

        #     # So now we have this stuff in a QueryResult...
        #     shoving: Bio.SearchIO._model.query.QueryResult = Bio.SearchIO.read(result, 'blast-xml')

        #     # print(shoving.hit_keys)
        #     for key in shoving.hit_keys:
        #         for name in hg19_names:
        #             if re.search(name, key):
        #                 # print(key)
        #                 print(shoving[key])
        #                 # Then use this to figure out which genes I actually need to use? Filter the known gene
        #                 # and UT database to a smaller set?
        #         # if key in 

        #     break
    # break

            # print(type(record))
            # print(record.alignments)

