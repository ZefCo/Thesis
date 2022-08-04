from requests import Response
import blat_api as ba
import ucsc_restapi as api
import pathlib
import pandas
import RQuery
import json

# Stopped at MGA_ANKLE1
# ENST00000389936_ENST00000394458
# which was having ENST Query errors, so skip it.
# Overall I have 8436 Blat/Enst/etc fusions
# out of 26809 - or 66% of total number of fusions
# and there are >40000

# Total there are 57613 entries in the UT database
# There are 40452 entries >=1000 nucleotides in sequence

# https://groups.google.com/a/soe.ucsc.edu/g/genome/c/U-w4b_ZS2j0?pli=1


def cheap_log_file_output(message = None):
    '''
    '''

    cheap_log_file = pathlib.Path.cwd() / "CheapLogFile.txt"
    cheap_log_file.touch(exist_ok=True)
    with open(cheap_log_file, 'a+') as log_file:
        log_file.write(message)



def main():
    '''
    '''
    output_headers = ['Hgene', 'Henst', 'Hstrand', 'Hchr', 'Tgene', 'Tenst', 'Tstrand', 'Tchr', 'tail_blockCount', 'tail_blockSizes', 'tail_tStarts', 'tail_txStart', 'tail_txEnd', 'tail_cdsStart', 'tail_cdsEnd', 'tail_exonCount', 'tail_exonStarts', 'tail_exonEnds', 'tail_name2', 'tail_cdsStartStat', 'tail_cdsEndStat', 'tail_exonFrames', 'head_blockCount', 'head_blockSizes', 'head_tStarts', 'head_txStart', 'head_txEnd', 'head_cdsStart', 'head_cdsEnd', 'head_exonCount', 'head_exonStarts', 'head_exonEnds', 'head_name2', 'head_cdsStartStat', 'head_cdsEndStat', 'head_exonFrames']
    
    min_length = 1000

    out_file = pathlib.Path.cwd().parent / "Data_Files" / f"UT_Blat_ENST_{min_length}.csv"
    if not out_file.is_file():
        pandas.DataFrame(columns=output_headers).to_csv(out_file, header = True, index = False)

    with open(pathlib.Path.cwd().parent / "Data_Files" / "UTData_cds.csv") as csvfile:
        utdata = pandas.read_csv(csvfile, header = 0)

    utdata = utdata[utdata["SeqLen"] >= min_length]

    rows, _ = utdata.shape
    utda_desiered = ["Hgene", "Henst", "Hstrand", "Hchr", "Tgene", "Tenst", "Tstrand", "Tchr"]
    enst_desiered = ['txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds','name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    blat_desiered = ["blockCount", "blockSizes", "tStarts"]

    # new_data = pandas.DataFrame()

    for row in range(rows):
        row_of_interest = utdata.iloc[row, :]
        sequence: str = row_of_interest["Seq"]

        row_of_interest: pandas.Series = row_of_interest[utda_desiered].squeeze()
        hgene, henst, hstrand, hchr, tgene, tenst, tstrand, tchr = row_of_interest["Hgene"], row_of_interest["Henst"], row_of_interest["Hstrand"], row_of_interest["Hchr"], row_of_interest["Tgene"], row_of_interest["Tenst"], row_of_interest["Tstrand"], row_of_interest["Tchr"]

        print(f"\n####\n{hgene}_{tgene}\n{henst}_{tenst}")

        try:
            blat: pandas.DataFrame = ba.blat_query(sequence)

        except json.decoder.JSONDecodeError as e:

            cheap_log_file_output(f"########\nBlat JSON Error: Will have to look in Blat API\n{henst}_{tenst}\n")
        except Exception as e:
            print("$$$$ Error at BLAT, adding to 'log'")

            cheap_log_file_output(f"########\nBlat Error\n{e}\n\t{type(e)}\n{henst}_{tenst}\n")

            continue

        if isinstance(blat, pandas.DataFrame):
            blat = blat[((blat["strand"] == hstrand) & (blat["tName"] == hchr)) | ((blat["strand"] == tstrand) & (blat["tName"] == tchr))]
            blat = blat[["strand", "tName", "blockCount", "blockSizes", "tStarts"]]

            # ENST identification
            bows, _ = blat.shape

            # # if we match to more then 2 things, just ditch it. This will still give me plenty of data to work with.
            # if bows >= 3:
            #     continue

            # if bows != 2:
            #     continue
            failed = False

            for bow in range(bows):
                blat_of_interest = blat.iloc[bow, :]

                enst_predict_chrom, enst_predict_strand = blat_of_interest["tName"], blat_of_interest["strand"]
                five_prime_end = -1 if enst_predict_strand in "-" else 0
                # This might need some explaination: UCSC Genome Browser reports everything Left to Right (LR), and infers the anit-sense/- strand from this.
                # So when reading the blat results for the - strand everything has to be read Right to Left. This means everything is backwards, which
                # provides a confusing picture. If you want the 5' end of a - strand you want to grab the LAST index and then ADD the block size
                # to get the end position, but this block is actually REVERSED (so you're grabbing the end, and going back to the begining). This means you 
                # get one of two functions to find the start end positions:

                # Honestly it seems a little weird, but if you think very carefully it makes sense. Just remember: the - strand information is BACKWARDS and they do this
                # so they dont have to save two strands of information. Everything is reported Left to Right, but the - strand is Right to Left.
                enst_predict_start, enst_predict_end = blat_of_interest["tStarts"][five_prime_end], blat_of_interest["tStarts"][five_prime_end] + blat_of_interest["blockSizes"][five_prime_end]

                
                print(f"Chrom: {enst_predict_chrom} # Start: {enst_predict_start} # End: {enst_predict_end} # Strand: {enst_predict_strand}")
                # print(blat_of_interest)

                enst_predict_url = api.ens_tracks(chrom = enst_predict_chrom, start = enst_predict_start, end = enst_predict_end)
                # print(enst_predict_url)

                enst_query: Response = RQuery.query(enst_predict_url)
                # print(f"ENST Response: {enst_query}\tType: {type(enst_query)}")
                # Find a way to deal with <Response [400]> cause it doens't show up in the 
                # if enst_query.
                
                try:
                    enst_query: pandas.DataFrame = api.convert2frame(enst_query)

                except AttributeError as e:
                    failed = True
                    print("$$$ Attribute Error at ENST Query")

                    cheap_log_file_output(f"########\nENST Attribute Error at ENST Query\n{enst_query}\n{enst_predict_url}\n")

                    continue

                except Exception as e:
                    failed = True
                    print("$$$$ Error at ENST Query, adding to 'log'")

                    cheap_log_file_output(f"########\nENST Error\n{e}\n\t{type(e)}\n{enst_predict_url}\n")

                    continue


                enst_query = enst_query[(enst_query["name"] == henst) | (enst_query["name"] == tenst)]
                # enst_name = enst
                # print(enst_query)
                # print(tuple(enst_query.index))
                # print(type(enst_query))

                # Yeah I probably could establish this before filtering the thing down to one, but this just seemd more direct: filter it to one thing and then go from there.
                try:
                    ht_name = enst_query.loc[tuple(enst_query.index)[0], "name"]
                except IndexError:
                    failed = True
                    print("\n$$$Index Error$$$")

                    cheap_log_file_output(f"########\nIndex Error\n{enst_query}\n{enst_predict_url}\n")

                    continue
                
                except Exception as e:
                    failed = True
                    print("\n$$$Error at 112$$$")

                    cheap_log_file_output(f"########\nNew at 112 Error\n{e}\n\t{type(e)}\n{enst_predict_url}\n")

                    continue


                local_blat: pandas.Series = blat_of_interest[blat_desiered].squeeze()
                local_enst: pandas.Series = enst_query[enst_desiered].squeeze()

                fusion_position = "head" if ht_name in henst else "tail"

                local_group: pandas.Series = pandas.concat([local_blat, local_enst], axis = 0)
                # print(local_group)

                new_index = blat_desiered + enst_desiered
                new_index = [f"{fusion_position}_{header}" for header in new_index]
                local_group: pandas.Series = local_group.set_axis(new_index)
                # print(local_group)

                row_of_interest = pandas.concat([row_of_interest, local_group])
                # print(row_of_interest)

            # This should adjust the order so that they are always in the same order
            if failed:
                continue
            else:
                try:
                    row_of_interest = row_of_interest[output_headers]
                except Exception as e:
                    print("$$$Error at resorting row of interst$$$")

                    cheap_log_file_output(f"########\nError at Resoriting Row of Interest\n{e}\n\t{type(e)}\n")

                    continue

                row_of_interest.to_frame().T.to_csv(out_file, header = None, index = None, mode = 'a')
                print(f"Finished {row}")

    #         try:
    #             new_data = pandas.concat([new_data, row_of_interest.to_frame()], axis = 1)
    #         except Exception as e:
    #             print("$$$$ Error at adding data to New Data, adding to 'log'")

    #             cheap_log_file = pathlib.Path.cwd() / "CheapLogFile.txt"
    #             cheap_log_file.touch(exist_ok=True)
    #             with open(cheap_log_file, 'a+') as log_file:
    #                 log_file.write(f"########\nENST Error\n{e}\n\t{type(e)}\n")

    #         _, new_rows = new_data.shape
    #         print(f"New Rows in New Data: {new_rows}")
    
    # # print(type(new_data))
    # # print(new_data.T)
    # new_data = new_data.T
    # new_data.to_csv(pathlib.Path.cwd().parent / "Data_Files" / f"UT_Blat_ENST_{min_length}.csv")
    #             # print(local_blat)
    #             # print(local_enst)
    #             # local_gene: pandas.Series = pandas.concat([local_blat, local_enst])
                
    #             # index_desiered = [f"{fusion_position}_{index}" for index in enst_desiered]
    #             # # print(enst_query.columns)
    #             # # enst_query = enst_query[enst_desiered]
    #             # enst_query = enst_query.set_axis(enst_desiered)
    #             # print(enst_query)

    #             # # if isinstance(query, requests.models.Response):
    #             # #     query: dict = json.loads(query.text)
    #             # #     britta = query['blat']

    #             #     # query = 



    #         # for col in blat:
    #             # value = blat.loc["0", col]
    #             # print(f"Column: {col}\t Type: {type(value)}")

            

    #     # if row == 0:
    #         # break

    # # print(blat)
    # # print(blat.loc["0", 'blockCount'])
    # # print(type(blat.loc['0', 'blockCount']))




if __name__ in '__main__':
    main()