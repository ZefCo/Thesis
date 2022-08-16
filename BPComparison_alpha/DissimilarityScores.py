import pandas
from CommonMethods import compliment_strand
import CommonMethods as CM
import pathlib
import ucsc_restapi as api
import RQuery



def main():
    '''
    '''
    score_length = 25
    headers_of_sets = ["HblockSizes", "HqStarts", "HtStarts",  "TblockSizes", "TqStarts", "TtStarts", "HexonStarts", "HexonEnds", "HexonFrames", "TexonStarts", "TexonEnds", "TexonFrames"]
    ut_be = CM.import_csv(pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "UTBE_1000.csv", header=0, col_sets=headers_of_sets)

    frame_headers = list(ut_be.columns)
    frame_headers = ["Classification"] + frame_headers
    ut_be["Classification"] = "Unknown"
    ut_be = ut_be[frame_headers]

    # Head 3' - : index = 0 of HtStarts, HblockSize
    # Head 3' + : index = HblockCount of HtStarts, HblockSize

    # Tail 5' - : index = TblockCount of TtStarts, TblockSize
    # Tail 5' + : index = 0 of TtStarts, TblockSize

    rows, _ = ut_be.shape
    # print(rows, type(rows))
    
    ut_bes = pandas.DataFrame()
    # print(ut_be)

    for row in range(0, 15):
        print(f"Doing for loop iteration {row}")
        row_of_interest: pandas.Series = ut_be.iloc[row, :].copy().squeeze()
        row_of_interest = classify_splice(row_of_interest)

        # type = string
        head_dir, tail_dir = row_of_interest["Hstrand"], row_of_interest["Tstrand"]
        # type = string
        head_chr, tail_chr = row_of_interest["Hchr"], row_of_interest["Tchr"]
        # type = int
        head_blockCount, tail_blockCount = row_of_interest["HblockCount"], row_of_interest["TblockCount"]
        # type = tuple(int)
        head_tStarts, tail_tStarts = row_of_interest["HtStarts"], row_of_interest["TtStarts"]
        # type = tuple(int)
        head_blockSize, tail_blockSize = row_of_interest["HblockSizes"], row_of_interest["TblockSizes"]

        # I could have written this as a seperate method, but I didn't, and it's because of how the head and tail get
        # turned around depending on weather they are in the sense (+) or antisense (-) direction. So instead of making
        # it simple yet abstract, I made it clunky yet direct. I think this is the right choice

        # I chose to do this as blockCount - 1 rather then slicing with [-1] because I felt it was a little safer. I could
        # grab the last thing in the block size and the Starts when needed, but I wanted to make explicet where the last
        # index comes from

        if head_dir in '+':
            head_start = head_tStarts[head_blockCount - 1]

            head_end = head_start + head_blockSize[head_blockCount - 1]

            head_index = ((head_blockSize[head_blockCount - 1] - score_length) if head_blockSize[head_blockCount - 1] > score_length else 0, head_blockSize[head_blockCount - 1])

        elif head_dir in '-':
            head_start = head_tStarts[0]

            head_end = head_start + head_blockSize[0]

            head_index = (0, score_length if score_length < head_blockSize[head_blockCount - 1] else head_blockSize[head_blockCount - 1])

        else:
            print("If somehow you get this message, you really messed up further upstream in the data aquisition: it can't match the direction of the gene to + or -")

        if tail_dir in '+':
            tail_start = tail_tStarts[0]

            tail_end = tail_start + tail_blockSize[0]

            tail_index = (0, score_length if score_length < tail_blockSize[tail_blockCount - 1] else tail_blockSize[tail_blockCount - 1])


        elif tail_dir in '-':
            tail_start = tail_tStarts[tail_blockCount - 1]

            tail_end = tail_start + tail_blockSize[tail_blockCount - 1]

            tail_index = ((tail_blockSize[tail_blockCount - 1] - score_length) if tail_blockSize[tail_blockCount - 1] > score_length else 0, tail_blockSize[tail_blockCount - 1])


        head_3prime_url = api.sequence(chrom = head_chr, start = head_start, end = head_end)
        # print(head_3prime_url)

        head_dna = api.convertRequest(RQuery.query(head_3prime_url))
        if head_dir in '-':
            # print("Finding Compliment")
            head_dna = compliment_strand(head_dna)
        
        head_dna = head_dna[head_index[0]: head_index[1]]
        
        # head_3prime_seq = 
        tail_5prime_url = api.sequence(chrom = tail_chr, start = tail_start, end = tail_end)
        # print(tail_5prime_url)

        tail_dna = api.convertRequest(RQuery.query(tail_5prime_url))
        if tail_dir in '-':
            # print("Finding Compliment")
            tail_dna = compliment_strand(tail_dna)

        tail_dna = tail_dna[tail_index[0]: tail_index[1]]

        print(f"########\nHead DNA:\t{head_dna}\nChr: {head_chr}\tDir: {head_dir}\tHead Index: {head_index}\nURL: {head_3prime_url}\nTail DNA:\t{tail_dna}\nChr: {tail_chr}\tDir: {tail_dir}\tTail Index: {tail_index}\nURL: {tail_5prime_url}")

        ut_bes = pandas.concat([ut_bes, row_of_interest.to_frame().T], axis = 0)

    # print(ut_bes)


def check_exons():
    '''
    This will check each exon for slippage
    '''




def classify_splice(fusion_data: pandas.Series, cis = "C-SAG", inter = "T-E", intra = "T-A", adjecent_def = 50_000, unknown = "Unknown") -> str:
    '''
    Classify the type of splice as:
    C-SAG ~ Cis-SAG ~ same chr, same dir, adjecent genes -> last exon/UTR is < 50 kb to start of exon/UTR of next gene
    T-A ~ Trans-IntrAgenic ~ same chr, same/diff dir, non-adjecent genes
    T-E ~ trans-IntErgenic ~ diff chr
    '''

    # type_check = ["Hchr", "Tchr", "Classification", "Hstrand", "Tstrand", "HcdsEnd", "TcdsStart"]
    # for tc in type_check:
    #     print(f"Index: {tc}\tValue: {fusion_data[tc]}\tType: {type(fusion_data[tc])}")

    # Type wise: everything but the "TcdsStart" and "HcdsEnd" are string. The two exceptions are int64

    if fusion_data["Hchr"] != fusion_data["Tchr"]:
        fusion_data["Classification"] = inter
    
    elif fusion_data["Hchr"] == fusion_data["Tchr"]:
        
        if fusion_data["Hstrand"] != fusion_data["Tstrand"]:
            fusion_data["Classification"] = intra
        
        elif fusion_data["Hstrand"] == fusion_data["Tstrand"]:
            
            if (fusion_data["Hstrand"] in "+") and (fusion_data["Tstrand"] in "+"):
            
                if fusion_data["HcdsEnd"] < fusion_data["TcdsStart"]:
            
                    if fusion_data["TcdsStart"] - fusion_data["HcdsEnd"] <= adjecent_def:
                        fusion_data["Classification"] = cis
                    
                    elif fusion_data["TcdsStart"] - fusion_data["HcdsEnd"] > adjecent_def:
                        fusion_data["Classification"] = intra

                    else:
                        fusion_data["Classification"] = unknown
                
                elif fusion_data["HcdsEnd"] > fusion_data["TcdsStart"]:
                    fusion_data["Classification"] = intra

                else:
                    fusion_data["Classification"] = unknown

            elif (fusion_data["Hstrand"] in "-") and (fusion_data["Tstrand"] in "-"):
                
                if fusion_data["HcdsEnd"] > fusion_data["TcdsStart"]:
                    
                    if fusion_data["HcdsEnd"] - fusion_data["TcdsStart"] <= adjecent_def:
                        fusion_data["Classification"] = cis

                    elif fusion_data["HcdsEnd"] - fusion_data["TcdsStart"] > adjecent_def:
                        fusion_data["Classification"] = intra
                    
                    else:
                        fusion_data["Classification"] = unknown
                
                elif fusion_data["HcdsEnd"] < fusion_data["TcdsStart"]:
                    fusion_data["Classification"] = intra

                else:
                    fusion_data["Classification"] = unknown
            
            else:
                fusion_data["Classification"] = "Unknown"
    
    else:
        fusion_data["Classification"] = "Unknown"

    return fusion_data



def sanity_type_check(name, thing):
    print(f"~~~~\nVariable: {name}\tValue: {thing}\tType: {type(thing)}\n")



# chr1:94,054,531-94,054,976           
if __name__ in '__main__':
    main()