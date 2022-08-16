from json import JSONDecodeError
from sqlite3 import DatabaseError
import blat_api as bpi
import ucsc_restapi as upi
# import RQuery
# import requests
import pandas
from GeneClass import Gene
from typing import Tuple


class FusionGene():
    def __init__(self, hgene: str, tgene: str, seq: str, henst: str, tenst: str, hstrand: str, tstrand: str, hchrm: str, tchrm: str, gdb: str = "hg19") -> None:
        # everything here, unless otherwise noted, is mostly likely a string
        self.hgene, self.tgene = hgene, tgene
        self.henst, self.tenst = henst, tenst
        self.hchrm, self.tchrm = hchrm, tchrm
        self.hstrand, self.tstrand = hstrand, tstrand
        self.seq = seq
        self.gdb = gdb

        self.gnames = f"{self.hgene}_{self.tgene}"
        self.enames = f"{self.henst}_{self.tenst}"

        self.error_attempts: int = 5


    def blat(self):
        '''
        Blats the Sequence. Will only keep the blockCount, blockSizes, tStarts

        This will then trigger the ENST Identification. Why? Because BLAT only tells the genomic location, not what you're looking at.
        '''
        filter_blat = lambda blat_frame, strand, chrom: blat_frame[(blat_frame["strand"] == strand) & (blat_frame["tName"] == chrom)]

        # blat_response = bpi.blat_query(self.seq)
        blat_response, blat_url = self.query_attempt(bpi.blat_query, self.seq)
        self.blat_url: str = blat_url
        # print(blat_response)
        if isinstance(blat_response, pandas.DataFrame):
            head_response: pandas.DataFrame = filter_blat(blat_response, self.hstrand, self.hchrm)
            tail_response: pandas.DataFrame = filter_blat(blat_response, self.tstrand, self.tchrm)
        
        else:
            print("!!!!\tBlat Failed\t!!!!")
            print("Printing to logs?")
            # This might not be the best way, but if it can't BLAT then I just want to escape out of this. It's not going to do anything
            # So yes, this is a bare return
            return

        if (head_response.shape[0] > 0) and (tail_response.shape[0] > 0):
            print(f"~~~~\tClean Blat\t~~~~")

            # print(f"Blat => Head: {head_response.shape[0]} - Tail: {tail_response.shape[0]}")

            head_enst, head_enst_url = self.enst_identify(enst = self.henst, blat_data = head_response)
            tail_enst, tail_enst_url = self.enst_identify(enst = self.tenst, blat_data = tail_response)

            self.head_enst_url, self.tail_enst_url = head_enst_url, tail_enst_url

            if (isinstance(head_enst, pandas.DataFrame) and (head_enst.shape[0] > 0)) and (isinstance(tail_enst, pandas.DataFrame) and (tail_enst.shape[0] > 0)):
                print(f"~~~~\tClean ENST\t~~~~")

                self.head_gene = Gene(self.hgene, chrm = self.hchrm, strand = self.hstrand, ename=self.henst,
                                    gname = head_enst["name2"],
                                    txStart = head_enst["txStart"],
                                    txEnd = head_enst["txEnd"],
                                    cdsStart = head_enst["cdsStart"],
                                    cdsEnd = head_enst["cdsEnd"],
                                    exonCount = head_enst["exonCount"],
                                    exonStarts = head_enst["exonStarts"],
                                    exonEnds = head_enst["exonEnds"],
                                    exonFrames = head_enst["exonFrames"])

                self.tail_gene = Gene(self.tgene, chrm = self.tchrm, strand = self.tstrand, ename=self.tenst,
                                    gname = tail_enst["name2"],
                                    txStart = tail_enst["txStart"],
                                    txEnd = tail_enst["txEnd"],
                                    cdsStart = tail_enst["cdsStart"],
                                    cdsEnd = tail_enst["cdsEnd"],
                                    exonCount = tail_enst["exonCount"],
                                    exonStarts = tail_enst["exonStarts"],
                                    exonEnds = tail_enst["exonEnds"],
                                    exonFrames = tail_enst["exonFrames"])

            else:
                head_rows = head_enst.shape[0] if isinstance(head_enst, pandas.DataFrame) else 0
                tail_rows = tail_enst.shape[0] if isinstance(tail_enst, pandas.DataFrame) else 0
                print(f"One did not ENST Identify Correctly\nHead Rows: {head_rows}\tTail Rows: {tail_rows}")

        # else:
        #     print(f"One does not identify correctly\nHead Rows: {head_response.shape[0]}\tTail Rows: {tail_response.shape[0]}")

        #     print("Printing to Logs?")
        # # exit()




    def classify(self):
        '''
        Classify if it's Trans Inter/Intra Genic, Cis SAG
        '''




    def enst_identify(self, blat_data: pandas.DataFrame, enst: str) -> Tuple[pandas.Series, str]:
        '''
        Pulls up the ENST track and matches it to the expected ENST name of the fusion.
        '''

        rows, _ = blat_data.shape

        for row in range(rows):
            row_of_interest = blat_data.iloc[row, :].copy()

            gstart, gend, chrm = row_of_interest["tStarts"][0], row_of_interest["tStarts"][0] + row_of_interest["blockSizes"][0], row_of_interest["tName"]

            # enst_url = upi.ens_tracks(chrom=chrm, start=gstart, end=gend)

            enst_response, enst_url = self.query_attempt(upi.ens_tracks, chrom=chrm, start=gstart, end=gend)
            # print(f"!!!!\t\tFresh out of Query Attempt:\t\t!!!!\n{enst_response}")

            # try:
            if isinstance(enst_response, pandas.DataFrame):
                if enst_response.shape[0] > 0:
                    enst_response = enst_response[enst_response["name"] == enst]

                    if enst_response.shape[0] == 1:
                        enst_response.squeeze()
                        break
                    elif enst_response.shape[0] > 1:
                        print(f"!!!!\t{enst} isn non unique\t!!!!")
                        break
                else:
                    enst_response = None
            # except Exception as e:
            #     print(f"Error: {e}\tType{type(e)}")
            #     print(enst_url)
            #     # print(enst_response.shape[0])

        return enst_response, enst_url


            # exit()



    def query_attempt(self, func, *args, **kwargs) -> Tuple[pandas.Series or pandas.DataFrame or Exception, str]:
        '''
        Returns the usable response from the UCSC Genome browser (typically a Pandas Series or DataFrame)
        along with the url used.
        '''
        # print(type(func))
        # print(*args)
        # print(kwargs)
        # exit()
        for attempt in range(self.error_attempts):
            try:
                ucsc_data, ucsc_url = func(*args, **kwargs)
                if attempt > 0:
                    print(f"Successfully finished after {attempt} failed attempts")
                
                break
                # print(ucsc_data)
                # exit()

            except JSONDecodeError as e:
                print("!!!!\tJSON Decoder Error\t!!!!")
                print("    \tretrying")
                
                ucsc_data = e
                ucsc_url = None

            except Exception as e:
                print("!!!! New Error trying to get UCSC Data !!!!")
                print(f"Error: {e}\tType: {type(e)}")

                ucsc_data = e
                
                break
                # exit()

        else:
            print("!!!!\tFailed to get UCSC Data, sending to log\t!!!!")

        return ucsc_data, ucsc_url


    def write_to_csv(self, outfile: str):
        '''
        Writes to the error or output file
        '''




def main():
    '''
    '''


if __name__ in '__main__':
    main()