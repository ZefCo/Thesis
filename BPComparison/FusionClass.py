from ast import Num
from json import JSONDecodeError
import blat_api as bpi
import ucsc_restapi as upi
# import RQuery
# import requests
import pandas
from GeneClass import Gene
from BlatClass import Blat
from ScoreingClasses import DissSimilarityScore
from typing import Tuple
import pathlib
import numpy as np


class FusionGene():
    def __init__(self, hgene: str, tgene: str, seq: str, henst: str, tenst: str, hstrand: str, tstrand: str, hchrm: str, tchrm: str, gdb: str = "hg19") -> None:
        # everything here, unless otherwise noted, is mostly likely a string
        self.gdb = gdb

        self.hgene, self.henst, self.hchrm, self.hstrand = hgene, henst, hchrm, hstrand
        self.tgene, self.tenst, self.tchrm, self.tstrand = tgene, tenst, tchrm, tstrand
        
        self.seq = seq

        self._clean_blat: bool = False
        self._clean_enst: bool = False

        self.classification: str = None
        # The reason for why so many types:
        # str means the distance between the two is meaningless (they are on different chromosomes)
        # complex means they are on different strands on the chromosome
        # int means they are on the same strand, same chromosome
        self.shortDistance: int or complex or str = None
        self.head2tailDistance: int or complex or str = None

        self.head_gene: Gene = None
        self.tail_gene: Gene = None

        self.head_blat: Blat = None
        self.tail_gene: Blat = None



    def blat(self):
        '''
        Blats the Sequence. Will only keep the blockCount, blockSizes, tStarts

        This will then trigger the ENST Identification. Why? Because BLAT only tells the genomic location, not what you're looking at.
        '''
        filter_blat = lambda blat_frame, strand, chrom: blat_frame[(blat_frame["strand"] == strand) & (blat_frame["tName"] == chrom)].copy()

        # blat_response = bpi.blat_query(self.seq)
        blat_response, blat_url = self.query_attempt(bpi.blat_query, self.seq)
        self.blat_url: str = blat_url
        # print(blat_response)
        if isinstance(blat_response, pandas.DataFrame):
            blat_response["name"] = "Unknown"
            head_response: pandas.DataFrame = filter_blat(blat_response, self.hstrand, self.hchrm)
            head_response.reset_index(drop = True, inplace = True)
            tail_response: pandas.DataFrame = filter_blat(blat_response, self.tstrand, self.tchrm)
            tail_response.reset_index(drop = True, inplace = True)
        
        else:
            print("!!!!\tBlat Failed\t!!!!")
            self._clean_blat = False


            # This might not be the best way, but if it can't BLAT then I just want to escape out of this. It's not going to do anything
            # So yes, this is a bare return
            return

        if (head_response.shape[0] > 0) and (tail_response.shape[0] > 0):
            print(f"~~~~\tClean Blat\t~~~~")
            self._clean_blat = True

            # print(f"Blat => Head: {head_response.shape[0]} - Tail: {tail_response.shape[0]}")

            head_enst, henstURL = self.enst_identify(enst = self.henst, blat_data = head_response)
            tail_enst, tenstURL = self.enst_identify(enst = self.tenst, blat_data = tail_response)

            self.henstURL, self.tenstURL = henstURL, tenstURL

            # if (isinstance(head_enst, pandas.Series) and (head_enst.shape[0] > 0)) and (isinstance(tail_enst, pandas.Series) and (tail_enst.shape[0] > 0)):
            if (isinstance(head_enst, pandas.Series)) and (isinstance(tail_enst, pandas.Series)):

                print(f"~~~~\tClean ENST\t~~~~")
                self._clean_enst = True

                head_response, tail_response = head_response[head_response["name"] == head_enst["name"]].squeeze(), tail_response[tail_response["name"] == tail_enst["name"]].squeeze()

                # Dumps all the blat data into the Blat class. Why am I doing that here? Because I didn't identify them until now. Now that they are fully identified it is worth it to 
                # put them into the Blat class
                self.head_blat = Blat(qName = self.hgene,
                                      tStart = head_response["tStart"],
                                      tEnd = head_response["tEnd"],
                                      blockCount = head_response["blockCount"],
                                      blockSizes = head_response["blockSizes"],
                                      tStarts = head_response["tStarts"])

                self.tail_blat = Blat(qName = self.tgene,
                                      tStart = tail_response["tStart"],
                                      tEnd = tail_response["tEnd"],
                                      blockCount = tail_response["blockCount"],
                                      blockSizes = tail_response["blockSizes"],
                                      tStarts = tail_response["tStarts"])

                # dumps all the ENST data needed into another class, the Gene class. The Gene class holds all the information about the gene. This way I can also contact UCSC in the gene class to get
                # additional data if need be.
                self.head_gene = Gene(self.hgene, # chrm = self.hchrm, strand = self.hstrand, ename=self.henst,
                                    gname = head_enst["name2"],
                                    txStart = head_enst["txStart"],
                                    txEnd = head_enst["txEnd"],
                                    cdsStart = head_enst["cdsStart"],
                                    cdsEnd = head_enst["cdsEnd"],
                                    exonCount = head_enst["exonCount"],
                                    exonStarts = head_enst["exonStarts"],
                                    exonEnds = head_enst["exonEnds"],
                                    exonFrames = head_enst["exonFrames"])

                self.tail_gene = Gene(self.tgene, # chrm = self.tchrm, strand = self.tstrand, ename=self.tenst,
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
                print("!!!!\tFailed to ENST\t!!!!")
                print(f"One did not ENST Identify Correctly\nHead Rows: {head_rows}\tTail Rows: {tail_rows}")
                self._clean_enst = False
         

        # else:
        #     print(f"One does not identify correctly\nHead Rows: {head_response.shape[0]}\tTail Rows: {tail_response.shape[0]}")

        #     print("Printing to Logs?")
        # # exit()




    def classify(self, adjacent_def: int = 50_000, unknown: str = "Unknown", cis: str = "C-SAG", inter: str = "T-E", intra: str = "T-A", tag_def: int = 100_000):
        '''
        Classify if it's Trans Inter/Intra Genic, Cis SAG. 

        Classify the type of splice as:
        C-SAG ~ Cis-SAG ~ same chr, same dir, adjecent genes -> last exon/UTR is < 50 kb to start of exon/UTR of next gene
        T-A ~ Trans-IntrAgenic ~ same chr, same/diff dir, non-adjecent genes
        T-E ~ trans-IntErgenic ~ diff chr
        
        Also look for TAD? The problem with TAD is it is defined in bases and therefore only really able to be found if it's TA. If it's TE, I don't
        spacially know where the chromosomes are.

        The last few clauses are redundant, since at initialization the classification is set to unknown, but I wanted to, for completion, include them.
        '''

        if self._clean_blat and self._clean_enst:

            if self.hchrm not in self.tchrm:
                self.classification = inter
                self.head2tailDistance = inter

            elif self.hchrm in self.tchrm:

                if self.hstrand not in self.tstrand:
                    self.classification = intra
                    
                    if self.hstrand in '-':
                        hposition = self.head_blat.tStarts[0]
                    elif self.hstrand in '+':
                        hposition = self.head_blat.tStarts[-1] + self.head_blat.blockSizes[-1]

                    if self.tstrand in '-':
                        tposition = self.tail_blat.tStarts[-1] + self.tail_blat.blockSizes[-1]
                    elif self.tstrand in '+':
                        tposition = self.tail_blat.tStarts[0]

                    self.head2tailDistance = complex(0, hposition - tposition)

                
                elif self.hstrand in self.tstrand:
                    if self.hstrand in '-':
                        hposition = self.head_blat.tStarts[0]
                    elif self.hstrand in '+':
                        hposition = self.head_blat.tStarts[-1] + self.head_blat.blockSizes[-1]

                    if self.tstrand in '-':
                        tposition = self.tail_blat.tStarts[-1] + self.tail_blat.blockSizes[-1]
                    elif self.tstrand in '+':
                        tposition = self.tail_blat.tStarts[0]

                    self.head2tailDistance = hposition - tposition

                    if abs(self.head2tailDistance) <= adjacent_def:
                        self.classification = f"{cis}"

                    elif abs(self.head2tailDistance) > adjacent_def:
                        self.classification = f"{intra}"

                    else:
                        self.classification = unknown
                
                else:
                    self.classification = unknown
            
            else:
                self.classification = unknown

        else:
            print("Cannot classify")
            self.classification = unknown


        

    def distance_measure(self):
        '''
        '''
        # distances, aistances = np.zeros((4, 1)), np.zeros((4, 1))
        distances, aistances = np.zeros((4)), np.zeros((4))

        if (self._clean_blat and self._clean_enst) and (self.hchrm in self.tchrm):
            head_positions: tuple = (self.head_blat.tStarts[0], self.head_blat.tStarts[-1] + self.head_blat.blockSizes[-1])
            tail_positions: tuple = (self.tail_blat.tStarts[0], self.tail_blat.tStarts[-1] + self.head_blat.blockSizes[-1])

            i = 0
            for hosition in head_positions:
                for tosition in tail_positions:
                    d = hosition - tosition
                    distances[i], aistances[i] = d, abs(d)
                    i += 1
            
            min_position = np.argmin(aistances)

            self.shortDistance = int(distances[min_position])

            if self.hstrand not in self.tstrand:
                self.shortDistance = complex(0, self.shortDistance)
            # print(f"Shortest Distance is {self.shortDistance}")

        elif (self._clean_blat and self._clean_enst) and (self.hchrm not in self.tchrm):
            self.shortDistance = "T-E"
        
        else:
            self.shortDistance = "Unknown"




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
                        enst_response = enst_response.squeeze()
                        blat_data.loc[row, "name"] = enst_response["name"]
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


    def find_junction(self, length: int = 10):
        '''
        To anyone who comes in after me to make adjusments: draw out what you're doing. UCSC Genome reports things in the + strand, so it's a little tricky

        I probably could clean this up a bit and turn it into some extra methods, but I'm having a hard time visualizing this without dealing with the personal
        shit going on, so I'm doing this very explicetly.

        Just need to grab the next expected exon and the slippage

        Looks like their is a rather major issue: some of these go right to the end of the exon, so I need to check that I can actually grab the thing I want.

        Move most of this over to the scoring class
        '''
        if self._clean_blat and self._clean_enst:
            head5primeExIn, tail3primeExIn = self._align_blat()

            if self.hstrand in "+":
                hEnd: int = self.head_blat.tStarts[-1] + self.head_blat.blockSizes[-1]
                hStart: int = hEnd - length
                h3prime_rev_dir: bool = True
                h5int_rev_dir: bool = False
                h5exo_rev_dir: bool = False
                
                hIntStart: int = hEnd
                hIntEnd: int = hEnd + 10

                head5primeExIn += 1

                hExoStart: int = self.head_gene.exonStarts[head5primeExIn]
                hExoEnd: int = hExoStart + length

            elif self.hstrand in "-":
                hStart: int = self.head_blat.tStarts[0]
                hEnd: int = hStart + length
                h3prime_rev_dir: bool = False
                h5int_rev_dir: bool = True
                h5exo_rev_dir: bool = True

                hIntEnd: int = hStart
                hIntStart: int = hStart - 10

                head5primeExIn -= 1

                hExoEnd: int = self.head_gene.exonEnds[head5primeExIn]
                hExoStart: int = hExoEnd - length

            if self.tstrand in "+":
                tStart: int = self.tail_blat.tStarts[0]
                tEnd: int = tStart + length
                t5prime_rev_dir: bool = False
                t3int_rev_dir: bool = True
                t3exo_rev_dir: bool = True

                tIntEnd: int = tStart
                tIntStart: int = tStart - 10

                tail3primeExIn -= 1

                tExoEnd: int = self.tail_gene.exonEnds[tail3primeExIn]
                tExoStart: int = tExoEnd - length
            
            elif self.tstrand in "-":
                tEnd: int = self.tail_blat.tStarts[-1] + self.tail_blat.blockSizes[-1]
                tStart: int = tEnd - length
                t5prime_rev_dir: bool = True
                t3int_rev_dir: bool = False
                t3exo_rev_dir: bool = False
                
                tIntEnd: int = tEnd + 10
                tIntStart: int = tEnd

                tail3primeExIn += 1

                tExoStart: int = self.tail_gene.exonStarts[tail3primeExIn]
                # chr1:94,495,984-94,495,993
                tExoEnd: int = tExoStart + length

            head3primeSeq, _ = upi.sequence(chrom = self.hchrm, start = hStart, end = hEnd, strand = self.hstrand, reverse_dir = h3prime_rev_dir)
            head5primeInt, _ = upi.sequence(chrom = self.hchrm, start = hIntStart, end = hIntEnd, strand = self.hstrand, reverse_dir = h5int_rev_dir)
            head5primeExo, _ = upi.sequence(chrom = self.hchrm, start = hExoStart, end = hExoEnd, strand = self.hstrand, reverse_dir = h5exo_rev_dir)

            tail5primeSeq, _ = upi.sequence(chrom = self.tchrm, start = tStart, end = tEnd, strand = self.tstrand, reverse_dir = t5prime_rev_dir)
            tail3primeInt, _ = upi.sequence(chrom = self.tchrm, start = tIntStart, end = tIntEnd, strand = self.tstrand, reverse_dir = t3int_rev_dir )
            tail3primeExo, _ = upi.sequence(chrom = self.tchrm, start = tExoStart, end = tExoEnd, strand = self.tstrand, reverse_dir = t3exo_rev_dir)

            # print(f"H3' Seq = {head3primeSeq}\n{url1}\nH5' Int = {head5primeInt}\n{url2}\nH5' Exo = {head5primeExo}\n{url3}\nT5' Seq = {tail5primeSeq}\n{url4}\nT3' Int = {tail3primeInt}\n{url5}\nT3' Exo = {tail3primeExo}\n{url6}")

            # print(f"H3' Seq = {head3primeSeq}\tT5' Seq = {tail5primeSeq}")
            # print(f"T3' Int = {tail3primeInt}\tH5' Int = {head5primeInt}")
            # print(f"T3' Exo = {tail3primeExo}\tH5' Exo = {head5primeExo}")

            # print(f"T5' Fus Seq = {tail5primeSeq}\tT3' Int Seq = {tail3primeInt}")
            # print(f"{urlT5}\n{urlT3}")

            toscore = {"H3_on_Tint3": [head3primeSeq, tail3primeInt],
                       "H3_on_Texo3": [head3primeSeq, tail3primeExo],
                       "T5_on_Hint5": [tail5primeSeq, head5primeInt],
                       "T5_on_Hexo5": [tail5primeSeq, head5primeExo]}
            
            self.disscores = DissSimilarityScore(print_score = True, **toscore)

            # for name, score in self.disscores:
            #     if not isinstance(score, np.ndarray):
            #         print(f"{name} = {score}")

            


    def query_attempt(self, func, *args, error_attempts: int = 5, **kwargs) -> Tuple[pandas.Series or pandas.DataFrame or Exception, str]:
        '''
        Returns the usable response from the UCSC Genome browser (typically a Pandas Series or DataFrame)
        along with the url used.
        '''
        # print(type(func))
        # print(*args)
        # print(kwargs)
        # exit()
        for attempt in range(error_attempts):
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


    def second_import(self, second_file: str or pathlib.Path = None, *args, **kwargs):
        '''
        For when things have already been identified and you just want to add details to the thing.
        '''


    def write_to_database(self, outfile: str or pathlib.Path, error_output: bool = False, *args, **kwargs):
        '''
        Writes to the error or output file

        Does not output None type, Bool type, or repetative names
        '''

        if isinstance(outfile, str):
            outfile: pathlib.Path = pathlib.Path(outfile)

        printable_row = pandas.Series(dtype=object)

        for key, value in vars(self).items():
            if ((value is not None) and (not isinstance(value, bool))):

                if isinstance(value, Gene):
                    indicator = "h" if value.name in self.hgene else "t"

                    for gey, galue in vars(value).items():

                        if (galue is not None) and (not isinstance(galue, bool)) and (gey not in "name"):
                            printable_row[f"{indicator}{gey}"] = galue
                            # print(f"{gey}: {galue}\nType: {type(galue)}")

                elif isinstance(value, Blat):
                    indicator = "h" if value.qName in self.hgene else "t"

                    for bey, balue in vars(value).items():

                        if(balue is not None) and (not isinstance(balue, bool)) and (bey not in "qName"):
                            printable_row[f"{indicator}{bey}"] = balue
                            # print(f"{bey}: {balue}\nType: {type(balue)}")

                elif isinstance(value, DissSimilarityScore):
                    for dey, dalue in vars(value).items():

                        if not isinstance(dalue, np.ndarray):
                            printable_row[f"{dey}"] = dalue

                else:
                    printable_row[f"{key}"] = value
                    # print(f"{key}: {value}\nType: {type(value)}")

        # print(printable_row)
        # print(printable_row.shape)

        if not outfile.is_file():
            pandas.DataFrame(columns=list(printable_row.index)).to_csv(outfile, header = True, index = False)


        printable_row.to_frame().T.to_csv(outfile, header = False, index = False, mode = 'a')



    def _align_blat(self) -> Tuple[int, int]:
        '''
        '''
        if self.hstrand in "+":
            h3prime = self.head_blat.tStarts[-1] + self.head_blat.blockSizes[-1]
            exon3primes = self.head_gene.exonEnds
        elif self.hstrand in "-":
            h3prime = self.head_blat.tStarts[0]
            exon3primes = self.head_gene.exonStarts

        if self.tstrand in "+":
            t5prime = self.tail_blat.tStarts[0]
            exon5primes = self.tail_gene.exonStarts
        elif self.tstrand in "-":
            t5prime = self.tail_blat.tStarts[-1] + self.tail_blat.blockSizes[-1]
            exon5primes = self.tail_gene.exonEnds

        hindex = self._local_align(h3prime, exon3primes)
        tindex = self._local_align(t5prime, exon5primes)

        # print(f"H Index: {hindex}")
        # print(f"T Index: {tindex}")

        return hindex, tindex


    def _local_align(self, blatPosition: int, exonPositions: tuple) -> int:
        '''

        '''
        # print(f"blatStart = {blatPosition}\nexonStarts = {exonStarts}")
        delta_blat = np.array(exonPositions)

        delta_blat = abs(delta_blat - blatPosition)

        min_delta = np.argmin(delta_blat)

        return min_delta


def main():
    '''
    '''


if __name__ in '__main__':
    main()