import pathlib
import numpy as np
import ucsc_restapi as upi
from FusionClass import FusionGene
import pandas

class DissSimilarityScore():
    '''
    Right now you just pass the fusion Class into this and call the parts you want. Maybe I should come up with a less specific way of doing this but... I'm not.
    '''
    def __init__(self, fusion: FusionGene, length: int = 10) -> None:

        # All of these are str type
        self.H3_Seq, self.H5_Int, self.H5_Exo, self.T5_Seq, self.T3_Int, self.T3_Exo = None, None, None, None, None, None
        self.sequences, self.scores = None, None

        self.fusion: FusionGene = fusion

        self._length: int = length
        self._weights()


    def find_junction(self):
        '''
        To anyone who comes in after me to make adjusments: draw out what you're doing. UCSC Genome reports things in the + strand, so it's a little tricky

        This establishes the sequences to be scored.

        '''
        head5primeIndex, tail3primeIndex = self.fusion._align_blat()
        # print(self.fusion.hstrand, head5primeIndex - 1, head5primeIndex, head5primeIndex + 1, self.fusion.head_gene.exonCount)
        # print(self.fusion.tstrand, tail3primeIndex - 1, tail3primeIndex, tail3primeIndex + 1, self.fusion.tail_gene.exonCount)

        if (((self.fusion.hstrand in "+") and ((head5primeIndex + 1) <= self.fusion.head_gene.exonCount)) or ((self.fusion.hstrand in "-") and ((head5primeIndex - 1) >= 0))):
            
            # print(tail3primeIndex, head5primeIndex)
            # print(self.fusion.tail_gene.exonCount, self.fusion.head_gene.exonCount)

            if self.fusion.hstrand in "+":
                hEnd: int = self.fusion.head_blat.tStarts[-1] + self.fusion.head_blat.blockSizes[-1]
                hStart: int = hEnd - self._length
                h3prime_rev_dir: bool = True
                h5int_rev_dir: bool = False
                h5exo_rev_dir: bool = False
                
                hIntStart: int = hEnd
                hIntEnd: int = hEnd + self._length

                if head5primeIndex + 1 < self.fusion.head_gene.exonCount:
                    head5primeIndex += 1

                hExoStart: int = self.fusion.head_gene.exonStarts[head5primeIndex]
                hExoEnd: int = hExoStart + self._length

            elif self.fusion.hstrand in "-":
                hStart: int = self.fusion.head_blat.tStarts[0]
                hEnd: int = hStart + self._length
                h3prime_rev_dir: bool = False
                h5int_rev_dir: bool = True
                h5exo_rev_dir: bool = True

                hIntEnd: int = hStart
                hIntStart: int = hStart - self._length

                if head5primeIndex - 1 > self.fusion.head_gene.exonCount:
                    head5primeIndex -= 1

                hExoEnd: int = self.fusion.head_gene.exonEnds[head5primeIndex]
                hExoStart: int = hExoEnd - self._length

            # all of these are str. The _ is also a str type, it's the URL, but the URL right now is not needed.
            self.H3_Seq, _ = upi.sequence(chrom = self.fusion.hchrm, start = hStart, end = hEnd, strand = self.fusion.hstrand, reverse_dir = h3prime_rev_dir)
            self.H5_Int, _ = upi.sequence(chrom = self.fusion.hchrm, start = hIntStart, end = hIntEnd, strand = self.fusion.hstrand, reverse_dir = h5int_rev_dir)
            self.H5_Exo, _ = upi.sequence(chrom = self.fusion.hchrm, start = hExoStart, end = hExoEnd, strand = self.fusion.hstrand, reverse_dir = h5exo_rev_dir)

        else:
            self.H3_Seq, self.H5_Int, self.H5_Exo = "", "", ""

            for _ in range(self._length):
                self.H3_Seq: str = f"{self.H3_Seq}x"
                self.H5_Int: str = f"{self.H5_Int}x"
                self.H5_Exo: str = f"{self.H5_Exo}x"


        if (((self.fusion.tstrand in "+") and ((tail3primeIndex - 1) >= 0)) or ((self.fusion.tstrand in "-") and ((tail3primeIndex + 1) <= self.fusion.tail_gene.exonCount))):
            if self.fusion.tstrand in "+":
                tStart: int = self.fusion.tail_blat.tStarts[0]
                tEnd: int = tStart + self._length
                t5prime_rev_dir: bool = False
                t3int_rev_dir: bool = True
                t3exo_rev_dir: bool = True

                tIntEnd: int = tStart
                tIntStart: int = tStart - self._length

                if tail3primeIndex -1 > self.fusion.tail_gene.exonCount:
                    tail3primeIndex -= 1

                tExoEnd: int = self.fusion.tail_gene.exonEnds[tail3primeIndex]
                tExoStart: int = tExoEnd - self._length
            
            elif self.fusion.tstrand in "-":
                tEnd: int = self.fusion.tail_blat.tStarts[-1] + self.fusion.tail_blat.blockSizes[-1]
                tStart: int = tEnd - self._length
                t5prime_rev_dir: bool = True
                t3int_rev_dir: bool = False
                t3exo_rev_dir: bool = False
                
                tIntEnd: int = tEnd + self._length
                tIntStart: int = tEnd

                if tail3primeIndex + 1 < self.fusion.tail_gene.exonCount:
                    tail3primeIndex += 1

                tExoStart: int = self.fusion.tail_gene.exonStarts[tail3primeIndex]
                # chr1:94,495,984-94,495,993
                tExoEnd: int = tExoStart + self._length

            self.T5_Seq, _ = upi.sequence(chrom = self.fusion.tchrm, start = tStart, end = tEnd, strand = self.fusion.tstrand, reverse_dir = t5prime_rev_dir)
            self.T3_Int, _ = upi.sequence(chrom = self.fusion.tchrm, start = tIntStart, end = tIntEnd, strand = self.fusion.tstrand, reverse_dir = t3int_rev_dir )
            self.T3_Exo, _ = upi.sequence(chrom = self.fusion.tchrm, start = tExoStart, end = tExoEnd, strand = self.fusion.tstrand, reverse_dir = t3exo_rev_dir)

        else:
            self.T5_Seq, self.T3_Int, self.T3_Exo = "", "", ""

            for _ in range(self._length):
                self.T5_Seq: str = f"{self.T5_Seq}x"
                self.T3_Int: str = f"{self.T3_Int}x"
                self.T3_Exo: str = f"{self.T3_Exo}x"

        self.sequences: dict = {"H3_Seq": self.H3_Seq,
                                "H5_Intron": self.H5_Int,
                                "H5_Exon": self.H5_Exo,
                                "T5_Seq": self.T5_Seq,
                                "T3_Intron": self.T3_Int,
                                "T3_Exon": self.T3_Exo}


    def score(self):
        '''
        '''
        none_check = lambda first, second: True if (first is not None) and (second is not None) else False

        if self.sequences is None:
            self.find_junction()

        # print(self.sequences)

        self.scores: dict = {}
        # print(self.sequences)

        for i, (keyi, seqi) in enumerate(self.sequences.items()):
            for j, (keyj, seqj) in enumerate(self.sequences.items()):
                
                key = f"{keyi}_{keyj}"
                reverse_key = f"{keyj}_{keyi}"

                if i == j:
                    pass
                
                elif none_check(seqi, seqj):

                    if reverse_key not in self.scores.keys():
                        local_g = self._little_g(seqi, seqj)
                        try:
                            local_score = np.sum(self._class__weight * local_g)
                        except ValueError as e:
                            local_score = 1
                            print(f"{keyi} -> {seqi}\t{keyj} -> {seqj}")
                        except Exception as e:
                            local_score = 1
                            print(f"!!!!\tNew Error\t!!!!")
                            print(f"    \tError: {e}\tType: {type(e)}")
                            print(f"{keyi} -> {seqi}\t{keyj} -> {seqj}")

                        self.scores[key] = local_score

                else:
                    self.scores[key] = 1



    def write_score(self, outfile: str or pathlib.Path):
        '''
        '''
        if isinstance(outfile, str):
            outfile: pathlib.Path = pathlib.Path(outfile)

        printable_row = pandas.Series(dtype=object)
        printable_row["Name"] = f"{self.fusion.hgene}_{self.fusion.tgene}"
        printable_row["ENST"] = f"{self.fusion.henst}_{self.fusion.tenst}"
        printable_row["HStrand"], printable_row["TStrand"] = self.fusion.hstrand, self.fusion.tstrand
        printable_row["HChr"], printable_row["TChr"] = self.fusion.hchrm, self.fusion.tchrm
        
        for key, value in self.sequences.items():
            printable_row[key] = value

        for key, value in self.scores.items():
            printable_row[key] = value

        printable_row["Classification"], printable_row["ShortDistance"], printable_row["Head2TailDistance"] = self.fusion.classification, self.fusion.shortDistance, self.fusion.head2tailDistance

        if not outfile.is_file():
            pandas.DataFrame(columns=list(printable_row.index)).to_csv(outfile, header = True, index = False)

        printable_row.to_frame().T.to_csv(outfile, index = False, header = False, mode = 'a')




    def _little_g(self, sequence1: str, sequence2: str):
        '''
        '''
        delta_len = len(sequence1) - len(sequence2)

        if delta_len > 0:
            for _ in range(delta_len):
                sequence2 = f"{sequence2}x"

        if delta_len < 0:
            for _ in range(abs(delta_len)):
                sequence1 = f"{sequence1}x"

        g_vec = np.array([0 if seq1 in sequence2[i] else 1 for i, seq1 in enumerate(sequence1)])

        return g_vec



    def _weights(self):
        '''
        Weights are caculated as Sum[1/2**i] from i = 1 to L. There is nothing saying you can't do this from 0 instead of 1, but the advantage of using 1 as your start index
        is that the score is scaled from 0 to 1. If you started from 1 the score is scaled from 0 to 2. No change in physical meaning, just adjusts your R value.
        '''
        self._class__weight = np.array([1/(2**(i + 1)) for i in range(self._length)])



def main():
    '''
    '''


if __name__ in "__main__":
    main()