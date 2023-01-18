import re
import ucsc_restapi as upi
import numpy as np
import pathlib
import pandas


class Gene():
    def __init__(self, name: str, ename: str = None, gname: str = None, ncibname: str = None,
                       chrm: str = None, 
                       strand: str = None, 
                       txStart: int = None, 
                       txEnd: int = None, 
                       cdsStart: int = None, 
                       cdsEnd: int = None, 
                       exonCount: int = None, 
                       exonStarts: tuple or list = None, 
                       exonEnds: tuple or list = None, 
                       exonFrames: tuple or list = None) -> None:

        '''
        For simplicity this is based on the ENST Track of UCSC genome browser
        Converts exonStarts	exonEnds exonFrames into numerical tuples
        '''

        self.name = name
        self.ename, self.gname, self.ncibname = ename, gname, ncibname
        self.chrm, self.strand = chrm, strand
        self.txStart, self.cdsStart, self.txEnd, self.cdsEnd = self._convert2int(txStart), self._convert2int(cdsStart), self._convert2int(txEnd), self._convert2int(cdsEnd)
        self.exonCount, self.exonStarts, self.exonEnds, self.exonFrames = self._convert2numerical(exonCount), self._convert2numerical(exonStarts), self._convert2numerical(exonEnds), self._convert2numerical(exonFrames)
        self.utr5_cords, self.cds_cords, self.utr3_cords, self.intron_cords = None, None, None, None
        self.utr5_seq, self.cds_seq, self.utr3_seq, self.intron_seq = None, None, None, None




    def sequence_breakdown(self):
        '''
        '''

        print(f"Breakdown for {self.name}")
        # print(f"\n####\n{self.name}\n####\n")
        # print(f"Chrom: {self.chrm}\tStrand: {self.strand}")
        if (self.utr3_cords is None) or (self.utr5_cords is None) or (self.cds_cords is None) or (self.intron_cords is None):
            self._sequence_coords()

        if len(self.utr5_cords) > 0:
            # print(f"UTR Seqs")
            self.utr5_seq = []
            for utr5 in self.utr5_cords:
                seq, _ = upi.sequence(chrom=self.chrm, start = utr5[0], end = utr5[1], strand = self.strand)
                self.utr5_seq.append(seq)
                # print(url)

        if len(self.cds_cords) > 0:
            # print(f"CDS Seqs")
            self.cds_seq = []
            for cds in self.cds_cords:
                seq, _ = upi.sequence(chrom = self.chrm, start = cds[0], end = cds[1], strand = self.strand)
                self.cds_seq.append(seq)
                # print(url)
        
        if len(self.utr3_cords) > 0:
            # print(f"UTR Seqs")
            self.utr3_seq = []
            for utr3 in self.utr3_cords:
                seq, _ = upi.sequence(chrom = self.chrm, start = utr3[0], end = utr3[1], strand = self.strand)
                self.utr3_seq.append(seq)
                # print(url)

        if len(self.intron_cords) > 0:
            # print(f"Intron Seqs")
            self.intron_seq = []
            for intron in self.intron_cords:
                seq, _ = upi.sequence(chrom = self.chrm, start = intron[0], end = intron[1], strand = self.strand)
                self.intron_seq.append(seq)
                # print(url)




    def write_sequences(self, outfolder: str or pathlib.Path):
        '''
        Writes the UTR, CDS, and Intron sequences into 4 differenct files
        '''

        print(f"Writing {self.name} Sequences")
        if isinstance(outfolder, str):
            outfolder: pathlib.Path = pathlib.Path(outfolder)

        if not outfolder.is_dir():
            outfolder.mkdir()

        printable_row = pandas.Series(dtype=object)
        # name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
        # This may print out None if there is nothing in there, but that's OK
        printable_row["Name"] = self.name
        printable_row["ename"], printable_row["gname"], printable_row["ncibname"] = self.ename, self.gname, self.ncibname
        printable_row["Chr"], printable_row["Strand"] = self.chrm, self.strand

        utr5_headers, cds_headers, utr3_headers, intron_headers = ["Name", "ename", "gname", "ncibname", "Chr", "Strand"], ["Name", "ename", "gname", "ncibname", "Chr", "Strand"], ["Name", "ename", "gname", "ncibname", "Chr", "Strand"], ["Name", "ename", "gname", "ncibname", "Chr", "Strand"]

        if isinstance(self.utr5_seq, list):
            for j, seq in enumerate(self.utr5_seq):
                header_seq, header_cor = f"UTR5Seq{j}", f"UTR5Cord{j}"
                utr5_headers.append(header_seq), utr5_headers.append(header_cor)

                printable_row[header_cor] = self.utr5_cords[j]
                printable_row[header_seq] = seq
            
            utr5_frame = printable_row[utr5_headers]
            utr5_file = outfolder / "UTR5_seq.csv"
            # print(utr5_frame)

            if not utr5_file.is_file():
                pandas.DataFrame(columns=list(printable_row.index)).to_csv(utr5_file, header = False, index = False)

            utr5_frame.to_frame().T.to_csv(utr5_file, index = False, header = False, mode = 'a')
            
        if isinstance(self.cds_seq, list):
            for j, seq in enumerate(self.cds_seq):
                header_seq, header_cor = f"CDSSeq{j}", f"CDSCord{j}"
                cds_headers.append(header_seq), cds_headers.append(header_cor)

                printable_row[header_cor] = self.cds_cords[j]
                printable_row[header_seq] = seq

            cds_frame = printable_row[cds_headers]
            cds_file = outfolder / "CDS_seq.csv"
            # print(cds_frame)

            if not cds_file.is_file():
                pandas.DataFrame(columns=list(printable_row.index)).to_csv(cds_file, header = False, index = False)

            cds_frame.to_frame().T.to_csv(cds_file, index = False, header = False, mode = 'a')

            
        if isinstance(self.utr3_seq, list):
            for j, seq in enumerate(self.utr3_seq):
                header_seq, header_cor = f"UTR3Seq{j}", f"UTR3Cord{j}"
                utr3_headers.append(header_seq), utr3_headers.append(header_cor)

                printable_row[header_cor] = self.utr3_cords[j]
                printable_row[header_seq] = seq

            utr3_frame = printable_row[utr3_headers]
            utr3_file = outfolder / "UTR3_seq.csv"
            # print(utr3_frame)

            if not utr3_file.is_file():
                pandas.DataFrame(columns=list(printable_row.index)).to_csv(utr3_file, header = False, index = False)

            utr3_frame.to_frame().T.to_csv(utr3_file, index = False, header = False, mode = 'a')


        if isinstance(self.intron_seq, list):
            for j, seq in enumerate(self.intron_seq):
                header_seq, header_cor = f"IntronSeq{j}", f"IntronCord{j}"
                intron_headers.append(header_seq), intron_headers.append(header_cor)

                printable_row[header_cor] = self.intron_cords[j]
                printable_row[header_seq] = seq

            intron_frame = printable_row[intron_headers]
            intron_file = outfolder / "Intron.csv"
            # print(intron_frame)

            if not intron_file.is_file():
                pandas.DataFrame(columns=list(printable_row.index)).to_csv(intron_file, header = False, index = False)

            intron_frame.to_frame().T.to_csv(intron_file, index = False, header = False, mode = 'a')


    def _sequence_coords(self):
        '''
        Gets the sequences for the UTR, CDS, and Introns
        '''

        delta, delta_prime = self._local_align()

        # print(delta, delta_prime)

        self.utr5_cords = []
        self.cds_cords = []
        self.utr3_cords = []
        self.intron_cords = []

        # if delta == delta_prime -> only UTR
        # else -> cds

        if delta_prime == delta:
            for j, start_cord in enumerate(self.exonStarts):
                self.utr5_cords.append((start_cord, self.exonEnds[j]))


        elif self.strand in "+":
            # print("#### Sense Strand ####")
            for j, start_cord in enumerate(self.exonStarts):
                if (j < delta) and (j < delta_prime):  # only in the utr region
                    self.utr5_cords.append((start_cord, self.exonEnds[j]))

                elif (j == delta) and (j < delta_prime):  # crossing over from the utr region to the cds region
                    self.utr5_cords.append((start_cord, self.cdsStart))
                    self.cds_cords.append((self.cdsStart, self.exonEnds[j]))
                
                elif (j > delta) and (j < delta_prime):  # only in the cds region
                    self.cds_cords.append((start_cord, self.exonEnds[j]))

                elif (j > delta) and (j == delta_prime):  # crossing over form the cds to the utr
                    self.cds_cords.append((start_cord, self.cdsEnd))
                    self.utr3_cords.append((self.cdsEnd, self.exonEnds[j]))

                elif (j > delta) and (j > delta_prime):  # only in the utr region
                    self.utr3_cords.append((start_cord, self.exonEnds[j]))

                elif (j == delta) and (j == delta_prime):  # this happens if everything is in the utr 5'
                    self.utr5_cords.append((start_cord, self.exonEnds[j]))

        
        elif self.strand in "-":
            # print("#### AntiSense Strand ####")
            for j, start_cord in enumerate(self.exonStarts):
                if (j < delta) and (j < delta_prime):  # only in the utr region
                    self.utr3_cords.append((start_cord, self.exonEnds[j]))

                elif (j == delta) and (j < delta_prime):  # crossing over from the utr region to the cds region
                    self.utr3_cords.append((start_cord, self.cdsStart))
                    self.cds_cords.append((self.cdsStart, self.exonEnds[j]))
                
                elif (j > delta) and (j < delta_prime):  # only in the cds region
                    self.cds_cords.append((start_cord, self.exonEnds[j]))

                elif (j > delta) and (j == delta_prime):  # crossing over form the cds to the utr
                    self.cds_cords.append((start_cord, self.cdsEnd))
                    self.utr5_cords.append((self.cdsEnd, self.exonEnds[j]))

                elif (j > delta) and (j > delta_prime):  # only in the utr region
                    self.utr5_cords.append((start_cord, self.exonEnds[j]))

                elif (j == delta) and (j == delta_prime):  # this happens if everything is in the utr 5'
                    self.utr5_cords.append((start_cord, self.exonEnds[j]))


        for i in range(len(self.exonEnds) - 1):
            self.intron_cords.append((self.exonEnds[i], self.exonStarts[i + 1])) 



    def _local_align(self):
        '''
        Finds the minimum difference between one set of coordinates and a specific location.
        '''
        delta_starts = np.array(self.exonStarts)
        delta_ends = np.array(self.exonEnds)

        delta, delta_prime = abs(delta_ends - self.cdsStart), abs(delta_starts - self.cdsEnd)
        delta, delta_prime = np.argmin(delta), np.argmin(delta_prime)

        return delta, delta_prime     



    def _convert2numerical(self, sequence: str):
        '''
        Because on some of these responses come back as lists, but their not type(list), their type(str). This converts
        them to a list. But it converts EVERYTHING in that column to a list for simplicity, so an entry of 1 is not in a
        list of 1. Trust me, in the long run this works out for the best.

        This was originally in the Blat API, but honestly I need to use it in several places so I moved it here. Easier access
        for everyone.

        OK technically this is a tuple, but list, tuple, set, tomato, potatoe
        
        Yes I know what the actually phrase is!
        '''
        if isinstance(sequence, str):
            # print(sequence, type(sequence))
            sequence = re.sub(r"\(|\)|\,$", "", sequence)
            # print(sequence, type(sequence))
            # exit()

            sequence = re.split(',', sequence)

            seqlist = []
            for strint in sequence:
                
                try:
                    strint = int(strint)

                except ValueError as e:
                    # print("!!! Value Error !!!")
                    # print(f"Error: {e}\tType: {type(e)}\nInput: {strint}")
                    strint = None

                except Exception as e:
                    # print("!!! Other Error !!!")
                    # print(f"Error: {e}\tType: {type(e)}\nInput: {strint}")
                    strint = None

                seqlist.append(strint)

            seqlist = tuple(seqlist)

        else:
            seqlist = sequence

        return seqlist


    def _convert2int(self, att: str) -> int:
        '''
        For when the thing is just a number.
        '''

        att = int(att)

        return att


    def _get_cds(self):
        '''
        Somehow we're going to get the cds sequence
        '''



def main():
    '''
    '''


if __name__ in '__main__':
    main()