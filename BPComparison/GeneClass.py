import re


class Gene():
    def __init__(self, name: str, ename: str = None, gname: str = None,
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
        self.ename, self.gname = ename, gname
        self.chrm, self.strand = chrm, strand
        self.txStart, self.cdsStart, self.txEnd, self.cdsEnd = txStart, cdsStart, txEnd, cdsEnd
        self.exonCount, self.exonStarts, self.exonEnds, self.exonFrames = self._convert2numerical(exonCount), self._convert2numerical(exonStarts), self._convert2numerical(exonEnds), self._convert2numerical(exonFrames)



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



def main():
    '''
    '''


if __name__ in '__main__':
    main()