


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
        '''

        self.name = name
        self.ename, self.gname = ename, gname
        self.chrm, self.strand = chrm, strand
        self.txStart, self.cdsStart, self.txEnd, self.cdsEnd = txStart, cdsStart, txEnd, cdsEnd
        self.exonCount, self.exonStarts, self.exonEnds, self.exonFrames = exonCount, exonStarts, exonEnds, exonFrames



def main():
    '''
    '''


if __name__ in '__main__':
    main()