
class Blat():
    def __init__(self, 
                 matches: int = None, misMatches: int = None, repMatches: int = None, nCount: int = None, 
                 qNumInsert: int = None, qBaseInsert: int = None, tNumInsert: int = None, tBaseInsert: int = None, 
                 strand: str = None, qName: str = None, qSize: int = None, qStart: int = None, qEnd: int = None, 
                 tName: str = None, tSize: int = None, tStart: int = None, tEnd: int = None, 
                 blockCount: int = None, blockSizes: tuple = None, qStarts: tuple = None, tStarts: tuple = None) -> None:

        '''
        Not sure what I'm going to do with this other then just dump a lot of the attributes here and worry about them later
        '''
        self.matches = matches
        self.misMatches = misMatches
        self.repMatches = repMatches
        self.nCount = nCount
        self.qNumInsert = qNumInsert
        self.qBaseInsert = qBaseInsert
        self.tNumInsert = tNumInsert
        self.tBaseInsert = tBaseInsert
        self.strand = strand
        self.qName = qName
        self.qSize = qSize
        self.qStart = qStart
        self.qEnd = qEnd
        self.tName = tName
        self.tSize = tSize
        self.tStart = tStart
        self.tEnd = tEnd
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.qStarts = qStarts
        self.tStarts = tStarts
        self.repMatches = repMatches




def main():
    pass


if __name__ in '__main__':
    main()