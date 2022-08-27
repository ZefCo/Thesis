import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import STAP
from rpy2.robjects.packages import importr


class RFun():
    def __init__(self) -> None:
        # Built in R Functions
        # self.function_name = robjects.r['function']
        # # Custom Functions
        # self.r = STAP(file_to_read, 'custom_function_name')
        ggplot = importr('gdata')

        self.read_csv = robjects.r["read.csv"]
        self.factor = robjects.r["factor"]
        self.levels = robjects.r["levels"]


    def create_factor(self, rframe: rpy2.robjects.vectors.DataFrame, column: str):
        '''
        '''
        col_index = list(rframe.colnames).index(column)
        rframe[col_index] = self.factor(rframe.rx2(column))

        return rframe


    def bar_and_whiskers(self, rframe: rpy2.robjects.vectors.DataFrame):
        '''
        '''

        



def main():
    '''
    '''

if __name__ in '__main__':
    main()