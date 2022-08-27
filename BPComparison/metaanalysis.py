import rpy2
import pandas
import pathlib
from RFunClass import RFun


# want min
# want max
# qaurtile of data
# Done by each column
# A lot of this can be done in ggplot

# import data to R
# apply catagories to Classification
# For each classification + global
    # y = column


def main(file_path = pathlib.Path or str):
    '''
    '''
    if isinstance(file_path, str):
        file_path = pathlib.Path(file_path)

    R = RFun()

    scoring_data: rpy2.robjects.vectors.DataFrame = R.read_csv(str(file_path), header=True)

    scoring_data = R.create_factor(scoring_data, "Classification")




if __name__ in '__main__':
    main(file_path = pathlib.Path.cwd().parent / "Data_Files" / "BPComp" / "Fusion_Scores_min100.csv")
