import pathlib
cwd = pathlib.Path.cwd()
import re



def the_most_import_method():
    '''
    Becuase it's the most important fucntion ever
    '''
    print(cwd)
    # hold: list = re.split(",", str)

if __name__ in "__main__":
    print(__name__)
    the_most_import_method()
