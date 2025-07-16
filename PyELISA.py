
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np

def main(args):

    data_file = args.data_file
    file_type = args.file_type.lower()
    
    validate_input(data_file, file_type)
    df = get_data(data_file, file_type)
    

def get_data(data_file, file_type):

    if file_type == 'csv':
        df = pd.read_csv(data_file)
    else:
        df = pd.read_csv(data_file, sep='\t')

    df['Dilution'] = pd.eval(df['Dilution']).astype(float)
    df['Dilution'] = np.log2( df['Dilution'] )

    return df
    
    
def validate_input(data_file, file_type):
    
    error_message = ''
    # VALIDATE FILE
    try:
        h = open(data_file)
        h.close()
    except:
        error_message += '\nERROR: The file could not be found. Please ensure the file is in the same directory as PyELISA or specify an absolute file path.\n'
        
    if file_type not in ['csv', 'tsv']:
        error_message += '\nERROR: Unrecognized file type. Valid file types of "tsv" or "csv"'
        
    if error_message:
        print(error_message)
        print('\nExiting prematurely...')
        exit()


def get_args(arguments):
    parser = argparse.ArgumentParser(description='Calculate endpoint titer (ET) from ELISA data', prog='PyELISA')
    
    parser.add_argument('data_file', help="""Your data file.""")
    parser.add_argument('-t', '--file_type', type=str, default='tsv', help="""Your file type (csv or tsv). Default=tsv""")

    args = parser.parse_args(arguments)
    
    return args


if __name__ == '__main__':
    import sys, argparse
    args = get_args(sys.argv[1:])
    main(args)