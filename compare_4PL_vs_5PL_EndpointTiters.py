
import pandas as pd

def main():
    
    get_endpoint_titers('Endpoint_Titer_Results_4PL.tsv')
    
def get_endpoint_titers(file):
    
    df = pd.read_csv(file, sep='\t')

    return df
    

if __name__ == '__main__':
    main()