
import pandas as pd

def main():
    
    df_4pl = get_endpoint_titers('Endpoint_Titer_Results_4PL.tsv')
    df_5pl = get_endpoint_titers('Endpoint_Titer_Results_5PL.tsv')
    
    # PREP DATA FOR PLOTTING
    df_5pl.drop('g', axis=1, inplace=True)
    df_4pl['regression_type'] = '4PL'
    df_5pl['regression_type'] = '5PL'
    df = pd.concat([df_4pl, df_5pl])
    df['xtick_label'] = df['Group'] + ', ' + df['Sample ID']
    
    
    
    
def get_endpoint_titers(file):
    
    df = pd.read_csv(file, sep='\t')

    return df
    

if __name__ == '__main__':
    main()