
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    
    df_4pl = get_endpoint_titers('Endpoint_Titer_Results_4PL.tsv')
    df_5pl = get_endpoint_titers('Endpoint_Titer_Results_5PL.tsv')
    
    # PREP DATA FOR PLOTTING
    df_5pl.drop('g', axis=1, inplace=True)
    df_4pl['regression_type'] = '4PL'
    df_5pl['regression_type'] = '5PL'
    df = pd.concat([df_4pl, df_5pl], ignore_index=True)
    df['xtick_label'] = df['Group'] + ', ' + df['Sample ID']
    
    plotting(df)
    
    
def plotting(df):
    
    sns.barplot(x='xtick_label', y='1 / Endpoint_Titer', data=df, hue='regression_type')
    plt.xticks(rotation=90)
    plt.xlabel('Group, Sample ID')
    fig = plt.gcf()
    fig.set_size_inches(7,5)
    plt.savefig('4PL_vs_5PL_EndpointTiter_Comparison.tif', bbox_inches='tight', dpi=300)
    plt.close()
    
    
def get_endpoint_titers(file):
    
    df = pd.read_csv(file, sep='\t')

    return df
    

if __name__ == '__main__':
    main()