
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
colors = sns.color_palette()


def main(args):

    data_file = args.data_file
    file_type = args.file_type.lower()
    
    validate_input(data_file, file_type)
    
    print('Gathering data...')
    df = get_data(data_file, file_type)
    
    category_cols = list(df.columns[:2])
    data_cols = list(df.columns[2:])

    df_long = df.melt(id_vars=category_cols, value_vars=data_cols, var_name='Individual', value_name='Absorbance')

    print('Fitting curves...')
    fit_df = fit_data(df_long)
    
    print('Making plots...')
    make_lineplots(df_long, 'Absorbance_vs_Dilution.jpg')
    make_lineplots_fitdata_v2(df_long, 'SigmoidFit_Absorbance_vs_Dilution.jpg', pd.DataFrame(fit_df))
    
    
def fit_data(df):
    
    fit_df = {'Groups':[],
            'Dilution':[],
            'Individual':[],
            'Absorbance':[]}
    categories = list(set(df['Groups']))
    individuals = list(set(df['Individual']))
    for cat in categories:
        for individual in individuals:
            single_df = df[ (df['Groups'] == cat) & (df['Individual'] == individual) ]
            popt, pcov = curve_fit(sigmoidal_fit, single_df['Dilution'], single_df['Absorbance'], maxfev=10000)

            num_points = 1000
            xvals = list(single_df['Dilution'])
            xrange = list(np.linspace(xvals[0], xvals[-1], num=num_points))
            yfit_vals = [sigmoidal_fit(x, *popt) for x in xrange]

            fit_df['Groups'] += [cat]*num_points
            fit_df['Individual'] += [individual]*num_points
            fit_df['Dilution'] += xrange
            fit_df['Absorbance'] += yfit_vals
            
    return fit_df
    

def sigmoidal_fit(x, a, b, c, d):
    fx = (a-d) / (1+((x/c)**b)) + d
    return fx
    
    
def make_lineplots(df, plot_name):
    
    df_copy = df.copy()
    df_copy['Dilution'] = np.log2( df_copy['Dilution'] )

    sns.set_style('darkgrid', rc={'xtick.bottom': True, 'ytick.left': True, 'xtick.color': 'black', 'ytick.color': 'black'})
    grid = sns.FacetGrid(df_copy, col='Individual', row='Groups', hue='Groups')

    grid.map(sns.scatterplot, 'Dilution', 'Absorbance')
    grid.map(sns.lineplot, 'Dilution', 'Absorbance')

    grid.set_titles(col_template="{col_name}", row_template="{row_name}")
    grid.fig.subplots_adjust(hspace=0.2) # Adjust spacing
    
    row_cats = []
    for cat in df_copy['Groups']:
        if cat not in row_cats:
            row_cats.append(cat)

    legend_elements = [Line2D([0], [0], marker='o', linestyle='None', markeredgecolor='None', label=cat, markerfacecolor=colors[i], markersize=10) for i, cat in enumerate(row_cats)]

    fig = plt.gcf()
    fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1,0.5), handletextpad=0.0, title='Groups')

    fig.set_size_inches(12, 8)

    plt.savefig(plot_name, bbox_inches='tight', dpi=600)
    plt.close()
    
    
def make_lineplots_fitdata_v2(df, plot_name, lines_df=None):
    
    df_copy = df.copy()
    df_copy['Dilution'] = np.log2( df_copy['Dilution'] )
    
    lines_df_copy = lines_df.copy()
    lines_df_copy['Dilution'] = np.log2( lines_df_copy['Dilution'] )
    
    df_copy['plot_type'] = 'scatter'
    lines_df_copy['plot_type'] = 'line'
    
    combined_df = pd.concat([df_copy, lines_df_copy], ignore_index=True)

    sns.set_style('darkgrid', rc={'xtick.bottom': True, 'ytick.left': True, 'xtick.color': 'black', 'ytick.color': 'black'})
    grid = sns.FacetGrid(df_copy, col='Individual', row='Groups', hue='Groups')
    
    grid.map(sns.scatterplot, 'Dilution', 'Absorbance', data=df_copy)

    row_cats = []
    col_cats = []
    for i, row_axes in enumerate(grid.axes):
        for j, ax in enumerate(row_axes):
            row_category = grid.row_names[i]
            col_category = grid.col_names[j]
            if row_category not in row_cats:
                row_cats.append(row_category)

            single_df = lines_df_copy[ (lines_df_copy['Groups'] == row_category) & (lines_df_copy['Individual'] == col_category)]
            sns.lineplot(x='Dilution', y='Absorbance', data=single_df, color='grey', ax=ax)
            ax.get_legend().remove()

    legend_elements = [Line2D([0], [0], marker='o', linestyle='None', markeredgecolor='None', label=cat, markerfacecolor=colors[i], markersize=10) for i, cat in enumerate(row_cats)]

    grid.set_titles(col_template="{col_name}", row_template="{row_name}")
    grid.fig.subplots_adjust(hspace=0.2) # Adjust spacing

    fig = plt.gcf()
    fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1,0.5), handletextpad=0.0, title='Groups')

    fig.set_size_inches(12, 8)
    plt.savefig(plot_name, bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_data(data_file, file_type):

    if file_type == 'csv':
        df = pd.read_csv(data_file)
    else:
        df = pd.read_csv(data_file, sep='\t')

    df['Dilution'] = pd.eval(df['Dilution']).astype(float)

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