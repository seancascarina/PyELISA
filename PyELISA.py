
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
    regression_type = args.regression_type.upper()
    threshold = args.threshold
    
    validate_input(data_file, file_type, regression_type, threshold)
    
    print('Gathering data...')
    df = get_data(data_file, file_type)
    
    category_cols = list(df.columns[:2])
    data_cols = list(df.columns[2:])

    df_long = df.melt(id_vars=category_cols, value_vars=data_cols, var_name='Individual', value_name='Absorbance')

    print('Fitting curves...')
    fit_df, et_df, categories, individuals = prep_containers(df_long)
    fit_df, et_df = fit_data(df_long, regression_type, threshold, fit_df, et_df, categories, individuals)

    print('Making plots...')
    make_lineplots(df_long, f'Absorbance_vs_Dilution_{regression_type}.jpg')
    make_lineplots_fitdata(df_long, f'SigmoidFit_Absorbance_vs_Dilution_{regression_type}.jpg', threshold, pd.DataFrame(fit_df))
    make_boxplot_endpoint_titers(et_df, regression_type)
    
    
def fit_data(df, regression_type, threshold, fit_df, et_df, categories, individuals):

    output = prep_output_file(regression_type)
    for cat in categories:
        for individual in individuals:
            single_df = df[ (df['Groups'] == cat) & (df['Individual'] == individual) ]

            num_points = 1000
            xvals = list(single_df['Dilution'])
            xrange = list(np.linspace(xvals[0], xvals[-1], num=num_points))
            
            if regression_type == '4PL':
                popt, pcov = curve_fit(fit_4PL, single_df['Dilution'], single_df['Absorbance'], maxfev=10000)
                yfit_vals = [fit_4PL(x, *popt) for x in xrange]
                a, b, c, d = popt
                if a > threshold:
                    print(f'\nWARNING: the sample "{individual}" from group "{cat}" has a curve-fit value at zero concentration that is higher than your threshold (threshold={threshold}), curve-fit absorbance at zero={a}. An endpoint titer could not be determined for this sample. Consider increasing the threshold or assaying more dilute samples to get lower absorbance readings.\n')
                    endpoint_titer = None
                else:
                    endpoint_titer = calc_endpoint_titer_4PL(a, b, c, d, threshold)
            else:
                popt, pcov = curve_fit(fit_5PL, single_df['Dilution'], single_df['Absorbance'], maxfev=10000)
                yfit_vals = [fit_5PL(x, *popt) for x in xrange]
                a, b, c, d, g = popt
                if a > threshold:
                    print(f'\nWARNING: the sample "{individual}" from group "{cat}" has a curve-fit value at zero concentration that is higher than your threshold (threshold={threshold}), curve-fit absorbance at zero={a}. An endpoint titer could not be determined for this sample. Consider increasing the threshold or assaying more dilute samples to get lower absorbance readings.\n')
                    endpoint_titer = None
                else:
                    endpoint_titer = calc_endpoint_titer_5PL(a, b, c, d, g, threshold)

            fit_df['Groups'] += [cat]*num_points
            fit_df['Individual'] += [individual]*num_points
            fit_df['Dilution'] += xrange
            fit_df['Absorbance'] += yfit_vals
            
            if endpoint_titer:
                et_reciprocal = 1 / endpoint_titer
                et_df['Endpoint titer'].append(1/et_reciprocal)
                et_df['Groups'].append(cat)
            else:
                et_reciprocal = None
                
            output.write('\t'.join([cat, individual, str(endpoint_titer), str(et_reciprocal)] + [str(x) for x in popt]) + '\n')
            
    output.close()
    
    return fit_df, et_df
    

def fit_5PL(x, a, b, c, d, g):
    fx = (a-d) / (1+((x/c)**b))**g + d
    return fx
    
    
def calc_endpoint_titer_5PL(a, b, c, d, g, threshold):
    endpoint_titer = c * (( (a-d) / (threshold-d) )**(1/g) - 1)**(1/b)
    return endpoint_titer
    

def fit_4PL(x, a, b, c, d):
    fx = (a-d) / (1+((x/c)**b)) + d
    return fx
    
    
def calc_endpoint_titer_4PL(a, b, c, d, threshold):
    endpoint_titer = c * (( (a-d) / (threshold-d) ) - 1)**(1/b)
    return endpoint_titer
    
    
def prep_containers(df):
    
    fit_df = {'Groups':[],
            'Dilution':[],
            'Individual':[],
            'Absorbance':[]}
    et_df = {'Endpoint titer':[],
            'Groups':[]}

    categories = []
    for cat in df['Groups']:
        if cat not in categories:
            categories.append(cat)
            
    individuals = list(set(df['Individual']))
    
    return fit_df, et_df, categories, individuals
    
    
def prep_output_file(regression_type):
    
    if regression_type == '4PL':
        header = '\t'.join(['Group', 'Sample ID', 'Endpoint Titer', '1 / Endpoint_Titer'] + list('abcd'))
    else:
        header = '\t'.join(['Group', 'Sample ID', 'Endpoint Titer', '1 / Endpoint_Titer'] + list('abcdg'))
        
    output = open('Endpoint_Titer_Results.tsv', 'w')
    output.write(header + '\n')
    
    return output
    
    
def make_boxplot_endpoint_titers(df, regression_type):
    
    sns.boxplot(x='Groups', y='Endpoint titer', data=df, hue='Groups', showfliers=False, palette=['0.7']*len(set(df['Groups'])))
    sns.stripplot(x='Groups', y='Endpoint titer', data=df, hue='Groups')
    plt.savefig(f'Endpoint_titer_Boxplots_{regression_type}.tif', bbox_inches='tight', dpi=600)
    plt.close()
    
    
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
    
    
def make_lineplots_fitdata(df, plot_name, threshold, lines_df=None):
    
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
            sns.lineplot(x=(min(single_df['Dilution']), max(single_df['Dilution'])), y=(threshold, threshold), linestyle='--', color='0.5', ax=ax)
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
    
    
def validate_input(data_file, file_type, regression_type, threshold):
    
    error_message = ''
    # VALIDATE FILE
    try:
        h = open(data_file)
        h.close()
    except:
        error_message += '\nERROR: The file could not be found. Please ensure the file is in the same directory as PyELISA or specify an absolute file path.\n'
        
    if file_type not in ['csv', 'tsv']:
        error_message += '\nERROR: Unrecognized file type. Valid file types are "tsv" or "csv" only.'
        
    if regression_type not in ['4PL', '5PL']:
        error_message += '\nERROR: Unrecognized value for regression_type parameter. Valid regression types are "4PL" or "5PL" only.'
        
    if not 0.0 <= threshold:
        error_message += '\nERROR: Threshold value for calculating endpoint titer should be > 0.'
        
    if error_message:
        print(error_message)
        print('\nExiting prematurely...')
        exit()


def get_args(arguments):
    parser = argparse.ArgumentParser(description='Calculate endpoint titer (ET) from ELISA data', prog='PyELISA')
    
    parser.add_argument('data_file', help="""Your data file.""")
    parser.add_argument('-f', '--file_type', type=str, default='csv', help="""The file type of the data_file provided. Options are "csv" or "tsv" Default=csv""")
    parser.add_argument('-r', '--regression_type', type=str, default='4PL', help="""The regression type used for curve fitting. Options are "4PL" or "5PL" (4-parameter logistic or 5-parameter logistic, respectively). Default=4PL""")
    parser.add_argument('-t', '--threshold', type=float, default=0.2, help="""The threshold to use for endpoint titer calculation. Default=0.2""")

    args = parser.parse_args(arguments)
    
    return args


if __name__ == '__main__':
    import sys, argparse
    args = get_args(sys.argv[1:])
    main(args)