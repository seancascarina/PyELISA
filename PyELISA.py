
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
colors = sns.color_palette()


def main(args):

    # EXTRACT COMMAND-LINE ARGUMENTS
    data_file = args.data_file
    file_type = args.file_type.lower()
    regression_type = args.regression_type.upper()
    threshold = args.threshold
    maxfev = args.maxfev
    
    validate_input(data_file, file_type, regression_type, threshold)
    
    print('Gathering data...')
    df = get_data(data_file, file_type)
    
    # SEPARATE GROUP AND SAMPLE_ID COLUMNS FROM DATA COLUMNS
    category_cols = list(df.columns[:2])
    data_cols = list(df.columns[2:])

    # CONVERT WIDE DICTIONARY TO LONG-FORM
    df_long = df.melt(id_vars=category_cols, value_vars=data_cols, var_name='Individual', value_name='Absorbance')

    print('Fitting curves...')
    fit_df, et_df, categories, individuals = prep_containers(df_long)
    fit_df, et_df = fit_data(df_long, regression_type, threshold, fit_df, et_df, categories, individuals, maxfev)

    print('Making plots...')
    make_lineplots(df_long, categories, f'Absorbance_vs_Dilution_{regression_type}.jpg')
    make_lineplots_fitdata(df_long, f'SigmoidFit_Absorbance_vs_Dilution_{regression_type}.jpg', threshold, pd.DataFrame(fit_df))
    make_boxplot_endpoint_titers(et_df, regression_type)
    
    
def fit_data(df, regression_type, threshold, fit_df, et_df, categories, individuals, maxfev):
    """
    Fit Absorbance vs. Dilution ELISA data with a 4-parameter logistic (4PL) or
    5-parameter logistic (5PL) regression model (user's choice).
    
    Returns:
        fit_df (dict) - Long-form dictionary with 1000 fit data points per 
            sample ID representing 4PL or 5PL curvefit values.
        et_df (dict) - Long-form dictionary with endpoint titer values
            for each group (for boxplot).
    """
    
    output = prep_output_file(regression_type)
    for cat in categories:
        for individual in individuals:
            
            # EXTRACT DATA FOR A SINGLE SAMPLE_ID ("individual")
            single_df = df[ (df['Groups'] == cat) & (df['Individual'] == individual) ]

            # DEFINE X-VALUES TO USE FOR CURVE FITTING
            num_points = 1000
            xvals = list(single_df['Dilution'])
            xrange = list(np.linspace(xvals[0], xvals[-1], num=num_points))
            
            # PERFORM CURVE-FITTING FOR 4PL
            if regression_type == '4PL':
                popt, pcov = curve_fit(fit_4PL, single_df['Dilution'], single_df['Absorbance'], maxfev=maxfev)
                yfit_vals = [fit_4PL(x, *popt) for x in xrange]
                a, b, c, d = popt   # EXTRACT OPTIMIZED REGRESSION PARAMETERS
                r_squared = calc_rsquared(single_df['Dilution'], single_df['Absorbance'], popt, regression_type)
                if a > threshold:
                    print(f'\nWARNING: the sample "{individual}" from group "{cat}" has a curve-fit value at zero concentration that is higher than your threshold (threshold={threshold}), curve-fit absorbance at zero={a}. An endpoint titer could not be determined for this sample. Consider increasing the threshold or assaying more dilute samples to get lower absorbance readings.\n')
                    endpoint_titer = None
                else:
                    endpoint_titer = calc_endpoint_titer_4PL(a, b, c, d, threshold)
            # OR PERFORM CURVE-FITTING FOR 5PL
            else:
                popt, pcov = curve_fit(fit_5PL, single_df['Dilution'], single_df['Absorbance'], maxfev=maxfev)
                yfit_vals = [fit_5PL(x, *popt) for x in xrange]
                a, b, c, d, g = popt   # EXTRACT OPTIMIZED REGRESSION PARAMETERS
                r_squared = calc_rsquared(single_df['Dilution'], single_df['Absorbance'], popt, regression_type)
                if a > threshold:
                    print(f'\nWARNING: the sample "{individual}" from group "{cat}" has a curve-fit value at zero concentration that is higher than your threshold (threshold={threshold}), curve-fit absorbance at zero={a}. An endpoint titer could not be determined for this sample. Consider increasing the threshold or assaying more dilute samples to get lower absorbance readings.\n')
                    endpoint_titer = None
                else:
                    endpoint_titer = calc_endpoint_titer_5PL(a, b, c, d, g, threshold)

            # STORE FIT DATA
            fit_df['Groups'] += [cat]*num_points
            fit_df['Individual'] += [individual]*num_points
            fit_df['Dilution'] += xrange
            fit_df['Absorbance'] += yfit_vals
            
            # STORE ENDPOINT TITER DATA
            if endpoint_titer:
                et_reciprocal = 1 / endpoint_titer
                et_df['Endpoint titer'].append(et_reciprocal)
                et_df['Groups'].append(cat)
            else:
                et_reciprocal = None
                
            output.write('\t'.join([cat, individual, str(endpoint_titer), str(et_reciprocal), str(r_squared)] + [str(x) for x in popt]) + '\n')
            
    output.close()
    
    return fit_df, et_df
    
    
def calc_rsquared(x_vals, y_vals, popt, regression_type):
    """
    Calculate the R-squared for each curve fit to assess goodnes of fit.
    
    Returns:
        r-squared (float)
    """

    # SUM OF SQUARED RESIDUALS
    if regression_type == '4PL':
        ss_res = sum([(fit_4PL(x_val, *popt) - y_vals.iloc[i])**2 for i, x_val in enumerate(x_vals)])
    else:
        ss_res = sum([(fit_5PL(x_val, *popt) - y_vals.iloc[i])**2 for i, x_val in enumerate(x_vals)])
    
    # TOTAL SUM OF SQUARES
    ss_tot = sum([(val - np.mean(y_vals))**2 for val in y_vals])
    
    r_squared = 1 - (ss_res / ss_tot)
    
    return r_squared
    
    
def fit_5PL(x, a, b, c, d, g):
    """
    Objective function for fitting a 5PL regression model.
    
    Returns:
        fx (float) - predicted absorbance value (y) for a given dilution (x).
    """
    
    fx = (a-d) / (1+((x/c)**b))**g + d
    return fx
    
    
def calc_endpoint_titer_5PL(a, b, c, d, g, threshold):
    """
    Calculate endpoint titer for using the 5PL regression model.
    
    Returns:
        endpoint_titer (float) - dilution at which the threshold crosses 
            the sigmoidal 5PL regression line (fitted curve).
    """
    
    endpoint_titer = c * (( (a-d) / (threshold-d) )**(1/g) - 1)**(1/b)
    return endpoint_titer
    

def fit_4PL(x, a, b, c, d):
    """
    Objective function for fitting a 4PL regression model.
    
    Returns:
        fx (float) - predicted absorbance value (y) for a given dilution (x).
    """
    
    fx = (a-d) / (1+((x/c)**b)) + d
    return fx
    
    
def calc_endpoint_titer_4PL(a, b, c, d, threshold):
    """
    Calculate endpoint titer for using the 4PL regression model.
    
    Returns:
        endpoint_titer (float) - dilution at which the threshold crosses 
            the sigmoidal 4PL regression line (fitted curve).
    """
    
    endpoint_titer = c * (( (a-d) / (threshold-d) ) - 1)**(1/b)
    return endpoint_titer
    
    
def prep_containers(df):
    """
    Prepare data containers for use in downstream functions.
    
    Returns:
        fit_df (dict) - Long-form dictionary that will contain 1000 fit data points per 
            sample ID representing 4PL or 5PL curvefit values.
        et_df (dict) - Long-form dictionary that will contain endpoint titer values
            for each group (for boxplot).
        categories (list) - List of all groups/categories in the user-provided data,
            with the original order preserved.
        individuals (list) - List of all individual sample IDs in the user-provided data,
            with the original order preserved.
    """
    
    fit_df = {'Groups':[],
            'Dilution':[],
            'Individual':[],
            'Absorbance':[]}
    et_df = {'Endpoint titer':[],
            'Groups':[]}

    # STORE LIST OF UNIQUE GROUPS/CATEGORIES WITH ORDER PRESERVED
    categories = []
    for cat in df['Groups']:
        if cat not in categories:
            categories.append(cat)

    # STORE LIST OF UNIQUE SAMPLE_IDs ("individual") WITH ORDER PRESERVED
    individuals = []
    for ind in df['Individual']:
        if ind not in individuals:
            individuals.append(ind)
    
    return fit_df, et_df, categories, individuals
    
    
def prep_output_file(regression_type):
    """
    Prepare output file to write endpoint titer data. Opens the file in write mode and 
    writes a header line specific to 4PL or 5PL regression.
    
    Returns:
        output (file handle) - handle of the output file to write endpoint titer data.
    """
    
    if regression_type == '4PL':
        header = '\t'.join(['Group', 'Sample ID', 'Endpoint Titer', '1 / Endpoint_Titer', 'Goodness of Fit (R-squared)'] + list('abcd'))
    else:
        header = '\t'.join(['Group', 'Sample ID', 'Endpoint Titer', '1 / Endpoint_Titer', 'Goodness of Fit (R-squared)'] + list('abcdg'))
        
    output = open(f'Endpoint_Titer_Results_{regression_type}.tsv', 'w')
    output.write(header + '\n')
    
    return output
    
    
def make_boxplot_endpoint_titers(df, regression_type):
    """
    Make boxplot with overlaid stripplot for the calculated endpoint titer values.
    
    Returns:
        None
    """
    
    sns.boxplot(x='Groups', y='Endpoint titer', data=df, hue='Groups', showfliers=False, palette=['0.7']*len(set(df['Groups'])))
    sns.stripplot(x='Groups', y='Endpoint titer', data=df, hue='Groups')
    plt.savefig(f'Endpoint_titer_Boxplots_{regression_type}.tif', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def make_lineplots(df, row_cats, plot_name):
    """
    Make figure with grid of simple line plots of Absorbance vs. Dilution for each unique sample ID.
    
    Returns:
        None
    """
    
    # MAKE COPY OF df AND CALCULATE LOG BASE 2 OF DILUTION
    df_copy = df.copy()
    df_copy['Dilution'] = np.log2( df_copy['Dilution'] )

    # SET STYLE TO BE SIMILAR TO ORIGINAL ELISA-R PLOTS
    sns.set_style('darkgrid', rc={'xtick.bottom': True, 'ytick.left': True, 'xtick.color': 'black', 'ytick.color': 'black'})
    
    # MAKE FACETGRID WITH SCATTERPLOTS AND LINEPLOTS
    grid = sns.FacetGrid(df_copy, col='Individual', row='Groups', hue='Groups')
    grid.map(sns.scatterplot, 'Dilution', 'Absorbance')
    grid.map(sns.lineplot, 'Dilution', 'Absorbance')

    # STYLIZE TITLES AND SUBPLOT SPACING
    grid.set_titles(col_template="{col_name}", row_template="{row_name}")
    grid.fig.subplots_adjust(hspace=0.2) # Adjust spacing

    # MAKE CUSTOM LEGEND ELEMENTS
    legend_elements = [Line2D([0], [0], marker='o', linestyle='None', markeredgecolor='None', label=cat, markerfacecolor=colors[i], markersize=10) for i, cat in enumerate(row_cats)]
    fig = plt.gcf()
    fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1,0.5), handletextpad=0.0, title='Groups')

    fig.set_size_inches(12, 8)
    plt.savefig(plot_name, bbox_inches='tight', dpi=600)
    plt.close()
    
    
def make_lineplots_fitdata(df, plot_name, threshold, lines_df=None):
    """
    Make figure with grid of curve-fit regression lines for each unique sample ID.
    
    Returns:
        None
    """
    
    # MAKE COPY OF df AND lines_df AND CALCULATE LOG BASE 2 OF DILUTION FOR EACH
    df_copy = df.copy()
    df_copy['Dilution'] = np.log2( df_copy['Dilution'] )
    
    lines_df_copy = lines_df.copy()
    lines_df_copy['Dilution'] = np.log2( lines_df_copy['Dilution'] )

    # SET STYLE TO BE SIMILAR TO ORIGINAL ELISA-R PLOTS
    sns.set_style('darkgrid', rc={'xtick.bottom': True, 'ytick.left': True, 'xtick.color': 'black', 'ytick.color': 'black'})
    
    # MAKE FACETGRID WITH SCATTERPLOTS AND LINEPLOTS
    grid = sns.FacetGrid(df_copy, col='Individual', row='Groups', hue='Groups')
    grid.map(sns.scatterplot, 'Dilution', 'Absorbance', data=df_copy)

    # LOOP OVER GRID TO OVERLAY UNIQUE CURVE FIT SIGMOID FOR EACH SUBPLOT
    row_cats = []
    col_cats = []
    for i, row_axes in enumerate(grid.axes):
        for j, ax in enumerate(row_axes):
            row_category = grid.row_names[i]
            col_category = grid.col_names[j]
            if row_category not in row_cats:
                row_cats.append(row_category)

            # EXTRACT AND PLOT CURVE-FIT VALUES SPECIFIC TO A SINGLE CONDITION+SAMPLE COMBINATION
            single_df = lines_df_copy[ (lines_df_copy['Groups'] == row_category) & (lines_df_copy['Individual'] == col_category)]
            sns.lineplot(x='Dilution', y='Absorbance', data=single_df, color='grey', ax=ax)
            sns.lineplot(x=(min(single_df['Dilution']), max(single_df['Dilution'])), y=(threshold, threshold), linestyle='--', color='0.5', ax=ax)
            ax.get_legend().remove()

    # STYLIZE TITLES AND SUBPLOT SPACING
    grid.set_titles(col_template="{col_name}", row_template="{row_name}")
    grid.fig.subplots_adjust(hspace=0.2) # Adjust spacing

    # MAKE CUSTOM LEGEND ELEMENTS
    legend_elements = [Line2D([0], [0], marker='o', linestyle='None', markeredgecolor='None', label=cat, markerfacecolor=colors[i], markersize=10) for i, cat in enumerate(row_cats)]
    fig = plt.gcf()
    fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1,0.5), handletextpad=0.0, title='Groups')

    fig.set_size_inches(12, 8)
    plt.savefig(plot_name, bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_data(data_file, file_type):
    """
    Read data from user-provided file and store as a dictionary.
    
    Returns:
        df (dict) - pandas DataFrame with the same structure as user-provided data.
    """
    
    if file_type == 'csv':
        df = pd.read_csv(data_file)
    else:
        df = pd.read_csv(data_file, sep='\t')

    df['Dilution'] = pd.eval(df['Dilution']).astype(float)

    return df
    
    
def validate_input(data_file, file_type, regression_type, threshold):
    """
    Run validation checks on command-line input from user.
    Checks that:
        1) the data_file could be successfully opened.
        2) the file_type is "csv" or "tsv" (comma-separated values or tab-separated values, respectively).
        3) the regression_type is "4PL" or "5PL" (4-parameter logistic or 5-parameter logistic, respectively).
        4) the threshold for calculating endpoint titer is >= 0.0.
        
    Returns:
        None
    """
    
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
        error_message += '\nERROR: Threshold value for calculating endpoint titer should be >= 0.'
        
    if error_message:
        print(error_message)
        print('\nExiting prematurely...')
        exit()


def get_args(arguments):
    """
    Gather and type-check command-line arguments.
    
    Returns:
        args (argparse namespace) - command line arguments used as variables.
    """
    
    parser = argparse.ArgumentParser(description='Calculate endpoint titer (ET) from ELISA data', prog='PyELISA')
    
    parser.add_argument('data_file', help="""Your data file.""")
    parser.add_argument('-f', '--file_type', type=str, default='csv', help="""The file type of the data_file provided. Options are "csv" or "tsv" Default=csv""")
    parser.add_argument('-r', '--regression_type', type=str, default='4PL', help="""The regression type used for curve fitting. Options are "4PL" or "5PL" (4-parameter logistic or 5-parameter logistic, respectively). Default=4PL""")
    parser.add_argument('-t', '--threshold', type=float, default=0.2, help="""The threshold to use for endpoint titer calculation. Default=0.2""")
    parser.add_argument('-x', '--maxfev', type=int, default=10000, help="""Maximum number of function evaluations during scipy.optimize curve fitting. Default=10000""")

    args = parser.parse_args(arguments)
    
    return args


if __name__ == '__main__':
    import sys, argparse
    args = get_args(sys.argv[1:])
    main(args)