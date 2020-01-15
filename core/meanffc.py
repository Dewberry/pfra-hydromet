import os
import logging
import warnings
import pathlib as pl
import scrapbook as sb
import papermill as pm
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from IPython.core.display import display
import matplotlib as mpl
import plotly.graph_objs as go
from plotly.offline import iplot
from matplotlib import pyplot as plt
logging.basicConfig(level=logging.ERROR)
plib = "pathlib.Path"


# --------------------------------------------------------------------------#

"""Functions called by PM_Sampler_Ops.ipynb, SSP_to_Mean_Curve.ipynb, 
   Stratified_Sampler.ipynb, and Make_Production_Run_List.ipynb. These 
   notebooks calculate the mean flow frequency curve, specify a range of 
   annual exceedance probabilities with corresponding weights, and assign
   those weights to discharge events.
"""

# --------------------------------------------------------------------------#


def make_directories(dir_lst: list, verbose: bool = True) -> None:
    """Check if each directory within the passed list exists and create any
       missing directories.
    """
    for directory in dir_lst:
        if not os.path.isdir(directory):
            os.makedirs(directory)
            if verbose:
                print('{0} - created'.format(str(directory)))
        else:
            if verbose:
                print('{0} - already exists'.format(str(directory)))
    return None


def list_ssp_files(path: plib, max_cl: float = 0.999,
                   verbose: bool = True) -> list:
    """Identify all .rpt files whose upper confidence limit is greater than
       0.5 and less than or equal to the maximum defined upper confidence
       limit.
    """
    ssp_results = []
    for file in pl.Path(path).glob('**/*.rpt'):
        filename = file.name
        assert '_' in filename, 'Filename does not include an "_" separating' \
                                ' the gage ID from the confidence limit, ' \
                                'should be "GageID_UpperConfidenceLimit"'
        split_stem = file.stem.split('_')
        assert len(split_stem) == 2, 'Filename contains more than two ' \
                                     'elements, should be ' \
                                     '"GageID_UpperConfidenceLimit"'
        cl = float(split_stem[1])
        if cl > 100.0:
            cl = cl/10.0
        if 50.0 < cl <= (max_cl*100):
            ssp_results.append(file)
            if verbose:
                print('{0} added to list'.format(filename))
        elif cl > (max_cl*100):
            if verbose:
                print('{0} not added to list, above maximum confidence'
                      ' limit'.format(filename))
        elif cl <= 50:
            if verbose:
                print('{0} not added to list, below minimum confidence '
                      'limit'.format(filename))
        else:
            if verbose:
                print('{0} not added to list, check naming '
                      'convention'.format(filename))
    assert len(ssp_results) > 0, 'No .rpt files identified in {0}'.format(path)
    return ssp_results


def make_ssp_table(ssp_results: list, version: str = '2_2') -> pd.DataFrame:
    """Create a table summarizing the SSP results where the index is the
       annual exceedance probability, each column is a confidence limit,
       and each cell is the corresponding discharge.
    """
    df = pd.DataFrame()
    assert_len = 'The .rpt files do not have the same number of annual ' \
                 'exceedance probabilities'
    for i, file in enumerate(ssp_results):
        if i == 0:
            df = GetFreqCurves(file, version=version)
        else:
            tmp = GetFreqCurves(file, version=version)
            cols = list(set(tmp.columns) - {'0.5', 'AEP'})
            tmp = tmp[cols].copy()
            assert df.shape[0] == tmp.shape[0], assert_len
            df = df.merge(tmp, left_index=True, right_index=True)
    df = df.set_index('AEP')
    df = df.reindex(sorted(df.columns), axis=1)
    cl_totals = [np.sum(df.iloc[:, i]) for i in range(df.shape[1])]
    assert np.all(np.diff(cl_totals) > 0), 'Q not increasing with CL as ' \
                                           'expected, check data'
    col_lst = [float(col) for col in df.columns]
    assert len(col_lst) == len(set(col_lst)), 'Duplicate columns'
    return df


def GetFreqCurves(f: plib, version: str = '2_2') -> pd.DataFrame:
    """Read the passed .rpt file and extract the annual exceedance
       probability, median flow frequency curve, and the upper and lower
       confidence limits. Note that the user must specify the HEC-SSP version
       used to create the .rpt files.
    """
    assert version in ['2_1', '2_2'], 'GetFreqCurve can only read .rpt files' \
                                      ' from versions 2.1 and 2.2 of HEC-SSP'
    if version == '2_1':
        line_breaks = [(1, 13), (14, 26), (27, 40), (41, 53), (54, 65)]
    else:
        line_breaks = [(1, 15), (16, 29), (34, 42), (46, 60), (61, 74)]
    read_results = False
    aep = []
    with open(f) as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if 'Upper Confidence Level:' in line:
                high = float(line.split(':')[1].replace('\n', ''))
            if 'Lower Confidence Level:' in line:
                low = float(line.split(':')[1].replace('\n', ''))
            if 'Frequency:' in line:
                aep.append(float(line.split(':')[1].replace('\n', '')) / 100.0)
            if 'Final Results' in line:
                read_results = True
            elif '<< Frequency Curve >>' in line and read_results:
                skiprows = i + 7
    assert (float(high) + float(low)) == 1.0, 'In {0} the upper and lower ' \
                                              'confidence limit values do not' \
                                              ' add to 1.0, check the user ' \
                                              'defined confidence limits in ' \
                                              'HEC-SSP'.format(f)
    cols = ['0.5', 'Variance', 'AEP', str(low), str(high)]
    df = pd.read_fwf(f, skiprows=skiprows, colspecs=line_breaks, names=cols)
    df = df[['0.5', str(low), str(high)]][0:len(aep)].copy()
    for col in df.columns:
        df[col] = df[col].apply(lambda x: float(x.replace(',', '')))
    df['AEP'] = aep
    return df


def monotonic_test(df: pd.DataFrame, adj_amount: float = 1.0) -> pd.DataFrame:
    """Test that the discharge increases with decreasing annual exceedance
       probability and adjust the discharge for the smallest annual exceedance
       probabilities if not.
    """
    no_adjust = True
    for col in df.columns:
        if np.diff(df[col]).max() >= 0:
            no_adjust = False
            maximum = df[col].max()
            idx = df[col].idxmax()
            diff = round(maximum - df.iloc[0][col], 1)
            adj_df_idx = df.loc[:idx][col].index
            num = len(adj_df_idx)
            val = np.arange(maximum, maximum + num * adj_amount, adj_amount)
            for i, v in enumerate(adj_df_idx):
                df.loc[v][col] = val[num - 1 - i]
            cl = float(col) * 100.0
            aep = df.iloc[0].name
            warnings.warn('Q not increasing with decreasing AEP for the {0}% '
                          'CL: difference of {1} cfs between {2} and {3}. '
                          'Adjusting Q'.format(cl, diff, aep, idx))
    if no_adjust:
        print('Discharge increases with decreasing annual exceedance '
              'probability for all confidence limits')
    return df


def zvar(cl: list) -> np.ndarray:
    """Used to calculate the standard normal z variate of the passed
       confidence limits or annual exceedance probabilities.
    """
    clz = np.array([norm.ppf((1 - clim)) for clim in cl])
    return clz


def binQ(df: pd.DataFrame) -> np.ndarray:
    """Determines the minimum and maximum discharge value for the passed
       dataframe and constructs an array of equally spaced discharge between
       these two values.
    """
    qmin = df.min().min()
    qmax = df.max().max()
    q = np.linspace(qmin, qmax, num=len(df))
    return q


def interp_AEP(df: pd.DataFrame, q: np.ndarray, clz: np.ndarray,
               aepz: np.ndarray, extrapolate: bool = True) -> pd.DataFrame:
    """Apply linear interpolation/extrapolation to calculate AEP(z) for
       each CL(z) and binned flow.
    """
    df1 = df.copy()
    df1['Q'] = q
    df1.set_index('Q', inplace=True)
    df1.columns = clz
    if not extrapolate:
        for cl in np.arange(len(clz)):
            q_min = df.iloc[:, cl].min()
            q_max = df.iloc[:, cl].max()
            f = interp1d(df.iloc[:, cl], aepz)
            for i in np.arange(len(q)):
                q_val = q[i]
                if q_val < q_min:
                    df1.iloc[i, cl] = aepz[-1]
                if (q_val >= q_min) & (q_val <= q_max):
                    df1.iloc[i, cl] = f(q_val)
                if q_val > q_max:
                    df1.iloc[i, cl] = aepz[0]
    if extrapolate:
        for cl in np.arange(len(clz)):
            f = interp1d(df.iloc[:, cl], aepz, fill_value='extrapolate')
            for i in np.arange(len(q)):
                q_val = q[i]
                df1.iloc[i, cl] = f(q_val)
    return df1


def zvar_inv(df: pd.DataFrame, cl: list) -> pd.DataFrame:
    """Calculate the the inverse of the standard normal Z variate for each
       annual exceedance probability and confidence limit.
    """
    df.columns = cl
    for clim in cl:
        df[clim] = 1 - norm.cdf(df[clim])
    return df


def mean_AEP(df: pd.DataFrame, exclude_tails: bool = True) -> list:
    """Calculate the mean (expected) value of the annual exceedance
       probability for each flow; the mean is equal to the area under the
       CDF, which is calculated using the trapezoidal rule. If exclude_tails
       is True than the integration only includes the area between the lower
       and upper confidence limits, else the integration includes the entire
       distribution, i.e. from 0.0 to the lower confidence limit and from the
       upper confidence limit to 1.0.
    """
    aepm_lst = []
    cl_arr = df.columns.values
    for logq in df.index:
        aep_arr = df.loc[logq].values
        aepm = trapz(aep_arr, x=cl_arr)
        if exclude_tails:
            scale_val = cl_arr[-1] - cl_arr[0]
            aepm = aepm/scale_val
        else:
            leftside = aep_arr[0] * cl_arr[0]
            rightside = aep_arr[-1] * cl_arr[0]
            aepm = aepm + leftside + rightside
        aepm_lst.append(aepm)
    return aepm_lst


def ffc_summary(standard_aep: list, standard_aepz: np.ndarray,
                add_ri: bool = False, verbose: bool = False) -> pd.DataFrame:
    """Initialize a summary table to store the mean and median flow frequency
       curves.
    """
    df = pd.DataFrame(data={'AEPz': standard_aepz}, index=standard_aep)
    df.index.name = 'AEP'
    if add_ri:
        df['RI'] = 1.0/df.index
    if verbose:
        display(df.head(2))
    return df


def interp_q(aepz: np.ndarray, q: np.ndarray):
    """Create a function for calculating the log flow given the normal z
       variate of the annual exceedance probability.
    """
    f = interp1d(aepz, q, fill_value='extrapolate')
    return f


def format_mean_curve(dic: dict, verbose: bool = True) -> pd.DataFrame:
    """Convert the mean flow frequency curve dictionary into a pandas
       dataframe.
    """
    df = pd.DataFrame.from_dict(dic)
    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'AEP'})
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.sort_values('AEP', ascending=False)
    df = df.set_index('AEP')
    if verbose:
        display(df.head(2))
    return df


def RI(min_val: int, max_val: int, nbin: int,
       verbose: bool = True) -> pd.DataFrame:
    """Divide the log of the range between the minimum and maximum
       recurrence interval into equal increments.
    """
    steps = 10 ** np.linspace(np.log10(min_val), np.log10(max_val), nbin + 1)
    ri = 1 / ((1 / steps[0:nbin] + 1 / steps[1:nbin + 1]) / 2.0)
    df = pd.DataFrame(data={'Center': ri})
    if verbose:
        display(df.head(2))
    return df


def AEP_weights(df: pd.DataFrame, min_val: float, max_val: float,
                verbose: bool = True) -> pd.DataFrame:
    """Given the center of the bin, calculate the weight of each binned
       interval, by taking the difference between the floor and the ceiling
       of the bin (left- and right-hand sides, respectively). The weight
       corresponds to the proportion of the annual exceedance probability
       within the binned interval. Note that the extra weight (events greater
       than the max_val) is added to the first bin.
    """
    n = df.shape[0]
    df['Ceiling'] = np.zeros(n)
    df['Floor'] = np.zeros(n)
    df['Ceiling'].iloc[0] = min_val
    for i in np.arange(n - 1):
        df.iloc[i]['Floor'] = (df.iloc[i]['Center'] +
                               df.iloc[i + 1]['Center'])/2
    df['Floor'].iloc[-1] = max_val
    for i in np.arange(1, n):
        df.iloc[i]['Ceiling'] = df.iloc[i - 1]['Floor']
    df1 = df.rdiv(1)
    df1['Weight'] = df1['Ceiling'] - df1['Floor']
    total = sum(df1['Weight'])
    df1['Weight'].iloc[0] += (0.5 - total)
    df1 = df1.set_index('Center')
    df1.index.name = 'AEP'
    if verbose:
        display(df1.head(2))
    return df1


def add_events(df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
    """Add an events column to the passed dataframe and set it as the index.
    """
    e = ['E{0}'.format(str(i).zfill(4)) for i in np.arange(1, df.shape[0] + 1)]
    df['Events'] = e
    df.reset_index(inplace=True)
    df = df.set_index('Events')
    df.index.name = ''
    if verbose:
        display(df.head(2))
    return df


# --------------------------------------------------------------------------#

"""Plotting functions called by SSP_to_Mean_Curve.ipynb. This notebook 
   calculates the mean flow frequency curve using Bulletin 17C confidence 
   limits calculated in HEC-SSP.
"""

# --------------------------------------------------------------------------#


def plot_ssp_transform(data: pd.DataFrame, aepz: np.ndarray,
                       clz: np.ndarray) -> None:
    """Plot the log10 discharge verses the standard normal z variate of the
       annual exceedance probability (note that each line represents a
       single confidence limit). Also plot the log10 discharge verse the
       standard normal z variate of the confidence limits (note that each
       line represents a single annual exceedance probability).
    """
    cmap = mpl.cm.viridis
    cmap1 = mpl.cm.viridis_r
    plt.figure(1, figsize=(24, 6))
    plt.clf()
    ax1 = plt.subplot2grid((1, 2), (0, 0),
                           xlabel='Annual Exceedance Probability, [z variate]',
                           ylabel='Discharge, [log(cfs)]',
                           title='Discharge vs. Annual Exceedance Probability'
                                 ' (Each Line is a Confidence Limit)')
    for i in np.arange(len(clz)):
        ax1.plot(aepz, data.iloc[:, i], linestyle="-", marker='.',
                 color=cmap(i / len(clz)))
    ax1 = plt.subplot2grid((1, 2), (0, 1),
                           xlabel='Confidence Limits, [z variate]',
                           ylabel='Discharge, [log(cfs)]',
                           title='Discharge vs. Confidence Limits (Each Line'
                                 ' is an Annual Exceedance Probability)')
    for i in np.arange(len(aepz)):
        ax1.plot(clz, data.iloc[i], linestyle="-", marker='.',
                 color=cmap1(i / len(aepz)))
    return None


def plot_ssp_interp(df: pd.DataFrame) -> None:
    """Plot the standard normal z variates of the annual exceedance probability
       verses the confidence limits (note that each line represents a single
       discharge).
    """
    cmap = mpl.cm.viridis
    plt.figure(2, figsize=(12, 6))
    plt.clf()
    ax1 = plt.subplot2grid((1, 1), (0, 0),
                           xlabel='Confidence Limits, [z variate]',
                           ylabel='Annual Exceedance Probability, [z variate]',
                           title='Annual Exceedance Probability vs. Confidence'
                                 ' Limits (Each Line is a Discharge)')
    clz = list(df.columns)
    for i, idx in enumerate(df.index):
        ax1.plot(clz, df.loc[idx].values, linestyle="-", marker='.',
                 color=cmap(i / df.shape[0]))
    return None


def plot_ssp_interptrans(df: pd.DataFrame) -> None:
    """Plot the annual exceedance probability verses the confidence limits
       (note that each line represents a single discharge).
    """
    cmap = mpl.cm.viridis
    plt.figure(3, figsize=(12, 6))
    plt.clf()
    ax1 = plt.subplot2grid((1, 1), (0, 0),
                           xlabel='Confidence Limits',
                           ylabel='Annual Exceedance Probability',
                           title='Annual Exceedance Probability vs. Confidence'
                                 ' Limits (Each Line is a Discharge)')
    cl = list(df.columns)
    for i, idx in enumerate(df.index):
        ax1.plot(cl, df.loc[idx].values, linestyle="-", marker='.',
                 color=cmap(i / df.shape[0]))
    return None


def plot_ssp_meanmed(aepz: np.ndarray, df: pd.DataFrame, aepmz: np.ndarray,
                     q: np.ndarray) -> None:
    """Plot the mean and median flow frequency curves using log10 discharge
       and the standard normal z variate of the annual exceedance probability.
    """
    plt.figure(4, figsize=(12, 6))
    plt.clf()
    ax1 = plt.subplot2grid((1, 1), (0, 0),
                           xlabel='Annual Exceedance Probability, [z variate]',
                           ylabel='Discharge, [log(cfs)]',
                           title='Mean and Median Flow Frequency Curve (Log10 '
                                 'Discharge, Z Variate AEP)')
    ax1.plot(aepz, df.iloc[:]['0.5'], linestyle='-', marker='.',
             label='Median', color='black')
    ax1.plot(aepmz, q, linestyle='-', marker='.', label='Mean', color='red')
    ax1.legend(loc='lower right', frameon=False)
    return None


def plot_ssp_meanmedffc(table: pd.DataFrame, gage_id: str) -> None:
    """Plot the mean and median flow frequency curves using both the annual
       exceedance probability and the recurrence interval on a log-log plot.
    """
    perc_aep = [aep*100 for aep in table.index]
    ri = [1.0/aep for aep in table.index]
    xlabels = ['Annual Exceedance Probability', 'Recurrence Interval']
    xunits = ['%', 'years']
    position = ['left', 'right']
    fig, ax = plt.subplots(1, 2, figsize=(24, 6))
    for i, idx in enumerate([perc_aep, ri]):
        ax[i].plot(idx, table['Q_Median_cfs'], linestyle='-',
                   marker='.', label='Median', color='black')
        ax[i].plot(idx, table['Q_Mean_cfs'], linestyle='-',
                   marker='.', label='Mean', color='red')
        ax[i].set_xscale('log')
        ax[i].set_yscale('log')
        ax[i].set_xlabel('{0}, [{1}]'.format(xlabels[i], xunits[i]))
        ax[i].set_ylabel('Discharge, [cfs]')
        ax[i].set_title('Discharge vs. {0} (Station ID: {1})'
                        ''.format(xlabels[i], gage_id))
        ax[i].grid(True, which='both')
        ax[i].legend(loc='lower {0}'.format(position[i]), frameon=True)
        if i == 0:
            format_xtick = ['{:}'.format(x) for x in ax[i].get_xticks()]
            ax[i].set_xticklabels(format_xtick)
    return None


def plotly_ssp_meanffc(table: pd.DataFrame, gage_id: str) -> None:
    """Interactively plot the mean flow frequency curve on a log-log plot.
    """
    perc_aep = [aep*100 for aep in table.index]
    ri = [round(1/aep) for aep in table.index]
    events = np.arange(1, len(ri))
    q = table['Q_Mean_cfs']
    trace = go.Scatter(x=perc_aep, y=q, name='RI', mode='markers',
                       text=['{0:,.0f}-year'.format(x) for x in ri],
                       hoverinfo='text+name')
    trace2 = go.Scatter(x=perc_aep, y=q, name='Q', mode='lines',
                        text=table['Q_Mean_cfs'
                                   ''].apply(lambda x: '{0:,.2f}'.format(x)),
                        hoverinfo='text+name',
                        line=dict(color='rgb(204, 0, 153)'))
    trace3 = go.Scatter(x=perc_aep, y=q, name='Event', mode='markers',
                        text=events, hoverinfo='text+name')
    data = [trace, trace2, trace3]
    layout = go.Layout(dict(title='Mean Discharge vs. Annual Exceedance '
                                  'Probability (Station ID: {0})'
                                  ''.format(gage_id),
                            xaxis=dict(title='Annual Exceedance Probability, '
                                             '[%]', type='log',
                                       autorange=True, tickmode='linear'),
                            yaxis=dict(title='Discharge, [cfs]', type='log',
                                       autorange=True),
                            legend=dict(orientation="h"),
                            font=dict(color='rgb(0,0,0)'),
                            paper_bgcolor='rgb(255,255,255)',
                            plot_bgcolor='rgb(255,255,255)', showlegend=False))
    fig = go.Figure(data=data, layout=layout)
    iplot(fig)
    return None
