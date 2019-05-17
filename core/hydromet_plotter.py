import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter

geoDF = 'GeoDataFrame'


#---------------------------------------------------------------------------#

'''Plotting functions called by EventTable.ipynb. This 
   notebook calculates excess rainfall by first randomly selecting a 
   precipitation recurrance interval and corresponding precipitation amount, 
   precipitation temporal distribution, and curve number for the area of 
   interest. The randomly selected precipitation data and curve number are 
   then used by the curve number approach to calculate the excess rainfall 
   amount for the corresponding recurrance interval. The procedure is 
   repeated for the specified number of events/recurrance intervals. 
'''

#---------------------------------------------------------------------------#


def plot_area_of_interest(geo_df: geoDF, select_data: str, 
                                                column: str) -> plt.subplots:
    '''Plots the column of the geodataframe with matplotlib.
    '''
    fig = geo_df.plot(column = column, categorical = True, figsize = (14, 18))
    fig.set_title('Area of Interest (ID: {})'.format(select_data))
    fig.grid()


def plot_rand_precip_data(df: pd.DataFrame, rand_data: list, duration: int,
                                                xlim: int=2) -> plt.subplots:
    '''Plots the precipitation amount selected from the lognormal 
       distribution with the return period on the x-axis and the amount of 
       precipitation on the y-axis, along with the expected values, lower 
       90% confidence limits, and upper 90% confidence limits.
    '''
    fig, ax = plt.subplots(figsize=(20, 6))
    ax.grid(True, which="both")
    ax.semilogx(df.index, df['Upper (90%)'], color='darkolivegreen', 
                        linewidth=3, label=r'Upper (90%) Confidence Limit')
    ax.semilogx(df.index, df['Expected Value'], color='darkred', linewidth=2, 
                                                    label='Expected Value')
    ax.semilogx(df.index, df['Lower (90%)'], color='darkblue', linewidth=3, 
                                    label=r'Lower (90%) Confidence Limit')
    for col in rand_data:
        ax.scatter(df.index, df[col], s=25, edgecolor='black', 
                    linewidth='1',  facecolor=np.random.rand(4,), label=col)
    def mil(x: float, pos: int) -> str:
        ''' Convert the passed x-value to a string.
        '''
        return '{}'.format(x)
    mil_formatter = FuncFormatter(mil)
    for axis in [ax.xaxis]:
        axis.set_major_formatter(mil_formatter)
    ax.set_xlabel('Return Period (years)', fontsize=18)
    ax.set_ylabel('P (inches)', fontsize=18)
    ax.legend()
    ax.set_xlim(xlim,)
    ax.set_ylim(0,)
    ax.set_title('Random Precipitation  \n{} Hour Duration'.format(duration),
                                                                fontsize=18)


def plot_deciles_by_quartile(curve_group: dict, qrank: list,
                qmap: dict, vol: int, reg: int, dur: int) -> plt.subplots:
    '''Plots the temporal distribution at each decile for each quartile. 
    '''
    fig, ax = plt.subplots(2,2, figsize=(24,10))
    for axi in ax.flat:
        axi.xaxis.set_major_locator(plt.MultipleLocator((
                                            curve_group['q1'].shape[0]-1)/6))
        axi.xaxis.set_minor_locator(plt.MultipleLocator(1))
    axis_num=[[0,0], [0,1], [1,0], [1,1]]
    for i, val in enumerate(qmap['map'].keys()):
        for col in curve_group[val].columns:
            plt.suptitle('Volume '+str(vol)+' Region '+str(reg)+' Duration '+str(dur),
                                        fontsize = 20, x  = 0.507, y = 1.02)
            ax[axis_num[i][0],axis_num[i][1]].plot(curve_group[val][col], 
                                                                label=col) 
            ax[axis_num[i][0],axis_num[i][1]].grid()
            ax[axis_num[i][0],axis_num[i][1]].set_title('Quartile {0}\n{1}%'
                    ' of Cases'.format(i+1, int(qrank[i]*100)), fontsize=16)
            ax[axis_num[i][0],axis_num[i][1]].legend(title='Deciles')
            ax[axis_num[i][0],axis_num[i][1]].set_xlabel('Time (hours)', 
                                                                fontsize=14)
            ax[axis_num[i][0],axis_num[i][1]].set_ylabel('Precip (% Total)', 
                                                                fontsize=14)
    plt.tight_layout()


def plot_decile_histogram(df: pd.DataFrame) -> plt.subplots:
    '''Plots a histogram of the randomly selected decile numbers within the
       passed dataframe.
    '''
    fig = df.hist(bins=20, figsize=(20,6), grid=False)


def plot_rainfall_and_excess(final_precip: pd.DataFrame, 
    cum_excess: pd.DataFrame, dur: int=24, iplot: bool=False) -> plt.subplot:
    '''Plots the cumulative rainfall and runoff for each randomly selected 
       event.
    '''
    fig, ax = plt.subplots(1, 2, figsize=(24,5))
    for axi in ax.flat:
        axi.xaxis.set_major_locator(plt.MultipleLocator(dur/6))
        axi.xaxis.set_minor_locator(plt.MultipleLocator(1))
    for col in final_precip.columns:
        ax[0].plot(final_precip[col])
        ax[1].plot(cum_excess[col]) 
    nevents = final_precip.shape[1]
    ax[0].set_title('{} Cumulative Rainfall' 
                                    'Events'.format(nevents), fontsize=18)
    ax[1].set_title('{} Cumulative Runoff' 
                                    'Events'.format(nevents), fontsize=18)
    ax[0].set_ylim(0, 1.1*final_precip.max().max())
    ax[1].set_ylim(0, 1.1*final_precip.max().max())
    ax[0].set_xlabel('Time (hours)', fontsize=18)
    ax[1].set_xlabel('Time (hours)', fontsize=18)
    ax[0].set_ylabel('Precip (inches)', fontsize=18)
    ax[1].set_ylabel('Runoff (inches)', fontsize=18)
    ax[0].grid(True)
    ax[1].grid(True)
    if iplot:
        plt.close(fig)
    return fig


def plot_curve_groups(reordered_group: dict, reordered_curves: pd.DataFrame, curve_test_df: pd.DataFrame, y_max: float, final: bool=True) -> plt.subplots:
    '''Plots the mean temporal distribution and the corresponding 
       individual temporal distributions of a curve group for each group as 
       separate plots.  
    '''
    c_df = reordered_curves.copy()
    nc_df = curve_test_df.copy()
    x = list(c_df.index.values)
    if not final: x += [c_df.shape[0]]
    for c in c_df.columns:
        fig, ax = plt.subplots(figsize=(30,8))
        for nc in reordered_group[c]:
            y = nc_df[nc].values
            if not final: y = np.insert(y, 0, 0)
            ax.plot(x, y, alpha=0.75, label = nc);
        y = c_df[c].values
        if not final: y = np.insert(y, 0, 0)
        ax.plot(x, y, color='black', linewidth='2', label='Mean Curve')    
        ax.grid()
        ax.set_xlabel('Duration, [hours]')
        ax.set_ylabel('Runoff, [inches]')
        ax.set_ylim(0, y_max*1.1)
        ax.set_title('Group {} Temporal Distribution'.format(c))
        ax.legend()


def plot_grouped_curves(final_curves: dict, y_max: float, 
                                        iplot: bool=False) -> plt.subplots:
    '''Plots the mean curve of each group of curves determined using the 
       convolution test as well as the curves that were not grouped. 
    '''
    fig, ax = plt.subplots(figsize=(30,8))
    for col in final_curves.columns:
        ax.plot(final_curves[col]);
    ax.grid()
    ax.set_xlabel('Duration, [hours]')
    ax.set_ylabel('Runoff, [inches]')
    ax.set_ylim(0, y_max*1.1)
    ax.set_title('{} Temporal Curves'.format(final_curves.shape[1]))
    if iplot:
        plt.close(fig)
    return fig


#---------------------------------------------------------------------------#