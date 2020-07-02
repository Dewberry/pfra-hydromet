import pathlib as pl
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy import interpolate, stats, special
from matplotlib import pyplot as plt
from cycler import cycler


#----------------------------------------------------------------------------------------------------------------------#
# Functions called by EventTable_Stratified.ipynb.
#----------------------------------------------------------------------------------------------------------------------#

def Q_SCS(R: np.array, CN: float, mu: float) -> float:
    """SCS-CN runoff formula.
    """
    S = 1000.0/CN-10.0
    return (R-mu*S)**2/(R-mu*S+S)


def Norm_Constant_GEV(x: np.ndarray, PMP: float) -> float:
    """Constant for distribution truncation at the PMP value.
    """ 
    return 1.0/stats.genextreme.cdf(PMP, x[2], x[0], x[1])


def Norm_Constant_LN(SD: float, mu: float, PMP: float) -> float:
    """Constant for distribution truncation at the PMP value. 
    """ 
    return 1.0/stats.lognorm.cdf(PMP, SD, scale = np.exp(mu))


def PDF_GEV(R: np.ndarray, x: np.ndarray, PMP: float) -> np.ndarray:
    """
    """
    return Norm_Constant_GEV(x, PMP)*stats.genextreme.pdf(R, x[2], x[0], x[1])


def CDF_GEV(R: np.ndarray, x: np.ndarray, PMP: float) -> float:
    """
    """
    return Norm_Constant_GEV(x, PMP)*stats.genextreme.cdf(R, x[2], x[0], x[1])


def PPF_GEV(P: np.ndarray, x: np.ndarray, PMP: float) -> np.ndarray:
    """
    """
    return stats.genextreme.ppf(P/Norm_Constant_GEV(x, PMP), x[2], x[0], x[1])


def GEV_Parameters(df: pd.DataFrame, GEV_Parameters: np.ndarray, bounds: tuple, ID: str, PMP: float) -> pd.DataFrame:
    """Function defines an objective function for finding the GEV parameters and then determines the best GEV parameters 
       that minimize the difference between the GEV and comparison data.
    """    
    def objective_func_GEV(x: np.ndarray) -> float: 
        """Calculates the sum of the squared residuals between the return interval and return interval calculated from 
           the GEV CDF with the differences normalized by the return interval. 
        """ 
        return sum(np.square((RI-1/(1-CDF_GEV(row[ID], x, PMP)))/RI) for RI, row in df.iterrows())
    solution = minimize(objective_func_GEV, GEV_Parameters, method='SLSQP', bounds=bounds, options={'disp': False})
    df_GEV_parameters = pd.DataFrame(data=solution.x, index=["mu", "sigma", "xi"], columns=["GEV {}".format(ID)])
    return df_GEV_parameters


def GEV_parameters_Fit(raw_precip: pd.DataFrame, ID: str, PMP: float) -> pd.DataFrame:
    """This function provides initial value for finding the GEV parameters and then finds the best GEV parameters using 
       the function GEV_parameters.
    """
    year = raw_precip.index.values
    weights = np.append(1/year[:-1]-1/year[1:], 1/year[-1])
    Avg = (weights*raw_precip[ID]).sum()
    GEV_parameters = np.array([Avg*0.8, 0.5, -0.25])
    bounds = ((Avg*0.7, Avg*1.0), (0.01, 1.1), (-0.5, 0.0))
    df_GEV_parameters = GEV_Parameters(raw_precip, GEV_parameters, bounds, ID, PMP)
    return df_GEV_parameters


def Avg_R_integrand(R: float, GEV_parameters: np.ndarray, PMP: float) -> float:
    """This function defines the integrand for calculating an average based on the GEV distribution.
    """
    return R*PDF_GEV(R, GEV_parameters, PMP)


def Avg_R(lower_bound: float, upper_bound: float, GEV_parameters: np.ndarray, PMP: float) -> float:
    """Calculates the average value of the GEV distribution of rainfall or runoff based on an upper and lower bound.
    """
    return quad(Avg_R_integrand, lower_bound, upper_bound, args=(GEV_parameters, PMP))

def GEV_RI(RI: np.ndarray, GEV_parameters: np.ndarray, PMP: float) -> np.ndarray:
    """Provides rainfall or runoff as a function of the return interval (RI).
    """
    return PPF_GEV(1-1.0/RI, GEV_parameters, PMP)


def objective_func_bound_GEV(RI_lower: float, RI_upper: float, RI_middle: float, GEV_parameters: np.ndarray, 
                             PMP: float) -> float:
    """Calculates the square of the error between the average rainfall or runoff calculated from the bin floor and 
       ceiling (given in terms of RI) and the rainfall or runoff of the return period of interest.
    """ 
    return np.square(Avg_R(GEV_RI(RI_lower, GEV_parameters, PMP), GEV_RI(RI_upper, GEV_parameters, PMP), GEV_parameters, 
                                  PMP)[0]/(1.0/RI_lower-1.0/RI_upper) - GEV_RI(RI_middle, GEV_parameters, PMP))


def bound_lower_GEV(RI_upper: float, RI_middle: float, GEV_parameters: np.ndarray, initial_value: float,
                    PMP: float) -> float:
    """Finds the rainfall or runoff bin floor given the bin ceiling and average return interval of the bin.
    """
    return minimize(objective_func_bound_GEV, initial_value, args = (RI_upper, RI_middle, GEV_parameters, PMP),
                    method='SLSQP', bounds=[(1.0, RI_upper*0.999)], options={'disp': False})


def PDF_QlS(Q: np.ndarray, S: float, mu: float, GEV_parameters: np.ndarray, PMP: float) -> np.ndarray:
    """This function provides the runoff PDF conditional on max potential retention, where Q=runoff, S=max potential 
       retention, and mu=initial abstraction parameter.
    """
    return (Q+2.0*S+np.sqrt(Q)*np.sqrt(Q+4.0*S))/(2.0*np.sqrt(Q)*np.sqrt(Q+4.0*S))*\
            PDF_GEV(1.0/2.0*(Q+np.sqrt(Q)*np.sqrt(Q+4.0*S)+2.0*S*mu), GEV_parameters, PMP)


def PDF_Q(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray,
          Delta_P: float, error_PQ: float) -> float: 
    """
    """
    return sum(Delta_P*PDF_QlS(Q, S_avg_partition, mu, GEV_parameters, PMP) for S_avg_partition in 
               partition_avg)/(1-error_PQ)

def Qzero_integrand(S: float, mu: float, alpha: float, beta: float, S_limit: float, GEV_parameters: np.ndarray,
                    PMP: float, error_PQ: float) -> float:
    """Defines the integrand for calculating the probability of zero runoff.
    """
    return (CDF_GEV(S*mu, GEV_parameters, PMP)-CDF_GEV(0, GEV_parameters, PMP))\
           *(1.0/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit)/(1-error_PQ)


def P_Qzero(mu: float, alpha: float, beta: float, S_limit: float, GEV_parameters: np.ndarray, PMP: float, 
            error_PQ: float) -> float:
    """Defines discrete probability of zero runoff (integrated).
    """
    return quad(Qzero_integrand, 0, S_limit, args =(mu, alpha, beta, S_limit, GEV_parameters, PMP, error_PQ)) 


def CDF_Q(Q: float, mu: float, alpha: float, beta: float, S_limit: float, GEV_parameters: np.ndarray, PMP: float, 
          partition_avg: np.ndarray, Delta_P: float, error_PQ: float) -> float:
    """Defines the cumulative distribution function for runoff. PDF PDF_Q(u) is integrated from zero to an arbitrary 
       runoff Q.
    """
    return quad(PDF_Q, 0.0, Q, args=(mu, GEV_parameters, PMP, partition_avg, Delta_P, error_PQ))[0]\
                +P_Qzero(mu, alpha, beta, S_limit, GEV_parameters, PMP, error_PQ)[0]

def Avg_Q_integrand(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, 
                    Delta_P: float, error_PQ: float) -> float:
    """
    """
    return Q*PDF_Q(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, error_PQ)


def Avg_Q(lower_bound: float, upper_bound: float, mu: float, GEV_parameters: np.ndarray, PMP: float, 
          partition_avg: np.ndarray, Delta_P: float, error_PQ: float) -> float:
    """
    """
    return quad(Avg_Q_integrand, lower_bound, upper_bound, args=(mu, GEV_parameters, PMP, partition_avg, Delta_P, 
                error_PQ))


def runoff_RI(RI: float, f_RI_Q: interpolate.interp1d) -> float:
    """Defines runoff as a function of the return interval (RI).
    """
    return f_RI_Q(RI)


def objective_func_bound_runoff_L(RI_lower: float, RI_upper: float, RI_middle: float, mu: float, 
                                  GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, 
                                  f_RI_Q: interpolate.interp1d, error_PQ: float) -> float:
    """Calculates the square of the error between the average runoff calculated from the bin floor and ceiling (given in 
       terms of RI) and the runoff of the return period of interest.
    """ 
    return np.square(Avg_Q(runoff_RI(RI_lower, f_RI_Q), runoff_RI(RI_upper,  f_RI_Q), mu, GEV_parameters, PMP, 
                      partition_avg, Delta_P, error_PQ)[0]/(1.0/RI_lower-1.0/RI_upper)-runoff_RI(RI_middle, f_RI_Q)) 


def Bound_L(RI_upper: float, RI_middle: float, mu: float, GEV_parameters: np.ndarray, PMP: float, 
            partition_avg: np.ndarray, Delta_P: float, initial_value: float, f_RI_Q: interpolate.interp1d, 
            error_PQ: float) -> float:
    """Calculates runoff bin floor given the bin ceiling and average return interval of the bin.
    """
    return minimize(objective_func_bound_runoff_L, initial_value, 
                    args = (RI_upper, RI_middle, mu, GEV_parameters, PMP, partition_avg, Delta_P, f_RI_Q, error_PQ),
                    method='SLSQP', bounds=[(1.0, RI_upper)], options={'disp': False})


def PDF_S(S: np.ndarray, alpha: float, beta: float, S_limit: float)-> float:
    """Defines the distribution of the max potential retention.
    """
    return (1.0/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit)


def S_avg_integrand(S: float, alpha: float, beta: float, S_limit: float) -> float:
    """Defines the integrand for finding the average value over each partition.
    """
    return S*PDF_S(S, alpha, beta, S_limit)


def S_avg_partition(alpha: float, beta: float, S_limit: float, lower_bound: float, upper_bound: float) -> float:
    """Defines the integration for the average value over each partition.
    """
    return quad(S_avg_integrand, lower_bound, upper_bound, args=(alpha, beta, S_limit))


def PDF_SlQ(S: np.ndarray, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, 
            Delta_P: float, alpha: float, beta: float, S_limit: float, error_PQ: float) -> float:
    """Defines the PDF of the max potential retention, S, conditional on runoff, Q.
    """
    return PDF_QlS(Q, S, mu, GEV_parameters, PMP)*PDF_S(S, alpha, beta, S_limit)/PDF_Q(Q, mu, GEV_parameters, PMP, 
                                                                                       partition_avg, Delta_P, error_PQ)


def CDF_SlQ(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, 
            Delta_P: float, alpha: float, beta: float, S_limit: float, error_PQ: float) -> float:
    """Defines the CDF of the max potential retention, S, conditional on runoff Q.
    """
    return quad(PDF_SlQ, 0, S, args=(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, 
                                     error_PQ))[0]


def Avg_SlQ_integrand(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, 
                      Delta_P: float, alpha: float, beta: float, S_limit: float, error_PQ: float) -> float:
    """Defines the integrand for calculating the average max potential retention.
    """
    return S*PDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, error_PQ)


def Avg_SlQ(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, 
            alpha: float, beta: float, S_limit: float, error_PQ: float, lower_bound: float, 
            upper_bound: float) -> float:
    """Derives the average values of the max potential retention by integrating between an upper and lower bound.
    """
    return quad(Avg_SlQ_integrand, lower_bound, upper_bound, args=(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, 
                                                                   alpha, beta, S_limit, error_PQ))[0]


def objective_func_median_S(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, 
                            partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, 
                            S_limit: float, error_PQ: float) -> float:
    """Calculates the square of the error between the the CDF value and the median value of 0.5.
    """ 
    return np.square(CDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, error_PQ)-0.5)


def Median_S(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, 
             alpha: float, beta: float, S_limit: float, error_PQ: float, bounds: list, Initial_Value: float) -> float:
    """
    """
    return minimize(objective_func_median_S, Initial_Value, 
                    args = (Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, error_PQ),
                    method='SLSQP', bounds=bounds, options={'disp': False})


def partition_S_avgs(n_partition: int, Delta_P: float, alpha: float, beta: float, S_limit: float) -> np.ndarray:
    """Calculates the average value of the max potential retention for n partitions of the distribution.
    """
    Bounds = np.linspace(0.0, 1.0, n_partition+1)
    Bounds_S = S_limit*stats.beta(alpha, beta).ppf(Bounds)
    Bounds_Lower = Bounds_S[:-1]
    Bounds_Upper = Bounds_S[1:]
    partition_avg = np.array([S_avg_partition(alpha, beta, S_limit, lower, upper)[0] 
                                  for lower, upper in zip(Bounds_Lower, Bounds_Upper)])/Delta_P
    return partition_avg


def weights_Rainfall(Return_Intervals: np.ndarray, GEV_parameters: np.ndarray, PMP: float, RI_upper_bound: float, 
                     NOAA_precip: pd.DataFrame, ID: str, CN: float, mu: float) -> pd.DataFrame:
    """Calculate the weights of the rainfall events. If the RI of interest are already in the mean curve tablee, RI values for the rainfall 
       are taken directly from the input data (NOAA_precip or mean precip curve) instead of being calculated from the fitted GEV.
    """
    Size = Return_Intervals.size
    Bin_Bounds_R_topdown = np.zeros(Size+1)
    Bin_Bounds_R_topdown[Size] = RI_upper_bound
    for i in range(0, Size):
        Bin_Bounds_R_topdown[Size-i-1] = bound_lower_GEV(Bin_Bounds_R_topdown[Size-i], Return_Intervals[Size-i-1], 
                                                         GEV_parameters, Return_Intervals[Size-i-1]*0.2, PMP).x[0] 
    lower_bound = GEV_RI(RI_upper_bound, GEV_parameters, PMP)
    Avg_PlusR = Avg_R(lower_bound, PMP, GEV_parameters, PMP)[0]/(1.0/RI_upper_bound)
    Prob_Plus = CDF_GEV(Avg_PlusR, GEV_parameters, PMP)  
    RI_index = np.append(Return_Intervals, 1.0/(1-Prob_Plus)).astype(int)
    weights_R_topdown = (1.0/Bin_Bounds_R_topdown[:-1]-1.0/Bin_Bounds_R_topdown[1:]).astype(float)
    weights_R_topdown = np.append(weights_R_topdown, 1.0/RI_upper_bound)
    data = np.vstack((Bin_Bounds_R_topdown, np.append(Bin_Bounds_R_topdown[1:], np.inf), weights_R_topdown)).T
    df_weights = pd.DataFrame(data=data, index=RI_index, columns=['Bin Floor', 'Bin Ceiling', 'Event Weight'])
    RI_data =NOAA_precip[NOAA_precip.index.isin(Return_Intervals)].index.to_numpy().astype(int)
    RI_index_calc = RI_index[np.isin(RI_index, RI_data, invert=True)]
    Precip_calculate = GEV_RI(RI_index_calc, GEV_parameters, PMP)
    df2 = pd.DataFrame(data = Precip_calculate, index = RI_index_calc, columns=[ID])  
    df_R_NOAA_E = NOAA_precip[NOAA_precip.index.isin(RI_data)].copy()
    df_precip = df_R_NOAA_E.append(df2)
    df_precip = pd.DataFrame(df_precip.iloc[:,0])#df_precip.drop('P_Median_in', axis=1)
    Q = Q_SCS(df_precip[ID].values, CN, mu)
    df_precip['Runoff'] = Q 
    return pd.concat([df_weights, df_precip], axis=1)


def runoff_GEV(mu: float, GEV_parameters, PMP: float, alpha: float, beta: float, S_limit: float, 
               partition_avg: np.ndarray, Delta_P: float, error_PQ: float, n_partitions_Q: int=40) -> tuple:
    """Calculates the values of runoff versus return period, fits a GEV distribution to the results, and returns 
       dataframes for both the GEV parameters and runoff as a function of the return interval.
    """
    Q_line = np.linspace(0.01, PMP - 0.01, n_partitions_Q+1)
    Return_PeriodQ = 1.0/(1-np.transpose([CDF_Q(Q, mu, alpha, beta, S_limit, GEV_parameters, PMP, partition_avg, 
                                                Delta_P, error_PQ) for Q in Q_line]))
    df_runoff = pd.DataFrame(Q_line, index=Return_PeriodQ, columns=['Runoff'])
    df_GEV_parameters_R = GEV_parameters_Fit(df_runoff, 'Runoff', PMP)
    return df_runoff, df_GEV_parameters_R


def runoff_weights(Return_Intervals: np.ndarray, RI_upper_bound: float, mu: float, GEV_Parameters_Runoff: pd.DataFrame,
                   GEV_Parameters_Rain: pd.DataFrame, PMP: float, partition_avg: np.ndarray, Delta_P: float, 
                   error_PQ: float) -> pd.DataFrame:
    """Calculate the weights of the runoff events assuming that a GEV distribution is the best fit to the runoff 
       distribution that was derived analytically (and implemented with numerical integration) based on the GEV PDF of 
       rainfall and the distribution of the max potential retention.
    """               
    Size = Return_Intervals.size
    Bin_Bounds = np.zeros(Size+1)
    Bin_Bounds[Size] = RI_upper_bound
    for i in range(0, Size):
        Bin_Bounds[Size-i-1] = bound_lower_GEV(Bin_Bounds[Size-i], Return_Intervals[Size-i-1], GEV_Parameters_Runoff, 
                                               np.array([Return_Intervals[Size-i-1]*0.2]), PMP).x[0] 
    lower_bound = GEV_RI(RI_upper_bound, GEV_Parameters_Runoff, PMP)
    Avg_Plus = Avg_Q(lower_bound, PMP, mu, GEV_Parameters_Rain, PMP, partition_avg, Delta_P, 
                     error_PQ)[0]/(1.0/RI_upper_bound)
    Prob_Plus = CDF_GEV(Avg_Plus, GEV_Parameters_Runoff, PMP)
    RI_index = np.append(Return_Intervals, 1.0/(1.0-Prob_Plus)).astype(int)
    Event_Amount = GEV_RI(RI_index, GEV_Parameters_Runoff, PMP)
    df_runoff = pd.DataFrame(data=Event_Amount, index=RI_index.astype(int), columns=['Runoff'])
    weights = (1.0/Bin_Bounds[:-1]-1.0/Bin_Bounds[1:]).astype(float)
    weights = np.append(weights, 1/RI_upper_bound)
    data = np.vstack((Bin_Bounds, np.append(Bin_Bounds[1:], np.inf), weights)).T
    df_weights = pd.DataFrame(data=data, index=RI_index.astype(int), columns=['Bin Floor', 'Bin Ceiling', 
                                                                              'Event Weight']) 
    return pd.concat([df_weights, df_runoff], axis=1)


def Scenarios_Avg_S_Median_S(df_weights_runoff: pd.DataFrame, mu: float, GEV_parameters: np.ndarray, PMP: float, 
                             partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, 
                             S_limit: float, error_PQ: float) -> pd.DataFrame:
    """Calculate median and average max potential retention scenarios for given runoff.
    """
    Runoff_Q = df_weights_runoff['Runoff'].values
    Return_Intervals_Q = df_weights_runoff.index.values.astype(int)
    Avg_S_list = [Avg_SlQ(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, error_PQ, 0.0, 
                          S_limit) for Q1 in Runoff_Q]
    R_Avg_S = [1.0/2.0*(Q+np.sqrt(Q)*np.sqrt(Q+4.0*S)+2.0*S*mu) for Q, S in zip(Runoff_Q, Avg_S_list)]
    Median_S_list = [Median_S(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, error_PQ, 
                              [(0.25, S_limit)], (0+S_limit)/3).x[0] for Q1 in Runoff_Q]
    R_Median_S = [1.0/2.0*(Q+np.sqrt(Q)*np.sqrt(Q+4.0*S)+2.0*S*mu) for Q, S in zip(Runoff_Q, Median_S_list)]
    new_data = np.vstack((Avg_S_list, R_Avg_S, Median_S_list, R_Median_S)).T
    df_SR1 = pd.DataFrame(data=new_data, index=Return_Intervals_Q, 
                          columns=['Avg. S', 'Rainfall (Avg. S)', 'Median S', 'Rainfall (Median S)']) 
    return pd.concat([df_weights_runoff, df_SR1], axis=1)


def Scenarios_low_and_high_S(df_runoff_SR1: pd.DataFrame, mu: float, GEV_parameters: np.ndarray, PMP: float, 
                             partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, 
                             S_limit: float, error_PQ: float) -> pd.DataFrame:
    """Calculate scenarios for high and low maximum potential retention.
    """
    weights_runoff = df_runoff_SR1['Event Weight'].values
    Runoff_Q = df_runoff_SR1['Runoff'].values
    Return_Intervals_Q = df_runoff_SR1.index.values.astype(int)
    Median_S_list = df_runoff_SR1['Median S'].values
    Avg_S_Lower50_list = [Avg_SlQ(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, error_PQ, 
                                  0.0, S1)/0.5 for Q1, S1 in zip(Runoff_Q, Median_S_list)]
    Avg_S_Upper50_list = [Avg_SlQ(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, error_PQ, 
                                  S1, S_limit)/0.5 for Q1, S1 in zip(Runoff_Q, Median_S_list)]
    R_Avg_S_Lower50 = [1.0/2.0*(Q+np.sqrt(Q)*np.sqrt(Q+4.0*S)+2.0*S*mu) for Q, S in zip(Runoff_Q, Avg_S_Lower50_list)]
    R_Avg_S_Upper50 = [1.0/2.0*(Q+np.sqrt(Q)*np.sqrt(Q+4.0*S)+2.0*S*mu) for Q, S in zip(Runoff_Q, Avg_S_Upper50_list)]
    new_data = np.vstack((weights_runoff*0.5, Runoff_Q, Avg_S_Lower50_list, R_Avg_S_Lower50, Avg_S_Upper50_list, 
                          R_Avg_S_Upper50)).T
    df_SR2 = pd.DataFrame(data=new_data, index=Return_Intervals_Q,
                          columns=['Event Weight', 'Runoff', 'Avg. S (Lower 50%)', 'Rainfall (Lower 50%)', 
                                   'Avg. S (Upper 50%)', 'Rainfall (Upper 50%)'])
    return df_SR2

def precip_hyetograph_nrcs(df : pd.DataFrame) -> pd.DataFrame:
    '''This function takes the dataframe precip table extracted from NOAA Atlas 14 and calculates
    the nested hyetograph for storm events classified by recurrence intervals. The function first 
    retreives the ratio of rainfall, and incremental intensity and then proceeds to get ratio, slope, and slope difference\
    and finally fits parabolic curvefrom 0 to 9 hours that passes through the ratios at 0
    , 6, and 9 hours. The function then fits curves for the remaining data until 12 hours.
    '''
    time_range = {'time':np.arange(start =0, stop = 241,step = 1)}
    ratio_to_24h = pd.DataFrame(time_range,columns = ['time']).set_index(['time'])

    dif = df.diff()
    dif.at['05m','value'] = df.at['05m','value']
    df['ratio'] = df/df.at['24h','value']
    i_val = {'05m': 12, '10m': 12, '15m': 12, '30m': 4, '60m': 2,
             '02h': 1, '03h': 1, '06h': 1/3, '12h': 1/6, '24h': 1/12}
    intensity_val = pd.DataFrame.from_dict(i_val,orient='index')
    df.insert(1,'increm_intensity',dif['value']*intensity_val[0],True)

    raw_rf = {'time':[0,6,9,10.5,11,11.5,11.75,11.875,11.917,12,12.083,12.125,
                      12.25,12.5,13,13.5,15,18,24]
              }
    raw_df = pd.DataFrame(raw_rf, columns = ['time'])
    temp_0 = 0.5 - df.sort_values('ratio',ascending=False)['ratio']*.5
    temp_12 = 0.5
    temp_24 = 1- temp_0.sort_values(0,ascending=False)
    raw_df['ratio'] = ""
    raw_df.loc[0:9,'ratio']= temp_0.values
    raw_df.loc[9:18,'ratio'] = temp_24.values
    raw_df.loc[9,'ratio'] = temp_12
    raw_df['slope_raw'] = raw_df['ratio'].diff() / raw_df['time'].diff()
    raw_df.loc[0,'slope_raw'] = 0
    raw_df['slope_dif'] = raw_df.loc[0:9]['slope_raw'].diff()

    df2 = raw_df.set_index(['time'])
    a = ((2/3)* df2.at[9.0,'ratio']-df2.at[6.0,'ratio'])/18
    b = (df2.at[6.0,'ratio']-36* a)/6
    low_12h = 4 * df.loc['24h','value']*(1/36 +2/9 * df.loc['06h','value']/df.loc['24h','value'])
    up_12h = 2/3 * df.loc['24h','value']*(5/6+2/3 * df.loc['06h','value']/df.loc['24h','value'])
    ##fix negatives
    if b < 0:
        0
    if a < 0:
        df2.at[9.0,'ratio']/81
    ##fix 0 slope
    if 18*a+b<0:
        df2.at[9.0,'ratio']/4.5
    if 18*a+b<0:
        (-1*b/18)
    a2 = (9/10.5* df2.at[10.5,'ratio'] - df2.at[9.0,'ratio'])/13.5
    b2 = (df2.at[9.0,'ratio'] -81 *a2) / 9
    ##check 2h rainfall
    up_2 = 2* df.loc['24h','value']*(0.5-(df2.at[11.5,'ratio']+ 3*df2.at[10.5,'ratio'])/4)+.01
    low_2 = 2* df.loc['24h','value']*(0.5-(3*df2.at[11.5,'ratio']+df2.at[10.5,'ratio'])/4)+.01
    if df.loc['02h','value']<low_2:
        test1 = low_2
    else:
        test1 = df.loc['02h','value']
    if df.loc['02h','value']> up_2:
        test2= up_2
    else:
        test2 = df.loc['02h','value']
    if test1 > test2:
        test3 = test1
    else:
        test3 = test2
    if test2 > test3:
        test4 = test2
    else:
        test4 = test3
    if test4>up_2:
        test_f = up_2
    else:
        test_f = test4
    a3 = 2*(df2.at[11.5,'ratio'] - 2*(0.5-0.5* test_f / df.loc['24h','value'])+ df2.at[10.5,'ratio'])
    b3 = df2.at[11.5,'ratio']- df2.at[10.5,'ratio'] - 22*a3
    c3 = (0.5 - 0.5 * test_f /df.loc['24h','value']) - 121* a3 - 11 * b3  
    ##write ratios
    ratio_to_24h['ratio'] = 0
    ratio_to_24h.loc[0:90,'ratio'] = a*np.power(ratio_to_24h.loc[0:90].index/10,2)+b*ratio_to_24h.loc[0:90].index/10
    ratio_to_24h.loc[91:105,'ratio'] = a2*np.power(ratio_to_24h.loc[91:105].index/10,2)+ b2*ratio_to_24h.loc[91:105].index/10
    ratio_to_24h.loc[106:115,'ratio'] = a3*np.power(ratio_to_24h.loc[106:115].index/10,2)+ b3*ratio_to_24h.loc[106:115].index/10 + c3
    ##extra work to get 11.6,11.7
    ratio_to_24h['slope'] = ratio_to_24h['ratio'].diff()/0.1                                                               
    if -0.867*ratio_to_24h.loc[115,'slope']+ 0.4337 < 0.399: 
        fac_116 = -0.867*ratio_to_24h.loc[115,'slope']+ 0.4337
    else:
        fac_116 = 0.399
    if -0.4917*ratio_to_24h.loc[115,'slope']+ 0.8182 < 0.799: 
        fac_117 = -0.4917*ratio_to_24h.loc[115,'slope']+ 0.8182
    else:
        fac_116 = 0.799
    ratio_to_24h.at[116,'ratio'] = df2.at[11.5,'ratio'] + fac_116 *(df2.at[11.75,'ratio']-df2.at[11.5,'ratio'])
    ratio_to_24h.at[117,'ratio'] = df2.at[11.5,'ratio'] + fac_117 *(df2.at[11.75,'ratio']-df2.at[11.5,'ratio'])                                                                                                                
    ratio_to_24h.at[118,'ratio'] = df2.at[11.75,'ratio'] + 0.4*(df2.at[11.875,'ratio']-df2.at[11.75,'ratio'])
    ratio_to_24h.at[119,'ratio'] = df2.at[11.875,'ratio'] + 0.6*(df2.at[11.917,'ratio']-df2.at[11.875,'ratio'])
    ratio_to_24h.loc[121:240,'ratio'] = 1 - ratio_to_24h.loc[0:119,'ratio'].sort_index(ascending=False).values
    ratio_to_24h.loc[120,'ratio'] = ratio_to_24h.at[121,'ratio'] - (df.at['05m','ratio'] +1/5*(df.at['10m','ratio']-df.at['05m','ratio']))
    ratio_to_24h.loc[0,'ratio'] = 0
    ratio_to_24h['slope'] = ratio_to_24h['ratio'].diff()/0.1
    ratio_to_24h.at[0,'slope'] = 0
    ratio_to_24h['t_step'] = ratio_to_24h.index*.1
    ratio_to_24h.index = ratio_to_24h.index*0.1
    return ratio_to_24h

def get_hyeto_input_data(temporal_precip_table_dir: str, event_or_quartile,
                                 display_print: bool=True) -> pd.DataFrame:
    '''Extracts the temporal distribution for precipitation frequency data for the specified duration
       from an Excel sheet and returns the dataframe with the data.  
    '''
    if type(event_or_quartile) == int:
        hyeto_precip = 'nrcs_hye_{}'.format(event_or_quartile)
        df = pd.read_excel(temporal_precip_table_dir, sheet_name= hyeto_precip, index_col=0)
        if display_print: print(display(df.head(2)))
        return df
    if type(event_or_quartile) == str:
        hyeto_precip = 'atlas_hye_{}'.format(event_or_quartile)
        df = pd.read_excel(temporal_precip_table_dir, sheet_name= hyeto_precip, index_col=0)
        weights_df = pd.read_excel(temporal_precip_table_dir, sheet_name= 'atlas_hye_weights', index_col=0)
        if display_print: print(display(df.head(2)))
        return df, weights_df

def hydro_out_to_dic(curve_df: pd.DataFrame, BCN: str) -> dict:
    '''This function takes the dataframe and adding additional data
    required for the dss file and json file creation.
    '''
    dic = {}
    df_dic = curve_df.to_dict()
    dates = list(curve_df.index)
    ordin = curve_df.index.name.title()
    events = {}
    for k, v in df_dic.items():
        if 'E' in k:
            events[k] = list(v.values())
    key ='H{0}'.format(str(24).zfill(2))
    val = {'time_idx_ordinate': ordin, 
            'run_duration_days': str(2),
            'time_idx': dates, 
            'pluvial_BC_units': 'inch/ts', 
            'BCName': {BCN: events}}         
    dic[key] = val
    return dic


def Rename_Final_Groups_Precip_Stratified(curve_weight: dict, dur: int) -> dict:
    '''Sorts the groups by their weight and then renames the groups so that
       the group with the largest weight is designed E0001 and the group with
       the next largest weight is designated E0002 (for the 6 hour duration). 
       The thounsands place is set to 0, 1, 2, 3 for the 6, 12, 24, and 96 
       hour durations, respectively. A dictionary mapping the original group
       names to the new group names is returned. 
    '''
    assert dur in [6, 12, 24, 96], "Naming convention not set for duration"
    rename_map = {}
    weights = curve_weight.values()
    dur_adj = {6:0, 12:1, 24:2, 96:3 }
    num = 1
    #for i in weights:
    for k, v in curve_weight.items():
       #if i==v:
        ID = 'E{0}{1}'.format(dur_adj[dur],str(num).zfill(3))
        rename_map[k] = ID 
        num+=1
    return rename_map   

#----------------------------------------------------------------------------------------------------------------------#
# Functions for calculating inputs to the mean precipitation curve calculation.
#----------------------------------------------------------------------------------------------------------------------#

def return_interval_data(raw_precip: pd.DataFrame, Return_Intervals_MC: np.ndarray, df_GEV_parameters: pd.DataFrame, 
                         PMP: float) -> pd.DataFrame:
    """Calculates the additional precipitation values for RI not in the original NOAA data. The additional precipitation 
       value are merged with the NOAA data. In addition for each RI, the parameters are calculated for a log-normal 
       distribution that represents the variability (uncertainty) of precipitation based on the 90-percent confidence 
       interval retrieved from the NOAA data.
    """
    Non_Exceedance_Prob = 1-1/Return_Intervals_MC 
    GEV_parameters_M = df_GEV_parameters['GEV Median'].values
    GEV_parameters_L = df_GEV_parameters['GEV Lower (90%)'].values
    GEV_parameters_U = df_GEV_parameters['GEV Upper (90%)'].values
    Precip_additional_M = PPF_GEV(Non_Exceedance_Prob, GEV_parameters_M, PMP)
    Precip_additional_L = PPF_GEV(Non_Exceedance_Prob, GEV_parameters_L, PMP)
    Precip_additional_U = PPF_GEV(Non_Exceedance_Prob, GEV_parameters_U, PMP)
    Precip_additional = np.vstack((Precip_additional_M, Precip_additional_L, Precip_additional_U)).T
    df1 = pd.DataFrame(data=Precip_additional, index=Return_Intervals_MC, 
                       columns=['Median', 'Lower (90%)', 'Upper (90%)']) 
    df2 = pd.concat([raw_precip, df1]).sort_index(kind='mergesort')
    df2['Log SD (Lower)'] = (np.log(df2['Median'].values) - np.log(df2['Lower (90%)'].values))/1.645
    df2['Log SD (Upper)'] = (np.log(df2['Upper (90%)'].values) - np.log(df2['Median'].values))/1.645
    df2['Max Log SD'] = np.maximum(df2['Log SD (Lower)'].values, df2['Log SD (Upper)'].values)
    median = df2['Median'].values
    mu_LN = np.log(median)
    SD = df2['Max Log SD'].values
    df2['mu LN'] = [mu_truncated_LN(SD1, PMP, median1, mu1).x[0] for median1, mu1, SD1 in zip(median, mu_LN, SD)]
    return df2

def mu_truncated_LN(sigma: float, PMP: float, median: float, Initial_Value: float) -> float:
    """Find the mu parameter when the median of the truncated (at the PMP) lognormal is equal to the true median value.
    """
    def objective_mu_LN(mu1: float, sigma: float, PMP: float, median: float) -> float:
        """
        """
        return np.square(median - np.exp(mu1-np.sqrt(2)*sigma*special.erfcinv(1/2*special.erfc((mu1-np.log(PMP))/
                        (np.sqrt(2)*sigma)))))
    return minimize(objective_mu_LN, Initial_Value, args = (sigma, PMP, median), method='SLSQP', 
                    bounds=[(0.0, Initial_Value*2)], options={'disp': False})


def mean_curve_input_table(CL: np.ndarray, return_interval_data: pd.DataFrame, PMP: float, 
                           outputs_dir: pl.WindowsPath) -> pd.DataFrame:
    """This function takes the return interval data and creates an input table of values for calculating the mean 
       curve. The function returns a table of precipitation values for the different AEP and confidence limits (CL) 
       based on a lognormal distribution that represents the variability of precipitation based on the 90-percent 
       confidence interval limits provided by NOAA Atlas 14.
    """
    mu_LN = return_interval_data['mu LN'].values
    SD = return_interval_data['Max Log SD'].values
    data = [stats.lognorm.ppf(CL/Norm_Constant_LN(SD1, mu1, PMP), SD1, scale = np.exp(mu1)) for mu1, SD1 in 
            zip(mu_LN, SD)]
    df_input = pd.DataFrame(data=data, columns = CL, 
                            index = 1/return_interval_data.index.values).sort_index(axis=0, ascending=True)
    df_input.index.name = 'AEP'
    df_input = df_input.drop([1])
    df_input.to_csv(outputs_dir)
    return df_input

#----------------------------------------------------------------------------------------------------------------------#
# Plotting Functions
#----------------------------------------------------------------------------------------------------------------------#

def plot_GEV_precip_curves(precip_data: pd.DataFrame, df_GEV_parameters: pd.DataFrame, PMP: float, 
                           Label1: str='') -> None:
    """This functions plots the GEV distributions and also associated GEV return frequency curves on top of the 
       precpitation curve data taken either from NOAA Atlas 14 or from the mean precipitation curve output.
    """
    color = ['r', 'k', 'k']
    _, ax = plt.subplots(1, 2, figsize=(10,4))
    for i, (_, columndata) in enumerate(df_GEV_parameters.iteritems()):
        Precip = np.linspace(PPF_GEV(1e-100, columndata.values, PMP), PPF_GEV(0.9999999, columndata.values, PMP), 1000)
        Return_Period = 1.0/(1-CDF_GEV(Precip, columndata.values, PMP))
        ax[0].plot(Precip, PDF_GEV(Precip, columndata.values, PMP), color[i] , lw=2, alpha=0.6)
        ax[1].plot(Return_Period, Precip, color[i], lw=2.5, alpha=0.6)
    for _, columndata in precip_data.iteritems():
        columndata.plot(style=['+-', 'o-', '.--'], logx=True)
    ax[0].set_xlabel(f'{Label1} [inches]')
    ax[0].set_ylabel('GEV PDF $p_R(R)$')
    ax[0].set_title(f'24-hour {Label1}')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Return Period [years]')
    ax[1].set_ylabel(f'{Label1} [inches]')
    ax[1].set_title(f'24-hour {Label1}')    
    return None

        
def plot_runoff_maxRetention_distributions(GEV_parameters_E: np.ndarray, PMP: float, fitted_cn: pd.DataFrame) -> None:
    """Plots the distribution of runoff conditional on the max potential retention and plots the distribution of the max 
       potential retention.
    """
    custom_cycler = cycler('color', ['0.1', '0.25', '0.4', '0.55']) + cycler('lw', [1, 1, 1, 1])
    S_limit = 1000.0/fitted_cn.iloc[0]['CN Lower Limit']-10.0
    alpha = fitted_cn.iloc[0]['alpha']
    beta = fitted_cn.iloc[0]['beta']
    mu = fitted_cn.iloc[0]['mu']
    Q = np.linspace(PPF_GEV(1e-100, GEV_parameters_E, PMP), PPF_GEV(0.99, GEV_parameters_E, PMP), 100)
    S = np.linspace(0.0, S_limit, 100)
    _, ax = plt.subplots(1, 2, figsize=(10,4))
    ax[0].set_prop_cycle(custom_cycler)
    SA = np.linspace(0.1, 3.5, 5)
    ax[0].plot(Q, np.transpose([PDF_QlS(Q, S1, mu, GEV_parameters_E, PMP) for S1 in SA]))
    ax[0].grid(linestyle='--')
    ax[0].set_ylim((0, 1.1))
    ax[0].set_xlabel('Runoff, Q [inches]')
    ax[0].set_ylabel('$p_Q(Q | S)$') 
    ax[0].set_title('Conditional Runoff Distribution')
    ax[1].set_prop_cycle(custom_cycler)
    ax[1].plot(S, (1.0/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit))
    ax[1].grid(linestyle='--')
    ax[1].set_xlabel('Max Potential Retention, S [inches]')
    ax[1].set_ylabel('$p_S(S)$') 
    ax[1].set_title('Max Potential Retention Distribution')
    plt.tight_layout()
    return None
    
def plot_runoff_distributions_final(GEV_parameters_Rain: np.ndarray, GEV_parameters_Runoff: np.ndarray, PMP: float, 
                                    fitted_cn: pd.DataFrame, partition_avg: np.ndarray, Delta_P: float, 
                                    error_PQ: float) -> None:
    """Plots the runoff distribution and the runoff return frequency curve in comparison to the original rainfall return 
       frequency curve.
    """
    mu = fitted_cn.iloc[0]['mu']
    Q1 = np.linspace(0.01, 6, 1000)
    Return_Period = np.geomspace(1.1, 10**7, 100000)
    Precip = GEV_RI(Return_Period, GEV_parameters_Rain, PMP)
    Runoff = GEV_RI(Return_Period, GEV_parameters_Runoff, PMP)
    _, ax = plt.subplots(1, 2, figsize=(10,4))
    ax[0].grid(linestyle='--')
    ax[0].set_xlabel('Runoff, Q [inches]')
    ax[0].set_ylabel('$p_Q(Q)$')
    ax[0].set_title('Runoff Distribution')
    ax[0].plot(Q1, PDF_Q(Q1, mu, GEV_parameters_Rain, PMP, partition_avg, Delta_P, error_PQ), lw = 1, color = '0.1')
    ax[1].set_xscale('log')
    ax[1].grid(linestyle='--')
    ax[1].set_ylim((0, PMP))
    ax[1].set_xlabel('Return Period [years]')
    ax[1].set_ylabel('Depth [inches]')
    ax[1].set_title('24-hour Event') 
    ax[1].plot(Return_Period, Runoff, 'b', lw=2, alpha=0.45, label = 'Runoff')
    ax[1].plot(Return_Period, Precip, 'r', lw=2, alpha=0.6, label='Rainfall')
    ax[1].legend()
    plt.tight_layout()
    return None
    
def plot_max_potential_retention_cond_runoff(GEV_parameters_E: np.ndarray, PMP: float, fitted_cn: pd.DataFrame, 
                                             partition_avg: np.ndarray, Delta_P: float, error_PQ: float) -> None:
    """Plots the distribution of the max potential retention conditional on different runoff values.
    """
    custom_cycler = cycler('color', ['0.1', '0.25', '0.4', '0.55']) + cycler('lw', [1, 1, 1, 1])
    S_limit = 1000.0/fitted_cn.iloc[0]['CN Lower Limit']-10
    alpha = fitted_cn.iloc[0]['alpha']
    beta = fitted_cn.iloc[0]['beta']
    mu = fitted_cn.iloc[0]['mu']
    S1 = np.linspace(0.01, S_limit, 1000)
    QA = np.linspace(0.5, PMP, 50)
    PDF_S = np.transpose([PDF_SlQ(S1, Q1, mu, GEV_parameters_E, PMP, partition_avg, Delta_P, alpha, beta, S_limit, 
                                  error_PQ) for Q1 in QA])
    _, ax = plt.subplots(1, 1, figsize=(6,4))
    ax.set_prop_cycle(custom_cycler)
    ax.plot(S1, PDF_S)
    ax.grid(linestyle = '--')
    ax.set_ylim((0, np.max(PDF_S)))
    ax.set_xlabel('Max Potential Retention, S [inches]')
    ax.set_ylabel('$p_S(S | Q)$')
    ax.set_title('Conditional Max Potential Retention Distribution')
    plt.tight_layout()
    return None