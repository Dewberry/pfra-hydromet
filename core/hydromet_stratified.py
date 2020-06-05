import os
import io
import re
import csv
import json
import time
import urllib
import shutil
import logging
import operator
import warnings
import pathlib as pl
from zipfile import ZipFile
from datetime import datetime
logging.basicConfig(level=logging.ERROR)
from collections import Counter, OrderedDict
from IPython.display import display, Markdown

import numpy as np
import pandas as pd
from scipy import stats, special
from nptyping import Array
from scipy.optimize import minimize
from scipy import interpolate, integrate
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import FormatStrFormatter

import fiona
import rasterio
from pyproj import Proj
import geopandas as gpd
from rasterio.mask import mask
from shapely.geometry import mapping
from cycler import cycler
from scipy.integrate import quad
from scipy import optimize


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
    return 1/stats.genextreme.cdf(PMP, x[2], x[0], x[1])


def Norm_Constant_LN(SD: float, mu: float, PMP: float) -> float:
    """Constant for distribution truncation at the PMP value. 
    """ 
    return 1/stats.lognorm.cdf(PMP, SD, scale = np.exp(mu))


def PDF_GEV(R: np.ndarray, x: np.ndarray, PMP: float) -> float:
    """
    """
    return Norm_Constant_GEV(x, PMP)*stats.genextreme.pdf(R, x[2], x[0], x[1])


def CDF_GEV(R: np.ndarray, x: np.ndarray, PMP: float) -> float:
    """
    """
    return Norm_Constant_GEV(x, PMP)*stats.genextreme.cdf(R, x[2], x[0], x[1])


def PPF_GEV(P: np.ndarray, x: np.ndarray, PMP: float) -> float:
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
    year = raw_precip.index.to_numpy()
    weights = np.append(1/year[:-1]-1/year[1:], 1/year[-1])
    Avg = (weights*raw_precip[ID]).sum()
    GEV_parameters = np.array([Avg*0.8, 0.5, -0.25])
    bounds = ((Avg*0.7, Avg*1.0), (0.01, 1.0), (-0.5, 0.0))
    df_GEV_parameters = GEV_Parameters(raw_precip, GEV_parameters, bounds, ID, PMP)
    return df_GEV_parameters


def Avg_R_integrand(R: float, GEV_parameters: np.ndarray, PMP: float) -> float:
    """This function defines the integrand for calculating an average based on the GEV distribution.
    """
    return R*PDF_GEV(R, GEV_parameters, PMP)


def Avg_R(lower_bound: float, upper_bound: float, GEV_parameters: np.ndarray, PMP: float) -> float:
    """Calculates the average value of the GEV distribution of rainfall based on an upper and lower bound.
    """
    return quad(Avg_R_integrand, lower_bound, upper_bound, args=(GEV_parameters, PMP))


def rainfall_RI(RI: float, GEV_parameters: np.ndarray, PMP: float) -> float:
    """Provides rainfall as a function of the return interval (RI).
    """
    return PPF_GEV(1-1/RI, GEV_parameters, PMP)


def objective_func_bound_rainfall(RI_lower: float, RI_upper: float, RI_middle: float, GEV_parameters: np.ndarray, 
                                  PMP: float) -> float:
    """Calculates the square of the error between the average rainfall calculated form the bin floor and ceiling (given 
       in terms of RI) and the rainfall of the return period of interest.
    """ 
    return np.square(Avg_R(rainfall_RI(RI_lower, GEV_parameters, PMP), rainfall_RI(RI_upper, GEV_parameters, PMP), 
                     GEV_parameters, PMP)[0]/(1.0/RI_lower-1.0/RI_upper) - rainfall_RI(RI_middle, GEV_parameters, PMP))


def Bound_rainfall_L(RI_upper: float, RI_middle: float, GEV_parameters: np.ndarray, initial_value: float, 
                     PMP: float) -> float:
    """Finds the rainfall bin floor given the bin ceiling and average return interval of the bin.
    """
    return minimize(objective_func_bound_rainfall, initial_value, args = (RI_upper, RI_middle, GEV_parameters, PMP),
                    method='SLSQP', bounds=[(1.0, RI_upper)], options={'disp': False})


def PDF_QlS(Q: float, S: float, mu: float, GEV_parameters: np.ndarray, PMP: float) -> float:
    """This function provides the runoff PDF conditional on Max. Potential Retention, where Q=Runoff, S=Max. Potential 
       Retention, and mu= initial abstraction parameter.
    """
    return (Q+2*S+np.sqrt(Q)*np.sqrt(Q+4*S))/(2*np.sqrt(Q)*np.sqrt(Q+4*S))*\
            PDF_GEV(1/2*(Q+np.sqrt(Q)*np.sqrt(Q+4*S)+2*S*mu), GEV_parameters, PMP)


#The runoff distribution consists of a continuous part plus a discrete probability of zero runoff
#The continuous runoff distribution
def PDF_Q(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, 
          Delta_P: float) -> float: 
    """
    """
    return sum(Delta_P*PDF_QlS(Q, S_avg_partition, mu, GEV_parameters, PMP) for S_avg_partition in partition_avg)


def Qzero_integrand(S: float, mu: float, alpha: float, beta: float, S_limit: float, GEV_parameters: np.ndarray,
                    PMP: float) -> float:
    """Defines the integrand for calculating the probability of zero runoff.
    """
    return (CDF_GEV(S*mu, GEV_parameters, PMP)-CDF_GEV(0, GEV_parameters, PMP))\
           *(1/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit)


def P_Qzero(mu: float, alpha: float, beta: float, S_limit: float, GEV_parameters: np.ndarray, PMP: float) -> float:
    """Defines discrete probability of zero runoff (integrated).
    """
    return quad(Qzero_integrand, 0, S_limit, args =(mu, alpha, beta, S_limit, GEV_parameters, PMP)) 


def CDF_Q(Q: float, mu: float, alpha: float, beta: float, S_limit: float, GEV_parameters: np.ndarray, PMP: float, 
          partition_avg: np.ndarray, Delta_P: float) -> float:
    """Defines the cumulative distribution function for runoff.
    """
    return quad(PDF_Q, 0, Q, args=(mu,  GEV_parameters, PMP, partition_avg, Delta_P))[0]\
                +P_Qzero(mu, alpha, beta, S_limit, GEV_parameters, PMP)[0]

def Avg_Q_integrand(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, 
                    Delta_P: float) -> float:
    """
    """
    return Q*PDF_Q(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P)


def Avg_Q(lower_bound: float, upper_bound: float, mu: float, GEV_parameters: np.ndarray, PMP: float, 
          partition_avg: np.ndarray, Delta_P: float) -> float:
    """
    """
    return quad(Avg_Q_integrand, lower_bound, upper_bound, args=(mu, GEV_parameters, PMP, partition_avg, Delta_P))


def runoff_RI(RI: float, tck_RI_Q: float) -> float:
    """Defines runoff as a function of the return interval (RI).
    """
    return interpolate.splev(RI, tck_RI_Q, der=0)


def objective_func_bound_runoff_L(RI_lower: float, RI_upper: float, RI_middle: float, mu: float, 
                                  GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, 
                                  tck_RI_Q: float) -> float:
    """Calculates the square of the error between the average runoff calculated form the bin floor and ceiling (given in 
       terms of RI) and the runoff of the return period of interest.
    """ 
    return (np.square(Avg_Q(runoff_RI(RI_lower, tck_RI_Q), runoff_RI(RI_upper,  tck_RI_Q), mu, GEV_parameters, PMP, 
                      partition_avg, Delta_P)[0]/(1.0/RI_lower-1.0/RI_upper)-runoff_RI(RI_middle, tck_RI_Q)))


def Bound_L(RI_upper: float, RI_middle: float, mu: float, GEV_parameters: np.ndarray, PMP: float, 
            partition_avg: np.ndarray, Delta_P: float, initial_value: np.ndarray, tck_RI_Q: float) -> float:
    """Calculates runoff bin floor given the bin ceiling and average return interval of the bin.
    """
    return minimize(objective_func_bound_runoff_L, initial_value, 
                    args = (RI_upper, RI_middle, mu, GEV_parameters, PMP, partition_avg, Delta_P, tck_RI_Q), 
                    method='SLSQP', bounds=[(1.0, RI_upper)], options={'disp': False})


def PDF_S(S: float, alpha: float, beta: float, S_limit: float)-> float:
    """Defines the distribution of the Max. Potential Retention.
    """
    return (1/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit)


def S_avg_integrand(S: float, alpha: float, beta: float, S_limit: float) -> float:
    """Defines the integrand for finding the average value over each partition.
    """
    return S*PDF_S(S, alpha, beta, S_limit)


def S_avg_partition(alpha: float, beta: float, S_limit: float, lower_bound: float, upper_bound: float) -> float:
    """Defines the integration for the average value over each partition.
    """
    return quad(S_avg_integrand, lower_bound, upper_bound, args=(alpha, beta, S_limit))


def PDF_SlQ(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, 
            Delta_P: float, alpha: float, beta: float, S_limit: float) -> float:
    """Defines the PDF of the max. potential retention, S, conditional on runoff, Q.
    """
    return PDF_QlS(Q, S, mu, GEV_parameters, PMP)*PDF_S(S, alpha, beta, S_limit)/PDF_Q(Q, mu,  GEV_parameters, PMP, 
                                                                                       partition_avg, Delta_P)


def CDF_SlQ(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, 
            Delta_P: float, alpha: float, beta: float, S_limit: float) -> float:
    """Defines the CDF of the max potential retention, S, conditional on runoff Q.
    """
    return quad(PDF_SlQ, 0, S, args=(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit))[0]


def Avg_S1Q_integrand(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, 
                      Delta_P: float, alpha: float, beta: float, S_limit: float) -> float:
    """Defines the integrand for calculating the average max potential retention.
    """
    return S*PDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit)


def Avg_SlQ(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, 
            alpha: float, beta: float, S_limit: float, lower_bound: float, upper_bound: float) -> float:
    """Derives the average values of the max potential retention by integrating between an upper and lower bound.
    """
    return quad(Avg_S1Q_integrand, lower_bound, upper_bound, args=(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, 
                                                                   alpha, beta, S_limit))


def objective_func_median_S(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, 
                            partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, 
                            S_limit: float) -> float:
    """Calculates the square of the error between the the CDF value and the median value of 0.5.
    """ 
    return np.square(CDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit)-0.5)


def Median_S(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, 
             alpha: float, beta: float, S_limit: float, bounds: list, Initial_Value: np.ndarray) -> float:
    """
    """
    print('Calculating Median S for Runoff = %s' %(Q))
    return minimize(objective_func_median_S, Initial_Value, \
                    args = (Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit),\
                    method='SLSQP',bounds=bounds, options={'disp': False})


def partition_S_avgs(n_partition: int, alpha: float, beta: float, S_limit: float) -> np.ndarray:
    """Calculates the average value of the max. potential retention for n partitions of the distribution.
    """
    Delta_P = 1.0/n_partition
    Bounds  = np.linspace(0, 1, n_partition+1)
    Bounds_S = S_limit*stats.beta(alpha, beta).ppf(Bounds)
    Bounds_Lower = Bounds_S[:-1]
    Bounds_Upper = Bounds_S[1:]
    partition_avg = np.transpose([S_avg_partition(alpha, beta, S_limit, lower, upper)[0] \
                              for lower, upper in zip(Bounds_Lower, Bounds_Upper)])/Delta_P
    return partition_avg


def weights_Rainfall(Return_Intervals: np.ndarray, GEV_parameters: np.ndarray, PMP: float, RI_upper_bound: float, 
                     NOAA_precip: pd.DataFrame, ID: str, RI_data: np.ndarray, CN: float, mu: float) -> pd.DataFrame:
    """Calculate the weights of the rainfall events. RI_data is the return intervals from which values for the rainfall 
       are taken directly from the input data (NOAA_precip) instead of being calculated from the fitted GEV.
    """
    Size = Return_Intervals.size
    Bin_Bounds_R_topdown = np.zeros(Size+1)
    Bin_Bounds_R_topdown[Size] = RI_upper_bound
    for i in range(0, Size):
        Bin_Bounds_R_topdown[Size-i-1] = Bound_rainfall_L(Bin_Bounds_R_topdown[Size-i], Return_Intervals[Size-i-1], 
                                                          GEV_parameters, Return_Intervals[Size-i-1]*0.2, PMP).x[0] 
    lower_bound = PPF_GEV(1-1.0/RI_upper_bound, GEV_parameters, PMP)
    Avg_PlusR = Avg_R(lower_bound, PMP, GEV_parameters, PMP)[0]/(1.0/RI_upper_bound)
    Prob_Plus = CDF_GEV(Avg_PlusR, GEV_parameters, PMP)  # Not sure why you are trying to find this value...
    RI_index = np.append(Return_Intervals, 1.0/(1-Prob_Plus)).astype(int)
    weights_R_topdown = (1.0/Bin_Bounds_R_topdown[:-1]-1.0/Bin_Bounds_R_topdown[1:]).astype(float)
    weights_R_topdown = np.append(weights_R_topdown, 1.0/RI_upper_bound)
    data = np.vstack((Bin_Bounds_R_topdown, np.append(Bin_Bounds_R_topdown[1:], np.inf), weights_R_topdown)).T
    df_weights = pd.DataFrame(data=data, index=RI_index, columns=['Bin Floor', 'Bin Ceiling', 'Event Weight']) 
    RI_index_calc = RI_index[np.isin(RI_index, RI_data, invert=True)]
    Precip_calculate = PPF_GEV(1-1/RI_index_calc, GEV_parameters, PMP)
    df2 = pd.DataFrame(data = Precip_calculate, index = RI_index_calc, columns=[ID])  
    df_R_NOAA_E = NOAA_precip[NOAA_precip.index.isin(RI_data)].copy()
    df_precip = df_R_NOAA_E.append(df2)
    df_precip = df_precip.drop('Median', axis=1)
    Q = Q_SCS(df_precip[ID].values, CN, mu)
    df_precip['Runoff'] = Q 
    return pd.concat([df_weights, df_precip], axis=1)


def runoff(Return_Intervals: np.ndarray, RI_upper_bound: float, mu: float, GEV_parameters: np.ndarray, PMP: float, 
           alpha: float, beta: float, S_limit: float, partition_avg: np.ndarray, Delta_P: float, error_PQ: float, 
           n_partitions_Q: int=40) -> tuple:
    """Calculate runoff and runoff probability weights and returns a dataframe of runoff values and runoff weights. The 
       function also returns funoff as a function of RI based on cubic spline interpolation and RI as a function of RI 
       based on a cubic spline interpolation. 
    """
    Q_line = np.linspace(.001, PMP - 0.001, n_partitions_Q+1)
    Return_PeriodQ = 1/(1-np.transpose([error_PQ + CDF_Q(Q, mu, alpha, beta, S_limit, GEV_parameters, PMP, partition_avg, Delta_P) for Q in Q_line]))
    #Define Runoff as a function of the return interval with a cublic spline interpolation
    tck_RI_Q = interpolate.splrep(Return_PeriodQ, Q_line)
    #Define return interval as a function of the runoffl with a cublic spline interpolation
    tck_Q_RI = interpolate.splrep(Q_line, Return_PeriodQ)
    Size = Return_Intervals.size
    Bin_Bounds = np.zeros(Size+1)
    Bin_Bounds[Size] = RI_upper_bound
    for i in range(0, Size):
        Bin_Bounds[Size-i-1] = Bound_L(Bin_Bounds[Size-i], Return_Intervals[Size-i-1], mu, GEV_parameters, PMP, 
                                       partition_avg, Delta_P, np.array([1.01]), tck_RI_Q ).x[0] 
        print('Bin Ceiling = {}, Bin Floor {}'.format(Bin_Bounds[Size-i], Bin_Bounds[Size-i-1]))
    #Calculate the Average/typical precipitation for the upper bin bound and PMP
    lower_bound = interpolate.splev(RI_upper_bound, tck_RI_Q, der=0)
    Avg_Q_Upper = Avg_Q(lower_bound, PMP, mu, GEV_parameters, PMP, partition_avg, Delta_P)[0]/(1/RI_upper_bound)
    Prob_Q_Upper = 1-(error_PQ + CDF_Q(Avg_Q_Upper, mu, alpha, beta, S_limit, GEV_parameters, PMP, partition_avg,
                                       Delta_P))
    Return_Intervals_Q = np.append(Return_Intervals, 1/Prob_Q_Upper)
    Runoff_Q = interpolate.splev(Return_Intervals_Q, tck_RI_Q, der=0)
    df_runoff = pd.DataFrame(data=Runoff_Q ,index=Return_Intervals_Q.astype(int), columns=['Runoff'])
    #weights of events   
    weights = (1.0/Bin_Bounds[:-1]-1.0/Bin_Bounds[1:]).astype(float)
    weights = np.append(weights, 1/RI_upper_bound)
    data = np.vstack((Bin_Bounds, np.append(Bin_Bounds[1:], np.inf), weights)).T
    df_weights =  pd.DataFrame(data=data, index=Return_Intervals_Q.astype(int), 
                               columns=['Bin Floor', 'Bin Ceiling','Event Weight']) 
    return tck_RI_Q, tck_Q_RI, pd.concat([df_weights, df_runoff], axis=1)


def Scenarios_Avg_S_Median_S(df_weights_runoff: pd.DataFrame, mu: float, GEV_parameters: np.ndarray, PMP: float, 
                             partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, 
                             S_limit: float) -> pd.DataFrame:
    """Calculate median and average max. potential retention scenarios for given runoff.
    """
    Runoff_Q = df_weights_runoff['Runoff'].to_numpy()
    Return_Intervals_Q = df_weights_runoff.index.to_numpy().astype(int)
    Avg_S_list = [Avg_SlQ(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit,0,S_limit)[0] for Q1 in Runoff_Q]
    R_Avg_S = [1/2*(Q+np.sqrt(Q)*np.sqrt(Q+4*S)+2*S*mu) for Q, S in zip(Runoff_Q, Avg_S_list) ]
    #Bounds and initial value for finding the Median max. potential retention
    bounds= [(0.25, S_limit)]
    Initial_Value = np.array([1.5])
    #Find the median value of the max. potential retention
    Median_S_list = [Median_S(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, bounds, Initial_Value).x[0] for Q1 in Runoff_Q]
    R_Median_S = [1/2*(Q+np.sqrt(Q)*np.sqrt(Q+4*S)+2*S*mu) for Q, S in zip(Runoff_Q, Median_S_list) ]
    #Save results as a dataframe
    new_data = np.vstack(( Avg_S_list, R_Avg_S, Median_S_list, R_Median_S)).T
    df_SR1 = pd.DataFrame(data=new_data, index=Return_Intervals_Q.astype(int), columns=['Avg. S', 'Rainfall', 'Median S', 'Rainfall']) 
    return pd.concat([df_weights_runoff, df_SR1], axis=1)


def Scenarios_low_and_high_S(df_runoff_SR1: pd.DataFrame, mu: float, GEV_parameters: np.ndarray, PMP: float, 
                             partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, 
                             S_limit: float) -> pd.DataFrame:
    """Calculate scenarios for high and low maximum potential retention.
    """
    weights_runoff = df_runoff_SR1['Event Weight'].to_numpy()
    Runoff_Q = df_runoff_SR1['Runoff'].to_numpy()
    Return_Intervals_Q = df_runoff_SR1.index.to_numpy().astype(int)
    Median_S_list = df_runoff_SR1['Median S'].to_numpy()
    Avg_S_Lower50_list = [Avg_SlQ(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha,beta,S_limit,0,S1)[0]/0.5 for Q1, S1 in zip(Runoff_Q, Median_S_list)]
    Avg_S_Upper50_list = [Avg_SlQ(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha,beta,S_limit,S1,S_limit)[0]/0.5 for Q1, S1 in zip(Runoff_Q, Median_S_list)]
    R_Avg_S_Lower50 = [1/2*(Q+np.sqrt(Q)*np.sqrt(Q+4*S)+2*S*mu) for Q, S in zip(Runoff_Q, Avg_S_Lower50_list) ]
    R_Avg_S_Upper50 = [1/2*(Q+np.sqrt(Q)*np.sqrt(Q+4*S)+2*S*mu) for Q, S in zip(Runoff_Q, Avg_S_Upper50_list) ]
    new_data = np.vstack(( weights_runoff*.5,    Runoff_Q ,  Avg_S_Lower50_list, R_Avg_S_Lower50, \
                      weights_runoff*.5,    Runoff_Q ,  Avg_S_Upper50_list, R_Avg_S_Upper50)).T
    df_SR2 = pd.DataFrame(data=new_data, index=Return_Intervals_Q.astype(int),\
                      columns=['Event Weight', 'Runoff','Avg. S (Lower 50%)', 'Rainfall', \
                               'Event Weight', 'Runoff','Avg. S (Upper 50%)', 'Rainfall'])
    return df_SR2


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


def mean_curve_input_table(CL: np.ndarray, return_interval_data: pd.DataFrame, PMP: float, outputs_dir: pl.WindowsPath, 
                           Project_Area: str, Pluvial_Model: str, BCN: str) -> pd.DataFrame:
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
    df_input.to_csv(outputs_dir/'Mean_Curve_Input_{0}_{1}_{2}.csv'.format(Project_Area, Pluvial_Model, BCN))
    return df_input


#----------------------------------------------------------------------------------------------------------------------#
# Plotting Functions
#----------------------------------------------------------------------------------------------------------------------#

def plot_GEV_precip_curves(precip_data: pd.DataFrame, df_GEV_parameters: pd.DataFrame, PMP: float) -> None:
    """This functions plots the GEV distributions and also associated GEV return frequency curves on top of the 
       precpitation curve data taken either from NOAA Atlas 14 or from the mean precipitation curve output.
    """
    color = ['r', 'k', 'k']
    _, ax = plt.subplots(1, 2, figsize=(10,4))
    for i, (_, columndata) in enumerate(df_GEV_parameters.iteritems()):
        Precip = np.linspace(PPF_GEV(1e-100, columndata.values, PMP), PPF_GEV(0.999999, columndata.values, PMP), 100)
        Return_Period = 1/(1-CDF_GEV(Precip, columndata.values, PMP))
        ax[0].plot(Precip, PDF_GEV(Precip, columndata.values, PMP), color[i] , lw=2, alpha=0.6, label='genextreme pdf')
        ax[1].plot(Return_Period, Precip, color[i], lw=2.5, alpha=0.6, label='genextreme pdf')
    for _, columndata in precip_data.iteritems():
        columndata.plot(style=['+-', 'o-', '.--'], logx=True)
    ax[0].set_xlabel('Rainfall [inches]')
    ax[0].set_ylabel('GEV PDF $p_R(R)$')
    ax[0].set_title('24-hour Rainfall')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Return Period [years]')
    ax[1].set_ylabel('Rainfall [inches]')
    ax[1].set_title('24-hour Rainfall')    
    return None

        
def plot_runoff_maxRetention_distributions(GEV_parameters_E: np.ndarray, PMP: float, fitted_cn: pd.DataFrame) -> None:
    """Plots the distribution of runoff conditional on the max potential retention and plots the distribution of the max 
       potential retention.
    """
    custom_cycler = cycler('color', ['.1', '.25', '.4', '.55']) + cycler('lw', [1, 1, 1, 1])
    S_limit = 1000.0/fitted_cn.iloc[0]['CN Lower Limit']-10.0
    alpha = fitted_cn.iloc[0]['alpha']
    beta = fitted_cn.iloc[0]['beta']
    mu = fitted_cn.iloc[0]['mu']
    Q = np.linspace(PPF_GEV(10**-100, GEV_parameters_E, PMP), PPF_GEV(0.99, GEV_parameters_E, PMP), 100)
    S = np.linspace(0, S_limit, 100)
    _, ax = plt.subplots(1, 2, figsize=(10,4))
    ax[0].set_prop_cycle(custom_cycler)
    SA = np.linspace(0.1, 3.5 , 5.0)
    ax[0].plot(Q, np.transpose([PDF_QlS(Q, S1, mu, GEV_parameters_E, PMP) for S1 in SA]))
    ax[0].grid(linestyle='--')
    ax[0].set_ylim((0, 1.1))
    ax[0].set_xlabel('Runoff, Q (inches)'),
    ax[0].set_ylabel('$p_Q(Q | S)$'), 
    ax[0].set_title('Conditional Runoff Distribution')
    ax[1].set_prop_cycle(custom_cycler)
    SA = np.linspace(0.1, 3.5 , 5.0)
    ax[1].plot(S, (1/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit))
    ax[1].grid(linestyle='--')
    ax[1].set_xlabel('Max. Potential Retention, S (inches)')
    ax[1].set_ylabel('$p_S(S)$') 
    ax[1].set_title('Max. Potential Retention Distribution')
    plt.tight_layout()
    plt.show()
    return None
    

def plot_runoff_distributions_final(GEV_parameters_E: np.ndarray, PMP: float, fitted_cn: pd.DataFrame, 
                                    partition_avg: np.ndarray, Delta_P: float, tck_RI_Q: float) -> None:
    """Plots the runoff distribution and the runoff return frequency curve in comparison to the original rainfall return 
       frequency curve.
    """
    mu = fitted_cn.iloc[0]['mu']
    Q1 = np.linspace(.01, 6, 1000)
    Return_PeriodQ = np.linspace(1, 10**5.4, 100000)
    Precip = np.linspace(PPF_GEV(10**-100, GEV_parameters_E, PMP), PPF_GEV(.999999, GEV_parameters_E, PMP), 100)
    Return_Period = 1/(1-CDF_GEV(Precip,  GEV_parameters_E , PMP))
    _, ax = plt.subplots(1, 2, figsize=(10,4))
    custom_cycler = cycler('color', ['.1', '.25', '.4', '.55']) + cycler('lw', [1, 1, 1, 1])
    ax[0].set_prop_cycle(custom_cycler)
    ax[0].grid(linestyle='--')
    ax[0].set_xlabel('Runoff, Q (inches)')
    ax[0].set_ylabel('$p_Q(Q)$')
    ax[0].set_title('Runoff Distribution')
    ax[0].plot(Q1, PDF_Q(Q1, mu, GEV_parameters_E, PMP, partition_avg, Delta_P))
    ax[1].set_xscale('log')
    ax[1].set_prop_cycle(custom_cycler)
    ax[1].grid(linestyle='--')
    ax[1].set_ylim((0, PMP))
    ax[1].set_xlabel('Runoff, Q (inches)')
    ax[1].set_ylabel('Rainfall (red) and Runoff (gray), inches$')
    ax[1].set_title('24-hour Event')
    ax[1].plot(Return_PeriodQ, interpolate.splev(Return_PeriodQ, tck_RI_Q, der=0), 'b', lw =2.9, alpha=.45)
    ax[1].plot(Return_Period, Precip, 'r-', lw=5, alpha=0.6, label='genextreme pdf')
    plt.tight_layout()
    plt.show()
    return None
    

def plot_max_potential_retention_cond_runoff(GEV_parameters_E: np.ndarray, PMP: float, fitted_cn: pd.DataFrame, 
                                             partition_avg: np.ndarray, Delta_P: float) -> None:
    """Plots the distribution of the max potential retention conditional on different runoff values.
    """
    S_limit = 1000/fitted_cn.iloc[0]['CN Lower Limit']-10
    alpha = fitted_cn.iloc[0]['alpha']
    beta = fitted_cn.iloc[0]['beta']
    mu = fitted_cn.iloc[0]['mu']
    S1 = np.linspace(.01, S_limit, 1000)
    _, ax = plt.subplots(1, 1, figsize=(10,7))
    custom_cycler = cycler('color', ['.1', '.25', '.4', '.55']) + cycler('lw', [1, 1, 1, 1])
    ax.set_prop_cycle(custom_cycler)
    QA = np.linspace(.5, 10 , 50)
    ax.plot(S1, np.transpose([PDF_SlQ(S1, Q1, mu, GEV_parameters_E, PMP, partition_avg, Delta_P, alpha, beta, S_limit) for Q1 in QA]))
    ax.grid(linestyle = '--')
    ax.set_ylim((0, .45))
    ax.set_xlabel('Max. Potential Retention, S (inches)')
    ax.set_ylabel('$p_S(S | Q)$')
    ax.set_title('Conditional Max. Potential Retention Distribution')
    plt.tight_layout()
    plt.show()
    return None