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
#import papermill as pm
#import scrapbook as sb
from zipfile import ZipFile
from datetime import datetime
logging.basicConfig(level=logging.ERROR)
from collections import Counter, OrderedDict
from IPython.display import display, Markdown

import numpy as np
import pandas as pd
from scipy import stats
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

#New additions from the previous hydromet.py file
from cycler import cycler
from scipy.integrate import quad
from scipy import optimize


#---------------------------------------------------------------------------#

'''Functions called by EventTable_Stratified.ipynb. 
'''

#---------------------------------------------------------------------------#

# SCS-CN runoff formula
def Q_SCS(R,CN,mu):
    S=1000/CN-10
    return (R-mu*S)**2/(R-mu*S+S)

#Define normalization constant to account for distribution being truncated at the PMP
def Norm_Constant_GEV(x: np.ndarray, PMP: float)-> float:
    ''''Constant for distribution truncation at the PMP value 
    ''' 
    return 1/stats.genextreme.cdf(PMP, x[2], x[0], x[1])

def Norm_Constant_LN(SD: float, mu: float, PMP: float)-> float:
    ''''Constant for distribution truncation at the PMP value 
    ''' 
    return 1/stats.lognorm.cdf(PMP, SD, scale = np.exp(mu))


def PDF_GEV(R, x: np.ndarray, PMP: float)-> float:
    return Norm_Constant_GEV(x, PMP)*stats.genextreme.pdf(R, x[2], x[0], x[1])

def CDF_GEV(R, x: np.ndarray, PMP: float)-> float:
    return Norm_Constant_GEV(x, PMP)*stats.genextreme.cdf(R, x[2], x[0], x[1])

def PPF_GEV(P, x: np.ndarray, PMP: float)-> float:
    return stats.genextreme.ppf(P/Norm_Constant_GEV(x, PMP) ,x[2], x[0], x[1])

def Fit_GEV_Parameters(df: pd.DataFrame, GEV_Parameters: np.ndarray, bounds, ID: str, PMP: float)->pd.DataFrame:
    def objective_func_GEV(x: np.ndarray)-> float: 
        ''''Calculates the sum of the squared residuals between the return interval
        and return interval calculated from the GEV CDF with the differences
        normalized by the return interval. 
        ''' 
        #Function provides the sqaured error between the return interal, RI, and GEV calculated RI, \
        #with the difference normalized by the return interval
        return sum(np.square( ( RI-1/(1-CDF_GEV(row[ID], x, PMP) ) )/RI ) \
               for RI, row in df.iterrows() )
    solution = minimize(objective_func_GEV, GEV_Parameters, method='SLSQP', bounds=bounds, options={ 'disp': True})
    df_GEV_parameters=pd.DataFrame(data=solution.x, index=["mu", "sigma", "xi"], columns=["GEV {}".format(ID)])
    return df_GEV_parameters


def Avg_R_integrand(R: float, GEV_parameters: np.ndarray, PMP: float)-> float:
    return R*PDF_GEV(R, GEV_parameters, PMP)

def Avg_R(lower_bound: float, upper_bound: float, GEV_parameters: np.ndarray, PMP: float)-> float:
    return quad(Avg_R_integrand, lower_bound, upper_bound, args=(GEV_parameters, PMP))

#Functions for the rainfall event weights

#Rainfall as a function of the return interval (RI)
def R_f(RI: float, GEV_parameters: np.ndarray, PMP: float)-> float:
        return PPF_GEV( 1-1/RI, GEV_parameters, PMP)

def objective_func_bound_rainfall(RI_lower: float, RI_upper: float, RI_middle: float, GEV_parameters: np.ndarray, PMP: float)->float:
    ''''Calculates the square of the error between the average rainfall
        calculated form the bin floor and ceiling (given in terms of RI)
        and the raifall of the return period of interest
     ''' 
    return (np.square(Avg_R( R_f(RI_lower, GEV_parameters, PMP), R_f(RI_upper, GEV_parameters, PMP), GEV_parameters, PMP)[0] \
                      /(1.0/RI_lower-1.0/RI_upper) - R_f(RI_middle, GEV_parameters, PMP) ) )

#Function for function the rainfall  floor given the bin ceiling and average return interval of the bin
def Bound_rainfall_L(RI_upper: float, RI_middle: float, GEV_parameters: np.ndarray, initial_value: np.ndarray, PMP: float)-> float:
    return minimize(objective_func_bound_rainfall, initial_value, \
                    args = (RI_upper, RI_middle, GEV_parameters, PMP),\
                    method='SLSQP',bounds= [(1.0, RI_upper)], options={ 'disp': False})

#Function for function the rainfall bin ceilving given the bin floor and average return interval of the bin
def Bound_rainfall_U(RI_lower: float, RI_middle: float, GEV_parameters: np.ndarray, initial_value: np.ndarray)-> float:
    return minimize(objective_func_bound_rainfall, initial_value, \
                    args = (RI_lower, RI_middle, GEV_parameters),\
                    method='SLSQP', bounds = [(1.0, RI_upper_bound)], options={ 'disp': False})


#Runoff Conditional on Max. Potential Retention
#Q=Runoff, S=Max. Potential Retention, and mu= initial abstraction parameter
def PDF_QlS(Q: float, S: float, mu: float, GEV_parameters: np.ndarray, PMP: float)-> float:
    return (Q+2*S+np.sqrt(Q)*np.sqrt(Q+4*S))/(2*np.sqrt(Q)*np.sqrt(Q+4*S))* \
            PDF_GEV(1/2*(Q+np.sqrt(Q)*np.sqrt(Q+4*S)+2*S*mu), GEV_parameters, PMP)

#The runoff distribution consists of a continuous part plus a discrete probablity of zero runoff
#The continuous runoff distribution
def PDF_Q(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float)-> float: 
    return sum(Delta_P*PDF_QlS(Q, S_avg_partition, mu, GEV_parameters, PMP) for S_avg_partition in partition_avg)

#The discrete probability of zero runoff (integrand)
def Qzero_integrand(S: float, mu: float, alpha: float, beta: float, S_limit: float, GEV_parameters: np.ndarray, PMP: float)-> float:
    return ( CDF_GEV(S*mu, GEV_parameters, PMP) \
            -CDF_GEV(0, GEV_parameters, PMP)) \
            *(1/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit)
#The discrete probabilility of zero runoff (integrated)
def P_Qzero(mu: float, alpha: float, beta: float, S_limit: float, GEV_parameters: np.ndarray, PMP: float)-> float:
    return quad(Qzero_integrand, 0, S_limit, args =(mu, alpha, beta, S_limit, GEV_parameters, PMP)) 

#Runoff cumulative distribution function
def CDF_Q(Q: float, mu: float, alpha: float, beta: float, S_limit: float, GEV_parameters: np.ndarray, PMP: float,  partition_avg: np.ndarray, Delta_P: float)-> float:
    return quad(PDF_Q, 0, Q, args=(mu,  GEV_parameters, PMP, partition_avg, Delta_P))[0]\
                +P_Qzero(mu, alpha, beta, S_limit, GEV_parameters, PMP)[0]

def Avg_Q_integrand(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float)-> float:
    return Q*PDF_Q(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P)

def Avg_Q(lower_bound: float, upper_bound: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float)-> float:
    return quad(Avg_Q_integrand, lower_bound, upper_bound, args=(mu, GEV_parameters, PMP, partition_avg, Delta_P))


#Find runoff event weights

#Runoff as a function of the return interval (RI)
def Q_f(RI: float, tck_RI_Q)-> float:
        return interpolate.splev(RI, tck_RI_Q, der=0)

def objective_func_bound_runoff_L(RI_lower: float, RI_upper: float, RI_middle: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, tck_RI_Q)-> float:
    ''''Calculates the square of the error between the average runoff
        calculated form the bin floor and ceiling (given in terms of RI)
        and the runoff of the return period of interest
     ''' 
    return (np.square(Avg_Q( Q_f(RI_lower, tck_RI_Q), Q_f(RI_upper,  tck_RI_Q), mu, \
                            GEV_parameters, PMP, partition_avg, Delta_P)[0]/(1.0/RI_lower -1.0/RI_upper) \
                      - Q_f(RI_middle, tck_RI_Q) ) )

#Function for function the runoff bin floor given the bin ceiling and average return interval of the bin
def Bound_L(RI_upper: float, RI_middle: float, mu: float, GEV_parameters: np.ndarray, PMP, partition_avg: np.ndarray, Delta_P: float, initial_value: np.ndarray, tck_RI_Q)-> float:
    return minimize(objective_func_bound_runoff_L, initial_value, \
                    args = (RI_upper, RI_middle, mu, GEV_parameters, PMP,\
                            partition_avg, Delta_P, tck_RI_Q),\
                    method='SLSQP',bounds= [(1.0, RI_upper)], options={ 'disp': False})

#Define distribution of Max. Potential Retention
def PDF_S(S: float, alpha: float, beta: float, S_limit: float)-> float:
    return (1/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit)

#Define the integrand for finding the average value over each patition
def S_avg_integrand(S: float, alpha: float, beta: float, S_limit: float)-> float:
    return S*PDF_S(S, alpha, beta, S_limit)

#Define the integration for the average value over each patition
def S_avg_partition(alpha: float, beta: float, S_limit: float, lower_bound: float, upper_bound: float)-> float:
    return quad(S_avg_integrand, lower_bound, upper_bound, args=(alpha, beta, S_limit))

#Define the distribution the maximum potential Retention Conditional on runoff
def PDF_SlQ(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, S_limit: float)-> float:
    return PDF_QlS(Q, S, mu, GEV_parameters, PMP)*PDF_S(S, alpha, beta, S_limit)/PDF_Q(Q, mu,  GEV_parameters, PMP, partition_avg, Delta_P)

def CDF_SlQ(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, S_limit: float)-> float:
    return quad(PDF_SlQ, 0, S, args=(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit))[0]

def Avg_S1Q_integrand(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, S_limit: float)-> float:
    return S*PDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit)

def Avg_SlQ(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, S_limit: float, lower_bound: float, upper_bound: float)-> float:
    return quad(Avg_S1Q_integrand, lower_bound, upper_bound, args=(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit))

#Find the Median Value of S, i.e., the max potential retention, given a value of runoff
def objective_func_median_S(S: float, Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, S_limit: float)-> float:
    ''''Calculates the square of the error between the the CDF value
        and the median value of 0.5
     ''' 
    return np.square(CDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit)-0.5)

def Median_S(Q: float, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, S_limit: float, bounds: list, Initial_Value: np.ndarray)-> float:
    print('Calculating Median S for Runoff = %s' %(Q))
    return minimize(objective_func_median_S, Initial_Value, \
                    args = (Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit),\
                    method='SLSQP',bounds=bounds, options={ 'disp': False})

#Calculate average value of the max. potential retention for n partitions of the distribution
def partition_S_avgs(n_partition: int, alpha: float, beta: float, S_limit: float)->np.ndarray:
    #Define Probability for each partition
    Delta_P=1/n_partition
    #Define the bounds in probability space
    Bounds  = np.linspace(0, 1, n_partition+1)
    #Define the bounds in terms of the max. potential retention
    Bounds_S = S_limit*stats.beta(alpha, beta).ppf(Bounds)
    #Define the lower and upper bounds for each partition
    Bounds_Lower = Bounds_S[:-1]
    Bounds_Upper = Bounds_S[1:]
    #Find the average value over each partition, which is the basis for the runoff distribution
    partition_avg = np.transpose([S_avg_partition(alpha, beta, S_limit, lower, upper)[0] \
                              for lower, upper in zip(Bounds_Lower, Bounds_Upper)])/Delta_P
    return partition_avg

#Calculate the weights of the rainfall events.
def weights_Rainfall(Return_Intervals: np.ndarray, GEV_parameters: np.ndarray, PMP: float, RI_upper_bound: float, NOAA_precip: pd.DataFrame, ID: str, RI_data: np.ndarray, CN: float, mu: float)->pd.DataFrame:
    '''RI_data is the return intervals from which values for the rainfall are taken directly from the input data (NOAA_precip) 
        instead of being calculated from the fitted GEV.
    '''
    Size = Return_Intervals.size #Length of Return Intervals array
    #Create array for bounds of bins for each event
    Bin_Bounds_R_topdown = np.zeros(Size+1)
    #Populate the top bound for a topdown calculation of the event weights.
    Bin_Bounds_R_topdown[Size] = RI_upper_bound
    for i in range(0, Size):
        Bin_Bounds_R_topdown[Size-i-1] = \
        Bound_rainfall_L(Bin_Bounds_R_topdown[Size-i], \
                         Return_Intervals[Size-i-1], \
                         GEV_parameters, \
                         np.array([ Return_Intervals[Size-i-1]*.2 ]), \
                         PMP).x[0] 
        print('Bin Ceiling = %s, Bin Average = %s, Bin Floor = %s' \
                %(Bin_Bounds_R_topdown[Size-i], Return_Intervals[Size-i-1], Bin_Bounds_R_topdown[Size-i-1]))
    #Calculate the Average/typical precipitation between the RI upper bound and PMP
    lower_bound = PPF_GEV( 1-1/ RI_upper_bound, GEV_parameters, PMP)
    Avg_PlusR =  Avg_R(lower_bound, PMP, GEV_parameters, PMP)[0]/(1/RI_upper_bound)
    Prob_Plus =  CDF_GEV(Avg_PlusR, GEV_parameters, PMP)
    #New RI index with additional event
    RI_index = np.append(Return_Intervals, 1/(1-Prob_Plus)).astype(int)
    #weights of events   
    weights_R_topdown = (1.0/Bin_Bounds_R_topdown[:-1]-1.0/Bin_Bounds_R_topdown[1:]).astype(float)
    weights_R_topdown = np.append(weights_R_topdown, 1/RI_upper_bound)
    data= np.vstack(( Bin_Bounds_R_topdown, \
                 np.append(Bin_Bounds_R_topdown[1:], np.inf),\
                 weights_R_topdown)).T
    #dataframe of weights
    df_weights =  pd.DataFrame(data=data, index=RI_index, \
                           columns=['Bin Floor', 'Bin Celing','Event Weight']) 
    #Based on GEV  calculate the  precipitation for each RI
  
    #Determine RI that are required where values are not provided in the input data
    RI_index_calc = RI_index[ np.isin(RI_index, RI_data, invert=True) ]
    Precip_calculate = PPF_GEV(1-1/RI_index_calc, GEV_parameters, PMP)
    
    df2=pd.DataFrame(data = Precip_calculate , index = RI_index_calc.astype(int), \
                 columns=[ID])  
    #dataframe of input values
    df_R_NOAA_E =  NOAA_precip[ NOAA_precip.index.isin(RI_data) ]    #NOT USED pd.DataFrame(NOAA_precip.loc[RI_NOAA_L : RI_NOAA_U][ID])
    #Combine input values and distribution based values
    df_precip = df_R_NOAA_E.append(df2)
    #Drop column of median values
    df_precip = df_precip.drop('Median',axis=1)
    #Calculate the runoff
    Q = Q_SCS( df_precip[ID].to_numpy(), CN, mu)
    #append Runoff to the table
    df_precip['Runoff'] = Q 
    return pd.concat([df_weights,  df_precip], axis=1)

#Calculate runoff and runoff weights
def runoff(Return_Intervals: np.ndarray, RI_upper_bound: float, mu: float, GEV_parameters, PMP: float, alpha: float, beta: float, S_limit: float, partition_avg: np.ndarray, Delta_P: float, error_PQ: float)->tuple:
    #Define Runoff as a function of RI based on cubic spline interpolation
    n_partitions_Q = 40 #30 was too little, so increased to 40
    #Determine
    Q_line = np.linspace(.001, PMP - 0.1, n_partitions_Q+1)
    Return_PeriodQ= 1/(1- np.transpose([error_PQ + CDF_Q(Q, mu, alpha, beta, S_limit, GEV_parameters, PMP, partition_avg, Delta_P) for Q in Q_line]))
    #Define Runoff as a function of the return interval with a cublic spline interpolation
    tck_RI_Q = interpolate.splrep(Return_PeriodQ, Q_line)
    #Define return interval as a function of the runoffl with a cublic spline interpolation
    tck_Q_RI = interpolate.splrep(Q_line, Return_PeriodQ)
    #Define runoff event bounds
    Size = Return_Intervals.size
    #Create array for bounds of bins for each event
    Bin_Bounds = np.zeros(Size+1)
    #Populate the top bound for a topdown calculation of the event weights.
    Bin_Bounds[Size] = RI_upper_bound
    for i in range(0, Size):
        Bin_Bounds[Size-i-1] = Bound_L(Bin_Bounds[Size-i], Return_Intervals[Size-i-1],\
                                       mu, GEV_parameters, PMP, partition_avg,\
                                       Delta_P, np.array([1.01]), tck_RI_Q ).x[0] 
        print('Bin Ceiling = %s, Bin Floor %s' %(Bin_Bounds[Size-i],Bin_Bounds[Size-i-1] ) )
    #Calculate the Average/typical precipitation for the upper bin bound and PMP
    lower_bound = interpolate.splev(RI_upper_bound, tck_RI_Q, der=0)
    Avg_Q_Upper = Avg_Q(lower_bound, PMP, mu, GEV_parameters,\
                        PMP, partition_avg, Delta_P)[0]/(1/RI_upper_bound)
    Prob_Q_Upper = 1- (error_PQ + CDF_Q( Avg_Q_Upper, mu, alpha, \
                                        beta, S_limit, GEV_parameters, \
                                        PMP, partition_avg, Delta_P))
    #Append event return interval to the list of return intervals
    Return_Intervals_Q = np.append(Return_Intervals,1/Prob_Q_Upper)
    #Calculate Runoff for Return Intervals
    Runoff_Q = interpolate.splev(Return_Intervals_Q, tck_RI_Q, der=0)
    #Dataframe of results
    df_runoff = pd.DataFrame(data=Runoff_Q ,index=Return_Intervals_Q.astype(int), columns=['Runoff'])
    #weights of events   
    weights = (1.0/Bin_Bounds[:-1]-1.0/Bin_Bounds[1:]).astype(float)
    weights = np.append(weights, 1/RI_upper_bound)
    data= np.vstack(( Bin_Bounds, \
                 np.append(Bin_Bounds[1:], np.inf),\
                 weights)).T
    #dataframe of weights
    df_weights =  pd.DataFrame(data=data, index=Return_Intervals_Q.astype(int), \
                           columns=['Bin Floor', 'Bin Celing','Event Weight']) 
    return tck_RI_Q, tck_Q_RI, pd.concat([ df_weights, df_runoff], axis=1)

#Calculate median and average max. potential retention scenarios for given runoff
def Scenarios_Avg_S_Median_S(df_weights_runoff: pd.DataFrame, mu: float, GEV_parameters, PMP: float, partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, S_limit: float)->pd.DataFrame:
    Runoff_Q = df_weights_runoff['Runoff'].to_numpy()
    Return_Intervals_Q = df_weights_runoff.index.to_numpy().astype(int)
    #Calculate Average value of the max. potential retention
    Avg_S_list = [Avg_SlQ(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit,0,S_limit)[0] for Q1 in Runoff_Q]
    #Calculate the Associated Average Rainfall
    R_Avg_S = [1/2*(Q+np.sqrt(Q)*np.sqrt(Q+4*S)+2*S*mu) for Q, S in zip(Runoff_Q, Avg_S_list) ]
    #Bounds and initial value for finding the Median max. potential retention
    bounds= [(.25, S_limit)]
    Initial_Value = np.array([1.5])
    #Find the median value of the max. potential retention
    Median_S_list = [Median_S(Q1, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, bounds, Initial_Value).x[0] for Q1 in Runoff_Q]
    #Calculate the associated rainfall for the median value of the max. potential retention
    R_Median_S = [1/2*(Q+np.sqrt(Q)*np.sqrt(Q+4*S)+2*S*mu) for Q, S in zip(Runoff_Q, Median_S_list) ]
    #Save results as a dataframe
    new_data = np.vstack(( Avg_S_list, R_Avg_S, Median_S_list, R_Median_S)).T
    df_SR1 = pd.DataFrame(data=new_data, index=Return_Intervals_Q.astype(int), columns=['Avg. S', 'Rainfall','Median S', 'Rainfall']) 
    return  pd.concat([df_weights_runoff, df_SR1], axis=1)

#Calculate scenarios for high and low maximum potential retention
def Scenarios_low_and_high_S(df_runoff_SR1: pd.DataFrame, mu: float, GEV_parameters: np.ndarray, PMP: float, partition_avg: np.ndarray, Delta_P: float, alpha: float, beta: float, S_limit: float)->pd.DataFrame:
    weights_runoff = df_runoff_SR1['Event Weight'].to_numpy()
    Runoff_Q =  df_runoff_SR1['Runoff'].to_numpy()
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


# Calculate additional values for the mean curve and merge with the raw precip data from NOAA
def Mean_Curve_RI_data(raw_precip: pd.DataFrame, Return_Intervals_MC: np.ndarray, df_GEV_parameters: pd.DataFrame, PMP: float)->pd.DataFrame:
    Non_Exceedance_Prob = 1-1/Return_Intervals_MC 
    GEV_parameters_M = df_GEV_parameters['GEV Median'].to_numpy()
    GEV_parameters_L = df_GEV_parameters['GEV Lower (90%)'].to_numpy()
    GEV_parameters_U = df_GEV_parameters['GEV Upper (90%)'].to_numpy()
    Precip_additional_M = PPF_GEV(Non_Exceedance_Prob, GEV_parameters_M, PMP)
    Precip_additional_L = PPF_GEV(Non_Exceedance_Prob, GEV_parameters_L, PMP)
    Precip_additional_U = PPF_GEV(Non_Exceedance_Prob, GEV_parameters_U, PMP)
    #
    Precip_additional = np.vstack((Precip_additional_M, \
                               Precip_additional_L,Precip_additional_U)).T
    df1=pd.DataFrame(data=Precip_additional, index= Return_Intervals_MC, \
                 columns=['Median', 'Lower (90%)', 'Upper (90%)'])  # 
    df2= pd.concat([raw_precip, df1]).sort_index(kind='mergesort')
    df2['Log SD (Lower)'] = (np.log(df2['Median'].to_numpy()) - np.log(df2['Lower (90%)'].to_numpy()))/1.645
    df2['Log SD (Upper)'] = (np.log(df2['Upper (90%)'].to_numpy()) - np.log(df2['Median'].to_numpy()) )/1.645
    df2['Max Log SD'] = np.maximum(df2['Log SD (Lower)'].to_numpy(), df2['Log SD (Upper)'].to_numpy())
    df2['mu LN'] = np.log(df2['Median'].to_numpy()) 
    return df2

##Find GEV parameters for NOAA Data

def GEV_parameters_Fit(raw_precip: pd.DataFrame, ID: str, PMP: float)->pd.DataFrame:
    #Extract return intervals (in years)
    year = raw_precip.index.to_numpy()
    #Weights  of different NOAA Atlas value from the 1 to 1000 year event
    weights = np.append(1/year[:-1]-1/year[1:], 1/year[-1])
    Avg   = ( weights*raw_precip[ID]).sum()
    GEV_parameters   = np.array([Avg*.8 , 0.5, -0.25])
    bounds   = (( Avg*0.7, Avg*1.0), (.01, 1), (-0.5, 0))
    df_GEV_parameters=Fit_GEV_Parameters(raw_precip, GEV_parameters, bounds, ID, PMP)
    return df_GEV_parameters

################################
##### Plotting Functions
################################

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

