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

#Define normalization constant to account for distribution being truncated at the PMP
def Norm_Constant(x: np.ndarray, PMP)-> float:
    ''''Constant for distribution truncation at the PMP value 
    ''' 
    return 1/stats.genextreme.cdf(PMP, x[2], x[0], x[1])

def PDF_GEV(R, x: np.ndarray, PMP)-> float:
    return Norm_Constant(x, PMP)*stats.genextreme.pdf(R ,x[2], x[0], x[1])

def CDF_GEV(R, x: np.ndarray, PMP)-> float:
    return Norm_Constant(x, PMP)*stats.genextreme.cdf(R ,x[2], x[0], x[1])

def PPF_GEV(P, x: np.ndarray, PMP)-> float:
    return stats.genextreme.ppf(P/Norm_Constant(x, PMP) ,x[2], x[0], x[1])

def Fit_GEV_Parameters(df: pd.DataFrame, GEV_Parameters: np.ndarray, bounds, ID: str, PMP: float):
    def objective_func_GEV(x: np.ndarray)-> float: 
        ''''Calculates the sum of the squared residuals between the return interval
        and return interval calculated from the GEV CDF with the differences
        normalized by the return interval. 
        ''' 
        #Function provides the sqaured error between the return interal, RI, and GEV calculated RI, \
        #with the difference normalized by the return interval
        return sum(np.square( ( RI-1/(1-CDF_GEV(row[ID], x, PMP) ) )/RI ) \
               for RI, row in df.iterrows() )
    solution = minimize(objective_func_GEV, GEV_Parameters, method='SLSQP',bounds=bounds, options={ 'disp': True})
    df_GEV_parameters=pd.DataFrame(data=solution.x, index=["mu", "sigma", "xi"], columns=["GEV {}".format(ID)])
    return df_GEV_parameters


def Avg_R_integrand(R, GEV_parameters, PMP):
    return R*PDF_GEV(R, GEV_parameters, PMP)

def Avg_R(lower_bound, upper_bound, GEV_parameters, PMP):
    return quad(Avg_R_integrand, lower_bound, upper_bound, args=(GEV_parameters, PMP))

#Functions for the rainfall event weights

#Rainfall as a function of the return interval (RI)
def R_f(RI, GEV_parameters, PMP)-> float:
        return PPF_GEV( 1-1/RI, GEV_parameters, PMP)

def objective_func_bound_rainfall(RI_lower, RI_upper, RI_middle, GEV_parameters, PMP):
    ''''Calculates the square of the error between the average rainfall
        calculated form the bin floor and ceiling (given in terms of RI)
        and the raifall of the return period of interest
     ''' 
    return (np.square(Avg_R( R_f(RI_lower, GEV_parameters, PMP), R_f(RI_upper, GEV_parameters, PMP), GEV_parameters, PMP)[0] \
                      /(1.0/RI_lower-1.0/RI_upper) - R_f(RI_middle, GEV_parameters, PMP) ) )

#Function for function the rainfall  floor given the bin ceiling and average return interval of the bin
def Bound_rainfall_L(RI_upper, RI_middle, GEV_parameters, initial_value, PMP)-> float:
    return minimize(objective_func_bound_rainfall, initial_value, \
                    args = (RI_upper, RI_middle, GEV_parameters, PMP),\
                    method='SLSQP',bounds= [(1.0, RI_upper)], options={ 'disp': False})

#Function for function the rainfall bin ceilving given the bin floor and average return interval of the bin
def Bound_rainfall_U(RI_lower, RI_middle, GEV_parameters, initial_value)-> float:
    return minimize(objective_func_bound_rainfall, initial_value, \
                    args = (RI_lower, RI_middle, GEV_parameters),\
                    method='SLSQP',bounds= [(1.0, RI_upper_bound)], options={ 'disp': False})


#Runoff Conditional on Max. Potential Retention
#Q=Runoff, S=Max. Potential Retention, and mu= initial abstraction parameter
def PDF_QlS(Q, S, mu, GEV_parameters, PMP):
    return (Q+2*S+np.sqrt(Q)*np.sqrt(Q+4*S))/(2*np.sqrt(Q)*np.sqrt(Q+4*S))* \
            PDF_GEV(1/2*(Q+np.sqrt(Q)*np.sqrt(Q+4*S)+2*S*mu), GEV_parameters, PMP)

#The runoff distribution consists of a continuous part plus a discrete probablity of zero runoff
#The continuous runoff distribution
def PDF_Q(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P): 
    return sum(Delta_P*PDF_QlS(Q, S_avg_partition, mu, GEV_parameters, PMP) for S_avg_partition in partition_avg)

#The discrete probability of zero runoff (integrand)
def Qzero_integrand(S, mu, alpha, beta, S_limit, GEV_parameters, PMP):
    return ( CDF_GEV(S*mu, GEV_parameters, PMP) \
            -CDF_GEV(0, GEV_parameters, PMP)) \
            *(1/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit)
#The discrete probabilility of zero runoff (integrated)
def P_Qzero(mu, alpha, beta, S_limit, GEV_parameters, PMP):
    return quad(Qzero_integrand, 0, S_limit, args =(mu, alpha, beta, S_limit, GEV_parameters, PMP)) 

#Runoff cumulative distribution function
def CDF_Q(Q, mu, alpha, beta, S_limit, GEV_parameters, PMP,  partition_avg, Delta_P):
    return quad(PDF_Q, 0, Q, args=(mu,  GEV_parameters, PMP, partition_avg, Delta_P))[0]\
                +P_Qzero(mu, alpha, beta, S_limit, GEV_parameters, PMP)[0]

def Avg_Q_integrand(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P):
    return Q*PDF_Q(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P)

def Avg_Q(lower_bound, upper_bound, mu, GEV_parameters, PMP, partition_avg, Delta_P):
    return quad(Avg_Q_integrand, lower_bound, upper_bound, args=(mu, GEV_parameters, PMP, partition_avg, Delta_P))


#Find runoff event weights

#Runoff as a function of the return interval (RI)
def Q_f(RI, tck_RI_Q)-> float:
        return interpolate.splev(RI, tck_RI_Q, der=0)

def objective_func_bound_runoff_L(RI_lower, RI_upper, RI_middle, mu, GEV_parameters, PMP, partition_avg, Delta_P, tck_RI_Q):
    ''''Calculates the square of the error between the average runoff
        calculated form the bin floor and ceiling (given in terms of RI)
        and the runoff of the return period of interest
     ''' 
    return (np.square(Avg_Q( Q_f(RI_lower, tck_RI_Q), Q_f(RI_upper,  tck_RI_Q), mu, \
                            GEV_parameters, PMP, partition_avg, Delta_P)[0]/(1.0/RI_lower -1.0/RI_upper) \
                      - Q_f(RI_middle, tck_RI_Q) ) )

#Function for function the runoff bin floor given the bin ceiling and average return interval of the bin
def Bound_L(RI_upper, RI_middle, mu, GEV_parameters, PMP, partition_avg, Delta_P, initial_value, tck_RI_Q)-> float:
    return minimize(objective_func_bound_runoff_L, initial_value, \
                    args = (RI_upper, RI_middle, mu, GEV_parameters, PMP,\
                            partition_avg, Delta_P, tck_RI_Q),\
                    method='SLSQP',bounds= [(1.0, RI_upper)], options={ 'disp': False})

#Define distribution of Max. Potential Retention
def PDF_S(S, alpha, beta, S_limit):
    return (1/S_limit)*stats.beta(alpha, beta).pdf(S/S_limit)

#Define the integrand for finding the average value over each patition
def S_avg_integrand(S, alpha, beta, S_limit) :
    return S*PDF_S(S, alpha, beta, S_limit)

#Define the integration for the average value over each patition
def S_avg_partition(alpha, beta, S_limit, lower_bound, upper_bound) :
    return quad(S_avg_integrand, lower_bound, upper_bound, args=(alpha, beta, S_limit))

#Define the distribution the maximum potential Retention Conditional on runoff
def PDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit):
    return PDF_QlS(Q, S, mu, GEV_parameters, PMP)*PDF_S(S, alpha, beta, S_limit)/PDF_Q(Q, mu,  GEV_parameters, PMP, partition_avg, Delta_P)

def CDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit):
    return quad(PDF_SlQ, 0, S, args=(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit))[0]

def Avg_S1Q_integrand(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit):
    return S*PDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit)

def Avg_SlQ(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, lower_bound, upper_bound):
    return quad(Avg_S1Q_integrand, lower_bound, upper_bound, args=(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit))

#Find the Median Value of S, i.e., the max potential retention, given a value of runoff
def objective_func_median_S(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit):
    ''''Calculates the square of the error between the the CDF value
        and the median value of 0.5
     ''' 
    return np.square(CDF_SlQ(S, Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit)-0.5)

def Median_S(Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit, bounds, Initial_Value):
    return minimize(objective_func_median_S, Initial_Value, \
                    args = (Q, mu, GEV_parameters, PMP, partition_avg, Delta_P, alpha, beta, S_limit),\
                    method='SLSQP',bounds=bounds, options={ 'disp': True})