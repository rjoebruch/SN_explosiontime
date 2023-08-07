#### Functions for fitting the exploision time of a SN 
from astropy.coordinates import SkyCoord 
from astropy.table import Table, vstack, hstack
from astropy.io import ascii 
import numpy as np
import matplotlib.pylab as plt
from numpy import random
import pandas as pd

import random
import math
import scipy 
from scipy.optimize import curve_fit
from scipy.special import gamma
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.stats import median_abs_deviation

import sfdmap 
import extinction
from PyAstronomy import pyasl
from astropy.coordinates import Distance


from astropy import cosmology

from iminuit import Minuit



import pprint



def rise_typeII(jd, a, t_exp, n):
    '''
    This function defines the type of rise we want to fit to the early light curve.

    f(x) = | a(t-T_exp)**n if t >= T_exp
           | 0             if t <= T_exp
    
    parameters
    ----------
    jd     [array] time interval
    a      [float] amplitude/normalisation scale
    T_exp  [float] estimated time of zero flux
    n      [float] fit exponent
    
    
    returns
    -------
    array
    
    '''
  
    hiuv = np.where( jd-t_exp > 0, a*((jd-t_exp)**n) , 0 )

    # rise   = a*((jd+T_exp)**n)
    # #a*np.sign(jd-T_exp)*(np.abs(jd-T_exp)**n)
    # # a*((jd-T_exp)**n)
    # condi0 = ( (jd + T_exp) <= 0 )

    # #condi0 = ((rise) >= 0 )
    # rise   = condi0*rise

    # return rise
    return hiuv





def rise_typeII_mag(jd, a, t_exp, n):
    '''
    This function defines the type of rise we want to fit to the early light curve.

    f(x) = | a(t-T_exp)**n if t >= T_exp
           | 0             if t <= T_exp
    
    parameters
    ----------
    jd     [array] time interval
    a      [float] amplitude/normalisation scale
    T_exp  [float] estimated time of zero flux
    n      [float] fit exponent
    
    
    returns
    -------
    array
    
    '''
  
    # hiuv = np.where( jd-t_exp > 0, a*((jd-t_exp)**n) , 0 ) 
    ''' can't have that conditionin the mag space because it's not a matter of zero flux!!
    '''

    # rise   = a*((jd+T_exp)**n)
    # #a*np.sign(jd-T_exp)*(np.abs(jd-T_exp)**n)
    # # a*((jd-T_exp)**n)
    # condi0 = ( (jd + T_exp) <= 0 )

    # #condi0 = ((rise) >= 0 )
    # rise   = condi0*rise

    # return rise
    return  a*(jd-t_exp)**n




def to_flux(m):
    '''
    this funtion converts magnitude to flux
    
    parameters
    ----------
    f.  column or array?
    
    returns
    -------
    '''
    
    return 10**(-m/2.5)



def error_onflux_from_mag(m,dm):
    '''
    this funtion converts magnitude to flux
    
    parameters
    ----------
    f.  column or array?
    
    returns
    -------
    '''
    
    return np.sqrt( (((-np.log10(10)/2.5)*(10**(-m/2.5)))**2) *  dm**2 )


