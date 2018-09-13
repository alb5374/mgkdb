#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post processing script containing:
    get_Qes(filepath):      input formatted_nrg file, return time and Qes vectors
    Qes_sat(Qes):           input Qes vector, return saturated Qes value
    get_scanlog(filepath):  input scan.log file, return scan parameter, growth 
                              rate, and frequency
    
@author: Austin Blackmon
"""

import numpy as np

def get_Qes(filepath):
    #generate arrays of time and Qes values
    time = np.genfromtxt(filepath, usecols=(0))
    Qes = np.genfromtxt(filepath, usecols=(6))
    
    #run Qes_sat function to find and return saturated Qes value
    Qes_saturated = Qes_sat(Qes)
    
    return(time,Qes_saturated)

def Qes_sat(Qes):
    last = len(Qes)
    step = int(len(Qes)*0.1)
    for i in range(0,last-step):
        #find saturation over range of 100 time steps
        varience = np.var(Qes[i:i+step])
        
        #test value for Qes_sat
        if varience < 0.05: 
            #find and return saturated Qes value
            varience = np.var(Qes[i:last])
            Qes_saturated = np.mean(Qes[i:last])
            return(Qes_saturated)
        
    return('No saturated state found')
        
def get_scanlog(filepath):
    #generate arrays of scan parameter, growth rate, and omega values
    scan_param = np.genfromtxt(filepath, usecols=(2))
    growth_rate = np.genfromtxt(filepath, usecols=(4))
    omega = np.genfromtxt(filepath, usecols=(5))
    
    #output arrays for scan_param, growth_rate, omega
    return(scan_param, growth_rate, omega)