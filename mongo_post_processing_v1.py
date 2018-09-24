#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post processing script containing:
    get_Qes(filepath):      input output dir and run extension, return saturated Qes
    find_params(filepath):              input 'parameters' filepath, return wanted values
    find_gamma(filepath):               input 'omega' filepath, return gamma value

    get_scanlog(filepath):  input scan.log file, return scan parameter, growth 
                               rate, and frequency
    
@author: Austin Blackmon
"""

from mongo_file_handling_v1 import *
import numpy as np

def get_nrg(out_dir, extension):
    #modified from IFSedge/get_nrg.py
    ncols=10
    time=np.empty(0,dtype='float')
    nrg0=np.empty((1,ncols))
    nrg1=np.empty((0,ncols),dtype='float')
    
    par = Parameters()
    par.Read_Pars(out_dir + '\\parameters_' + extension)
    pars = par.pardict 
    
    nspec = 0
    for key in pars:
        if key.find('name') != -1:
            nspec = nspec + 1    
    
    if nspec<=2:
        nrg2=np.empty((0,10),dtype='float')
    if nspec<=3:
        nrg2=np.empty((0,10),dtype='float')
        nrg3=np.empty((0,10),dtype='float')
        
    f=open(out_dir + '\\nrg_' + extension,'r')
    nrg_in=f.read()

    nrg_in_lines=nrg_in.split('\n')
    for j in range(len(nrg_in_lines)):
        if nrg_in_lines[j] and j % (nspec+1) == 0:
            time=np.append(time,nrg_in_lines[j])
        elif nrg_in_lines[j] and j % (nspec+1) == 1:
            nline=nrg_in_lines[j].split()
            for i in range(ncols):
                nrg0[0,i]=nline[i]
            nrg1=np.append(nrg1,nrg0,axis=0)
        elif nspec>=2 and nrg_in_lines[j] and j % (nspec+1) ==2:
            nline=nrg_in_lines[j].split()
            for i in range(ncols):
                nrg0[0,i]=nline[i]
            nrg2=np.append(nrg2,nrg0,axis=0)
        elif nspec==3 and nrg_in_lines[j] and j % (nspec+1) ==3:
            nline=nrg_in_lines[j].split()
            for i in range(ncols):
                nrg0[0,i]=nline[i]
            nrg3=np.append(nrg3,nrg0,axis=0)

    if nspec==1:
        return time,nrg1
    elif nspec==2:
        return time,nrg1,nrg2
    else:
        return time,nrg1,nrg2,nrg3

def get_Qes(out_dir, extension):
    #generate arrays of time and Qes values
    nrg = get_nrg(out_dir, extension)
    time = nrg[0]
    Qes = nrg[1][:,-4]
    
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

def find_params(filepath):
    #if parameters not defined, set to 0
    kx, ky, omn, omt = 0, 0, 0, 0
    
    #grab parameters dictionary from ParIO.py - Parameters()
    par = Parameters()
    par.Read_Pars(filepath)
    pars = par.pardict 
    
    #if parameters defined in dict, grab their values
    if 'kymin' in pars:
        ky = pars['kymin']
    if 'kx_center' in pars:
        kx = pars['kx_center']
    if 'omn1' in pars:
        omn = pars['omn1']
    if 'omt1' in pars:
        omt = pars['omt1']
        
    return(kx, ky, omn, omt)
    
def find_gamma(filepath):
    #if gamma not found, set to 0
    gamma = 0
    
    #read scan.log
    scan = open(filepath, 'r') 
    scan = scan.read()
    
    #grab only gamma from 'omega' file and return value
    gamma = float(scan.split()[1])
    return(gamma)
         
def get_scanlog(filepath):
    #generate arrays of scan parameter, growth rate, and omega values
    scan_param = np.genfromtxt(filepath, usecols=(2))
    growth_rate = np.genfromtxt(filepath, usecols=(4))
    omega = np.genfromtxt(filepath, usecols=(5))
    
    #output arrays for scan_param, growth_rate, omega
    return(scan_param, growth_rate, omega)