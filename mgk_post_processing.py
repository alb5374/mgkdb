#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post processing script containing:
    get_nrg(out_dir, suffix):               input put
    get_Qes(filepath):                      input 'nrg' filepath return saturated Qes for
                                               nonlinear runs
    find_params(filepath):                  input 'parameters' filepath, return wanted values
    find_omega(filepath):                   input 'omega' filepath, return omega, gamma values
                                               for linear runs
    get_scanlog(filepath):                  input scan.log file, return scan parameter, growth 
                                               rate, and frequency for linear runs
    get_omega_from_field(out_dir, suffix):  input output directory and run suffix, return
                                               omega, gamma values for nonlinear runs
                                               ***in development***
    get_quasilinear(filepath):              ***in development***
@author: Austin Blackmon
"""

from mgk_file_handling import *
import numpy as np
import optparse as op
import matplotlib.pyplot as plt
from fieldlib import *
from ParIO import * 
from finite_differences import *

def get_nrg(out_dir, suffix):
    #modified from IFSedge/get_nrg.py
    ncols=10
    time=np.empty(0,dtype='float')
    nrg0=np.empty((1,ncols))
    nrg1=np.empty((0,ncols),dtype='float')
    
    par = Parameters()
    par.Read_Pars(out_dir + '\\parameters_' + suffix)
    pars = par.pardict 
    
    nspec = pars['n_spec']
#    for key in pars:
#        if key.find('name') != -1:
#            nspec = nspec + 1    
    
    if nspec<=2:
        nrg2=np.empty((0,10),dtype='float')
    if nspec<=3:
        nrg2=np.empty((0,10),dtype='float')
        nrg3=np.empty((0,10),dtype='float')
        
    f=open(out_dir + '\\nrg_' + suffix,'r')
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

def get_Qes(out_dir, suffix):
    '''
    input output dir and run suffix, return saturated Qes
    '''
    
    #generate arrays of time and Qes values
    nrg = get_nrg(out_dir, suffix)
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
    
def find_omega(filepath):
    #if gamma not found, set to 0
    gamma = 0
    omega = 0
    
    #read scan.log
    scan = open(filepath, 'r') 
    scan = scan.read()
    
    #grab only gamma from 'omega' file and return value
    gamma = float(scan.split()[1])
    omega = float(scan.split()[2])
    return(gamma,omega)
         
def get_scanlog(filepath):
    #generate arrays of scan parameter, growth rate, and omega values
    scan_param = np.genfromtxt(filepath, usecols=(2))
    growth_rate = np.genfromtxt(filepath, usecols=(4))
    omega = np.genfromtxt(filepath, usecols=(5))
    
    #output arrays for scan_param, growth_rate, omega
    return(scan_param, growth_rate, omega)
    
def get_quasilinear(filepath):
        #if parameters not defined, set to 0
    gamma, kx, ky  = 0, 0, 0
    
    #grab parameters dictionary from ParIO.py - Parameters()
    par = Parameters()
    par.Read_Pars(filepath)
    pars = par.pardict 
    
    #if parameters defined in dict, grab their values
    if 'kymin' in pars:
        ky = pars['kymin']
    if 'kx_center' in pars:
        kx = pars['kx_center']   
    gamma = find_gamma(filepath)
    quasi_gamma =  gamma / (kx**2 + ky**2)
        
def get_omega_from_field(out_dir, suffix):
    calc_from_apar=0
    par = Parameters()
    par.Read_Pars(out_dir+'\\parameters_'+suffix)
    pars = par.pardict
    if pars['n_spec'] == 1:
        time, nrgi = get_nrg(out_dir, suffix)
    elif pars['n_spec'] == 2:
        time, nrgi, nrge = get_nrg(suffix)
    elif pars['n_spec'] == 3:
        time, nrgi, nrge, nrg2 = get_nrg(suffix)
    else:
        sys.exit("n_spec must be 1,2,3.")
    
    tstart = 24.0
    tend = 25.0

    
    
    istart = np.argmin(abs(np.array(field.tfld)-tstart))
    print("istart,start_time",istart,field.tfld[istart])
    iend = np.argmin(abs(np.array(field.tfld)-tend))
    print("iend,end_time",iend,field.tfld[iend])
    
    
    #field.set_time(field.tfld[-1],len(field.tfld)-1)
    field.set_time(field.tfld[-1])
    imax = np.unravel_index(np.argmax(abs(field.phi()[:,0,:])),(field.nz,field.nx))
    phi = np.empty(0,dtype='complex128')
    if pars['n_fields'] > 1:
        imaxa = np.unravel_index(np.argmax(abs(field.apar()[:,0,:])),(field.nz,field.nx))
        apar = np.empty(0,dtype='complex128')
    
    time = np.empty(0)
    for i in range(istart,iend):
        field.set_time(field.tfld[i])
        phi = np.append(phi,field.phi()[imax[0],0,imax[1]])
        if pars['n_fields'] > 1:
            apar = np.append(apar,field.apar()[imaxa[0],0,imaxa[1]])
        time = np.append(time,field.tfld[i])
         

    if len(phi) < 2.0:
        output_zeros = True
        omega = 0.0+0.0J
    else:
        output_zeros = False
        if calc_from_apar:
            if pars['n_fields'] < 2:
                stop
            omega = np.log(apar/np.roll(apar,1))
            dt = time - np.roll(time,1)
            omega /= dt
            omega = np.delete(omega,0)
            time = np.delete(time,0)
        else:
            omega = np.log(phi/np.roll(phi,1))
            dt = time - np.roll(time,1)
            omega /= dt
            omega = np.delete(omega,0)
            time = np.delete(time,0)
    
    gam_avg = np.average(np.real(omega))
    om_avg = np.average(np.imag(omega))
    
    
    if output_zeros:
        f=open(out_dir+'\\omega_'+suffix,'w')
        f.write(str(pars['kymin'])+'    '+str(0.0)+'    '+str(0.0)+'\n')
        f.close()
    else:
        plt.plot(time,np.real(omega),label='gamma')
        plt.plot(time,np.imag(omega),label='omega')
        plt.xlabel('t(a/cs)')
        plt.ylabel('omega(cs/a)')
        plt.legend(loc='upper left')
        plt.show()
    
        f=open(out_dir+'\\omega_'+suffix,'w')
        f.write(str(pars['kymin'])+'    '+str(gam_avg)+'    '+str(om_avg)+'\n')
        f.close()

        