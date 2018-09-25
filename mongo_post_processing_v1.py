#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post processing script containing:
    get_Qes(filepath):      input output dir and run suffix, return saturated Qes
    find_params(filepath):              input 'parameters' filepath, return wanted values
    find_gamma(filepath):               input 'omega' filepath, return gamma value

    get_scanlog(filepath):  input scan.log file, return scan parameter, growth 
                               rate, and frequency
    
@author: Austin Blackmon
"""

from mongo_file_handling_v1 import *
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
    
def get_quasilinear(filepath):
        #if parameters not defined, set to 0
    kx, ky  = 0, 0
    
    #grab parameters dictionary from ParIO.py - Parameters()
    par = Parameters()
    par.Read_Pars(filepath)
    pars = par.pardict 
    
    #if parameters defined in dict, grab their values
    if 'kymin' in pars:
        ky = pars['kymin']
    if 'kx_center' in pars:
        kx = pars['kx_center']
        
def get_omega_from_field(out_dir, suffix):   
#    parser=op.OptionParser(description='Calculates growth rate and frequncy from field file.')
#    parser.add_option('--apar','-a',action='store_const',const=1,help = 'Calculate from Apar instead of phi.')
#    options,args=parser.parse_args()
#    if len(args)!=1:
#        exit("""
#    Please include run number as argument (e.g., 0001)."
#        \n""")
#    calc_from_apar=options.apar
#    suffix = args[0]
    
    par = Parameters()
    par.Read_Pars(out_dir+'parameters_'+suffix)
    pars = par.pardict
    if pars['n_spec'] == 1:
        time, nrge = get_nrg(out_dir, suffix)
    elif pars['n_spec'] == 2:
        time, nrgi, nrge = get_nrg(suffix)
    elif pars['n_spec'] == 3:
        time, nrgi, nrge, nrg2 = get_nrg(suffix)
    else:
        sys.exit("n_spec must be 1,2,3.")
    
    print("Calculating growth rate from phi.")
    plt.semilogy(time,nrgi[:,0])
    plt.xlabel('time')
    plt.show()
    
    tstart = float(raw_input("Enter start time: "))
    tend = float(raw_input("Enter end time: "))
    
    field = fieldfile(out_dir+'field_'+suffix,pars)
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
    
    print("imax",imax)
    
    time = np.empty(0)
    for i in range(istart,iend):
        #field.set_time(field.tfld[i],i)
        field.set_time(field.tfld[i])
        phi = np.append(phi,field.phi()[imax[0],0,imax[1]])
        if pars['n_fields'] > 1:
            apar = np.append(apar,field.apar()[imaxa[0],0,imaxa[1]])
        time = np.append(time,field.tfld[i])
         
    #plt.semilogy(time,np.abs(phi))
    #plt.semilogy(time,np.abs(apar))
    #plt.show()
    #omega = fd_d1_o4(np.log(phi),time)
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
            print('omega',omega)
            omega = np.delete(omega,0)
            time = np.delete(time,0)
    
    gam_avg = np.average(np.real(omega))
    om_avg = np.average(np.imag(omega))
    print("Gamma:",gam_avg)
    print("Omega:",om_avg)
    
    
    if output_zeros:
        f=open(out_dir+'omega_'+suffix,'w')
        f.write(str(pars['kymin'])+'    '+str(0.0)+'    '+str(0.0)+'\n')
        f.close()
    else:
        plt.plot(time,np.real(omega),label='gamma')
        plt.plot(time,np.imag(omega),label='omega')
        plt.xlabel('t(a/cs)')
        plt.ylabel('omega(cs/a)')
        plt.legend(loc='upper left')
        plt.show()
    
        selection = int(float(raw_input('How to proceed:\n1. Accept calculation \n2. Manually enter gamma and omega \n3. Don\'t output anything\n')))
        if selection == 1:
            f=open(out_dir+'omega_'+suffix,'w')
            f.write(str(pars['kymin'])+'    '+str(gam_avg)+'    '+str(om_avg)+'\n')
            f.close()
        elif selection == 2:
            gam_avg = float(raw_input('Enter gamma: '))
            om_avg = float(raw_input('Enter omega: '))
            f=open(out_dir+'omega_'+suffix,'w')
            f.write(str(pars['kymin'])+'    '+str(gam_avg)+'    '+str(om_avg)+'\n')
            f.close()
        
