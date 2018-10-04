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
    get_quasilinear(filepath):              ***in development***
    get_omega_from_field(out_dir, suffix):  input output directory and run suffix, return
                                               omega, gamma values for nonlinear runs
                                               ***in development***
    plot_linear(out_dir,scan_param,freq):   input output directory, scan parameter, and desired
                                               frequency, saves plot
                                               possible scan_param: 'kx', 'ky', 'TiTe', 'omn', 'omt'
                                               possible freq: 'gamma', 'omega'
@author: Austin Blackmon
"""

from mgk_file_handling import *
import numpy as np
import optparse as op
import matplotlib.pyplot as plt
from fieldlib import *
from ParIO import * 
from finite_differences import *
from sys import path
from sys import exit
import os

def get_nspec(out_dir,suffix):
    #grab parameters dictionary from ParIO.py - Parameters()
    par = Parameters()
    par.Read_Pars(out_dir + '/parameters' + suffix)
    pars = par.pardict 
    
    #find 'n_spec' value in parameters dictionary
    nspec = pars['n_spec']
    
    return(nspec)
    
def get_nrg(out_dir, suffix):
    #modified from IFSedge/get_nrg.py
    
    #initializations
    ncols=10
    time=np.empty(0,dtype='float')
    nrg0=np.empty((1,ncols))
    nrg1=np.empty((0,ncols),dtype='float')
    
    #grab 'n_spec' from 'parameters'
    nspec = get_nspec(out_dir,suffix)
    
    #separate initializations for different 'n_spec' values
    if nspec<=2:
        nrg2=np.empty((0,10),dtype='float')
    if nspec<=3:
        nrg2=np.empty((0,10),dtype='float')
        nrg3=np.empty((0,10),dtype='float')
    
    
    #open 'nrg' file
    f=open(out_dir + '/nrg' + suffix,'r')
    nrg_in=f.read()

    #format 'nrg' file for reading
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

    #return 'time' and 'nrgx' arrays
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
    #if parameters not defined, set to 'None'
    kx, ky, omn, omt = 'None', 'None', 'None', 'None'
    
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
        
    #return k
    return(kx, ky, omn, omt)
    
def find_omega(filepath):
    #if gamma/omega not found, set to 'None'
    gamma, omega = 'None', 'None'
    
    ### RUN CALC_GR FUNCTION IF NOT FOUND -  ###
    
    #read scan.log
    scan = open(filepath, 'r') 
    scan = scan.read()
    
    #grab gamma and omega from 'omega' file and return values
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
    par.Read_Pars(out_dir+'/parameters'+suffix)
    pars = par.pardict
    
    #find 'n_spec' value in parameters dictionary
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
    iend = np.argmin(abs(np.array(field.tfld)-tend))    
    
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
        f=open(out_dir+'/omega'+suffix,'w')
        f.write(str(pars['kymin'])+'    '+str(0.0)+'    '+str(0.0)+'\n')
        f.close()
    else:
        plt.plot(time,np.real(omega),label='gamma')
        plt.plot(time,np.imag(omega),label='omega')
        plt.xlabel('t(a/cs)')
        plt.ylabel('omega(cs/a)')
        plt.legend(loc='upper left')
        plt.show()
    
        f=open(out_dir+'/omega'+suffix,'w')
        f.write(str(pars['kymin'])+'    '+str(gam_avg)+'    '+str(om_avg)+'\n')
        f.close()

def plot_linear(out_dir,scan_param,freq):
    #style and margin adjustments
    plt.gcf().subplots_adjust(bottom=0.2,left=0.2)
    
    #check scan_param input, set xlabel
    if scan_param == 'kx':
        xlabel = r'$k_x$'
    elif scan_param =='ky':
        xlabel = r'$k_y$'
    elif scan_param == 'TiTe':
        xlabel = r'$T_i/T_e$'   
    elif scan_param == 'omn':
        xlabel = r'$\omega_n$'   
    elif scan_param == 'omt':
        xlabel = r'$\omega_T$'
        
    #check freq input, set ylabel
    if freq == 'gamma':
        ylabel = r'$\gamma$'
        column = (4)
    elif freq == 'omega':
        ylabel = r'$\omega$'
        column = (5)
        
    #formatting
    titlesize=22
    axissize=22
    plt.figure(figsize=(10,10))
    
    #grab scan_param column and freq column from 'scan.log'
    x0 = np.genfromtxt(out_dir +'/scan.log', usecols=(2), skip_header=1)
    y0 = np.genfromtxt(out_dir +'/scan.log', usecols=column, skip_header=1) 
    
    #plot
    plt.plot(x0,y0,color='#990099',label=out_dir,marker='*',ms='14',ls='-')
    
    #axis, title labels
    plt.title(out_dir,y=1.02,fontsize=titlesize)
    plt.xlabel(xlabel, fontsize=axissize)
    plt.xticks(color='k', size=22)
    plt.ylabel(ylabel, fontsize=axissize)
    plt.yticks(color='k', size=22)
    
    #legend location
    plt.legend(loc='best',numpoints=1,fontsize=14)
    
    #save and close figure
    plt.savefig(out_dir + '/' + scan_param + '_vs_' + frequency +'.png')
    plt.savefig(out_dir + '/' + scan_param + '_vs_' + frequency +'.svg')
    plt.close()
    
    ### ADD PLOT MULTIPLE RUNS - AUTOMATE LABELING, COLORING, ETC ###

def calc_gamma(out_dir,suffix,ncols=10):
    """dens,upar,tpar,tperp,Ges,Gem,Qes,Qem,Pes,Pem"""
    
    #grab 'n_spec' from 'parameters'
    nspec = get_nspec(out_dir,suffix)
    
    if nspec == 2:
        time,nrgi,nrge=get_nrg(out_dir,suffix)
    elif nspec == 3:
        time,nrgi,nrg2,nrge=get_nrg(out_dir,suffix)
    else:
        time,nrgi = get_nrg(out_dir,suffix)

    if nrgi[-1,0]/nrgi[-2,0] < 1.0e-10:
        nrgi=np.delete(nrgi,-1,0)
        if nspec > 1:
            nrge=np.delete(nrge,-1,0)
        time=np.delete(time,-1,0)
################# AUTOMATE FINDING START TIME - INDEX AFTER D/DT >>>> 1 ####################
    start_time=time[-1]-2.0
    start_index=np.argmin(abs(time-start_time))
    ntime=len(time)-start_index

    dlogdt=np.zeros((ntime,1))
    for i in range(ntime-1):
        i0=i+start_index
        dlogdt[i,0]=0.5*(nrgi[i0+1,0]-nrgi[i0,0])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,0]+nrgi[i0,0]))
    avg_gr=np.zeros(10)

    for i in range(1):
        avg_gr[i]=np.sum(dlogdt[:-1,i])/len(dlogdt[:-1,i])

    for i in range(ntime-1):
        if nspec > 1:
            i0=i+start_index
            dlogdt[i,0]=0.5*(nrge[i0+1,0]-nrge[i0,0])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,0]+nrge[i0,0]))    
            
    for i in range(1,10):
        avg_gr[i]=np.sum(dlogdt[:-1,i-5])/len(dlogdt[:-1,i-5])
 
    
    momname=list()
    momname.append('ni')
    momname.append('Tpar_i  ')
    momname.append('Tperp_i ')
    momname.append('Qes_i   ')
    momname.append('Qem_i   ')
    momname.append('ne      ')
    momname.append('Tpar_e  ')
    momname.append('Tperp_e ')
    momname.append('Qes_e   ')
    momname.append('Qem_e  ')

    #print avg_gr
    #print "Select growth rate to keep:"
    #print "Average Growth Rates:"
    #print "0:ni      ",avg_gr[0]
    #print "1:Tpar_i  ",avg_gr[1]
    #print "2:Tperp_i ",avg_gr[2]
    #print "3:Qes_i   ",avg_gr[3]
    #print "4:Qem_i   ",avg_gr[4]
    #print "5:ne      ",avg_gr[5]
    #print "6:Tpar_e  ",avg_gr[6]
    #print "7:Tperp_e ",avg_gr[7]
    #print "8:Qes_e   ",avg_gr[8]
    #print "9:Qem_e  ",avg_gr[9]
    #print "-1:none"

    fit =  np.e**(2.0*avg_gr[0]*time[start_index:])*(nrgi[-1,0]/np.e**(2.0*avg_gr[0]*time[-1]))
    err = abs(np.sum(nrgi[start_index:,0]-fit[:])/np.sum(nrgi[start_index:,0]))

    if err > 1.0e-2:
        plt.semilogy(time,nrgi[:,0],'-x')
        plt.semilogy(time[start_index-500:],np.e**(2.0*avg_gr[0]*time[start_index-500:])*(nrgi[-1,0]/np.e**(2.0*avg_gr[0]*time[-1])),'--',color='green')
        plt.show()
        test = raw_input("Accept calculation? (y=yes)")
        if test=='y':
            return avg_gr[0]
        else:
            return calc_gamma(out_dir, suffix)
    else:
        return avg_gr[0]
def get_suffixes(out_dir):
    suffixes = []
    
    #scan files in GENE output directory, find all run suffixes, return as list
    files = next(os.walk(out_dir))[2]
    for count, name in enumerate(files, start=0):
        if name.startswith('parameters_'):
            suffix = name.split('_',1)[1]
            suffix = '_' + suffix
            suffixes.append(suffix)
        elif name.lower().startswith('parameters.dat'):
            suffixes = ['.dat']                
    return suffixes

def scan_info(out_dir, suffix):
    par = Parameters()
    par.Read_Pars('parameters')
    pars = par.pardict
#    edge_opt = pars['edge_opt']
#    print("edge_opt = ",edge_opt)
#    numscan_tot = pars['scan_dims']
    numscan = get_suffixes(out_dir)
    numscan_tot = len(numscan)
    scan_info = np.zeros((numscan_tot,14),dtype='float64')
    
    for i in range(numscan_tot):
        par0 = Parameters()
        print ("Analyzing ",suffix)
        if os.path.isfile(out_dir + '/parameters' + suffix):
            par0.Read_Pars(out_dir + '/parameters' + suffix)
            pars0 = par0.pardict
            nspec = pars0['n_spec']
            print(pars0['kymin'])
            scan_info[i,0] = pars0['kymin']
            if 'x0' in pars0:
                scan_info[i,1] = pars0['x0']
            else:
                scan_info[i,1] = pars['x0']
            if 'kx_center' in pars0:
                scan_info[i,2] = pars0['kx_center']
            else:
                scan_info[i,2] = 0.0
            if 'n0_global' in pars0:
                scan_info[i,3] = pars0['n0_global']
            else:
                scan_info[i,3] = np.nan
        else:
            par0.Read_Pars(out_dir + '/parameters' + suffix)
            pars0 = par0.pardict
            nspec = pars0['n_spec']
            scan_info[i,0] = float(str(pars0['kymin']).split()[0])
            scan_info[i,1] = float(str(pars0['x0']).split()[0])
            if 'kx_center' in pars0:
                scan_info[i,2] = float(str(pars0['kx_center']).split()[0])
            else:
                scan_info[i,2] = 0.0
            if 'n0_global' in pars0:
                scan_info[i,3] = pars0['n0_global']
            else:
                scan_info[i,3] = 0.0
        if os.path.isfile(out_dir + '/omega' + suffix):
            omega0 = np.genfromtxt(out_dir + '/omega' + suffix)
            if omega0.any() and omega0[1] != 0.0:
                scan_info[i,4]=omega0[1]
                scan_info[i,5]=omega0[2]
            elif True:
                scan_info[i,4]=calc_gamma(out_dir, suffix)
                scan_info[i,5]= 0.0
                np.savetxt(out_dir + '/omega' + suffix,[scan_info[i,0],scan_info[i,4],np.nan])
            else:
                scan_info[i,4]=np.nan
                scan_info[i,5]=np.nan
        elif calc_grs and os.path.isfile(out_dir + '/nrg' + suffix):
            scan_info[i,4]=calc_gamma(out_dir, suffix)
            scan_info[i,5] = 0.0
            np.savetxt(out_dir + 'omega' + suffix,[scan_info[i,0],scan_info[i,4],np.nan])
        else:
            scan_info[i,4]=np.nan
            scan_info[i,5]=np.nan
        
        if os.path.isfile(out_dir + '/field' + suffix):
            field = fieldfile(out_dir + '/field' + suffix,pars0)
            #field.set_time(field.tfld[-1],len(field.tfld)-1)
            field.set_time(field.tfld[-1])
            fntot = field.nz*field.nx
    
            dz = float(2.0)/float(field.nz)
            zgrid = np.arange(fntot)/float(fntot-1)*(2*field.nx-dz)-field.nx
            zgrid0 = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0
            phi = np.zeros(fntot,dtype='complex128')
            apar = np.zeros(fntot,dtype='complex128')
            phikx = field.phi()[:,0,:]
            aparkx = field.phi()[:,0,:]
            phase_fac = -np.e**(-2.0*np.pi*(0.0+1.0J)*pars0['n0_global']*pars0['q0'])
            for j in range(field.nx/2+1):
                phi[(j+field.nx/2)*field.nz:(j+field.nx/2+1)*field.nz]=field.phi()[:,0,j]*phase_fac**j
                if j < field.nx/2:
                    phi[(field.nx/2-j-1)*field.nz : (field.nx/2-j)*field.nz ]=field.phi()[:,0,-1-j]*phase_fac**(-(j+1))
                if pars0['n_fields']>1:
                    apar[(j+field.nx/2)*field.nz:(j+field.nx/2+1)*field.nz]=field.apar()[:,0,j]*phase_fac**j
                    if j < field.nx/2:
                        apar[(field.nx/2-j-1)*field.nz : (field.nx/2-j)*field.nz ]=field.apar()[:,0,-1-j]*phase_fac**(-(j+1))
        
            zavg=np.sum(np.abs(phi)*np.abs(zgrid))/np.sum(np.abs(phi))
            scan_info[i,6] = zavg
            cfunc,zed,corr_len=my_corr_func_complex(phi,phi,zgrid,show_plot=False)
            scan_info[i,7] = corr_len
            parity_factor_apar = np.abs(np.sum(apar))/np.sum(np.abs(apar))
            scan_info[i,8] = parity_factor_apar
            parity_factor_phi = np.abs(np.sum(phi))/np.sum(np.abs(phi))
            scan_info[i,9] = parity_factor_phi
    
            #KBM test with E||
            gpars,geometry = read_geometry_local(pars0['magn_geometry'][1:-1]+'_' + suffix)
            jacxB = geometry['gjacobian']*geometry['gBfield']
            if scan_info[i,5] == scan_info[i,5]:
                omega_complex = (scan_info[i,5]*(0.0+1.0J) + scan_info[i,4])
                gradphi = fd_d1_o4(phi,zgrid)
                for j in range(pars0['nx0']):
                    gradphi[pars0['nz0']*j:pars0['nz0']*(j+1)] = gradphi[pars0['nz0']*j:pars0['nz0']*(j+1)]/jacxB[:]/np.pi
                #plt.plot(gradphi)
                #plt.plot(omega_complex*apar)
                #plt.show()
                diff = np.sum(np.abs(gradphi + omega_complex*apar))
                phi_cont = np.sum(np.abs(gradphi))
                apar_cont = np.sum(np.abs(omega_complex*apar))
                scan_info[i,11] = diff/(phi_cont+apar_cont)
                #print 'diff/abs',scan_info[i,10]
            else:
                scan_info[i,11] = np.nan
            phi0 = np.empty(np.shape(phikx),dtype = 'complex') 
            apar0 = np.empty(np.shape(aparkx),dtype = 'complex') 
            #print "Shape of phikx",np.shape(phikx)
            #print "Shape of zgrid",np.shape(zgrid)
            #dummy = raw_input("Press any key:\n")
            #if edge_opt > 0.0: 
            #    for ix in range(len(phikx[0,:])):
            #        phi0[:,ix] = remove_edge_opt_complex(phikx[:,ix],edge_opt)
            #        apar0[:,ix] = remove_edge_opt_complex(aparkx[:,ix],edge_opt)
            phi0 = phikx
            apar0 = aparkx
            #Calculate <gamma_HB> / gamma
            geomfile = pars0['magn_geometry'][1:-1]+'_' + suffix
            print("geomfile",geomfile)
            zgrid_pp, Btheta_R, prefactor = get_abs_psi_prime(geomfile,'../rbsProfs',pars['x0'])
            rbs = np.genfromtxt('../rbsProfs')
            ind_rbs_x0 = np.argmin(abs(rbs[:,0]-pars['x0'])) 
            gamma_HB_norm_x0 = rbs[ind_rbs_x0,9]
            ind_z0 = np.argmin(abs(zgrid_pp)) 
            prefactor_norm = prefactor/prefactor[ind_z0]
            gamma_HB_theta = abs(gamma_HB_norm_x0*prefactor_norm)
            #plt.plot(zgrid_pp,gamma_HB_theta)
            #plt.xlabel(r'$z/\pi$',size=18)
            #plt.ylabel(r'$\gamma_{HB}(c_s/a)$',size=18)
            #plt.show()
            gamma_HB_sum = 0.0
            phi_sum = 0.0
            for ix in range(len(phi0[0,:])):
                gamma_HB_sum += np.sum(abs(phi0[:,ix])**2*gamma_HB_theta*geometry['gjacobian'])
                phi_sum += np.sum(abs(phi0[:,ix])**2*geometry['gjacobian'])
            gamma_HB_avg = gamma_HB_sum / phi_sum
            scan_info[i,12] = gamma_HB_avg
            #gamma_HB_theta = abs(gamma_HB_norm_x0*prefactor_norm)
            #theta0 = (scan_info[i,2]/(pars0['shat']*pars0['kymin']*np.pi))
            #if theta0 > 1.0:
            #    theta0 -= 2.0
            #ind_theta0 = np.argmin(abs(zgrid_pp-theta0)) 
            #print 'kx_center',scan_info[i,2]
            #print 'shat',pars0['shat']
            #print 'kymin',pars0['kymin']
            #print 'theta0',theta0
            #print 'ind_theta0',ind_theta0
            n_info[i,13] = np.min(gamma_HB_theta)
        else:
            scan_info[i,6] = np.nan
            scan_info[i,7] = np.nan
            scan_info[i,8] = np.nan
            scan_info[i,9] = np.nan
            scan_info[i,11] = np.nan
            scan_info[i,12] = np.nan
            scan_info[i,13] = np.nan
    
        if os.path.isfile(out_dir + 'nrg' + suffix):
            if nspec==1:
                tn,nrg1=get_nrg0('_' + suffix,nspec=nspec)
                scan_info[i,10]=nrg1[-1,7]/abs(nrg1[-1,6])
            elif nspec==2:
                tn,nrg1,nrg2=get_nrg0('_' + suffix,nspec=nspec)
                scan_info[i,10]=nrg2[-1,7]/(abs(nrg2[-1,6])+abs(nrg1[-1,6]))
            elif nspec==3:
                tn,nrg1,nrg2,nrg3=get_nrg0('_' + suffix,nspec=nspec)
                scan_info[i,10]=nrg3[-1,7]/(abs(nrg3[-1,6])+abs(nrg1[-1,6]))
            else:
                sys.exit("Not ready for nspec>2")
        else:
            scan_info[i,10] = np.nan
    
    print(np.shape(scan_info))
    
    f=out_dir + '/scan_info.dat'
    head = '1.kymin 2.x0 3.kx_center 4.n0_global 5.gamma(cs/a) 6.omega(cs/a) 7.<z> 8.lambda_z 9.parity(apar) 10.parity(phi) 11.QEM/QES 12.Epar cancelation 13.gamma_HB_avg 14.gamma_HB_min'
    np.savetxt(f,scan_info, header = head)