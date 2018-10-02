#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File handling script for formatting output files, getting file lists, and
reading and writing to database containing:
    get_file_list(out_dir,begin):       input GENE output directory and base filepath 
                                           (nrg, energy, etc), return full list of files 
                                           in directory
    get_suffixes(out_dir):            input GENE output directory, return list of run 
                                           suffixes in the directory
    gridfs_put(filepath):               input filepath, upload  file to database, and 
                                           return object_id of uploaded file
    gridfs_read(db_file):               input database filename, return contents of file
    upload_to_mongo   
    isLinear
@author: Austin Blackmon
"""

from mgk_post_processing import *
from ParIO import * 
import numpy as np
import os
import gridfs
from pymongo import MongoClient

def get_file_list(out_dir,begin):
    files_list = []
    
    #unwanted filetype suffixes for general list
    bad_ext = ('.ps','.png')
    
    #scan files in GENE output directory, ignoring files in '/in_par', and return list
    for dirpath, dirnames, files in os.walk(out_dir):
        for count, name in enumerate(files, start=0):
            if name.lower().startswith(begin) and name.find('in_par') == -1  \
            and name.endswith(bad_ext) == False:
                files = os.path.join(dirpath, name)
                files_list.append(files)
                
    return files_list     


def get_suffixes(out_dir):
    suffixes = []
    
    #scan files in GENE output directory, find all run suffixes, return as list
    for dirpath, dirnames, files in os.walk(out_dir):
        for count, name in enumerate(files, start=0):
            if name.lower().startswith('codemods'):
                extension = name.split('_',1)[1]
                suffixes.append(extension)
    return suffixes

def gridfs_put(filepath):
    #set directory and filepath
    file = open(filepath, 'rb')

    #connect to 'ETG' database
    db = MongoClient().ETG

    #upload file to 'fs.files' collection
    fs = gridfs.GridFS(db)
    dbfile = fs.put(file, encoding='UTF-8', filepath=filepath)
    file.close()
    
    #grab '_id' for uploaded file
    object_id = str(dbfile)
    return(object_id)
    
def gridfs_read(db_file):
    #connect to 'ETG' database
    db = MongoClient().ETG
    
    #open 'filepath'
    fs = gridfs.GridFS(db)
    file = fs.find_one({"filepath": db_file})
    contents = file.read()
    return(contents)
    

def upload_to_mongo(out_dir, user, linear, confidence, input_heat, keywords):
    #generate file lists for relevant files
    autopar_files = get_file_list(out_dir, 'autopar')
    codemods_files = get_file_list(out_dir, 'codemods')
    energy_files = get_file_list(out_dir, 'energy')
    field_files = get_file_list(out_dir, 'field')
    geneerr_files = get_file_list(out_dir, 'geneerr')
    mom_files = get_file_list(out_dir, 'mom')
    nrg_files = get_file_list(out_dir, 'nrg')
    omega_files = get_file_list(out_dir, 'omega')
    parameters_files = get_file_list(out_dir, 'parameters')
    s_alpha_files = get_file_list(out_dir, 's_alpha')
    scanlog_files = get_file_list(out_dir, 'scan')
    vsp_files = get_file_list(out_dir, 'vsp')

    #list of all GENE output files to be uploaded
    output_files = autopar_files + codemods_files + energy_files + \
        field_files + geneerr_files + mom_files + nrg_files +  omega_files + \
        parameters_files + s_alpha_files + scanlog_files + vsp_files
    
    #upload all GENE output files to database
    object_ids = [] 
    for files in output_files:
        object_ids.append(gridfs_put(files)) 
            
    #map from_output files to object ids
    for i in range(0,len(output_files)):
        object_ids[i] = object_ids[i] + '    ' +  output_files[i]
        
    ### USE DICTIONARY ###
    scan_id = 'None'
    scanlog_id = 'None'
    codemods_id = 'None'
    submit_id = 'None'
    parameters_id = 'None'
    eqdisk_id = 'None'
    efit_id = 'None'
    autopar_id = 'None'
    energy_id = 'None'
    field_id = 'None'
    mom_id = 'None'
    nrg_id = 'None'
    vsp_id = 'None'
    omega_id = 'None'
    s_alpha_id = 'None'
    
    #for linear runs
    if linear:
        #connect to 'ETG' database, 'LinearRuns' collection
        runs = MongoClient().ETG.LinearRuns
        suffixes = get_suffixes(out_dir)
        for suffix in suffixes:
            for line in object_ids:
                if line.find('codemods_' + suffix) != -1:
                    codemods_id = line.split()[0]
                if line.find('parameters_' + suffix) != -1:
                    parameters_id = line.split()[0]
                if line.find('autopar_' + suffix) != -1:
                    autopar_id = line.split()[0]
                if line.find('energy_' + suffix) != -1:
                    energy_id = line.split()[0]
                if line.find('nrg_' + suffix) != -1:
                    nrg_id = line.split()[0]
                if line.find('omega_' + suffix) != -1:
                    omega_id = line.split()[0]
                if line.find('scan.log') != -1:
                    scanlog_id = line.split()[0]
                if line.find('s_alpha_' + suffix) != -1:
                    s_alpha_id = line.split()[0]
    
            #find relevant parameters from in/output
            gamma = find_omega(out_dir + '\\omega_' + suffix)[0]
            omega = find_omega(out_dir + '\\omega_' + suffix)[1]
            params = find_params(out_dir + '\\parameters_' + suffix)
            kx = params[0]
            ky = params[1]
            omn = params[2]
            omt = params[3]
            
            
            #document format for linear runs  
            run_data = {"user": user,
                        "run_collection_name": out_dir,
                        "run_suffix": '_' + suffix,
                        "keywords": keywords,
                        "confidence": confidence,
                        "codemods_id": codemods_id,
                        "submitcmd_id": submit_id,
                        "parameters_id": parameters_id,
                        "eqdisk_id": eqdisk_id,
                        "efit_id": efit_id,
                        "autopar_id": autopar_id,
                        "energy_id": energy_id,
                        "field_id": field_id,
                        "mom_id": mom_id,
                        "nrg_id": nrg_id,
                        "omega_id": omega_id,
                        "scanlog_id": scanlog_id,
                        "vsp_id": vsp_id,
                        "gamma": gamma,
                        "omega": omega,
                        "ky": ky,
                        "kx": kx,
                        "omt": omt,
                        "omn": omn
                        }
            #insert run_data into database
            runs.insert_one(run_data).inserted_id
    
    #for nonlinear runs
    if not linear:
        #connect to 'ETG' database, 'NonlinRuns' collection
        runs = MongoClient().ETG.NonlinRuns
        suffixes = get_suffixes(out_dir)
        for suffix in suffixes:
            for line in object_ids:
                if line.find('codemods_' + suffix) != -1:
                    codemods_id = line.split()[0]
                if line.find('parameters_' + suffix) != -1:
                    parameters_id = line.split()[0]
                if line.find('autopar_' + suffix) != -1:
                    autopar_id = line.split()[0]
                if line.find('energy_' + suffix) != -1:
                    energy_id = line.split()[0]
                if line.find('nrg_' + suffix) != -1:
                    nrg_id = line.split()[0]
                if line.find('omega_' + suffix) != -1:
                    omega_id = line.split()[0]
                if line.find('scan.log') != -1:
                    scanlog_id = line.split()[0]
                if line.find('s_alpha_' + suffix) != -1:
                    s_alpha_id = line.split()[0]
                if line.find('field_' + suffix) != -1:
                    field_id = line.split()[0]
                if line.find('mom_' + suffix) != -1:
                    mom_id = line.split()[0]
                if line.find('vsp_' + suffix) != -1:
                    vsp_id = line.split()[0]
                    
            #find relevant parameters from in/output
            Qes = get_Qes(out_dir, suffix)
            params = find_params(out_dir + '\\parameters_' + suffix)
            kx = params[0]
            ky = params[1]
            omn = params[2]  # check n_spec for suffix
            omt = params[3]
            
            #document format for nonlinear runs
            run_data = {"user": user,
                        "run_collection_name": out_dir,
                        "run_suffix": '_' + suffix,
                        "keywords": keywords,
                        "confidence": confidence,                       
                        "codemods_id": codemods_id,
                        "submitcmd_id": submit_id,
                        "parameters_id": parameters_id,
                        "eqdisk_id": eqdisk_id,
                        "efit_id": efit_id,
                        "autopar_id": autopar_id,
                        "energy_id": energy_id,
                        "field_id": field_id,
                        "mom_id": mom_id,
                        "nrg_id": nrg_id,
                        "omega_id":omega_id,
                        "scanlog_id": scanlog_id,
                        "vsp_id": vsp_id,
                        "Qes" : Qes,
                        "ky" : ky,
                        "kx" : kx,
                        "omt" : omt,
                        "omn" : omn                            
                       }
            
            #insert run_data into database
            runs.insert_one(run_data).inserted_id                

def isLinear(name):
    if os.path.isfile(name + '\\parameters'):
        par = Parameters()
        par.Read_Pars(name + '\\parameters')
        pars = par.pardict
        linear = not pars['nonlinear']
        return(linear)
    elif name.find('linear') != -1:
        linear = True 
        return(linear)
    elif name.find('nonlin') != -1:
        linear = False
        return(linear)