#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Base script to handle uploading of GENE folders to database

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
    
    #connect to 'ETG' database, 'Runs' collection
    runs = MongoClient().ETG.Runs
    
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
    formatted_energy_id = 'None'
    field_id = 'None'
    mom_id = 'None'
    nrg_id = 'None'
    formatted_nrg_id = 'None'
    vsp_id = 'None'
    omega_id = 'None'
    s_alpha_id = 'None'
    diagnostics_id = 'None'
    
    #for linear runs
    if linear:
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
                        "run_suffix": out_dir + '_' + suffix,
                        "keywords": keywords,
                        "codemods": codemods_id,
                        "submitcmd": submit_id,
                        "confidence": confidence,
                        "parameters": parameters_id,
                        "eqdisk": eqdisk_id,
                        "efit": efit_id,
                        "autopar": autopar_id,
                        "energy": energy_id,
                        "field": field_id,
                        "mom": mom_id,
                        "nrg": nrg_id,
                        "omega": omega_id,
                        "scanlog": scanlog_id,
                        "vsp": vsp_id,
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
                        "codemods": codemods_id,
                        "submitcmd": submit_id,
                        "confidence": confidence,
                        "parameters": parameters_id,
                        "eqdisk": eqdisk_id,
                        "efit": efit_id,
                        "autopar": autopar_id,
                        "energy": energy_id,
                        "field": field_id,
                        "mom": mom_id,
                        "nrg": nrg_id,
                        "omega":omega_id,
                        "scanlog": scanlog_id,
                        "vsp": vsp_id,
                        "Qes" : Qes,
                        "ky" : ky,
                        "kx" : kx,
                        "omt" : omt,
                        "omn" : omn                            
                       }
            
            #insert run_data into database
            runs.insert_one(run_data).inserted_id                
