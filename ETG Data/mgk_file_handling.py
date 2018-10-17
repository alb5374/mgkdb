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
from pymongo import MongoClient
import os
import gridfs
import re

def get_file_list(out_dir,begin):
    files_list = []
    
    #unwanted filetype suffixes for general list
    bad_ext = ('.ps','.png')
    
    #scan files in GENE output directory, ignoring files in '/in_par', and return list
    files = next(os.walk(out_dir))[2]
    for count, name in enumerate(files, start=0):
        if name.startswith(begin) and name.endswith(bad_ext) == False and not os.path.isdir('in_par'):
            file = out_dir + '/' + name
            files_list.append(file)
    return files_list     


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
    
def gridfs_read(filepath):
    #connect to 'ETG' database
    db = MongoClient().ETG
    
    #open 'filepath'
    fs = gridfs.GridFS(db)
    file = fs.find_one({"filepath": filepath})
    contents = file.read()
    return(contents)
    
def isLinear(name):
    if os.path.isfile(name + '/parameters'):
        par = Parameters()
        par.Read_Pars(name + '/parameters')
        pars = par.pardict
        linear = not pars['nonlinear']
        return(linear)
    elif name.find('linear') != -1:
        linear = True 
        return(linear)
    elif name.find('nonlin') != -1:
        linear = False
        return(linear)
        
def isUploaded(out_dir,runs):
    inDb = runs.find({ "run_collection_name": out_dir })
    for run in inDb:
        runIn = run["run_collection_name"]
        return(runIn == out_dir)

def update_mongo(out_dir,runs_coll):
    fields = ["user", "run_collection_name" ,"run_suffix" ,"keywords", "confidence",                        "codemods_id", "submitcmd_id", "parameters_id", "eqdisk_id", "efit_id", "autopar_id", "energy_id", "nrg_id", "omega_id", "scanlog_id", "scaninfo_id", "Qes", "ky", "kx", "omt", "omn", "gamma"]
    update_null = input('Would you like to update only "Null" fields?  (y/n)')
    for field in fields:
        if update_null == 'y' or update_null == 'Y':
            key = 'None'
        key = 'None'
        runs_coll.find({field: key })
        

def get_object_ids(out_dir):
        #generate file lists for relevant files
        autopar_files = get_file_list(out_dir, 'autopar')
        codemods_files = get_file_list(out_dir, 'codemods')
        energy_files = get_file_list(out_dir, 'energy')
        geneerr_files = get_file_list(out_dir, 'geneerr')
        nrg_files = get_file_list(out_dir, 'nrg')
        omega_files = get_file_list(out_dir, 'omega')
        parameters_files = get_file_list(out_dir, 'parameters')
        s_alpha_files = get_file_list(out_dir, 's_alpha')
        scanlog_files = get_file_list(out_dir, 'scan.log')
        scan_info_files = get_file_list(out_dir, 'scan_info.dat')

        #list of all GENE output files to be uploaded
        output_files = autopar_files + codemods_files + energy_files + geneerr_files +            nrg_files +  omega_files + parameters_files + s_alpha_files + scanlog_files +            scan_info_files
        
        #upload all GENE output files to database
        object_ids = [] 
        for files in output_files:
            object_ids.append(gridfs_put(files)) 
                
        #map from_output files to object ids
        for i in range(0,len(output_files)):
            object_ids[i] = object_ids[i] + '    ' +  output_files[i]
        return object_ids

def upload_to_mongo(out_dir, user, linear, confidence, input_heat, keywords):
    #initialize files dictionary
    files_dict =  {'scan_id': 'None', 
                   'scanlog_id': 'None', 
                   'scaninfo_id': 'None', 
                   'codemods_id': 'None', 
                   'submit_id': 'None', 
                   'parameters_id': 'None', 
                   'eqdisk_id': 'None', 
                   'efit_id': 'None', 
                   'autopar_id': 'None', 
                   'energy_id': 'None', 
                   'nrg_id': 'None', 
                   'omega_id': 'None', 
                   's_alpha_id': 'None'
                  }

    #for linear runs
    if linear:
        #connect to 'ETG' database, 'LinearRuns' collection
        runs_coll = MongoClient().ETG.LinearRuns
        if isUploaded(out_dir, runs_coll):
            update = input('Folder already uploaded, update the the runs? (y/n)')
            if update == 'y' or update == 'Y':
                update_mongo(out_dir, runs_coll)
        else:
            scan_info(out_dir)
            object_ids = get_object_ids(out_dir)            
            suffixes = get_suffixes(out_dir)
            for suffix in suffixes:
                for line in object_ids:
                    if line.find('codemods' + suffix) != -1:
                        files_dict['codemods_id'] = line.split()[0]
                    if line.find('parameters' + suffix) != -1:
                        files_dict['parameters_id'] = line.split()[0]
                    if line.find('autopar' + suffix) != -1:
                        files_dict['autopar_id'] = line.split()[0]
                    if line.find('energy' + suffix) != -1:
                        files_dict['energy_id'] = line.split()[0]
                    if line.find('nrg' + suffix) != -1:
                        files_dict['nrg_id'] = line.split()[0]
                    if line.find('omega' + suffix) != -1:
                        files_dict['omega_id'] = line.split()[0]
                    if line.find('scan.log') != -1:
                        files_dict['scanlog_id'] = line.split()[0]
                    if line.find('scan_info.dat') != -1:
                        files_dict['scaninfo_id'] = line.split()[0]
                    if line.find('s_alpha' + suffix) != -1:
                        files_dict['s_alpha_id'] = line.split()[0]
        
                #find relevant parameters from in/output
                gamma = find_omega(out_dir + '/omega' + suffix)[0]
                omega = find_omega(out_dir + '/omega' + suffix)[1]
                params = find_params(out_dir + '/parameters' + suffix)
                kx = params[0]
                ky = params[1]
                omn = params[2]
                omt = params[3]
                
                
                #document format for linear runs  
                run_data_dict = {"user": user,
                                 "run_collection_name": out_dir,
                                 "run_suffix": '' + suffix,
                                 "keywords": keywords,
                                 "confidence": confidence,
                                 "gamma": gamma,
                                 "omega": omega,
                                 "ky": ky,
                                 "kx": kx,
                                 "omt": omt,
                                 "omn": omn
                                 }
                
                run_data =  {**run_data_dict, **files_dict}
                
                #insert run_data into database
                runs_coll.insert_one(run_data).inserted_id
            print('Run collection \'' + out_dir + '\' uploaded succesfully.')
    
    #for nonlinear runs
    if not linear:
        #connect to 'ETG' database, 'NonlinRuns' collection
        runs_coll = MongoClient().ETG.NonlinRuns
        if isUploaded(out_dir,runs_coll):
            print('Folder already uploaded')
        else:
            object_ids = get_object_ids(out_dir)
            
            suffixes = get_suffixes(out_dir)
            for suffix in suffixes:
                for line in object_ids:                   
                    if line.find('codemods' + suffix) != -1:
                        files_dict['codemods_id'] = line.split()[0]
                    if line.find('parameters' + suffix) != -1:
                        files_dict['parameters_id'] = line.split()[0]
                    if line.find('autopar' + suffix) != -1:
                        files_dict['autopar_id'] = line.split()[0]
                    if line.find('energy' + suffix) != -1:
                        files_dict['energy_id'] = line.split()[0]
                    if line.find('nrg' + suffix) != -1:
                        files_dict['nrg_id'] = line.split()[0]
                    if line.find('omega' + suffix) != -1:
                        files_dict['omega_id'] = line.split()[0]
                    if line.find('scan.log') != -1:
                        files_dict['scanlog_id'] = line.split()[0]
                    if line.find('scan_info.dat') != -1:
                        files_dict['scaninfo_id'] = line.split()[0]
                    if line.find('s_alpha' + suffix) != -1:
                        files_dict['s_alpha_id'] = line.split()[0]
                        
                #find relevant parameters from in/output
                Qes = get_Qes(out_dir, suffix)
                params = find_params(out_dir + '/parameters' + suffix)
                kx = params[0]
                ky = params[1]
                omn = params[2]  # check n_spec for suffix
                omt = params[3]
                
                #document format for nonlinear runs
                run_data_dict = {"user": user,
                                 "run_collection_name": out_dir,
                                 "run_suffix": '' + suffix,
                                 "keywords": keywords,
                                 "confidence": confidence,                       
                                 "Qes" : Qes,
                                 "ky" : ky,
                                 "kx" : kx,
                                 "omt" : omt,
                                 "omn" : omn                            
                                }
                
                run_data =  {**run_data_dict, **files_dict}
                
                #insert run_data into database
                runs_coll.insert_one(run_data).inserted_id
            print('Run collection \'' + out_dir + '\' uploaded succesfully.')

def upload_big(out_dir, linear):
    
    field_files = get_file_list(out_dir, 'field')
    mom_files = get_file_list(out_dir, 'mom')
    vsp_files = get_file_list(out_dir, 'vsp')
    
    output_files = field_files + mom_files + vsp_files

    object_ids = [] 
    for files in output_files:
        object_ids.append(gridfs_put(files)) 
            
    #map from_output files to object ids
    for i in range(0,len(output_files)):
        object_ids[i] = object_ids[i] + '    ' +  output_files[i]
    
    if not linear:
        object_ids = get_object_ids(out_dir)
        field_id = 'None'
        mom_id = 'None'
        vsp_id = 'None'
        suffixes = get_suffixes(out_dir)
        for suffix in suffixes:
            for line in object_ids:
                if line.find('field' + suffix) != -1:
                    field_id = line.split()[0]
                if line.find('mom' + suffix) != -1:
                    mom_id = line.split()[0]
                if line.find('vsp' + suffix) != -1:
                    vsp_id = line.split()[0]        

            run_data = {"field_id": field_id,
                        "mom_id": mom_id,
                        "vsp_id": vsp_id,
                       }
            
            runs.update(run_data).inserted_id