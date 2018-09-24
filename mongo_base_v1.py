#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Base script to handle uploading of GENE folders to database

@author: Austin Blackmon
"""


from mongo_post_processing_v1 import *
from mongo_file_handling_v1 import *
from ParIO import * 
import numpy as np


#######################################################################

user = 'A. Blackmon'

out_dir = '78697.51105_hager_nonlin_r95'
linear = 0

#out_dir = 'salpha_lin_kyscan'
#linear = 1

confidence = '5'  ### 1-10, 1: little confidence, 10: well checked ###
input_heat = ''
if linear:
    lin = 'linear'
else:
    lin = 'nonlin'
    
# Please enter any relevant keywords, i.e., ETG, ITG, pedestal, core #
keywords = 'ETG, pedestal, GENE, ' + lin

######################################################################

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

#format 'energy' and 'nrg' files and get formatted file lists
for files in energy_files:
    format_energy(files)
for files in nrg_files:
    format_nrg(files)
formatted_energy_files = get_file_list(out_dir, 'formatted_energy')
formatted_nrg_files = get_file_list(out_dir, 'formatted_nrg')

#list of all GENE output files to be uploaded
output_files = autopar_files + codemods_files + energy_files + \
    formatted_energy_files + field_files + geneerr_files + \
    mom_files + nrg_files + formatted_nrg_files + omega_files + \
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
### USE DICTIONARY ###


#for linear runs
if linear:
    extensions = get_extensions(out_dir)
    for extension in extensions:
        for line in object_ids:
            if line.find('codemods_' + extension) != -1:
                codemods_id = line.split()[0]
            if line.find('parameters_' + extension) != -1:
                parameters_id = line.split()[0]
            if line.find('autopar_' + extension) != -1:
                autopar_id = line.split()[0]
            if line.find('energy_' + extension) != -1 \
            and line.find('formatted') == -1:
                energy_id = line.split()[0]
            if line.find('formatted_energy_' + extension) != -1:
                formatted_energy_id = line.split()[0]
            if line.find('nrg_' + extension) != -1 \
            and line.find('formatted') == -1 :
                nrg_id = line.split()[0]
            if line.find('formatted_nrg_' + extension) != -1:
                formatted_nrg_id = line.split()[0]
            if line.find('omega_' + extension) != -1:
                omega_id = line.split()[0]
            if line.find('scan.log') != -1:
                scanlog_id = line.split()[0]
            if line.find('s_alpha_' + extension) != -1:
                s_alpha_id = line.split()[0]

        #find relevant parameters from in/output
        gamma = find_gamma(out_dir + '\\omega_' + extension)
        params = find_params(out_dir + '\\parameters_' + extension)
        kx = params[0]
        ky = params[1]
        omn = params[2]  # check n_spec for suffix
        omt = params[3]
        
        
        #document format for linear runs  
        run_data = {"user": user,
                    "run_collection_name": out_dir,
                    "run_extension": out_dir + '_' + extension,
                    "keywords": keywords,
                    "codemods": codemods_id,
                    "submitcmd": submit_id,
                    "confidence": confidence,
                    "parameters": parameters_id,
                    "eqdisk": eqdisk_id,
                    "efit": efit_id,
                    "autopar": autopar_id,
                    "energy": energy_id,
                    "formatted_energy": formatted_energy_id,
                    "field": field_id,
                    "mom": mom_id,
                    "nrg": nrg_id,
                    "formatted_nrg": formatted_nrg_id,
                    "omega":omega_id,
                    "scanlog": scanlog_id,
                    "vsp": vsp_id,
                    "gamma": gamma,
                    "ky" : ky,
                    "kx" : kx,
                    "omt" : omt,
                    "omn" : omn
                    }
        #insert run_data into database
        runs.insert_one(run_data).inserted_id
                
                
                
#if linear:
#    extensions = get_extensions(out_dir)
#    for extension in extensions:
#        for line in object_ids:

########### WHY IS ONLY 1 NONLIN RUN BEING UPLOADED?!?!?!?! #############

#for nonlinear runs
if not linear:
    extensions = get_extensions(out_dir)
    for extension in extensions:
        for line in object_ids:
            if line.find('codemods_' + extension) != -1:
                codemods_id = line.split()[0]
            if line.find('parameters_' + extension) != -1:
                parameters_id = line.split()[0]
            if line.find('autopar_' + extension) != -1:
                autopar_id = line.split()[0]
            if line.find('energy_' + extension) != -1 \
            and line.find('formatted') == -1:
                energy_id = line.split()[0]
            if line.find('formatted_energy_' + extension) != -1:
                formatted_energy_id = line.split()[0]
            if line.find('nrg_' + extension) != -1 \
            and line.find('formatted') == -1 :
                nrg_id = line.split()[0]
            if line.find('formatted_nrg_' + extension) != -1:
                formatted_nrg_id = line.split()[0]
            if line.find('omega_' + extension) != -1:
                omega_id = line.split()[0]
            if line.find('scan.log') != -1:
                scanlog_id = line.split()[0]
            if line.find('s_alpha_' + extension) != -1:
                s_alpha_id = line.split()[0]
            if line.find('field_' + extension) != -1:
                field_id = line.split()[0]
            if line.find('mom_' + extension) != -1:
                mom_id = line.split()[0]
            if line.find('vsp_' + extension) != -1:
                vsp_id = line.split()[0]
                
        #find relevant parameters from in/output
        Qes = get_Qes(out_dir, extension)
        params = find_params(out_dir + '\\parameters_' + extension)
        kx = params[0]
        ky = params[1]
        omn = params[2]  # check n_spec for suffix
        omt = params[3]
        
        #document format for nonlinear runs
        run_data = {"user": user,
                    "run_collection_name": out_dir,
                    "run_extension": '_' + extension,
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
        print(extension,line)
########################

# HOW TO INTEGRATE ADIOS??!?!

#########################
                