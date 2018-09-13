#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Base script to handle uploading of GENE folders to database

@author: Austin Blackmon
"""


from mongo_file_handling_v1 import *
import numpy as np

#############################
user = 'A. Blackmon'
out_dir = '78697.51105_hager_nonlin_r95'
linear = True
scan = True
confidence = '5'
keywords = 'ETG, pedestal, GENE'
input_heat = ''
#############################

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

#format 'energy' and 'nrg' files
for files in energy_files:
    format_energy(files)

for files in nrg_files:
    format_nrg(files)

#replace energy and nrg file lists with formatted file lists
formatted_energy_files = get_file_list(out_dir, 'formatted_energy')
formatted_nrg_files = get_file_list(out_dir, 'formatted_nrg')

#list of all GENE output files to be uploaded
output_files = autopar_files + codemods_files + energy_files + \
    field_files + geneerr_files + mom_files + nrg_files + \
    omega_files + parameters_files + s_alpha_files + scanlog_files +\
    vsp_files

#upload all GENE output files to database
object_ids = [] 

for files in output_files:
    object_ids.append(gridfs_put(files))
        
object_id = np.chararray((len(output_files) - 1,len(output_files) - 1))
for i in range(0,len(output_files) - 1):
    object_id[0,i] = output_files[i]

#connect to 'ETG' database, 'runs' collection
runs = MongoClient().ETG.Runs

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
scanlog_id = 'None'
diagnostics_id = 'None'

if scan:
    numruns = len(codemods_files)
    for i in range(0,numruns - 1):
        for line in object_ids:
            if line.find('codemods') != -1 \
            and line.endswith(str(i).zfill(4)):
                codemods_id = line.split('    ')[0]
            if line.find('parameters') != -1 \
            and line.endswith(str(i).zfill(4)):
                parameters_id = line.split('    ')[0]
            if line.find('autopar') != -1 \
            and line.endswith(str(i).zfill(4)):
                autopar_id = line.split('    ')[0]
            if line.find('energy') != -1 \
            and line.endswith(str(i).zfill(4)):
                energy_id = line.split('    ')[0]
            if line.find('formatted_energy') != -1 \
            and line.endswith(str(i).zfill(4)):
                energy_id = line.split('      ')[0]
            if line.find('nrg') != -1 \
            and line.endswith(str(i).zfill(4)):
                nrg_id = line.split('    ')[0]
            if line.find('formatted_nrg') != -1 \
            and line.endswith(str(i).zfill(4)):
                nrg_id = line.split('    ')[0]
            if line.find('omega') != -1 \
            and line.endswith(str(i).zfill(4)):
                omega_id = line.split('    ')[0]
            if line.find('scan.log') != -1:
                scanlog_id = line.split('    ')[0]
            if line.find('s_alpha') != -1 \
            and line.endswith(str(i).zfill(4)):
                s_alpha_id = line.split('    ')[0]
                
        #document format for 'runs' collection  
            if codemods_id != 'None':              
                run_data = {"user": user,
                            "runname": out_dir + '_' + str(i).zfill(4),
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
                            "vsp": vsp_id,
                            "omega":omega_id,
                            "scanlog": scanlog_id,
                            "diagnostics": diagnostics_id
                           }
                runs.insert_one(run_data).inserted_id


if not linear:
    numruns = len(codemods_files)
    for i in range(0,numruns - 1):
        for line in object_ids:
            if line.find('codemods') != -1 \
            and not line.endswith(filetypes):           
                codemods_id = line.split('    ')[0]
            if line.find('parameters') != -1 \
            and not line.endswith(filetypes):
                parameters_id = line.split('    ')[0]
            if line.find('autopar') != -1 \
            and not line.endswith(filetypes):
                autopar_id = line.split('    ')[0]
            if line.find('energy') != -1 \
            and not line.endswith(filetypes):
                energy_id = line.split('    ')[0]
            if line.find('field') != -1 \
            and not line.endswith(filetypes):
                field_id = line.split('    ')[0]
            if line.find('mom') != -1 \
            and not line.endswith(filetypes):
                mom_id = line.split('    ')[0]
            if line.find('nrg') != -1 \
            and not line.endswith(filetypes):
                nrg_id = line.split('    ')[0]
            if line.find('vsp') != -1 \
            and not line.endswith(filetypes):
                vsp_id = line.split('    ')[0]

        #document format for 'runs' collection
            if codemods_id != 'None':
                run_data = {"user": user,
                            "filename": out_dir + '_' + str(i+1).zfill(4),
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
                            "vsp": vsp_id,
                            "omega":omega_id,
                            "scanlog": scanlog_id,
                            "diagnostics": diagnostics_id
                           }
                
                runs.insert_one(run_data).inserted_id