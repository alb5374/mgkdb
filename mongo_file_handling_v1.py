#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
File handling script for formatting output files, getting file lists, and
reading and writing to database containing:
    format_energy(filepath):            input filepath, return formatted GENE energy 
                                          output formatted energy file (header removed)
    format_nrg(filepath):               input filepath, return formatted GENE nrg output
                                          files ### currently only for n_spec=1 ###
    get_file_list(out_dir,begin):       input GENE output directory and base filepath 
                                          (nrg, energy, etc), return full list of files 
                                          in directory
    get_formatted_list(out_dir,begin):  input GENE output directory and base filepath 
                                          (nrg, energy) of formatted files, return full 
                                          list of formatted files in directory
    gridfs_put(filepath):               input filepath, upload  file to database, and 
                                          return object_id of uploaded file
    gridfs_read(db_file):              input database filename, return contents of file
        

@author: Austin Blackmon
"""

from pymongo import MongoClient
import gridfs
import os
import numpy as np

def format_energy(filepath):
    #open files to fix columns
    energyFile = open(filepath, 'r')
    filepath_split = filepath.split('\\')
    energyFormattedFile = open(filepath_split[0] + '\\formatted_' + filepath_split[1],'w')
    
    #edit energy(...) file and write to energy(...)_formatted
    with energyFile as f:
        for count, line in enumerate(f, start=1):

            #grab time from odd number lines
            if count > 15:
                energyFormattedFile.write(line)

    #close files to fix columns
    energyFile.close()
    energyFormattedFile.close()

def format_nrg(filepath):
    #open files to fix columns
    nrgFile = open(filepath, 'r')
    filepath_split = filepath.split('\\')
    nrgFormattedFile = open(filepath_split[0] + '\\formatted_' + filepath_split[-1],'w')
    
    #edit nrg(...) file and write to nrg(...)_formatted
    with nrgFile as f:
        for count, line in enumerate(f, start=1):

            #grab time from odd number lines
            if count % 2 == 1:
                time = line.replace('\n','')

            #add time to lines of data on even numbered lines
            else:
                line_data = line
                line = time + line_data
                line = line.replace('     ', '')
                nrgFormattedFile.write(line)

    #close files to fix columns
    nrgFile.close()
    nrgFormattedFile.close()

def get_file_list(out_dir,begin):
    files_list = []
    count = 0
    filetypes = ('.ps','.png')
    
    #scan files in GENE output directory
    for dirpath, dirnames, files in os.walk(out_dir):
        for name in files:
            if name.lower().startswith(begin) and name.find('in_par') == -1  \
            and name.endswith(filetypes) == False:
                files = os.path.join(dirpath, name)
                files_list.append(files)
                count = count + 1
    return files_list
  
def get_formatted_list(out_dir,begin):
    files_list = []
    count = 0
    
    #scan files
    for dirpath, dirnames, files in os.walk(out_dir):
        for name in files:
            if name.lower().startswith(begin):
                files = os.path.join(dirpath, name)
                files_list.append(files)
                count = count + 1
    return files_list          

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
    

  
    
    #loop through all runs, only use 'etg' keyword, check what run is 
    #    closest to base run based on certain parameter, least squared fit