#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
File handling script for formatting output files, getting file lists, and
reading and writing to database containing:
    get_file_list(out_dir,begin):       input GENE output directory and base filepath 
                                           (nrg, energy, etc), return full list of files 
                                           in directory
    get_extensions(out_dir):            input GENE output directory, return list of run 
                                           extensions in the directory
    gridfs_put(filepath):               input filepath, upload  file to database, and 
                                           return object_id of uploaded file
    gridfs_read(db_file):               input database filename, return contents of file
        

@author: Austin Blackmon
"""

from pymongo import MongoClient
import gridfs
import os
import numpy as np
from ParIO import * 

def get_file_list(out_dir,begin):
    files_list = []
    
    #unwanted filetype extensions for general list
    bad_ext = ('.ps','.png')
    
    #scan files in GENE output directory, ignoring files in '/in_par', and return list
    for dirpath, dirnames, files in os.walk(out_dir):
        for count, name in enumerate(files, start=0):
            if name.lower().startswith(begin) and name.find('in_par') == -1  \
            and name.endswith(bad_ext) == False:
                files = os.path.join(dirpath, name)
                files_list.append(files)
                
    return files_list     


def get_extensions(out_dir):
    extensions = []
    
    #scan files in GENE output directory, find all run extensions, return as list
    for dirpath, dirnames, files in os.walk(out_dir):
        for count, name in enumerate(files, start=0):
            if name.lower().startswith('codemods'):
                extension = name.split('_',1)[1]
                extensions.append(extension)
    return extensions

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