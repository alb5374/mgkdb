# -*- coding: utf-8 -*-
"""
Main script to handle uploading GENE runs to the MGK database
Required fields:    user
                    output_folder
                    multiple_runs (True or False)
                    
Optional fields:    confidence
                    input_heat
                    keywords
                    
@author: Austin Blackmon
"""

from mgk_file_handling import *
from ParIO import *
import os


########################################################################

user = 'A. Blackmon'

output_folder = '.'     ### Set as '.' for current directory ###
multiple_runs = True    ### Automate scanning through a directory of numerous runs ###

if not multiple_runs:
    confidence = '5'     ### 1-10, 1: little confidence, 10: well checked ###
else:
    confidence = 'None'  ### Set if  same for all runs, else set as 'None' ###
    
input_heat = 'None'      ### Set if input heat is known, else set as 'None' ###
    
### enter any relevant keywords, i.e., ETG, ITG, pedestal, core ###
keywords = 'ETG, pedestal, GENE, '

#######################################################################

#scan through a directory for more than one run
if multiple_runs:    
    #scan through directory for run directories
    for dirpath, dirnames, files in os.walk(output_folder):
        for count, name in enumerate(dirnames, start=0):
            #make   list of directories
            folder = os.path.join(name)
            #check if run is linear or nonlinear
            linear = isLinear(name)
            if linear:
                lin = 'linear'
            else:
                lin = 'nonlin'
            #add linear/nonlin to keywords
            keywords_lin = keywords + lin
            
            #send run list to upload_to_mongo to be uploaded
            upload_to_mongo(folder, user, linear, confidence, input_heat, keywords_lin)

#submit a single run
else: 
    #check if run is linear or nonlinear
    linear = isLinear(output_folder)
    if linear:
        lin = 'linear'
    else:
        lin = 'nonlin'
    #add linear/nonlin to keywords
    keywords_lin = keywords + lin
    
    #send run to upload_to_mongo to be uploaded
    upload_to_mongo(output_folder, user, linear, confidence, input_heat, keywords_lin)