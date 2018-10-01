# -*- coding: utf-8 -*-
"""
Main script to handle uploading GENE runs to the MGK database
Required fields:    user
                    output_folder
                    multiple_runs (True or False)
                    linear (True or False)
                    
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

output_folder = '.'
multiple_runs = True
#linear = True   ######## CHECK PARAMS FILE INSTEAD ########

if not multiple_runs:
    confidence = '5'  ### 1-10, 1: little confidence, 10: well checked ###
else:
    confidence = 'None'  ### Set if same for all runs ###
    
#if linear:
#    lin = 'linear'
#else:
#    lin = 'nonlin'
#    
input_heat = 'None'
    
### enter any relevant keywords, i.e., ETG, ITG, pedestal, core ###
keywords = 'ETG, pedestal, GENE, ' #+ lin

#######################################################################

if multiple_runs:
    folder_list = []
    for dirpath, dirnames, files in os.walk(output_folder):
        
        for count, name in enumerate(dirnames, start=0):

            folder = os.path.join(name)
            folder_list.append(folder)
            
            linear = isLinear(name)
                  
            for folder in folder_list:
                upload_to_mongo(folder, user, linear, confidence, input_heat, keywords)
else:
    upload_to_mongo(output_folder, user, linear, confidence, input_heat, keywords)