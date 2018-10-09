# **MGKDB**
MGKDB is a tool for storing, accessing, and learning from a wealth of fusion results.  The main focus of MGKDB will be in the pedestal region of tokamaks, with a large portion of data coming from simulation results.  Currently, GENE is the only code compatible with the database, though compatibility with more codes will be added in the future.

MGKDB is written for Python 3.6, using Anaconda 1.6.  As of yet, no other versions of Python or Anaconda have been tested.  Any other scritps necessary for running MGKDB are included in this repository.

## **Instructions**
1. Clone the repository.
2. In mgk_uploader.py, you will find the code listed below, sectioned off by ####...####.  This section is where you edit variables to fit your needs.  
	* ** *Required* ** fields are: ```user```, ```output_folder```, and ```multiple_runs```.  
	* *Desired* fields are: ```input_heat```, ```confidence```, and ```keywords```.  
		* ```input_heat```, in MW, is for simulations based on experiments with known input heat
		* ```confidence``` will allow MGKDB to guage the checks that went into setting up your simulation.  Low values are for simulations that were quickly thrown together, with little or no prior checks performed.  High values are for simulations for which a wide array of numerical and physical checks were performed
		* ```keywords``` are very helpful to provide MGKDB with metadata and allows for smart searching through the database.  Please include as many as possible!
3. MGKDB will automate finding important parameters and upload them to the database.  A message will display if your run was uploaded successfully.  

*Thank you for contributing!*
```python 
########################################################################

### REQUIRED ###
user = 'Your Name'
output_folder = '.'     ### Set as '.' for current directory ###
multiple_runs = True    ### Automate scanning through a directory of numerous runs ###
#################

### DESIRED ###
if not multiple_runs:
    confidence = ''     ### 1-10, 1: little confidence, 10: well checked ###
else:
    confidence = 'None'  ### Set if  same for all runs, else set as 'None' ###
input_heat = 'None'      ### Set if input heat is known, else set as 'None' ###
keywords = 'ETG, pedestal, GENE'  ### enter any relevant keywords, i.e., ETG, ITG, pedestal, core ###
###############

########################################################################
```
## Future of MGKDB
* Integration into GENE diagnostic tool
* Intuitive interface for user access
* 



	