#!/bin/python


'''
    SCRIPT TO PRINT OUT THE ASYMPTOTIC FORECAST RMSES & ROIS OF A PARTICULAR SETUP
    
    COMMAND LINE INPUTS:
        1) NUMBER OF TRIALS
        2) PATH TO DIRECTORY CONTAINING ALL TRIALS

    THIS SCRIPT ASSUMES THAT EACH TRIAL HAS A SUB-DIRECTORY "long" CONTAINING LONG DA EXPERIMENT RUNS

    THE ASYMPTOTIC RMSE AND ROI FOR EVERY TRIAL WILL BE PRINTED OUT.
'''


import numpy as np
from netCDF4 import Dataset as ncopen
from sys import argv
from os.path import exists


# Read user inputs from command line
num_trials = int(argv[1])
path_to_trial_parent_dir = argv[2]


# Generate list of long DA experiment directories
list_trial_dirs = [ '%s/Trial%03d/long' % (path_to_trial_parent_dir, itrial) for itrial in range( 1, num_trials+1) ]



# Compute forecast RMSEs
# ----------------------
# Init array to hold results
forecast_rmse_list = np.zeros( num_trials , dtype='f')

# Loop over all trial directories
for itrial in range(num_trials):

    mydir=list_trial_dirs[itrial]

    # Determine filenames
    fcst_fname = mydir+"/forecast.nc"
    true_fname = mydir+"/true_state.nc"

    # If either files do not exist, skip over this trial
    if ( ( not exists(fcst_fname) ) or ( not exists(true_fname) ) ):
        print( fcst_fname )
        forecast_rmse_list[itrial] = -1.111
        continue

    # If both files exist, load them
    fcst_file = ncopen( fcst_fname, 'r')
    true_file = ncopen( true_fname, 'r')


    # Handling exception case where DA experiment broke down halfway
    fcst_timedim = fcst_file.variables['state'][:].shape[0]
    true_timedim = true_file.variables['state'][:].shape[0]
    if ( fcst_timedim != true_timedim ):
        forecast_rmse_list[itrial] = -2.222
        continue
        

    # If none of the exception cases are triggered, compute the RMSE
    fcst_rmse = np.mean( 
        np.sqrt(
            np.mean(
                (
                    np.mean( fcst_file.variables['state'][500:,:,:], axis=1 )
                    - true_file.variables['state'][500:,0,:]
                )**2,
                axis=-1
            )
        ), axis=0 )

    # Handling exception case where weird values are obtained
    if ( fcst_rmse == np.nan or fcst_rmse == np.inf or np.ma.is_masked(fcst_rmse) or fcst_rmse > 100):
        fcst_rmse = -8.888

    forecast_rmse_list[itrial] = fcst_rmse

    # Free-up file handles
    fcst_file.close()
    true_file.close()

# --- End of loop over trials



# Print out forecast RMSEs
rmse_string = ""
for itrial in range( num_trials ):
    substring = "  %10.4f  " % forecast_rmse_list[itrial]
    rmse_string += substring
print( rmse_string )