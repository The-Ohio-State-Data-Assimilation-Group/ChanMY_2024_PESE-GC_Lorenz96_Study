#!/bin/python

# =============================================================================
# SCRIPT TO SELECT THE ROI WITH THE BEST PERFORMANCE AND LOAD RELEVANT INFO
# =============================================================================

# Load useful modules
import numpy as np
from sys import argv
from os.path import exists as check_exist


# Load list of all directories to compare prior RMSEs
roi_list = argv[1]
roi_list = roi_list.split(",")[:-1]
dir_list = ['halfROI'+roi_string for roi_string in roi_list ]
nROI = len(dir_list)

# Generate arrays to hold important metrics
analysis_rmse_arr = np.zeros(nROI)
forecast_cr_arr = np.zeros(nROI)
analysis_bias_arr = np.zeros(nROI)
analysis_stability_arr = np.zeros(nROI)


# Read in the analysis RMSEs from each ROI tuning directory
for dd in range(nROI):

    # Read directory name
    dirname = dir_list[dd]

    # Generate textfile name
    txt_fname = "%s/useful_stats.txt" % dirname

    # if textfile exists, load stats. Otherwise, set to nan
    if( check_exist( txt_fname ) ):
        f = open( txt_fname, "r" )
        all_lines = f.readlines()
        f.close()

        # Load posterior stats
        analysis_line_splitted = all_lines[1].split()
        analysis_rmse_arr[dd] = float( analysis_line_splitted[1] )
        analysis_bias_arr[dd] = float( analysis_line_splitted[5] )

        # Load rms(dxdt)
        analysis_stability_arr[dd] = float( all_lines[2].split()[1] )        

        # Load forecasted CR
        forecast_cr_arr[dd] = float( all_lines[1].split()[3] )
    
    else:

        # Set values to NaN if file does not exist
        analysis_rmse_arr[dd] = np.nan
        forecast_cr_arr[dd] =   np.nan
        analysis_bias_arr[dd] = np.nan
        analysis_stability_arr[dd] = np.nan

    # ---- End of file-existence checker
   

# --- End of loop over each ROI tuning directory

    

# Overwrite bad RMSE values
mask = (
    (
        np.isnan( analysis_rmse_arr ) + np.isinf( analysis_rmse_arr ) + ( analysis_rmse_arr > 9999)
    ) > 0
)
analysis_rmse_arr[mask] = 9999



# Determine best-performing ROI 
best_ind = np.argmin( analysis_rmse_arr )


# Print out useful information
best_rmse = analysis_rmse_arr[best_ind]
best_cr = forecast_cr_arr[best_ind]
best_bias = analysis_bias_arr[best_ind]
best_stability = analysis_stability_arr[best_ind]
best_roi = float(dir_list[best_ind].replace("halfROI",""))

print( "ROI:%f,RMSE:%f,CR:%f,BIAS:%f,DXDT:%f" % (best_roi, best_rmse, best_cr, best_bias, best_stability ) )

