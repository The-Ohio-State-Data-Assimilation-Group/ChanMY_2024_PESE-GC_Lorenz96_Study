#!/bin/python

# =============================================================================
# SCRIPT TO SELECT DIRECTORY WITH THE BEST FORECAST RMSES
# =============================================================================

# Load useful modules
import numpy as np
from sys import argv


# Load list of all directories to compare prior RMSEs
roi_list = argv[1]
roi_list = roi_list.split(",")[:-1]
dir_list = ['halfROI'+roi_string for roi_string in roi_list ]

# Generate list to hold forecast RMSE
fcst_rmse_list = []


# Read in the prior RMSEs from each ROI tuning directory
for dirname in dir_list:

    # Load text file & parse prior rmse
    f = open( "%s/useful_stats.txt" % dirname, "r" )
    all_lines = f.readlines()
    prior_rmse = float( all_lines[0].split()[1] )
    fcst_rmse_list.append( prior_rmse )
    f.close()

# --- End of loop over each ROI tuning directory

    


# Save directory names and the RMSE values into a text file
f = open('ROI_RMSE_values.txt','w')
string = ""
for dd in range(len(dir_list)):
    string += "ROI: %9.3f    RMSE: %f\n" % ( float(roi_list[dd]), fcst_rmse_list[dd] )
# --- End of loop over directory list
f.write( string )
f.close()




# Determine best-performing directory and print out
best_rmse_ind = np.argmin( fcst_rmse_list )
print(dir_list[best_rmse_ind])

